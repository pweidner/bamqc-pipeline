#!/usr/bin/env python3
import argparse, gzip, sys, os
import numpy as np
import pandas as pd

def read_counts(path):
    """Read per-window counts.

    Supported formats:
    - legacy: chrom, start, end, counts
    - strand-aware: chrom, start, end, watson_count, crick_count, counts
    """
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", header=None, comment="#")
    if df.shape[1] < 4:
        raise ValueError(f"{path}: expected at least 4 columns (chrom, start, end, counts). Got {df.shape[1]}")

    if df.shape[1] == 6:
        df.columns = ["chrom", "start", "end", "watson_count", "crick_count", "counts"]
    else:
        df.columns = [*(
            ["chrom", "start", "end"] + [f"c{i}" for i in range(3, df.shape[1] - 1)]
        ), "counts"]

    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df["counts"] = pd.to_numeric(df["counts"], errors="coerce").fillna(0)
    for col in ("watson_count", "crick_count"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    df["binsize"] = (df["end"] - df["start"]).astype(np.int64)
    return df

def entropy(counts: np.ndarray) -> float:
    tot = counts.sum()
    if tot <= 0: return np.nan
    p = counts / tot
    # ignore zeros via masked array to avoid log(0)
    p = p[p > 0]
    return float(-(p * np.log(p)).sum())

def spikiness(counts: np.ndarray) -> float:
    tot = counts.sum()
    if tot <= 0: return np.nan
    return float(np.abs(np.diff(counts)).sum() / tot)

def gini_index(x: np.ndarray) -> float:
    # standard Gini on non-negative
    x = np.asarray(x, dtype=float)
    if x.size == 0: return np.nan
    if np.any(x < 0): 
        # shift if negative (shouldn't happen for counts)
        x = x - x.min()
    if x.sum() == 0: return 0.0
    x = np.sort(x)
    n = x.size
    cum = np.cumsum(x)
    # Gini = 1 + 1/n - 2 * sum_i ((n+1-i)*x_i)/(n*sum x)
    g = 1 + 1.0/n - 2.0 * (cum.sum() / (n * x.sum()))
    return float(g)

def pct_at_least(counts: np.ndarray, thresh: int) -> float:
    if counts.size == 0: return np.nan
    return float((counts >= thresh).mean())

def fold80_penalty(counts: np.ndarray) -> float:
    """
    Approximate Fold-80 penalty used in exome QC:
    mean_coverage / coverage_at_20th_percentile.
    Larger => worse uniformity. (Picard reports similar metric.)
    """
    if counts.size == 0: return np.nan
    mean_cov = counts.mean()
    if mean_cov == 0: return np.nan
    cov20 = np.percentile(counts, 20.0)
    cov20 = max(cov20, 1e-12)
    return float(mean_cov / cov20)

def cv(counts: np.ndarray) -> float:
    if counts.size == 0: return np.nan
    mean = counts.mean()
    if mean == 0: return np.nan
    return float(counts.std(ddof=0) / mean)

def mad(counts: np.ndarray) -> float:
    if counts.size == 0: return np.nan
    med = np.median(counts)
    return float(np.median(np.abs(counts - med)))

def background_estimate(df: pd.DataFrame) -> float:
    """Estimate strand background from Watson/Crick-biased bins.

    This follows the breakpointR idea: classify bins by W-C balance, estimate
    Crick leakage in WW bins and Watson leakage in CC bins, then average the
    available estimates. Watson is the minus strand, Crick the plus strand.
    """
    if not {"watson_count", "crick_count"}.issubset(df.columns):
        return np.nan

    w = df["watson_count"].to_numpy(dtype=float)
    c = df["crick_count"].to_numpy(dtype=float)
    total = w + c
    informative = total >= 10
    if not np.any(informative):
        return np.nan

    ratio = np.full(total.shape, np.nan, dtype=float)
    ratio[informative] = (w[informative] - c[informative]) / total[informative]
    ww = informative & (ratio > 0.8)
    cc = informative & (ratio < -0.8)

    estimates = []
    if np.any(ww):
        estimates.append((c[ww].sum() + ww.sum()) / (w[ww].sum() + ww.sum()))
    if np.any(cc):
        estimates.append((w[cc].sum() + cc.sum()) / (c[cc].sum() + cc.sum()))

    return float(np.mean(estimates)) if estimates else np.nan


def read_gc_table(gc_path):
    """
    Expect a table aligned to the same windows with a GC column in [0,1] or [0,100].
    Accepts outputs like `bedtools nuc` (column 'GC') or custom tsv with 'gc'/'GC'.
    Must be sorted/identical order to counts, or provide keys to merge on.
    For robustness, we merge on chrom/start/end if those columns exist.
    """
    op = gzip.open if gc_path.endswith(".gz") else open
    df = pd.read_csv(op(gc_path, "rt"), sep="\t", comment="#")
    # Try to locate GC column
    gc_col = None
    for cand in ["GC","gc","gc_fraction","gc_content"]:
        if cand in df.columns:
            gc_col = cand
            break
    if gc_col is None:
        # bedtools nuc header often has '%GC' too
        for cand in ["%GC","pct_gc","PCT_GC"]:
            if cand in df.columns:
                gc_col = cand
                # convert to fraction later
                break
    if gc_col is None:
        raise ValueError(f"Could not find a GC column in {gc_path}.")
    # Normalize GC to 0..1
    gc = df[gc_col].astype(float).values
    if gc.max() > 1.0 + 1e-6:
        gc = gc / 100.0
    # If chrom/start/end present, keep them for merge; else assume same order.
    keep_cols = [c for c in ["chrom","start","end"] if c in df.columns]
    out = pd.DataFrame({"gc": gc})
    if keep_cols:
        out[keep_cols] = df[keep_cols].values
    return out

def main():
    ap = argparse.ArgumentParser(description="QC metrics from per-window counts")
    ap.add_argument("--counts", required=True, help="Per-window counts TSV/bed.gz (chrom,start,end,counts)")
    ap.add_argument("--sample", required=True, help="Library/sample name")
    ap.add_argument("--out",    required=True, help="Output TSV (one row)")
    # Coverage thresholds for PCT>=X (comma-separated integers)
    ap.add_argument("--pct-thresholds", default="1,10,30")
    # Optional extras (used if provided)
    ap.add_argument("--gc-table", help="TSV with GC per window (e.g. bedtools nuc). Must align with the same windows; can be merged via chrom/start/end if present.")
    ap.add_argument("--preseq", help="preseq lc_extrap output for this Library (optional complexity metrics)")
    args = ap.parse_args()

    df = read_counts(args.counts)
    counts = df["counts"].values.astype(float)
    binsize_mean = float(df["binsize"].mean()) if "binsize" in df else np.nan
    bin_background = background_estimate(df)

    # Basic summary
    total = float(counts.sum())
    avg_read_count = float(counts.mean())

    # Aneufinder-style metrics
    ent = entropy(counts)
    spike = spikiness(counts)

    # Variability metrics
    gini = gini_index(counts)
    coeff_var = cv(counts)
    abs_mad = mad(counts)
    sd_cov = float(counts.std(ddof=0))

    # PCT>=X coverage thresholds
    thresh = [int(x) for x in args.pct_thresholds.split(",") if x.strip() != ""]
    pct_metrics = { f"pct_ge_{t}x": pct_at_least(counts, t) for t in thresh }

    # Fold-80 penalty
    fold80 = fold80_penalty(counts)

    # GC-bias (optional): Pearson r between counts (or normalized) and GC fraction
    gc_r = np.nan
    if args.gc_table:
        gc_df = read_gc_table(args.gc_table)
        # Merge if keys present, else assume same order
        if set(["chrom","start","end"]).issubset(gc_df.columns) and set(["chrom","start","end"]).issubset(df.columns):
            m = pd.merge(df[["chrom","start","end","counts"]], gc_df, on=["chrom","start","end"], how="inner")
            if m.shape[0] > 0:
                y = m["counts"].astype(float).values
                x = m["gc"].astype(float).values
            elif gc_df.shape[0] == df.shape[0]:
                # Chromosome prefixes can differ after BAM/reference harmonization.
                # The workflow preserves window order, so row order is a safe fallback.
                y = counts
                x = gc_df["gc"].astype(float).values
            else:
                raise ValueError("GC table keys do not match counts table and row counts differ.")
        else:
            # same row order
            if gc_df.shape[0] != df.shape[0]:
                raise ValueError("GC table row count does not match counts table and no keys to merge on.")
            y = counts
            x = gc_df["gc"].astype(float).values
        # normalize counts by binsize to mitigate variable windows (if any)
        if "binsize" in df.columns:
            y = y / df["binsize"].astype(float).values
        # drop NaNs
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() > 2:
            gc_r = float(np.corrcoef(x[mask], y[mask])[0,1])

    # Optional: preseq complexity.
    complexity_at_observed = np.nan
    complexity_saturation  = np.nan
    preseq_status = "missing"
    if args.preseq and os.path.exists(args.preseq):
        try:
            p_df = pd.read_csv(args.preseq, sep="\t", comment="#")
            if {"TOTAL_READS", "EXPECTED_DISTINCT"}.issubset(set(p_df.columns)):
                total_col = "TOTAL_READS"
                distinct_col = "EXPECTED_DISTINCT"
                numeric = p_df[[total_col, distinct_col]].apply(pd.to_numeric, errors="coerce").dropna()
            else:
                p_df2 = pd.read_csv(args.preseq, sep="\t", comment="#", header=None)
                if p_df2.shape[1] >= 2:
                    numeric = p_df2.iloc[:, [0, 1]].apply(pd.to_numeric, errors="coerce").dropna()
                    numeric.columns = ["TOTAL_READS", "EXPECTED_DISTINCT"]
                    total_col = "TOTAL_READS"
                    distinct_col = "EXPECTED_DISTINCT"
                else:
                    numeric = pd.DataFrame(columns=["TOTAL_READS", "EXPECTED_DISTINCT"])

            if numeric.empty:
                preseq_status = "empty"
            else:
                idx = (numeric[total_col].astype(float) - total).abs().idxmin()
                row = numeric.loc[idx]
                tot = float(row[total_col])
                distinct = float(row[distinct_col])
                if tot > 0:
                    complexity_at_observed = distinct
                    complexity_saturation = distinct / tot
                    preseq_status = "ok"
                else:
                    preseq_status = "stub"
        except Exception as e:
            preseq_status = "error"
            print(f"[qc_from_counts] Warning reading preseq file {args.preseq}: {e}", file=sys.stderr)


    # Assemble one-row output
    out = {
        "Library": args.sample,
        "n.bins": int(df.shape[0]),
        "avg.binsize": binsize_mean,
        "total.read.count": total,
        "avg.read.count": avg_read_count,
        "spikiness": spike,
        "entropy": ent,
        "coverage_gini": gini,
        "coverage_cv": coeff_var,
        "coverage_mad": abs_mad,
        "coverage_sd": sd_cov,
        "fold80_penalty": fold80,
        "gc_pearson_r": gc_r,
        "background": bin_background,
        "preseq_distinct_at_observed": complexity_at_observed,
        "preseq_saturation": complexity_saturation,
        "preseq_curve_status": preseq_status,
    }
    out.update(pct_metrics)

    # Write TSV (single line)
    odir = os.path.dirname(args.out)
    if odir and not os.path.exists(odir):
        os.makedirs(odir, exist_ok=True)
    pd.DataFrame([out]).to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()
