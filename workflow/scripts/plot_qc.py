#!/usr/bin/env python3
import argparse, os, sys, textwrap
import pandas as pd
import numpy as np

# NEW: force non-interactive backend
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

PLOT_ALIASES = {
    "#Mapped": "alf_mapped_n",
    "MappedFraction": "alf_mapped_frac",
    "MappedProperFraction": "alf_mapped_proper_pair_frac",
    "DuplicateFraction": "alf_duplicate_frac",
    "UnmappedFraction": "alf_unmapped_frac",
    "ErrorRate": "alf_error_frac",
    "MedianInsertSize": "alf_insert_size_med",
    "MedianCoverage": "alf_coverage_med",
    "SDCoverage": "alf_coverage_sd",
    "FractionCovered": "alf_covered_frac",
    "entropy": "bin_entropy",
    "spikiness": "bin_spikiness",
    "gc_r": "bin_gc_r",
    "total.read.count": "bin_total_read_count",
}

def add_plot_aliases(df):
    """Expose stable plotting aliases while keeping final_qc.tsv prefixed."""
    for alias, source in PLOT_ALIASES.items():
        if alias not in df.columns and source in df.columns:
            df[alias] = df[source]
    return df

def nice_num(x):
    if pd.isna(x): return "NA"
    try:
        x = float(x)
    except Exception:
        return str(x)
    if abs(x) >= 1e6:
        return f"{x/1e6:.2f}M"
    if abs(x) >= 1e3:
        return f"{x/1e3:.2f}k"
    return f"{x:.2f}"

def add_title(fig, title):
    fig.suptitle(title, fontsize=16, fontweight="bold", y=0.99)

def page_summary(fig, row):
    ax = fig.add_axes([0.05, 0.12, 0.9, 0.8])
    ax.axis("off")

    left = [
        ("Sample", row.get("Sample", "NA")),
        ("Library", row.get("Library", "NA")),
        ("Reads mapped", nice_num(row.get("#Mapped", np.nan))),
        ("Mapped fraction", f"{100*float(row.get('MappedFraction', np.nan)):.2f}%"
                           if pd.notna(row.get("MappedFraction")) else "NA"),
        ("Proper pairs", f"{100*float(row.get('MappedProperFraction', np.nan)):.2f}%"
                         if pd.notna(row.get("MappedProperFraction")) else "NA"),
        ("Duplicates", f"{100*float(row.get('DuplicateFraction', np.nan)):.2f}%"
                       if pd.notna(row.get("DuplicateFraction")) else "NA"),
        ("Error rate", f"{100*float(row.get('ErrorRate', np.nan)):.3f}%"
                       if pd.notna(row.get("ErrorRate")) else "NA"),
    ]
    right = [
        ("Median insert size", nice_num(row.get("MedianInsertSize", np.nan))),
        ("Median coverage", nice_num(row.get("MedianCoverage", np.nan))),
        ("SD coverage", nice_num(row.get("SDCoverage", np.nan))),
        ("Fraction covered", f"{100*float(row.get('FractionCovered', np.nan)):.2f}%"
                             if pd.notna(row.get("FractionCovered")) else "NA"),
        ("Entropy", f"{float(row.get('entropy', np.nan)):.3f}"
                    if pd.notna(row.get("entropy")) else "NA"),
        ("Spikiness", f"{float(row.get('spikiness', np.nan)):.4f}"
                      if pd.notna(row.get("spikiness")) else "NA"),
        ("GC correlation", f"{float(row.get('gc_r', np.nan)):.4f}"
                           if pd.notna(row.get("gc_r")) else "NA"),
        ("Preseq saturation", f"{100*float(row.get('preseq_saturation', np.nan)):.2f}%"
                              if pd.notna(row.get("preseq_saturation")) else "NA"),
    ]

    text_lines = []
    w = 38
    text_lines.append("SUMMARY\n" + "—"*w)
    for k,v in left:
        text_lines.append(f"{k:<22} : {v}")
    text_lines.append("")
    text_lines.append("COVERAGE / COMPLEXITY\n" + "—"*w)
    for k,v in right:
        text_lines.append(f"{k:<22} : {v}")

    ax.text(0.02, 0.95, "\n".join(text_lines), va="top", ha="left",
            family="DejaVu Sans Mono", fontsize=11)

def page_bars(fig, row):
    ax = fig.add_axes([0.1, 0.18, 0.85, 0.75])
    labels = [
        "Mapped", "ProperPair", "Duplicate", "Unmapped",
        "ErrorRate", "FractionCovered"
    ]
    keys = [
        "MappedFraction", "MappedProperFraction", "DuplicateFraction",
        "UnmappedFraction", "ErrorRate", "FractionCovered"
    ]
    vals = []
    for k in keys:
        v = row.get(k, np.nan)
        vals.append(np.nan if pd.isna(v) else float(v))

    # Convert rates to %
    percent_idx = [0,1,2,3,5]
    y = []
    for i,v in enumerate(vals):
        if pd.isna(v): y.append(np.nan)
        else: y.append(100*v if i in percent_idx else 100*v)  # keep all in %
    bars = ax.bar(labels, y)
    ax.set_ylabel("Percent")
    ax.set_ylim(0, 100)
    ax.grid(axis='y', alpha=0.2)
    for rect, val in zip(bars, y):
        if pd.notna(val):
            ax.text(rect.get_x()+rect.get_width()/2, rect.get_height()+2,
                    f"{val:.1f}%", ha="center", va="bottom", fontsize=9)
    ax.set_title("Key Fractions", fontweight="bold")

def page_counts_qc(fig, row):
    ax = fig.add_axes([0.1, 0.18, 0.85, 0.75])
    labels = ["Entropy", "Spikiness", "GC r"]
    vals = [
        row.get("entropy", np.nan),
        row.get("spikiness", np.nan),
        row.get("gc_r", np.nan),
    ]
    # handle nan
    heights = [np.nan if pd.isna(v) else float(v) for v in vals]
    bars = ax.bar(labels, [0 if pd.isna(v) else v for v in heights])
    ax.set_ylabel("Value")
    ax.grid(axis='y', alpha=0.2)
    for rect, val in zip(bars, heights):
        if pd.notna(val):
            ax.text(rect.get_x()+rect.get_width()/2, rect.get_height()+max(0.02, 0.02*rect.get_height()),
                    f"{val:.3g}", ha="center", va="bottom", fontsize=9)
    ax.set_title("Counts-based QC", fontweight="bold")

def page_preseq(fig, sample, preseq_dir):
    ax = fig.add_axes([0.12, 0.18, 0.84, 0.75])
    pfile = os.path.join(preseq_dir, f"{sample}.lc.tsv")
    if not os.path.exists(pfile):
        ax.text(0.5, 0.5, "No preseq curve found.", ha="center", va="center")
        ax.axis("off")
        return
    df = pd.read_csv(pfile, sep="\t", comment="#", header=None, names=["total", "distinct"])
    if df.empty:
        ax.text(0.5, 0.5, "Empty preseq file.", ha="center", va="center")
        ax.axis("off")
        return
    ax.plot(df["total"], df["distinct"], lw=2)
    ax.set_xlabel("Total reads (mapped primary)")
    ax.set_ylabel("Distinct reads (pred.)")
    ax.set_title("Library Complexity (preseq lc_extrap)", fontweight="bold")
    ax.grid(alpha=0.2)

# --- ADD: small helpers for cohort plots ---
def _col_ok(df, col):
    return (col in df.columns) and (df[col].notna().any())

def _jitter(n, width=0.10, rng=None):
    if rng is None:
        rng = np.random.default_rng(42)
    return rng.uniform(-width, width, size=n)

def add_panel_violin_box_jitter(ax, data, title, ylabel=None, log=False):
    """Draw violin + box + jitter for a 1D array-like."""
    vals = np.array(pd.to_numeric(pd.Series(data), errors="coerce"), dtype=float)
    vals = vals[~np.isnan(vals)]
    if vals.size == 0:
        ax.text(0.5, 0.5, f"No data for {title}", ha="center", va="center")
        ax.set_axis_off()
        return

    # violin
    vp = ax.violinplot(vals, positions=[1], showmeans=False, showmedians=False, widths=0.8)
    # box
    ax.boxplot(vals, positions=[1], widths=0.25, manage_ticks=False)
    # jitter
    x = 1 + _jitter(len(vals), 0.10)
    ax.scatter(x, vals, s=9, alpha=0.4)

    if log:
        ax.set_yscale("log")
    ax.set_xlim(0.4, 1.6)
    ax.set_xticks([1])
    ax.set_xticklabels([title])
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.grid(axis="y", alpha=0.2)


def make_run_summary_pdf(df, out_pdf, preseq_dir=None):
    """
    Build a single multi-panel PDF summarizing all libraries in final_qc.tsv.
    Plots what’s available; silently skips missing columns.
    """
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    # Collect candidate metrics -> (column_name, pretty_title, ylabel, logscale?)
    # We only plot ones that exist in df.
    candidates = [
        ("#Mapped",               "Mapped reads",                "reads",        True),
        ("MappedFraction",        "Mapped fraction",             "fraction",     False),
        ("MappedProperFraction",  "Proper pairs",                "fraction",     False),
        ("DuplicateFraction",     "Duplicates",                  "fraction",     False),
        ("UnmappedFraction",      "Unmapped",                    "fraction",     False),
        ("ErrorRate",             "Error rate",                  "fraction",     False),
        ("MedianInsertSize",      "Median insert size",          "bp",           False),
        ("MedianCoverage",        "Median coverage",             "x",            True),
        ("SDCoverage",            "Coverage SD",                 "x",            True),
        ("FractionCovered",       "Fraction covered",            "fraction",     False),
        ("entropy",               "Entropy",                     "",             False),
        ("spikiness",             "Spikiness",                   "",             False),
        ("gc_r",                  "GC correlation",              "",             False),
        ("preseq_saturation",     "Preseq saturation",           "fraction",     False),
    ]
    metrics = [(c,t,y,l) for (c,t,y,l) in candidates if _col_ok(df, c)]

    # If nothing, make a stub page
    with PdfPages(out_pdf) as pdf:
        # Cover / header page
        fig = plt.figure(figsize=(11, 8.5))
        fig.suptitle("Run-level QC Summary", fontsize=18, fontweight="bold", y=0.98)
        ax = fig.add_axes([0.05, 0.12, 0.9, 0.8]); ax.axis("off")
        n_lib = df.shape[0]
        lib_names = df["Library"].tolist() if "Library" in df.columns else [f"S{i+1}" for i in range(n_lib)]
        lines = [
            f"Total libraries: {n_lib}",
            "",
            "This PDF summarizes per-library QC metrics across the run.",
            "Panels show violin + box + jitter to visualize distribution and outliers.",
            "",
            "Included metrics:",
        ]
        for _, title, _, _ in metrics:
            lines.append(f"  • {title}")
        ax.text(0.02, 0.95, "\n".join(lines), va="top", ha="left", fontsize=12)
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        if not metrics:
            # no data to plot
            return

        # Make pages with up to 6 panels each
        per_page = 6
        pages = [metrics[i:i+per_page] for i in range(0, len(metrics), per_page)]

        for page_metrics in pages:
            fig = plt.figure(figsize=(11, 8.5))
            fig.suptitle("Run-level QC – Distributions", fontsize=16, fontweight="bold", y=0.98)
            gs = GridSpec(3, 2, figure=fig, left=0.07, right=0.98, top=0.92, bottom=0.08, hspace=0.35, wspace=0.25)

            for idx, (col, title, ylabel, logscale) in enumerate(page_metrics):
                r, c = divmod(idx, 2)
                ax = fig.add_subplot(gs[r, c])
                add_panel_violin_box_jitter(ax, df[col], title, ylabel=ylabel, log=logscale)

            # Optional: a tiny outlier table on the last axis if space allows and Library present
            # (We skip to keep things clean. Easy to add later.)
            pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # OPTIONAL: overview scatter (proper vs duplicate, etc.) if both cols exist
        combos = [
            ("MappedProperFraction", "DuplicateFraction", "Proper vs Duplicate"),
            ("MedianCoverage", "DuplicateFraction", "Coverage vs Duplicate"),
            ("entropy", "spikiness", "Entropy vs Spikiness"),
        ]
        any_combo = False
        for xcol, ycol, title in combos:
            if _col_ok(df, xcol) and _col_ok(df, ycol):
                any_combo = True
                fig = plt.figure(figsize=(11, 8.5))
                ax = fig.add_axes([0.1, 0.15, 0.85, 0.78])
                ax.scatter(df[xcol], df[ycol], s=18, alpha=0.6)
                ax.set_xlabel(xcol); ax.set_ylabel(ycol)
                ax.set_title(title, fontweight="bold")
                ax.grid(alpha=0.2)
                pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # OPTIONAL: preseq snapshot – fraction saturated histogram if present
        if _col_ok(df, "preseq_saturation"):
            fig = plt.figure(figsize=(11, 8.5))
            ax = fig.add_axes([0.1, 0.15, 0.85, 0.78])
            vals = pd.to_numeric(df["preseq_saturation"], errors="coerce").dropna()
            if len(vals) > 0:
                ax.hist(vals, bins=20)
                ax.set_title("Preseq saturation – distribution", fontweight="bold")
                ax.set_xlabel("Saturation (fraction)"); ax.set_ylabel("Libraries")
                ax.grid(axis="y", alpha=0.2)
            else:
                ax.text(0.5, 0.5, "No preseq_saturation values", ha="center", va="center")
                ax.axis("off")
            pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


def write_summary_pdf(df_full, outdir, summary_name="summary.pdf", preseq_dir=None):
    os.makedirs(outdir, exist_ok=True)
    out_pdf = os.path.join(outdir, summary_name)
    make_run_summary_pdf(df_full, out_pdf, preseq_dir=preseq_dir)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", required=True, help="final_qc.tsv")
    ap.add_argument("--preseq_dir", required=True, help="dir with {sample}.lc.tsv")
    ap.add_argument("--outdir", required=True, help="output directory for PDFs")
    ap.add_argument("--only", default=None, help="Plot only this Library/sample (exact match to 'Library' column).")
    # NEW:
    ap.add_argument("--make-summary", action="store_true",
                    help="If set, generate run-level summary.pdf.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df_full = add_plot_aliases(pd.read_csv(args.final, sep="\t"))

    # NEW: only write summary when explicitly requested
    if args.make_summary:
        write_summary_pdf(df_full, args.outdir, summary_name="summary.pdf", preseq_dir=args.preseq_dir)

    # If a specific library is requested, generate that per-sample PDF
    df = df_full
    if args.only is not None:
        df = df.loc[df["Library"] == args.only].copy()
        if df.empty:
            print(f"[plot_qc] WARNING: --only {args.only} not found in final_qc.tsv", file=sys.stderr)
            return 0

    if "Library" not in df.columns:
        print("[plot_qc] ERROR: 'Library' column not found in final_qc.", file=sys.stderr)
        sys.exit(1)

    for _, row in df.iterrows():
        lib = row["Library"]
        pdf_path = os.path.join(args.outdir, f"{lib}.qc.pdf")
        with PdfPages(pdf_path) as pdf:
            fig = plt.figure(figsize=(11.0, 8.5))
            add_title(fig, f"QC Summary – {lib}")
            page_summary(fig, row)
            pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

            fig = plt.figure(figsize=(11.0, 8.5))
            add_title(fig, f"Key Fractions – {lib}")
            page_bars(fig, row)
            pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

            fig = plt.figure(figsize=(11.0, 8.5))
            add_title(fig, f"Counts-based QC – {lib}")
            page_counts_qc(fig, row)
            pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

            fig = plt.figure(figsize=(11.0, 8.5))
            add_title(fig, f"Preseq – {lib}")
            page_preseq(fig, lib, args.preseq_dir)
            pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    return 0

if __name__ == "__main__":
    sys.exit(main())
