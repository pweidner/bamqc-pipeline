#!/usr/bin/env python3
import argparse
import os
import re
import sys
import pandas as pd

NUMERIC_COLUMNS = {
    "mapped", "suppl", "dupl", "mapq", "read2", "good", "pass1", "medbin",
    "nb_p", "nb_r", "nb_a",
}


def snake(name):
    name = str(name).strip().replace(".", "_")
    name = re.sub(r"[^A-Za-z0-9_]+", "_", name)
    name = re.sub(r"__+", "_", name)
    return name.strip("_").lower()


def strip_bam(name):
    name = str(name)
    return name[:-4] if name.endswith(".bam") else name


def strip_configured_ext(name, bam_ext):
    name = str(name)
    suffixes = [bam_ext]
    if str(bam_ext).endswith(".bam"):
        suffixes.append(str(bam_ext)[:-4])
    for suffix in sorted(set(s for s in suffixes if s), key=len, reverse=True):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return name


def value_keys(value, bam_ext):
    if pd.isna(value):
        return []
    raw = str(value).strip()
    if not raw:
        return []
    base = os.path.basename(raw)
    keys = [
        raw,
        os.path.abspath(raw) if os.path.isabs(raw) else raw,
        base,
        strip_bam(base),
        strip_configured_ext(base, bam_ext),
    ]
    if bam_ext and not base.endswith(bam_ext):
        with_ext = base + bam_ext
        keys.extend([with_ext, strip_bam(with_ext), strip_configured_ext(with_ext, bam_ext)])
    return list(dict.fromkeys(k for k in keys if k))


def add_key(mapping, ambiguous, key, library):
    if not key or pd.isna(library):
        return
    library = str(library)
    if key in mapping and mapping[key] != library:
        ambiguous.add(key)
        return
    mapping[key] = library


def library_from_sample_cell(sample, cell, bam_ext):
    if pd.isna(sample) or pd.isna(cell):
        return pd.NA
    sample = str(sample).strip()
    cell = str(cell).strip()
    if not sample or not cell or sample == "nan" or cell == "nan":
        return pd.NA
    ext_no_bam = strip_bam(str(bam_ext or ""))
    base = strip_configured_ext(os.path.basename(cell), bam_ext)
    if ext_no_bam and not base.endswith(ext_no_bam):
        base = base + ext_no_bam
    return sample + "_" + base


def build_library_lookup(libmap, bam_ext):
    key_to_library = {}
    ambiguous = set()
    for _, row in libmap.iterrows():
        library = row.get("Library")
        for col in ("bam", "cell"):
            if col in libmap.columns:
                for key in value_keys(row[col], bam_ext):
                    add_key(key_to_library, ambiguous, key, library)
    for key in ambiguous:
        key_to_library.pop(key, None)
    return key_to_library, ambiguous


def read_counts_info(path):
    df = pd.read_csv(path, sep=r"\s+", comment="#", engine="python")
    df.columns = [snake(c) for c in df.columns]
    df["counts_info_path"] = str(path)
    return df


def resolve_library(row, key_to_library, bam_ext):
    for col in ("bam", "cell"):
        if col in row:
            for key in value_keys(row[col], bam_ext):
                if key in key_to_library:
                    return key_to_library[key]
    return library_from_sample_cell(row.get("sample"), row.get("cell"), bam_ext)


def main():
    ap = argparse.ArgumentParser(
        description="Aggregate mosaicatcher/scTRIP *.info_raw or *.info tables and normalize them to BamQC Library IDs."
    )
    ap.add_argument("counts_info", nargs="*", help="Sample-level counts/*.info_raw or counts/*.info files")
    ap.add_argument("--libmap", required=True, help="BamQC metadata/library_map.tsv")
    ap.add_argument("--bam-ext", default=".sort.mdup.bam", help="Configured BAM suffix")
    ap.add_argument("--out", required=True, help="Run-level output TSV")
    args = ap.parse_args()

    frames = []
    for path in args.counts_info:
        try:
            df = read_counts_info(path)
        except Exception as exc:
            print(f"[aggregate_sctrip_counts_info] skipped {path}: {exc}", file=sys.stderr)
            continue
        if df.empty:
            print(f"[aggregate_sctrip_counts_info] skipped empty table {path}", file=sys.stderr)
            continue
        frames.append(df)

    raw = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["Library"])
    libmap = pd.read_csv(args.libmap, sep="\t")
    key_to_library, ambiguous = build_library_lookup(libmap, args.bam_ext)
    if ambiguous:
        print(f"[aggregate_sctrip_counts_info] ignored {len(ambiguous)} ambiguous cell/BAM keys", file=sys.stderr)

    if not raw.empty:
        raw.insert(0, "Library", raw.apply(lambda row: resolve_library(row, key_to_library, args.bam_ext), axis=1))
        for col in raw.columns:
            if col in NUMERIC_COLUMNS:
                raw[col] = pd.to_numeric(raw[col], errors="coerce")

        raw = raw.drop_duplicates()
        matched = raw["Library"].notna() & (raw["Library"].astype(str).str.len() > 0)
        duplicated = raw.loc[matched, "Library"].duplicated(keep=False)
        if duplicated.any():
            ndup = int(duplicated.sum())
            print(
                f"[aggregate_sctrip_counts_info] warning: {ndup} rows share duplicated Library IDs; keeping first per Library",
                file=sys.stderr,
            )
            matched_rows = raw.loc[matched].drop_duplicates(subset=["Library"], keep="first")
            raw = pd.concat([matched_rows, raw.loc[~matched]], ignore_index=True)

        unmatched_n = int((raw["Library"].isna() | (raw["Library"].astype(str).str.len() == 0)).sum())
        if unmatched_n:
            print(f"[aggregate_sctrip_counts_info] warning: {unmatched_n} rows could not be mapped to Library", file=sys.stderr)

    rename = {}
    for col in raw.columns:
        if col == "Library":
            continue
        rename[col] = col if col.startswith("sctrip_") else "sctrip_" + snake(col)
    raw = raw.rename(columns=rename)

    preferred = ["Library", "sctrip_sample", "sctrip_cell"]
    ordered = [c for c in preferred if c in raw.columns] + [c for c in raw.columns if c not in preferred]
    raw = raw[ordered] if ordered else pd.DataFrame(columns=["Library"])

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    raw.to_csv(args.out, sep="\t", index=False)
    print(f"[aggregate_sctrip_counts_info] wrote {raw.shape[0]} rows and {raw.shape[1]} columns to {args.out}")


if __name__ == "__main__":
    main()
