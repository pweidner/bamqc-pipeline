# Output Contract

`bamqc-pipeline` writes one row per discovered library.

## Primary Tables

- `final_qc.tsv`: the single final QC table. It always includes core Alfred,
  binned-coverage, and preseq metrics; includes optional HIER scTRIP/mosaicatcher
  counts info metrics; and includes `ash_*` columns when Ashley outputs are
  available or Ashley computation is enabled.
- `sctrip_counts_info.tsv`: optional run-level aggregation of HIER
  `counts/*.info_raw` tables, falling back to `counts/*.info` when raw tables are
  absent for a sample.
- `alignment_summary_metrics.tsv`: parsed Alfred `ME` metrics before final
  column normalization.
- `binned/{Library}.bins.tsv.gz`: per-window `chrom`, `start`, `end`,
  Watson/minus count from non-read2 reads, Crick/plus count from non-read2
  reads, and total filtered read count. Older four-column count files remain
  readable by the QC script but cannot produce `bin_background`.
- `metadata/library_map.tsv`: mapping from input BAM basename to the stable
  `Library` identifier used throughout the workflow.

## Column Naming

- `Library` and `Sample` are the join keys.
- `alf_*` columns come from Alfred.
- `bin_*` columns come from genome-window read counts; `bin_background` is a
  breakpointR-inspired single-value strand background estimate from
  Watson/Crick-biased non-read2 bin counts.
- `preseq_*` columns come from preseq library-complexity curves.
- `sctrip_*` columns come from optional HIER scTRIP/mosaicatcher counts info
  tables.
- `ash_*` columns come from existing or newly computed Ashley QC outputs.

Plotting scripts consume the prefixed schema above and include legacy aliases
only for compatibility with older output tables.
