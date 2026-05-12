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
- `metadata/library_map.tsv`: mapping from input BAM basename to the stable
  `Library` identifier used throughout the workflow.

## Column Naming

- `Library` and `Sample` are the join keys.
- `alf_*` columns come from Alfred.
- `bin_*` columns come from genome-window read counts.
- `preseq_*` columns come from preseq library-complexity curves.
- `sctrip_*` columns come from optional HIER scTRIP/mosaicatcher counts info
  tables.
- `ash_*` columns come from existing or newly computed Ashley QC outputs.

Plotting scripts consume the prefixed schema above and include legacy aliases
only for compatibility with older output tables.
