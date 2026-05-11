# Output Contract

`bamqc-pipeline` writes one row per discovered library.

## Primary Tables

- `final_qc.tsv`: core Alfred, binned-coverage, and preseq metrics.
- `final_qc_with_ashleys.tsv`: `final_qc.tsv` plus Ashley QC predictions and
  features when `ashleys.enabled: true`.
- `alignment_summary_metrics.tsv`: parsed Alfred `ME` metrics before final
  column normalization.
- `metadata/library_map.tsv`: mapping from input BAM basename to the stable
  `Library` identifier used throughout the workflow.

## Column Naming

- `Library` and `Sample` are the join keys.
- `alf_*` columns come from Alfred.
- `bin_*` columns come from genome-window read counts.
- `preseq_*` columns come from preseq library-complexity curves.
- `ash_*` columns come from Ashley QC.

Plotting scripts consume the prefixed schema above and include legacy aliases
only for compatibility with older output tables.
