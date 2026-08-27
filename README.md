# bamqc-pipeline

`bamqc-pipeline` is a Snakemake workflow for per-library BAM QC in Strand-seq and mosaicatcher-style projects. It discovers BAMs, computes alignment and coverage metrics, estimates library complexity, optionally folds in scTRIP/mosaicatcher count tables and Ashley QC, and writes one consistent `final_qc.tsv` for downstream filtering and plotting.

## What It Produces

The main deliverable is:

```text
final_qc.tsv
```

Each row is one discovered library/cell. Columns are grouped by prefix so the source of each metric stays obvious:

- `alf_*`: alignment, mapping, duplicate, insert-size, MAPQ, and error metrics from Alfred.
- `bin_*`: fixed-window coverage shape metrics from bedtools-derived bins, including `bin_background` from Watson/Crick-biased bins.
- `preseq_*`: library-complexity estimates from preseq.
- `sctrip_*`: optional HIER-mode scTRIP/mosaicatcher `counts/*.info_raw` or fallback `counts/*.info` metrics.
- `ash_*`: optional Ashley labels/features, copied from existing outputs or generated when Ashley is enabled.

## Input Layouts

The workflow auto-detects one of two layouts from `data_location`:

```text
FLAT
/path/to/run/*.sort.mdup.bam

HIER
/path/to/run/<Sample>/bam/*.sort.mdup.bam
/path/to/run/<Sample>/counts/*.info_raw   # optional, preferred
/path/to/run/<Sample>/counts/*.info       # optional fallback
```

`bam_ext` controls the BAM suffix. HIER mode creates stable `Library` identifiers by combining the sample folder and BAM basename, while retaining `Sample` for grouping.

## Installation

Clone the workflow into a path without spaces. Ashley QC is only needed when `ashleys.enabled: true` and usable Ashley outputs are not already present.

```bash
cd
mkdir -p work
cd work
git clone https://github.com/pweidner/bamqc-pipeline.git
git clone https://github.com/friendsofstrandseq/ashleys-qc.git
cd bamqc-pipeline
```

Create or activate a Snakemake environment. The rule-specific tools are installed from `workflow/envs/*.yaml` during execution.

```bash
mamba create -n snakemake snakemake=7.32.0 -y
conda activate snakemake
```

## Configuration

Tracked defaults live in `config/config.yaml`. For cluster runs, put machine- and run-specific values in ignored `config/config.local.yaml`; it is loaded automatically when present.

```yaml
ref: hg38
reference_path: /ref/dir               # contains hg38.fa and hg38.fa.fai
data_location: /path/to/input          # FLAT or HIER input root
output_location: /path/to/output
window: 200000
chromosomes: "1-22,X,Y"                # or "all" for every reference contig
plot: true

bam_ext: ".sort.mdup.bam"
tmp_dir: /tmp
large_bam_threshold_mb: 800

ashleys:
  enabled: true
  bin: /abs/path/ashleys-qc/bin/ashleys.py
  model_path: /abs/path/ashleys-qc/models/svc_default.pkl
  win_sizes: [5000000,2000000,1000000,800000,600000,400000,200000]
  threads: 32
  mem_mb: 200000
  conda_env: envs/ashleys.yaml
```

BAM-dependent jobs use `large_bam_threshold_mb` to choose between small and large resource profiles. By default, BAMs over 800 MB use the `medium` partition, more memory, and a longer runtime for `alfred_qc`, `coverage_counts`, and `preseq_lc`.

## Run

```bash
cp config/config.local.example.yaml config/config.local.yaml
# edit config/config.local.yaml for the run
snakemake --profile workflow/profiles --keep-going
```

When reusing an existing output directory from an older `bamqc-pipeline` version, force regeneration of binned counts or use a fresh `output_location` to populate `bin_background`; older `binned/*.bins.tsv.gz` files do not contain Watson/Crick counts.

```bash
snakemake --profile workflow/profiles --keep-going --forcerun coverage_counts qc_from_counts make_core_qc finalize_qc plot_run_summary_R
```

## Output Structure

```text
output_location/
├── final_qc.tsv                         # single final table; Ashley columns included when available
├── alignment_summary_metrics.tsv        # parsed Alfred summary across libraries
├── sctrip_counts_info.tsv               # optional HIER counts/*.info_raw or counts/*.info aggregation
│
├── metadata/
│   └── library_map.tsv                  # cell/BAM to Library/Sample mapping
│
├── stats-by-lib/
│   └── {Library}.qc.tsv.gz              # per-library Alfred output
│
├── binned/
│   └── {Library}.bins.tsv.gz            # chrom/start/end/Watson/Crick read1/total bin counts
│
├── qc-from-bins/
│   └── {Library}.counts_qc.tsv          # bin-derived metrics before final merge
│
├── preseq/
│   └── {Library}.lc.tsv                 # library complexity curves
│
├── ashleys/
│   ├── features.tsv                     # merged or computed Ashley features
│   └── prediction/
│       └── prediction.tsv               # merged labels or predictions
│
└── plots/
    ├── per-lib-qc/{Library}.qc.pdf      # optional per-library PDF
    └── run_summary.pdf                  # run-level QC relationship plots
```

## Output Metrics (what they mean and how to read them)

This pipeline produces per-library QC summaries in these primary tables:

- **`final_qc.tsv`** — the single final QC table. It always includes core **Alfred**, **bin-wise coverage**, and **preseq** metrics; adds optional HIER **scTRIP/mosaicatcher counts info**; and adds `ash_*` columns when Ashley outputs are available or Ashley computation is enabled.
- **`sctrip_counts_info.tsv`** — optional run-level aggregation of HIER `counts/*.info_raw` tables, falling back to `counts/*.info` when raw tables are absent.

All non-identifier columns are prefixed by their producing tool to make provenance explicit.

---

## 1. Identifiers (no prefix)

| Column | Description |
|------|-------------|
| `Library` | Unique library identifier used throughout the pipeline (e.g. `DRUG-CDTR-P1AZA-C_A5573_L1_i301.sort.mdup`). |
| `Sample` | Sample / condition identifier grouping multiple libraries (e.g. `DRUG-CDTR-P1AZA-C`). |

---

## 2. Alfred alignment and BAM QC (`alf_*`)

Metrics derived from **`alfred qc`**, summarizing mapping, alignment accuracy, and coverage statistics.

### Read filtering & mapping

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `alf_qcfail_n` | QC-failed reads | High values indicate poor read quality. |
| `alf_qcfail_frac` | Fraction QC-failed | >0.05 often indicates a problematic library. |
| `alf_duplicate_marked_n` | Duplicate reads | |
| `alf_duplicate_frac` | Duplicate fraction | **High = low complexity or over-sequencing.** |
| `alf_unmapped_n` | Unmapped reads | |
| `alf_unmapped_frac` | Fraction unmapped | High values may indicate contamination or wrong reference. |
| `alf_mapped_n` | Mapped reads | |
| `alf_mapped_frac` | Fraction mapped | Healthy libraries are typically high (>0.8). |

### Read balance & orientation

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `alf_mapped_read1_n`, `alf_mapped_read2_n` | Read1 / Read2 mapped counts | |
| `alf_mapped2_vs_mapped1_ratio` | Read2 / Read1 ratio | **~1.0 expected** for paired-end data. |
| `alf_mapped_forward_frac` | Fraction forward strand | |
| `alf_mapped_reverse_frac` | Fraction reverse strand | **~0.5 / 0.5 expected** unless protocol-biased. |

### Alignment types

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `alf_secondary_alignments_frac` | Secondary alignments | Elevated values indicate multi-mapping / repeats. |
| `alf_supplementary_alignments_frac` | Supplementary alignments | Can indicate SVs, chimeras, or mapping artifacts. |
| `alf_spliced_alignments_frac` | Spliced alignments | Typically low for DNA; higher for RNA or odd mapping. |

### Pairing & concordance

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `alf_mapped_pairs_frac` | Fraction of mapped read pairs | |
| `alf_mapped_same_chr_frac` | Pairs on same chromosome | Low values may indicate discordant mapping. |
| `alf_mapped_proper_pair_frac` | Properly paired reads | **High = good library structure.** |

### Alignment accuracy & errors

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `alf_match_frac` | Matched bases / aligned bases | **Closer to 1 = higher accuracy.** |
| `alf_mismatch_frac` | Mismatched base fraction | |
| `alf_deletion_frac`, `alf_insertion_frac` | Indel rates | Elevated values can indicate mapping or chemistry artifacts. |
| `alf_error_frac` | Aggregate alignment error rate | Lower is better. |

### Clipping & context

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `alf_soft_clip_frac` | Soft-clipped bases | High = adapters, short inserts, or mapping difficulty. |
| `alf_hard_clip_frac` | Hard-clipped bases | Usually near zero. |
| `alf_homopolymer_context_del/ins` | Indels in homopolymers | Elevated values indicate systematic indel artifacts. |

### Read length, coverage & MAPQ

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `alf_read_length_med` | Median read length | |
| `alf_insert_size_med` | Median insert size | |
| `alf_mapq_med` | Median mapping quality | **Higher = more confident mapping.** |
| `alf_coverage_med` | Median coverage | |
| `alf_coverage_sd` | Coverage standard deviation | High = uneven coverage. |
| `alf_covered_frac` | Fraction of reference covered | Low = sparse library. |

---

## 3. Bin-wise coverage metrics (`bin_*`)

Computed from fixed-size genome windows using `bedtools coverage -counts` and summarized in `qc_from_counts.py`. Bin files store Watson-strand (`-`) and Crick-strand (`+`) counts from non-read2 reads plus total filtered read counts per window, so the same scheduled job supports coverage-shape QC and Strand-seq background estimation.

### Basic bin descriptors

| Column | Meaning |
|------|--------|
| `bin_n_bins` | Number of windows used. |
| `bin_avg_binsize` | Mean window size (bp). |
| `bin_total_read_count` | Total reads across all bins. |
| `bin_avg_read_count` | Mean reads per bin. |

### Coverage uniformity & signal shape

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `bin_entropy` | Shannon entropy of bin counts | **Higher = more even coverage.** |
| `bin_spikiness` | Local coverage jaggedness | **Higher = noisier / uneven signal.** |
| `bin_gini` | Gini index of coverage | **0 = uniform, higher = uneven.** |
| `bin_cv` | Coefficient of variation | |
| `bin_mad` | Median absolute deviation | Robust variability measure. |
| `bin_sd` | Standard deviation | |
| `bin_background` | breakpointR-inspired strand background estimate from Watson/Crick-biased bins | **Lower = cleaner strand separation.** Useful with high `sctrip_good` or `bin_total_read_count` and low `bin_spikiness` when selecting strong libraries. |

### Uniformity & GC bias

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `bin_fold80` | Fold-80 penalty | **~1 ideal**, higher = worse uniformity. |
| `bin_gc_r` | GC-coverage correlation | Large magnitude = GC bias. |

### Depth thresholds

| Column | Meaning |
|------|--------|
| `bin_pct_ge_1x` | Fraction of bins with >=1 read. |
| `bin_pct_ge_10x` | Fraction with >=10 reads. |
| `bin_pct_ge_30x` | Fraction with >=30 reads. |

---

## 4. Library complexity (preseq) (`preseq_*`)

Derived from **`preseq lc_extrap`**, estimating how many *unique* DNA fragments are present.

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `preseq_distinct_at_observed` | Expected number of distinct fragments at observed depth | Higher = more complex library. |
| `preseq_saturation` | Distinct / total reads at observed depth | **0 = highly duplicated**, **1 = highly complex**. |
| `preseq_curve_status` | Whether the preseq curve was parsed as `ok`, `missing`, `empty`, `stub`, or `error`. | Non-`ok` values should be checked in logs. |

**Rule of thumb**
- `preseq_saturation ~ 1` -> sequencing deeper will still yield new information
- `preseq_saturation ~ 0` -> sequencing deeper mostly yields duplicates

---

## 5. scTRIP / mosaicatcher counts info (`sctrip_*`)

Available in HIER mode when per-sample `counts/*.info_raw` or `counts/*.info` files exist. Raw tables are preferred because they retain all cells.

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `sctrip_sample` | mosaicatcher sample label from the counts info table | Should match `Sample`. |
| `sctrip_cell` | Cell/library name from the counts info table | Used with `Sample`/BAM path to join into `Library`. |
| `sctrip_mapped` | Total reads seen by scTRIP/mosaicatcher counting | Global input depth before filters. |
| `sctrip_suppl` | Supplementary, secondary, or QC-failed reads filtered out | High values suggest alignment artifacts or poor read quality. |
| `sctrip_dupl` | PCR duplicate reads filtered out | High values indicate low library complexity. |
| `sctrip_mapq` | Reads filtered out for low mapping quality | High values suggest ambiguous mapping. |
| `sctrip_read2` | Read 2 records filtered out by the counter | Expected for read-1-oriented counting logic. |
| `sctrip_good` | Reads used for counting | Useful direct comparison to BamQC `bin_total_read_count`. |
| `sctrip_pass1` | mosaicatcher coverage pass flag | `0` means downstream NB fields should be treated cautiously. |
| `sctrip_nb_p`, `sctrip_nb_r`, `sctrip_nb_a` | Negative-binomial parameters estimated by mosaicatcher | Present mainly for compatibility and deeper troubleshooting. |
| `sctrip_bam` | BAM path recorded in the counts info table | Helps audit joins. |
| `sctrip_counts_info_path` | Source `*.info_raw` or `*.info` file | Helps trace run-level aggregation. |

---

## 6. Ashley QC predictions (`ash_*`)

Merged from existing Ashley labels/predictions when present. If `ashleys.enabled: true` and labels/predictions are missing, they are generated with **ashleys-qc**. If `ashleys.enabled: false`, missing Ashley predictions are not generated.

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `ash_label` | Predicted QC class (model-specific) | |
| `ash_prob` | Prediction confidence | Values near 0.5 = ambiguous. |
| `ash_cell` | BAM identifier used by Ashley | |
| `ash_sample` | Sample label used by Ashley | |

---

## 7. Ashley feature vectors (`ash_*`)

Multi-scale **Watson-strand bin features** and read category fractions.

### Window-bin distributions

For each window size (`5mb`, `2mb`, `1mb`, `0_8mb`, `0_6mb`, `0_4mb`, `0_2mb`):

| Column pattern | Meaning |
|--------------|--------|
| `ash_w10_*` … `ash_w100_*` | Fraction of windows falling into Watson% bins (0-10%, …, 90-100%). |
| `ash_total_*` | Total fraction across bins (~1.0 if normalized). |

**Interpretation**
- Smooth, balanced distributions indicate stable coverage.
- Skewed distributions can indicate CNVs, strand imbalance, or technical artifacts.

### Mapping & "good read" fractions

| Column | Meaning |
|------|--------|
| `ash_p_unmap` | Fraction unmapped. |
| `ash_p_map` | Fraction mapped. |
| `ash_p_supp` | Fraction supplementary. |
| `ash_p_dup` | Fraction duplicate. |
| `ash_p_mq` | Fraction passing MAPQ filter. |
| `ash_p_read2` | Fraction read2. |
| `ash_p_good` | Fraction of reads passing all Ashley filters (usable signal). |

---

## Practical Interpretation Summary

- **Low complexity** -> high `alf_duplicate_frac`, low `preseq_saturation`
- **Uneven coverage** -> high `bin_spikiness`, `bin_gini`, `bin_fold80`
- **High strand background** -> high `bin_background`; strongest libraries combine low background, low spikiness, and high `sctrip_good` / `bin_total_read_count`
- **GC bias** -> large `|bin_gc_r|`
- **Mapping problems** -> low `alf_mapped_frac`, low `alf_mapq_med`
- **scTRIP count failure** -> `sctrip_pass1 == 0` or very low `sctrip_good` -> inspect mosaicatcher counts logs
- **Ashley disagreement** -> `ash_prob low` -> inspect manually

---



## Notes
- **Entropy** and **spikiness** reflect coverage evenness (low entropy or high spikiness = uneven).
- **Fold80 penalty** follows the Picard metric (ideal = 1, higher = less uniform).
- **Preseq** metrics allow extrapolation of unique reads vs sequencing depth.
- **scTRIP/mosaicatcher counts info** is included only for HIER runs that contain per-sample `counts/*.info_raw` or fallback `counts/*.info` tables.
- **Ashley's QC** integrates pretrained classification of Strand-seq libraries by coverage pattern and W->C balance.

---

## 8. Citations

- **Alfred** — Rausch *et al.*, *Genome Res* (2019)
- **preseq** — Daley & Smith, *Bioinformatics* (2013)
- **ASHLEYS** - Gros *et al.*, *Bioinformatics* (2021)
- **bedtools** — Quinlan & Hall, *Bioinformatics* (2010)
- **Snakemake** — Köster & Rahmann, *Bioinformatics* (2012)

---

## 9. Roadmap

- contamination checks
- More GC/coverage plots (Lorenz, violin)
- Optional HTML dashboard

---

Happy QC-ing.
