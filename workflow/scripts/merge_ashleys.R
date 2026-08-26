#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

strip_bam <- function(x) {
  x <- basename(x)
  x <- str_replace(x, "\\.sort\\.mdup\\.bam$", "")
  x <- str_replace(x, "\\.sort\\.mdup$", "")
  x <- str_replace(x, "\\.mdup\\.bam$", "")
  x <- str_replace(x, "\\.mdup$", "")
  x <- str_replace(x, "\\.sort\\.bam$", "")
  x <- str_replace(x, "\\.sort$", "")
  x <- str_replace(x, "\\.bam$", "")
  x
}

read_tsv_quiet <- function(path, ...) {
  readr::read_tsv(path, guess_max = 1e6, show_col_types = FALSE, progress = FALSE, ...)
}

parse_cli_args <- function(argv) {
  if (length(argv) %% 2 != 0) {
    stop("Expected command-line arguments as --key value pairs.")
  }
  args <- list()
  if (length(argv) > 0) {
    for (i in seq(1, length(argv), by = 2)) {
      key <- sub("^--", "", argv[[i]])
      args[[key]] <- argv[[i + 1]]
    }
  }
  missing <- setdiff(c("final", "pred", "out"), names(args))
  if (length(missing) > 0) {
    stop("Missing required arguments: ", paste(paste0("--", missing), collapse = ", "))
  }
  if (is.null(args$feat)) {
    args$feat <- NULL
  }
  args
}

# Find a column in df by case-insensitive match
find_col_ci <- function(df, candidates) {
  nm <- names(df)
  low <- tolower(nm)
  for (cand in candidates) {
    hit <- which(low == tolower(cand))
    if (length(hit) > 0) return(nm[hit[1]])
  }
  return(NULL)
}

# Remove a literal prefix if present (vectorized)
remove_literal_prefix <- function(x, prefix) {
  ok <- !is.na(prefix) & startsWith(x, prefix)
  x2 <- x
  x2[ok] <- substr(x[ok], nchar(prefix[ok]) + 1, nchar(x[ok]))
  x2
}

# Normalize Ashley prediction keys. Historical tables may lack Sample or rows.
normalize_pred_keys <- function(pred) {
  cell_col <- find_col_ci(pred, c("cell", "sample_name", "library", "Library"))
  sample_col <- find_col_ci(pred, c("sample", "Sample", "sample_id"))
  pred_col <- find_col_ci(pred, c("prediction", "label", "ash_label"))
  prob_col <- find_col_ci(pred, c("probability", "prob", "ash_prob"))

  if (is.null(cell_col)) {
    stop("prediction.tsv missing required cell column; columns are: ", paste(names(pred), collapse=", "))
  }

  out <- tibble(
    cell = as.character(pred[[cell_col]]),
    sample = if (!is.null(sample_col)) as.character(pred[[sample_col]]) else NA_character_,
    prediction = if (!is.null(pred_col)) as.character(pred[[pred_col]]) else NA_character_,
    probability = if (!is.null(prob_col)) suppressWarnings(as.numeric(pred[[prob_col]])) else NA_real_
  ) %>%
    filter(!is.na(.data$cell), .data$cell != "") %>%
    mutate(
      Base = strip_bam(.data$cell),
      Library_key = if_else(
        !is.na(.data$sample) & .data$sample != "",
        paste0(.data$sample, "_", .data$Base, ".sort.mdup"),
        NA_character_
      )
    )

  if (is.null(pred_col)) {
    warning("[merge_ashleys] prediction.tsv lacks a prediction/label column; ash_label will be NA.")
  }
  if (is.null(prob_col)) {
    warning("[merge_ashleys] prediction.tsv lacks a probability/prob column; ash_prob will be NA.")
  }
  if (is.null(sample_col) && nrow(out) > 0) {
    warning("[merge_ashleys] prediction.tsv lacks sample column; joining predictions by unique cell basename.")
  }

  out
}

# Normalize Ashley features keys → Base
normalize_feat_keys <- function(feat) {
  nm <- names(feat)
  key <- dplyr::case_when(
    "sample_name" %in% nm ~ "sample_name",
    "cell"        %in% nm ~ "cell",
    TRUE ~ NA_character_
  )
  if (is.na(key)) {
    warning("[merge_ashleys] features.tsv lacks 'sample_name' or 'cell', will skip merging features.")
    return(NULL)
  }
  feat %>% mutate(Base = strip_bam(.data[[key]]))
}

# Normalize Ashley feature column names into safe snake-ish keys (but keep meaning)
normalize_ash_feature_name <- function(nm) {
  nm2 <- tolower(nm)
  nm2 <- str_replace_all(nm2, "\\.", "_")
  # collapse x.0mb -> xmb for 5.0 / 2.0 / 1.0 specifically
  nm2 <- str_replace_all(nm2, "5_0mb", "5mb")
  nm2 <- str_replace_all(nm2, "2_0mb", "2mb")
  nm2 <- str_replace_all(nm2, "1_0mb", "1mb")
  # for decimals like 0.8mb -> 0_8mb etc already handled by '.'->'_'
  nm2 <- str_replace_all(nm2, "__+", "_")
  nm2 <- str_replace_all(nm2, "^_+|_+$", "")
  nm2
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

final_path <- normalizePath(args$final, mustWork = TRUE)
pred_path  <- normalizePath(args$pred,  mustWork = TRUE)
feat_path  <- if (!is.null(args$feat)) normalizePath(args$feat, mustWork = FALSE) else NULL

message("[merge_ashleys] Reading inputs ...")
final <- read_tsv_quiet(final_path)
pred  <- read_tsv_quiet(pred_path)

# ---- Prepare keys on final ----
lib_col  <- find_col_ci(final, c("Library"))
samp_col <- find_col_ci(final, c("Sample"))
if (is.null(lib_col))  stop("[merge_ashleys] final_qc.tsv missing 'Library' column")
if (is.null(samp_col)) stop("[merge_ashleys] final_qc.tsv missing 'Sample' column")

final <- final %>%
  mutate(
    Library_raw = .data[[lib_col]],
    Sample_raw  = .data[[samp_col]]
  )

# Base = remove "<Sample>_" prefix from Library, then strip suffixes
prefix <- paste0(final$Sample_raw, "_")
lib_no_prefix <- remove_literal_prefix(final$Library_raw, prefix)
final <- final %>% mutate(Base = strip_bam(lib_no_prefix))

# Build Library_key matching Ashley normalization
final <- final %>% mutate(Library_key = paste0(.data$Sample_raw, "_", .data$Base, ".sort.mdup"))

# ---- Join predictions ----
message("[merge_ashleys] Joining Ashley predictions ...")
pred2 <- normalize_pred_keys(pred)

if (nrow(pred2) == 0) {
  merged <- final %>%
    mutate(
      ash_cell = NA_character_,
      ash_sample = NA_character_,
      ash_label = NA_character_,
      ash_prob = NA_real_
    )
  message("[merge_ashleys] No Ashley prediction rows found; keeping ash_* prediction columns as NA.")
} else if (any(!is.na(pred2$Library_key))) {
  pred_min <- pred2 %>%
    filter(!is.na(.data$Library_key)) %>%
    distinct(.data$Library_key, .keep_all = TRUE) %>%
    transmute(
      Library_key = .data$Library_key,
      ash_cell    = .data$cell,
      ash_sample  = .data$sample,
      ash_label   = .data$prediction,
      ash_prob    = .data$probability
    )

  merged <- final %>% left_join(pred_min, by = "Library_key")
  message("[merge_ashleys] Matched predictions: ", sum(!is.na(merged$ash_cell)), " / ", nrow(merged))
} else {
  if (any(duplicated(final$Base))) {
    stop("[merge_ashleys] prediction.tsv has no sample column and final table has duplicated cell basenames; cannot join unambiguously.")
  }
  if (any(duplicated(pred2$Base))) {
    stop("[merge_ashleys] prediction.tsv has no sample column and duplicated cell basenames; cannot join unambiguously.")
  }

  pred_min <- pred2 %>%
    transmute(
      Base = .data$Base,
      ash_cell = .data$cell,
      ash_sample = NA_character_,
      ash_label = .data$prediction,
      ash_prob = .data$probability
    )

  merged <- final %>% left_join(pred_min, by = "Base")
  message("[merge_ashleys] Matched predictions by basename: ", sum(!is.na(merged$ash_cell)), " / ", nrow(merged))
}

# ---- Join features (optional) ----
if (!is.null(feat_path) && file.exists(feat_path) && file.info(feat_path)$size > 0) {
  message("[merge_ashleys] Joining Ashley features ...")
  feat <- read_tsv_quiet(feat_path)
  feat <- normalize_feat_keys(feat)

  if (!is.null(feat)) {
    # Prefix all feature columns except join key; normalize names
    feat_cols <- setdiff(names(feat), c("Base"))
    feat_pref <- feat %>%
      mutate(Base = strip_bam(.data$Base)) %>%
      rename_with(.cols = all_of(feat_cols), .fn = ~ paste0("ash_", normalize_ash_feature_name(.x)))

    if (any(duplicated(feat_pref$Base))) {
      warning("[merge_ashleys] features.tsv has duplicated cell basenames; keeping the first feature row per basename.")
      feat_pref <- feat_pref %>% distinct(.data$Base, .keep_all = TRUE)
    }

    merged <- merged %>% left_join(feat_pref, by = "Base")
    message("[merge_ashleys] Feature columns added: ", length(setdiff(names(feat_pref), "Base")))
  } else {
    message("[merge_ashleys] Skipping features merge (no usable key).")
  }
} else {
  message("[merge_ashleys] No features file supplied or empty; skipping features merge.")
}

# Drop helper columns if you don't want them in final output
# Comment these out if you prefer to keep debug columns.
merged <- merged %>% select(-Library_raw, -Sample_raw, -Base, -Library_key)

# Atomic write
out_path <- normalizePath(args$out, mustWork = FALSE)
out_dir  <- dirname(out_path)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("[merge_ashleys] Writing to: ", out_path)
tmp_path <- paste0(out_path, ".tmp")
readr::write_tsv(merged, tmp_path)

if (!file.exists(tmp_path) || file.size(tmp_path) <= 0) {
  stop("[merge_ashleys] ERROR: tmp file missing/empty: ", tmp_path)
}
ok <- file.rename(tmp_path, out_path)
if (!isTRUE(ok)) {
  file.copy(tmp_path, out_path, overwrite = TRUE)
  unlink(tmp_path)
}
if (!file.exists(out_path) || file.size(out_path) <= 0) {
  stop("[merge_ashleys] ERROR: final file missing/empty: ", out_path)
}
message("[merge_ashleys] Done.")
