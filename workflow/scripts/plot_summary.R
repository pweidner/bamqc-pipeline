#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

lab_si <- function(accuracy = NULL) {
  scales::label_number(scale_cut = scales::cut_si(""), accuracy = accuracy)
}

args <- list()
argv <- commandArgs(trailingOnly = TRUE)
for (i in seq(1, length(argv), by = 2)) {
  key <- gsub("^--", "", argv[i])
  val <- if ((i + 1) <= length(argv)) argv[i + 1] else TRUE
  args[[key]] <- val
}

if (is.null(args$final) || is.null(args$out)) {
  stop("Usage: plot_summary.R --final final_qc.tsv --out summary.pdf [--count-col 'bin_total_read_count']")
}

final_path <- args$final
out_pdf <- args$out
count_col <- if (!is.null(args[["count-col"]])) args[["count-col"]] else "bin_total_read_count"

df <- suppressMessages(
  readr::read_tsv(final_path, progress = FALSE, guess_max = 1e6, show_col_types = FALSE)
)

as_num <- function(x) suppressWarnings(as.numeric(x))
finite <- function(x) is.finite(as_num(x))
median_finite <- function(x) {
  x <- as_num(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

get_col_any <- function(candidates, default = NA_real_) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0) df[[hit[1]]] else rep(default, nrow(df))
}

get_chr_any <- function(candidates, default = NA_character_) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0) as.character(df[[hit[1]]]) else rep(default, nrow(df))
}

rank01 <- function(x, high_good = TRUE) {
  x <- as_num(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  vals <- x[ok]
  if (length(unique(vals)) < 2) {
    out[ok] <- 0.5
  } else {
    out[ok] <- (rank(vals, ties.method = "average") - 1) / (length(vals) - 1)
  }
  if (!high_good) out[ok] <- 1 - out[ok]
  out
}

safe_spearman <- function(x, y) {
  x <- as_num(x); y <- as_num(y)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

rho_label <- function(x, y) {
  r <- safe_spearman(x, y)
  if (is.na(r)) "" else paste0("Spearman rho = ", sprintf("%.2f", r))
}

ink <- "#202124"
muted <- "#6b7280"
axis_grey <- "#3f3f46"
paper <- "#ffffff"
panel_line <- "#bfc4c9"
blue <- "#2f6f9f"
teal <- "#1b8a78"
vermilion <- "#b54736"
gold <- "#c7972e"
slate <- "#475569"

seq_pal <- c("#f4f1ec", "#d8d2c4", "#9bb7bd", "#5a8ca3", "#2f5f85", "#223b5a")
div_pal <- c("#2f6f9f", "#f7f7f7", "#b54736")

theme_bamqc <- theme_classic(base_size = 9) +
  theme(
    plot.background = element_rect(fill = paper, color = NA),
    panel.background = element_rect(fill = paper, color = NA),
    axis.line = element_line(color = axis_grey, linewidth = 0.3),
    axis.ticks = element_line(color = axis_grey, linewidth = 0.25),
    axis.text = element_text(color = ink, size = 8),
    axis.title = element_text(color = ink, size = 9),
    plot.title = element_text(face = "bold", color = ink, size = 10, margin = margin(b = 3)),
    plot.subtitle = element_text(color = muted, size = 8, margin = margin(b = 4)),
    legend.title = element_text(size = 8, color = ink),
    legend.text = element_text(size = 7, color = ink),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(4, "mm"),
    plot.tag = element_text(face = "bold", size = 11, color = ink),
    plot.margin = margin(5, 7, 5, 5)
  )

empty_panel <- function(title, subtitle = "Metric not available in this run") {
  ggplot() +
    annotate("text", x = 0, y = 0, label = subtitle, size = 3, color = muted) +
    xlim(-1, 1) + ylim(-1, 1) +
    labs(title = title, x = NULL, y = NULL) +
    theme_bamqc +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
}

Library <- get_chr_any(c("Library"), default = seq_len(nrow(df)))
Sample <- get_chr_any(c("Sample", "sctrip_sample"), default = "All libraries")

TotalReadCount <- as_num(get_col_any(c(count_col, "bin_total_read_count", "total.read.count")))
SctripGood <- as_num(get_col_any(c("sctrip_good")))
UsableReads <- ifelse(is.finite(SctripGood), SctripGood, TotalReadCount)
DuplicateFraction <- as_num(get_col_any(c("alf_duplicate_frac", "DuplicateFraction")))
Entropy <- as_num(get_col_any(c("bin_entropy", "entropy")))
Spikiness <- as_num(get_col_any(c("bin_spikiness", "spikiness")))
MappedFraction <- as_num(get_col_any(c("alf_mapped_frac", "MappedFraction")))
MappedProperFraction <- as_num(get_col_any(c("alf_mapped_proper_pair_frac", "MappedProperFraction")))
GCR <- as_num(get_col_any(c("bin_gc_r", "gc_pearson_r")))
PreseqSat <- as_num(get_col_any(c("preseq_saturation")))
Background <- as_num(get_col_any(c("bin_background", "background")))
Gini <- as_num(get_col_any(c("bin_gini", "coverage_gini")))
Fold80 <- as_num(get_col_any(c("bin_fold80", "fold80_penalty")))
CoverageCV <- as_num(get_col_any(c("bin_cv", "coverage_cv")))

score_mat <- cbind(
  usable_reads = rank01(log10(pmax(UsableReads, 1)), high_good = TRUE),
  entropy = rank01(Entropy, high_good = TRUE),
  mapped = rank01(MappedFraction, high_good = TRUE),
  preseq = rank01(PreseqSat, high_good = TRUE),
  spikiness = rank01(Spikiness, high_good = FALSE),
  duplicates = rank01(DuplicateFraction, high_good = FALSE),
  background = rank01(Background, high_good = FALSE),
  gc_bias = rank01(abs(GCR), high_good = FALSE)
)
QualityScore <- rowMeans(score_mat, na.rm = TRUE)
QualityScore[!is.finite(QualityScore)] <- NA_real_

qc <- tibble(
  Library = Library,
  Sample = Sample,
  UsableReads = UsableReads,
  TotalReadCount = TotalReadCount,
  SctripGood = SctripGood,
  DuplicateFraction = DuplicateFraction,
  Entropy = Entropy,
  Spikiness = Spikiness,
  MappedFraction = MappedFraction,
  MappedProperFraction = MappedProperFraction,
  GCR = GCR,
  AbsGCR = abs(GCR),
  PreseqSat = PreseqSat,
  Background = Background,
  Gini = Gini,
  Fold80 = Fold80,
  CoverageCV = CoverageCV,
  QualityScore = QualityScore
)

make_depth_spike <- function() {
  dat <- qc %>% filter(is.finite(UsableReads), UsableReads > 0, is.finite(Spikiness))
  if (nrow(dat) == 0) return(empty_panel("Usable reads and bin spikiness"))
  color_var <- if (any(is.finite(dat$Background))) "Background" else "QualityScore"
  g <- ggplot(dat, aes(UsableReads, Spikiness, color = .data[[color_var]])) +
    geom_point(alpha = 0.72, size = 1.35, stroke = 0, na.rm = TRUE) +
    scale_x_continuous(trans = "log10", labels = lab_si()) +
    labs(
      title = "Usable reads and bin spikiness",
      subtitle = rho_label(log10(dat$UsableReads), dat$Spikiness),
      x = "Usable reads per cell",
      y = "bin_spikiness",
      color = if (color_var == "Background") "bin_background" else "QC score"
    ) +
    theme_bamqc
  if (color_var == "Background") {
    g + scale_color_gradient(low = blue, high = vermilion, na.value = "#d4d4d8")
  } else {
    g + scale_color_gradient(low = "#d4d4d8", high = teal, na.value = "#d4d4d8")
  }
}

make_background_spike <- function() {
  dat <- qc %>% filter(is.finite(Background), is.finite(Spikiness), is.finite(UsableReads), UsableReads > 0)
  if (nrow(dat) == 0) return(empty_panel("Background and coverage shape"))
  g <- ggplot(dat, aes(Background, Spikiness, color = UsableReads)) +
    geom_point(alpha = 0.72, size = 1.35, stroke = 0, na.rm = TRUE) +
    labs(
      title = "Background and coverage shape",
      subtitle = rho_label(dat$Background, dat$Spikiness),
      x = "bin_background",
      y = "bin_spikiness",
      color = "Usable reads"
    ) +
    scale_color_gradientn(colors = seq_pal, labels = lab_si(), na.value = "#d4d4d8") +
    theme_bamqc
  if (nrow(dat) > 2) g <- g + geom_smooth(method = "lm", se = FALSE, linewidth = 0.35, color = ink, na.rm = TRUE)
  g
}

make_dup_preseq <- function() {
  dat <- qc %>% filter(is.finite(DuplicateFraction), is.finite(PreseqSat))
  if (nrow(dat) == 0) return(empty_panel("Duplication and library complexity"))
  g <- ggplot(dat, aes(DuplicateFraction, PreseqSat, color = MappedFraction)) +
    geom_point(alpha = 0.72, size = 1.35, stroke = 0, na.rm = TRUE) +
    scale_x_continuous(labels = label_percent(accuracy = 1)) +
    scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0, NA)) +
    scale_color_gradient(low = gold, high = blue, na.value = "#d4d4d8", labels = label_percent(accuracy = 1)) +
    labs(
      title = "Duplication and library complexity",
      subtitle = rho_label(dat$DuplicateFraction, dat$PreseqSat),
      x = "Duplicate fraction",
      y = "preseq_saturation",
      color = "Mapped"
    ) +
    theme_bamqc
  if (nrow(dat) > 2) g <- g + geom_smooth(method = "lm", se = FALSE, linewidth = 0.35, color = ink, na.rm = TRUE)
  g
}

make_quality_rank <- function() {
  dat <- qc %>% filter(is.finite(QualityScore)) %>% arrange(desc(QualityScore)) %>% mutate(Rank = row_number())
  if (nrow(dat) == 0) return(empty_panel("QC score ranking"))
  ggplot(dat, aes(Rank, QualityScore, color = Background)) +
    geom_line(color = panel_line, linewidth = 0.25, na.rm = TRUE) +
    geom_point(alpha = 0.78, size = 1.3, stroke = 0, na.rm = TRUE) +
    scale_y_continuous(labels = label_number(accuracy = 0.01), limits = c(0, 1)) +
    scale_color_gradient(low = blue, high = vermilion, na.value = slate) +
    labs(
      title = "Ranked multi-metric QC score",
      subtitle = "Ranks high depth, high complexity, low noise, low background",
      x = "Library rank",
      y = "QC score",
      color = "bin_background"
    ) +
    theme_bamqc
}

make_gc_uniformity <- function() {
  dat <- qc %>% filter(is.finite(AbsGCR), is.finite(CoverageCV))
  if (nrow(dat) == 0) return(empty_panel("GC bias and coverage variability"))
  ggplot(dat, aes(AbsGCR, CoverageCV, color = Spikiness)) +
    geom_point(alpha = 0.72, size = 1.35, stroke = 0, na.rm = TRUE) +
    scale_color_gradient(low = teal, high = vermilion, na.value = "#d4d4d8") +
    labs(
      title = "GC bias and coverage variability",
      subtitle = rho_label(dat$AbsGCR, dat$CoverageCV),
      x = "abs(bin_gc_r)",
      y = "bin_cv",
      color = "bin_spikiness"
    ) +
    theme_bamqc
}

make_background_hist <- function() {
  dat <- qc %>% filter(is.finite(Background))
  if (nrow(dat) == 0) return(empty_panel("Background distribution"))
  ggplot(dat, aes(Background)) +
    geom_histogram(bins = 35, fill = blue, alpha = 0.88, color = paper, linewidth = 0.15, na.rm = TRUE) +
    labs(
      title = "Background distribution",
      subtitle = paste0("Median = ", sprintf("%.4f", median_finite(dat$Background))),
      x = "bin_background",
      y = "Libraries"
    ) +
    theme_bamqc
}

make_sample_summary <- function() {
  dat <- qc %>%
    filter(!is.na(Sample), Sample != "") %>%
    group_by(Sample) %>%
    summarize(
      n = dplyr::n(),
      median_score = median_finite(QualityScore),
      median_usable = median_finite(UsableReads),
      median_background = median_finite(Background),
      .groups = "drop"
    ) %>%
    filter(is.finite(median_usable), median_usable > 0, is.finite(median_score))
  if (nrow(dat) == 0 || length(unique(dat$Sample)) < 2) return(empty_panel("Sample-level overview"))
  ggplot(dat, aes(median_usable, median_score, size = n, fill = median_background)) +
    geom_point(shape = 21, color = ink, alpha = 0.78, stroke = 0.25, na.rm = TRUE) +
    scale_x_continuous(trans = "log10", labels = lab_si()) +
    scale_y_continuous(limits = c(0, 1), labels = label_number(accuracy = 0.01)) +
    scale_size_continuous(range = c(1.8, 7), breaks = pretty_breaks(3)) +
    scale_fill_gradient(low = blue, high = vermilion, na.value = "#d4d4d8") +
    labs(
      title = "Sample-level overview",
      subtitle = "Point size is library count",
      x = "Median usable reads",
      y = "Median QC score",
      size = "Libraries",
      fill = "Median background"
    ) +
    theme_bamqc
}

make_correlation_heatmap <- function() {
  metrics <- data.frame(
    `usable reads` = log10(pmax(qc$UsableReads, 1)),
    `bin_background` = qc$Background,
    `bin_spikiness` = qc$Spikiness,
    `bin_entropy` = qc$Entropy,
    `duplicate fraction` = qc$DuplicateFraction,
    `preseq_saturation` = qc$PreseqSat,
    `mapped fraction` = qc$MappedFraction,
    `abs GC correlation` = qc$AbsGCR,
    `bin_cv` = qc$CoverageCV,
    check.names = FALSE
  )
  keep <- vapply(metrics, function(x) {
    x <- as_num(x)
    sum(is.finite(x)) >= 3 && length(unique(x[is.finite(x)])) > 1
  }, logical(1))
  metrics <- metrics[, keep, drop = FALSE]
  if (ncol(metrics) < 2) return(empty_panel("Spearman metric correlations"))
  cmat <- suppressWarnings(cor(metrics, use = "pairwise.complete.obs", method = "spearman"))
  cdat <- as.data.frame(as.table(cmat), stringsAsFactors = FALSE)
  names(cdat) <- c("MetricX", "MetricY", "rho")
  levels <- colnames(metrics)
  cdat$MetricX <- factor(cdat$MetricX, levels = levels)
  cdat$MetricY <- factor(cdat$MetricY, levels = rev(levels))
  ggplot(cdat, aes(MetricX, MetricY, fill = rho)) +
    geom_tile(color = paper, linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 2.15, color = ink, na.rm = TRUE) +
    scale_fill_gradient2(low = div_pal[1], mid = div_pal[2], high = div_pal[3], midpoint = 0, limits = c(-1, 1), na.value = "#e5e7eb") +
    labs(title = "Spearman metric correlations", x = NULL, y = NULL, fill = "rho") +
    theme_bamqc +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
}

page1 <- (make_depth_spike() | make_background_spike()) /
  (make_dup_preseq() | make_quality_rank()) +
  plot_annotation(
    title = "bamqc-pipeline run summary",
    subtitle = paste0(nrow(qc), " libraries; lower bin_background, lower spikiness, and higher usable reads mark stronger libraries"),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, color = ink),
      plot.subtitle = element_text(size = 9, color = muted)
    )
  )

page2 <- (make_correlation_heatmap() | make_gc_uniformity()) /
  (make_background_hist() | make_sample_summary()) +
  plot_annotation(
    title = "QC metric relationships",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 14, color = ink))
  )

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
pdf(out_pdf, width = 10.5, height = 8.0, onefile = TRUE, useDingbats = FALSE)
print(page1)
print(page2)
invisible(dev.off())
