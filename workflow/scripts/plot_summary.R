#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# --- helper: SI axis labels (replacement for deprecated label_number_si) ---
lab_si <- function(accuracy = NULL) {
  scales::label_number(scale_cut = scales::cut_si(""), accuracy = accuracy)
}

# optional but nice
have_met <- requireNamespace("MetBrewer", quietly = TRUE)
have_vir <- requireNamespace("viridisLite", quietly = TRUE)

args <- list()
argv <- commandArgs(trailingOnly = TRUE)
for (i in seq(1, length(argv), by = 2)) {
  key <- gsub("^--", "", argv[i])
  val <- if ((i + 1) <= length(argv)) argv[i + 1] else TRUE
  args[[key]] <- val
}

# Required args
if (is.null(args$final) || is.null(args$out)) {
  stop("Usage: plot_summary.R --final final_qc.tsv --out summary.pdf [--count-col 'bin_total_read_count']")
}

final_path <- args$final
out_pdf    <- args$out
count_col  <- if (!is.null(args[["count-col"]])) args[["count-col"]] else "bin_total_read_count"

# quiet readr column-spec chatter
df <- suppressMessages(
  readr::read_tsv(final_path, progress = FALSE, guess_max = 1e6, show_col_types = FALSE)
)

get_col_any <- function(candidates, default = NA_real_) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0) df[[hit[1]]] else rep(default, nrow(df))
}

# pick palette
pal_disc <- if (have_met) MetBrewer::met.brewer("Hokusai3", 9, direction = -1) else hue_pal()(9)
pal_fill <- pal_disc[6]
pal_line <- pal_disc[8]
pal_grad <- if (have_met) {
  MetBrewer::met.brewer("Hokusai3", 9)
} else if (have_vir) {
  viridisLite::viridis(9)
} else {
  colorRampPalette(hue_pal()(3))(9)
}

theme_run <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# Normalise current prefixed output names, with legacy fallbacks for older runs.
TotalReadCount <- get_col_any(c(count_col, "bin_total_read_count", "total.read.count"))
DuplicateFraction <- get_col_any(c("alf_duplicate_frac", "DuplicateFraction"))
Entropy    <- get_col_any(c("bin_entropy", "entropy"))
Spikiness  <- get_col_any(c("bin_spikiness", "spikiness"))

MappedFraction       <- get_col_any(c("alf_mapped_frac", "MappedFraction"))
MappedProperFraction <- get_col_any(c("alf_mapped_proper_pair_frac", "MappedProperFraction"))
MedianCoverage       <- get_col_any(c("alf_coverage_med", "MedianCoverage"))
MedianInsertSize     <- get_col_any(c("alf_insert_size_med", "MedianInsertSize"))
ErrorRate            <- get_col_any(c("alf_error_frac", "ErrorRate"))
GCR                  <- get_col_any(c("bin_gc_r", "gc_pearson_r"))
PreseqSat            <- get_col_any(c("preseq_saturation"))

# ---------- PAGE 1 ----------
# A) jitter+violin+box of TotalReadCount
df_plot <- tibble(group = "All libraries", y = as.numeric(TotalReadCount)) %>% filter(!is.na(y))

p_total <-
  ggplot(df_plot, aes(x = group, y = y)) +
  geom_violin(fill = pal_fill, alpha = 0.65, color = NA, width = 0.85, na.rm = TRUE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", na.rm = TRUE) +
  geom_jitter(width = 0.08, height = 0, size = 0.9, alpha = 0.35, na.rm = TRUE) +
  scale_y_continuous(labels = lab_si()) +
  labs(x = NULL, y = "Total number of reads per cell", title = "Total read count per cell") +
  theme_run + theme(axis.text.x = element_blank())

# B) duplicate rate histogram
p_dup <-
  ggplot(tibble(x = as.numeric(DuplicateFraction)) %>% filter(!is.na(x)),
         aes(x = x)) +
  geom_histogram(bins = 40, fill = pal_fill, alpha = 0.85, color = NA) +
  scale_x_continuous(labels = label_percent(accuracy = 1)) +
  labs(x = "Duplicate rate", y = "count", title = "Duplicate rate") +
  theme_run

# C) spikiness vs entropy colored by total reads
p_se <-
  ggplot(tibble(Entropy = as.numeric(Entropy),
                Spikiness = as.numeric(Spikiness),
                Reads = as.numeric(TotalReadCount)) %>%
           filter(!is.na(Entropy), !is.na(Spikiness), !is.na(Reads)),
         aes(Entropy, Spikiness, color = Reads)) +
  geom_point(alpha = 0.35, size = 1.2) +
  scale_color_gradientn(colors = pal_grad, labels = lab_si()) +
  labs(x = "Entropy", y = "Spikiness",
       color = "Total reads",
       title = "Spikiness vs Entropy") +
  theme_run

page1 <- (p_total | p_dup) / p_se +
  plot_annotation(title = "Run overview", theme = theme(plot.title = element_text(face = "bold", size = 14)))

# ---------- PAGE 2 ----------
plots2 <- list()

# helper builders that no-op gracefully when data missing
make_hist <- function(vec, title, xlab, to_percent = FALSE, bins = 40) {
  dat <- tibble(x = as.numeric(vec)) %>% filter(!is.na(x))
  if (nrow(dat) == 0) return(ggplot() + theme_void() + ggtitle(paste("No data:", title)))
  g <- ggplot(dat, aes(x)) +
    geom_histogram(bins = bins, fill = pal_fill, alpha = 0.85, color = NA) +
    labs(title = title, x = xlab, y = "count") + theme_run
  if (to_percent) g <- g + scale_x_continuous(labels = label_percent())
  g
}
make_jvb <- function(vec, title, ylab, log10y = FALSE) {
  dat <- tibble(group = "All libraries", y = as.numeric(vec)) %>% filter(!is.na(y))
  if (nrow(dat) == 0) return(ggplot() + theme_void() + ggtitle(paste("No data:", title)))
  g <- ggplot(dat, aes(group, y)) +
    geom_violin(fill = pal_fill, alpha = 0.65, color = NA, width = 0.85) +
    geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white") +
    geom_jitter(width = 0.08, height = 0, size = 0.9, alpha = 0.35) +
    labs(title = title, x = NULL, y = ylab) +
    theme_run + theme(axis.text.x = element_blank())
  if (log10y) g <- g + scale_y_continuous(trans = "log10", labels = lab_si())
  g
}
make_scatter <- function(xv, yv, title, xlab, ylab) {
  dat <- tibble(x = as.numeric(xv), y = as.numeric(yv)) %>% filter(!is.na(x), !is.na(y))
  if (nrow(dat) == 0) return(ggplot() + theme_void() + ggtitle(paste("No data:", title)))
  ggplot(dat, aes(x, y)) +
    geom_point(alpha = 0.35, size = 1.2, color = pal_line) +
    labs(title = title, x = xlab, y = ylab) + theme_run
}

plots2 <- c(
  list(make_hist(MappedFraction, "Mapped fraction", "fraction", to_percent = TRUE)),
  list(make_hist(MappedProperFraction, "Proper pairs", "fraction", to_percent = TRUE)),
  list(make_hist(ErrorRate, "Error rate", "fraction", to_percent = TRUE)),
  list(make_jvb(MedianCoverage, "Median coverage", "x (log10)", log10y = TRUE)),
  list(make_jvb(MedianInsertSize, "Median insert size", "bp", log10y = FALSE)),
  list(make_hist(GCR, "GC correlation", "Pearson r")),
  list(make_hist(PreseqSat, "Preseq saturation", "fraction", to_percent = TRUE))
)

# keep only first 6 for a clean grid; extra will start another row
page2 <- wrap_plots(plots2[1:min(6, length(plots2))], ncol = 3) +
  plot_annotation(title = "Additional run metrics")

# write multi-page pdf
pdf(out_pdf, width = 11, height = 8.5, onefile = TRUE)
print(page1)
print(page2)
# if more than 6 extra panels exist, put them on a 3rd page
if (length(plots2) > 6) {
  print(wrap_plots(plots2[7:length(plots2)], ncol = 3) +
          plot_annotation(title = "Additional run metrics (cont.)"))
}
invisible(dev.off())
