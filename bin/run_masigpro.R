#!/usr/bin/env Rscript
# bin/run_masigpro.R
# Pseudobulk maSigPro DEA per ROI, restricted to a configurable sample subset.

suppressPackageStartupMessages({
  library(optparse)
  library(maSigPro)
  library(tidyverse)
  library(edgeR)
})

option_list <- list(
  make_option("--delineation_dir", type = "character"),
  make_option("--counts_dir", type = "character"),
  make_option("--sample_metadata", type = "character",
              help = "CSV with columns: sample_id,condition,time"),
  make_option("--gene_list", type = "character", default = "NO_FILE",
              help = "Optional newline-delimited gene list file(s), comma-separated"),
  make_option("--samples", type = "character", default = "ALL",
              help = "Comma-separated sample_ids to include, or ALL"),
  make_option("--rois", type = "character", default = "Dorsal_gm,Medial_gm,Ventral_gm",
              help = "Comma-separated ROI names matching delineation CSV values"),
  make_option("--outdir", type = "character", default = "masigpro_results")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Sample metadata + optional restriction
# -----------------------------
sample_info <- read_csv(opt$sample_metadata, show_col_types = FALSE)

if (opt$samples != "ALL") {
  keep_samples <- strsplit(opt$samples, ",")[[1]]
  sample_info <- sample_info %>% filter(sample_id %in% keep_samples)
}
stopifnot(nrow(sample_info) > 0)

# -----------------------------
# Gene list(s) of interest (optional — if NO_FILE, use all genes)
# -----------------------------
use_gene_filter <- opt$gene_list != "NO_FILE"
if (use_gene_filter) {
  gene_files <- strsplit(opt$gene_list, ",")[[1]]
  genes_of_interest <- unique(unlist(lapply(gene_files, readLines)))
  cat("Loaded", length(genes_of_interest), "genes of interest\n")
}

roi_list <- strsplit(opt$rois, ",")[[1]]

for (roi in roi_list) {
  cat("Running ROI:", roi, "\n")

  pb_list <- list()
  meta_list <- list()

  for (i in seq_len(nrow(sample_info))) {
    sample <- sample_info$sample_id[i]
    cat("  Sample:", sample, "\n")

    mn_file <- file.path(opt$delineation_dir, paste0(sample, "_manual_delineation.csv"))
    if (!file.exists(mn_file)) { warning("Missing delineation file for ", sample); next }
    mn_data <- read.csv(mn_file)
    mn_barcodes <- sub(".*_", "", mn_data$X[mn_data$manual_delineation == roi])

    counts_file <- file.path(opt$counts_dir, paste0(sample, "_counts.csv"))
    if (!file.exists(counts_file)) { warning("Missing counts file for ", sample); next }
    counts_data <- read.csv(counts_file, row.names = 1, check.names = FALSE)
    colnames(counts_data) <- sub(".*_", "", colnames(counts_data))

    mn_barcodes <- intersect(mn_barcodes, colnames(counts_data))
    if (length(mn_barcodes) == 0) { warning("No cells in ROI for ", sample); next }

    mn_counts <- counts_data[, mn_barcodes, drop = FALSE]
    pb <- matrix(rowSums(mn_counts), ncol = 1)
    rownames(pb) <- rownames(mn_counts)
    colnames(pb) <- sample
    pb_list[[sample]] <- pb

    meta_list[[sample]] <- tibble(
      sample = sample,
      condition = sample_info$condition[i],
      time = sample_info$time[i]
    )
  }

  if (length(pb_list) < 2) {
    warning("Not enough samples with data for ROI ", roi, " — skipping")
    next
  }

  all_genes <- Reduce(union, lapply(pb_list, rownames))
  pb_list_filled <- lapply(pb_list, function(pb) {
    missing <- setdiff(all_genes, rownames(pb))
    if (length(missing) > 0) {
      zeros <- matrix(0, nrow = length(missing), ncol = 1)
      rownames(zeros) <- missing; colnames(zeros) <- colnames(pb)
      pb <- rbind(pb, zeros)
    }
    pb[all_genes, , drop = FALSE]
  })

  expr_roi <- Reduce(cbind, pb_list_filled)
  meta_roi <- bind_rows(meta_list)
  meta_roi <- meta_roi[match(colnames(expr_roi), meta_roi$sample), ]

  # edgeR filter + TMM normalization
  dge <- DGEList(counts = expr_roi)
  group <- interaction(meta_roi$condition, meta_roi$time)
  dge <- dge[filterByExpr(dge, group = group), , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")

  if (use_gene_filter) {
    dge <- dge[intersect(rownames(dge), genes_of_interest), ]
  }
  cat("  Testing", nrow(dge), "genes\n")

  expr_round <- round(cpm(dge, log = FALSE))
  roi_dir <- file.path(opt$outdir, paste0("roi_", roi))
  dir.create(roi_dir, recursive = TRUE, showWarnings = FALSE)

  write.table(t(rownames(dge)),
              file = file.path(roi_dir, paste0("initial_filtered_genes_", roi, ".txt")),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  # -----------------------------
  # maSigPro design — generalized to whatever conditions/timepoints are present
  # -----------------------------
  meta_roi$rep_group <- as.numeric(as.factor(paste(meta_roi$condition, meta_roi$time)))
  conditions_present <- unique(meta_roi$condition)

  edesign <- data.frame(Time = meta_roi$time, Replicate = meta_roi$rep_group)
  for (cond in setdiff(conditions_present, "healthy")) {
    edesign[[cond]] <- as.numeric(meta_roi$condition %in% c(cond, "healthy"))
  }
  rownames(edesign) <- meta_roi$sample

  design.matrix <- make.design.matrix(edesign, degree = 2)
  rownames(design.matrix$dis) <- rownames(edesign)

  fit <- p.vector(expr_round, design.matrix, counts = TRUE, Q = 0.05)
  tfit <- T.fit(fit)
  sigs <- get.siggenes(tfit, rsq = 0.4, vars = "groups")

  for (contrast_name in names(sigs$sig.genes)) {
    contrast <- sigs$sig.genes[[contrast_name]]

    pdf(file.path(roi_dir, paste0("gene_profiles_", roi, "_", contrast_name, ".pdf")),
        width = 24, height = 12)
    cluster <- see.genes(contrast, show.fit = TRUE, dis = design.matrix$dis,
                          cluster.method = "hclust", cluster.data = 1, k = 9)
    dev.off()

    write.table(cluster$cut,
                file = file.path(roi_dir, paste0("output_cluster_", roi, "_", contrast_name, ".txt")),
                sep = "\t", quote = FALSE)
    write.table(contrast$sig.profiles,
                file = file.path(roi_dir, paste0("output_sig_profiles_", roi, "_", contrast_name, ".txt")),
                sep = "\t", quote = FALSE)
    write.table(contrast$sig.pvalues,
                file = file.path(roi_dir, paste0("output_sig_pvalues_", roi, "_", contrast_name, ".txt")),
                sep = "\t", quote = FALSE)
  }
}
