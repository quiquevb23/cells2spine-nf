#!/usr/bin/env Rscript

# ==========================================================
# CLI MaSigPro R PIPELINE FOR NEXTFLOW
# ==========================================================

suppressPackageStartupMessages({
  library(optparse)
  library(maSigPro)
  library(tidyverse)
  library(edgeR)
})

# Define CLI options
option_list <- list(
  make_option(c("-d", "--delineation_dir"), type = "character", default = NULL,
              help = "Path to directory containing *_manual_delineation.csv files", metavar = "character"),
  make_option(c("-c", "--counts_dir"), type = "character", default = NULL,
              help = "Path to directory containing *_counts.csv files", metavar = "character"),
  make_option(c("-m", "--sample_metadata"), type = "character", default = NULL,
              help = "Path to sample metadata CSV file (cols: sample_id, condition, time)", metavar = "character"),
  make_option(c("-g", "--gene_list"), type = "character", default = "NO_FILE",
              help = "Path to optional txt file or directory with gene lists", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = "masigpro_results",
              help = "Output directory [default= %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$delineation_dir) || is.null(opt$counts_dir) || is.null(opt$sample_metadata)) {
  print_help(opt_parser)
  stop("Error: --delineation_dir, --counts_dir, and --sample_metadata are required arguments.", call. = FALSE)
}

# -----------------------------
# 1. Read Dynamic Metadata
# -----------------------------
sample_info <- read.csv(opt$sample_metadata, stringsAsFactors = FALSE)

required_cols <- c("sample_id", "condition", "time")
if (!all(required_cols %in% colnames(sample_info))) {
  stop(paste("Error: --sample_metadata CSV must contain columns:", paste(required_cols, collapse = ", ")))
}

cat("Loaded metadata for", nrow(sample_info), "samples.\n")

# -----------------------------
# 2. Parse Genes of Interest (Optional)
# -----------------------------
genes_of_interest <- NULL
if (opt$gene_list != "NO_FILE" && file.exists(opt$gene_list)) {
  if (dir.exists(opt$gene_list)) {
    gene_files <- list.files(opt$gene_list, full.names = TRUE)
    gene_sets <- lapply(gene_files, readLines)
    genes_of_interest <- unique(unlist(gene_sets))
  } else {
    genes_of_interest <- unique(readLines(opt$gene_list))
  }
  cat("Loaded", length(genes_of_interest), "genes of interest.\n")
} else {
  cat("No gene filter provided. Running maSigPro on all expressed genes.\n")
}

# -----------------------------
# 3. Detect ROIs Dynamically
# -----------------------------
delineation_files <- list.files(opt$delineation_dir, pattern = "*_manual_delineation.csv$", full.names = TRUE)
if (length(delineation_files) == 0) {
  stop("No delineation files found matching pattern '*_manual_delineation.csv'")
}

first_del <- read.csv(delineation_files[1])
roi_list <- unique(first_del$manual_delineation)
roi_list <- roi_list[!is.na(roi_list) & roi_list != ""]

cat("Detected ROIs:", paste(roi_list, collapse = ", "), "\n")

# -----------------------------
# 4. Iterate over ROIs
# -----------------------------
for (roi in roi_list) {
  cat("\n==========================================\n")
  cat("Running ROI:", roi, "\n")
  cat("==========================================\n")
  
  pb_list <- list()
  meta_list <- list()
  
  for (i in seq_len(nrow(sample_info))) {
    sample <- sample_info$sample_id[i]
    
    mn_file <- file.path(opt$delineation_dir, paste0(sample, "_manual_delineation.csv"))
    counts_file <- file.path(opt$counts_dir, paste0(sample, "_counts.csv"))
    
    if (!file.exists(mn_file) || !file.exists(counts_file)) {
      warning(paste("Missing files for sample", sample, "- skipping."))
      next
    }
    
    mn_data <- read.csv(mn_file)
    mn_barcodes <- mn_data$X[mn_data$manual_delineation == roi]
    mn_barcodes <- sub(".*_", "", mn_barcodes)
    
    counts_data <- read.csv(counts_file, row.names = 1, check.names = FALSE)
    colnames(counts_data) <- sub(".*_", "", colnames(counts_data))
    
    mn_barcodes <- intersect(mn_barcodes, colnames(counts_data))
    
    if (length(mn_barcodes) == 0) {
      next
    }
    
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
  
  if (length(pb_list) == 0) {
    cat("No samples found containing cells for ROI:", roi, "\n")
    next
  }
  
  # Merge Pseudobulk
  all_genes <- Reduce(union, lapply(pb_list, rownames))
  
  pb_list_filled <- lapply(pb_list, function(pb) {
    missing <- setdiff(all_genes, rownames(pb))
    if (length(missing) > 0) {
      zeros <- matrix(0, nrow = length(missing), ncol = 1)
      rownames(zeros) <- missing
      colnames(zeros) <- colnames(pb)
      pb <- rbind(pb, zeros)
    }
    pb[all_genes, , drop = FALSE]
  })
  
  expr_roi <- Reduce(cbind, pb_list_filled)
  meta_roi <- bind_rows(meta_list)
  meta_roi <- meta_roi[match(colnames(expr_roi), meta_roi$sample), ]
  
  # Filtering with edgeR
  dge <- DGEList(counts = expr_roi)
  group <- interaction(meta_roi$condition, meta_roi$time)
  keep <- filterByExpr(dge, group = group)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Apply optional gene set filter
  if (!is.null(genes_of_interest)) {
    genes_keep <- intersect(rownames(dge), genes_of_interest)
    dge <- dge[genes_keep, ]
  }
  
  cat("Testing", nrow(dge), "genes in ROI:", roi, "\n")
  if (nrow(dge) < 3) {
    cat("Skipping maSigPro fit for ROI", roi, "- insufficient genes after filtering.\n")
    next
  }
  
  expr_norm <- cpm(dge, log = FALSE)
  expr_round <- round(expr_norm)
  
  roi_dir <- file.path(opt$outdir, paste0("roi_", roi))
  dir.create(roi_dir, recursive = TRUE, showWarnings = FALSE)
  
  write.table(t(rownames(dge)),
              file = file.path(roi_dir, paste0("initial_filtered_genes_", roi, ".txt")),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # MaSigPro Modeling
  meta_roi$rep_group <- as.numeric(as.factor(paste(meta_roi$condition, meta_roi$time)))
  
  edesign <- data.frame(
    Time = meta_roi$time,
    Replicate = meta_roi$rep_group,
    SCI = as.numeric(meta_roi$condition %in% c("SCI", "healthy")),
    bPAC = as.numeric(meta_roi$condition %in% c("bPAC", "healthy"))
  )
  rownames(edesign) <- meta_roi$sample
  
  design.matrix <- make.design.matrix(edesign, degree = 2)
  rownames(design.matrix$dis) <- rownames(edesign)
  
  fit <- p.vector(expr_round, design.matrix, counts = TRUE, Q = 0.05)
  tfit <- T.fit(fit)
  sigs <- get.siggenes(tfit, rsq = 0.4, vars = "groups")
  
  # Render PDF & Outputs
  if (!is.null(sigs$sig.genes$bPACvsSCI)) {
    pdf(file.path(roi_dir, paste0("gene_profiles_", roi, ".pdf")), width = 24, height = 12)
    cluster <- see.genes(
      sigs$sig.genes$bPACvsSCI,
      show.fit = TRUE,
      dis = design.matrix$dis,
      cluster.method = "hclust",
      cluster.data = 1,
      k = 9
    )
    dev.off()
    
    write.table(cluster$cut,
                file = file.path(roi_dir, paste0("output_cluster_", roi, ".txt")),
                sep = "\t", quote = FALSE)
    write.table(sigs$sig.genes$bPACvsSCI$sig.profiles,
                file = file.path(roi_dir, paste0("output_sig_profiles_bPACvsSCI_", roi, ".txt")),
                sep = "\t", quote = FALSE)
    write.table(sigs$sig.genes$bPACvsSCI$sig.pvalues,
                file = file.path(roi_dir, paste0("output_sig_pvalues_bPACvsSCI_", roi, ".txt")),
                sep = "\t", quote = FALSE)
  }
}

cat("\nmaSigPro execution complete! Outputs saved to:", opt$outdir, "\n")
