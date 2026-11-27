#!/usr/bin/env Rscript
# bin/masigpro_expression.R
# Time-course gene expression per area with maSigPro (paper: "Differential gene
# expression profile analysis").
#   1. sum spot counts per area and sample (pseudobulk)
#   2. filterByExpr (condition x time groups) + TMM on all retained genes
#   3. optionally restrict to gene lists, keeping the full-transcriptome
#      library sizes / norm factors
#   4. maSigPro on rounded TMM-CPM: negative-binomial GLM, degree 2,
#      p.vector Q (BH), T.fit backward (alfa), get.siggenes R2 > rsq, vars = "groups"
#   5. cluster each contrast's significant genes into temporal profiles
# Design (series, shared baseline time zero, reference series) is auto-detected
# from sample_metadata (sample, condition, time, group); see build_design().

suppressPackageStartupMessages({
  library(optparse)
  library(edgeR)
})

.this_dir <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else "."
})
source(file.path(.this_dir, "masigpro_utils.R"))

option_list <- list(
  make_option("--counts_dir", type = "character",
              help = "Directory with <sample>_counts.csv (genes x spots)"),
  make_option("--delineation_dir", type = "character", default = "NO_DELINEATION"),
  make_option("--sample_metadata", type = "character",
              help = "CSV: sample, condition, time, group"),
  make_option("--baseline_group", type = "character", default = "",
              help = "Group used as the shared time zero of every series (e.g. SHAM)"),
  make_option("--group_order", type = "character", default = "",
              help = "Comma-separated series order; the first is the maSigPro reference"),
  make_option("--areas", type = "character", default = "ALL_SPOTS",
              help = "Comma-separated delineation labels, or ALL_SPOTS"),
  make_option("--gene_list", type = "character", default = "NO_FILE",
              help = "Comma-separated gene list files (one gene per line); union is tested"),
  make_option("--degree", type = "integer", default = 2),
  make_option("--Q", type = "double", default = 0.05),
  make_option("--alfa", type = "double", default = 0.05),
  make_option("--rsq", type = "double", default = 0.4),
  make_option("--k", type = "integer", default = 9, help = "Clusters per contrast"),
  make_option("--outdir", type = "character", default = "masigpro_expression")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

meta  <- read_meta(opt$sample_metadata)
areas <- split_arg(opt$areas)
if (length(areas) && opt$delineation_dir == "NO_DELINEATION") {
  warning("areas given but no delineation_dir: using whole-sample pseudobulk")
  areas <- character(0)
}
group_order <- split_arg(opt$group_order)

# Whole-cohort check first, so a 2-timepoint config exits cleanly.
d_all <- build_design(meta, opt$baseline_group, group_order, opt$degree)
if (!is.null(d_all$skip)) {
  msg <- paste("maSigPro expression skipped:", d_all$skip)
  cat(msg, "\n"); writeLines(msg, file.path(opt$outdir, "SKIPPED.txt"))
  quit(status = 0)
}

genes_of_interest <- NULL
gene_files <- split_arg(opt$gene_list)
if (length(gene_files)) {
  genes_of_interest <- unique(trimws(unlist(lapply(gene_files, readLines))))
  genes_of_interest <- genes_of_interest[nzchar(genes_of_interest)]
  cat("Loaded", length(genes_of_interest), "genes of interest from", length(gene_files), "file(s)\n")
}

read_counts <- function(f) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::fread(f, data.table = FALSE, check.names = FALSE)
  } else {
    dt <- read.csv(f, check.names = FALSE)
  }
  genes <- as.character(dt[[1]])
  dt <- dt[, -1, drop = FALSE]
  list(genes = genes, df = dt)
}

# -----------------------------------------------------------------------------
# 1. Pseudobulk per area x sample (one pass over each counts file)
# -----------------------------------------------------------------------------
area_keys <- if (length(areas)) areas else "ALL_SPOTS"
pb <- setNames(lapply(area_keys, function(a) list()), area_keys)
n_spots <- list()

for (s in meta$sample) {
  f <- file.path(opt$counts_dir, paste0(s, "_counts.csv"))
  if (!file.exists(f)) { warning("missing counts for ", s); next }
  cat("Reading", s, "\n")
  cc <- read_counts(f)
  sp <- spots_by_area(areas, opt$delineation_dir, s, colnames(cc$df))
  for (a in names(sp)) {
    bc <- sp[[a]]
    n_spots[[length(n_spots) + 1]] <- data.frame(area = a, sample = s, n_spots = length(bc))
    if (length(bc) == 0) { warning("no spots in ", a, " for ", s); next }
    v <- rowSums(as.matrix(cc$df[, bc, drop = FALSE]))
    names(v) <- cc$genes
    if (anyDuplicated(names(v))) v <- tapply(v, names(v), sum)
    pb[[a]][[s]] <- v
  }
  rm(cc); invisible(gc())
}
if (length(n_spots))
  write.table(do.call(rbind, n_spots), file.path(opt$outdir, "spots_per_area_sample.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------------------------------------------------------
# 2-5. Per area
# -----------------------------------------------------------------------------
summary_rows <- list()
for (a in area_keys) {
  cat("\n=== Area:", a, "===\n")
  adir <- file.path(opt$outdir, a)
  dir.create(adir, recursive = TRUE, showWarnings = FALSE)
  if (length(pb[[a]]) < 2) { warning("fewer than 2 samples with spots in ", a); next }

  genes <- Reduce(union, lapply(pb[[a]], names))
  counts <- sapply(pb[[a]], function(v) { x <- v[genes]; x[is.na(x)] <- 0; x })
  rownames(counts) <- genes
  m <- meta[match(colnames(counts), meta$sample), ]

  d <- build_design(m, opt$baseline_group, group_order, opt$degree)
  if (!is.null(d$skip)) {
    msg <- paste("skipped:", d$skip); cat(" ", msg, "\n")
    writeLines(msg, file.path(adir, "SKIPPED.txt")); next
  }
  describe_design(d, m)
  write.table(cbind(sample = rownames(d$edesign), d$edesign), file.path(adir, "edesign.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  dge <- DGEList(counts = counts)
  grp <- factor(paste(m$condition, m$time, sep = "_"))
  dge <- dge[filterByExpr(dge, group = grp), , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  cpm_all <- cpm(dge, log = FALSE)          # lib sizes + TMM from the full filtered set
  tested <- if (is.null(genes_of_interest)) rownames(cpm_all) else intersect(rownames(cpm_all), genes_of_interest)
  cat("  genes after filterByExpr:", nrow(dge), "| tested:", length(tested), "\n")

  write.table(cbind(gene = rownames(counts), counts), file.path(adir, "pseudobulk_counts.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(cbind(gene = rownames(cpm_all), round(cpm_all, 3)), file.path(adir, "tmm_cpm.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(tested, file.path(adir, "tested_genes.txt"))
  if (length(tested) < 2) { warning("fewer than 2 genes to test in ", a); next }

  expr <- round(cpm_all[tested, , drop = FALSE])
  res <- run_masigpro(expr, d, counts = TRUE, Q = opt$Q, alfa = opt$alfa, rsq = opt$rsq)
  sig <- write_masigpro_tables(res, adir, item = "gene")

  # temporal clusters per contrast
  if (!is.null(res$sigs)) for (cn in names(res$sigs$sig.genes)) {
    ct <- res$sigs$sig.genes[[cn]]
    n <- n_sig(ct)
    if (n < 3) next
    k <- min(opt$k, n - 1)
    pdf(file.path(adir, paste0("gene_profiles_", cn, ".pdf")), width = 24, height = 12)
    cl <- tryCatch(see.genes(ct, show.fit = TRUE, dis = d$design$dis, edesign = d$edesign,
                             cluster.method = "hclust", cluster.data = 1, k = k),
                   error = function(e) { warning("see.genes ", a, " ", cn, ": ", conditionMessage(e)); NULL })
    dev.off()
    if (!is.null(cl))
      write.table(data.frame(gene = names(cl$cut), cluster = cl$cut),
                  file.path(adir, paste0("clusters_", cn, ".tsv")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
  }

  n_by <- if (nrow(sig)) table(sig$contrast) else integer(0)
  summary_rows[[a]] <- data.frame(
    area = a, n_samples = ncol(counts), n_filtered = nrow(dge), n_tested = length(tested),
    n_pvector = if (is.null(res$fit$i)) 0L else res$fit$i,
    n_sig_any = length(unique(sig$feature)),
    per_contrast = if (length(n_by)) paste(names(n_by), as.integer(n_by), sep = "=", collapse = ";") else "")
}

if (length(summary_rows)) {
  s <- do.call(rbind, summary_rows)
  write.table(s, file.path(opt$outdir, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n"); print(s)
}
