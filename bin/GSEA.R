#!/usr/bin/env Rscript

# Activate renv project environment so packages are found
if (file.exists("renv.lock") && requireNamespace("renv", quietly = TRUE)) {
  renv::activate()
}

library(optparse)
library(clusterProfiler)
library(org.Rn.eg.db)
library(dplyr)
library(readr)
library(stringr)
library(tools)
library(ReactomePA)

option_list <- list(
  make_option(c("--deg_dir"), type="character", default=NULL,
              help="Directory containing DEG CSV files", metavar="character"),
  make_option(c("--output_dir"), type="character", default=NULL,
              help="Output directory for ORA results", metavar="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$deg_dir) || is.null(opt$output_dir)) {
  stop("❌ Both --deg_dir and --output_dir must be specified")
}

cat("📂 Using DEG directory:", opt$deg_dir, "\n")
cat("📂 Output directory:", opt$output_dir, "\n")

# List area subfolders
area_dirs <- list.dirs(opt$deg_dir, recursive = FALSE)

if (length(area_dirs) == 0) {
  # No subfolders, fallback to deg_dir itself
  area_dirs <- opt$deg_dir
}


# Find CSV files matching pattern
#csv_files <- list.files(opt$deg_dir, pattern = "\\.csv$", full.names = TRUE)

#if (length(csv_files) == 0) {
#  cat("⚠️  No matching DEG CSV files found in:", opt$deg_dir, "\n")
#  quit(status = 1)
#}

if (!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)

run_gsea_go <- function(ranked_genes, category, output_file) {
  gsea_result <- gseGO(
    geneList = ranked_genes,
    OrgDb = org.Rn.eg.db,
    keyType = "SYMBOL",
    ont = category,
    nPerm = 3000,
    minGSSize = 5,
    maxGSSize = 1000,
    pvalueCutoff = 0.1,
    verbose = TRUE
  )
  
  if (!is.null(gsea_result) && nrow(as.data.frame(gsea_result)) > 0) {
    write.csv(as.data.frame(gsea_result), output_file, row.names = FALSE)
  }
}


# Shared Entrez conversion for KEGG and Reactome
prepare_entrez_ranking <- function(ranked_genes) {
  # Map SYMBOL → ENTREZID and remove duplicates
  entrez_ids <- bitr(names(ranked_genes),
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Rn.eg.db) %>%
                distinct(ENTREZID, .keep_all = TRUE)

  # Merge ranking with Entrez IDs
  merged <- merge(
    data.frame(SYMBOL = names(ranked_genes), logFC = ranked_genes),
    entrez_ids,
    by = "SYMBOL"
  )
  
  # Create named numeric vector: names = Entrez IDs
  ranked_entrez <- setNames(as.numeric(merged$logFC), merged$ENTREZID)
  ranked_entrez <- ranked_entrez[order(ranked_entrez, decreasing = TRUE)]
  
  return(ranked_entrez)
}

run_gsea_kegg <- function(ranked_entrez, output_file) {
  gsea_result <- gseKEGG(
    geneList     = ranked_entrez,
    organism     = "rno",
    nPerm        = 3000,
    minGSSize    = 5,
    maxGSSize    = 1000,
    pvalueCutoff = 0.1,
    verbose      = TRUE
  )

  if (!is.null(gsea_result) && nrow(as.data.frame(gsea_result)) > 0) {
    gsea_result <- setReadable(gsea_result, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
    write.csv(as.data.frame(gsea_result), output_file, row.names = FALSE)

  }
}

run_gsea_reactome <- function(ranked_entrez, output_file) {
  gsea_result <- gsePathway(
    geneList     = ranked_entrez,
    organism     = "rat",
    nPerm        = 3000,
    minGSSize    = 5,
    maxGSSize    = 1000,
    pvalueCutoff = 0.1,
    verbose      = TRUE
  )

  if (!is.null(gsea_result) && nrow(as.data.frame(gsea_result)) > 0) {
    gsea_result <- setReadable(gsea_result, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
    write.csv(as.data.frame(gsea_result), output_file, row.names = FALSE)
  }
}

for (area_dir in area_dirs) {
  area_name <- basename(area_dir)
  cat("📂 Processing area:", area_name, "\n")

  csv_files <- list.files(area_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) {
    cat("⚠️  No DEG CSVs found in", area_dir, "\n")
    next
  }

  # Create area-specific output folder
  area_output_dir <- file.path(opt$output_dir, area_name)
  if (!dir.exists(area_output_dir)) dir.create(area_output_dir, recursive = TRUE)

  # ---- Then your existing csv_file loop goes here ----
  for (csv_file in csv_files) {
    cat("Processing cell type:", csv_file, "\n")
    deg_df <- read_csv(csv_file, show_col_types = FALSE)

    # Normalize column names
    if ("Z" %in% names(deg_df)) {
      deg_df <- deg_df %>% rename(logFC = Z)
    }
    if ("padj" %in% names(deg_df)) {
      deg_df <- deg_df %>% rename(adj.P.Val = padj)
    }

    required_cols <- c("gene", "logFC", "adj.P.Val")
    if (!all(required_cols %in% names(deg_df))) {
      cat("⚠️  Skipping file (missing required columns):", csv_file, "\n")
      next
    }

    # Prepare named vector for GSEA: names = gene symbols, values = logFC
    ranked_genes <- deg_df %>%
      filter(!is.na(gene), !is.na(logFC)) %>%
      distinct(gene, .keep_all = TRUE) %>%
      arrange(desc(logFC)) %>%
      {setNames(.$logFC, .$gene)}

    base_filename <- file_path_sans_ext(basename(csv_file))
  
    # Run GSEA for GO categories
    for (category in c("BP", "CC", "MF")) {
      out_file <- file.path(area_output_dir, paste0(base_filename, "_GSEA_GO_", category, ".csv"))
      run_gsea_go(ranked_genes, category, out_file)
      }

    # Prepare Entrez once for both KEGG & Reactome
    ranked_entrez <- prepare_entrez_ranking(ranked_genes)
    
    # KEGG
    run_gsea_kegg(ranked_entrez, file.path(area_output_dir, paste0(base_filename, "_GSEA_KEGG.csv")))

    # Reactome
    run_gsea_reactome(ranked_entrez, file.path(area_output_dir, paste0(base_filename, "_GSEA_Reactome.csv")))
  }
}

cat("✅ GSEA analysis with GO, KEGG, and Reactome completed. Results saved in:", opt$output_dir)





#for (csv_file in csv_files) {
#  deg_df <- read_csv(csv_file, show_col_types = FALSE)

#  required_cols <- c("gene", "logFC")
#  if (!all(required_cols %in% names(deg_df))) {
#    cat("⚠️  Skipping file (missing required columns):", csv_file, "\n")
#    next
#  }
  
  # Prepare named vector for GSEA: names = gene symbols, values = logFC
#  ranked_genes <- deg_df %>%
#    filter(!is.na(gene), !is.na(logFC)) %>%
#    distinct(gene, .keep_all = TRUE) %>%
#    arrange(desc(logFC)) %>%
#    {setNames(.$logFC, .$gene)}
  
#  base_filename <- file_path_sans_ext(basename(csv_file))
  
  # Run GSEA for GO categories
#  for (category in c("BP", "CC", "MF")) {
#    out_file <- file.path(opt$output_dir, paste0(base_filename, "_GSEA_GO_", category, ".csv"))
#    run_gsea_go(ranked_genes, category, out_file)
#  }

  
  # Prepare Entrez once for both KEGG & Reactome
#  ranked_entrez <- prepare_entrez_ranking(ranked_genes)
  
  # KEGG
#  run_gsea_kegg(ranked_entrez, file.path(opt$output_dir, paste0(base_filename, "_GSEA_KEGG.csv")))
  
  # Reactome
#  run_gsea_reactome(ranked_entrez, file.path(opt$output_dir, paste0(base_filename, "_GSEA_Reactome.csv")))
#}

#cat("✅ GSEA analysis with GO, KEGG, and Reactome completed. Results saved in:", opt$output_dir, "\n")
