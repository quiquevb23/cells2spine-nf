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

convert_to_entrez <- function(genes) {
  suppressMessages({
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)
  }) %>%
    distinct(ENTREZID, .keep_all = TRUE) %>%
    pull(ENTREZID) %>%
    unique()
}

run_enrichment <- function(gene_list, category = NULL, enrichment_type, output_file) {
  if (length(gene_list) == 0) return(NULL)
  
  enrich_result <- NULL
  
  if (enrichment_type == "GO") {
    enrich_result <- enrichGO(
      gene          = gene_list,
      OrgDb         = org.Rn.eg.db,
      keyType       = "SYMBOL",
      ont           = category,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05
    )
    
  } else if (enrichment_type == "KEGG") {
    enrich_result <- enrichKEGG(
      gene          = gene_list,     # already Entrez IDs
      organism      = "rno",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05
    )
    
  } else if (enrichment_type == "Reactome") {
    enrich_result <- enrichPathway(
      gene          = gene_list,     # already Entrez IDs
      organism      = "rat",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      readable      = TRUE
    )
  }
  
  if (!is.null(enrich_result) && nrow(as.data.frame(enrich_result)) > 0) {
    write.csv(as.data.frame(enrich_result), output_file, row.names = FALSE)
  }
}

# ===== Main loop =====

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
    deg_df <- read_csv(csv_file, show_col_types = FALSE)


    # Normalize column names (to bypass the requirement)
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

    sig_genes <- deg_df %>%
      filter(!is.na(gene), !is.na(adj.P.Val), adj.P.Val < 0.05)

    up_genes   <- sig_genes %>% filter(logFC > 0)  %>% pull(gene) %>% unique()
    down_genes <- sig_genes %>% filter(logFC < 0)  %>% pull(gene) %>% unique()

    up_entrez   <- convert_to_entrez(up_genes)
    down_entrez <- convert_to_entrez(down_genes)
  
    base_filename <- file_path_sans_ext(basename(csv_file))

    # GO: BP, CC, MF
    for (category in c("BP", "CC", "MF")) {
      run_enrichment(up_genes,   category, "GO", file.path(area_output_dir, paste0(base_filename, "_ORA_GO_UP_", category, ".csv")))
      run_enrichment(down_genes, category, "GO", file.path(area_output_dir, paste0(base_filename, "_ORA_GO_DOWN_", category, ".csv")))
    }

    # KEGG
    run_enrichment(up_entrez,   NULL, "KEGG", file.path(area_output_dir, paste0(base_filename, "_ORA_KEGG_UP.csv")))
    run_enrichment(down_entrez, NULL, "KEGG", file.path(area_output_dir, paste0(base_filename, "_ORA_KEGG_DOWN.csv")))

    # Reactome
    run_enrichment(up_entrez,   NULL, "Reactome", file.path(area_output_dir, paste0(base_filename, "_ORA_Reactome_UP.csv")))
    run_enrichment(down_entrez, NULL, "Reactome", file.path(area_output_dir, paste0(base_filename, "_ORA_Reactome_DOWN.csv")))
  }
}

#csv_files <- list.files(opt$deg_dir, pattern = "\\.csv$", full.names = TRUE)

#if (length(csv_files) == 0) {
#  cat("⚠️  No matching ORA CSV files found in:", opt$deg_dir, "\n")
#}

#if (!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)

#for (csv_file in csv_files) {
#  deg_df <- read_csv(csv_file, show_col_types = FALSE)
  
#  required_cols <- c("gene", "logFC", "FDR")
#  if (!all(required_cols %in% names(deg_df))) {
#    cat("⚠️  Skipping file (missing required columns):", csv_file, "\n")
#    next
#  }
  
#  sig_genes <- deg_df %>%
#    filter(!is.na(gene), !is.na(FDR), FDR < 0.05)
  
#  up_genes   <- sig_genes %>% filter(logFC > 0)  %>% pull(gene) %>% unique()
#  down_genes <- sig_genes %>% filter(logFC < 0)  %>% pull(gene) %>% unique()
  
#  up_entrez   <- convert_to_entrez(up_genes)
#  down_entrez <- convert_to_entrez(down_genes)
  
#  base_filename <- file_path_sans_ext(basename(csv_file))
  
  # GO: BP, CC, MF
#  for (category in c("BP", "CC", "MF")) {
#    run_enrichment(up_genes,   category, "GO", file.path(opt$output_dir, paste0(base_filename, "_ORA_GO_UP_", category, ".csv")))
#    run_enrichment(down_genes, category, "GO", file.path(opt$output_dir, paste0(base_filename, "_ORA_GO_DOWN_", category, ".csv")))
#  }
  
#  # KEGG
#  run_enrichment(up_entrez,   NULL, "KEGG", file.path(opt$output_dir, paste0(base_filename, "_ORA_KEGG_UP.csv")))
#  run_enrichment(down_entrez, NULL, "KEGG", file.path(opt$output_dir, paste0(base_filename, "_ORA_KEGG_DOWN.csv")))
  
  # Reactome
#  run_enrichment(up_entrez,   NULL, "Reactome", file.path(opt$output_dir, paste0(base_filename, "_ORA_Reactome_UP.csv")))
#  run_enrichment(down_entrez, NULL, "Reactome", file.path(opt$output_dir, paste0(base_filename, "_ORA_Reactome_DOWN.csv")))
#}

#cat("✅ ORA analysis with GO, KEGG, and Reactome completed. Results saved in:", opt$output_dir, "\n")


