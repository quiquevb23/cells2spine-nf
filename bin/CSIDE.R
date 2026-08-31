#!/usr/bin/env Rscript

# ---------------------------
# Parse arguments
# ---------------------------
# Load the .h5ad file from cell2location
args <- commandArgs(trailingOnly = TRUE)

ref_level <- args[1]
spatial_input <- args[2]     # Here we have the .rds dir for each sample after RCTD
output_base_dir <- args[3]   # Main base directory for outputs like references
ref_label <- args[4]
masked_celltypes <- args[5]
renv_project <- args[6]

# position 7: comma-separated conditions string
conditions <- strsplit(args[7], ",")[[1]]

# positions 8 onward: samples vector
#samples <- args[11:length(args)]

# everything after that are samples, "::", and conditions_map
rest <- args[8:length(args)]
split_idx <- which(rest == "::")

if(length(split_idx) == 0) stop("[Rscript] '::' sentinel missing for conditions_map")
samples <- rest[1:(split_idx-1)]
conditions_map <- rest[(split_idx+1):length(rest)]

# 2️⃣ check if "::conditions" sentinel is present
cond_sentinel_idx <- which(conditions_map == "::conditions")
if(length(cond_sentinel_idx) != 1) stop("[Rscript] '::conditions' sentinel missing for test/ref")
# extract test/ref
test_cond <- conditions_map[cond_sentinel_idx + 1]
ref_cond  <- conditions_map[cond_sentinel_idx + 2]

# remove sentinel + test/ref from conditions_map
conditions_map <- conditions_map[1:(cond_sentinel_idx-1)]

# Turn "Spatial_1:healthy" into a named vector
condition_lookup <- setNames(
  sub(".*:", "", conditions_map),  # condition
  sub(":.*", "", conditions_map)   # sample
)

# ---------------------------
# Debug prints
# ---------------------------
cat("[DEBUG] samples:\n")
print(samples)
cat("[DEBUG] condition_lookup:\n")
print(condition_lookup)
cat("[DEBUG] test_cond =", test_cond, " | ref_cond =", ref_cond, "\n")

#########LOOKS LIKE THIS#########
#condition_lookup <- c(
#  "Spatial_1" = "healthy",
#  "Spatial_2" = "injured10",
#  "Spatial_3" = "injured30",
#  "Spatial_4" = "treated10"
#)

#max_cores <- as.integer(args[4])

# Define CSIDE_OUTPUT_DIR: RUN_NAME + CSIDE
cside_output_dir <- file.path(output_base_dir, "CSIDE")

# create the directory if it does not exist
if (!dir.exists(cside_output_dir)) {
  dir.create(cside_output_dir, recursive = TRUE)
}

# ---------------------------
# Load renv environment
# ---------------------------

library(renv)
renv::restore(prompt = FALSE)


# ---------------------------
# Now load project libraries
# ---------------------------
suppressPackageStartupMessages({
  library(spacexr)
  library(Matrix)
  library(limma)
})

# 3. Apply function patches
#source("./cside_fixes.R")
# This function is a patch that solves bugs in CSIDE.population.inference

# ---------------------------
# Helper functions
# ---------------------------

# Dir for input data (spatial samples and references)
data_dir = file.path(output_base_dir, "data")
cell2loc_main_dir = file.path(output_base_dir, "cell2location_map")

# ---- Gene filtering (analogous to cell2location.filter_genes) ----
filter_genes <- function(counts,
                         cell_count_cutoff = 5,
                         cell_percentage_cutoff2 = 0.03,
                         nonz_mean_cutoff = 1.12) {

    # Ensure dgCMatrix for fast row/col operations
    if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")

    n_cells <- ncol(counts)
    n_genes <- nrow(counts)

    # Number of cells with non-zero counts per gene (fast for sparse)
    cell_counts <- Matrix::rowSums(counts > 0)

    # Fraction of cells expressing the gene
    cell_fractions <- cell_counts / n_cells

    # Equivalent percentage cutoff expressed as an absolute cell count
    percent_cutoff_cells <- ceiling(cell_percentage_cutoff2 * n_cells)

    # Compute mean of non-zero counts per gene efficiently:
    # total counts per gene divided by number of non-zero cells (guarding zeros)
    total_counts_per_gene <- Matrix::rowSums(counts)  # sum across cells
    nonz_means <- numeric(n_genes)
    nonz_idx <- which(cell_counts > 0)
    if (length(nonz_idx) > 0) {
        nonz_means[nonz_idx] <- total_counts_per_gene[nonz_idx] / cell_counts[nonz_idx]
    }

    # Apply logic:
    # - Exclude genes with cell_counts < cell_count_cutoff
    # - Include genes with cell_counts >= percent_cutoff_cells (i.e. >= percentage)
    # - For genes with cell_count_cutoff <= cell_counts < percent_cutoff_cells,
    #   include only if nonz_means >= nonz_mean_cutoff

    keep <- logical(n_genes)

    # Genes that are expressed in at least the percentage cutoff -> always keep
    if (percent_cutoff_cells <= 1) {
        # if percent cutoff resolves to 0 or 1 cell, include those logically by fraction check
        keep[cell_fractions >= cell_percentage_cutoff2] <- TRUE
    } else {
        keep[cell_counts >= percent_cutoff_cells] <- TRUE
    }

    # Genes below absolute cutoff -> always excluded (already FALSE)
    # Genes in the intermediate range:
    intermediate_idx <- which((cell_counts >= cell_count_cutoff) &
                              (cell_counts < percent_cutoff_cells))
    if (length(intermediate_idx) > 0) {
        keep[intermediate_idx] <- nonz_means[intermediate_idx] >= nonz_mean_cutoff
    }

    # Informative message
    cat(sprintf("Filtering genes: keeping %d of %d genes (abs_cutoff=%d, pct_cutoff=%d cells, nonz_mean_cutoff=%.3g)\n",
                sum(keep), n_genes, cell_count_cutoff, percent_cutoff_cells, nonz_mean_cutoff))

    return(keep)  # logical vector
}

# Function to load reference data
load_reference_data <- function(label, data_dir) {
    label_dir <- file.path(data_dir, paste0("reference_", label))

    counts_file <- file.path(label_dir, paste0("counts.mtx"))
    cell_types_file <- file.path(label_dir, paste0("cell_types.csv"))
    genes_file <- file.path(label_dir, paste0("genes.csv"))
    barcodes_file <-  file.path(label_dir, paste0("cells.csv"))
    
    # Load row (genes) and column (cell) names, skipping the first line
    genes <- read.csv(genes_file, header=FALSE, stringsAsFactors=FALSE)[,1]
    barcodes <- read.csv(barcodes_file, header=FALSE, stringsAsFactors=FALSE)[,1]

    # Load counts matrix 
    counts <- readMM(counts_file)
    counts <- as(counts, "CsparseMatrix")  # Convert to sparse format
    
    # Ensure correct dimensions before assigning names
    if (length(genes) != nrow(counts)) {
        stop("Error: Number of genes does not match row count in counts matrix!")
    }
    if (length(barcodes) != ncol(counts)) {
        stop("Error: Number of barcodes does not match column count in counts matrix!")
    }
    
    # Handle duplicate gene names
    if (any(duplicated(genes))) {
        warning("Duplicate gene names detected! Making names unique.")
        genes <- make.unique(genes)  # Appends ".1", ".2", etc. to duplicates
    }

    # Handle duplicate gene names
    if (any(duplicated(barcodes))) {
        warning("Duplicate gene names detected! Making names unique.")
        barcodes <- make.unique(barcodes)  # Appends ".1", ".2", etc. to duplicates
    }

    # Assign row and column names
    rownames(counts) <- genes
    colnames(counts) <- barcodes

    # Load cell types metadata (No header → Manually define column names)
    cell_types_df <- read.csv(cell_types_file, header=FALSE, stringsAsFactors=FALSE)
    colnames(cell_types_df) <- c("barcodes", "cell_type")  # Manually assign names

    # Ensure unique barcode names in metadata
    if (any(duplicated(cell_types_df$barcodes))) {
        warning("Duplicate barcodes in metadata detected! Making them unique.")
        cell_types_df$barcodes <- make.unique(cell_types_df$barcodes)
    }

    # Ensure barcodes match the counts matrix
    cell_types_df <- cell_types_df[cell_types_df$barcodes %in% colnames(counts), ]

    # ❗ Apply masking here
    if (masked_celltypes != "" && masked_celltypes != "NA") {
        to_mask <- unlist(strsplit(masked_celltypes, ","))
        cat("Masking cell types:", paste(to_mask, collapse=", "), "\n")
        cell_types_df <- cell_types_df[!cell_types_df$cell_type %in% to_mask, ]
    }

    cat("Number of genes before filtering:", length(rownames(counts)), "\n")
    # ---- NEW: Apply gene filtering like cell2location for the reference ----
    #keep_genes <- filter_genes(counts)

    # do not filter genes as cell2location does
    keep_genes <- rep(TRUE, nrow(counts))

    cat("Number of genes kept after filtering:", sum(keep_genes), "\n")

    counts <- counts[keep_genes, , drop = FALSE]

    cell_types <- factor(cell_types_df$cell_type)
    names(cell_types) <- cell_types_df$barcodes

    cat("Cell Type Counts:\n")
    cell_type_counts <- table(cell_types)
    print(cell_type_counts)

    # Remove cell types with fewer than 25 barcodes (let's skip this step as we did not do it in cell2loc)
    valid_cell_types <- names(cell_type_counts[cell_type_counts >= 0])
    removed_cell_types <- setdiff(names(cell_type_counts), valid_cell_types)
    if (length(removed_cell_types) > 0) {
        cat("Removing cell types with fewer than 25 cells:", paste(removed_cell_types, collapse = ", "), "\n")
    }

    cell_types_df <- cell_types_df[cell_types_df$cell_type %in% valid_cell_types, ]

    # Update counts matrix and cell types
    counts <- counts[, cell_types_df$barcodes, drop=FALSE]
    cell_types <- factor(cell_types_df$cell_type)
    names(cell_types) <- cell_types_df$barcodes
    
    cat("Cell Type Counts after filtering:\n")
    print(table(cell_types))

    nUMI <- colSums(counts) # Sum over rows (genes) to get per-cell nUMI

    return(Reference(counts, cell_types, nUMI))
}

#HELPER: subsetting region of RCTD object

subset_RCTD_replicates <- function(RCTD_reps, barcodes) {
  new_reps <- RCTD_reps
  for (i in seq_along(new_reps@RCTD.reps)) {
    rep_obj <- new_reps@RCTD.reps[[i]]
    keep <- colnames(rep_obj@spatialRNA@counts) %in% barcodes
    
    # Subset spatial RNA
    rep_obj@spatialRNA@counts <- rep_obj@spatialRNA@counts[, keep, drop = FALSE]
    rep_obj@spatialRNA@nUMI   <- rep_obj@spatialRNA@nUMI[keep]
    rep_obj@spatialRNA@coords <- rep_obj@spatialRNA@coords[keep, , drop = FALSE]
    
    new_reps@RCTD.reps[[i]] <- rep_obj
  }
  return(new_reps)
}

cat("Loading sc references")

# Load reference data for both labels (dynamically)
references <- lapply(conditions, function(cond) {
  load_reference_data(cond, data_dir)
})
names(references) <- conditions

# Print summaries
for (cond in conditions) {
  cat(sprintf("%s reference cell types:\n", cond))
  print(table(references[[cond]]@cell_types))
}

#reference_healthy <- load_reference_data("healthy", data_dir)
#reference_injured <- load_reference_data("injured", data_dir)

#cat("Healthy reference cell types:\n")
#print(table(reference_healthy@cell_types))

#cat("Injured reference cell types:\n")
#print(table(reference_injured@cell_types))



# Function to load spatial data (counts are from cell2loc) 
load_spatial_data <- function(sample_name, sample_dir, cell2loc_dir) {
  counts_file <- file.path(cell2loc_dir, paste0(sample_name, "_counts.csv"))
  coords_file <- file.path(sample_dir, paste0(sample_name, "_coordinates.csv"))
  
  counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
  
  # Load coordinates (barcodes are in the first column)
  coords <- read.csv(coords_file, row.names = 1, check.names = FALSE)
  
  nUMI <- colSums(counts)
  
  return(SpatialRNA(coords, counts, nUMI))
}


###COMBINE REFERENCES

# Combine counts matrices
#combined_counts <- cbind(reference_healthy@counts, reference_injured@counts)

# Combine cell type labels
#combined_cell_types <- factor(c(as.character(reference_healthy@cell_types),
#                                as.character(reference_injured@cell_types)))
#names(combined_cell_types) <- c(names(reference_healthy@cell_types),
#                                names(reference_injured@cell_types))

# Combine nUMI
#combined_nUMI <- c(reference_healthy@nUMI, reference_injured@nUMI)

# Create combined reference
#reference_combined <- Reference(combined_counts, combined_cell_types, combined_nUMI)


### LOAD SPATIAL REPS

# Spatial samples are in data_dir + {sample}

pucks <- lapply(samples, function(sample) {
  sample_dir <- file.path(data_dir, sample)
  cell2loc_dir <- file.path(cell2loc_main_dir, sample)
  load_spatial_data(sample, sample_dir, cell2loc_dir)
})
names(pucks) <- samples



# Obtain region annotations which are in data_dir + {sample} (make a dataframe)
region_annotations <- lapply(samples, function(sample) {
  sample_dir <- file.path(data_dir, sample)
  file <- file.path(sample_dir, paste0(sample, "_manual_delineation.csv"))
  reg <- read.csv(file, header = TRUE, check.names = FALSE)

  rownames(reg) <- reg[, 1]
  reg <- reg[, -1, drop = FALSE]
  return(reg)
})
names(region_annotations) <- samples


#----------------OLD WAY ---------------
# List only the CSV files matching "Spatial_<number>_counts.csv"
#sample_files <- list.files(tmp_dir, pattern = "^Spatial_\\d+_counts\\.csv$")
# Extract numbers and order by them
#sample_nums <- as.integer(sub("Spatial_(\\d+)_counts\\.csv", "\\1", sample_files))
#sample_files <- sample_files[order(sample_nums)]

#sample_names <- sub("_counts\\.csv$", "", sample_files)

#pucks <- lapply(sample_names, function(sample) load_spatial_data(sample, tmp_dir))
#names(pucks) <- sample_names

# Example: first 2 replicates healthy (0), last 2 injured (1)
#replicate_names <- sample_names
#group_ids <- c(0,0,1,1)

#exvar_list <- lapply(sample_names, function(sample) {
#  if (sample %in% c("Spatial_1", "Spatial_2")) {
#    val <- 0  # healthy
#  } else {
#    val <- 1  # injured
#  }
#  barcodes <- colnames(pucks[[sample]]@counts)   # get barcodes
  # Create a simple region variable (for example)
#  exvar <- rep(val, length(barcodes))
#  names(exvar) <- barcodes
#  return(exvar)
#})
#names(exvar_list) <- sample_names  # make sure names match replicate_names

#print(names(exvar_list))

#group_ids <- rep(1, length(pucks))  # same group for CSIDE, or adjust if needed

# ---------------------------
# Build RCTD replicates
# ---------------------------
#myRCTD.reps <- create.RCTD.replicates(
#  pucks,
#  reference_combined,
#  replicate_names,
#  group_ids = group_ids,
#  max_cores = 4
#)

#myRCTD.reps <- run.RCTD.replicates(myRCTD.reps)
#saveRDS(myRCTD.reps, file.path(cside_output_dir, "myRCTD_reps_afterRCTD.rds"))

#cat("RCTD replicates created:\n")
#print(names(myRCTD.reps@r))   # safe check

#cell_types <- levels(myRCTD.reps@RCTD.reps[[1]]@reference@cell_types) #this needs to be assigned for CSIDE

## LOAD previous saved rds ##
#myRCTD.reps <- readRDS(file.path(cside_output_dir,'myRCTDde_reps_afterRCTD.rds'))

# ---------------------------
# Run RCTD independently per sample with condition-specific reference
# ---------------------------

#RCTD_list <- list()
#for (i in seq_along(pucks)) {
#  sample_name <- replicate_names[i]
#  condition <- group_ids[i]
  
  # choose reference
#  ref <- if (condition == 0) reference_healthy else reference_injured
  
  # build and run RCTD object
#  rctd <- create.RCTD(
#    pucks[[i]],
#    ref,
#    max_cores = 4,
#  )
  
#  rctd <- run.RCTD(rctd, doublet_mode = "full")
#  RCTD_list[[sample_name]] <- rctd
#}
#----------------------------OLD------------



# ---------------------------GOOD WAY TO DO RCTD
# Run RCTD independently per sample with condition-specific reference
# Using condition_lookup
# ---------------------------

RCTD_list <- list()

for (sample_name in names(pucks)) {
  condition <- condition_lookup[[sample_name]]  # find condition for sample
  ref <- references[[condition]]                # pick reference dynamically

  rctd <- create.RCTD(
    pucks[[sample_name]],
    ref,
    CELL_MIN_INSTANCE = 1,
    max_cores = 4
  )

  rctd <- run.RCTD(rctd, doublet_mode = "full")
  RCTD_list[[sample_name]] <- rctd
}

replicate_names <- vapply(names(RCTD_list),
                          function(s) condition_lookup[[s]],
                          character(1))

# ---------------------------
# Combine into a replicate object
# ---------------------------
myRCTD.reps <- merge_RCTD_objects(
  RCTD_list,
  replicate_names = replicate_names,
)

saveRDS(myRCTD.reps, file.path(cside_output_dir, "myRCTD_reps_afterRCTD.rds"))

# ---------------------------
# Load RCTD object
# ---------------------------

myRCTD.reps <- readRDS(file.path(cside_output_dir, "myRCTD_reps_afterRCTD.rds"))



#for (i in seq_along(myRCTD.reps@RCTD.reps)) {
#  myRCTD.reps@RCTD.reps[[i]]@config[["doublet_mode"]] <- "full"
  
#  myRCTD.reps@RCTD.reps[[i]] <- run.CSIDE.intercept(
#    myRCTD.reps@RCTD.reps[[i]],
#    cell_type_threshold = 0,
#    doublet_mode = FALSE,
#    weight_threshold = 0.1   # or higher
#  )
#}
#saveRDS(myRCTD.reps, file.path(cside_output_dir, "myRCTD_reps_afterIntercept_ct_thr_0.rds"))


#myRCTD.reps <- readRDS(file.path(cside_output_dir, "myRCTD_reps_afterIntercept_ct_thr_0.rds"))


# ---- after loading myRCTD.reps from RDS ----

# 1) Inspect what you just loaded (quick check)
cat("Inspecting loaded myRCTD.reps:\n")
print(class(myRCTD.reps))
cat("Replicate count:", length(myRCTD.reps@RCTD.reps), "\n")
cat("Available slots in a replicate (example):\n")
print(slotNames(myRCTD.reps@RCTD.reps[[1]]))

# 2) sample names (derive from myRCTD.reps) and ensure they are set
sample_names <- names(myRCTD.reps@RCTD.reps)
if (is.null(sample_names) || any(sample_names == "")) {
  sample_names <- paste0("Sample_", seq_along(myRCTD.reps@RCTD.reps))
  names(myRCTD.reps@RCTD.reps) <- sample_names
}
samples <- sample_names   # keep your previous variable name for compatibility
cat("Samples recovered:", paste(samples, collapse = ", "), "\n")

# 3) RCTD_list alias (so old code referencing RCTD_list keeps working)
RCTD_list <- myRCTD.reps@RCTD.reps

# 5) Ensure condition_lookup exists; if not, recreate from conditions_map (CLI parsing earlier)
if (!exists("condition_lookup")) {
  if (exists("conditions_map")) {
    condition_lookup <- setNames(
      sub(".*:", "", conditions_map),
      sub(":.*", "", conditions_map)
    )
  } else {
    stop("condition_lookup missing and cannot reconstruct: ensure 'conditions_map' arg was passed.")
  }
}
# Quick check
missing_map <- samples[is.na(vapply(samples, function(s) condition_lookup[[s]], character(1)))]
if (length(missing_map) > 0) stop("No condition mapping for samples: ", paste(missing_map, collapse = ", "))

##############
# 6) Build meta.design.matrix (one textual 'condition' column and binary columns per condition)
#meta.design.matrix <- data.frame(
#  condition = vapply(samples, function(s) condition_lookup[[s]], character(1)),
#  stringsAsFactors = FALSE
#)
#rownames(meta.design.matrix) <- samples

# Explicitly set 'healthy' as the reference
#meta.design.matrix$condition <- factor(meta.design.matrix$condition, levels = c("healthy", setdiff(conditions, "healthy")))

# Assuming "healthy" is the label for samples 1 & 2 and "injured" for 3 & 4
#meta.design.matrix$condition <- factor(meta.design.matrix$condition,
#                                       levels = c("healthy", "injured"))
#cat("Condition factor levels (reference first):", levels(meta.design.matrix$condition), "\n")
###############

# Get all unique conditions
all_conditions <- unique(vapply(samples, function(s) condition_lookup[[s]], character(1)))

# Optionally, set a reference (e.g., "healthy") if it exists
ref_cond <- if ("healthy" %in% all_conditions) "healthy" else all_conditions[1]

meta.design.matrix <- data.frame(
  condition = factor(vapply(samples, function(s) condition_lookup[[s]], character(1)),
                     levels = c(ref_cond, setdiff(all_conditions, ref_cond))),
  stringsAsFactors = FALSE
)
rownames(meta.design.matrix) <- samples

cat("Condition factor levels (reference first):", levels(meta.design.matrix$condition), "\n")
##############

for (cond in conditions) {
  meta.design.matrix[[cond]] <- as.numeric(meta.design.matrix$condition == cond)
}
cat("meta.design.matrix built. Rownames:", paste(rownames(meta.design.matrix), collapse = ", "), "\n")

# 7) Sanity checks on de_results slots that downstream code expects
cat("Checking presence of de_results slots in replicates:\n")
check_slots <- sapply(myRCTD.reps@RCTD.reps, function(r) {
  dr <- r@de_results
  c(has_gene_fits = !is.null(dr$gene_fits),
    has_results_weights = !is.null(r@results$weights))
})
print(check_slots)

# 8) Safe writing helper (converts sparse matrices to dense before writing)
safe_write_matrix_csv <- function(mat, path, ...) {
  if (inherits(mat, "dgCMatrix") || inherits(mat, "CsparseMatrix") || inherits(mat, "TsparseMatrix")) {
    mat <- as.matrix(mat)
  }
  tryCatch({
    write.csv(mat, path, quote = FALSE, ...)
  }, error = function(e) {
    warning("Failed to write ", path, ": ", e$message)
  })
}



# --------------------------- review this part
# Save just weights 
# ---------------------------
for (i in seq_along(myRCTD.reps@RCTD.reps)) {
  rep <- myRCTD.reps@RCTD.reps[[i]]
  sample_name <- names(RCTD_list)[i]

  outdir <- file.path(output_base_dir, "RCTD", sample_name)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  # weights
  if (!is.null(rep@results$weights)) {
    safe_write_matrix_csv(rep@results$weights, file.path(outdir, paste0(sample_name, "_RCTD_weights.csv")))
  } else warning("No weights found for sample ", sample_name)
}



# After running CSIDE with only intercept, we get specific gene expression values for
# each cell type, we can run region-specific (subset data) Wilcox test to get cell type
# specific gene differences: Better do Z-test computable from CSIDE output

library(dplyr)

#-----------------------------
# Function 1: Average across replicates
#-----------------------------

run_region_dea_avg <- function(myRCTD.reps, regions, region_annotations, design_matrix, outdir) {
  results_list <- list()

  for (reg in regions) {
    regdir <- file.path(outdir, reg)
    dir.create(regdir, showWarnings = FALSE, recursive = TRUE)

    cat("Running DEA for region:", reg, "\n")

    # Collect results from each replicate
    reg_results <- list()
    for (i in seq_along(myRCTD.reps@RCTD.reps)) {
      sample_name <- names(region_annotations)[i]
      region_barcodes <- names(region_annotations[[sample_name]])[region_annotations[[sample_name]] == reg]

      subset_reps <- subset_RCTD_replicates(myRCTD.reps, region_barcodes)
      rep <- subset_reps@RCTD.reps[[i]]

      mean_mat <- rep@de_results$gene_fits$mean_val   # genes x cell types
      se_mat   <- rep@de_results$gene_fits$s_mat      # genes x cell types

      reg_results[[i]] <- list(mean = mean_mat, se = se_mat)
    }

    # 1. Intersect genes and cell types across replicates
    common_genes <- Reduce(intersect, lapply(reg_results, function(x) rownames(x$mean)))
    common_celltypes <- Reduce(intersect, lapply(reg_results, function(x) colnames(x$mean)))

    if (length(common_genes) == 0 || length(common_celltypes) == 0) {
      warning(paste("No common genes or cell types in region:", reg))
      next
    }

    # Subset all matrices to common genes/celltypes
    for (i in seq_along(reg_results)) {
      colnames(reg_results[[i]]$se) <- colnames(reg_results[[i]]$mean) #Important all should have same naming
      reg_results[[i]]$mean <- reg_results[[i]]$mean[common_genes, common_celltypes, drop = FALSE]
      reg_results[[i]]$se   <- reg_results[[i]]$se[common_genes, common_celltypes, drop = FALSE]
    }

    #####CHANGED TO MAKE IT DYNAMIC AND USE INVERSE-VARIANCE WEIGHTED
    cond1 <- conditions[1]
    cond2 <- conditions[2]

    idx1 <- which(design_matrix$condition == cond1)
    idx2 <- which(design_matrix$condition == cond2)

    mean_cond1 <- matrix(NA, nrow = length(common_genes), ncol = length(common_celltypes),
                         dimnames = list(common_genes, common_celltypes))
    mean_cond2 <- matrix(NA, nrow = length(common_genes), ncol = length(common_celltypes),
                         dimnames = list(common_genes, common_celltypes))
    se_cond1 <- matrix(NA, nrow = length(common_genes), ncol = length(common_celltypes),
                       dimnames = list(common_genes, common_celltypes))
    se_cond2 <- matrix(NA, nrow = length(common_genes), ncol = length(common_celltypes),
                       dimnames = list(common_genes, common_celltypes))

    # Compute inverse-variance weighted mean & SE per gene/celltype
    for (g in common_genes) {
      for (ct in common_celltypes) {
        # condition 1
        means1 <- sapply(idx1, function(i) reg_results[[i]]$mean[g, ct])
        ses1   <- sapply(idx1, function(i) reg_results[[i]]$se[g, ct])
        w1 <- 1 / (ses1^2)
        mean_cond1[g, ct] <- sum(w1 * means1) / sum(w1)
        se_cond1[g, ct]   <- sqrt(1 / sum(w1))

        # condition 2
        means2 <- sapply(idx2, function(i) reg_results[[i]]$mean[g, ct])
        ses2   <- sapply(idx2, function(i) reg_results[[i]]$se[g, ct])
        w2 <- 1 / (ses2^2)
        mean_cond2[g, ct] <- sum(w2 * means2) / sum(w2)
        se_cond2[g, ct]   <- sqrt(1 / sum(w2))
      }
    }

    # Z-test
    Z <- (mean_cond2 - mean_cond1) / sqrt(se_cond2^2 + se_cond1^2)
    pvals <- 2 * (1 - pnorm(abs(Z)))

    # Save CSV per cell type
    for (ct in common_celltypes) {
      df <- data.frame(
        gene = rownames(Z),
        Z = Z[, ct],
        pval = pvals[, ct]
      )
      df$padj <- p.adjust(df$pval, method = "fdr")

    ####OLD HARD-CODED WAY
    #healthy_idx <- which(design_matrix$injured == 0)
    #injured_idx <- which(design_matrix$injured == 1)

    #mean_healthy <- Reduce("+", lapply(healthy_idx, function(i) reg_results[[i]]$mean)) / length(healthy_idx)
    #mean_injured <- Reduce("+", lapply(injured_idx, function(i) reg_results[[i]]$mean)) / length(injured_idx)

    #se_healthy <- Reduce("+", lapply(healthy_idx, function(i) reg_results[[i]]$se)) / length(healthy_idx)
    #se_injured <- Reduce("+", lapply(injured_idx, function(i) reg_results[[i]]$se)) / length(injured_idx)

    #Z <- (mean_injured - mean_healthy) / sqrt(se_injured^2 + se_healthy^2)
    #pvals <- 2 * (1 - pnorm(abs(Z)))

    # Loop over cell types and save each as separate CSV
    #for (ct in colnames(Z)) {
    #  df <- data.frame(
    #    gene = rownames(Z),
    #    Z = Z[, ct],
    #    pval = pvals[, ct]
    #  )
    #  df$padj <- p.adjust(df$pval, method = "fdr")

      outfile <- file.path(regdir, paste0("DEA_", reg, "_", ct, "_Z_CSIDE.csv"))
      write.csv(df, outfile, row.names = FALSE)
    }
  }
}

#-----------------------------
# Function 2: Meta-analysis (Stouffer’s Z)
#-----------------------------
run_region_dea_meta <- function(myRCTD.reps, regions, region_annotations, design_matrix, outdir) {
  for (reg in regions) {
    regdir <- file.path(outdir, reg)
    dir.create(regdir, showWarnings = FALSE, recursive = TRUE)

    cat("Running DEA (meta-analysis) for region:", reg, "\n")

    # Collect results from each replicate
    reg_results <- list()
    for (i in seq_along(myRCTD.reps@RCTD.reps)) {
      sample_name <- names(region_annotations)[i]
      region_barcodes <- names(region_annotations[[sample_name]])[region_annotations[[sample_name]] == reg]

      subset_reps <- subset_RCTD_replicates(myRCTD.reps, region_barcodes)
      rep <- subset_reps@RCTD.reps[[i]]

      mean_mat <- rep@de_results$gene_fits$mean_val
      se_mat   <- rep@de_results$gene_fits$s_mat

      reg_results[[i]] <- list(
        mean = mean_mat,
        se   = se_mat,
        cond = ifelse(design_matrix$injured[i] == 0, "healthy", "injured")
      )
    }
    # Intersections
    common_genes <- Reduce(intersect, lapply(reg_results, function(x) rownames(x$mean)))
    common_celltypes <- Reduce(intersect, lapply(reg_results, function(x) colnames(x$mean)))

    if (length(common_genes) == 0 || length(common_celltypes) == 0) {
      warning(paste("No common genes/celltypes for region:", reg))
      next
    }

    # Subset all matrices to common genes/celltypes
    for (i in seq_along(reg_results)) {
      colnames(reg_results[[i]]$se) <- colnames(reg_results[[i]]$mean) #Important all should have same naming
      reg_results[[i]]$mean <- reg_results[[i]]$mean[common_genes, common_celltypes, drop = FALSE]
      reg_results[[i]]$se   <- reg_results[[i]]$se[common_genes, common_celltypes, drop = FALSE]
    }

    # Pair up healthy vs injured for each replicate and compute Z
    # (assumes equal number of reps per condition)
    all_Zs <- list()
    for (i in seq_along(replicate_Zs)) {
      for (j in seq_along(replicate_Zs)) {
        if (replicate_Zs[[i]]$cond == "healthy" && replicate_Zs[[j]]$cond == "injured") {
          m1 <- replicate_Zs[[i]]$mean
          m2 <- replicate_Zs[[j]]$mean
          s1 <- replicate_Zs[[i]]$se
          s2 <- replicate_Zs[[j]]$se

          Z <- (m2 - m1) / sqrt(s1^2 + s2^2)
          all_Zs[[length(all_Zs)+1]] <- Z
        }
      }
    }

    if (length(all_Zs) == 0) {
      warning(paste("No valid healthy–injured pairs in region:", reg))
      next
    }

    # Convert list of Z matrices to array
    Z_array <- simplify2array(all_Zs)  # genes x celltypes x replicates
    
    # Meta-analysis (Stouffer’s method, unweighted)
    Z_meta <- apply(Z_array, c(1, 2), function(zvec) {
      zvec <- na.omit(zvec)
      if (length(zvec) == 0) return(NA)
      sum(zvec) / sqrt(length(zvec))
    })

    pvals <- 2 * (1 - pnorm(abs(Z_meta)))

    # Loop over cell types and save each as separate CSV
    for (ct in colnames(Z_meta)) {
      df <- data.frame(
        gene = rownames(Z_meta),
        Z = Z_meta[, ct],
        pval = pvals[, ct]
      )
      df$padj <- p.adjust(df$pval, method = "fdr")

      outfile <- file.path(regdir, paste0("DEA_", reg, "_", ct, "_meta_CSIDE.csv"))
      write.csv(df, outfile, row.names = FALSE)
    }
  }
}


#Select regions and ensure drop NA and empty strings
regions <- unique(unlist(region_annotations))
regions <- regions[!is.na(regions) & trimws(regions) != ""]


#####
## FUNCTION: LIMMA
#####


#######################################################

run_limma_dea_long <- function(myRCTD.reps, regions, region_annotations, design_matrix, outdir, pseudocount = 1e-3, already_log = TRUE) {
  if (!requireNamespace("limma", quietly = TRUE)) stop("Please install 'limma'")
  if (!requireNamespace("tools", quietly = TRUE)) stop("Please install 'tools'")

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  results_list <- list()

  # Samples (names of replicates) - ensure present
  all_samples <- names(myRCTD.reps@RCTD.reps)
  if (is.null(all_samples)) {
    all_samples <- paste0("Sample_", seq_along(myRCTD.reps@RCTD.reps))
    names(myRCTD.reps@RCTD.reps) <- all_samples
  }

  for (reg in regions) {
    cat("\n[run_limma_dea] === region:", reg, "===\n")
    regdir <- normalizePath(file.path(outdir, reg), mustWork = FALSE)
    dir.create(regdir, showWarnings = FALSE, recursive = TRUE)

    # collect per-sample mean matrices for this region only when available
    reg_results <- list()

    for (i in seq_along(myRCTD.reps@RCTD.reps)) {
      sample_name <- names(region_annotations)[i]

      rep <- myRCTD.reps@RCTD.reps[[i]]

      mean_mat <- rep@de_results$gene_fits$mean_val   # genes x cell types
      se_mat   <- rep@de_results$gene_fits$s_mat      # genes x cell types

      reg_results[[i]] <- list(sample = sample_name, mean = mean_mat, se = se_mat)
    }

    if (length(reg_results) == 0) {
      warning("[run_limma_dea] No valid samples with gene_fits for region: ", reg)
      next
    }


    # Intersect genes and cell types across VALID samples only
    common_genes <- Reduce(intersect, lapply(reg_results, function(x) rownames(x$mean)))
    common_celltypes <- Reduce(intersect, lapply(reg_results, function(x) colnames(x$mean)))

    if (length(common_genes) == 0 || length(common_celltypes) == 0) {
      warning("[run_limma_dea] No common genes or cell types for region: ", reg)
      next
    }

    cat("[run_limma_dea] Using", length(reg_results), "samples for region", reg,
        " — common genes:", length(common_genes),
        " common celltypes:", length(common_celltypes), "\n")

    # For each cell type, build genes x samples matrix, log-transform and run limma
    for (ct in common_celltypes) {
      cat("  - celltype:", ct, " ... building expression matrix\n")

      # Build matrix genes x samples (columns in same order as sample_names_valid)
      expr_mat <- sapply(reg_results, function(x) {
        # keep only common_genes and the celltype column
        # if celltype missing in any sample, this would have been filtered by common_celltypes
        x$mean[common_genes, ct]
      })
      # sapply returns matrix genes x samples; set colnames using sample field
      #colnames(expr_mat) <- vapply(reg_results, function(x) x$sample, character(1))
      colnames(expr_mat) <- vapply(reg_results, function(x) x$sample, character(1))
      rownames(expr_mat) <- common_genes

      # Remove any columns with all NA or all zero (defensive)
      keep_cols <- colSums(is.finite(expr_mat)) > 0
      if (any(!keep_cols)) {
        warning("    some samples produced only NA/Inf for celltype ", ct, " — dropping them")
        expr_mat <- expr_mat[, keep_cols, drop = FALSE]
      }
      if (ncol(expr_mat) == 0) {
        warning("    no valid columns remain for ", ct, " in region ", reg, " — skipping")
        next
      }

      # --- 🔧 NEW BLOCK: handle NA genes ---
      # Drop genes (rows) that are entirely NA
      expr_mat <- expr_mat[rowSums(is.na(expr_mat)) < ncol(expr_mat), , drop = FALSE]

      # Option 1: drop genes that have *any* NA (strict)
      expr_mat <- expr_mat[complete.cases(expr_mat), , drop = FALSE]

      if (nrow(expr_mat) == 0) {
        warning("    no valid genes remain for celltype ", ct, " in region ", reg, " — skipping")
        next
      }

      # Match design_matrix rows to expr_mat columns
      samples_for_test <- colnames(expr_mat)
      if (!all(samples_for_test %in% rownames(design_matrix))) {
        missing_samples <- setdiff(samples_for_test, rownames(design_matrix))
        stop("[run_limma_dea] design_matrix missing rows for samples: ", paste(missing_samples, collapse = ", "))
      }

      # Subset design matrix to the samples present and preserve order
      design_sub <- design_matrix[samples_for_test, , drop = FALSE]

      # --- 🔧 Normalization block ---
      # If estimates are already in log scale (RCTD), convert back to linear first
      expr_mat <- exp(expr_mat)

      expr_log2 <- log2(expr_mat + pseudocount)

      expr_norm <- limma::normalizeBetweenArrays(expr_log2, method = "quantile")

      expr_use <- expr_norm 

      # CSIDE estimates are already in log scale
      #expr_log <- log2(expr_mat + pseudocount)

      #expr_use <- expr_mat  # no transformation

      # Build design: use intercept model (~ condition) rather than 0+ to simplify coef picking
      # expect design_sub to have column 'condition' (text)
      if (!("condition" %in% colnames(design_sub))) {
        stop("[run_limma_dea] design_matrix must contain a 'condition' column")
      }
      
      #design <- stats::model.matrix(~ condition, data = design_sub)

      # Ensure strict alignment
      design_sub <- design_sub[colnames(expr_use), , drop = FALSE]

      design <- stats::model.matrix(~ condition, data = design_sub)

      cat("\n--- LIMMA DEBUG ---\n")
      cat("Region:", reg, "\n")
      cat("Celltype:", ct, "\n")
      cat("expr_use dim:", dim(expr_use), "\n")
      cat("design dim:", dim(design), "\n")
      cat("expr colnames:\n")
      print(colnames(expr_use))
      cat("design rownames:\n")
      print(rownames(design))
      cat("All names match:",
          all(colnames(expr_use) == rownames(design)), "\n")
      cat("-------------------\n")

      # Fit limma
      fit <- limma::lmFit(expr_use, design)

      # Check residual degrees of freedom: need at least 1 residual df to run eBayes
      # residual df = ncol(expr_log) - qr(design)$rank
      resid_df <- ncol(expr_use) - qr(design)$rank
      if (resid_df <= 0) {
        warning(sprintf("[run_limma_dea] Not enough residual degrees of freedom for region=%s celltype=%s (n=%d, rank=%d). Skipping.",
                        reg, ct, ncol(expr_use), qr(design)$rank))
        next
      }

      # Decide which coefficient to test: if model has >1 columns, pick the last (condition effect)
      design_colnames <- colnames(design)
      if (length(design_colnames) < 2) {
        warning("[run_limma_dea] Design has <2 columns for region=", reg, " ct=", ct, " — skipping")
        next
      }

      # If there are exactly 2 columns (Intercept + conditionX), coef=2 is the condition effect
      coef_index <- 2
      coef_name <- design_colnames[coef_index]

      # empirical Bayes
      fit2 <- limma::eBayes(fit, trend = TRUE)

      # topTable for the specified coefficient
      tt <- limma::topTable(fit2, coef = coef_index, number = Inf, sort.by = "P", adjust.method = "BH")
      tt$gene <- rownames(tt)
      tt <- tt[, c("gene", "logFC", "AveExpr", "P.Value", "adj.P.Val")]

      # Save file
      outfn <- normalizePath(file.path(regdir, paste0("DEA_", reg, "_", tools::file_path_sans_ext(ct), "_limma.csv")), mustWork = FALSE)
      write.csv(tt, outfn, row.names = FALSE)
      cat("    -> wrote (confirmed):", outfn, file.exists(outfn), "\n")

      results_list[[paste(reg, ct, sep = "_")]] <- tt
    } # end celltype loop
  } # end region loop

  return(results_list)
}


# Run averaging method
#dea_avg_results <- run_region_dea_avg(myRCTD.reps, regions, region_annotations, meta.design.matrix, cside_output_dir)

# Run meta-analysis method
#dea_meta_results <- run_region_dea_meta(myRCTD.reps, regions, region_annotations, meta.design.matrix, cside_output_dir)

# Run limma 
#dea_limma_results <- run_limma_dea_long(myRCTD.reps, regions, region_annotations, meta.design.matrix, cside_output_dir)

subset_RCTD_replicates_complete <- function(rep_obj, barcodes) {
  # Defensive checks
  if (!inherits(rep_obj, "RCTD")) {
    stop("Input must be a single RCTD replicate object, not the entire myRCTD.reps container.")
  }

  # Filter barcodes that exist in object
  keep <- colnames(rep_obj@spatialRNA@counts) %in% barcodes
  if (sum(keep) == 0) {
    warning("No matching barcodes found in this replicate; returning original object.")
    return(rep_obj)
  }

  # 1. Subset spatial RNA
  rep_obj@spatialRNA@counts <- rep_obj@spatialRNA@counts[, keep, drop = FALSE]
  rep_obj@spatialRNA@nUMI   <- rep_obj@spatialRNA@nUMI[keep]
  rep_obj@spatialRNA@coords <- rep_obj@spatialRNA@coords[keep, , drop = FALSE]

  # 2. Subset results weights if available
  if (!is.null(rep_obj@results$weights)) {
    rep_obj@results$weights <- rep_obj@results$weights[keep, , drop = FALSE]
  }

  # 6. Return subset replicate
  return(rep_obj)
}

#FUNCTION TO take rep_region, compute per celltype metrics and apply filtering rule

#Min of spots need to be 0.05 percent, not absolute number

filter_rep_region_celltypes <- function(
  rep_region,
  weight_threshold = 0.1,
  min_spot_fraction = 0.05,  # <- now fraction, not absolute number
  max_weight_threshold = 0.2,
  median_weight_threshold = 0.05,
  verbose = TRUE
) {
  if (is.null(rep_region@results$weights)) {
    stop("No weights found in rep_region@results$weights. Adjust accessor.")
  }

  wmat <- rep_region@results$weights
  barcodes_here <- colnames(rep_region@spatialRNA@counts)

  # Subset to barcodes present in rep_region
  if (!is.null(rownames(wmat))) {
    keep_bar <- rownames(wmat) %in% barcodes_here
    wmat <- wmat[keep_bar, , drop = FALSE]
  }

  n_total_spots <- nrow(wmat)
  if (n_total_spots == 0) {
    warning("No barcodes found in weights matrix. Returning original rep_region.")
    return(list(rep_region = rep_region, metrics = NULL))
  }

  # Compute minimum number of spots = fraction * total
  min_spots_per_sample <- ceiling(min_spot_fraction * n_total_spots)

  # Handle sparse vs dense
  is_sparse <- inherits(wmat, "sparseMatrix")

  if (is_sparse) {
    n_spots_vec <- Matrix::colSums(wmat >= weight_threshold)
    max_vec <- apply(wmat, 2, max)
    median_vec <- apply(wmat, 2, median)
    mean_vec <- apply(wmat, 2, mean)
  } else {
    n_spots_vec <- colSums(wmat >= weight_threshold, na.rm = TRUE)
    max_vec <- apply(wmat, 2, max, na.rm = TRUE)
    median_vec <- apply(wmat, 2, median, na.rm = TRUE)
    mean_vec <- colMeans(wmat, na.rm = TRUE)
  }

  metrics_df <- data.frame(
    cell_type = colnames(wmat),
    n_spots = as.integer(n_spots_vec),
    max_weight = as.numeric(max_vec),
    median_weight = as.numeric(median_vec),
    mean_weight = as.numeric(mean_vec),
    stringsAsFactors = FALSE
  )

  # Decide which cell types to keep
  metrics_df$keep <- with(metrics_df,
                          (n_spots >= min_spots_per_sample) &
                          ((max_weight >= max_weight_threshold) |
                           (median_weight >= median_weight_threshold))
  )

  keep_ct <- metrics_df$cell_type[metrics_df$keep]

  # Verbose output
  if (verbose) {
    cat(sprintf("  Total spots = %d | Min spots cutoff (%.1f%%) = %d\n",
                n_total_spots, 100 * min_spot_fraction, min_spots_per_sample))
    cat(sprintf("  Filter summary: total CTs=%d, kept=%d\n",
                nrow(metrics_df), sum(metrics_df$keep)))
    print(metrics_df[order(-metrics_df$n_spots), ])
  }

  # Subset de_results$gene_fits if present
  if (!is.null(rep_region@de_results$gene_fits)) {
    gf <- rep_region@de_results$gene_fits$mean_val
    s_mat <- rep_region@de_results$gene_fits$s_mat

    if (!is.null(gf)) {
      ct_present <- colnames(gf)
      keep_ct_final <- intersect(ct_present, keep_ct)

      if (length(keep_ct_final) == 0) {
        warning("No cell types left after filtering.")
        rep_region@de_results$gene_fits <- NULL
      } else {
        rep_region@de_results$gene_fits$mean_val <- gf[, keep_ct_final, drop = FALSE]
        if (!is.null(s_mat)) {
          keep_ct_smat <- intersect(colnames(s_mat), keep_ct_final)
          if (length(keep_ct_smat) == 0) {
            warning("No matching cell types found in s_mat after filtering; removing s_mat.")
            rep_region@de_results$gene_fits$s_mat <- NULL
          } else {
            rep_region@de_results$gene_fits$s_mat <- s_mat[, keep_ct_smat, drop = FALSE]
          }
        }
      }
    }
  }

  return(list(rep_region = rep_region, metrics = metrics_df))
}

#' Identify condition-unique cell types and save top-N gene metadata
#'
#' @param region_reps RCTD replicates object for the region
#' @param meta_design_matrix Design matrix with sample names and condition column
#' @param gene_expr_dir Directory where per-celltype gene expression CSVs are stored
#' @param region Name of the region
#' @param outdir Directory to save the metadata CSV
#' @return Data frame of condition-unique cell types and top-N gene CSVs (also saved as CSV)
get_condition_unique_ct_topN <- function(region_reps,
                                         meta_design_matrix,
                                         gene_expr_dir,
                                         region,
                                         outdir) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install 'dplyr'")
  
  # Initialize metadata table
  topN_meta <- data.frame(
    region = character(),
    cell_type = character(),
    condition = character(),
    csv_file = character(),
    stringsAsFactors = FALSE
  )
  
  # Samples in region_reps
  samples <- names(region_reps@RCTD.reps)
  
  # Iterate over cell types observed in any replicate
  all_cts <- unique(unlist(lapply(region_reps@RCTD.reps, function(r) {
    if (!is.null(r@de_results$gene_fits)) colnames(r@de_results$gene_fits$mean_val) else NULL
  })))
  
  for (ct in all_cts) {
    # Presence of this CT across replicates
    ct_presence <- sapply(samples, function(samp) {
      r <- region_reps@RCTD.reps[[samp]]
      !is.null(r@de_results$gene_fits) && ct %in% colnames(r@de_results$gene_fits$mean_val)
    })
    
    if (any(ct_presence) && !all(ct_presence)) {
      # Conditions where this CT is present
      condition_for_ct <- unique(meta_design_matrix$condition[match(samples[ct_presence], rownames(meta_design_matrix))])
      
      # Collect CSVs for this CT across replicates where present
      csv_files <- file.path(gene_expr_dir, paste0(samples[ct_presence], "_", ct, "_gene_expr.csv"))
      
      # Add to metadata table
      topN_meta <- rbind(topN_meta, data.frame(
        region = region,
        cell_type = ct,
        condition = paste(condition_for_ct, collapse = ";"),
        csv_file = paste(csv_files, collapse = ";"),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Save metadata CSV
  meta_outfile <- file.path(outdir, paste0("condition_unique_CT_topN_", region, ".csv"))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  write.csv(topN_meta, meta_outfile, row.names = FALSE)
  
  message("[get_condition_unique_ct_topN] Saved metadata for region ", region, " to ", meta_outfile)
  
  return(topN_meta)
}


##RUNNING LIMMA ON SUBSETS OF RCTD

for (reg in regions) {
  #If region is 'Unknown' continue
  if (reg == "Unknown") {
    next
  }

  cat("\n=== Processing region:", reg, "===\n")
  #Copy of object for slicing
  region_reps <- myRCTD.reps 

  for (i in seq_along(myRCTD.reps@RCTD.reps)) {
    rep <- myRCTD.reps@RCTD.reps[[i]]
    sample_name <- names(myRCTD.reps@RCTD.reps)[i]

    # Identify barcodes belonging to this region
    reg_ann_df <- region_annotations[[sample_name]]

    # Select region-specific barcodes
    region_barcodes <- rownames(reg_ann_df)[reg_ann_df[,] == reg]

    cat("Sample:", sample_name,
        "| total RCTD barcodes:", ncol(rep@spatialRNA@counts),
        "| region-specific:", length(region_barcodes), "\n")

    if (length(region_barcodes) == 0) {
      cat("  → Skipping: no barcodes in this region\n")
      next
    }

    # Subset replicate to region
    rep_region <- subset_RCTD_replicates_complete(rep, region_barcodes)

    # Filter unsupported cell types BEFORE running CSIDE
    # Run CSIDE intercept on subset
    # celltype_threshold = 0: we include all celltypes 
    # weight_threshold = 0.1: per-spot per cell-type cutoff to be included for analysis
    
    res <- filter_rep_region_celltypes(
      rep_region,
      weight_threshold = 0.1,        # same as your per-spot cutoff
      min_spot_fraction = 0.05,     # keep cell types that are in >=5% spots
      max_weight_threshold = 0.2,   # optional: drop celltypes with very small max weight
      median_weight_threshold = 0.05,
      verbose = TRUE
    )
    
    # 2. Extract the filtered object (list)
    rep_region_filtered <- res$rep_region
    metrics_df <- res$metrics
    ## Identify kept cell types
    kept_cts <- metrics_df$cell_type[metrics_df$keep]

    # Filter reference inside the replicate (there is no ref in the replicate)
    #ref_info <- rep_region_filtered@reference@cell_type_info$info
    #ref_means <- rep_region_filtered@reference@cell_type_info$mean_profile

    #rep_region_filtered@reference@cell_type_info$info <-
    #  ref_info[ref_info$cell_type %in% kept_cts, , drop = FALSE]

    #rep_region_filtered@reference@cell_type_info$mean_profile <-
    #  ref_means[, kept_cts, drop = FALSE]

    #cat("Filtered reference to", length(kept_cts), "cell types:", paste(kept_cts, collapse=", "), "\n")

    # (optional) Inspect metrics to see which cell types were dropped, run CSIDE only with present celltypes
    head(res$metrics)

    rep_region <- run.CSIDE.intercept(
      rep_region_filtered,
      cell_types = kept_cts,
      cell_type_threshold = 0,
      doublet_mode = FALSE,
      weight_threshold = 0.1
    )

    # Remove genes that did not converge for particular celltypes (adds NAs)
    
    expr_mat <- rep_region@de_results$gene_fits$mean_val
    con_ok   <- rep_region@de_results$gene_fits$con_mat
    common_celltypes <- intersect(colnames(expr_mat), colnames(con_ok))

    # 2. Subset both matrices to these common cell types
    expr_mat_aligned <- expr_mat[, common_celltypes, drop = FALSE]
    con_ok_aligned <- con_ok[, common_celltypes, drop = FALSE]

    # 3. Mask non-converged entries
    expr_mat_filtered <- expr_mat_aligned
    expr_mat_filtered[!con_ok_aligned] <- NA

    # Assign back to rep_region
    rep_region@de_results$gene_fits$mean_val <- expr_mat_filtered

    # --- Save filtered region object ---
    region_dir <- file.path(cside_output_dir, reg)
    dir.create(region_dir, recursive = TRUE, showWarnings = FALSE)

    saveRDS(
      rep_region,
      file.path(region_dir, paste0("myRCTD_reps_afterIntercept_filtered_", reg, "_", sample_name, ".rds"))
    )

    # --- Save filtering metrics ---
    if (!is.null(metrics_df)) {
      write.csv(
        metrics_df,
        file.path(region_dir, paste0("filter_metrics_", sample_name, "_", reg, ".csv")),
        row.names = FALSE
      )
    }

    # --- Save per-celltype gene expression tables ---
    outdir <- file.path(output_base_dir, "RCTD", sample_name)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    gene_expr_dir <- file.path(outdir, "gene_expr_ct", reg)
    dir.create(gene_expr_dir, recursive = TRUE, showWarnings = FALSE)

    if (!is.null(rep_region@de_results$gene_fits)) {
      gene_expr_list <- rep_region@de_results$gene_fits$mean_val
      #New lines:    
      kept_cts <- metrics_df$cell_type[metrics_df$keep]
      ct_to_save <- intersect(colnames(gene_expr_list), kept_cts)
      ########
      #for (ct in colnames(gene_expr_list)) {
      for (ct in ct_to_save) {
        expr_mat <- gene_expr_list[, ct]
        #expr_cside_linear <- exp(expr_mat)
        # don't exponentiate, treat them as they are
        expr_cside_linear <- expr_mat
        outfile <- file.path(gene_expr_dir, paste0(sample_name, "_", ct, "_gene_expr.csv"))
        safe_write_matrix_csv(expr_cside_linear, outfile)
      }
    }

    # Example call
    gene_expr_dir <- file.path(output_base_dir, "RCTD", sample_name, "gene_expr_ct", reg) 
    topN_meta <- get_condition_unique_ct_topN(
      region_reps = region_reps,
      meta_design_matrix = meta.design.matrix,
      gene_expr_dir = gene_expr_dir,
      region = reg,
      outdir = file.path(cside_output_dir, reg)
    )

    # Save filtered replicate into list for DEA
    region_reps@RCTD.reps[[i]] <- rep_region
  }

  # Run limma once per region
  dea_limma_results <- run_limma_dea_long(
    region_reps,
    regions = reg,
    region_annotations,
    design_matrix = meta.design.matrix,
    outdir = normalizePath(cside_output_dir, mustWork = FALSE)
  )
}





###########################################

##NEW APROACH: OPTION WITH <3 REPS:

#cs_intercept_results <- list()
#for (s in samples) {
#  cs <- run.CSIDE.intercept(RCTD_list[[s]], cell_type_threshold = 0, doublet_mode = FALSE, weight_threshold = 0.1)
  # extract a data.frame with mean and se per gene x celltype
#  cs_intercept_results[[s]] <- extract_mean_se(cs) # you need to implement/inspect how to extract; e.g. cs@results...
#}

# Now compute group-level mean and SE from per-replicate means+SEs using inverse-variance weighting:
#group_meta <- function(sample_names) {
  # For each gene x celltype:
#  genes_ct <- cs_intercept_results[[sample_names[1]]]$feature # example
#  out <- data.frame(feature = genes_ct, mean = NA_real_, se = NA_real_, stringsAsFactors = FALSE)
#  for (i in seq_along(genes_ct)) {
    # collect per-sample estimates and SEs
#    ests <- sapply(sample_names, function(s) cs_intercept_results[[s]]$mean[i])
#    ses  <- sapply(sample_names, function(s) cs_intercept_results[[s]]$se[i])
#    w <- 1 / (ses^2)
#    mu_hat <- sum(w * ests) / sum(w)
#    se_hat <- sqrt(1 / sum(w))
#    out$mean[i] <- mu_hat
#    out$se[i] <- se_hat
#  }
#  out
#}

#group1 <- group_meta(samples_by_cond[[1]])
#group2 <- group_meta(samples_by_cond[[2]])

# Z-test per gene x celltype:
#z <- (group1$mean - group2$mean) / sqrt(group1$se^2 + group2$se^2)
#pval <- 2 * pnorm(-abs(z))




##The code below dont work
#region_specific_dea <- list()

#for (sample in sample_names) {
#  delineation <- region_annotations[[sample]]
#  barcodes <- names(delineation)
#  regions <- unique(delineation)
  
#  for (reg in regions) {
#    cat("Running DEA for", sample, "region:", reg, "\n")
    
    # Select barcodes for this region
#    region_barcodes <- barcodes[delineation == reg]
    
    # Subset the CSIDE object
#    subset_reps <- subset_RCTD_replicates(myRCTD.reps, region_barcodes)
    
    # Run CSIDE DEA for this region
    #cside_dea_reg <- CSIDE.population.inference(myRCTD.reps, fdr = 0.01,
    #                                           meta = TRUE,
    #                                           meta.design.matrix = meta.design.matrix,
    #                                           meta.test_var = 'injured'
    #                                           )

    # Run Z Test 


    # Store results as csv
#    region_specific_dea[[paste0(sample, "_", reg)]] <- de_res
#  }
#}

#saveRDS(region_specific_dea, file.path(cside_output_dir, "region_specific_DEA_results.rds"))



# ---------------------------
# Build meta design matrix
# ---------------------------
#meta.design.matrix <- data.frame(
#  intercept = 1,
#  condition = myRCTD.reps@group_ids
#)

#rownames(meta.design.matrix) <- replicate_names


# --- Step 0: sanitize cell_type_info in all replicates ---
#sanitize_cell_type_info <- function(rctd) {
#  mat <- rctd@cell_type_info$info[[1]]
#  mat[is.na(mat)] <- 0                  # replace NAs
#  zero_rows <- rowSums(mat) == 0
#  mat[zero_rows, ] <- 1e-6              # small non-zero to avoid division by zero
#  rctd@cell_type_info$info[[1]] <- mat
#  return(rctd)
#}

#myRCTD.reps@RCTD.reps <- lapply(myRCTD.reps@RCTD.reps, sanitize_cell_type_info)

#For now it does now work ->

# ---------------------------
# Run CSIDE meta regression
# ---------------------------
#myRCTD.reps <- CSIDE.population.inference(
#  myRCTD.reps,
#  fdr = 0.1,
#  log_fc_thresh = 0.4,
#  meta = TRUE,
#  meta.design.matrix = meta.design.matrix,
#  meta.test_var = "condition",
#  CT.PROP = 0 #lower threshold to include more genes
#)

#saveRDS(myRCTD.reps, file.path(cside_output_dir, "myRCTD_reps_pop_inference.rds"))

# ---------------------------
# Save results
# ---------------------------
#sig_genes <- myRCTD.reps@population_sig_gene_list
#de_results <- myRCTD.reps@population_de_results

#saveRDS(sig_genes, file.path(cside_output_dir, "sig_genes.rds"))
#saveRDS(de_results, file.path(cside_output_dir, "de_results.rds"))

#cat("✅ CSIDE pipeline complete. Results saved to", cside_output_dir, "\n")


# Dummy pixel-level covariates (to satisfy full-rank check)
#X.replicates <- lapply(seq_along(myRCTD.reps@RCTD.reps), function(i) {
#  barcodes <- colnames(myRCTD.reps@RCTD.reps[[i]]@spatialRNA@counts)
#  matrix(1, nrow = length(barcodes), ncol = 1, dimnames = list(barcodes, "dummy"))
#})

# 0 = control, 1 = injured, per replicate
#explanatory.variable.replicates <- lapply(seq_along(myRCTD.reps@RCTD.reps), function(i) {
#  val <- if (i <= 2) 0 else 1
  # a single-column matrix with one row is fine
#  matrix(val, nrow = 1, ncol = 1, dimnames = list("replicate", "injured"))
#})

#RCTD_controls <- myRCTD.reps  # fake object to satisfy the internal check

#i <- count_cell_types(
#  RCTD_controls@RCTD.reps[[1]],
#  colnames(RCTD_controls@RCTD.reps[[1]]@spatialRNA@counts),
#  cell_types,
#  cell_type_threshold = 25,
#  doublet_mode = T,
#  weight_threshold = 0.75
#)

#print(i)


# Run population-level C-SIDE across replicates
#RCTD_controls <- run.CSIDE.replicates(
#  RCTD_controls,
#  de_mode = "general",
#  X.replicates = X.replicates,
#  explanatory.variable.replicates = explanatory.variable.replicates,
#  cell_types = cell_types,
#  cell_type_threshold = 10,
#  gene_threshold = 5e-05,
#  doublet_mode = TRUE,
#  fdr = 0.25,
#  population_de = TRUE                               # optional: run meta-analysis immediately
#)

#saveRDS(RCTD_controls,file.path(cside_output_dir,'myRCTDde_reps_afterCSIDE.rds'))

#RCTD_controls <- CSIDE.population.inference(RCTD_controls, fdr = 0.25)

#saveRDS(RCTD_controls,file.path(cside_output_dir,'myRCTDde_reps_pop_inference.rds'))

