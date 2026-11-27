# bin/masigpro_utils.R
# Shared helpers for masigpro_expression.R and masigpro_abundance.R.
# Sourced from the calling script's directory (both live in /usr/local/bin in the image).

suppressPackageStartupMessages({
  library(maSigPro)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || identical(a, "")) b else a

split_arg <- function(x) {
  if (is.null(x) || !nzchar(x) || x %in% c("NONE", "NO_FILE", "ALL_SPOTS")) return(character(0))
  trimws(strsplit(x, ",")[[1]])
}

# -----------------------------------------------------------------------------
# Sample metadata: CSV with columns sample, condition, time, group
#   group : the time series a sample belongs to (e.g. SCI, bPAC, SHAM)
#   the baseline group (e.g. SHAM / healthy) is the shared time-zero of every series
# -----------------------------------------------------------------------------
read_meta <- function(path) {
  meta <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("sample", "condition", "time", "group")
  miss <- setdiff(need, colnames(meta))
  if (length(miss)) stop("sample_metadata is missing column(s): ", paste(miss, collapse = ", "),
                         " (expected ", paste(need, collapse = ", "), ")")
  meta$sample <- as.character(meta$sample)
  meta$group  <- as.character(meta$group)
  meta$time   <- as.numeric(meta$time)
  if (anyNA(meta$time)) stop("sample_metadata: 'time' must be numeric for every sample")
  if (anyDuplicated(meta$sample)) stop("sample_metadata: duplicated sample ids")
  meta
}

# -----------------------------------------------------------------------------
# Auto-detected maSigPro design.
#   - series groups = every group except the baseline, ordered by group_order
#   - baseline samples get 1 in every series column (shared time zero, as in
#     maSigPro's edesignCT example); replicates = same series + time
#   - first series is the reference: get.siggenes(vars = "groups") returns
#     "<ref>" (its own time profile) and "<other>vs<ref>" for each other series
#   - degree = min(requested, n distinct times - 1); needs >= 3 distinct times
# Returns list(edesign, design, series, degree) or list(skip = "<reason>").
# -----------------------------------------------------------------------------
build_design <- function(meta, baseline_group = "", group_order = character(0), degree = 2) {
  has_base <- nzchar(baseline_group) && baseline_group %in% meta$group
  if (nzchar(baseline_group) && !has_base)
    warning("baseline_group '", baseline_group, "' not among the sample groups: no shared time zero")

  series <- unique(meta$group[meta$group != baseline_group])
  series <- c(intersect(group_order, series), setdiff(series, group_order))
  if (length(series) == 0) return(list(skip = "no non-baseline group in the samples"))

  times <- sort(unique(meta$time))
  if (length(times) < 3)
    return(list(skip = paste0("only ", length(times), " distinct timepoint(s) (",
                              paste(times, collapse = ", "), "); maSigPro needs >= 3")))
  if (has_base) {
    bt <- unique(meta$time[meta$group == baseline_group])
    if (length(bt) != 1 || bt != min(meta$time))
      warning("baseline samples are expected at a single, earliest timepoint; found: ",
              paste(bt, collapse = ", "))
  }
  for (s in series) {
    nt <- length(unique(meta$time[meta$group %in% c(s, baseline_group)]))
    if (nt < 2) return(list(skip = paste0("series '", s, "' has only ", nt, " timepoint(s)")))
  }
  deg <- min(degree, length(times) - 1)

  rep_key <- ifelse(meta$group == baseline_group, paste("BASELINE", meta$time),
                    paste(meta$group, meta$time))
  edesign <- data.frame(Time = meta$time,
                        Replicate = as.integer(factor(rep_key, levels = unique(rep_key))))
  for (s in series) {
    edesign[[make.names(s)]] <- as.numeric(meta$group == s | (has_base & meta$group == baseline_group))
  }
  rownames(edesign) <- meta$sample

  design <- make.design.matrix(edesign, degree = deg)
  rownames(design$dis) <- rownames(edesign)
  list(edesign = edesign, design = design, series = make.names(series), degree = deg)
}

describe_design <- function(d, meta) {
  cat("  maSigPro design: degree", d$degree, "| reference series:", d$series[1],
      "| other series:", if (length(d$series) > 1) paste(d$series[-1], collapse = ", ") else "none", "\n")
  print(cbind(d$edesign, condition = meta$condition[match(rownames(d$edesign), meta$sample)]))
}

# -----------------------------------------------------------------------------
# Two-step maSigPro. counts = TRUE -> negative-binomial GLM, else gaussian.
# Returns list(fit, tfit, sigs) with NULLs where a step had nothing to work on.
# -----------------------------------------------------------------------------
run_masigpro <- function(mat, d, counts, Q = 0.05, alfa = 0.05, rsq = 0.4) {
  mat <- as.matrix(mat[, rownames(d$design$dis), drop = FALSE])
  min_obs <- min(6, ncol(mat))
  fit <- p.vector(mat, d$design, Q = Q, MT.adjust = "BH", counts = counts, min.obs = min_obs)
  if (is.null(fit$i) || fit$i == 0) return(list(fit = fit, tfit = NULL, sigs = NULL))
  tfit <- T.fit(fit, step.method = "backward", alfa = alfa)
  sigs <- tryCatch(get.siggenes(tfit, rsq = rsq, vars = "groups"),
                   error = function(e) { warning("get.siggenes: ", conditionMessage(e)); NULL })
  list(fit = fit, tfit = tfit, sigs = sigs)
}

n_sig <- function(contrast) {
  sp <- contrast$sig.profiles
  if (is.null(sp) || length(sp) == 0) 0L else NROW(sp)
}

# Write the tables every maSigPro run produces; returns a data.frame
# (contrast, feature) of significant features.
write_masigpro_tables <- function(res, outdir, item = "gene") {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  pv <- data.frame(feature = rownames(res$fit$p.vector) %||% names(res$fit$p.vector),
                   p_value = as.numeric(res$fit$p.vector),
                   p_adj_BH = as.numeric(res$fit$p.adjusted))
  write.table(pv, file.path(outdir, "pvector_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  if (is.null(res$tfit)) {
    cat("  no ", item, "s pass p.vector at the chosen Q\n", sep = "")
    return(data.frame(contrast = character(0), feature = character(0)))
  }
  sol <- as.data.frame(res$tfit$sol)
  write.table(cbind(feature = rownames(sol), sol), file.path(outdir, "tfit_sol.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  if (is.null(res$sigs)) return(data.frame(contrast = character(0), feature = character(0)))
  out <- list()
  for (cn in names(res$sigs$sig.genes)) {
    ct <- res$sigs$sig.genes[[cn]]
    n <- n_sig(ct)
    cat("  contrast ", cn, ": ", n, " significant ", item, "s\n", sep = "")
    if (n == 0) next
    write.table(cbind(feature = rownames(ct$sig.profiles), ct$sig.profiles),
                file.path(outdir, paste0("sig_profiles_", cn, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(cbind(feature = rownames(ct$sig.pvalues), ct$sig.pvalues),
                file.path(outdir, paste0("sig_pvalues_", cn, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    out[[cn]] <- data.frame(contrast = cn, feature = rownames(ct$sig.profiles))
  }
  if (length(out)) do.call(rbind, out) else data.frame(contrast = character(0), feature = character(0))
}

# -----------------------------------------------------------------------------
# Delineation / barcode helpers
# -----------------------------------------------------------------------------
read_delineation <- function(delineation_dir, sample) {
  f <- file.path(delineation_dir, paste0(sample, "_manual_delineation.csv"))
  if (!file.exists(f)) { warning("missing delineation for ", sample); return(NULL) }
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  setNames(as.character(d$manual_delineation), as.character(d[[1]]))
}

# Map delineation barcodes onto data barcodes. Exact match first; if most
# don't match, fall back to the part after the last "_" (sample prefixes).
align_barcodes <- function(delin_bc, data_bc) {
  if (mean(delin_bc %in% data_bc) >= 0.5) return(delin_bc)
  strip <- function(x) sub(".*_", "", x)
  key <- strip(data_bc)
  if (anyDuplicated(key)) stop("barcodes are not unique after stripping prefixes")
  data_bc[match(strip(delin_bc), key)]
}

# spots (data barcodes) per area for one sample; areas = character(0) -> ALL_SPOTS
spots_by_area <- function(areas, delineation_dir, sample, data_bc) {
  if (length(areas) == 0) return(list(ALL_SPOTS = data_bc))
  lab <- read_delineation(delineation_dir, sample)
  if (is.null(lab)) return(list())
  mapped <- align_barcodes(names(lab), data_bc)
  res <- lapply(areas, function(a) unique(na.omit(mapped[lab == a])))
  names(res) <- areas
  res
}
