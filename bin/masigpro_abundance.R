#!/usr/bin/env Rscript
# bin/masigpro_abundance.R
# Time-course cell-type abundance per area with maSigPro on CLR compositions
# (paper: "Differential cell-type abundance analysis"). Works on RCTD or
# cell2location spot weights.
#   1. spot level: weights / spot total, entries < --spot_min set to 0, re-normalised
#   2. area level per sample: mean per-spot proportion of each cell type
#      (sum of spot proportions / n spots: ROI-size independent)
#   3. drop cell types with mean proportion <= --ct_min across samples in the
#      area, re-close, CLR with --pseudocount
#   4. maSigPro (gaussian), degree 2, p.vector Q (BH), T.fit backward (alfa),
#      get.siggenes R2 > rsq, vars = "groups"
#   5. log2FC vs baseline per group x time = (mean CLR - baseline mean CLR) / ln 2,
#      heatmap and CLR trajectories of significant cell types
# Design is auto-detected from sample_metadata (sample, condition, time, group).

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
})

.this_dir <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else "."
})
source(file.path(.this_dir, "masigpro_utils.R"))

option_list <- list(
  make_option("--method", type = "character", default = "rctd",
              help = "rctd or cell2loc (sets the default weights file pattern)"),
  make_option("--weights_dir", type = "character",
              help = "RCTD/ (RCTD/<s>/<s>_RCTD_weights.csv) or cell2location_map/ (<s>/<s>_cell2loc_weights.csv)"),
  make_option("--weights_pattern", type = "character", default = "",
              help = "Override: path under weights_dir with {sample} placeholders"),
  make_option("--delineation_dir", type = "character", default = "NO_DELINEATION"),
  make_option("--sample_metadata", type = "character"),
  make_option("--baseline_group", type = "character", default = ""),
  make_option("--group_order", type = "character", default = ""),
  make_option("--areas", type = "character", default = "ALL_SPOTS"),
  make_option("--spot_min", type = "double", default = 0.05),
  make_option("--ct_min", type = "double", default = 0.02),
  make_option("--pseudocount", type = "double", default = 1e-4),
  make_option("--degree", type = "integer", default = 2),
  make_option("--Q", type = "double", default = 0.05),
  make_option("--alfa", type = "double", default = 0.05),
  make_option("--rsq", type = "double", default = 0.4),
  make_option("--outdir", type = "character", default = "masigpro_abundance")
)
opt <- parse_args(OptionParser(option_list = option_list))
opt$method <- tolower(opt$method)
pattern <- if (nzchar(opt$weights_pattern)) opt$weights_pattern else switch(opt$method,
  rctd     = "{sample}/{sample}_RCTD_weights.csv",
  cell2loc = "{sample}/{sample}_cell2loc_weights.csv",
  stop("--method must be rctd or cell2loc"))
outdir <- file.path(opt$outdir, opt$method)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta  <- read_meta(opt$sample_metadata)
areas <- split_arg(opt$areas)
if (length(areas) && opt$delineation_dir == "NO_DELINEATION") {
  warning("areas given but no delineation_dir: using whole-sample composition")
  areas <- character(0)
}
group_order <- split_arg(opt$group_order)
base <- opt$baseline_group

d_all <- build_design(meta, base, group_order, opt$degree)
if (!is.null(d_all$skip)) {
  msg <- paste("maSigPro abundance skipped:", d_all$skip)
  cat(msg, "\n"); writeLines(msg, file.path(outdir, "SKIPPED.txt"))
  quit(status = 0)
}

# -----------------------------------------------------------------------------
# 1-2. Per-sample area compositions
# -----------------------------------------------------------------------------
area_keys <- if (length(areas)) areas else "ALL_SPOTS"
comp <- setNames(lapply(area_keys, function(a) list()), area_keys)
spot_log <- list()

for (s in meta$sample) {
  f <- file.path(opt$weights_dir, gsub("{sample}", s, pattern, fixed = TRUE))
  if (!file.exists(f)) { warning("missing weights for ", s, ": ", f); next }
  w <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
  w[is.na(w) | w < 0] <- 0
  rs <- rowSums(w)
  w <- w[rs > 0, , drop = FALSE] / rs[rs > 0]
  w[w < opt$spot_min] <- 0
  rs <- rowSums(w)
  w <- w[rs > 0, , drop = FALSE] / rs[rs > 0]

  sp <- spots_by_area(areas, opt$delineation_dir, s, rownames(w))
  for (a in names(sp)) {
    bc <- sp[[a]]
    spot_log[[length(spot_log) + 1]] <- data.frame(area = a, sample = s, n_spots = length(bc))
    if (length(bc) == 0) { warning("no spots in ", a, " for ", s); next }
    comp[[a]][[s]] <- colSums(w[bc, , drop = FALSE]) / length(bc)
  }
}
if (length(spot_log))
  write.table(do.call(rbind, spot_log), file.path(outdir, "spots_per_area_sample.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

clr <- function(P, eps) { L <- log(P + eps); sweep(L, 2, colMeans(L)) }

# group x time label used for log2FC columns / plots
meta$series_time <- ifelse(meta$group == base, base, paste0(meta$group, "_", meta$time))

summary_rows <- list(); lfc_all <- list(); sig_all <- list(); clr_all <- list()
for (a in area_keys) {
  cat("\n=== Area:", a, "===\n")
  adir <- file.path(outdir, a)
  dir.create(adir, recursive = TRUE, showWarnings = FALSE)
  if (length(comp[[a]]) < 2) { warning("fewer than 2 samples with spots in ", a); next }

  cts <- Reduce(union, lapply(comp[[a]], names))
  P <- sapply(comp[[a]], function(v) { x <- v[cts]; x[is.na(x)] <- 0; x })
  rownames(P) <- cts
  m <- meta[match(colnames(P), meta$sample), ]

  # 3. filter, re-close, CLR
  mean_prop <- rowMeans(P)
  keep <- mean_prop > opt$ct_min
  write.table(data.frame(celltype = cts, mean_proportion = mean_prop, retained = keep),
              file.path(adir, "celltype_filter.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(cbind(celltype = rownames(P), P), file.path(adir, "proportions.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  if (sum(keep) < 2) { warning("fewer than 2 cell types above ct_min in ", a); next }
  Pk <- P[keep, , drop = FALSE]
  Pk <- sweep(Pk, 2, colSums(Pk), "/")
  C <- clr(Pk, opt$pseudocount)
  write.table(cbind(celltype = rownames(C), C), file.path(adir, "clr.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  clr_long <- data.frame(area = a, celltype = rep(rownames(C), ncol(C)),
                         sample = rep(colnames(C), each = nrow(C)), clr = as.vector(C))
  clr_long <- merge(clr_long, m[, c("sample", "condition", "time", "group")], by = "sample")
  clr_all[[a]] <- clr_long

  # 4. maSigPro
  d <- build_design(m, base, group_order, opt$degree)
  sig <- data.frame(contrast = character(0), feature = character(0))
  if (!is.null(d$skip)) {
    msg <- paste("skipped:", d$skip); cat(" ", msg, "\n")
    writeLines(msg, file.path(adir, "SKIPPED.txt"))
  } else {
    describe_design(d, m)
    res <- run_masigpro(C, d, counts = FALSE, Q = opt$Q, alfa = opt$alfa, rsq = opt$rsq)
    sig <- write_masigpro_tables(res, adir, item = "cell type")
  }
  if (nrow(sig)) sig_all[[a]] <- cbind(area = a, sig)

  # 5. log2FC vs baseline
  if (base %in% m$group) {
    bmean <- rowMeans(C[, m$group == base, drop = FALSE])
    st <- unique(m$series_time[m$group != base])
    lfc <- sapply(st, function(x) (rowMeans(C[, m$series_time == x, drop = FALSE]) - bmean) / log(2))
    lfc <- matrix(lfc, nrow = nrow(C), dimnames = list(rownames(C), st))
    lfc_all[[a]] <- data.frame(area = a, celltype = rep(rownames(lfc), ncol(lfc)),
                               series_time = rep(colnames(lfc), each = nrow(lfc)),
                               log2fc = as.vector(lfc))
  }

  summary_rows[[a]] <- data.frame(
    area = a, n_samples = ncol(P), n_celltypes = length(cts), n_retained = sum(keep),
    n_sig_any = length(unique(sig$feature)),
    sig_vs_reference = paste(unique(sig$feature[grepl("vs", sig$contrast)]), collapse = ";"))
}

# -----------------------------------------------------------------------------
# Combined tables and plots
# -----------------------------------------------------------------------------
if (length(summary_rows)) {
  s <- do.call(rbind, summary_rows)
  write.table(s, file.path(outdir, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n"); print(s)
}
sig_df <- if (length(sig_all)) do.call(rbind, sig_all) else NULL
if (!is.null(sig_df))
  write.table(sig_df, file.path(outdir, "significant_all_areas.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
clr_df <- if (length(clr_all)) do.call(rbind, clr_all) else NULL
if (!is.null(clr_df))
  write.table(clr_df, file.path(outdir, "clr_all_areas_long.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Heatmap: rows = every cell type seen, x marks cell types not retained in an
# area, * marks significant "<series>vs<reference>" trajectories.
if (length(lfc_all)) {
  lfc_df <- do.call(rbind, lfc_all)
  write.table(lfc_df, file.path(outdir, "log2fc_vs_baseline.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  all_ct <- sort(unique(unlist(lapply(comp, function(l) unlist(lapply(l, names))))))
  nb <- meta[meta$group != base, ]
  g_lv <- c(intersect(group_order, nb$group), setdiff(sort(unique(nb$group)), group_order))
  st_levels <- unique(nb$series_time[order(match(nb$group, g_lv), nb$time)])
  grid <- expand.grid(area = names(lfc_all), celltype = all_ct, series_time = st_levels,
                      stringsAsFactors = FALSE)
  hm <- merge(grid, lfc_df, all.x = TRUE)
  hm$label <- ifelse(is.na(hm$log2fc), "x", "")
  vs <- if (!is.null(sig_df)) unique(sig_df[grepl("vs", sig_df$contrast), c("area", "feature")]) else NULL
  if (!is.null(vs) && nrow(vs)) {
    star <- paste(hm$area, hm$celltype) %in% paste(vs$area, vs$feature)
    hm$label[star & !is.na(hm$log2fc)] <- "*"
  }
  hm$fill <- pmax(pmin(hm$log2fc, 3), -3)
  hm$series_time <- factor(hm$series_time, levels = st_levels)
  hm$celltype <- factor(hm$celltype, levels = rev(all_ct))

  p <- ggplot(hm, aes(series_time, celltype, fill = fill)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = label), size = 3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-3, 3), na.value = "grey85",
                         name = paste0("log2FC CLR\nvs ", base)) +
    facet_wrap(~ area, nrow = 1) +
    labs(x = NULL, y = NULL,
         title = paste0(toupper(opt$method), ": cell-type abundance vs ", base),
         subtitle = "* significant series-vs-reference trajectory (maSigPro); x = not retained in area") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  hw <- 2 + 1.6 * length(lfc_all) * max(1, length(st_levels)) / 2
  hh <- 2 + 0.22 * length(all_ct)
  for (ext in c("pdf", "png"))
    ggsave(file.path(outdir, paste0("heatmap_log2fc_vs_baseline.", ext)), p,
           width = hw, height = hh, dpi = 200, limitsize = FALSE)
}

# Trajectories of significant cell types: CLR over time per series, mean +/- SE,
# individual samples overlaid; baseline samples drawn as time zero of every series.
if (!is.null(sig_df) && !is.null(clr_df)) {
  series <- setdiff(unique(meta$group), base)
  pdf(file.path(outdir, "trajectories_significant.pdf"), width = 9, height = 6)
  for (a in unique(sig_df$area)) {
    for (ct in unique(sig_df$feature[sig_df$area == a])) {
      x <- clr_df[clr_df$area == a & clr_df$celltype == ct, ]
      xb <- x[x$group == base, ]
      x <- x[x$group != base, ]
      if (nrow(xb)) x <- rbind(x, do.call(rbind, lapply(series, function(g) transform(xb, group = g))))
      agg <- do.call(rbind, lapply(split(x, list(x$group, x$time), drop = TRUE), function(z)
        data.frame(group = z$group[1], time = z$time[1], mean = mean(z$clr),
                   se = if (nrow(z) > 1) sd(z$clr) / sqrt(nrow(z)) else 0)))
      cons <- paste(sig_df$contrast[sig_df$area == a & sig_df$feature == ct], collapse = ", ")
      p <- ggplot() +
        geom_ribbon(data = agg, aes(time, ymin = mean - se, ymax = mean + se, fill = group), alpha = 0.2) +
        geom_line(data = agg, aes(time, mean, colour = group), linewidth = 1) +
        geom_point(data = x, aes(time, clr, colour = group), position = position_jitter(width = 0.3, height = 0),
                   alpha = 0.7) +
        labs(title = paste(ct, "-", a), subtitle = paste("significant in:", cons),
             x = "time", y = "CLR abundance") +
        theme_bw()
      print(p)
    }
  }
  dev.off()
}
