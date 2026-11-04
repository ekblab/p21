## p21 | IPA upstream regulators and biological functions
## Interactive ComplexHeatmap exploration of the ORIGINAL IPA Comparison exports.
## Run in RStudio from project root: source("R/08_IPA_upstream_heatmap.R")
## Does not run IPA or modify the original input tables.

## 1 | Files and dependencies -----------------------------------------------------
required <- c("here", "ComplexHeatmap", "circlize", "grid", "ekbSeq")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
here::i_am("R/08_IPA_upstream_heatmap.R")
## Keep the ORIGINAL IPA comparison exports where they are.
input_dir <- here::here("Results", "07_IPA")
## Figures produced by script 08 belong to its own Results folder.
figure_dir <- here::here("Results", "08_IPA_visualisations", "figures")

## 'KO_minus_WT' = six formal difference-of-differences (default).
## 'irradiated_vs_control' = twelve within-genotype irradiation responses.
comparison <- "KO_minus_WT"

## Export PNG + SVG after showing both heatmaps in the RStudio Plots pane.
export_figures <- TRUE

read_ipa_pair <- function(prefix, comparison) {
  files <- file.path(input_dir, paste0(prefix, "_comparison_", comparison,
                                       c("_zscores.txt", "_bh_pval.txt")))
  if (any(!file.exists(files))) {
    stop("Missing original IPA export(s):\n", paste(files[!file.exists(files)], collapse = "\n"))
  }
  reader <- function(path) {
    x <- utils::read.delim(path, skip = 2L, check.names = FALSE, quote = "",
                           comment.char = "", na.strings = c("N/A", "NA", ""),
                           stringsAsFactors = FALSE, fileEncoding = "UTF-8")
    if (ncol(x) < 2L) stop("Unexpected IPA table: ", basename(path))
    x
  }
  zs <- reader(files[1L])
  ps <- reader(files[2L])
  if (!identical(names(zs), names(ps)) || !identical(names(zs)[1L],
                  if (prefix == "upstream_regulators") "Upstream Regulators" else "Diseases and Bio Functions")) {
    stop("IPA header mismatch: ", prefix)
  }
  if (!setequal(zs[[1L]], ps[[1L]])) stop("IPA feature names differ between z-score and BH exports")
  list(z = zs, p = ps)
}

ipa_columns <- function(comparison) {
  if (comparison == "KO_minus_WT") {
    c(paste0("HT1080_p21KOminusWT_", c("7GY", "UVA", "UVB")),
      paste0("MCF10A_p21KOminusWT_", c("7GY", "UVA", "UVB")))
  } else {
    unlist(lapply(c("HT1080", "MCF10A"), function(cell)
      unlist(lapply(c("WT", "p21KO"), function(genotype)
        paste(cell, genotype, c("7GY", "UVA", "UVB"), sep = "_")))))
  }
}

## Join by feature name, NEVER by export row position. Ambiguous selected
## names are an error: IPA does not supply a unique stable ID in these exports.
select_ipa <- function(pair, features, comparison) {
  if (anyDuplicated(features)) stop("Selected features must be unique")
  a <- pair$z
  b <- pair$p
  present <- features[features %in% a[[1L]]]
  not_present <- setdiff(features, present)
  if (length(not_present)) message("Not in this IPA export: ", paste(not_present, collapse = ", "))
  if (!length(present)) stop("No selected features found")
  if (anyDuplicated(a[[1L]][a[[1L]] %in% present]) ||
      anyDuplicated(b[[1L]][b[[1L]] %in% present]))
    stop("Ambiguous duplicate IPA feature name(s) among selected features")
  col_order <- ipa_columns(comparison)
  if (!setequal(col_order, names(a)[-1L])) stop("Unexpected contrast columns")
  z <- as.matrix(a[match(present, a[[1L]]), col_order, drop = FALSE])
  p <- as.matrix(b[match(present, b[[1L]]), col_order, drop = FALSE])
  storage.mode(z) <- "numeric"
  storage.mode(p) <- "numeric"
  rownames(z) <- rownames(p) <- present
  if (any(is.infinite(z)) || any(is.infinite(p))) stop("Infinite IPA values in input")
  list(z = z, neg_log10_bh = p)
}

## The same exact blue–white–red colors/breaks as module Figures 13/14.
cap <- 6
colors <- c("#12386F", "#306FB3", "#94BDE2", "#E0ECF8", "#FFFFFF",
            "#FAE1DD", "#E89A94", "#C94644", "#941E2C")
breaks <- cap * c(-1, -.75, -.40, -.12, 0, .12, .40, .75, 1)
col_fun <- circlize::colorRamp2(breaks, colors)

contrast_labels <- if (comparison == "KO_minus_WT") {
  paste(rep(c("HT1080", "MCF10A"), each = 3),
        rep(c("7 Gy", "UVA", "UVB"), 2), sep = " | ")
} else {
  paste(rep(c("HT1080", "MCF10A"), each = 6),
        rep(rep(c("WT", "KO"), each = 3), 2),
        rep(c("7 Gy", "UVA", "UVB"), 4), sep = " | ")
}

## 2 | Upstream regulators --------------------------------------------------------
## Deliberately broad, hypothesis-driven candidates; edit this vector to explore.
regulators <- c(
  "TP53", "CDKN1A", "CDKN2A", "RB1", "RBL1", "RBL2",
  "E2F", "E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "E2F6", "FOXM1", "MYC", "MYCN",
  "TGFB1", "TGFB2", "SMAD3", "SMAD4", "SMAD7", "STAT3", "STAT1", "CEBPB", "TNF",
  "VEGF", "YAP1", "HGF", "EGFR", "KRAS", "FOS", "JUN",
  "PPARD", "PPARGC1A", "TFEB", "ATF4", "NUPR1", "NFE2L2", "HSF1", "RICTOR", "DELE1", "XBP1"
)
upstream <- select_ipa(read_ipa_pair("upstream_regulators", comparison), regulators, comparison)
z <- upstream$z
log10_bh <- upstream$neg_log10_bh       # IPA stores -log10(BH-adjusted p)

## A regulator with no predicted z-score in ANY contrast cannot be clustered.
all_missing <- rowSums(is.finite(z)) == 0L
if (any(all_missing)) {
  message("No z-score in any contrast (omitting from heatmap): ",
          paste(rownames(z)[all_missing], collapse = ", "))
  z <- z[!all_missing, , drop = FALSE]
  log10_bh <- log10_bh[!all_missing, , drop = FALSE]
}

row_labels <- rownames(z)
row_labels[row_labels == "TP53"] <- "TP53 (p53)"
row_labels[row_labels == "CDKN1A"] <- "CDKN1A (p21)"
row_labels[row_labels == "CDKN2A"] <- "CDKN2A (p16/p14ARF locus)"

## Clustering ONLY: replace missing values by within-row means.
## Original NA predictions remain NA in 'z' and are shown in grey.
## All-NA rows were removed above; no biological zero is manufactured.
z_for_clustering <- z
for (i in seq_len(nrow(z_for_clustering))) {
  absent <- !is.finite(z_for_clustering[i, ])
  if (any(absent)) {
    z_for_clustering[i, absent] <- mean(z_for_clustering[i, !absent])
  }
}
row_dendrogram <- stats::as.dendrogram(stats::hclust(stats::dist(z_for_clustering)))

## Edit this Heatmap() call directly in RStudio:
##   row_split = 3    : cut the dendrogram into three groups (your preferred view)
##   row_split = NULL : no displayed split
##   cluster_rows = FALSE : preserve the manual 'regulators' ordering
##   column clustering/splits are disabled to preserve biological column order.
## row_km > 1 is not used: k-means on the original matrix fails with IPA NA values.
ht_upstream <- ComplexHeatmap::Heatmap(
  z,
  name = "IPA activation\nz-score",
  col = col_fun,
  na_col = "#E8E8E8",  # distinguish 'not predicted' (grey) from z = 0 (white)
  cluster_rows = row_dendrogram,
  cluster_columns = FALSE,
  row_km = 1,
  column_km = 1,
  row_split = 3,
  column_split = NULL,
  row_labels = row_labels,
  column_labels = contrast_labels,
  row_names_side = "right",
  row_names_gp = grid::gpar(fontsize = 8.5),
  column_names_gp = grid::gpar(fontsize = 9),
  column_names_rot = 45,
  column_title = "IPA upstream regulators \n p21 interaction (KO - WT)",
  column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
  rect_gp = grid::gpar(col = "#F2F2F2", lwd = 0.25),
  border = TRUE,
  heatmap_legend_param = list(at = c(-6, -3, 0, 3, 6))
)

## Viewport rendering: view figure in the RStudio Plots pane.
ComplexHeatmap::draw(ht_upstream, heatmap_legend_side = "right")

## 3 | IPA Diseases & Bio Functions ----------------------------------------------
## Functions were chosen A PRIORI for the exploratory figure as interpretable
## cell survival/death, cell-cycle, migration and invasion phenotypes.
## These are IPA predictions from the six FORMAL p21 interactions, not assays.
functions <- c(
  "Cell survival", "Cell viability", "Apoptosis", "Cell death of tumor cell lines",
  "Cell proliferation of tumor cell lines", "Cell cycle progression",
  "Migration of cells", "Invasion of cells"
)
function_data <- select_ipa(read_ipa_pair("functions", "KO_minus_WT"),
                            functions, "KO_minus_WT")
z_functions <- function_data$z
p_functions <- function_data$neg_log10_bh

## Cell labels show the original (unclipped) score; * means BH-adjusted p < 0.05.
## The color scale is clipped at +/- 6, consistent with the regulator heatmap.
function_cell <- function(j, i, x, y, width, height, fill) {
  value <- z_functions[i, j]
  bh <- p_functions[i, j]
  label <- if (is.finite(value)) sprintf("%+.1f", value) else "—"
  if (is.finite(value) && is.finite(bh) && bh > -log10(0.05)) label <- paste0(label, "*")
  txt_color <- if (is.finite(value) && abs(value) >= 3.8) "white" else "#233245"
  grid::grid.text(label, x, y, gp = grid::gpar(col = txt_color, fontsize = 9))
}

## Separate, editable ComplexHeatmap() call for the functions heatmap.
ht_functions <- ComplexHeatmap::Heatmap(
  z_functions,
  name = "IPA activation\nz-score",
  col = col_fun,
  na_col = "#E8E8E8",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_km = 1,
  column_km = 1,
  row_split = NULL,
  column_split = NULL,
  column_labels = paste(rep(c("HT1080", "MCF10A"), each = 3),
                        rep(c("7 Gy", "UVA", "UVB"), 2), sep = " | "),
  row_names_side = "left",
  row_names_gp = grid::gpar(fontsize = 9),
  column_names_gp = grid::gpar(fontsize = 9),
  column_names_rot = 45,
  column_title = "IPA biological functions \n p21 interaction (KO - WT)",
  column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
  rect_gp = grid::gpar(col = "#F2F2F2", lwd = 0.25),
  border = TRUE,
  cell_fun = function_cell,
  heatmap_legend_param = list(at = c(-6, -3, 0, 3, 6))
)
ComplexHeatmap::draw(ht_functions, heatmap_legend_side = "right")
message("* in the function heatmap: BH-adjusted p < 0.05; grey: no IPA z-score.")

## 4 | Optional PNG and editable SVG via ekbSeq -------------------------------
if (export_figures) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  suffix <- if (comparison == "KO_minus_WT") "interactions" else "primary_responses"
  ekbSeq::export_plot_dual(
    file.path(figure_dir, paste0("IPA_upstream_regulators_", suffix)),
    plot_expr = function() ComplexHeatmap::draw(ht_upstream, heatmap_legend_side = "right"),
    width = 6, height = 11, dpi = 300, rasterize = FALSE,
    grid_capture = TRUE, svg_backend = "svglite"
  )
  ekbSeq::export_plot_dual(
    file.path(figure_dir, "IPA_functions_interactions"),
    plot_expr = function() ComplexHeatmap::draw(ht_functions, heatmap_legend_side = "right"),
    width = 6, height = 6.5, dpi = 300, rasterize = FALSE,
    grid_capture = TRUE, svg_backend = "svglite"
  )
}
