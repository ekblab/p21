## ************************************************************************************************
## p21 bulk RNA-seq
## 04 - Formal p21-dependent irradiation responses (24 h)
##
## Biological question:
##   Does p21 loss alter the transcriptional response to 7 Gy, UVA or UVB
##   within HT1080 or MCF10A cells?
##
## Primary model (A): ~ processing + group
## Interaction (difference in differences):
##   (p21-negative irradiated - p21-negative control)
##     - (WT irradiated - WT control)
##
## Analysis decisions:
##   - Reuse the fitted 91-sample DESeq2 Model A; do NOT refit a new model.
##   - Six formal Wald contrasts: two cell lines x three irradiations.
##   - Positive interaction LFC = MORE POSITIVE treatment response after p21 loss;
##     negative interaction LFC = MORE NEGATIVE response. Neither implies
##     induction/repression without inspecting both genotype-specific responses.
##   - Significance: DESeq2 BH-adjusted padj < 0.05; NO log2FC cutoff.
##   - Unshrunken interaction LFC for Wald tests/volcano plots; ashr-shrunken
##     interaction LFC for biological interpretation (NOT a second test).
##   - Validate the interaction against fitted GLM coefficients AND the original
##     treatment contrasts (excluding zero-count-pair special cases in DESeq2).
##   - No module integration, GO/GSEA, IPA or DREAM analysis in this script.
##
## Output (all below Results/04_p21_interactions/):
##   tables/{unshrunken,shrunken}/: six annotated gene-level CSVs each
##   tables/interaction_summary.csv: six-row result overview
##   tables/interaction_contrast_definitions.csv: explicit four-group contrasts
##   tables/7Gy_UVA_interaction_concordance.csv: cross-radiation comparison
##   tables/selected_gene_responses.csv: curated WT/KO response examples
##   figures/: interaction counts, volcano atlas, 7 Gy/UVA concordance,
##             selected-gene WT/KO dot plots (all PNG/SVG)
##   The single Data/processed_data.rds is read-only; no additional RDS files.
## ************************************************************************************************

rm(list = ls())
gc()

## ******
## Project setup, dependencies and paths
if (!requireNamespace("here", quietly = TRUE)) stop("Install 'here'.")
here::i_am("R/04_p21_interactions.R")
source(here::here("R", "00_setup.R"))
source(here::here("R", "plot_style.R"))
for (pkg in c("DESeq2", "SummarizedExperiment", "ashr", "ggplot2", "patchwork",
              "EnhancedVolcano", "ekbSeq", "scales")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing R package: ", pkg)
}

alpha <- 0.05
root <- here::here("Results", "04_p21_interactions")
tab_raw_dir <- file.path(root, "tables", "unshrunken")
tab_shr_dir <- file.path(root, "tables", "shrunken")
fig_dir <- file.path(root, "figures")
primary_raw_dir <- here::here("Results", "02_primary_DE", "01_PRIMARY_model_A", "tables", "unshrunken")
primary_shr_dir <- here::here("Results", "02_primary_DE", "01_PRIMARY_model_A", "tables", "shrunken")
for (path in c(tab_raw_dir, tab_shr_dir, fig_dir))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

## ************************************************************************************************
## Load the fitted primary model and verify the analysis context
processed <- readRDS(here::here("Data", "processed_data.rds"))
dds <- processed$dds
annotation <- processed$annotation
meta <- as.data.frame(SummarizedExperiment::colData(dds))
stopifnot(ncol(dds) == 91L, nrow(dds) == 61806L,
          identical(rownames(dds), annotation$ENSEMBL), !anyDuplicated(annotation$ENSEMBL),
          identical(deparse(DESeq2::design(dds)), deparse(~ processing + group)),
          all(c("cell_line", "genotype", "radiation", "processing", "group") %in% names(meta)),
          length(DESeq2::resultsNames(dds)) == ncol(stats::model.matrix(DESeq2::design(dds), meta)))
check_design(meta, ~ processing + group)

## ******
## Define six biological comparisons and their four component groups
plan <- expand.grid(radiation = c("7gy", "UVA", "UVB"),
                    cell_line = c("HT1080", "MCF10A"), stringsAsFactors = FALSE)
plan$id <- with(plan, paste(cell_line, radiation, sep = "_"))
plan$KO_treated <- with(plan, mapply(group_name, cell_line, "p21_neg", radiation))
plan$KO_control <- with(plan, mapply(group_name, cell_line, "p21_neg", "control"))
plan$WT_treated <- with(plan, mapply(group_name, cell_line, "wt", radiation))
plan$WT_control <- with(plan, mapply(group_name, cell_line, "wt", "control"))
stopifnot(nrow(plan) == 6L, all(unlist(plan[c("KO_treated", "KO_control",
                                                  "WT_treated", "WT_control")]) %in% levels(meta$group)))
utils::write.csv(plan, file.path(root, "tables", "interaction_contrast_definitions.csv"), row.names = FALSE)

## ******
## Read the corresponding WT/KO response tables from script 02.
## These are ONLY for annotation and validation, never separate significance tests
## of the interaction. All interaction p-values come directly from DESeq2.
read_primary <- function(line, genotype, radiation, method = c("unshrunken", "shrunken")) {
  method <- match.arg(method)
  id <- paste(line, genotype, radiation, sep = "_")
  filename <- if (method == "unshrunken") paste0(id, "_vs_control_unshrunken.csv") else
    paste0(id, "_vs_control_ashr.csv")
  path <- file.path(if (method == "unshrunken") primary_raw_dir else primary_shr_dir, filename)
  if (!file.exists(path)) stop("Missing primary Model-A CSV: ", path)
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("ENSEMBL", "log2FoldChange", "padj") %in% names(x)) ||
      anyNA(x$ENSEMBL) || anyDuplicated(x$ENSEMBL)) stop("Invalid primary CSV: ", path)
  expected <- if (method == "unshrunken") "Wald_unshrunken" else "ashr_shrunken"
  if ("lfc_method" %in% names(x) && any(x$lfc_method != expected, na.rm = TRUE))
    stop("Unexpected LFC estimator in: ", path)
  x
}

## Annotate a formal interaction with the separate treatment responses.
## NOTE: an ashr-shrunken interaction is not the subtraction of two independently
## shrunken treatment LFCs. Only the *unshrunken* estimates obey this identity.
annotate_interaction <- function(res, wt, ko, method, tested_significance = NULL) {
  tab <- result_table(res, annotation)
  get <- function(src, field) src[[field]][match(tab$ENSEMBL, src$ENSEMBL)]
  if (anyNA(match(tab$ENSEMBL, wt$ENSEMBL)) || anyNA(match(tab$ENSEMBL, ko$ENSEMBL)))
    stop("Response tables and interaction results have incompatible gene IDs.")
  tab$WT_response_log2FC <- get(wt, "log2FoldChange")
  tab$p21neg_response_log2FC <- get(ko, "log2FoldChange")
  tab$WT_response_padj <- get(wt, "padj")
  tab$p21neg_response_padj <- get(ko, "padj")
  tab$lfc_method <- method
  tab$significant <- if (is.null(tested_significance)) !is.na(tab$padj) & tab$padj < alpha else tested_significance
  ## Direction refers to the *unshrunken interaction*, not to expression level.
  tab$interaction_direction <- ifelse(tab$significant & tab$log2FoldChange > 0,
    "more_positive_in_p21neg", ifelse(tab$significant & tab$log2FoldChange < 0,
    "more_negative_in_p21neg", "not_significant"))
  tab[order(tab$padj, na.last = TRUE), , drop = FALSE]
}

## ************************************************************************************************
## Formal interaction contrasts from the existing DESeq2 Wald fit
##
## Build contrast coefficients through the *model matrix*, not by mixing
## observed processing distributions. The +1/-1 group weights sum to zero;
## processing cancels exactly under the additive primary model.
## Filtering below applies ONLY to QC of the exported response estimates;
## it does NOT filter tested genes or change the formal interaction test.
expressed <- rowSums(DESeq2::counts(dds, normalized = TRUE) >= 10) >= 5
raw_counts <- DESeq2::counts(dds, normalized = FALSE)
model_coef <- stats::coef(dds)  # fitted log2-scale coefficients, not shrunken
stopifnot(identical(rownames(model_coef), rownames(dds)),
          identical(colnames(model_coef), DESeq2::resultsNames(dds)))
summary_rows <- vector("list", nrow(plan))
volcano_panels <- vector("list", nrow(plan))
## Keep only columns needed for the two additional figures (no re-reading six
## large CSV files and no intermediate RDS objects).
viz_tables <- list()

for (i in seq_len(nrow(plan))) {
  p <- plan[i, ]
  message(sprintf("Interaction %d/%d: %s (p21-negative response minus WT response)",
                  i, nrow(plan), p$id))
  weights <- c(setNames( 1, p$KO_treated), setNames(-1, p$KO_control),
               setNames(-1, p$WT_treated), setNames( 1, p$WT_control))
  contrast <- group_contrast(dds, weights)  # defined in R/project_helpers.R
  stopifnot(length(contrast) == length(DESeq2::resultsNames(dds)),
            isTRUE(all.equal(unname(contrast[1]), 0)),
            identical(names(contrast), DESeq2::resultsNames(dds)))

  ## One joint Wald test of the interaction (includes the coefficient covariance).
  ## Comparing 'significant in KO' vs 'not significant in WT' is NOT this test.
  res <- DESeq2::results(dds, contrast = unname(contrast), alpha = alpha)
  shrunk <- DESeq2::lfcShrink(dds, res = res, type = "ashr", quiet = TRUE)
  stopifnot(identical(rownames(res), rownames(shrunk)),
            isTRUE(all.equal(res$pvalue, shrunk$pvalue, check.attributes = FALSE)),
            isTRUE(all.equal(res$padj, shrunk$padj, check.attributes = FALSE)))

  wt_raw <- read_primary(p$cell_line, "wt", p$radiation, "unshrunken")
  ko_raw <- read_primary(p$cell_line, "p21_neg", p$radiation, "unshrunken")
  wt_shr <- read_primary(p$cell_line, "wt", p$radiation, "shrunken")
  ko_shr <- read_primary(p$cell_line, "p21_neg", p$radiation, "shrunken")

  ## **** test 1: verify the four-group contrast against the fitted GLM.
  ## No reference CSV is involved in this check. Every supported estimate
  ## must equal c' %*% beta, using the ORIGINAL unshrunken coefficients.
  expected_model <- as.vector(model_coef %*% unname(contrast))

  ## **** test 2: validate against script 02's exported treatment contrasts.
  ## DESeq2::results(contrast = c("group", treated, control)) forces LFC=0 if
  ## BOTH groups have zero counts for a gene, even if it is expressed elsewhere.
  ## The four-group interaction need not have the same zeroing behavior. Thus
  ## the identity KO_response - WT_response is checked only where BOTH pairs
  ## (KO treated/control AND WT treated/control) contain at least one count.
  ## The global expression filter alone is not sufficient for this exception.
  idx_wt <- match(rownames(res), wt_raw$ENSEMBL)
  idx_ko <- match(rownames(res), ko_raw$ENSEMBL)
  if (anyNA(idx_wt) || anyNA(idx_ko)) stop("Missing reference genes for ", p$id)
  expected_csv <- ko_raw$log2FoldChange[idx_ko] - wt_raw$log2FoldChange[idx_wt]
  ko_pair <- meta$group %in% c(p$KO_treated, p$KO_control)
  wt_pair <- meta$group %in% c(p$WT_treated, p$WT_control)
  valid_pair <- rowSums(raw_counts[, ko_pair, drop = FALSE]) > 0 &
                rowSums(raw_counts[, wt_pair, drop = FALSE]) > 0
  comparable <- expressed & valid_pair & is.finite(res$log2FoldChange) &
                is.finite(expected_csv) & is.finite(expected_model)
  if (sum(comparable) < 1000L)
    stop("Too few well-expressed genes with supported WT/KO comparisons: ", p$id)

  error_model <- abs(res$log2FoldChange - expected_model)
  error_csv <- abs(res$log2FoldChange - expected_csv)
  max_model_error <- max(error_model[comparable])
  max_error <- max(error_csv[comparable])
  n_zero_pair <- sum(expressed & !valid_pair)
  message(sprintf("  LFC check: %d genes tested; %d expressed genes lack counts in a WT/KO pair; GLM max error %.3g; primary CSV max error %.3g",
                  sum(comparable), n_zero_pair, max_model_error, max_error))

  ## Do NOT suppress a genuine discrepancy. Report an example so mismatched
  ## primary CSVs, coefficient ordering or stale fitted data can be diagnosed.
  if (max_model_error > 1e-5 || max_error > 1e-5) {
    worst <- which.max(ifelse(comparable, pmax(error_model, error_csv), -Inf))
    stop(sprintf(paste0("Inconsistent interaction contrast for %s: GLM max %.4g, ",
                        "primary CSV max %.4g; example %s (interaction %.4f, ",
                        "model %.4f, KO-WT CSV %.4f). Check that scripts 02/04 ",
                        "use the same processed_data.rds and Model-A results."),
                 p$id, max_model_error, max_error, rownames(res)[worst],
                 res$log2FoldChange[worst], expected_model[worst], expected_csv[worst]))
  }

  raw_tab <- annotate_interaction(res, wt_raw, ko_raw, "Wald_unshrunken")
  shr_tab <- annotate_interaction(shrunk, wt_shr, ko_shr, "ashr_shrunken")
  if (p$radiation %in% c("7gy", "UVA")) {
    viz_tables[[p$id]] <- shr_tab[, c("ENSEMBL", "SYMBOL", "baseMean", "padj",
      "log2FoldChange", "WT_response_log2FC", "p21neg_response_log2FC")]
  }
  ## Preserve the original Wald direction in both exports: shrinkage does not
  ## establish independent significance or define the direction of the test.
  shr_tab$interaction_direction <- raw_tab$interaction_direction[match(shr_tab$ENSEMBL, raw_tab$ENSEMBL)]
  stopifnot(identical(raw_tab$ENSEMBL, shr_tab$ENSEMBL),
            identical(raw_tab$significant, shr_tab$significant))
  utils::write.csv(raw_tab, file.path(tab_raw_dir, paste0(p$id, "_p21_interaction_unshrunken.csv")), row.names = FALSE)
  utils::write.csv(shr_tab, file.path(tab_shr_dir, paste0(p$id, "_p21_interaction_ashr.csv")), row.names = FALSE)

  ## One summary row per biological contrast. Direction of an interaction is
  ## interpreted relative to WT, not as up- vs downregulated expression.
  n_pos <- sum(raw_tab$interaction_direction == "more_positive_in_p21neg")
  n_neg <- sum(raw_tab$interaction_direction == "more_negative_in_p21neg")
  summary_rows[[i]] <- data.frame(cell_line = p$cell_line, radiation = p$radiation,
    contrast_id = p$id, n_genes = nrow(raw_tab), n_pvalue = sum(!is.na(raw_tab$pvalue)),
    n_padj = sum(!is.na(raw_tab$padj)), n_interaction_FDR05 = n_pos + n_neg,
    n_more_positive_p21neg = n_pos, n_more_negative_p21neg = n_neg,
    n_LFC_checked = sum(comparable), n_expressed_zero_pair = n_zero_pair,
    max_abs_LFC_model_error = max_model_error, max_abs_LFC_check_error = max_error)

  ## ******
  ## Six-panel volcano atlas (the original unshrunken interaction Wald effects)
  ## Axes use each panel's actual data range: no percentile clipping/sampling.
  v <- raw_tab[is.finite(raw_tab$log2FoldChange) & !is.na(raw_tab$padj), , drop = FALSE]
  if (!nrow(v)) stop("No plottable interactions: ", p$id)
  missing_symbol <- is.na(v$SYMBOL) | !nzchar(v$SYMBOL)
  v$SYMBOL[missing_symbol] <- v$ENSEMBL[missing_symbol]
  if (any(v$padj < 0)) stop("Negative adjusted p-value in ", p$id)
  if (any(v$padj == 0)) {
    message(p$id, ": replacing ", sum(v$padj == 0), " underflowed padj in plot ONLY.")
    v$padj[v$padj == 0] <- .Machine$double.xmin
  }
  x_range <- range(v$log2FoldChange)
  x_margin <- max(0.2, diff(x_range) * 0.055)
  y_max <- max(-log10(v$padj))
  volcano_color <- ifelse(v$interaction_direction == "more_positive_in_p21neg",
    npg_colors[["red"]], ifelse(v$interaction_direction == "more_negative_in_p21neg",
    npg_colors[["navy"]], "#B8B8B8"))
  names(volcano_color) <- ifelse(v$interaction_direction == "more_positive_in_p21neg",
    "Positive interaction", ifelse(v$interaction_direction == "more_negative_in_p21neg",
    "Negative interaction", "Not significant"))
  irradiation_label <- if (p$radiation == "7gy") "7 Gy" else p$radiation
  volcano_panels[[i]] <- ekbSeq::plot_volcano(
    v, nlabel = 10, title = paste(p$cell_line, irradiation_label, sep = " | "),
    pointsize = 0.85, labSize = 3.0, plot.signif = TRUE, pThres = alpha, lfcThres = 0,
    xlim = x_range + c(-x_margin, x_margin), ylim = c(-0.2, y_max * 1.16 + 0.75),
    colCustom = volcano_color, colAlpha = 0.60, caption = "", subtitle = "", legendPosition = "none",
    axisLabSize = 9, titleLabSize = 10, drawConnectors = TRUE,
    widthConnectors = 0.25, max.overlaps = 10, maxoverlapsConnectors = 10,
    cutoffLineType = "dashed", cutoffLineCol = "grey65", cutoffLineWidth = 0.35,
    xlab = expression(Delta*log[2]~"response (p21-negative - WT)"),
    ylab = expression(-log[10]~"adjusted p-value"),
    border = "full", borderWidth = 0.5, borderColour = "black"
  ) + p21_theme(9) + ggplot2::theme(
    panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.45),
    axis.line = ggplot2::element_blank(), legend.position = "none",
    plot.title = ggplot2::element_text(face = "bold", size = 10),
    plot.margin = ggplot2::margin(7, 9, 7, 9))
}

## ************************************************************************************************
## Six-contrast summary and interpretation guardrails
summary_interaction <- do.call(rbind, summary_rows)
stopifnot(nrow(summary_interaction) == 6L,
          all(summary_interaction$n_interaction_FDR05 ==
              summary_interaction$n_more_positive_p21neg + summary_interaction$n_more_negative_p21neg))
utils::write.csv(summary_interaction, file.path(root, "tables", "interaction_summary.csv"), row.names = FALSE)
print(summary_interaction[, c("cell_line", "radiation", "n_interaction_FDR05",
                              "n_more_positive_p21neg", "n_more_negative_p21neg")], row.names = FALSE)

## ************************************************************************************************
## Figure 1: stacked count of genes with formal p21 x irradiation interactions
## Stacks report the sign of the *difference in treatment effects*.
bar_data <- rbind(
  data.frame(summary_interaction, direction = "More negative in p21−",
             count = summary_interaction$n_more_negative_p21neg),
  data.frame(summary_interaction, direction = "More positive in p21−",
             count = summary_interaction$n_more_positive_p21neg))
bar_data$direction <- factor(bar_data$direction,
                             levels = c("More negative in p21−", "More positive in p21−"))
bar_data$x <- match(bar_data$contrast_id, summary_interaction$contrast_id)
bar_data$x[bar_data$cell_line == "MCF10A"] <- bar_data$x[bar_data$cell_line == "MCF10A"] + 0.6
bar_width <- 0.70
max_n <- max(summary_interaction$n_interaction_FDR05)
y_top <- max(5, ceiling(max_n * 1.22 / 5) * 5)
locations <- unique(bar_data[, c("contrast_id", "cell_line", "radiation", "x")])
locations <- locations[match(summary_interaction$contrast_id, locations$contrast_id), ]
locations$total <- summary_interaction$n_interaction_FDR05

p_counts <- ggplot2::ggplot(bar_data,
    ggplot2::aes(x = x, y = count, fill = direction)) +
  ggplot2::geom_col(width = bar_width, color = "black", linewidth = 0.34) +
  ggplot2::geom_text(data = locations,
     ggplot2::aes(x = x, y = total, label = scales::comma(total)),
     inherit.aes = FALSE, vjust = -0.55, fontface = "bold", size = 3.35) +
  ggplot2::geom_text(data = locations,
     ggplot2::aes(x = x, y = -0.045 * y_top,
                  label = ifelse(radiation == "7gy", "7 Gy", radiation)),
     inherit.aes = FALSE, size = 3.55) +
  ggplot2::annotate("segment", x = 0.65, xend = 3.35, y = -0.115*y_top, yend = -0.115*y_top) +
  ggplot2::annotate("segment", x = 4.25, xend = 6.95, y = -0.115*y_top, yend = -0.115*y_top) +
  ggplot2::annotate("text", x = 2, y = -0.15*y_top, label = "HT1080", size = 4, fontface = "bold") +
  ggplot2::annotate("text", x = 5.6, y = -0.15*y_top, label = "MCF10A", size = 4, fontface = "bold") +
  ggplot2::scale_fill_manual(values = c("More negative in p21−" = npg_colors[["navy"]],
                                          "More positive in p21−" = npg_colors[["red"]]), name = NULL) +
  ggplot2::scale_x_continuous(limits = c(0.55, 7.05), breaks = NULL,
                               expand = ggplot2::expansion(mult = 0)) +
  ggplot2::scale_y_continuous(labels = scales::label_comma(),
                               breaks = scales::breaks_pretty(n = 6),
                               expand = ggplot2::expansion(mult = 0)) +
  ggplot2::coord_cartesian(ylim = c(0, y_top), clip = "off", expand = FALSE) +
  ggplot2::labs(title = "p21-dependent transcriptional responses to irradiation",
    subtitle = "DESeq2 Model A | difference-in-differences Wald tests | BH-adjusted p < 0.05\nRed/blue compare signed irradiation responses (KO - WT), not absolute response strength",
    caption = NULL,  # Interpretation is specified in the subtitle; avoids overlap with grouped x-axis labels.
    x = NULL, y = "Genes with significant interaction") +
  p21_theme(11) + ggplot2::theme(
    panel.border = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.65),
    axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.65),
    axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
    legend.position = "top", legend.justification = "left",
    plot.title = ggplot2::element_text(face = "bold", size = 16),
    plot.subtitle = ggplot2::element_text(size = 10.5, color = "grey30"),
    plot.caption = ggplot2::element_text(size = 9, color = "grey35"),
    plot.margin = ggplot2::margin(t = 7, r = 15, b = 95, l = 7))
ekbSeq::export_plot_dual(file.path(fig_dir, "p21_interaction_counts"),
                         p_counts, width = 12.5, height = 7.0, rasterize = FALSE)

## ************************************************************************************************
## Figure 2: six independent volcano panels (two cell lines x three treatments)
## Every panel has independent real axis limits; no LFC threshold or clipping.
## Use a real, shared key below the panels. Sign refers to *the difference
## in response LFC*, never necessarily the observed absolute expression.
interaction_key <- data.frame(
  x = c(1, 5, 9), y = 1,
  label = c("Higher response LFC in p21 KO", "Higher response LFC in WT", "Not significant"),
  color = c(npg_colors[["red"]], npg_colors[["navy"]], "#B8B8B8"))
p_key <- ggplot2::ggplot(interaction_key, ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_point(ggplot2::aes(color = color), size = 3) +
  ggplot2::geom_text(ggplot2::aes(x = x + 0.18, label = label), hjust = 0, size = 3.05) +
  ggplot2::scale_color_identity() +
  ggplot2::scale_x_continuous(limits = c(0.6, 13), expand = c(0, 0)) +
  ggplot2::coord_cartesian(ylim = c(0.65, 1.35), clip = "off") +
  ggplot2::theme_void() +
  ggplot2::theme(plot.margin = ggplot2::margin(0, 5, 0, 5))

volcano_atlas <- (patchwork::wrap_plots(volcano_panels, ncol = 3) / p_key) +
  patchwork::plot_layout(heights = c(21, 1)) +
  patchwork::plot_annotation(
    title = "p21-dependent irradiation responses: formal interactions",
    subtitle = "KO response - WT response | Wald interaction LFC | BH-adjusted p < 0.05",
    caption = "Red: more positive response in KO; blue: more positive response in WT. This need not mean absolute expression is higher.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16),
      plot.subtitle = ggplot2::element_text(size = 10.5, color = "grey30"),
      plot.caption = ggplot2::element_text(size = 9, color = "grey35")))
ekbSeq::export_plot_dual(file.path(fig_dir, "p21_interaction_volcano_atlas"),
                         volcano_atlas, width = 18, height = 13.4, rasterize = FALSE)

## ************************************************************************************************
## Figure 3: concordance of p21 interactions after 7 Gy and UVA
##
## These are *directly ashr-shrunken interaction* estimates, NOT differences
## of separately shrunken treatment responses. Correlations are descriptive:
## even strongly concordant genes do not by themselves prove shared regulators.
## All genes with finite estimates and nonzero mean count are shown. Red points
## require a formal BH-significant interaction in BOTH irradiation conditions.
concordance_panels <- vector("list", 2L)
concordance_rows <- vector("list", 2L)

for (j in seq_along(c("HT1080", "MCF10A"))) {
  cell <- c("HT1080", "MCF10A")[[j]]
  a <- viz_tables[[paste0(cell, "_7gy")]]
  b <- viz_tables[[paste0(cell, "_UVA")]]
  idx <- match(a$ENSEMBL, b$ENSEMBL)
  if (anyNA(idx) || anyDuplicated(a$ENSEMBL)) stop("Mismatched interaction genes: ", cell)

  paired <- data.frame(ENSEMBL = a$ENSEMBL, baseMean = a$baseMean,
                       Gy7 = a$log2FoldChange, UVA = b$log2FoldChange[idx],
                       sig_7gy = !is.na(a$padj) & a$padj < alpha,
                       sig_UVA = !is.na(b$padj[idx]) & b$padj[idx] < alpha)
  valid <- is.finite(paired$Gy7) & is.finite(paired$UVA) &
           is.finite(paired$baseMean) & paired$baseMean > 0
  paired <- paired[valid, , drop = FALSE]
  paired$status <- ifelse(paired$sig_7gy & paired$sig_UVA, "Both conditions",
                          ifelse(paired$sig_7gy | paired$sig_UVA, "One condition", "Neither"))
  shared <- paired$status == "Both conditions"
  same_sign <- sum(sign(paired$Gy7[shared]) == sign(paired$UVA[shared]))
  r <- if (nrow(paired) > 2L) stats::cor(paired$Gy7, paired$UVA) else NA_real_
  concordance_rows[[j]] <- data.frame(cell_line = cell, n_genes = nrow(paired),
    pearson_r = r, n_significant_7gy = sum(paired$sig_7gy),
    n_significant_UVA = sum(paired$sig_UVA), n_shared_FDR05 = sum(shared),
    n_shared_same_sign = same_sign,
    pct_shared_same_sign = if (sum(shared)) 100 * same_sign / sum(shared) else NA_real_)

  ## Identical x/y ranges within each cell-line panel allow comparison with the
  ## identity line. Ranges are chosen from the REAL extrema (no quantile caps).
  xy_max <- max(abs(c(paired$Gy7, paired$UVA))) * 1.07
  scatter_colors <- c("Neither" = "#B8B8B8", "One condition" = npg_colors[["navy"]],
                      "Both conditions" = npg_colors[["red"]])
  paired$status <- factor(paired$status, levels = names(scatter_colors))
  concordance_panels[[j]] <- ggplot2::ggplot(paired,
      ggplot2::aes(x = Gy7, y = UVA)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "grey50",
                         linewidth = 0.45, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, color = "grey75") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, color = "grey75") +
    ggplot2::geom_point(data = paired[paired$status == "Neither", , drop = FALSE],
                        ggplot2::aes(color = status), size = 0.45, alpha = 0.16) +
    ggplot2::geom_point(data = paired[paired$status == "One condition", , drop = FALSE],
                        ggplot2::aes(color = status), size = 0.75, alpha = 0.55) +
    ggplot2::geom_point(data = paired[paired$status == "Both conditions", , drop = FALSE],
                        ggplot2::aes(color = status), size = 1.12, alpha = 0.92) +
    ggplot2::scale_color_manual(values = scatter_colors, name = "Significant interaction") +
    ggplot2::scale_x_continuous(limits = c(-xy_max, xy_max)) +
    ggplot2::scale_y_continuous(limits = c(-xy_max, xy_max)) +
    ggplot2::coord_fixed(ratio = 1, expand = FALSE) +
    ggplot2::labs(title = cell,
      subtitle = sprintf("Pearson r = %.3f  |  Shared FDR < 0.05: %s  |  Same sign: %s/%s",
        r, scales::comma(sum(shared)), same_sign, sum(shared)),
      x = "7 Gy interaction (ashr log2FC)", y = "UVA interaction (ashr log2FC)") +
    p21_theme(10) + ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey30"),
      legend.position = "bottom", legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 9))
  message(sprintf("7 Gy / UVA concordance: %s | r = %.3f | %d shared FDR < 0.05",
                  cell, r, sum(shared)))
}

concordance_summary <- do.call(rbind, concordance_rows)
utils::write.csv(concordance_summary, file.path(root, "tables", "7Gy_UVA_interaction_concordance.csv"),
                 row.names = FALSE)
p_concordance <- patchwork::wrap_plots(concordance_panels, ncol = 2, guides = "collect") +
  patchwork::plot_annotation(
    title = "Concordance of p21-dependent irradiation responses",
    subtitle = "Comparison of directly ashr-shrunken interactions after 7 Gy and UVA",
    caption = "Each point is a gene with nonzero mean expression. Red: BH FDR < 0.05 in both contrasts; blue: one. Dashed line: y = x.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey30"),
      plot.caption = ggplot2::element_text(size = 9, color = "grey35")))
ekbSeq::export_plot_dual(file.path(fig_dir, "p21_7Gy_UVA_interaction_concordance"),
                         p_concordance, width = 13.8, height = 7.6, rasterize = FALSE)

## ************************************************************************************************
## Figure 4: selected-gene WT vs p21-negative irradiation responses
##
## The points show the SEPARATELY ashr-shrunken treatment/control log2FCs
## from script 02. Connecting lines highlight the difference in response;
## statistical evidence comes from the DIRECT four-group Wald interaction.
## Candidates are examples chosen a priori for illustration, not a new gene
## selection/enrichment test. Export their values and interaction FDR.
examples <- list(
  HT1080_UVA = c("MCM2", "MCM7", "KIF11", "KNL1", "FOXM1"),
  MCF10A_7gy = c("IL6", "INHBA", "CCN2", "AXL", "SLC7A11"),
  MCF10A_UVA = c("IL6", "INHBA", "CCN2", "AXL", "SLC7A11"))
selected_rows <- list()

for (id in names(examples)) {
  tab <- viz_tables[[id]]
  selected <- lapply(examples[[id]], function(symbol) {
    hit <- which(tab$SYMBOL == symbol)
    if (!length(hit)) { warning(id, ": missing candidate ", symbol); return(NULL) }
    if (length(hit) > 1L) {
      hit <- hit[which.max(tab$baseMean[hit])]
      warning(id, ": multiple features for ", symbol, "; selected highest baseMean")
    }
    tab[hit, , drop = FALSE]
  })
  selected <- do.call(rbind, selected)
  if (is.null(selected) || !nrow(selected)) stop("No curated genes found: ", id)
  selected$contrast_id <- id
  selected$gene_order <- match(selected$SYMBOL, examples[[id]])
  selected_rows[[id]] <- selected
}
selected_genes <- do.call(rbind, selected_rows)
selected_genes$interaction_FDR05 <- !is.na(selected_genes$padj) & selected_genes$padj < alpha
utils::write.csv(selected_genes[, c("contrast_id", "ENSEMBL", "SYMBOL", "gene_order",
  "WT_response_log2FC", "p21neg_response_log2FC", "log2FoldChange", "padj", "interaction_FDR05")],
  file.path(root, "tables", "selected_gene_responses.csv"), row.names = FALSE)

## Keep a common x scale for both MCF10A panels; use a tighter HT1080 scale
## so its smaller log2FC differences remain legible. No points are clipped.
mcf_lim <- max(abs(c(selected_rows$MCF10A_7gy$WT_response_log2FC,
  selected_rows$MCF10A_7gy$p21neg_response_log2FC,
  selected_rows$MCF10A_UVA$WT_response_log2FC,
  selected_rows$MCF10A_UVA$p21neg_response_log2FC)), na.rm = TRUE) * 1.12
ht_lim <- max(abs(c(selected_rows$HT1080_UVA$WT_response_log2FC,
  selected_rows$HT1080_UVA$p21neg_response_log2FC)), na.rm = TRUE) * 1.16

dot_panels <- lapply(names(examples), function(id) {
  tab <- selected_rows[[id]]
  tab$gene <- factor(tab$SYMBOL, levels = rev(examples[[id]]))
  dots <- rbind(data.frame(gene = tab$gene, genotype = "WT",
                           log2FC = tab$WT_response_log2FC),
                data.frame(gene = tab$gene, genotype = "p21-KO",
                           log2FC = tab$p21neg_response_log2FC))
  dots$genotype <- factor(dots$genotype, levels = c("WT", "p21-KO"))
  caption <- if (id == "HT1080_UVA") "HT1080 | UVA" else
    paste0("MCF10A | ", if (grepl("7gy$", id)) "7 Gy" else "UVA")
  bound <- if (id == "HT1080_UVA") ht_lim else mcf_lim
  ggplot2::ggplot(tab, ggplot2::aes(y = gene)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.4) +
    ggplot2::geom_segment(ggplot2::aes(x = WT_response_log2FC,
      xend = p21neg_response_log2FC, yend = gene), linewidth = 0.8, color = "grey55") +
    ggplot2::geom_point(data = dots, ggplot2::aes(x = log2FC, color = genotype),
      inherit.aes = TRUE, size = 2.6) +
    ggplot2::scale_color_manual(values = c("WT" = npg_colors[["navy"]],
                                           "p21-KO" = npg_colors[["red"]]), name = NULL) +
    ggplot2::scale_x_continuous(limits = c(-bound, bound)) +
    ggplot2::labs(title = caption, x = "Irradiation vs control (ashr log2FC)", y = NULL) +
    p21_theme(10) + ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      axis.text.y = ggplot2::element_text(face = "italic", size = 10),
      legend.position = "bottom")
})
p_dots <- patchwork::wrap_plots(dot_panels, ncol = 3, guides = "collect") +
  patchwork::plot_annotation(
    title = "Examples of p21-modulated irradiation responses",
    subtitle = "Separate WT and KO treatment effects (ashr log2FC); lines join the same gene",
    caption = "Illustrative genes selected for biological relevance. Formal interaction FDRs are in tables/selected_gene_responses.csv; point differences are not ashr interaction estimates.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey30"),
      plot.caption = ggplot2::element_text(size = 9, color = "grey35")))
ekbSeq::export_plot_dual(file.path(fig_dir, "p21_selected_genes_WT_vs_KO"),
                         p_dots, width = 15.5, height = 6.4, rasterize = FALSE)

## ************************************************************************************************
## Summary and downstream boundaries
## Script 05 will integrate THESE six tested effects with the existing L2
## modules. Lack of a significant interaction does NOT establish p21 independence.
## BH FDR is controlled within each contrast, not jointly across all six contrasts.
message("Interaction results: ", root)
message("Done: 6 direct Wald interactions; 6 ashr tables; 4 figures (counts, volcano, scatter, dot plot).")