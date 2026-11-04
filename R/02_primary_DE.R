## ************************************************************************************************
## p21 bulk RNA-seq
## 02 - Primary differential expression after irradiation (24 h)
##
## Biological question:
##   Which genes respond to 7 Gy, UVA and UVB in HT1080 and MCF10A,
##   separately for WT and p21-negative cells?
##
## Primary model (A):  ~ processing + group
## Sensitivity model (B): ~ processing + sequencing_run + group
##
## Analysis decisions:
##   - Use the 91 samples and endogenous counts from Data/processed_data.rds.
##   - Compare irradiation with the matched control within cell line/genotype.
##   - DEGs: DESeq2 BH-adjusted padj < 0.05; NO mandatory log2FC cutoff.
##   - Wald log2FC (unshrunken): DE tests and volcano plots.
##   - ashr-shrunken log2FC: interpretation and downstream response analyses.
##   - CSV name and lfc_method column explicitly identify the LFC estimator.
##   - Only Model A generates biological result tables and figures.
##   - Model B is a compact, separate robustness check, not a second DEG set.
##   - One processed RDS already exists; all results are CSV and PNG/SVG.
##
## Output (all below Results/02_primary_DE/):
##   01_PRIMARY_model_A/tables/{unshrunken,shrunken}/: 12 labeled CSVs each
##   01_PRIMARY_model_A/tables/DEG_summary.csv: one 12-row overview
##   01_PRIMARY_model_A/figures/: one-panel grouped/stacked DEG graph and volcano atlas
##   02_SENSITIVITY_model_B/: one comparison CSV and one-panel A/B diagnostic graph
## ************************************************************************************************

rm(list = ls())
gc()

## ******
## Project setup, dependencies and paths
if (!requireNamespace("here", quietly = TRUE)) stop("Install the 'here' package.")
here::i_am("R/02_primary_DE.R")
source(here::here("R", "00_setup.R"))
source(here::here("R", "plot_style.R"))
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install 'ggplot2'.")
if (!requireNamespace("patchwork", quietly = TRUE)) stop("Install 'patchwork' for the volcano atlas.")
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) stop("Install 'EnhancedVolcano' for ekbSeq::plot_volcano().")

alpha <- 0.05
primary_dir <- here::here("Results", "02_primary_DE", "01_PRIMARY_model_A")
sensitivity_dir <- here::here("Results", "02_primary_DE", "02_SENSITIVITY_model_B")
raw_dir <- file.path(primary_dir, "tables", "unshrunken")
shrunk_dir <- file.path(primary_dir, "tables", "shrunken")
fig_dir <- file.path(primary_dir, "figures")
sens_tab_dir <- file.path(sensitivity_dir, "tables")
sens_fig_dir <- file.path(sensitivity_dir, "figures")
for (path in c(raw_dir, shrunk_dir, fig_dir, sens_tab_dir, sens_fig_dir))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

## ************************************************************************************************
## Load the single processed object and check the primary model
processed <- readRDS(here::here("Data", "processed_data.rds"))
dds <- processed$dds
annotation <- processed$annotation
meta <- as.data.frame(SummarizedExperiment::colData(dds))

stopifnot(ncol(dds) == 91L, nrow(dds) == 61806L)
stopifnot(identical(rownames(dds), annotation$ENSEMBL), !anyDuplicated(annotation$ENSEMBL))
stopifnot(all(c("cell_line", "genotype", "radiation", "processing", "sequencing_run", "group") %in% names(meta)))
stopifnot(identical(deparse(DESeq2::design(dds)), deparse(~ processing + group)))
stopifnot(length(DESeq2::resultsNames(dds)) > 0L)
check_design(meta, ~ processing + group)
check_design(meta, ~ processing + sequencing_run + group)

## ************************************************************************************************
## Define the twelve irradiation-versus-control comparisons
## Group contrasts hold processing fixed; the processing term cancels.
## No unadjusted comparison of observed EKB vs UKD group means is used.
contrast_plan <- expand.grid(radiation = c("7gy", "UVA", "UVB"),
                             genotype = c("wt", "p21_neg"),
                             cell_line = c("HT1080", "MCF10A"), stringsAsFactors = FALSE)
contrast_plan$treatment <- with(contrast_plan, group_name(cell_line, genotype, radiation))
contrast_plan$control <- with(contrast_plan, group_name(cell_line, genotype, "control"))
contrast_plan$contrast_id <- with(contrast_plan, paste(cell_line, genotype, radiation, sep = "_"))
stopifnot(nrow(contrast_plan) == 12L,
          all(c(contrast_plan$treatment, contrast_plan$control) %in% levels(meta$group)))

## ******
## Small helpers used only by this script
annotate_de <- function(res, annotation) {
  tab <- result_table(res, annotation)
  tab$significant <- !is.na(tab$padj) & tab$padj < alpha
  tab$direction <- ifelse(tab$significant & tab$log2FoldChange > 0, "up",
                          ifelse(tab$significant & tab$log2FoldChange < 0, "down", "not significant"))
  tab[order(tab$padj, na.last = TRUE), , drop = FALSE]
}

deg_counts <- function(tab) {
  c(n_total = nrow(tab), n_pvalue = sum(!is.na(tab$pvalue)),
    n_padj = sum(!is.na(tab$padj)), n_DEG = sum(tab$significant),
    n_up = sum(tab$direction == "up"), n_down = sum(tab$direction == "down"))
}

## ************************************************************************************************
## Model A: primary DESeq2 results and ashr shrinkage
##
## DESeq2::results() supplies unshrunken Wald LFC, SE and the BH-adjusted padj.
## lfcShrink(type = "ashr", res = ...) changes the LFC/SE only; the original
## Wald p-values and padj remain the significance basis. We never re-test a
## shrunken LFC or impose an arbitrary magnitude cutoff for DEG classification.
summary_rows <- vector("list", nrow(contrast_plan))
volcano_panels <- vector("list", nrow(contrast_plan))

for (i in seq_len(nrow(contrast_plan))) {
  info <- contrast_plan[i, ]
  contrast <- c("group", info$treatment, info$control)
  message("Primary contrast ", i, "/12: ", info$contrast_id)

  res <- DESeq2::results(dds, contrast = contrast, alpha = alpha)
  shrunk <- DESeq2::lfcShrink(dds, res = res, type = "ashr", quiet = TRUE)
  raw_tab <- annotate_de(res, annotation)
  shrunk_tab <- annotate_de(shrunk, annotation)
  shrunk_tab <- shrunk_tab[match(raw_tab$ENSEMBL, shrunk_tab$ENSEMBL), , drop = FALSE]

  ## Both tables refer to the same Wald significance test.
  stopifnot(identical(raw_tab$ENSEMBL, shrunk_tab$ENSEMBL),
            isTRUE(all.equal(raw_tab$pvalue, shrunk_tab$pvalue, check.attributes = FALSE)),
            isTRUE(all.equal(raw_tab$padj, shrunk_tab$padj, check.attributes = FALSE)))

  ## Explicit provenance in both the filenames AND CSV columns. The padj
  ## remains the original Wald/DESeq2 result even in the ashr LFC table.
  raw_tab$lfc_method <- "Wald_unshrunken"
  shrunk_tab$lfc_method <- "ashr_shrunken"
  utils::write.csv(raw_tab, file.path(raw_dir, paste0(info$contrast_id, "_vs_control_unshrunken.csv")), row.names = FALSE)
  utils::write.csv(shrunk_tab, file.path(shrunk_dir, paste0(info$contrast_id, "_vs_control_ashr.csv")), row.names = FALSE)

  counts <- deg_counts(raw_tab)
  summary_rows[[i]] <- data.frame(contrast_id = info$contrast_id,
    cell_line = info$cell_line, genotype = info$genotype, radiation = info$radiation,
    n_total = counts[["n_total"]], n_pvalue = counts[["n_pvalue"]],
    n_padj = counts[["n_padj"]], n_DEG = counts[["n_DEG"]],
    n_up = counts[["n_up"]], n_down = counts[["n_down"]])

  ## Publication-ready volcano via ekbSeq::plot_volcano(), one plot per contrast.
  ## Use ORIGINAL, unshrunken Wald log2FC and BH-adjusted padj. No clipping,
  ## percentile trimming, subsampling or modification of real data values.
  ## Each panel receives independent x/y limits based on its actual extrema.
  v <- raw_tab[is.finite(raw_tab$log2FoldChange) & !is.na(raw_tab$padj), , drop = FALSE]
  if (!nrow(v)) stop("No plottable volcano data for ", info$contrast_id)
  v$SYMBOL[is.na(v$SYMBOL) | !nzchar(v$SYMBOL)] <- v$ENSEMBL[is.na(v$SYMBOL) | !nzchar(v$SYMBOL)]

  ## Exact zeros can arise from floating-point underflow for extremely small
  ## adjusted p-values. Substitute the smallest positive double ONLY in the
  ## plotting copy so -log10(padj) remains finite; CSV results are unchanged.
  if (any(v$padj < 0, na.rm = TRUE)) stop("Negative adjusted p-value encountered.")
  if (any(v$padj == 0)) {
    message(info$contrast_id, ": plotting ", sum(v$padj == 0), " underflowed padj value(s) at machine precision.")
    v$padj[v$padj == 0] <- .Machine$double.xmin
  }

  x_range <- range(v$log2FoldChange)
  x_margin <- max(0.2, diff(x_range) * 0.055)
  x_limits <- x_range + c(-x_margin, x_margin)
  y_max <- max(-log10(v$padj))
  y_limits <- c(-0.2, y_max * 1.16 + 0.75) # room for extreme points and gene labels

  ## Custom colors are determined from ORIGINAL Wald padj and LFC direction.
  ## Important: lfcThres = 0 does NOT introduce an additional DEG cutoff.
  volcano_color <- ifelse(v$direction == "up", npg_colors[["red"]],
                          ifelse(v$direction == "down", npg_colors[["navy"]], "#B8B8B8"))
  names(volcano_color) <- ifelse(v$direction == "up", "Upregulated",
                                ifelse(v$direction == "down", "Downregulated", "Not significant"))

  genotype_label <- if (info$genotype == "wt") "WT" else "p21-negative"
  radiation_label <- if (info$radiation == "7gy") "7 Gy" else info$radiation
  volcano_panels[[i]] <- ekbSeq::plot_volcano(
    v, nlabel = 2, title = paste(info$cell_line, genotype_label, radiation_label, sep = " | "),
    pointsize = 0.8, labSize = 2.5, plot.signif = TRUE, pThres = alpha, lfcThres = 0,
    xlim = x_limits, ylim = y_limits, colCustom = volcano_color,
    colAlpha = 0.60, caption = "", subtitle = "", legendPosition = "none",
    axisLabSize = 9, titleLabSize = 10, drawConnectors = TRUE,
    widthConnectors = 0.25, max.overlaps = 10, maxoverlapsConnectors = 10,
    cutoffLineType = "dashed", cutoffLineCol = "grey65", cutoffLineWidth = 0.35,
    xlab = expression(log[2]~"fold change"), ylab = expression(-log[10]~"adjusted p-value"),
    border = "full", borderWidth = 0.5, borderColour = "black"
  ) + p21_theme(9) + ggplot2::theme(
    panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.45),
    axis.line = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", size = 10, margin = ggplot2::margin(b = 3)),
    plot.margin = ggplot2::margin(7, 9, 7, 9),
    legend.position = "none"
  )

}

## ******
## One concise DEG overview for interpretation and the analysis report
summary_de <- do.call(rbind, summary_rows)
summary_de$pct_DEG_of_pvalue_valid <- round(100 * summary_de$n_DEG / summary_de$n_pvalue, 1)
stopifnot(all(summary_de$n_DEG == summary_de$n_up + summary_de$n_down))
utils::write.csv(summary_de, file.path(primary_dir, "tables", "DEG_summary.csv"), row.names = FALSE)
print(summary_de[, c("cell_line", "genotype", "radiation", "n_DEG", "n_up", "n_down")], row.names = FALSE)

## ************************************************************************************************
## GraphPad/Western-blot-style grouped, stacked DEG counts
##
## Numeric x positions give exact bar spacing: paired bars TOUCH (no dodge gap),
## while conditions and cell lines are separated by deliberate whitespace.
## The second and third label tiers are placed below the x-axis using
## coord_cartesian(clip = "off"); only the left and bottom axes remain visible.
## WT bars are opaque; p21-negative bars use alpha = 0.7 in BOTH figures.
make_deg_positions <- function(sensitivity = FALSE, bar_width = 0.62) {
  rows <- list()
  x_left <- 0
  i <- 0L
  for (cell in c("HT1080", "MCF10A")) {
    for (radiation in c("7gy", "UVA", "UVB")) {
      for (genotype in c("wt", "p21_neg")) {
        for (model in if (sensitivity) c("A", "B") else "A") {
          i <- i + 1L
          rows[[i]] <- data.frame(cell_line = cell, radiation = radiation,
                                  genotype = genotype, model = model, x = x_left + bar_width / 2)
          x_left <- x_left + bar_width
        }
        ## Separate WT from p21− in the sensitivity plot; A and B still touch.
        if (sensitivity && genotype == "wt") x_left <- x_left + 0.22
      }
      if (radiation != "UVB") x_left <- x_left + 0.80
    }
    if (cell == "HT1080") x_left <- x_left + 1.35
  }
  do.call(rbind, rows)
}

make_group_brackets <- function(positions, keys, y, label_fn, bar_width = 0.62) {
  groups <- unique(positions[keys])
  do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
    sub <- positions
    for (key in keys) sub <- sub[sub[[key]] == groups[[key]][i], , drop = FALSE]
    data.frame(x = min(sub$x) - bar_width / 2 + 0.035,
               xend = max(sub$x) + bar_width / 2 - 0.035,
               xmid = mean(range(sub$x)), y = y, label = label_fn(groups[i, , drop = FALSE]))
  }))
}

## One plotting function is used for the primary (12 bars) and sensitivity
## (24 bars) figures. It does not change DESeq2 results or exported tables.
plot_deg_bars <- function(counts_wide, sensitivity = FALSE) {
  bar_width <- 0.62
  positions <- make_deg_positions(sensitivity, bar_width)
  keys <- c("cell_line", "radiation", "genotype", "model")
  if (sensitivity) {
    input <- do.call(rbind, lapply(c("A", "B"), function(model) {
      data.frame(counts_wide[, c("cell_line", "radiation", "genotype")], model = model,
                 n_DEG = counts_wide[[paste0("n_DEG_", model)]],
                 n_up = counts_wide[[paste0("n_up_", model)]],
                 n_down = counts_wide[[paste0("n_down_", model)]])
    }))
  } else {
    input <- data.frame(counts_wide[, c("cell_line", "radiation", "genotype")], model = "A",
                        n_DEG = counts_wide$n_DEG, n_up = counts_wide$n_up,
                        n_down = counts_wide$n_down)
  }
  totals <- merge(input, positions, by = keys, sort = FALSE)
  stopifnot(nrow(totals) == nrow(positions), all(totals$n_up + totals$n_down == totals$n_DEG))
  bars <- rbind(data.frame(totals, direction = "Upregulated", n = totals$n_up),
                data.frame(totals, direction = "Downregulated", n = totals$n_down))
  bars$direction <- factor(bars$direction, levels = c("Downregulated", "Upregulated"))
  y_top <- ceiling(max(totals$n_DEG) * 1.17 / 500) * 500

  ## Multi-level annotation like a Western blot: bar labels, brackets for
  ## irradiation, and a larger bracket to identify each cell line.
  if (sensitivity) {
    bar_labels <- transform(positions, label = model, y = -0.036 * y_top)
    genotype_lines <- make_group_brackets(positions, c("cell_line", "radiation", "genotype"),
                                         -0.085 * y_top,
                                         function(r) if (r$genotype == "wt") "WT" else "p21−")
    genotype_lines$label_y <- -0.119 * y_top
    condition_y <- -0.163 * y_top
    condition_label_y <- -0.196 * y_top
    cell_line_y <- -0.245 * y_top
    cell_line_label_y <- -0.279 * y_top
  } else {
    bar_labels <- transform(positions, label = ifelse(genotype == "wt", "WT", "p21−"),
                            y = -0.036 * y_top)
    genotype_lines <- NULL
    condition_y <- -0.088 * y_top
    condition_label_y <- -0.124 * y_top
    cell_line_y <- -0.170 * y_top
    cell_line_label_y <- -0.206 * y_top
  }
  conditions <- make_group_brackets(positions, c("cell_line", "radiation"), condition_y,
                                    function(r) if (r$radiation == "7gy") "7 Gy" else r$radiation)
  conditions$label_y <- condition_label_y
  cell_lines <- make_group_brackets(positions, "cell_line", cell_line_y, function(r) r$cell_line)
  cell_lines$label_y <- cell_line_label_y

  base <- ggplot2::ggplot(bars, ggplot2::aes(x = x, y = n, fill = direction, alpha = genotype)) +
    ggplot2::geom_col(width = bar_width, color = "black", linewidth = 0.32) +
    ggplot2::geom_text(data = totals,
                       ggplot2::aes(x = x, y = n_DEG, label = scales::comma(n_DEG)),
                       inherit.aes = FALSE, vjust = -0.55,
                       size = if (sensitivity) 2.65 else 3.15, fontface = "bold") +
    ggplot2::geom_text(data = bar_labels, ggplot2::aes(x = x, y = y, label = label),
                       inherit.aes = FALSE, size = if (sensitivity) 3.0 else 3.35,
                       fontface = "bold") +
    ggplot2::geom_segment(data = conditions,
                          ggplot2::aes(x = x, xend = xend, y = y, yend = y),
                          inherit.aes = FALSE, linewidth = 0.42, color = "black") +
    ggplot2::geom_text(data = conditions,
                       ggplot2::aes(x = xmid, y = label_y, label = label),
                       inherit.aes = FALSE, size = 3.45) +
    ggplot2::geom_segment(data = cell_lines,
                          ggplot2::aes(x = x, xend = xend, y = y, yend = y),
                          inherit.aes = FALSE, linewidth = 0.50, color = "black") +
    ggplot2::geom_text(data = cell_lines,
                       ggplot2::aes(x = xmid, y = label_y, label = label),
                       inherit.aes = FALSE, size = 4.0, fontface = "bold") +
    ggplot2::scale_fill_manual(values = c(Downregulated = npg_colors[["navy"]],
                                            Upregulated = npg_colors[["red"]]), name = NULL) +
    ggplot2::scale_alpha_manual(values = c(wt = 1, p21_neg = 0.7), guide = "none") +
    ggplot2::scale_x_continuous(limits = range(positions$x) + c(-0.75, 0.75),
                                breaks = NULL, expand = ggplot2::expansion(mult = 0)) +
    ggplot2::scale_y_continuous(labels = scales::label_comma(),
                                breaks = scales::breaks_pretty(n = 6),
                                expand = ggplot2::expansion(mult = 0)) +
    ggplot2::coord_cartesian(ylim = c(0, y_top), clip = "off", expand = FALSE) +
    ggplot2::labs(
      title = if (sensitivity) "Sensitivity analysis: sequencing-run adjustment" else
        "Transcriptional responses to irradiation",
      subtitle = if (sensitivity)
        "Model A: processing + group  |  Model B: processing + sequencing run + group" else
        "DESeq2 Model A  |  BH-adjusted p < 0.05  |  no log2FC cutoff",
      x = NULL, y = "Differentially expressed genes") +
    p21_theme(11) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.65),
      axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.65),
      axis.ticks.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
      legend.position = "top", legend.justification = "left",
      legend.text = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(face = "bold", size = 17),
      plot.subtitle = ggplot2::element_text(size = 10.5, color = "grey30"),
      plot.margin = ggplot2::margin(t = 8, r = 18, b = if (sensitivity) 154 else 115, l = 9))

  if (sensitivity) {
    base <- base +
      ggplot2::geom_segment(data = genotype_lines,
                            ggplot2::aes(x = x, xend = xend, y = y, yend = y),
                            inherit.aes = FALSE, linewidth = 0.35, color = "black") +
      ggplot2::geom_text(data = genotype_lines,
                         ggplot2::aes(x = xmid, y = label_y, label = label),
                         inherit.aes = FALSE, size = 3.25, fontface = "bold") +
      ggplot2::labs(caption = "Model A is primary; Model B is a sensitivity analysis only.") +
      ggplot2::theme(plot.caption = ggplot2::element_text(size = 9, color = "grey35"))
  }
  base
}

## ************************************************************************************************
## Primary figure 1: 12 bars in ONE graph (no facets), WT and p21− touching.
## Stacked up/down segments and black outlines; alpha = 0.7 for p21−.
p_deg <- plot_deg_bars(summary_de)
ekbSeq::export_plot_dual(file.path(fig_dir, "DEG_counts_by_condition"), p_deg,
                         width = 13.8, height = 7.7, rasterize = FALSE)

## ************************************************************************************************
## Primary figure 2: volcano atlas from twelve ekbSeq::plot_volcano() objects
## Each panel is scaled to its OWN real minimum/maximum (no percentile caps).
## All statistical values remain unchanged in the exported CSVs.
volcano_atlas <- patchwork::wrap_plots(volcano_panels, ncol = 3) +
  patchwork::plot_annotation(
    title = "Irradiation response: primary differential expression",
    subtitle = paste("DESeq2 Model A | unshrunken log2FC | BH-adjusted p < 0.05 | no log2FC cutoff",
                     "| red: up; blue: down; grey: not significant"),
    caption = "Independent axis ranges per panel; no log2FC or adjusted p-value truncation.",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 17),
                           plot.subtitle = ggplot2::element_text(size = 10, color = "grey30"),
                           plot.caption = ggplot2::element_text(size = 9, color = "grey35")))
ekbSeq::export_plot_dual(file.path(fig_dir, "irradiation_12_panel_volcano_atlas"), volcano_atlas,
                         width = 19, height = 18, rasterize = FALSE)

## ******
## Reporting atlases: two cell-line-specific six-panel figures.
## Use precisely the same Wald/DESeq2 panels as the existing 12-panel atlas.
## Rows: WT / p21−. Columns: 7 Gy / UVA / UVB. No tests are repeated.
## A separate, visible color key avoids adding 6 redundant panel legends.
make_volcano_key <- function(labels) {
  stopifnot(length(labels) == 3L)
  key <- data.frame(x = c(1, 5, 9), y = 1, label = labels,
    color = c(npg_colors[["red"]], npg_colors[["navy"]], "#B8B8B8"))
  ggplot2::ggplot(key, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(color = color), size = 3.1) +
    ggplot2::geom_text(ggplot2::aes(x = x + 0.17, label = label),
                       hjust = 0, size = 3.5) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_x_continuous(limits = c(0.6, 13), expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = c(0.65, 1.35), clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 5, 0, 5))
}

primary_key <- make_volcano_key(c("Upregulated", "Downregulated", "Not significant"))
for (cell in c("HT1080", "MCF10A")) {
  idx <- which(contrast_plan$cell_line == cell)
  stopifnot(length(idx) == 6L,
            identical(contrast_plan$radiation[idx], rep(c("7gy", "UVA", "UVB"), 2)),
            identical(contrast_plan$genotype[idx], rep(c("wt", "p21_neg"), each = 3)))
  six <- (patchwork::wrap_plots(volcano_panels[idx], ncol = 3) / primary_key) +
    patchwork::plot_layout(heights = c(22, 1)) +
    patchwork::plot_annotation(
      title = paste(cell, "| irradiation versus control"),
      subtitle = "Top: WT | Bottom: p21 KO | Wald log2FC | BH FDR < 0.05 | no LFC cutoff",
      caption = "Panel-specific observed axis limits; colors encode significant response direction.",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 16),
        plot.subtitle = ggplot2::element_text(size = 10, color = "grey30"),
        plot.caption = ggplot2::element_text(size = 9, color = "grey35")))
  ekbSeq::export_plot_dual(file.path(fig_dir, paste0("irradiation_volcano_", cell, "_6panel")),
                            six, width = 15.6, height = 10.7, rasterize = FALSE)
}

## ************************************************************************************************
## Model B: a SINGLE sensitivity analysis, not a second biological DEG set
## The run coefficient is weakly anchored (2 vs 3 HT1080 WT, 7-Gy UKD
## replicates); therefore we do not treat B-only DEGs as a new signature.
## We export ONE compact table comparing A/B across the twelve contrasts.
message("Fitting Model B for sensitivity analysis only ...")
dds_B <- DESeq2::DESeqDataSetFromMatrix(
  countData = DESeq2::counts(dds, normalized = FALSE), colData = meta,
  design = ~ processing + sequencing_run + group)
dds_B <- DESeq2::DESeq(dds_B, quiet = TRUE, parallel = FALSE)

sensitivity_rows <- vector("list", nrow(contrast_plan))
for (i in seq_len(nrow(contrast_plan))) {
  info <- contrast_plan[i, ]
  contrast <- c("group", info$treatment, info$control)
  a <- DESeq2::results(dds, contrast = contrast, alpha = alpha)
  b <- DESeq2::results(dds_B, contrast = contrast, alpha = alpha)
  stopifnot(identical(rownames(a), rownames(b)))

  sig_a <- !is.na(a$padj) & a$padj < alpha
  sig_b <- !is.na(b$padj) & b$padj < alpha
  comparable <- is.finite(a$log2FoldChange) & is.finite(b$log2FoldChange)
  overlap <- sum(sig_a & sig_b)
  union <- sum(sig_a | sig_b)
  sensitivity_rows[[i]] <- data.frame(contrast_id = info$contrast_id,
    cell_line = info$cell_line, genotype = info$genotype, radiation = info$radiation,
    n_DEG_A = sum(sig_a), n_DEG_B = sum(sig_b),
    n_up_A = sum(sig_a & a$log2FoldChange > 0, na.rm = TRUE),
    n_down_A = sum(sig_a & a$log2FoldChange < 0, na.rm = TRUE),
    n_up_B = sum(sig_b & b$log2FoldChange > 0, na.rm = TRUE),
    n_down_B = sum(sig_b & b$log2FoldChange < 0, na.rm = TRUE),
    n_shared = overlap,
    n_A_only = sum(sig_a & !sig_b), n_B_only = sum(!sig_a & sig_b),
    jaccard = if (union) overlap / union else NA_real_,
    pearson_log2FC = if (sum(comparable) > 2L) stats::cor(a$log2FoldChange[comparable], b$log2FoldChange[comparable]) else NA_real_,
    median_abs_delta_log2FC = if (any(comparable)) stats::median(abs(a$log2FoldChange[comparable] -
                                                   b$log2FoldChange[comparable])) else NA_real_,
    n_padj_A = sum(!is.na(a$padj)), n_padj_B = sum(!is.na(b$padj)))
}
sensitivity <- do.call(rbind, sensitivity_rows)
utils::write.csv(sensitivity, file.path(sens_tab_dir, "model_A_vs_B_summary.csv"), row.names = FALSE)

## Sensitivity figure: ONE graph with touching A/B bars inside each
## genotype and irradiation; independent of downstream biological interpretation.
stopifnot(all(sensitivity$n_DEG_A == sensitivity$n_up_A + sensitivity$n_down_A),
          all(sensitivity$n_DEG_B == sensitivity$n_up_B + sensitivity$n_down_B))
p_sens <- plot_deg_bars(sensitivity, sensitivity = TRUE)
ekbSeq::export_plot_dual(file.path(sens_fig_dir, "DEG_counts_A_vs_B"), p_sens,
                         width = 18.5, height = 8.3, rasterize = FALSE)

## ************************************************************************************************
## Summary and boundaries for downstream analyses
## Read primary shrunken CSVs for response PCA/heatmap, unshrunken CSVs for
## Wald-based reporting. The sensitivity directory must not become a second
## biological input, and no downstream analyses need Model B result objects.
message("Primary results: ", primary_dir)
message("Sensitivity check: ", sensitivity_dir)
message("Done: 12 primary contrasts, 12 ashr-shrunken tables, one DEG summary.")