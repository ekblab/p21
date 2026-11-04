## ************************************************************************************************
## p21 bulk RNA-seq
## 05 - Integration of irradiation-response modules and formal p21 interactions (24 h)
##
## Biological questions:
##   Which previously defined response modules carry p21-modulated genes?
##   Do formal p21 interactions agree with the observed WT/KO response patterns?
##   Are some modules candidates for a later, biologically justified merger?
##
## Analysis decisions:
##   - KEEP the eight L2 modules from script 03 unchanged. No re-clustering.
##   - Read only the primary ashr responses (script 03) and directly ashr-
##     shrunken four-group interactions (script 04; Model A).
##   - BH-adjusted interaction padj < 0.05; no log2FC threshold.
##   - Module-level interaction effects are DESCRIPTIVE averages across genes.
##     Directly shrunken interaction LFCs are NOT KO_shrunken - WT_shrunken.
##   - Exploratory one-sided Fisher tests ask if significant interaction genes
##     are enriched against OTHER TESTABLE L2-module genes (not the genome).
##     Adjust 8 x 6 Fisher p-values together using BH; this is a secondary,
##     conditional analysis because modules were inferred from these data.
##   - GO-BP remains in script 03; no new GO/GSEA, IPA or DREAM here.
##   - No automatic module merging or final renaming.
##
## Output (all below Results/05_module_integration/):
##   tables/module_integration_summary.csv: 48 module x contrast summaries
##   tables/significant_module_interaction_genes.csv: significant gene-level hits
##   tables/module_response_profile_similarity.csv: 28 module-pair correlations
##   tables/module_WT_aligned_response_change.csv: interpretable WT/KO changes
##   figures/module_response_pairs.{png,svg}: general WT/KO response programs
##   figures/module_response_strength_change.{png,svg}: WT-oriented strengthening/weakening
##   tables/module_pair_subsampling_stability.csv: descriptive robustness of four pairs
##   One processed RDS remains untouched; this script creates no RDS files.
## ************************************************************************************************

rm(list = ls())
gc()

## ******
## Project setup, dependencies and paths
if (!requireNamespace("here", quietly = TRUE)) stop("Install the 'here' package.")
here::i_am("R/05_module_integration.R")
source(here::here("R", "00_setup.R"))
source(here::here("R", "plot_style.R"))
for (pkg in c("ggplot2", "ekbSeq", "scales")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing R package: ", pkg)
}

alpha <- 0.05
root <- here::here("Results", "05_module_integration")
tab_dir <- file.path(root, "tables")
fig_dir <- file.path(root, "figures")
mod_dir <- here::here("Results", "03_response_architecture", "02_MODULE_L2")
response_path <- here::here("Results", "03_response_architecture", "response_matrix_DEG_union_ashr.csv")
interaction_dir <- here::here("Results", "04_p21_interactions", "tables", "shrunken")
for (d in c(tab_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

status <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))
  flush.console()
}
plot_theme <- function(size = 10) {
  if (exists("p21_theme", mode = "function")) p21_theme(size) else ggplot2::theme_classic(base_size = size)
}

## ************************************************************************************************
## Read fixed L2 module assignments and original ashr response profiles
## Row alignment is ALWAYS by ENSEMBL; never by CSV row position or gene symbol.
status("1/5 | Reading the fixed L2 modules and twelve irradiation-response profiles")
mem_path <- file.path(mod_dir, "module_membership.csv")
required <- c(mem_path, response_path)
if (any(!file.exists(required))) stop("Missing script-03 outputs:\n", paste(required[!file.exists(required)], collapse = "\n"))
mem <- utils::read.csv(mem_path, stringsAsFactors = FALSE)
resp <- utils::read.csv(response_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("ENSEMBL", "SYMBOL", "module") %in% names(mem)) ||
    !all(c("ENSEMBL", "SYMBOL") %in% names(resp))) stop("Missing mandatory module/response columns.")
if (anyNA(mem$ENSEMBL) || anyDuplicated(mem$ENSEMBL) ||
    anyNA(resp$ENSEMBL) || anyDuplicated(resp$ENSEMBL) ||
    !setequal(mem$ENSEMBL, resp$ENSEMBL)) stop("Module membership and response matrix have incompatible gene IDs.")
module_ids <- sprintf("M%02d", 1:8)
if (!setequal(unique(mem$module), module_ids)) stop("Expected the original eight L2 modules (M01-M08).")
mem <- mem[order(match(mem$module, module_ids)), , drop = FALSE]
resp <- resp[match(mem$ENSEMBL, resp$ENSEMBL), , drop = FALSE]
stopifnot(identical(mem$ENSEMBL, resp$ENSEMBL))

## Preserve the exact twelve-column order from script 03, with paired WT/KO.
profiles <- expand.grid(genotype = c("wt", "p21_neg"), radiation = c("7gy", "UVA", "UVB"),
                        cell_line = c("HT1080", "MCF10A"), stringsAsFactors = FALSE)
profiles$id <- with(profiles, paste(cell_line, genotype, radiation, sep = "_"))
if (!setequal(names(resp)[-(1:2)], profiles$id)) stop("Unexpected response columns; rerun script 03.")
response <- as.matrix(resp[, profiles$id, drop = FALSE])
storage.mode(response) <- "double"
if (any(!is.finite(response))) stop("The module response matrix contains non-finite log2FC values.")
rownames(response) <- mem$ENSEMBL

## The previous descriptive names remain metadata, not definitive annotation.
labels_path <- file.path(mod_dir, "module_descriptive_labels.csv")
module_labels <- data.frame(module = module_ids, descriptive_pattern = NA_character_)
if (file.exists(labels_path)) {
  old <- utils::read.csv(labels_path, stringsAsFactors = FALSE)
  if (all(c("module", "descriptive_pattern") %in% names(old)))
    module_labels$descriptive_pattern <- old$descriptive_pattern[match(module_ids, old$module)]
}

## ************************************************************************************************
## Read all six directly shrunken interaction contrasts from script 04
## The 'padj' in these ashr tables is still the ORIGINAL DESeq2 Wald FDR.
status("2/5 | Reading six directly shrunken interaction tables")
plan <- expand.grid(radiation = c("7gy", "UVA", "UVB"),
                    cell_line = c("HT1080", "MCF10A"), stringsAsFactors = FALSE)
plan$id <- with(plan, paste(cell_line, radiation, sep = "_"))
paths <- file.path(interaction_dir, paste0(plan$id, "_p21_interaction_ashr.csv"))
if (any(!file.exists(paths))) stop("Missing script-04 outputs:\n", paste(paths[!file.exists(paths)], collapse = "\n"))

read_interaction <- function(path, genes) {
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  needed <- c("ENSEMBL", "log2FoldChange", "padj", "baseMean", "interaction_direction",
              "WT_response_log2FC", "p21neg_response_log2FC")
  if (!all(needed %in% names(x)) || anyNA(x$ENSEMBL) || anyDuplicated(x$ENSEMBL))
    stop("Invalid interaction result table: ", path)
  if (!"lfc_method" %in% names(x) ||
      anyNA(x$lfc_method) || any(x$lfc_method != "ashr_shrunken"))
    stop("Expected DIRECTLY ashr-shrunken interaction estimates: ", path)
  idx <- match(genes, x$ENSEMBL)
  if (anyNA(idx)) stop("Script 04 is missing L2-module genes: ", path)
  x[idx, needed, drop = FALSE]
}
tabs <- lapply(paths, read_interaction, genes = mem$ENSEMBL)
stopifnot(all(vapply(tabs, nrow, integer(1)) == nrow(mem)))

## These response columns are independently shrunken (script 02), whereas
## log2FoldChange in these tables is the DIRECTLY shrunken four-group contrast.
## In general, they must NOT be subtracted to reconstruct the latter.
## CSV vectors are unnamed; matrix columns carry ENSEMBL names. Gene alignment
## was validated above, so compare values without their names (not less strictly).
check_response <- function(csv_lfc, matrix_lfc, id, genotype) {
  ref <- unname(matrix_lfc)
  if (isTRUE(all.equal(csv_lfc, ref, tolerance = 1e-6))) return(invisible(TRUE))
  delta <- abs(csv_lfc - ref)
  delta[!is.finite(delta)] <- Inf
  worst <- which.max(delta)
  stop(sprintf("Response LFC mismatch (%s, %s), ENSEMBL %s: script 04 = %.8g; script 03 = %.8g; abs difference = %.6g",
               id, genotype, mem$ENSEMBL[worst], csv_lfc[worst], ref[worst], delta[worst]))
}
for (j in seq_len(nrow(plan))) {
  p <- plan[j, ]
  wt <- paste(p$cell_line, "wt", p$radiation, sep = "_")
  ko <- paste(p$cell_line, "p21_neg", p$radiation, sep = "_")
  check_response(tabs[[j]]$WT_response_log2FC, response[, wt], p$id, "WT")
  check_response(tabs[[j]]$p21neg_response_log2FC, response[, ko], p$id, "p21-KO")
  if (any(tabs[[j]]$baseMean < 0, na.rm = TRUE)) stop("Negative baseMean: ", p$id)
  valid_dir <- c("more_positive_in_p21neg", "more_negative_in_p21neg", "not_significant")
  if (!all(tabs[[j]]$interaction_direction %in% valid_dir)) stop("Unrecognized interaction direction: ", p$id)
  flag <- is.finite(tabs[[j]]$padj) & tabs[[j]]$padj < alpha
  if (!identical(flag, tabs[[j]]$interaction_direction != "not_significant"))
    stop("Script-04 significance/direction mismatch for ", p$id)
}

## ************************************************************************************************
## Module-level summaries and exploratory enrichment of significant interactions
## Test universe = L2 genes having a testable Wald adjusted p-value for the
## respective interaction. This controls for differing module testability.
## Fisher tests are ONE-SIDED enrichment tests, not tests of a mean LFC.
status("3/5 | Summarizing interaction direction, strength and module enrichment")
rows <- vector("list", length(module_ids) * nrow(plan))
values <- vector("list", nrow(plan))
hits <- list()
index <- 0L
for (j in seq_len(nrow(plan))) {
  x <- tabs[[j]]
  tested <- is.finite(x$padj)
  sig <- tested & x$padj < alpha
  lfc <- x$log2FoldChange
  if (any(!is.finite(lfc))) stop("Non-finite ashr interaction LFC in ", plan$id[j])
  values[[j]] <- lfc

  for (m in module_ids) {
    index <- index + 1L
    in_mod <- mem$module == m
    a <- sum(in_mod & sig)
    b <- sum(in_mod & tested & !sig)
    c <- sum(!in_mod & sig)
    d <- sum(!in_mod & tested & !sig)
    if (a + b == 0L || c + d == 0L) stop("No testable genes in module/background: ", m)
    test <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE),
                               alternative = "greater")
    p <- plan[j, ]
    wt_name <- paste(p$cell_line, "wt", p$radiation, sep = "_")
    ko_name <- paste(p$cell_line, "p21_neg", p$radiation, sep = "_")
    mean_wt <- mean(response[in_mod, wt_name])
    mean_ko <- mean(response[in_mod, ko_name])
    rows[[index]] <- data.frame(
      module = m, n_module_genes = sum(in_mod),
      descriptive_pattern = module_labels$descriptive_pattern[match(m, module_ids)],
      contrast_id = p$id, cell_line = p$cell_line, radiation = p$radiation,
      mean_WT_response_ashr = mean_wt,
      mean_p21neg_response_ashr = mean_ko,
      descriptive_KO_minus_WT = mean_ko - mean_wt,
      mean_direct_interaction_ashr = mean(lfc[in_mod]),
      median_direct_interaction_ashr = stats::median(lfc[in_mod]),
      n_tested = a + b, n_significant = a,
      pct_significant_tested = 100 * a / (a + b),
      n_more_positive_p21neg = sum(in_mod & sig & x$interaction_direction == "more_positive_in_p21neg"),
      n_more_negative_p21neg = sum(in_mod & sig & x$interaction_direction == "more_negative_in_p21neg"),
      pct_significant_other_modules = 100 * c / (c + d),
      enrichment_odds_ratio = unname(test$estimate), enrichment_p = test$p.value)
  }

  ## Export only significant module-member genes, not six duplicated full CSVs.
  ix <- which(sig)
  if (length(ix)) hits[[j]] <- data.frame(
    ENSEMBL = mem$ENSEMBL[ix], SYMBOL = mem$SYMBOL[ix], module = mem$module[ix],
    contrast_id = plan$id[j], direct_interaction_ashr = lfc[ix],
    interaction_padj = x$padj[ix], interaction_direction = x$interaction_direction[ix],
    WT_response_ashr = x$WT_response_log2FC[ix],
    p21neg_response_ashr = x$p21neg_response_log2FC[ix])
}

summary_tab <- do.call(rbind, rows)
summary_tab$enrichment_FDR_BH48 <- stats::p.adjust(summary_tab$enrichment_p, method = "BH")
summary_tab$significantly_enriched <- summary_tab$enrichment_FDR_BH48 < alpha
stopifnot(nrow(summary_tab) == 48L,
          all(summary_tab$n_significant == summary_tab$n_more_positive_p21neg +
                summary_tab$n_more_negative_p21neg))
utils::write.csv(summary_tab, file.path(tab_dir, "module_integration_summary.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, hits), file.path(tab_dir, "significant_module_interaction_genes.csv"), row.names = FALSE)
message("Significant module-by-condition enrichments (exploratory BH < 0.05): ",
        sum(summary_tab$significantly_enriched), "/48")

## ******
## Diagnostic only: pairwise similarity of the eight fixed mean response patterns.
## A high correlation does NOT automatically justify merging two gene modules.
centers <- t(vapply(module_ids, function(m) colMeans(response[mem$module == m, , drop = FALSE]),
                    numeric(ncol(response))))
rownames(centers) <- module_ids
pair_idx <- utils::combn(module_ids, 2)
pair_similarity <- data.frame(module_A = pair_idx[1, ], module_B = pair_idx[2, ],
  pearson_mean_response = vapply(seq_len(ncol(pair_idx)), function(i)
    stats::cor(centers[pair_idx[1, i], ], centers[pair_idx[2, i], ]), numeric(1)))
utils::write.csv(pair_similarity, file.path(tab_dir, "module_response_profile_similarity.csv"), row.names = FALSE)

## ************************************************************************************************
## Prepare descriptive means for the two complementary module heatmaps
## Original WT/KO treatment responses and direct p21 interactions have DIFFERENT
## estimators. The second heatmap is a clearly labeled DESCRIPTIVE comparison
## of WT/KO response means, not a directly shrunken interaction estimator.
response_dat <- do.call(rbind, lapply(seq_len(nrow(profiles)), function(i) {
  p <- profiles[i, ]
  data.frame(module = module_ids, cell_line = p$cell_line, radiation = p$radiation,
             genotype = p$genotype, profile_id = p$id, mean_LFC = centers[, p$id])
}))
response_dat$module <- factor(response_dat$module, levels = rev(module_ids))
response_dat$cell_line <- factor(response_dat$cell_line, levels = c("HT1080", "MCF10A"))
response_dat$profile_id <- factor(response_dat$profile_id, levels = profiles$id)
summary_tab$module <- factor(summary_tab$module, levels = rev(module_ids))
summary_tab$cell_line <- factor(summary_tab$cell_line, levels = c("HT1080", "MCF10A"))
summary_tab$radiation <- factor(summary_tab$radiation, levels = c("7gy", "UVA", "UVB"))

blue_red <- c("#3C5488", "#F8F8F8", "#E64B35")
response_tick <- with(profiles, paste(ifelse(radiation == "7gy", "7 Gy", radiation),
                                      ifelse(genotype == "wt", "WT", "KO"), sep = "\n"))
names(response_tick) <- profiles$id

## ************************************************************************************************
## Publication-style heatmaps: paired response modules
## These are two DIFFERENT descriptive quantities; separate color limits are
## intentional. They follow the layout of the original script-03 heatmap.
## Strong negative correlations between module mean profiles motivate the
## pairwise row order; they do NOT imply a shared pathway or justify merging.
status("4/5 | Exporting two concise heatmaps with paired module order")
pair_order <- c("M01", "M08", "M02", "M06", "M03", "M05", "M04", "M07")
pair_names <- module_labels$descriptive_pattern[match(pair_order, module_ids)]
pair_names[is.na(pair_names) | !nzchar(pair_names)] <- "Response pattern"
row_labels <- setNames(paste(pair_order, "·", pair_names), pair_order)
row_levels <- rev(unname(row_labels))

response_pair <- response_dat
response_pair$module <- factor(row_labels[as.character(response_pair$module)], levels = row_levels)
interaction_pair <- summary_tab
interaction_pair$module <- factor(row_labels[as.character(interaction_pair$module)], levels = row_levels)

## Tiles show original WT/KO responses or a WT-oriented descriptive
## comparison. The latter is NOT the directly shrunken interaction.
## Tile labels improve readability of modest module means.
tile_label <- function(x) ifelse(abs(x) < 0.005, "0.00", sprintf("%+.2f", x))
response_pair$value_label <- tile_label(response_pair$mean_LFC)

pair_theme <- plot_theme(11) + ggplot2::theme(
  panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
  panel.spacing.x = grid::unit(0.22, "cm"),
  strip.background = ggplot2::element_blank(),
  strip.text = ggplot2::element_text(face = "bold", size = 12),
  axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
  axis.text.x = ggplot2::element_text(size = 10, color = "black"),
  axis.text.y = ggplot2::element_text(size = 10, color = "black"),
  legend.position = "right", legend.title = ggplot2::element_text(size = 10),
  plot.title = ggplot2::element_text(face = "bold", size = 16),
  plot.subtitle = ggplot2::element_text(size = 10, color = "grey30"),
  plot.caption = ggplot2::element_text(size = 9, color = "grey35"),
  plot.margin = ggplot2::margin(10, 12, 10, 12))

## A: the general responses (same data as script 03, simply reordered).
p_response_pairs <- ggplot2::ggplot(response_pair,
  ggplot2::aes(x = profile_id, y = module, fill = mean_LFC)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.95) +
  ggplot2::geom_text(ggplot2::aes(label = value_label), size = 3.1) +
  ggplot2::geom_hline(yintercept = c(2.5, 4.5, 6.5), color = "grey60", linewidth = 0.40) +
  ggplot2::facet_grid(~ cell_line, scales = "free_x", space = "free_x") +
  ggplot2::scale_x_discrete(labels = response_tick, drop = TRUE) +
  ggplot2::scale_fill_gradient2(low = blue_red[1], mid = "white", high = blue_red[3],
    midpoint = 0, limits = c(-2, 2), oob = scales::squish, name = "Mean ashr\nlog2FC") +
  ggplot2::labs(title = "Radiation-response modules: opposing expression programs",
    subtitle = "Original WT and p21-KO treatment effects | modules grouped by opposing response pattern",
    x = NULL, y = NULL,
    caption = "Means are descriptive; anticorrelated module profiles do not establish distinct mechanisms.") +
  pair_theme

## B: HOW the original radiation response changes after p21 loss.
## Orient the descriptive KO-WT difference to the SIGN of the WT module mean:
##   WT induction:  KO - WT;  WT repression: WT - KO.
## Positive = stronger response in the WT direction; negative = weaker or reversed.
## This is NOT a directly shrunken interaction and NOT a module-wide test.
## If the KO mean changes sign, flag it explicitly; an inversion is not simply
## a weaker response. Very weak WT means (<0.10) are flagged in the TABLE
## because their sign is an uncertain reference (not another hypothesis test).
interaction_pair$wt_aligned_change <- with(interaction_pair,
  sign(mean_WT_response_ashr) * (mean_p21neg_response_ashr - mean_WT_response_ashr))
interaction_pair$direction_reversal <- with(interaction_pair,
  mean_WT_response_ashr * mean_p21neg_response_ashr < 0)
interaction_pair$weak_WT_reference <- abs(interaction_pair$mean_WT_response_ashr) < 0.10
interaction_pair$value_label <- paste0(tile_label(interaction_pair$wt_aligned_change),
  ifelse(interaction_pair$direction_reversal, "†", ""))

## Save the interpretive quantities alongside the existing formal statistics.
## No filtering, re-fitting, or change to Fisher tests or Wald FDR values.
response_strength <- interaction_pair[, c("module", "contrast_id", "cell_line", "radiation",
  "mean_WT_response_ashr", "mean_p21neg_response_ashr", "wt_aligned_change",
  "direction_reversal", "weak_WT_reference")]
response_strength$module <- as.character(summary_tab$module) # retain M01-M08 IDs, not plotting labels
utils::write.csv(response_strength, file.path(tab_dir, "module_WT_aligned_response_change.csv"),
                 row.names = FALSE)

max_effect <- max(abs(interaction_pair$wt_aligned_change), na.rm = TRUE)
color_limit <- max(0.05, ceiling(max_effect / 0.05) * 0.05)
p_interaction_pairs <- ggplot2::ggplot(interaction_pair,
  ggplot2::aes(x = radiation, y = module, fill = wt_aligned_change)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.95) +
  ggplot2::geom_text(ggplot2::aes(label = value_label), size = 3.6) +
  ggplot2::geom_hline(yintercept = c(2.5, 4.5, 6.5), color = "grey60", linewidth = 0.40) +
  ggplot2::facet_grid(~ cell_line, scales = "free_x", space = "free_x") +
  ggplot2::scale_x_discrete(labels = c("7gy" = "7 Gy", "UVA" = "UVA", "UVB" = "UVB")) +
  ggplot2::scale_fill_gradient2(low = blue_red[1], mid = "white", high = blue_red[3],
    midpoint = 0, limits = c(-color_limit, color_limit), name = "Change in WT-aligned\nmean response (log2FC)") +
  ggplot2::labs(title = "p21 loss changes the strength of radiation-response programs",
    subtitle = "Blue: attenuated WT response  |  Red: enhanced WT response  |  descriptive module means",
    x = NULL, y = NULL,
    caption = paste("WT-aligned change = (KO mean − WT mean) × sign(WT mean).",
                    "† indicates a direction reversal; weak WT reference (|mean| < 0.10) flagged in CSV.",
                    "Formal gene-level Wald FDR and exploratory module enrichment are reported separately.")) +
  pair_theme

ekbSeq::export_plot_dual(file.path(fig_dir, "module_response_pairs"), p_response_pairs,
                         width = 19.5, height = 7.5, rasterize = FALSE)
ekbSeq::export_plot_dual(file.path(fig_dir, "module_response_strength_change"), p_interaction_pairs,
                         width = 15.2, height = 7.5, rasterize = FALSE)

## ******
## Conditional half-gene resampling: is each inverse pair driven by many genes?
## This assesses stability of DESCRIPTIVE centroid correlations under subsampling
## *within fixed modules*. It is not a p-value, not an independent validation,
## and does not address biological-replicate variance or clustering selection.
status("5/5 | Checking descriptive stability of four inverse response pairs")
set.seed(20261012)
pair_list <- list(c("M01", "M08"), c("M02", "M06"),
                  c("M03", "M05"), c("M04", "M07"))
n_split <- 200L
pair_rows <- lapply(pair_list, function(ab) {
  a <- which(mem$module == ab[1]); b <- which(mem$module == ab[2])
  r_full <- stats::cor(centers[ab[1], ], centers[ab[2], ])
  r_half <- replicate(n_split, {
    aa <- sample(a, max(2L, floor(length(a) / 2)))
    bb <- sample(b, max(2L, floor(length(b) / 2)))
    stats::cor(colMeans(response[aa, , drop = FALSE]),
               colMeans(response[bb, , drop = FALSE]))
  })
  data.frame(module_A = ab[1], module_B = ab[2],
    pearson_full = r_full, median_half_gene_r = stats::median(r_half),
    q05_half_gene_r = unname(stats::quantile(r_half, 0.05)),
    q95_half_gene_r = unname(stats::quantile(r_half, 0.95)),
    fraction_half_gene_r_below_minus_0_9 = mean(r_half < -0.9),
    n_resamples = n_split)
})
utils::write.csv(do.call(rbind, pair_rows),
                 file.path(tab_dir, "module_pair_subsampling_stability.csv"), row.names = FALSE)

## ************************************************************************************************
## Summary and boundaries
## Correlated genes, module selection from the same data and differences in
## testability limit the inferential interpretation of Fisher enrichment.
## An enriched module is NOT necessarily uniformly p21-dependent, and a module
## without enrichment is NOT proven p21-independent. Biological module names
## and any merging remain deferred until these patterns are reviewed.
status("Complete |: 8 unchanged modules, 6 contrasts, 48 summaries")
message("Results: ", root)
