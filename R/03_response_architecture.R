## p21 / 03 | Primary response analysis and GO enrichment
## Run directly to execute the COMPLETE workflow (this file + both next stages):
##   source("R/03_response_architecture.R")
## For selective re-evaluation or replotting, source the dedicated stage.
##
## Stage 1: twelve irradiation responses, k=8 L2 modules, GO-BP ORA,
## Wang semantic redundancy reduction, and processing-matched diagnostics.
## All biological statistics and semantically reduced GO terms are exported
## before the two downstream stages run. No flags or run-mode switches.

if (!requireNamespace("here", quietly = TRUE)) stop("Missing package: here")
source(here::here("R", "03_shared.R"), local = TRUE)
p21_check(c("ggplot2", "clusterProfiler", "org.Hs.eg.db", "GOSemSim",
            "scales", "ekbSeq", "DESeq2", "SummarizedExperiment"))

## ******
## Project setup, dependencies and paths
source(here::here("R", "00_setup.R"))
source(here::here("R", "plot_style.R"))
set.seed(21)
alpha <- 0.05
k <- 8L              # chosen resolution after comparing k = 4..12
n_pattern <- 5L      # genes ranked by module-profile similarity
n_biological <- 5L   # additional curated / GO-supported genes
root <- p21_result_root
primary <- here::here("Results", "02_primary_DE", "01_PRIMARY_model_A", "tables")
pca_dir <- p21_pca_dir
mod_dir <- p21_module_dir
sens_dir <- p21_sensitivity_dir
for (d in c(pca_dir, mod_dir, sens_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
## ******
## Small helpers used only by this script
## Reuse the common plot style and export PNG/SVG through ekbSeq.
plot_theme <- function(size = 11) {
  if (exists("p21_theme", mode = "function")) p21_theme(size) else ggplot2::theme_classic(base_size = size)
}
export <- function(path, p, w, h) p21_export(path, p, w, h)

## Timestamped messages remain visible when sourced in RStudio.
status <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))
  flush.console()
}

## ************************************************************************************************
## Assemble twelve primary irradiation-response profiles
## All responses originate from Model A and use ashr-shrunken log2FC.
## The DEG union determines the module input; no extra effect-size cutoff.
status("1/5 | Reading twelve ashr tables and the processed object")
profiles <- expand.grid(genotype = c("wt", "p21_neg"), radiation = c("7gy", "UVA", "UVB"),
                        cell_line = c("HT1080", "MCF10A"), stringsAsFactors = FALSE)
profiles$id <- with(profiles, paste(cell_line, genotype, radiation, sep = "_"))
profiles$pair <- with(profiles, paste(cell_line, radiation, sep = "_"))
profiles$label <- paste(ifelse(profiles$radiation == "7gy", "7 Gy", profiles$radiation),
                        ifelse(profiles$genotype == "wt", "WT", "p21−"), sep = "\n")
paths <- file.path(primary, "shrunken", paste0(profiles$id, "_vs_control_ashr.csv"))
if (any(!file.exists(paths))) stop("Missing Model-A ashr CSV files:\n", paste(paths[!file.exists(paths)], collapse = "\n"))
read_result <- function(path) {
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("ENSEMBL", "SYMBOL", "padj", "log2FoldChange") %in% names(x))) stop("Invalid input: ", path)
  if (anyNA(x$ENSEMBL) || anyDuplicated(x$ENSEMBL)) stop("Non-unique ENSEMBL IDs: ", path)
  if ("lfc_method" %in% names(x) && any(x$lfc_method != "ashr_shrunken", na.rm = TRUE)) stop("Not ashr: ", path)
  x
}
## ******
## Join the twelve gene-wise results by stable ENSEMBL IDs
tabs <- lapply(paths, read_result)
ids <- tabs[[1]]$ENSEMBL
if (!all(vapply(tabs, function(t) setequal(t$ENSEMBL, ids), logical(1)))) stop("Incompatible gene universes.")
mat_field <- function(field) vapply(tabs, function(t) as.numeric(t[[field]][match(ids, t$ENSEMBL)]), numeric(length(ids)))
lfc <- mat_field("log2FoldChange")
padj <- mat_field("padj")
rownames(lfc) <- rownames(padj) <- ids
colnames(lfc) <- colnames(padj) <- profiles$id
symbol <- tabs[[1]]$SYMBOL[match(ids, tabs[[1]]$ENSEMBL)]

## ******
## Define the expression-filtered set and the union of responsive genes
processed <- readRDS(here::here("Data", "processed_data.rds"))
if (is.null(processed$dds) || is.null(processed$annotation)) stop("Missing dds / annotation in processed object.")
anno <- as.data.frame(processed$annotation)
if (!all(c("ENSEMBL", "ENTREZID") %in% names(anno))) stop("processed$annotation needs ENSEMBL and ENTREZID for GO-BP.")
anno <- anno[match(ids, anno$ENSEMBL), , drop = FALSE]
if (anyNA(anno$ENSEMBL)) stop("Missing gene annotation entries.")
entrez <- as.character(anno$ENTREZID)
if ("GENENAME" %in% names(anno)) gene_name <- as.character(anno$GENENAME) else gene_name <- rep("", length(ids))
cts <- DESeq2::counts(processed$dds, normalized = TRUE)
if (!setequal(rownames(cts), ids)) stop("Counts differ from script 02 gene universe.")
expressed <- rowSums(cts[match(ids, rownames(cts)), , drop = FALSE] >= 10) >= 5
complete <- rowSums(is.finite(lfc)) == ncol(lfc)
sig_union <- rowSums(is.finite(padj) & padj < alpha) > 0
nonzero <- rowSums(lfc^2, na.rm = TRUE) > 1e-14
mod_ids <- ids[sig_union & complete & nonzero]
exp_ids <- ids[expressed & complete]
if (length(mod_ids) < 100 || length(exp_ids) < 100) stop("Unexpectedly small response gene sets.")
mod <- lfc[mod_ids, , drop = FALSE]
utils::write.csv(data.frame(ENSEMBL = mod_ids, SYMBOL = symbol[match(mod_ids, ids)], mod, check.names = FALSE),
                 file.path(root, "response_matrix_DEG_union_ashr.csv"), row.names = FALSE)
utils::write.csv(data.frame(metric = c("DEG_union", "adequately_expressed", "L2_modules"),
                            n_genes = c(sum(sig_union), length(exp_ids), nrow(mod))),
                 file.path(root, "analysis_inventory.csv"), row.names = FALSE)

## ************************************************************************************************
## Response PCA: shared projection for both cell lines
## Each response (not each biological sample) is one PCA observation.
## Gene-center, but do NOT variance-scale the twelve response vectors.
## DEG union: PC1/2, PC1/3, PC2/3; expressed genes: PC1/2 for robustness.
status("2/5 | Response PCAs and PC1-PC4 gene loadings")
rad_col <- c("7gy" = "#3C5488", UVA = "#E64B35", UVB = "#00A087")
run_pca <- function(x, label, stem, extra_planes = FALSE) {
  pc <- stats::prcomp(t(x), center = TRUE, scale. = FALSE)
  pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  coord <- cbind(profiles, as.data.frame(pc$x[, paste0("PC", 1:4), drop = FALSE]))
  utils::write.csv(coord, file.path(pca_dir, paste0(stem, "_coordinates.csv")), row.names = FALSE)
  load <- data.frame(ENSEMBL = rownames(pc$rotation), SYMBOL = symbol[match(rownames(pc$rotation), ids)],
                     pc$rotation[, paste0("PC", 1:4), drop = FALSE], check.names = FALSE)
  utils::write.csv(load, file.path(pca_dir, paste0(stem, "_gene_loadings_PC1_PC4.csv")), row.names = FALSE)
  utils::write.csv(data.frame(PC = seq_along(pct), variance_percent = pct),
                   file.path(pca_dir, paste0(stem, "_variance.csv")), row.names = FALSE)
  coord$point_type <- factor(paste(coord$cell_line, coord$genotype, sep = "_"),
                             levels = c("HT1080_wt", "HT1080_p21_neg", "MCF10A_wt", "MCF10A_p21_neg"))
  wt <- coord[coord$genotype == "wt", ]; ko <- coord[coord$genotype == "p21_neg", ]
  planes <- if (extra_planes) list(c(1, 2), c(1, 3), c(2, 3)) else list(c(1, 2))
  for (plane in planes) {
    xp <- paste0("PC", plane[1]); yp <- paste0("PC", plane[2])
    lines <- data.frame(x = wt[[xp]], y = wt[[yp]], xend = ko[[xp]][match(wt$pair, ko$pair)],
                        yend = ko[[yp]][match(wt$pair, ko$pair)], radiation = wt$radiation)
    p <- ggplot2::ggplot(coord, ggplot2::aes(x = .data[[xp]], y = .data[[yp]])) +
      ggplot2::geom_segment(data = lines, ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = radiation),
                            inherit.aes = FALSE, linewidth = 0.75, alpha = 0.65,
                            arrow = grid::arrow(length = grid::unit(0.10, "cm"))) +
      ggplot2::geom_point(ggplot2::aes(color = radiation, shape = point_type), size = 3.8, stroke = 1.0) +
      ggplot2::scale_color_manual(values = rad_col, breaks = c("7gy", "UVA", "UVB"),
                                  labels = c("7 Gy", "UVA", "UVB"), name = "Irradiation") +
      ggplot2::scale_shape_manual(values = c(1, 16, 2, 17),
          labels = c("HT1080 WT", "HT1080 p21−", "MCF10A WT", "MCF10A p21−"), name = "Cell line / genotype") +
      ggplot2::coord_fixed() + plot_theme(12) + ggplot2::theme(legend.position = "bottom", legend.box = "vertical") +
      ggplot2::labs(title = paste0("Irradiation response PCA | ", label),
        subtitle = "One shared PCA; WT → p21− transitions within each irradiation and cell line",
        x = sprintf("%s (%.1f%%)", xp, pct[plane[1]]), y = sprintf("%s (%.1f%%)", yp, pct[plane[2]]),
        caption = sprintf("%s genes | Model-A ashr log2FC; gene-centered, not variance-scaled", nrow(x)))
    export(file.path(pca_dir, paste0(stem, "_", xp, "_", yp)), p, 9.0, 7.0)
  }
  pc
}
## Both PCAs use the same twelve responses; only the gene universe changes.
pca_deg <- run_pca(mod, "DEG union (primary)", "PCA_DEG_union", TRUE)
run_pca(lfc[exp_ids, , drop = FALSE], "all adequately expressed genes (robustness)", "PCA_all_expressed")

## ************************************************************************************************
## L2-normalized gene modules
## L2 removes overall vector magnitude during clustering, but preserves sign
## and relative response pattern. All output heatmaps use ORIGINAL ashr log2FC.
## Eight modules are exploratory clusters, not discrete validated programs.
status("3/5 | Clustering ", nrow(mod), " responsive genes into ", k, " L2 modules")
mag <- sqrt(rowSums(mod^2))
normalized <- sweep(mod, 1, mag, "/")
set.seed(21)
fit <- stats::kmeans(normalized, centers = k, nstart = 40, iter.max = 150)
cent <- t(vapply(seq_len(k), function(i) colMeans(mod[fit$cluster == i, , drop = FALSE]), numeric(ncol(mod))))
ord <- order(rowMeans(cent), decreasing = TRUE)
module_ids <- sprintf("M%02d", seq_len(k))
module <- module_ids[match(fit$cluster, ord)]
cent <- cent[ord, , drop = FALSE]
rownames(cent) <- module_ids
colnames(cent) <- profiles$id  # keep explicit condition IDs in the exported mean-response CSV
names(module) <- rownames(mod)
unit_cent <- t(vapply(module_ids, function(m) {
  v <- colMeans(normalized[module == m, , drop = FALSE]); v / sqrt(sum(v^2))
}, numeric(ncol(mod))))
cosine <- rowSums(normalized * unit_cent[match(module, module_ids), , drop = FALSE])
mem <- data.frame(ENSEMBL = mod_ids, SYMBOL = symbol[match(mod_ids, ids)], ENTREZID = entrez[match(mod_ids, ids)],
                  GENENAME = gene_name[match(mod_ids, ids)], module = unname(module),
                  profile_cosine = cosine, response_L2_magnitude = mag, stringsAsFactors = FALSE)
utils::write.csv(mem, file.path(mod_dir, "module_membership.csv"), row.names = FALSE)
## ******
## Relate gene modules to PC1–PC4 through their squared loadings.
## Squared loading contributions sum to one within each PC.
loadings <- pca_deg$rotation[mod_ids, paste0("PC", 1:4), drop = FALSE]
contrib <- do.call(rbind, lapply(module_ids, function(m) {
  i <- module == m; frac <- colSums(loadings[i, , drop = FALSE]^2)
  data.frame(module = m, n_genes = sum(i), PC = names(frac),
             fraction_squared_loadings = as.numeric(frac),
             enrichment_over_module_size = as.numeric(frac) / mean(i))
}))
utils::write.csv(contrib, file.path(mod_dir, "module_PCA_contribution_PC1_PC4.csv"), row.names = FALSE)
utils::write.csv(data.frame(module = module_ids, n_genes = as.integer(table(factor(module, levels = module_ids))),
                            cent, check.names = FALSE), file.path(mod_dir, "module_mean_ashr_LFC.csv"), row.names = FALSE)

## ******
## Descriptive response labels (not inferred GO or mechanism).
## Retain stable module IDs in every exported table.
desc <- c("HT1080 UVA-associated induction", "MCF10A 7 Gy/UVA induction", "HT1080 UVB induction",
          "MCF10A UVB induction", "HT1080 UVB repression", "MCF10A 7 Gy/UVA repression",
          "MCF10A UVB repression", "HT1080 UVA-associated repression")
module_info <- data.frame(module = module_ids, descriptive_pattern = desc,
                          n_genes = vapply(module_ids, function(m) sum(module == m), integer(1)))
utils::write.csv(module_info, file.path(mod_dir, "module_descriptive_labels.csv"), row.names = FALSE)
## Module-level WT → p21− differences are descriptive, NOT interaction tests.
delta_rows <- list()
for (m in module_ids) {
  for (cl in c("HT1080", "MCF10A")) {
    for (rad in c("7gy", "UVA", "UVB")) {
      ix <- which(profiles$cell_line == cl & profiles$radiation == rad)
      wt_idx <- ix[profiles$genotype[ix] == "wt"]
      ko_idx <- ix[profiles$genotype[ix] == "p21_neg"]
      delta_rows[[length(delta_rows) + 1L]] <- data.frame(
        module = m, cell_line = cl, irradiation = rad,
        WT = cent[m, wt_idx], p21_negative = cent[m, ko_idx],
        KO_minus_WT = cent[m, ko_idx] - cent[m, wt_idx]
      )
    }
  }
}
delta <- do.call(rbind, delta_rows)
utils::write.csv(delta, file.path(mod_dir, "module_genotype_deltas_descriptive.csv"), row.names = FALSE)

## ************************************************************************************************
## GO Biological Process: module over-representation (ORA)
## Universe: genes in the responsive DEG union with an ENTREZ ID, not the genome.
## BH correction is applied within each module; modules are not independent.
##
## Runtime: simplify() can be expensive if it receives all GO terms. We keep
## all enrichment results, but simplify ONLY FDR-significant terms for display.
## Semantic data are cached and progress/timing is printed for each module.
status("4/5 | GO-BP ORA; ", length(module_ids), " modules, with per-module progress")
valid_id <- function(x) !is.na(x) & nzchar(x) & x != "NA"
universe <- unique(mem$ENTREZID[valid_id(mem$ENTREZID)])
go_all <- list()
go_top <- list()
go_nonredundant <- list()  # Complete significant set after semantic redundancy reduction
go_semantic_success <- setNames(rep(FALSE, length(module_ids)), module_ids)
go_members <- setNames(vector("list", k), module_ids)
sem_data <- NULL                             # cached ontology graph for GO-term similarity
sem_data_attempted <- FALSE                 # do not retry if unavailable

for (i in seq_along(module_ids)) {
  m <- module_ids[i]
  keys <- unique(mem$ENTREZID[mem$module == m & valid_id(mem$ENTREZID)])
  status(sprintf("GO-BP %d/%d | %s | %d ENTREZ genes | enrichGO started", i, k, m, length(keys)))

  t0 <- proc.time()[["elapsed"]]
  go <- tryCatch(
    clusterProfiler::enrichGO(gene = keys, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID", ont = "BP", universe = universe, pAdjustMethod = "BH",
      pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 15, maxGSSize = 500, readable = FALSE),
    error = function(e) { warning("GO-BP ", m, ": ", conditionMessage(e)); NULL }
  )
  tab <- if (is.null(go)) data.frame() else as.data.frame(go)
  n_sig <- if (nrow(tab)) sum(is.finite(tab$p.adjust) & tab$p.adjust < alpha) else 0L
  status(sprintf("GO-BP %d/%d | %s | enrichGO done in %.1fs | %d terms, %d FDR<0.05",
                 i, k, m, proc.time()[["elapsed"]] - t0, nrow(tab), n_sig))
  if (!nrow(tab)) {
    go_semantic_success[[m]] <- TRUE
    next
  }

  tab$module <- m
  go_all[[m]] <- tab                        # retain ALL tested terms in the CSV
  sig <- tab[is.finite(tab$p.adjust) & tab$p.adjust < alpha, , drop = FALSE]
  go_semantic_success[[m]] <- nrow(sig) <= 1L

  ## Semantic reduction supplies a separate, auditable representative term set; full ORA is retained.
  ## Shared semantic data and FDR filtering limit expensive pairwise comparisons.
  if (nrow(sig) > 1L) {
    if (!sem_data_attempted) {
      status("GO-BP | Preparing shared GO semantic data (once)...")
      sem_data <- tryCatch(
        GOSemSim::godata("org.Hs.eg.db", ont = "BP", computeIC = FALSE),
        error = function(e) { warning("GO semantic data: ", conditionMessage(e)); NULL }
      )
      sem_data_attempted <- TRUE
    }
    if (!is.null(sem_data)) {
      status(sprintf("GO-BP %d/%d | %s | simplifying %d significant terms...", i, k, m, nrow(sig)))
      t1 <- proc.time()[["elapsed"]]
      go_sig <- go
      go_sig@result <- go@result[go@result$ID %in% sig$ID, , drop = FALSE]
      reduced <- tryCatch(
        clusterProfiler::simplify(go_sig, cutoff = 0.75, by = "p.adjust",
                                  select_fun = min, semData = sem_data),
        error = function(e) { warning("GO-BP simplify ", m, ": ", conditionMessage(e)); NULL }
      )
      if (!is.null(reduced)) {
        sig <- sig[sig$ID %in% as.data.frame(reduced)$ID, , drop = FALSE]
        go_semantic_success[[m]] <- TRUE
      }
      status(sprintf("GO-BP %d/%d | %s | simplify done in %.1fs | %d representative terms",
                     i, k, m, proc.time()[["elapsed"]] - t1, nrow(sig)))
    }
  }

  sig <- sig[order(sig$p.adjust, -sig$Count), , drop = FALSE]
  go_nonredundant[[m]] <- sig
  go_top[[m]] <- head(sig, 3L)
  go_members[[m]] <- if (nrow(sig)) {
    unlist(strsplit(paste(head(sig$geneID, 5L), collapse = "/"), "/", fixed = TRUE))
  } else character()
  status(sprintf("GO-BP %d/%d | %s finished | %d GO terms selected for figure", i, k, m, nrow(go_top[[m]])))
}
status("GO-BP | All eight modules completed; writing enrichment tables")
combined <- if (length(go_all)) do.call(rbind, go_all) else data.frame()
utils::write.csv(combined, file.path(mod_dir, "GO_BP_ORA_all_terms.csv"), row.names = FALSE)
chosen <- if (length(go_top) && any(lengths(go_top) > 0)) do.call(rbind, go_top[lengths(go_top) > 0]) else data.frame()
utils::write.csv(chosen, file.path(mod_dir, "GO_BP_ORA_top3_per_module.csv"), row.names = FALSE)
## Preserve all FDR-significant GO terms after Wang semantic simplification.
## Theme summaries cannot reliably use the original raw GO list (heavy redundancy).
nonredundant <- if (length(go_nonredundant)) do.call(rbind, go_nonredundant) else data.frame()
utils::write.csv(nonredundant, file.path(mod_dir, "GO_BP_ORA_nonredundant_terms.csv"), row.names = FALSE)
utils::write.csv(data.frame(module = module_ids,
  wang_simplification_success = unname(go_semantic_success)),
  file.path(mod_dir, "GO_BP_semantic_reduction_status.csv"), row.names = FALSE)
if (any(!go_semantic_success)) stop(
  "GO semantic reduction incomplete for: ",
  paste(names(go_semantic_success)[!go_semantic_success], collapse = ", "),
  ". See GO_BP_semantic_reduction_status.csv; do not classify unreduced terms.")
utils::write.csv(data.frame(module = module_ids, n_significant_GO_BP = vapply(module_ids, function(m)
  if (is.null(go_all[[m]])) 0L else sum(go_all[[m]]$p.adjust < 0.05, na.rm = TRUE), integer(1)),
  n_top_terms = vapply(module_ids, function(m) if (is.null(go_top[[m]])) 0L else nrow(go_top[[m]]), integer(1))),
  file.path(mod_dir, "GO_BP_ORA_summary.csv"), row.names = FALSE)

## ****
## Reader-friendly annotation: one row per module, retaining descriptive
## RESPONSE patterns separately from enriched GO functions. GO enrichment
## alone cannot determine up-/downregulation of a biological process.
module_go <- module_info
module_go$n_significant_GO_BP <- vapply(module_ids, function(m) {
  x <- go_all[[m]]
  if (is.null(x)) 0L else sum(is.finite(x$p.adjust) & x$p.adjust < alpha)
}, integer(1))
for (j in seq_len(3L)) {
  module_go[[paste0("GO_BP_", j)]] <- vapply(module_ids, function(m) {
    x <- go_top[[m]]
    if (is.null(x) || nrow(x) < j) NA_character_ else as.character(x$Description[j])
  }, character(1))
  module_go[[paste0("GO_FDR_", j)]] <- vapply(module_ids, function(m) {
    x <- go_top[[m]]
    if (is.null(x) || nrow(x) < j) NA_real_ else as.numeric(x$p.adjust[j])
  }, numeric(1))
}
utils::write.csv(module_go, file.path(mod_dir, "module_GO_BP_overview.csv"), row.names = FALSE)

## ************************************************************************************************
## Select representative genes for the module heatmap
## Five genes best matching the module profile plus five curated or GO-supported
## genes from the SAME module. Export the reason behind every selection.
status("5/5 | Selecting representative and biologically informative genes")
curated <- c("CDKN1A", "MDM2", "GADD45A", "DDB2", "BTG2", "ATM", "BAX", "IL6", "IL1B",
             "INHBA", "SERPINB2", "CD44", "JUN", "FOS", "CXCL8", "NFKBIA", "TNFAIP3",
             "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "GINS2", "DBF4",
             "PCNA", "CDC6", "CDK2", "CCNB1", "RRM2", "DDIT3", "GADD45B", "NQO1", "HMOX1")
rep_list <- lapply(module_ids, function(m) {
  x <- mem[mem$module == m, , drop = FALSE]
  x <- x[order(-x$profile_cosine, -x$response_L2_magnitude), ]
  chosen_stat <- head(x, n_pattern); chosen_stat$selection_reason <- "pattern_representative"
  x <- x[!x$ENSEMBL %in% chosen_stat$ENSEMBL, , drop = FALSE]
  x$curated_rank <- match(x$SYMBOL, curated)
  x$GO_support <- as.integer(x$ENTREZID %in% go_members[[m]])
  x$has_symbol <- valid_id(x$SYMBOL) & !grepl("^ENSG", x$SYMBOL)
  x <- x[order(is.na(x$curated_rank), x$curated_rank, -x$GO_support,
               -x$has_symbol, -x$response_L2_magnitude, -x$profile_cosine), , drop = FALSE]
  chosen_bio <- head(x, n_biological)
  chosen_bio$selection_reason <- ifelse(!is.na(chosen_bio$curated_rank), "curated_known_gene",
                              ifelse(chosen_bio$GO_support > 0, "enriched_GO_term_gene", "annotated_high_effect"))
  if (nrow(chosen_bio) < n_biological) stop("Not enough genes in module ", m)
  rbind(chosen_stat[, c(names(mem), "selection_reason")],
        chosen_bio[, c(names(mem), "selection_reason")])
})
representative <- do.call(rbind, rep_list)
representative$gene_label <- ifelse(valid_id(representative$SYMBOL), representative$SYMBOL, representative$ENSEMBL)
representative$gene_label <- make.unique(representative$gene_label)
utils::write.csv(representative, file.path(mod_dir, "representative_genes_with_selection_reason.csv"), row.names = FALSE)
## Localized Excel-compatible export: GENENAME can contain commas.
## Preserve the standard quoted CSV for R and other machine-readable tools.
stopifnot(nrow(representative) == 80L,
          all(table(factor(representative$module, levels = module_ids)) == 10L),
          !anyDuplicated(representative$ENSEMBL))
utils::write.csv2(representative,
  file.path(mod_dir, "representative_genes_with_selection_reason_excel_DE.csv"), row.names = FALSE)

## Final module heatmaps and GO-theme figures are generated in
## R/03_final_visualisations.R from the exported tables.

## ************************************************************************************************
## Sensitivity analysis: processing-matched irradiation contrasts
## EKB: UVB vs EKB controls (four genotype/cell-line contexts).
## UKD: 7 Gy/UVA vs UKD controls (three contexts; MCF10A WT unavailable).
## Refit local raw counts with local dispersion/size-factor estimation.
## This does NOT test a processing × irradiation interaction; HT1080 WT 7 Gy
## also spans sequencing runs and is therefore not fully run-matched.
status("Sensitivity | Refitting within-site DESeq2 models for ten sensitivity contrasts")
meta <- as.data.frame(SummarizedExperiment::colData(processed$dds))
if (!all(c("cell_line", "genotype", "radiation", "processing", "sequencing_run", "group") %in% names(meta)))
  stop("Missing sample metadata for matched-control sensitivity analysis.")
sens_rows <- list()
for (site in c("EKB", "UKD")) {
  status("Processing sensitivity | ", site, " | fitting local DESeq2 model...")
  keep <- as.character(meta$processing) == site
  if (site == "UKD") keep <- keep & !(meta$cell_line == "MCF10A" & meta$genotype == "wt")
  site_meta <- meta[keep, , drop = FALSE]
  site_meta$local_group <- droplevels(factor(as.character(site_meta$group)))
  if (qr(stats::model.matrix(~ local_group, site_meta))$rank != nlevels(site_meta$local_group))
    stop("Rank-deficient local model for processing site ", site)
  local_dds <- DESeq2::DESeqDataSetFromMatrix(countData = DESeq2::counts(processed$dds)[, keep, drop = FALSE],
    colData = site_meta, design = ~ local_group)
  local_dds <- DESeq2::DESeq(local_dds, quiet = TRUE)
  status("Processing sensitivity | ", site, " | model fitted; comparing treatment effects")
  subset_plan <- profiles[profiles$radiation %in% if (site == "EKB") "UVB" else c("7gy", "UVA"), ]
  if (site == "UKD") subset_plan <- subset_plan[!(subset_plan$cell_line == "MCF10A" & subset_plan$genotype == "wt"), ]
  for (j in seq_len(nrow(subset_plan))) {
    r <- subset_plan[j, ]
    treatment <- group_name(r$cell_line, r$genotype, r$radiation)
    control <- group_name(r$cell_line, r$genotype, "control")
    if (!all(c(treatment, control) %in% levels(site_meta$local_group))) next
    local <- DESeq2::results(local_dds, contrast = c("local_group", treatment, control), alpha = alpha)
    reference_path <- file.path(primary, "unshrunken", paste0(r$id, "_vs_control_unshrunken.csv"))
    if (!file.exists(reference_path)) stop("Missing reference unshrunken file: ", reference_path)
    ref <- utils::read.csv(reference_path, stringsAsFactors = FALSE)
    ref_lfc <- ref$log2FoldChange[match(rownames(local), ref$ENSEMBL)]
    ok <- is.finite(local$log2FoldChange) & is.finite(ref_lfc) & expressed[match(rownames(local), ids)]
    sample_treat <- site_meta[as.character(site_meta$group) == treatment, , drop = FALSE]
    sample_ctrl <- site_meta[as.character(site_meta$group) == control, , drop = FALSE]
    matched_runs <- intersect(unique(as.character(sample_treat$sequencing_run)), unique(as.character(sample_ctrl$sequencing_run)))
    sens_rows[[length(sens_rows) + 1]] <- data.frame(contrast = r$id, processing = site,
      n_treatment = nrow(sample_treat), n_control = nrow(sample_ctrl),
      treatment_runs = paste(unique(as.character(sample_treat$sequencing_run)), collapse = ";"),
      control_runs = paste(unique(as.character(sample_ctrl$sequencing_run)), collapse = ";"),
      all_samples_run_matched = setequal(unique(as.character(sample_treat$sequencing_run)),
                                        unique(as.character(sample_ctrl$sequencing_run))),
      overlapping_runs = paste(matched_runs, collapse = ";"),
      n_comparable_genes = sum(ok), Pearson_r = if (sum(ok) > 2) stats::cor(local$log2FoldChange[ok], ref_lfc[ok]) else NA_real_,
      median_abs_delta_log2FC = stats::median(abs(local$log2FoldChange[ok] - ref_lfc[ok])),
      pct_same_direction = 100 * mean(sign(local$log2FoldChange[ok]) == sign(ref_lfc[ok])))
  }
}
sensitivity <- do.call(rbind, sens_rows)
utils::write.csv(sensitivity, file.path(sens_dir, "processing_matched_vs_primary_LFC_summary.csv"), row.names = FALSE)

## ******
## Compact sensitivity figure: every matched-site contrast is shown once.
## The comparison uses unshrunken Wald LFC; refitting changes size factors
## and dispersion estimates as well as the control composition. A local
## refit is NOT evidence for a processing × treatment interaction.
sensitivity$cell_line <- profiles$cell_line[match(sensitivity$contrast, profiles$id)]
sensitivity$genotype <- profiles$genotype[match(sensitivity$contrast, profiles$id)]
sensitivity$radiation <- profiles$radiation[match(sensitivity$contrast, profiles$id)]
sensitivity <- sensitivity[order(match(sensitivity$cell_line, c("HT1080", "MCF10A")),
                                 match(sensitivity$radiation, c("7gy", "UVA", "UVB")),
                                 match(sensitivity$genotype, c("wt", "p21_neg"))), , drop = FALSE]
sensitivity$label <- with(sensitivity, paste(cell_line, ifelse(genotype == "wt", "WT", "p21−"),
                                             ifelse(radiation == "7gy", "7 Gy", radiation), sep = " | "))
sensitivity$label <- factor(sensitivity$label, levels = rev(sensitivity$label))
sensitivity$run_context <- ifelse(sensitivity$all_samples_run_matched, "Run-matched", "Mixed runs")

## Each contrast appears once: the x-axis quantifies LFC deviation, while
## the Pearson coefficient is printed next to its point (not encoded as color).
p_sensitivity <- ggplot2::ggplot(sensitivity,
                                ggplot2::aes(x = median_abs_delta_log2FC, y = label)) +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = median_abs_delta_log2FC, yend = label),
                        color = "grey78", linewidth = 0.6) +
  ggplot2::geom_point(ggplot2::aes(color = processing, shape = run_context), size = 3.2, stroke = 0.85) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("r = %.3f", Pearson_r)),
                     hjust = -0.15, size = 3.2, color = "grey20", show.legend = FALSE) +
  ggplot2::scale_color_manual(values = c(EKB = "#3C5488", UKD = "#00A087"), name = "Processing") +
  ggplot2::scale_shape_manual(values = c("Run-matched" = 16, "Mixed runs" = 2), name = "Sequencing runs") +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.30)),
                              labels = scales::label_number(accuracy = 0.02)) +
  plot_theme(11) +
  ggplot2::labs(title = "Processing-matched control sensitivity",
                subtitle = sprintf("Primary Model A versus processing-specific refits | %s expressed genes/contrast",
                                   scales::comma(stats::median(sensitivity$n_comparable_genes))),
                x = "Median absolute difference in unshrunken log2FC", y = NULL,
                caption = "Values next to points are Pearson r. HT1080 WT 7 Gy spans two sequencing runs; MCF10A WT UVA/7 Gy cannot be checked.") +
  ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
                 legend.position = "bottom", legend.box = "horizontal",
                 axis.line = ggplot2::element_line(color = "black", linewidth = 0.55),
                 axis.text.y = ggplot2::element_text(size = 10),
                 plot.caption = ggplot2::element_text(size = 9, color = "grey35"))
export(file.path(sens_dir, "processing_matched_vs_primary_LFC_diagnostic"), p_sensitivity, 11.8, 7.0)
writeLines(c("Sensitivity only: contrasts refit within processing site, using only matching-site controls.",
 "EKB: four UVB contrasts; UKD: 7 Gy/UVA in HT1080 WT, HT1080 p21-negative and MCF10A p21-negative.",
 "MCF10A WT 7 Gy/UVA is excluded because there are no UKD WT control samples.",
 "HT1080 WT 7 Gy treatment samples span sequencing runs; therefore this is NOT fully sequencing-run matched.",
 "Refitting changes size-factor and dispersion estimation as well as the control composition.",
 "Model A remains the biological primary analysis. Do not promote this sensitivity subset to a new DEG set."),
 file.path(sens_dir, "INTERPRETATION_NOTES.txt"))
## ************************************************************************************************
## Summary: keep Model A as the sole biological analysis; processing-matched
## refits are diagnostics only. Formal interaction contrasts follow in script 04.
status("Complete | ", nrow(mod), " response genes; eight L2 modules; outputs in ", root)

## Run the next two stages only after stage 1 completed successfully.
## Each can also be rerun independently without recalculating these results.
message("Stage 1/3 complete: primary response analysis + GO semantic reduction")
source(here::here("R", "03_cluster_evaluation.R"), local = TRUE)
source(here::here("R", "03_final_visualisations.R"), local = TRUE)
message("Script 03 complete: analysis, cluster evaluation and final figures/audit")
