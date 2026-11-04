## ************************************************************************************************
## p21 bulk RNA-seq
## 06 - GO functional landscape: irradiation-response modules and p21 interactions
##
## Biological questions:
##   Which GO Biological Processes characterize the eight L2 response modules?
##   Which GO processes are positively/negatively enriched across the six
##   formal p21 x irradiation interactions (KO response minus WT response)?
##   Which broad functional themes connect the two levels (with recorded GO IDs)?
##
## Analysis decisions:
##   - REUSE script 03's full GO-BP ORA; do not recompute module enrichment.
##   - GSEA ranks ALL testable, mapped genes by the SIGNED UN-SHRUNKEN Wald
##     statistic from each of script 04's six formal interaction tests.
##   - Collapse duplicate ENTREZ IDs by highest baseMean (not maximum z).
##   - BH FDR < 0.05 is computed separately within each GSEA contrast.
##   - Use Wang GO-semantic redundancy reduction for FIGURES only;
##     export ALL tested ORA/GSEA terms. No statistical tests on term overlap.
##   - Shared landscape shows six THEMES with separate representative GO IDs
##     for ORA and GSEA (no claim that the GO IDs necessarily match).
##   - A shared GO term is descriptive concordance, not independent validation:
##     underlying RNA-seq data overlap, and ORA/GSEA have different nulls.
##   - IPA/DREAM are separate later analyses; no new fit, clusters or RDS.
##
## Output (all below Results/06_functional_GO/):
##   tables/GO_BP_GSEA_all_terms.csv: full six-contrast GSEA results
##   tables/GO_BP_ORA_display_terms.csv: curated-by-rule module terms for plot
##   tables/GO_BP_GSEA_display_terms.csv: GSEA terms for plot
##   tables/functional_landscape_terms.csv: six-theme matrix with GO-ID provenance
##   tables/functional_landscape_theme_definitions.csv: ORA/GSEA GO IDs per theme
##   tables/GO_module_leading_edge_overlap.csv: overlapping ENTREZ IDs/symbols
##   figures/GO_BP_module_ORA_dotplot.{png,svg}
##   figures/GO_BP_interaction_GSEA_heatmap.{png,svg}
##   figures/GO_BP_functional_landscape.{png,svg}
##   Existing six-contrast GO GSEA can be reused for display-only revisions.
## ************************************************************************************************

rm(list = ls())
gc()

## ******
## Project setup, dependencies and paths
if (!requireNamespace("here", quietly = TRUE)) stop("Install the 'here' package.")
here::i_am("R/06_functional_GO.R")
source(here::here("R", "00_setup.R"))
source(here::here("R", "plot_style.R"))
for (pkg in c("ggplot2", "scales", "clusterProfiler", "org.Hs.eg.db", "fgsea",
              "GOSemSim", "ekbSeq", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing R package: ", pkg)
}

alpha <- 0.05
## fgsea reports p-values below this bound as 1e-10. Setting this to 0
## resolves more extreme p-values, at the cost of potentially longer runs.
## The default is retained here for consistency with the initial analysis.
gsea_eps <- 1e-10
## Reuse the already computed complete GSEA CSV when updating figures only.
## Set FALSE if input contrasts, annotation, or GSEA parameters have changed.
reuse_existing_gsea <- TRUE
semantic_cutoff <- 0.75   # Wang GO semantic similarity: representative display only
pair_order <- c("M01", "M08", "M02", "M06", "M03", "M05", "M04", "M07")
root <- here::here("Results", "06_functional_GO")
tab_dir <- file.path(root, "tables")
fig_dir <- file.path(root, "figures")
mod_dir <- here::here("Results", "03_response_architecture", "02_MODULE_L2")
raw_dir <- here::here("Results", "04_p21_interactions", "tables", "unshrunken")
for (d in c(tab_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
status <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))
  flush.console()
}
plot_theme <- function(size = 10) {
  if (exists("p21_theme", mode = "function")) p21_theme(size) else ggplot2::theme_classic(base_size = size)
}
export <- function(name, p, w, h) ekbSeq::export_plot_dual(
  file.path(fig_dir, name), p, width = w, height = h, rasterize = FALSE)

## ******
## Small, deterministic helpers: preserve IDs and representativeness.
valid_entrez <- function(x) !is.na(x) & nzchar(x) & x != "NA"
split_ids <- function(x) {
  if (length(x) == 0L || is.na(x) || !nzchar(x)) return(character())
  unique(strsplit(as.character(x), "/", fixed = TRUE)[[1]])
}
## Ensure portable, readable labels. GO ID is always kept as the actual key.
short_name <- function(x, n = 57L) ifelse(nchar(x) > n, paste0(substr(x, 1L, n - 3L), "..."), x)

## ******
## Semantic GO pruning for FIGURES only. Never remove terms from full results.
## GOSemSim compares GO ontology meaning rather than overlapping gene lists.
## Limit candidate pools to avoid expensive all-vs-all GO similarity matrices.
status("Preparing GO-BP semantic similarity (Wang; reused throughout)")
sem_data <- GOSemSim::godata(ont = "BP", computeIC = FALSE)

select_semantic <- function(tab, n = 3L, anchors = character(), max_candidates = 65L) {
  if (!nrow(tab)) return(tab)
  stopifnot(all(c("ID", "p.adjust") %in% names(tab)))
  tab <- tab[order(tab$p.adjust, tab$ID), , drop = FALSE]
  tab <- tab[!duplicated(tab$ID), , drop = FALSE]
  favored <- tab[tab$ID %in% anchors, , drop = FALSE]
  if (nrow(favored)) favored <- favored[order(match(favored$ID, anchors)), , drop = FALSE]
  rest <- tab[!tab$ID %in% favored$ID, , drop = FALSE]
  candidates <- rbind(favored, head(rest, max_candidates))
  if (nrow(candidates) == 1L) return(candidates)
  sim <- GOSemSim::mgoSim(candidates$ID, candidates$ID,
                          semData = sem_data, measure = "Wang", combine = NULL)
  selected <- integer()
  for (k in seq_len(nrow(candidates))) {
    ## NA similarities are not evidence of redundancy; keep the term.
    related <- if (length(selected)) sim[k, selected] > semantic_cutoff else FALSE
    if (!any(related, na.rm = TRUE)) selected <- c(selected, k)
    if (length(selected) >= n) break
  }
  candidates[selected, , drop = FALSE]
}

## Explicit biological categories only for the FIGURE, not new hypothesis tests.
## The theme view may use different representative GO IDs on the ORA/GSEA sides;
## both exact IDs are recorded and never presented as the same statistical test.
theme_levels <- c("Cell cycle & DNA replication", "Mitochondrial energy", 
                  "Inflammation & cell communication", "Autophagy & lysosome",
                  "Cytoskeleton & adhesion", "Cilium organization")
classify_theme <- function(description) {
  d <- tolower(as.character(description))
  out <- rep(NA_character_, length(d))
  rules <- list(
    "Cell cycle & DNA replication" = "dna replication|chromosome|chromatid|cell cycle|mitotic|dna metabolic process",
    "Mitochondrial energy" = "mitochondri|respiration|respiratory|electron transport chain|atp synthesis|proton motive",
    "Inflammation & cell communication" = "inflamm|cytokine|chemotax|leukocyte|angiogen|cell communication|cell adhesion",
    "Autophagy & lysosome" = "autophag|lysosom|vacuol|vesicul",
    "Cytoskeleton & adhesion" = "actin|cytoskeleton|gtpase|cell.junction|focal adhesion",
    "Cilium organization" = "cili(um|a|ary)|ciliogenesis"
  )
  for (theme in names(rules)) {
    take <- is.na(out) & !is.na(d) & grepl(rules[[theme]], d, perl = TRUE)
    out[take] <- theme
  }
  out
}

## One GO ID per method/theme across conditions; reproducibly selected based
## on the number of significant contexts and minimum FDR (not maximum NES).
representative_id <- function(tab, theme, anchors = character()) {
  x <- tab[tab$theme == theme & !is.na(tab$ID) & is.finite(tab$p.adjust), , drop = FALSE]
  if (!nrow(x)) return(NA_character_)
  x$sig <- x$p.adjust < alpha
  counts <- aggregate(x$sig, list(ID = x$ID), sum)
  names(counts)[2] <- "n_significant"
  best <- aggregate(x$p.adjust, list(ID = x$ID), min)
  names(best)[2] <- "min_FDR"
  candidates <- merge(counts, best, by = "ID")
  candidates <- candidates[order(-candidates$n_significant, candidates$min_FDR, candidates$ID), ]
  ## Anchor preference is limited to IDs that are significantly represented.
  hit_anchor <- intersect(anchors, candidates$ID[candidates$n_significant > 0L])
  if (length(hit_anchor)) return(hit_anchor[1])
  candidates$ID[1]
}

## ************************************************************************************************
## Reuse the full module GO-BP ORA from script 03
## Background = DEG-union genes with ENTREZ annotation (script-03 definition).
## Do NOT interpret GO-ORA enrichment as process activation or repression.
status("1/5 | Reading full module GO-BP ORA and fixed L2 membership")
ora_path <- file.path(mod_dir, "GO_BP_ORA_all_terms.csv")
mem_path <- file.path(mod_dir, "module_membership.csv")
for (path in c(ora_path, mem_path)) if (!file.exists(path)) stop("Missing script-03 output: ", path)
ora <- utils::read.csv(ora_path, stringsAsFactors = FALSE)
mem <- utils::read.csv(mem_path, stringsAsFactors = FALSE)
if (!all(c("ID", "Description", "p.adjust", "geneID", "GeneRatio", "Count", "module") %in% names(ora)) ||
    !all(c("ENSEMBL", "module") %in% names(mem)) || anyDuplicated(mem$ENSEMBL) ||
    !setequal(unique(mem$module), pair_order)) stop("Unexpected GO-ORA or membership columns.")
ora <- ora[ora$module %in% pair_order & !is.na(ora$ID), , drop = FALSE]
ora_sig <- ora[is.finite(ora$p.adjust) & ora$p.adjust < alpha, , drop = FALSE]
## Up to three representative processes per module; selected by semantically
## distinct GO IDs. Reviewed biological anchors take precedence only when
## enriched (FDR < 0.05); full enrichment results remain unchanged.
focus_go <- list(
  M01 = c("GO:0006914", "GO:0007040"),
  M08 = c("GO:0006260", "GO:0051276"),
  M02 = c("GO:0006954", "GO:0001525"),
  M06 = c("GO:0140053", "GO:0032543"),
  M03 = character(), M05 = c("GO:0006695"),
  M04 = c("GO:0007264", "GO:0007015"),
  M07 = c("GO:0060271", "GO:0044782")
)
ora_shown <- do.call(rbind, lapply(pair_order, function(m) {
  x <- ora_sig[ora_sig$module == m, , drop = FALSE]
  if (!nrow(x)) return(NULL)
  status("Module ORA | ", m, " | semantic selection among ", nrow(x), " significant GO terms")
  z <- select_semantic(x, n = 3L, anchors = focus_go[[m]], max_candidates = 60L)
  z$display_selection <- ifelse(z$ID %in% focus_go[[m]], "reviewed_anchor", "Wang_nonredundant_FDR")
  z
}))
if (is.null(ora_shown)) ora_shown <- ora_sig[FALSE, , drop = FALSE]
ora_shown <- as.data.frame(ora_shown)
if (nrow(ora_shown)) {
  ratio <- strsplit(as.character(ora_shown$GeneRatio), "/", fixed = TRUE)
  ora_shown$gene_ratio_percent <- vapply(ratio, function(z) if (length(z) == 2L)
    100 * as.numeric(z[1]) / as.numeric(z[2]) else NA_real_, numeric(1))
  ora_shown$minus_log10_FDR <- -log10(pmax(ora_shown$p.adjust, .Machine$double.xmin))
} else {
  ora_shown$gene_ratio_percent <- numeric(0)
  ora_shown$minus_log10_FDR <- numeric(0)
}
utils::write.csv(ora_shown, file.path(tab_dir, "GO_BP_ORA_display_terms.csv"), row.names = FALSE)
status("Module ORA | ", nrow(ora_sig), " significant results; ",
       nrow(ora_shown), " representative GO-module pairs selected")

## ************************************************************************************************
## Six GO-BP GSEAs of formal p21 x radiation Wald interaction statistics
## GSEA is performed on ALL testable genes (not the significant subset and not
## the DEG union). Wald statistic is signed KO-response minus WT-response.
## ENTREZ duplicates retain the most abundant feature (highest baseMean),
## not the greatest |stat|; this avoids selecting extreme Wald values.
status("2/5 | Loading or computing six GO-BP interaction GSEAs")
plan <- expand.grid(radiation = c("7gy", "UVA", "UVB"),
                    cell_line = c("HT1080", "MCF10A"), stringsAsFactors = FALSE)
plan$id <- with(plan, paste(cell_line, radiation, sep = "_"))
raw_paths <- file.path(raw_dir, paste0(plan$id, "_p21_interaction_unshrunken.csv"))
if (any(!file.exists(raw_paths))) stop("Missing script-04 unshrunken files:\n",
                                       paste(raw_paths[!file.exists(raw_paths)], collapse = "\n"))
full_gsea_path <- file.path(tab_dir, "GO_BP_GSEA_all_terms.csv")
annotation_map <- NULL
using_cached_gsea <- reuse_existing_gsea && file.exists(full_gsea_path)
if (using_cached_gsea) {
  ## Reuse a completed six-contrast analysis for FIGURE-ONLY revisions.
  ## The user must set reuse_existing_gsea = FALSE when changing the input data.
  status("GO-GSEA | reusing existing complete results; no rerun of six enrichments")
  gsea <- utils::read.csv(full_gsea_path, stringsAsFactors = FALSE)
  stopifnot(all(c("ID", "Description", "NES", "p.adjust", "core_enrichment",
                  "contrast_id") %in% names(gsea)),
            setequal(unique(gsea$contrast_id), plan$id))
  first <- utils::read.csv(raw_paths[1], stringsAsFactors = FALSE)
  annotation_map <- first[, c("ENSEMBL", "ENTREZID", "SYMBOL"), drop = FALSE]
} else {
  all_gsea <- list()
  annotation_map <- NULL
  for (i in seq_len(nrow(plan))) {
    id <- plan$id[i]
    status(sprintf("GO-GSEA %d/6 | %s | reading Wald statistics", i, id))
    x <- utils::read.csv(raw_paths[i], stringsAsFactors = FALSE)
    needed <- c("ENSEMBL", "ENTREZID", "SYMBOL", "baseMean", "stat", "padj", "lfc_method")
    if (!all(needed %in% names(x)) || anyNA(x$ENSEMBL) || anyDuplicated(x$ENSEMBL) ||
        anyNA(x$lfc_method) || any(x$lfc_method != "Wald_unshrunken"))
      stop("Unexpected Wald result format: ", raw_paths[i])
    if (is.null(annotation_map)) annotation_map <- x[, c("ENSEMBL", "ENTREZID", "SYMBOL"), drop = FALSE]
    keys <- as.character(x$ENTREZID)
    keep <- valid_entrez(keys) & is.finite(x$stat)
    ranked <- data.frame(ENTREZID = keys[keep], stat = x$stat[keep],
                         baseMean = x$baseMean[keep], ENSEMBL = x$ENSEMBL[keep],
                         stringsAsFactors = FALSE)
    ranked <- ranked[order(ranked$ENTREZID, -ranked$baseMean, ranked$ENSEMBL), , drop = FALSE]
    ranked <- ranked[!duplicated(ranked$ENTREZID), , drop = FALSE]
    gene_list <- ranked$stat
    names(gene_list) <- ranked$ENTREZID
    ## Tie-breaking by ENTREZ ID is deterministic, without perturbing Wald scores.
    ## Equal statistics are legitimate; do not add random jitter to the ranks.
    gene_list <- gene_list[order(-gene_list, names(gene_list), method = "radix")]
    n_tied <- length(gene_list) - length(unique(unname(gene_list)))
    if (n_tied) status(sprintf("GO-GSEA %d/6 | %s | %d tied Wald ranks (%.2f%%); deterministic ENTREZ ordering",
                               i, id, n_tied, 100 * n_tied / length(gene_list)))
    if (length(gene_list) < 1000L) stop("Too few uniquely mapped Wald genes: ", id)
    status(sprintf("GO-GSEA %d/6 | %s | %d ranked ENTREZ genes | gseGO started",
                   i, id, length(gene_list)))
    t0 <- proc.time()[["elapsed"]]
    ## pvalueCutoff=1 retains nonsignificant terms for an honest complete table.
    ## fgsea multilevel permutation/adaptive precision is handled by the package.
    ## fgsea uses resampling to estimate p-values: fix the RNG seed per contrast.
    set.seed(20261012L + i)
    fit <- tryCatch(clusterProfiler::gseGO(
      geneList = gene_list, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID",
      ont = "BP", exponent = 1, minGSSize = 15, maxGSSize = 500,
      eps = gsea_eps, pvalueCutoff = 1, pAdjustMethod = "BH", by = "fgsea", verbose = FALSE),
      error = function(e) stop("GO-GSEA failed for ", id, ": ", conditionMessage(e)))
    z <- as.data.frame(fit)
    if (nrow(z)) {
      ## gseGO / gseaResult uses 'setSize', whereas enrichGO / enrichResult
      ## uses 'Count'. Normalize this metadata field; never change statistics.
      essential <- c("ID", "Description", "NES", "p.adjust", "core_enrichment")
      if (!all(essential %in% names(z)))
        stop("Unexpected GSEA output columns for ", id, ": ", paste(names(z), collapse = ", "),
             ". Required: ", paste(essential, collapse = ", "))
      if (!"setSize" %in% names(z)) {
        if ("Count" %in% names(z)) z$setSize <- z$Count
        else if ("size" %in% names(z)) z$setSize <- z$size
        else stop("Missing gene-set size in gseGO result for ", id,
                  "; returned columns: ", paste(names(z), collapse = ", "))
      }
      if (!"Count" %in% names(z)) z$Count <- z$setSize
      z$contrast_id <- id
      z$cell_line <- plan$cell_line[i]
      z$radiation <- plan$radiation[i]
      all_gsea[[id]] <- z
    }
    status(sprintf("GO-GSEA %d/6 | %s | %.1fs | %d tested terms, %d FDR<0.05",
      i, id, proc.time()[["elapsed"]] - t0, nrow(z),
      if (nrow(z)) sum(is.finite(z$p.adjust) & z$p.adjust < alpha) else 0L))
  }
  
  ## Keep column structure when none of the contrasts returns a term.
  gsea <- if (length(all_gsea)) do.call(rbind, all_gsea) else data.frame(
    ID = character(), Description = character(), NES = numeric(), p.adjust = numeric(),
    core_enrichment = character(), Count = integer(), contrast_id = character(),
    cell_line = character(), radiation = character())

}
if (!using_cached_gsea)
  utils::write.csv(gsea, full_gsea_path, row.names = FALSE)
gsea_sig <- gsea[is.finite(gsea$p.adjust) & gsea$p.adjust < alpha, , drop = FALSE]

## ************************************************************************************************
## Figure A - compact GO-ORA results organized BY MODULE, not a sparse matrix
## Horizontal axis shows -log10 FDR (capped for display); dots show fraction of
## annotated module genes. GO biology is functional membership, not activation.
status("3/5 | Plotting semantically nonredundant module GO-ORA")
fig_theme <- plot_theme(10.5) + ggplot2::theme(
  panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
  axis.line = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
  strip.background = ggplot2::element_blank(),
  plot.title = ggplot2::element_text(face = "bold", size = 16),
  plot.subtitle = ggplot2::element_text(color = "grey35", size = 10),
  plot.caption = ggplot2::element_text(color = "grey40", size = 8.5))

if (nrow(ora_shown)) {
  ## Show only modules with significant terms. M03 has no significant
  ## GO-BP term and is documented in the caption, NOT as an empty facet.
  ## Empty facets combined with free_y/free_y space can break grid viewports.
  op <- ora_shown
  op$module <- droplevels(factor(op$module, levels = pair_order))
  op$term_key <- paste(op$module, op$ID)
  op$term_key <- factor(op$term_key, levels = rev(unique(op$term_key)))
  term_labels <- setNames(op$Description, as.character(op$term_key))
  p_ora <- ggplot2::ggplot(op, ggplot2::aes(x = pmin(minus_log10_FDR, 30), y = term_key)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = pmin(minus_log10_FDR, 30),
                                      yend = term_key), color = "grey85", linewidth = 0.45) +
    ggplot2::geom_point(ggplot2::aes(size = gene_ratio_percent), color = "#3C5488") +
    ggplot2::facet_grid(module ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
    ggplot2::scale_y_discrete(labels = function(z) short_name(term_labels[z], 61L)) +
    ggplot2::scale_x_continuous(limits = c(0, 31), breaks = c(0, 10, 20, 30),
                                expand = ggplot2::expansion(mult = 0)) +
    ggplot2::scale_size_continuous(range = c(2.6, 6.3),
                                   name = "Genes in GO term\n(% of mapped module)") +
    ggplot2::labs(title = "Biological processes within radiation-response modules",
      subtitle = "Up to 3 semantically distinct GO-BP terms per module (Wang similarity)",
      x = expression(paste(-log[10], " adjusted p-value (display cap: 30)")), y = NULL,
      caption = "Response direction is shown separately in script 03/05. M03: no GO-BP term at FDR < 0.05.") +
    fig_theme + ggplot2::theme(
      strip.text.y = ggplot2::element_text(angle = 0, face = "bold", size = 10),
      panel.spacing.y = grid::unit(0.22, "lines"),
      axis.text.y = ggplot2::element_text(size = 9.3, color = "black"),
      legend.position = "right")
  export("GO_BP_module_ORA_dotplot", p_ora, 14.5, 10.0)
} else status("No significant module ORA terms; dotplot not generated")

## ************************************************************************************************
## Figure B - GO-GSEA heatmap: semantic pruning across the SIX interaction tests
## Identical GO IDs are used across the six columns. NES sign refers to the
## SIGNED KO minus WT contrast, not stronger/weaker absolute response.
status("4/5 | Semantic GSEA term reduction and publication-style heatmap")
selected_gsea <- character()
if (nrow(gsea_sig)) {
  gs <- gsea_sig[order(gsea_sig$p.adjust, gsea_sig$ID), , drop = FALSE]
  counts <- table(gs$ID)
  gs$n_contrasts_significant <- as.integer(counts[gs$ID])
  gs <- gs[order(-gs$n_contrasts_significant, gs$p.adjust, gs$ID), , drop = FALSE]
  gs <- gs[!duplicated(gs$ID), , drop = FALSE]
  ## Reviewed functionally informative anchors only if already significant.
  anchors <- c("GO:0006260", "GO:0007059", "GO:0006119", "GO:0022900",
               "GO:0006954", "GO:0007264", "GO:0060271")
  chosen <- select_semantic(gs, n = 12L,
                            anchors = intersect(anchors, gs$ID), max_candidates = 120L)
  selected_gsea <- chosen$ID
  chosen$display_selection <- ifelse(chosen$ID %in% anchors,
                                     "reviewed_anchor", "Wang_nonredundant_multicontrast")
  chosen$theme <- classify_theme(chosen$Description)
  utils::write.csv(chosen, file.path(tab_dir, "GO_BP_GSEA_display_terms.csv"), row.names = FALSE)
  status("GSEA | ", nrow(gsea_sig), " significant rows across six contrasts; ",
         length(selected_gsea), " semantically distinct GO IDs displayed")

  term_data <- chosen[, c("ID", "Description", "theme")]
  term_data$theme[is.na(term_data$theme)] <- "Other GO processes"
  theme_order <- c(theme_levels, "Other GO processes")
  term_data <- term_data[order(match(term_data$theme, theme_order),
                               match(term_data$ID, selected_gsea)), , drop = FALSE]
  id_order <- term_data$ID
  show <- expand.grid(ID = id_order, contrast_id = plan$id, stringsAsFactors = FALSE)
  idx <- match(paste(show$ID, show$contrast_id), paste(gsea$ID, gsea$contrast_id))
  show$NES <- gsea$NES[idx]
  show$FDR <- gsea$p.adjust[idx]
  show$significant <- is.finite(show$FDR) & show$FDR < alpha
  show$label <- ifelse(show$significant, "*", "")
  show$cell_line <- sub("_.*$", "", show$contrast_id)
  show$radiation <- sub("^.*_", "", show$contrast_id)
  show$radiation <- factor(show$radiation, levels = c("7gy", "UVA", "UVB"),
                            labels = c("7 Gy", "UVA", "UVB"))
  show$ID <- factor(show$ID, levels = rev(id_order))
  descriptions <- setNames(term_data$Description, term_data$ID)

  p_gsea <- ggplot2::ggplot(show, ggplot2::aes(x = radiation, y = ID, fill = NES)) +
    ggplot2::geom_tile(ggplot2::aes(alpha = ifelse(significant, "FDR < 0.05", "other")),
                       color = "white", linewidth = 0.9) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 4.2) +
    ggplot2::facet_grid(. ~ cell_line, scales = "free_x", space = "free_x") +
    ggplot2::scale_y_discrete(labels = function(z) short_name(descriptions[z], 57L)) +
    ggplot2::scale_alpha_manual(values = c("FDR < 0.05" = 1, "other" = 0.28),
                                guide = "none") +
    ggplot2::scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35",
                                  midpoint = 0, na.value = "grey92", name = "NES",
                                  limits = c(-2.5, 2.5), oob = scales::squish) +
    ggplot2::labs(title = "Coordinated p21-dependent transcriptional programs",
      subtitle = "Six formal interactions | semantically distinct GO-BP terms | * FDR < 0.05",
      x = NULL, y = NULL,
      caption = "Red: GO genes rank toward more positive KO-minus-WT interactions; blue: more negative. NES is NOT activation or response strength.") +
    fig_theme + ggplot2::theme(
      strip.text.x = ggplot2::element_text(face = "bold", size = 11),
      axis.text.y = ggplot2::element_text(size = 9.5, color = "black"),
      legend.position = "right")
  export("GO_BP_interaction_GSEA_heatmap", p_gsea, 15, max(7.0, length(id_order) * 0.45 + 2.5))
} else {
  utils::write.csv(gsea_sig, file.path(tab_dir, "GO_BP_GSEA_display_terms.csv"), row.names = FALSE)
  status("No significant GSEA GO terms; heatmap omitted")
}

## ************************************************************************************************
## Figure C - compact thematic synthesis: SAME biological themes, with a
## transparent independent representative GO ID for ORA and GSEA.
## These two IDs may DIFFER. Do not interpret them as an exact GO-term match;
## there is no pooled p-value and no new statistical gene-set test.
status("5/5 | Building six-theme functional landscape and provenance table")
ora_theme <- ora_sig
ora_theme$theme <- classify_theme(ora_theme$Description)
gsea_theme <- gsea
gsea_theme$theme <- classify_theme(gsea_theme$Description)

ora_anchor <- list(
  "Cell cycle & DNA replication" = c("GO:0006260"),
  "Mitochondrial energy" = c("GO:0140053", "GO:0032543"),
  "Inflammation & cell communication" = c("GO:0006954"),
  "Autophagy & lysosome" = c("GO:0006914"),
  "Cytoskeleton & adhesion" = c("GO:0007264"),
  "Cilium organization" = c("GO:0060271"))
gsea_anchor <- list(
  "Cell cycle & DNA replication" = c("GO:0006260"),
  "Mitochondrial energy" = c("GO:0006119", "GO:0022900"),
  "Inflammation & cell communication" = c("GO:0006954"),
  "Autophagy & lysosome" = c("GO:0006914"),
  "Cytoskeleton & adhesion" = c("GO:0007264"),
  "Cilium organization" = c("GO:0060271"))

chosen_theme <- data.frame(theme = theme_levels, ORA_ID = NA_character_,
                           GSEA_ID = NA_character_, stringsAsFactors = FALSE)
for (i in seq_along(theme_levels)) {
  theme <- theme_levels[i]
  chosen_theme$ORA_ID[i] <- representative_id(ora_theme, theme, ora_anchor[[theme]])
  chosen_theme$GSEA_ID[i] <- representative_id(gsea_theme, theme, gsea_anchor[[theme]])
}
lookup <- rbind(ora[, c("ID", "Description")], gsea[, c("ID", "Description")])
lookup <- lookup[!duplicated(lookup$ID), , drop = FALSE]
name_lookup <- setNames(lookup$Description, lookup$ID)
chosen_theme$ORA_term <- unname(name_lookup[chosen_theme$ORA_ID])
chosen_theme$GSEA_term <- unname(name_lookup[chosen_theme$GSEA_ID])
utils::write.csv(chosen_theme, file.path(tab_dir, "functional_landscape_theme_definitions.csv"),
                 row.names = FALSE)

ora_grid <- expand.grid(theme = theme_levels, module = pair_order, stringsAsFactors = FALSE)
ora_grid$GO_ID <- chosen_theme$ORA_ID[match(ora_grid$theme, chosen_theme$theme)]
ora_grid$Description <- chosen_theme$ORA_term[match(ora_grid$theme, chosen_theme$theme)]
idx <- match(paste(ora_grid$GO_ID, ora_grid$module), paste(ora$ID, ora$module))
ora_grid$FDR <- ora$p.adjust[idx]
ora_grid$significant <- is.finite(ora_grid$FDR) & ora_grid$FDR < alpha
ora_grid$score <- ifelse(ora_grid$significant,
                         -log10(pmax(ora_grid$FDR, .Machine$double.xmin)), NA_real_)

gsea_grid <- expand.grid(theme = theme_levels, contrast_id = plan$id, stringsAsFactors = FALSE)
gsea_grid$GO_ID <- chosen_theme$GSEA_ID[match(gsea_grid$theme, chosen_theme$theme)]
gsea_grid$Description <- chosen_theme$GSEA_term[match(gsea_grid$theme, chosen_theme$theme)]
idx <- match(paste(gsea_grid$GO_ID, gsea_grid$contrast_id),
             paste(gsea$ID, gsea$contrast_id))
gsea_grid$FDR <- gsea$p.adjust[idx]
gsea_grid$score <- gsea$NES[idx]
gsea_grid$significant <- is.finite(gsea_grid$FDR) & gsea_grid$FDR < alpha

utils::write.csv(rbind(
  data.frame(analysis = "module_ORA", GO_ID = ora_grid$GO_ID,
             Description = ora_grid$Description, theme = ora_grid$theme,
             condition = ora_grid$module, FDR = ora_grid$FDR,
             statistic = ora_grid$score, significant = ora_grid$significant),
  data.frame(analysis = "interaction_GSEA", GO_ID = gsea_grid$GO_ID,
             Description = gsea_grid$Description, theme = gsea_grid$theme,
             condition = gsea_grid$contrast_id, FDR = gsea_grid$FDR,
             statistic = gsea_grid$score, significant = gsea_grid$significant)),
  file.path(tab_dir, "functional_landscape_terms.csv"), row.names = FALSE)

ora_grid$theme <- factor(ora_grid$theme, levels = rev(theme_levels))
gsea_grid$theme <- factor(gsea_grid$theme, levels = rev(theme_levels))
ora_grid$module <- factor(ora_grid$module, levels = pair_order)
gsea_grid$contrast_id <- factor(gsea_grid$contrast_id, levels = plan$id)

p_left <- ggplot2::ggplot(subset(ora_grid, significant),
                          ggplot2::aes(x = module, y = theme, color = pmin(score, 30))) +
  ggplot2::geom_point(size = 4.2) +
  ggplot2::scale_color_gradient(low = "#9ECAE1", high = "#3C5488",
                                limits = c(0, 30), oob = scales::squish,
                                name = expression(-log[10]~"ORA FDR")) +
  ggplot2::scale_y_discrete(drop = FALSE) +
  ggplot2::labs(title = "General radiation response", subtitle = "GO-ORA by module",
                x = NULL, y = NULL) + fig_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9),
                 legend.position = "bottom")
gsea_grid$label <- ifelse(gsea_grid$significant, "*", "")
p_right <- ggplot2::ggplot(gsea_grid,
                           ggplot2::aes(x = contrast_id, y = theme, fill = score)) +
  ggplot2::geom_tile(ggplot2::aes(alpha = ifelse(significant, "FDR < 0.05", "other")),
                     color = "white", linewidth = 1.0) +
  ggplot2::geom_text(ggplot2::aes(label = label), size = 4.1) +
  ggplot2::scale_alpha_manual(values = c("FDR < 0.05" = 1, "other" = 0.30),
                              guide = "none") +
  ggplot2::scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35",
                                midpoint = 0, limits = c(-2.5, 2.5),
                                oob = scales::squish, na.value = "grey93", name = "GSEA NES") +
  ggplot2::scale_y_discrete(drop = FALSE) +
  ggplot2::scale_x_discrete(labels = function(z)
    gsub(" 7gy$", " 7 Gy", gsub("_", " ", z))) +
  ggplot2::labs(title = "p21 modulation", subtitle = "GO-GSEA by formal interaction",
                x = NULL, y = NULL) + fig_theme +
  ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(angle = 25, hjust = 1, size = 8.5),
                 legend.position = "bottom")

landscape <- p_left + p_right + patchwork::plot_layout(widths = c(1.05, 1.30)) +
  patchwork::plot_annotation(
    title = "Functional themes of irradiation response and p21 modulation",
    subtitle = "Six descriptive themes | selected representative GO terms per method (IDs can differ)",
    caption = "Stars: GSEA FDR < 0.05. ORA and GSEA are different analyses; exact GO IDs, tests and FDRs are listed in the accompanying CSV.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      plot.caption = ggplot2::element_text(size = 9)))
export("GO_BP_functional_landscape", landscape, 18, 7.3)

## ******
## Preserve a compact exact-GO overlap table for selected representative IDs.
## This is optional provenance, not an additional biological analysis.
annotation_map$ENTREZID <- as.character(annotation_map$ENTREZID)
annotation_map <- annotation_map[valid_entrez(annotation_map$ENTREZID), , drop = FALSE]
annotation_map <- annotation_map[!duplicated(annotation_map$ENTREZID), , drop = FALSE]
symbol_map <- setNames(annotation_map$SYMBOL, annotation_map$ENTREZID)
selected_ids <- unique(na.omit(c(chosen_theme$ORA_ID, chosen_theme$GSEA_ID)))
overlap_rows <- list()
for (id in selected_ids) {
  a <- ora_sig[ora_sig$ID == id, , drop = FALSE]
  b <- gsea_sig[gsea_sig$ID == id, , drop = FALSE]
  if (!nrow(a) || !nrow(b)) next
  for (i in seq_len(nrow(a))) for (j in seq_len(nrow(b))) {
    shared <- intersect(split_ids(a$geneID[i]), split_ids(b$core_enrichment[j]))
    genes <- unname(symbol_map[shared])
    genes[is.na(genes) | !nzchar(genes)] <- shared[is.na(genes) | !nzchar(genes)]
    overlap_rows[[length(overlap_rows) + 1L]] <- data.frame(
      GO_ID = id, Description = a$Description[i], module = a$module[i],
      contrast_id = b$contrast_id[j], ORA_FDR = a$p.adjust[i],
      GSEA_FDR = b$p.adjust[j], GSEA_NES = b$NES[j],
      n_shared = length(shared), shared_ENTREZID = paste(shared, collapse = ";"),
      shared_SYMBOL = paste(genes, collapse = ";"))
  }
}
empty <- data.frame(GO_ID = character(), Description = character(), module = character(),
  contrast_id = character(), ORA_FDR = numeric(), GSEA_FDR = numeric(),
  GSEA_NES = numeric(), n_shared = integer(),
  shared_ENTREZID = character(), shared_SYMBOL = character())
overlap_tab <- if (length(overlap_rows)) do.call(rbind, overlap_rows) else empty
utils::write.csv(overlap_tab, file.path(tab_dir, "GO_module_leading_edge_overlap.csv"),
                 row.names = FALSE)
status("Complete | full ORA/GSEA preserved; representative figures and provenance: ", root)

