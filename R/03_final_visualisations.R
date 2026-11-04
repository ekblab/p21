## p21 / 03 | Final module visualisations and GO biological-category audit
## Standalone: source("R/03_final_visualisations.R")
## Inputs: exported module/response CSVs and the fully semantically reduced
## GO-BP table generated in stage 1. No DESeq2, no k-means, no GO enrichment.
## All plots: ekbSeq::export_plot_dual (PNG + SVG).
## Audit: mandatory GO_BP_theme_category_audit.xlsx + CSVs with every
## GO ID, term, rule, override, multiple matches, module & plotted status.
## Categories are curated EXPLORATORY biological themes, not GO-slim classes.

if (!requireNamespace("here", quietly = TRUE)) stop("Missing package: here")
source(here::here("R", "03_shared.R"), local = TRUE)
p21_check(c("ggplot2", "scales", "patchwork", "ekbSeq", "openxlsx"))

## ************************************************************************************************
## Report figures from already exported module tables (no re-analysis)
## Run this script directly to refresh final figures and the category audit.
## It is also sourced automatically after the complete stage-1 + stage-2 workflow.
## Assignment comes from the EXISTING k-means CSV; figures use ORIGINAL
## ashr log2FC. Styling does not alter any modeled estimates.
## ************************************************************************************************


figure_dir <- p21_module_dir
response_file <- file.path(p21_result_root, "response_matrix_DEG_union_ashr.csv")
module_order <- p21_module_order
module_colors <- p21_module_colors
condition_ids <- c(
  "HT1080_wt_7gy", "HT1080_p21_neg_7gy", "HT1080_wt_UVA", "HT1080_p21_neg_UVA",
  "HT1080_wt_UVB", "HT1080_p21_neg_UVB", "MCF10A_wt_7gy", "MCF10A_p21_neg_7gy",
  "MCF10A_wt_UVA", "MCF10A_p21_neg_UVA", "MCF10A_wt_UVB", "MCF10A_p21_neg_UVB")

read_module_csv <- function(path, required) {
  if (!file.exists(path)) stop("Missing module input: ", path)
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Unexpected columns in ", path, ": ", paste(missing, collapse = ", "))
  x
}

## Same two-tier header for both output figures (matches the schematic).
## Top: cell line; middle: irradiation; bottom: WT / KO.
add_condition_header <- function(top) {
  cells <- data.frame(x = c(3.5, 9.5), x0 = c(.55, 6.55), x1 = c(6.45, 12.45),
                      label = c("HT1080", "MCF10A"))
  treatments <- data.frame(x = c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5),
                           x0 = seq(.56, 10.56, 2), x1 = seq(2.44, 12.44, 2),
                           label = rep(c("7 Gy", "UVA", "UVB"), 2))
  genotypes <- data.frame(x = 1:12, label = rep(c("WT", "KO"), 6))
  list(
    ggplot2::geom_segment(data = cells,
      ggplot2::aes(x = x0, xend = x1, y = top + 2.17, yend = top + 2.17),
      inherit.aes = FALSE, linewidth = .50, color = "#222222"),
    ggplot2::geom_text(data = cells,
      ggplot2::aes(x = x, y = top + 2.48, label = label),
      inherit.aes = FALSE, fontface = "bold", size = 4.15),
    ggplot2::geom_segment(data = treatments,
      ggplot2::aes(x = x0, xend = x1, y = top + 1.32, yend = top + 1.32),
      inherit.aes = FALSE, linewidth = .33, color = "#444444"),
    ggplot2::geom_text(data = treatments,
      ggplot2::aes(x = x, y = top + 1.68, label = label),
      inherit.aes = FALSE, size = 3.25),
    ggplot2::geom_text(data = genotypes,
      ggplot2::aes(x = x, y = top + .60, label = label),
      inherit.aes = FALSE, size = 2.82, color = "#222222"),
    ggplot2::geom_vline(xintercept = 6.5, color = "#333333", linewidth = .75),
    ggplot2::geom_vline(xintercept = c(2.5, 4.5, 8.5, 10.5),
                        color = "grey90", linewidth = .25)
  )
}

module_fill <- function(cap, title) {
  ## Higher contrast at moderate values than an evenly spaced large-range gradient.
  breaks <- c(-cap, -.75*cap, -.40*cap, -.12*cap, 0,
              .12*cap, .40*cap, .75*cap, cap)
  ggplot2::scale_fill_gradientn(
    colors = c("#12386F", "#306FB3", "#94BDE2", "#E0ECF8", "#FFFFFF",
               "#FAE1DD", "#E89A94", "#C94644", "#941E2C"),
    values = scales::rescale(breaks, from = c(-cap, cap)),
    limits = c(-cap, cap), oob = scales::squish,
    name = title, breaks = c(-cap, 0, cap))
}

module_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) + ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    legend.position = "bottom", legend.box = "horizontal",
    plot.title = ggplot2::element_text(face = "bold", size = 15),
    plot.subtitle = ggplot2::element_text(size = 10, color = "grey30"),
    plot.caption = ggplot2::element_text(size = 8, color = "grey35"),
    plot.margin = ggplot2::margin(10, 16, 10, 10))
}

export_module_plot <- function(plot, filename, width, height) {
  p21_export(file.path(figure_dir, filename), plot, width = width, height = height)
}

plot_module_means <- function(order, name) {
  df <- read_module_csv(file.path(figure_dir, "module_mean_ashr_LFC.csv"), "module")
  if (anyDuplicated(df$module) || !setequal(df$module, module_order))
    stop("Module means must contain M01-M08 exactly once")
  if (all(condition_ids %in% names(df))) {
    mat <- as.matrix(df[match(order, df$module), condition_ids, drop = FALSE])
  } else {
    ## Older script-03 exports can lack explicit condition column names.
    ## Use the separately exported *named* module/genotype/irradiation table
    ## rather than assuming any ambiguous column order.
    dd <- read_module_csv(file.path(figure_dir, "module_genotype_deltas_descriptive.csv"),
      c("module", "cell_line", "irradiation", "WT", "p21_negative"))
    keys <- paste(dd$module, dd$cell_line, dd$irradiation)
    if (anyDuplicated(keys)) stop("Duplicate module/cell line/irradiation rows")
    condition_line <- rep(c("HT1080", "MCF10A"), each = 6)
    condition_rad <- rep(rep(c("7gy", "UVA", "UVB"), each = 2), 2)
    condition_is_KO <- rep(c(FALSE, TRUE), 6)
    mat <- matrix(NA_real_, nrow = 8, ncol = 12)
    for (i in seq_along(order)) {
      for (j in seq_along(condition_ids)) {
        row <- match(paste(order[i], condition_line[j], condition_rad[j]), keys)
        if (is.na(row)) stop("Incomplete module genotype-delta export")
        mat[i, j] <- if (condition_is_KO[j]) dd$p21_negative[row] else dd$WT[row]
      }
    }
    message("Original means CSV lacks named columns; using annotated genotype-delta table.")
  }
  storage.mode(mat) <- "numeric"
  if (any(!is.finite(mat))) stop("Nonfinite module means")
  row_y <- if (identical(order, module_order)) {
    c(13.25, 12.00, 9.80, 8.55, 6.35, 5.10, 2.90, 1.65)
  } else {
    rev(seq(1.65, by = 1.52, length.out = 8))
  }
  heat <- data.frame(x = rep(1:12, times = 8), y = rep(row_y, each = 12),
                     log2FC = as.vector(t(mat)))
  labs <- data.frame(module = order, y = row_y,
                     text = paste(order, "·", c("HT1080 UVA induction", "HT1080 UVA repression",
                                                "MCF10A 7 Gy/UVA induction", "MCF10A 7 Gy/UVA repression",
                                                "HT1080 UVB induction", "HT1080 UVB repression",
                                                "MCF10A UVB induction", "MCF10A UVB repression")[match(order, module_order)]))
  ## Strictly label profiles rather than claiming pathway or protein activity.
  if (!identical(order, module_order)) {
    labels <- read_module_csv(file.path(figure_dir, "module_descriptive_labels.csv"),
                              c("module", "descriptive_pattern"))
    labs$text <- paste(labs$module, "·", labels$descriptive_pattern[match(labs$module, labels$module)])
  }
  labs$color <- unname(module_colors[labs$module])
  top <- max(row_y) + .47
  cap <- max(1.0, ceiling(max(abs(mat)) * 10) / 10)
  p <- ggplot2::ggplot(heat, ggplot2::aes(x = x, y = y, fill = log2FC)) +
    ggplot2::geom_tile(width = .96, height = .98, color = "white", linewidth = .35) +
    {if (identical(order, module_order))
      ggplot2::geom_segment(data = data.frame(y = c(10.90, 7.45, 4.00)),
        ggplot2::aes(x = -4.1, xend = 12.45, y = y, yend = y), inherit.aes = FALSE,
        color = "grey82", linewidth = .30)} +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.2f", log2FC)),
                       size = 2.55, color = "#222222") +
    ggplot2::geom_rect(data = labs,
      ggplot2::aes(xmin = -4.05, xmax = -3.89, ymin = y - .48, ymax = y + .48),
      inherit.aes = FALSE, fill = labs$color, color = NA) +
    ggplot2::geom_text(data = labs,
      ggplot2::aes(x = -3.75, y = y, label = text),
      inherit.aes = FALSE, hjust = 0, size = 2.75, color = "#222222") +
    add_condition_header(top) + module_fill(cap, "Mean ashr log2FC") +
    ggplot2::scale_x_continuous(limits = c(-4.20, 12.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(.25, top + 2.75), expand = c(0, 0)) +
    module_theme(10) +
    ggplot2::labs(title = "Mean irradiation responses across L2 modules",
      subtitle = "Original ashr log2FC | 8 exploratory response-pattern modules",
      caption = "Blue: negative; red: positive response vs control. Means are descriptive, not interaction tests.")
  export_module_plot(p, name, 17.2, 7.5)
}

select_module_genes <- function() {
  genes <- read_module_csv(file.path(figure_dir, "representative_genes_with_selection_reason.csv"),
    c("ENSEMBL", "SYMBOL", "module", "profile_cosine", "response_L2_magnitude", "selection_reason"))
  if (anyNA(genes$ENSEMBL) || anyDuplicated(genes$ENSEMBL)) stop("Missing/duplicate ENSEMBL in shortlist")
  if (!identical(as.integer(table(factor(genes$module, levels = sprintf("M%02d", 1:8)))), rep(10L, 8)))
    stop("Expected exactly ten preselected genes per module")
  ## Fixed shortlist: selected for module-pattern similarity AND biological relevance.
  choices <- list(
    M01 = c("CDKN1A", "MDM2", "BTG2"),
    M08 = c("MCM7", "TFDP1", "MCM2"),
    M02 = c("IL6", "INHBA", "IL1B"),
    M06 = c("PIK3CA", "PRKAG2", "PDK4"),
    M03 = c("MUC13", "LCE1C", "CEACAM1"),
    M05 = c("EPHB1", "NTM", "PTPRD"),
    M04 = c("ANKRD1", "CCL26", "TNFAIP3"),
    M07 = c("PPM1L", "TENM2", "WDPCP"))
  keys <- do.call(rbind, lapply(names(choices), function(mod)
    data.frame(module = mod, SYMBOL = choices[[mod]], rank = 1:3)))
  if (anyDuplicated(paste(genes$module, genes$SYMBOL)))
    stop("Nonunique module/SYMBOL pairs in the original shortlist")
  idx <- match(paste(keys$module, keys$SYMBOL), paste(genes$module, genes$SYMBOL))
  if (anyNA(idx)) stop("Missing candidate in original module shortlist: ",
                       paste(keys$SYMBOL[is.na(idx)], collapse = ", "))
  selected <- genes[idx, , drop = FALSE]
  selected$rank <- keys$rank
  if (anyDuplicated(selected$ENSEMBL) || any(!is.finite(selected$profile_cosine)) ||
      any(selected$profile_cosine < .75)) stop("Selected genes fail module/profile audit")
  utils::write.csv(selected, file.path(figure_dir, "representative_genes_selected_3per_module.csv"), row.names = FALSE)
  utils::write.csv2(selected, file.path(figure_dir, "representative_genes_selected_3per_module_excel_DE.csv"), row.names = FALSE)
  utils::write.csv2(genes, file.path(figure_dir, "representative_genes_with_selection_reason_excel_DE.csv"), row.names = FALSE)
  selected
}

plot_module_genes <- function() {
  genes <- select_module_genes()
  resp <- read_module_csv(response_file, c("ENSEMBL", condition_ids))
  if (anyNA(resp$ENSEMBL) || anyDuplicated(resp$ENSEMBL)) stop("Invalid response-matrix ENSEMBL IDs")
  idx <- match(genes$ENSEMBL, resp$ENSEMBL)
  if (anyNA(idx)) stop("Selected ENSEMBL IDs absent from response matrix")
  mat <- as.matrix(resp[idx, condition_ids, drop = FALSE]); storage.mode(mat) <- "numeric"
  if (any(!is.finite(mat))) stop("Nonfinite selected gene effects")

  block <- match(genes$module, module_order)
  module_base <- c(29.4, 26.0, 21.4, 18.0, 13.4, 10.0, 5.4, 2.0)
  row_y <- module_base[block] + 3 - genes$rank
  df <- data.frame(x = rep(1:12, times = nrow(genes)), y = rep(row_y, each = 12),
                   log2FC = as.vector(t(mat)))
  row_labels <- data.frame(y = row_y, symbol = genes$SYMBOL)
  mlabel <- data.frame(module = module_order,
                       y = module_base + 1,
                       ymin = module_base - .42,
                       ymax = module_base + 2.42)
  mlabel$color <- unname(module_colors[mlabel$module])
  top <- max(row_y) + .49
  cap <- max(1.5, as.numeric(stats::quantile(abs(mat), .98, na.rm = TRUE)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = log2FC)) +
    ggplot2::geom_tile(width = .96, height = .91, color = "white", linewidth = .18) +
    ggplot2::geom_segment(data = data.frame(y = c(24.7, 16.7, 8.7)),
      ggplot2::aes(x = -2.75, xend = 12.45, y = y, yend = y), inherit.aes = FALSE,
      color = "grey82", linewidth = .36) +
    ggplot2::geom_rect(data = mlabel,
      ggplot2::aes(xmin = -1.72, xmax = -1.52, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = mlabel$color, color = NA) +
    ggplot2::geom_text(data = mlabel,
      ggplot2::aes(x = -2.05, y = y, label = module),
      inherit.aes = FALSE, hjust = 1, fontface = "bold", size = 3.65) +
    ggplot2::geom_text(data = row_labels,
      ggplot2::aes(x = .32, y = y, label = symbol),
      inherit.aes = FALSE, hjust = 1, size = 3.10, color = "#222222") +
    add_condition_header(top) + module_fill(cap, "Gene ashr log2FC") +
    ggplot2::scale_x_continuous(limits = c(-3.05, 12.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(.35, top + 2.75), expand = c(0, 0)) +
    module_theme(10) +
    ggplot2::labs(title = "Representative genes across radiation-response modules",
      subtitle = "Three biologically informative, pattern-consistent genes per module",
      caption = "Original ashr log2FC, not clustering-normalized values. Module colors denote membership.")
  export_module_plot(p, "representative_gene_heatmap_3per_module", 14.5, 11.5)
}

message("Updating report module figures from existing CSV results...")
plot_module_means(module_order, "module_mean_response_heatmap_ordered")
plot_module_means(sprintf("M%02d", 1:8), "module_mean_response_heatmap_header")
plot_module_genes()
message("Module report figures saved to: ", figure_dir)

## ************************************************************************************************
## Exploratory GO biological-theme characterization from the COMPLETE,
## semantically reduced, FDR-significant module ORA sets.
##
## Summary proportions are percentages of representative GO TERM counts,
## NOT percentages of genes, pathway activity, or ontology-independent events.
## A transparent, mutually exclusive keyword dictionary assigns each term to
## one high-level theme for visualization; unclassified terms remain "Other".
## The manually reviewed categories are exploratory, not GO-slim annotations.
## Table outputs allow auditing every term/category mapping.
## Categories are manual interpretive bins, NOT formal GO hierarchy partitions.
##
## Main figure A: theme-aware GO-BP dotplot (representatives from leading themes).
## Main figure B: theme-composition donuts with term counts (M03 = none).
## Original ORA estimates, significance and modules are NEVER modified.
## ************************************************************************************************

go_theme_patterns <- c(
  "Immune / inflammatory" = "immune|immun|leukocyt|lymphocyt|\\bT[ -]cell\\b|\\bB[ -]cell\\b|myeloid|mononuclear|macrophag|phagocyt|antigen|interferon|interleukin|cytokine(?!sis)|chemokine|inflamm|defense response|antiviral|viral (?:process|transcription)|virus|bacteri|lipopolysaccharide|pathogen|NF-kappaB|hematopoietic|hematopoiesis|host.mediate|symbiont|wound heal",
  "Cell cycle / genome" = "cell cycle|cell division|cytokinesis|chromosom|chromatid|chromatin|heterochromatin|kinetochor|centromer|centrosom|nucleosom|histone|epigenet|spindle|replication|\\bDNA\\b|DNA[- ]|recombin|mitotic|meiotic|nuclear division|metaphase|anaphase|G0.to.G1|G1.to.S|cyclin.dependent|anaphase.promoting|protein.DNA complex",
  "Mitochondria / energy" = "mitochondri|oxidative phosphorylation|electron transport chain|tricarboxylic|aerobic respiration|cellular respiration|ATP synthesis|energy derivation|precursor metabolites and energy",
  "RNA / ribosome" = "\\b(?:[a-z]?RNA|rRNA|ncRNA|snRNA)\\b|ribosom|ribonucleo|splic|snRNP|spliceosom|translation",
  "Cilia" = "cili(?:um|a|ary|ogenesis)|axonem",
  "Lipid / sterol" = "cholesterol|sterol|fatty acid|lipid|terpenoid|alcohol biosynthe|alcohol metabol|lipoprotein|glycolipid|liposaccharide",
  "Vesicles / autophagy" = "autophag|lysosom|vacuol|endosom|endocyt|vesicl|Golgi|endomembrane|exocyt|membrane (?:docking|fusion|tether)|organelle (?:fusion|membrane fusion|localization by membrane tethering)|secretion|export from cell|protein localization to vacuole",
  "Ion transport / homeostasis" = "\\b(?:sodium|potassium|calcium|cation|metal ion|iron ion|anion|chloride|zinc ion|proton|pH|ion)\\b.*(?:transport|homeostasis|release)|(?:transport|homeostasis|release|regulation of pH).*\\b(?:sodium|potassium|calcium|cation|metal ion|iron ion|anion|chloride|zinc ion|proton|pH|ions?)\\b",
  "Growth / cell death" = "apoptot|apoptos|cell death|cell killing|cell growth|\\bgrowth\\b(?! factor)|proliferat|maintenance of cell number|cell population|regulation of anatomical structure size|regulation of cellular component size",
  "Vascular / hemostasis" = "angiogen|vascul|vascular|blood vessel|blood pressure|coagulation|hemostasis|haemostasis|hemodynamic|endotheli|tube morphogen",
  "Adhesion / cytoskeleton" = "actin|myosin|myofibril|cytoskel|adhesion|extracellular matrix|extracellular structure|collagen|integrin|migration|motility|locomot|chemotax|junction|myofibril|microtubul|filament|cell polarity|stress fib|Rho protein|muscle contraction|fiber organization",
  "Development / morphology" = "develop|morphogen|differenti|tissue|epithel|skin|epiderm|muscle|cell projection|neuron|neural|axon|synap|neurogen|gliogen|organismal|reproduct|pregnancy|hair cycle|anatomical structure|maturation|stem cell population",
  "Signaling / stress" = "signal|MAPK|ERK|small GTPase|phosphorylat|growth factor|cell activation|cellular response|response to|stress|reactive oxygen|oxidative stress|xenobiotic|mechanical stimulus|Wnt|NF.kappaB|hormone secretion|hormone processing|endocrine",
  "Other metabolism" = "metabol|catabol|biosynthe|glycosyl|proteoly|peptid|amino acid|carbohydrat|glycoprotein|glycosamino|monocarboxylate|organic acid|nucleobase|nucleotide|nucleoside|folic acid|aldehyde|methylation|protein processing"
)

## GO-ID overrides are independently curated exceptions to broad text matching.
## Applied to the same GO identifier in every module; they are NOT model data.
## Review GO-ID exceptions here; rerun stage 3 to regenerate the audit and figures.
go_theme_overrides <- data.frame(
  ID = c("GO:0032467", "GO:0000910", "GO:0051302", "GO:0051781", "GO:0002562", "GO:0002200", "GO:0016445", "GO:0003158", "GO:0072012", "GO:2000108", "GO:1903509", "GO:0030239", "GO:0071560", "GO:0008334", "GO:0060986", "GO:0035239", "GO:2000036", "GO:0019827", "GO:0071396", "GO:0006913", "GO:0046822", "GO:0034504", "GO:0006606", "GO:0006999", "GO:0006885", "GO:0007565", "GO:0019953", "GO:0044703", "GO:0044706", "GO:0007009", "GO:0031579", "GO:0140694", "GO:0000079"),
  theme = c("Cell cycle / genome", "Cell cycle / genome", "Cell cycle / genome", "Cell cycle / genome", "Immune / inflammatory", "Immune / inflammatory", "Immune / inflammatory", "Vascular / hemostasis", "Vascular / hemostasis", "Growth / cell death", "Lipid / sterol", "Adhesion / cytoskeleton", "Signaling / stress", "RNA / ribosome", "Signaling / stress", "Development / morphology", "Development / morphology", "Development / morphology", "Signaling / stress", "Other / unassigned", "Other / unassigned", "Other / unassigned", "Other / unassigned", "Other / unassigned", "Ion transport / homeostasis", "Other / unassigned", "Other / unassigned", "Other / unassigned", "Other / unassigned", "Other / unassigned", "Other / unassigned", "Other / unassigned", "Cell cycle / genome"),
  rationale = c("Cytokinesis is cell division, not cytokine signaling", "Cytokinesis is cell division, not cytokine signaling", "Cell division is explicitly cell-cycle-related", "Positive regulation of cell division", "Somatic immune-receptor diversification", "Somatic immune-receptor diversification", "Somatic immunoglobulin diversification", "Endothelium development", "Glomerulus vasculature development", "Specific regulation of leukocyte apoptosis; cell-death endpoint", "Liposaccharide metabolism is lipid-associated", "Myofibril assembly is contractile cytoskeleton", "Response to TGF-beta ligand", "Histone mRNA metabolic process is post-transcriptional RNA metabolism", "Endocrine hormone release is an endocrine signaling function", "Tube morphogenesis is not specifically vascular without a tissue qualifier", "Stem cell population maintenance denotes a differentiation-state process", "Stem cell population maintenance denotes a differentiation-state process", "Cellular response to lipid is a stimulus response, not lipid metabolism", "Nucleocytoplasmic transport is not itself cell-cycle regulation", "Regulation of nucleocytoplasmic transport is not itself cell-cycle regulation", "Nuclear import is not itself cell-cycle regulation", "Nuclear import is not itself cell-cycle regulation", "Nuclear pore organization is not itself cell-cycle regulation", "pH regulation is ionic homeostasis", "Pregnancy is a broad multicellular term not a specific cellular module process", "Sexual reproduction is not itself cell-cycle regulation", "Multi-organism reproductive process is broad", "Multi-organism multicellular process is broad", "Plasma membrane organization is broader than endomembrane/vesicles", "Membrane raft organization not vesicular trafficking", "Membraneless organelle assembly nonspecific", "Cyclin-dependent kinase activity"),
  stringsAsFactors = FALSE
)
stopifnot(!anyDuplicated(go_theme_overrides$ID),
          all(go_theme_overrides$theme %in%
                c(names(go_theme_patterns), "Other / unassigned")))

## Each semantic-reduced GO ID appears in at most one module row, but the same
## biological term can recur across modules. Apply a single stable assignment
## to a given GO ID, record all competing keyword matches for transparency.
audit_go_theme <- function(term, go_id) {
  if (anyNA(term) || any(!nzchar(trimws(term))) ||
      anyNA(go_id) || any(!nzchar(go_id)))
    stop("GO descriptions and IDs must not be empty")
  if (length(term) != length(go_id)) stop("Mismatched GO terms and IDs")
  result <- lapply(seq_along(term), function(i) {
    description <- term[[i]]
    hits <- lapply(go_theme_patterns, function(pattern) {
      mm <- regexpr(pattern, description, ignore.case = TRUE, perl = TRUE)
      if (mm[1L] < 0L) return(NA_character_)
      substr(description, mm[1L], mm[1L] + attr(mm, "match.length")[1L] - 1L)
    })
    hit_idx <- which(!vapply(hits, is.na, logical(1)))
    priority <- if (length(hit_idx)) hit_idx[1L] else NA_integer_
    category <- if (is.na(priority)) "Other / unassigned" else
      names(go_theme_patterns)[priority]
    method <- if (is.na(priority)) "unclassified" else "first_matching_specific_regex"
    reason <- ""
    override <- match(go_id[[i]], go_theme_overrides$ID)
    if (!is.na(override)) {
      category <- go_theme_overrides$theme[[override]]
      method <- "curated_GO_ID_override"
      reason <- go_theme_overrides$rationale[[override]]
    }
    data.frame(
      theme = category,
      rule_priority = if (identical(method, "curated_GO_ID_override")) NA_integer_ else priority,
      matched_text = if (is.na(priority) || identical(method, "curated_GO_ID_override")) "" else hits[[priority]],
      winning_regex = if (is.na(priority) || identical(method, "curated_GO_ID_override")) "" else unname(go_theme_patterns[priority]),
      all_matching_themes = paste(names(go_theme_patterns)[hit_idx], collapse = "; "),
      all_matching_text = paste(unlist(hits[hit_idx], use.names = FALSE), collapse = "; "),
      n_matching_themes = length(hit_idx),
      classification_method = method,
      override_rationale = reason,
      needs_review = category == "Other / unassigned" ||
        (length(hit_idx) > 1L && !identical(method, "curated_GO_ID_override")),
      audit_flag = if (identical(method, "curated_GO_ID_override")) "CURATED_OVERRIDE" else
        if (!length(hit_idx)) "UNASSIGNED" else if (length(hit_idx) > 1L)
          "MULTIPLE_MATCHES_PRIORITY_USED" else "SINGLE_MATCH",
      stringsAsFactors = FALSE
    )
  })
  audit <- do.call(rbind, result)
  rownames(audit) <- NULL
  audit
}

plot_module_go_views <- function() {
  all_file <- file.path(figure_dir, "GO_BP_ORA_all_terms.csv")
  reduced_file <- file.path(figure_dir, "GO_BP_ORA_nonredundant_terms.csv")
  reduction_status_file <- file.path(figure_dir, "GO_BP_semantic_reduction_status.csv")
  if (!file.exists(all_file)) stop("Missing complete module GO-BP table: ", all_file)
  all <- read_module_csv(all_file, c("ID", "Description", "module", "p.adjust", "GeneRatio", "Count"))
  all <- all[all$module %in% module_order & is.finite(all$p.adjust) &
               all$p.adjust < 0.05, , drop = FALSE]
  if (anyDuplicated(paste(all$ID, all$module))) stop("Duplicate GO/module entry")

  ## The full analysis writes the Wang-simplified term sets. Do not quietly
  ## construct "biological proportions" from thousands of redundant raw terms.
  if (!file.exists(reduction_status_file)) stop(
    "GO theme plots require semantic-reduction status: ", reduction_status_file)
  audit <- read_module_csv(reduction_status_file, c("module", "wang_simplification_success"))
  if (!setequal(audit$module, module_order) || anyNA(audit$wang_simplification_success) ||
      any(!audit$wang_simplification_success)) {
    stop("GO theme plots require complete semantic reduction. See: ", reduction_status_file)
  }
  if (!file.exists(reduced_file)) {
    stop("Missing semantically reduced ORA table: ", reduced_file)
  }
  reduced <- read_module_csv(reduced_file,
                             c("ID", "Description", "module", "p.adjust", "GeneRatio", "Count"))
  reduced <- reduced[reduced$module %in% module_order & is.finite(reduced$p.adjust) &
                       reduced$p.adjust < 0.05, , drop = FALSE]
  if (anyDuplicated(paste(reduced$ID, reduced$module))) stop("Nonunique GO/module in reduced table")
  if (!all(paste(reduced$ID, reduced$module) %in% paste(all$ID, all$module)))
    stop("Reduced GO terms must be a subset of significant full GO results")
  if (!nrow(reduced)) stop("No FDR-significant reduced GO terms")
  ## All occurrences of a GO ID must have the same GO name across modules.
  id_to_label <- split(as.character(reduced$Description), as.character(reduced$ID))
  if (any(vapply(id_to_label, function(z) length(unique(z)) > 1L, logical(1))))
    stop("Same GO ID has inconsistent descriptions across modules")


  ## The category assignment itself does not depend on enrichment significance:
  ## it depends ONLY on the description and the ordered regex dictionary.
  classification <- audit_go_theme(as.character(reduced$Description), as.character(reduced$ID))
  for (field in names(classification)) reduced[[field]] <- classification[[field]]
  reduced$minus_log10_FDR <- -log10(pmax(reduced$p.adjust, .Machine$double.xmin))
  reduced$gene_ratio_percent <- vapply(strsplit(as.character(reduced$GeneRatio), "/", fixed = TRUE),
    function(x) if (length(x) == 2L) 100 * as.numeric(x[1L]) / as.numeric(x[2L]) else NA_real_, numeric(1))
  if (any(!is.finite(reduced$gene_ratio_percent))) stop("Invalid GO GeneRatio")

  ## Separate editable/readable record of the rule ORDER and full regexes.
  ## It is descriptive metadata, NOT an additional ontology analysis.
  dictionary <- data.frame(
    rule_priority = seq_along(go_theme_patterns),
    biological_theme = names(go_theme_patterns),
    regex_pattern = unname(go_theme_patterns),
    match_mode = "case-insensitive PCRE on GO Description",
    assignment_logic = "ordered specific regex, followed by curated GO-ID override",
    stringsAsFactors = FALSE
  )
  utils::write.csv(dictionary, file.path(figure_dir, "GO_BP_theme_keyword_dictionary.csv"),
                   row.names = FALSE)

  ## Curated GO-ID exceptions accompany the machine-readable regex dictionary.
  utils::write.csv(go_theme_overrides,
    file.path(figure_dir, "GO_BP_theme_curated_overrides.csv"), row.names = FALSE)

  ## All themes count terms after per-module semantic redundancy reduction.
  tally <- as.data.frame(table(module = factor(reduced$module, levels = module_order),
                               theme = reduced$theme), stringsAsFactors = FALSE)
  tally <- tally[tally$Freq > 0L, , drop = FALSE]
  names(tally)[names(tally) == "Freq"] <- "n_nonredundant_GO_terms"
  count_total <- stats::aggregate(n_nonredundant_GO_terms ~ module, tally, sum)
  tally$n_module_terms <- count_total$n_nonredundant_GO_terms[match(tally$module, count_total$module)]
  tally$percent_of_terms <- 100 * tally$n_nonredundant_GO_terms / tally$n_module_terms
  tally <- tally[order(tally$module, -tally$n_nonredundant_GO_terms, tally$theme), , drop = FALSE]
  utils::write.csv(tally, file.path(figure_dir, "GO_BP_module_theme_composition.csv"), row.names = FALSE)

  if (nrow(reduced) != sum(tally$n_nonredundant_GO_terms) ||
      any(!is.finite(tally$percent_of_terms))) {
    stop("GO-theme audit totals do not match the reduced GO-term table")
  }

  ## Dominant labels are descriptive; do not call a module "mainly X" when
  ## the largest category is tiny or an unassigned bucket.
  dominant <- do.call(rbind, lapply(module_order, function(m) {
    tab <- tally[as.character(tally$module) == m, , drop = FALSE]
    if (!nrow(tab)) return(data.frame(module = m, n_reduced_terms = 0L,
      dominant_theme = "No significant GO-BP enrichment", fraction_percent = NA_real_,
      status = "No FDR-significant GO-BP terms", stringsAsFactors = FALSE))
    top <- tab[1L, ]
    runner_up <- if (nrow(tab) > 1L) tab$percent_of_terms[2L] else 0
    qualifies <- top$percent_of_terms >= 40 &&
      (top$percent_of_terms - runner_up) >= 15 &&
      top$theme != "Other / unassigned"
    data.frame(module = m, n_reduced_terms = as.integer(top$n_module_terms),
       dominant_theme = as.character(top$theme), fraction_percent = top$percent_of_terms,
       status = if (qualifies) "Predominant theme (descriptive)" else "Mixed / uncertain themes",
       stringsAsFactors = FALSE)
  }))
  utils::write.csv(dominant, file.path(figure_dir, "GO_BP_module_dominant_themes.csv"), row.names = FALSE)
  if (any(dominant$status == "Mixed / uncertain themes"))
    message("GO theme caveat: some modules have mixed / weakly classified ontology themes")

  ## Representative terms for the category composition shown in panel B.
  ## Represent up to four largest *biological* categories per module.
  ## If the leading category makes up >= 45% of reduced terms, include a
  ## second, distinct GO term from that category (maximum four total).
  ## Curated, predeclared GO IDs take priority only when enriched in a
  ## corresponding module; otherwise choose specific, FDR-supported terms.
  ## Generic "regulation of X" labels are deprioritized within each category
  ## where more specific representative processes are available; adjusted
  ## p-value and term size break ties. This is a DISPLAY selection, not a new
  ## enrichment test and not a ranking solely by statistical significance.
  ## These GO IDs were checked against the existing reduced ORA results.
  ## They favor recognizable module-level biology and shared M02/M04 terms.
  ## A preferred ID cannot be shown unless it is significant AND belongs to
  ## one of the four most frequent categories in the module.
  preferred_GO_ids <- c(
    "GO:0006914", "GO:0016050", # autophagy and vesicle organization
    "GO:0051276", "GO:0006260", "GO:0031507", # chromosome/DNA/chromatin
    "GO:0006954", "GO:0098609", "GO:0048729", # shared M02/M04 processes
    "GO:0007169", "GO:0007015", # growth-factor signaling, actin organization
    "GO:0032543", "GO:0006119", # mitochondrial translation, OXPHOS
    "GO:0033489", "GO:0034329", "GO:0050808", "GO:0006814", # sterol/adhesion/synapse (M05)
    "GO:0060271", "GO:0044782", "GO:0044282" # cilia/catabolism (M07)
  )
  preferred_by_module <- list(
    M01 = c("GO:0006914", "GO:0016050"),
    M08 = c("GO:0051276", "GO:0006260"),
    M02 = c("GO:0006954", "GO:0048729", "GO:0007169", "GO:0098609"),
    M06 = c("GO:0032543", "GO:0006119"),
    M05 = c("GO:0034329", "GO:0050808", "GO:0033489", "GO:0006814"),
    M04 = c("GO:0007015", "GO:0048729", "GO:0006954"),
    M07 = c("GO:0060271", "GO:0044282")
  )
  choose_theme_terms <- function(block, n, module) {
    if (!nrow(block)) return(block)
    priority <- unique(c(preferred_by_module[[module]], preferred_GO_ids))
    block$preferred <- match(block$ID, priority)
    block$generic_label <- grepl(
      "^(positive |negative )?regulation of |^(cellular|biological) process$",
      block$Description, ignore.case = TRUE)
    block <- block[order(is.na(block$preferred), block$preferred,
                         block$generic_label, block$p.adjust, -block$Count), , drop = FALSE]
    block$generic_label <- NULL
    block$preferred <- NULL
    head(block, n)
  }
  selected <- do.call(rbind, lapply(module_order, function(m) {
    counts <- tally[as.character(tally$module) == m &
                      tally$theme != "Other / unassigned", , drop = FALSE]
    if (!nrow(counts)) return(NULL)
    ## Do not promote small incidental themes as the module's main biology.
    ## Show a theme only if it represents >=10% of reduced terms or >=10 terms.
    ## This is a DISPLAY criterion; all themes remain visible in panel B/audit.
    counts <- counts[counts$percent_of_terms >= 10 |
                     counts$n_nonredundant_GO_terms >= 10L, , drop = FALSE]
    if (!nrow(counts)) return(NULL)
    leading <- head(counts, 4L)
    blocks <- lapply(seq_len(nrow(leading)), function(i) {
      theme <- as.character(leading$theme[i])
      n <- if (i == 1L && leading$percent_of_terms[i] >= 45 &&
                 leading$n_nonredundant_GO_terms[i] >= 2L) 2L else 1L
      block <- reduced[reduced$module == m & reduced$theme == theme, , drop = FALSE]
      chosen <- choose_theme_terms(block, n, m)
      chosen$theme_n_terms <- leading$n_nonredundant_GO_terms[i]
      chosen$theme_percent <- leading$percent_of_terms[i]
      chosen$selection_rank_in_theme <- seq_len(nrow(chosen))
      chosen
    })
    head(do.call(rbind, blocks), 4L)
  }))
  if (is.null(selected) || !nrow(selected)) stop("No GO terms chosen for display")
  ## A term can be relevant to multiple modules; one axis label is sufficient.
  ## Other significant module associations for that GO term appear as dots.
  selected <- selected[!duplicated(selected$ID), , drop = FALSE]
  utils::write.csv(selected, file.path(figure_dir, "GO_BP_theme_representatives.csv"),
                   row.names = FALSE)

  ## Enrich the EXISTING full term-assignment audit with exact plot membership.
  ## The same GO ID can be significant in multiple modules; only the chosen
  ## home module's GO-ID/module pair is marked selected_for_dotplot.
  audit_key <- paste(reduced$module, reduced$ID, sep = "||")
  selected_key <- paste(selected$module, selected$ID, sep = "||")
  if (anyDuplicated(audit_key) || anyDuplicated(selected_key))
    stop("Duplicate module/GO-ID audit key")
  reduced$selected_for_dotplot <- audit_key %in% selected_key
  selected_ix <- match(audit_key, selected_key)
  reduced$dotplot_theme_rank <- as.integer(selected$selection_rank_in_theme[selected_ix])
  reduced$category_method <- reduced$classification_method
  reduced <- reduced[order(match(reduced$module, module_order),
                           reduced$theme, reduced$p.adjust, reduced$ID), , drop = FALSE]
  ## Canonical comma-delimited CSV remains unchanged for R / reproducibility.
  utils::write.csv(reduced, file.path(figure_dir, "GO_BP_theme_term_assignments.csv"),
                   row.names = FALSE, na = "")

  ## Excel-friendly audit: NEVER store multiple category matches in a single
  ## semicolon-joined field. Explicit, separately named match columns mean
  ## that a term matching two or three themes remains unambiguous in Excel.
  split_audit_matches <- function(value) {
    if (is.na(value) || !nzchar(value)) character() else
      strsplit(value, "; ", fixed = TRUE)[[1L]]
  }
  n_match_slots <- max(1L, max(reduced$n_matching_themes, na.rm = TRUE))
  audit_excel <- reduced[, c("module", "ID", "Description", "theme",
    "rule_priority", "matched_text", "n_matching_themes", "audit_flag",
    "classification_method", "override_rationale", "needs_review",
    "p.adjust", "Count", "GeneRatio", "gene_ratio_percent",
    "selected_for_dotplot", "dotplot_theme_rank"), drop = FALSE]
  for (slot in seq_len(n_match_slots)) {
    audit_excel[[paste0("matched_theme_", slot)]] <- vapply(
      reduced$all_matching_themes, function(value) {
        x <- split_audit_matches(value)
        if (length(x) >= slot) x[[slot]] else ""
      }, character(1))
    audit_excel[[paste0("matched_keyword_", slot)]] <- vapply(
      reduced$all_matching_text, function(value) {
        x <- split_audit_matches(value)
        if (length(x) >= slot) x[[slot]] else ""
      }, character(1))
  }
  ## CSV2 follows German Excel defaults: semicolon-separated fields,
  ## decimal commas, and quoted text containing delimiters.
  utils::write.csv2(audit_excel,
    file.path(figure_dir, "GO_BP_theme_term_assignments_excel_DE.csv"),
    row.names = FALSE, na = "", fileEncoding = "UTF-8")

  ## Normalized rule-level audit: ONE ROW PER MATCHED TERM X RULE.
  ## More than one row does not mean more than one assigned pie category;
  ## only is_assigned_theme = TRUE defines the winner for each GO term.
  rule_matches <- lapply(seq_len(nrow(reduced)), function(i) {
    themes <- split_audit_matches(reduced$all_matching_themes[[i]])
    literal_hits <- split_audit_matches(reduced$all_matching_text[[i]])
    if (!length(themes)) return(NULL)
    if (length(themes) != length(literal_hits))
      stop("Mismatched audit matches for GO term ", reduced$ID[[i]])
    data.frame(
      module = reduced$module[[i]], ID = reduced$ID[[i]],
      Description = reduced$Description[[i]],
      matching_theme = themes,
      rule_priority = match(themes, names(go_theme_patterns)),
      matched_keyword = literal_hits,
      is_assigned_theme = themes == reduced$theme[[i]],
      assigned_theme = reduced$theme[[i]],
      stringsAsFactors = FALSE)
  })
  rule_matches <- do.call(rbind, rule_matches[!vapply(rule_matches, is.null, logical(1))])
  if (is.null(rule_matches)) {
    rule_matches <- data.frame(module = character(), ID = character(),
      Description = character(), matching_theme = character(),
      rule_priority = integer(), matched_keyword = character(),
      is_assigned_theme = logical(), assigned_theme = character())
  }
  stopifnot(nrow(rule_matches) == sum(reduced$n_matching_themes))
  utils::write.csv(rule_matches,
    file.path(figure_dir, "GO_BP_theme_rule_matches_long.csv"), row.names = FALSE)

  ## Mandatory: one auditable Excel workbook, not optional/silent.
  audit_xlsx <- file.path(figure_dir, "GO_BP_theme_category_audit.xlsx")
  ekbSeq::export_excel(
    list(Term_assignments = audit_excel,
         Module_composition = tally,
         All_rule_matches = rule_matches,
         Keyword_rules = dictionary,
         Curated_overrides = go_theme_overrides),
    file = audit_xlsx, styled = TRUE)
  if (!file.exists(audit_xlsx) || file.info(audit_xlsx)$size == 0L)
    stop("GO audit XLSX missing or empty after export: ", audit_xlsx)
  flagged <- sum(reduced$needs_review)
  message("GO category audit: ", nrow(reduced), " module/GO-term assignments; ",
          flagged, " flagged for review; curated GO-ID overrides are already resolved. See ",
          file.path(figure_dir, "GO_BP_theme_term_assignments.csv"))

  edges <- all[all$ID %in% selected$ID, , drop = FALSE]
  edges$theme <- selected$theme[match(edges$ID, selected$ID)]
  edges$gene_ratio_percent <- vapply(strsplit(as.character(edges$GeneRatio), "/", fixed = TRUE),
     function(x) if (length(x) == 2L) 100 * as.numeric(x[1L]) / as.numeric(x[2L]) else NA_real_, numeric(1))
  edges$minus_log10_FDR <- pmin(25, -log10(pmax(edges$p.adjust, .Machine$double.xmin)))
  selected <- selected[order(match(selected$module, module_order)), , drop = FALSE]
  selected$term <- as.character(selected$Description)
  dup <- duplicated(selected$term) | duplicated(selected$term, fromLast = TRUE)
  selected$term[dup] <- paste0(selected$term[dup], " (", selected$ID[dup], ")")
  edges$term <- selected$term[match(edges$ID, selected$ID)]
  edges$module <- factor(edges$module, levels = module_order)
  edges$term <- factor(edges$term, levels = rev(selected$term))

  p_dot <- ggplot2::ggplot(edges,
    ggplot2::aes(x = module, y = term, size = gene_ratio_percent, fill = minus_log10_FDR)) +
    ggplot2::geom_point(shape = 21, color = "#243B55", stroke = .30) +
    ggplot2::scale_fill_gradient(low = "#DCEAF8", high = "#12386F",
      name = "-log10(BH FDR)\n(capped at 25)") +
    ggplot2::scale_size_area(max_size = 7.6, name = "Gene ratio (%)") +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::labs(title = "A  |  Representative GO biological processes",
      subtitle = "Representative processes from up to four leading themes/module; BH FDR < 0.05",
      x = "Response module", y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "#EDF0F4", linewidth = .28),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold", size = 9),
      axis.text.y = ggplot2::element_text(size = 8.8),
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      legend.position = "bottom")
  export_module_plot(p_dot, "module_GO_BP_dotplot", 12.5, 9.5)

  ## Retain the earlier optional module-GO bipartite network as an alternative.
  module_nodes <- data.frame(module = module_order,
    y = seq(100, 0, length.out = length(module_order)), stringsAsFactors = FALSE)
  selected$term_y <- NA_real_
  for (m in module_order) {
    ix <- which(selected$module == m)
    if (length(ix)) selected$term_y[ix] <- module_nodes$y[match(m, module_nodes$module)] +
      seq(4, -4, length.out = length(ix))
  }
  selected$n_modules <- vapply(selected$ID, function(id)
    length(unique(as.character(edges$module[edges$ID == id]))), integer(1))
  edges$mod_y <- module_nodes$y[match(as.character(edges$module), module_nodes$module)]
  edges$term_y <- selected$term_y[match(as.character(edges$term), selected$term)]
  module_nodes$color <- unname(module_colors[module_nodes$module])
  p_network <- ggplot2::ggplot() +
    ggplot2::geom_curve(data = edges,
      ggplot2::aes(x = .05, y = mod_y, xend = .97, yend = term_y,
                   color = module, linewidth = minus_log10_FDR),
      curvature = .11, alpha = .60, show.legend = FALSE) +
    ggplot2::geom_point(data = module_nodes,
      ggplot2::aes(x = 0, y = y, fill = module), shape = 21,
      size = 7.6, color = "#374151", stroke = .6, show.legend = FALSE) +
    ggplot2::geom_text(data = module_nodes,
      ggplot2::aes(x = 0, y = y, label = module), fontface = "bold", size = 3.5) +
    ggplot2::geom_point(data = selected,
      ggplot2::aes(x = 1, y = term_y),
      color = ifelse(selected$n_modules > 1L, "#A52A39", "#34495E"),
      size = 3.9, show.legend = FALSE) +
    ggplot2::geom_text(data = selected,
      ggplot2::aes(x = 1.075, y = term_y, label = term),
      hjust = 0, size = 3.2, color = "#243141") +
    ggplot2::scale_color_manual(values = module_colors, guide = "none") +
    ggplot2::scale_fill_manual(values = module_colors, guide = "none") +
    ggplot2::scale_linewidth_continuous(range = c(.32, 2.0), guide = "none") +
    ggplot2::coord_cartesian(xlim = c(-.22, 2.9), ylim = c(-9, 109), clip = "off") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.margin = ggplot2::margin(15, 20, 18, 8)) +
    ggplot2::labs(title = "Module-GO association network",
      subtitle = "Edges: significant GO-BP enrichment (BH FDR < 0.05)",
      caption = "Theme-selected GO representatives; shared terms can connect multiple modules. No causal network.")
  export_module_plot(p_network, "module_GO_BP_network", 14.5, 10.5)

  ## Each miniature pie is a distribution of nonredundant GO-TERM COUNTS,
  ## not a quantitative share of responding genes. M03 has no pie at all.
  theme_levels <- c(names(go_theme_patterns), "Other / unassigned")
  theme_colors <- c(
    "Cell cycle / genome" = "#3C5488", "Mitochondria / energy" = "#7E6148",
    "RNA / ribosome" = "#8491B4", "Cilia" = "#00A087",
    "Lipid / sterol" = "#E6A742", "Vesicles / autophagy" = "#8F63A7",
    "Immune / inflammatory" = "#E64B35", "Vascular / hemostasis" = "#C95677",
    "Adhesion / cytoskeleton" = "#4DBBD5", "Development / morphology" = "#91D1C2",
    "Signaling / stress" = "#F39B7F", "Other metabolism" = "#B79D75",
    "Growth / cell death" = "#D685A2",
    "Ion transport / homeostasis" = "#55A6A4",
    "Other / unassigned" = "#C9CDD2")
  stopifnot(setequal(names(theme_colors), theme_levels))
  tally$module <- factor(as.character(tally$module), levels = module_order)
  tally$theme <- factor(tally$theme, levels = theme_levels)
  tally <- tally[order(tally$module, tally$theme), , drop = FALSE]
  tally$ymax <- ave(tally$n_nonredundant_GO_terms, tally$module,
                    FUN = function(x) cumsum(x) / sum(x))
  tally$ymin <- tally$ymax - tally$percent_of_terms / 100
  tally$ymid <- (tally$ymin + tally$ymax) / 2

  ## Label sufficiently wide slices with the absolute NUMBER OF REDUCED
  ## GO-BP TERMS in the corresponding biological category. Unlabeled small
  ## slices are still represented in the legend and complete audit CSV.
  tally$slice_label <- ifelse(tally$percent_of_terms >= 9,
                              as.character(tally$n_nonredundant_GO_terms), "")
  col_rgb <- grDevices::col2rgb(unname(theme_colors[as.character(tally$theme)]))
  luminance <- colSums(col_rgb * c(0.2126, 0.7152, 0.0722)) / 255
  tally$label_color <- ifelse(luminance < 0.55, "white", "#1B2838")
  pie_totals <- dominant[dominant$n_reduced_terms > 0L,
                         c("module", "n_reduced_terms"), drop = FALSE]
  pie_totals$module <- factor(pie_totals$module, levels = module_order)
  pie_totals$center_label <- paste0("n=", pie_totals$n_reduced_terms)

  pie <- ggplot2::ggplot(tally) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = .34, xmax = 1.0, ymin = ymin, ymax = ymax, fill = theme),
      color = "white", linewidth = .26) +
    ggplot2::geom_text(
      data = tally[nzchar(tally$slice_label), , drop = FALSE],
      ggplot2::aes(x = .70, y = ymid, label = slice_label, color = label_color),
      fontface = "bold", size = 3, show.legend = FALSE) +
    ggplot2::geom_text(
      data = pie_totals,
      ggplot2::aes(x = 0.06, y = .5, label = center_label),
      inherit.aes = FALSE, fontface = "bold", size = 2.85, color = "#25364A") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::facet_wrap(~module, ncol = 2, drop = FALSE) +
    ggplot2::scale_fill_manual(values = theme_colors, drop = FALSE,
                               name = "Biological theme") +
    ggplot2::scale_color_identity(guide = "none") +
    ggplot2::xlim(0, 1.15) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 10),
                   legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 7),
                   plot.title = ggplot2::element_text(face = "bold", size = 13)) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2, byrow = TRUE)) +
    ggplot2::labs(title = "B  |  Functional themes",
      subtitle = "Counts of reduced GO-BP terms; center = module total")
  ## A genuinely empty module must remain empty; it has zero significant terms.
  if ("M03" %in% module_order) {
    pie <- pie + ggplot2::geom_text(
      data = data.frame(module = factor("M03", levels = module_order),
                        x = .6, y = .5, label = "n=0\nNo enriched\nGO-BP terms"),
      ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE,
      size = 2.9, color = "grey40")
  }
  export_module_plot(pie, "module_GO_BP_theme_pies", 9.0, 12.0)
  combo <- (p_dot + pie + patchwork::plot_layout(widths = c(1.7, 1.3))) +
    patchwork::plot_annotation(title = "Biological characterization of radiation-response modules",
      subtitle = "Selected representative GO-BP terms and all semantically reduced significant GO-BP terms",
      caption = "Pie slices count semantically reduced GO-BP terms (FDR < 0.05), not genes or pathway activity. Keyword categories include documented GO-ID overrides; small slices are unlabeled.")
  export_module_plot(combo, "module_GO_BP_dotplot_theme_pies", 18.0, 12.4)
  message("GO biological themes: term-level assignments + summary exported under ", figure_dir)
  invisible(list(dotplot = p_dot, pies = pie, dominant = dominant))
}

plot_module_go_views()
