#!/usr/bin/env Rscript
# 09_DREAM_analysis.R
# p21 irradiation study | targeted DREAM analysis
# Input: published Engeland et al. (2022) DREAM targets; six KO-WT interaction
# contrasts (Wald + ashr); twelve ashr irradiation-versus-control contrasts.
# All DESeq2 results are read as-is. No model is refitted.
#
# Run at the p21 project root (see README_DREAM_09.md for options).
# GSEA: all genes with finite ORIGINAL Wald statistics; BH across six contexts.
# Effects and biological response plots: ashr-shrunken log2FC.

# ANALYSIS ROADMAP
# Question: Are the published DREAM targets enriched among genes whose response
# to irradiation changes when p21 is lost? Evaluate each cell line / radiation
# separately; then describe the WT and KO responses of the same target genes.
#
# Effect definition for every interaction (log2 scale):
#   delta_p21 = (KO_irradiated - KO_control) - (WT_irradiated - WT_control)
# Positive delta: more positive / less negative irradiation response in KO.
# Negative delta: less positive / more negative irradiation response in KO.
# This is an adjusted difference-in-differences estimate from the upstream model;
# here we import those contrasts and DO NOT refit or re-estimate the model.
#
# Distinct evidence layers must not be conflated:
#   (1) GSEA: enrichment of a predefined gene set in signed Wald rankings.
#   (2) Gene-level interaction: adjusted differential irradiation response.
#   (3) WT/KO response: within-genotype irradiation-vs-control effect.
#   (4) GO composition: descriptive annotations of GSEA leading-edge genes.
# None is a direct assay of DREAM complex binding or transcriptional activity.
options(stringsAsFactors = FALSE)
# fgseaMultilevel may use random sampling; fix RNG for repeatable results.
# Identical output also depends on fgsea/R versions and identical input files.
set.seed(21)

# ---- 00 | Runtime and input configuration -----------------------------------

# CLI entry point: keep inputs overridable, but use explicit project defaults.
# --skip-go changes descriptive GO annotations only; it does NOT skip GSEA.
parse_options <- function(args) {
  opts <- list(root = ".", targets = NULL, interactions = NULL,
               responses = NULL, output = NULL, skip_go = FALSE)
  if ("--help" %in% args) {
    cat("Usage: Rscript R/09_DREAM_analysis.R [--root PATH]",
        "[--targets CSV/XLSX] [--interactions DIR]",
        "[--responses DIR]",
        "[--output DIR] [--skip-go]\n")
    quit(status = 0)
  }
  if ("--skip-go" %in% args) {
    opts$skip_go <- TRUE
    args <- args[args != "--skip-go"]
  }
  if (length(args) %% 2L != 0L) stop("Options require a value; see --help.")
  if (!length(args)) return(opts)
  for (i in seq(1L, length(args), by = 2L)) {
    key <- sub("^--", "", args[i]); val <- args[i + 1L]
    if (!startsWith(args[i], "--") || !key %in% names(opts) || key == "skip_go") {
      stop("Unknown option: ", args[i])
    }
    opts[[key]] <- val
  }
  opts
}

# Fail before analysis if mandatory packages are absent. The GO database is
# only necessary when GO-based leading-edge categories are requested.
check_packages <- function(skip_go) {
  packages <- c("fgsea", "ggplot2", "patchwork", "ekbSeq", "openxlsx", "readxl", "ggsci")
  if (!skip_go) packages <- c(packages, "AnnotationDbi", "org.Hs.eg.db")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Required packages unavailable: ", paste(missing, collapse = ", "),
         "\nInstall CRAN/Bioconductor packages and ekbSeq 1.2.x first.",
         "\nUse --skip-go to omit optional GO-based composition annotations.")
  }
  invisible(TRUE)
}

resolve_path <- function(root, path) {
  if (is.null(path)) return(NULL)
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) return(path)
  file.path(root, path)
}

# p21.Rproj is a convenience check, not a hard precondition: analysts may
# supply absolute paths from a different working directory.
get_project_root <- function(root) {
  path <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(path, "p21.Rproj"))) return(path)
  if (basename(path) == "R" && file.exists(file.path(dirname(path), "p21.Rproj"))) {
    return(dirname(path))
  }
  message("No p21.Rproj at root: ", path, " (explicit input paths still supported).")
  path
}

opts <- parse_options(commandArgs(trailingOnly = TRUE))
root <- get_project_root(opts$root)
check_packages(opts$skip_go)

# Use the regular p21 results directories, not the temporary upload archives.
# Raw Wald and ashr-shrunken interaction tables are produced by the same step;
# genotype-specific irradiation responses are read from primary DE results.
# Override a path only when your LOCAL project structure differs.
targets_file <- resolve_path(root, if (is.null(opts$targets))
  "Data/Engeland_DREAM_targets_updated.xlsx" else opts$targets)
interaction_dir <- resolve_path(root, if (is.null(opts$interactions))
  "Results/04_p21_interactions" else opts$interactions)
response_dir <- resolve_path(root, if (is.null(opts$responses))
  "Results/02_primary_DE" else opts$responses)

if (!file.exists(targets_file)) stop("Missing DREAM target table: ", targets_file)
if (!dir.exists(interaction_dir)) stop("Missing interaction results directory: ", interaction_dir)
if (!dir.exists(response_dir)) stop("Missing primary DE results directory: ", response_dir)
out_dir <- if (!is.null(opts$output)) resolve_path(root, opts$output) else
  file.path(root, "Results", "09_DREAM_analysis")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 01 | Read tables, assert unique Ensembl identifiers ---------------------

# Read the original CSV exports produced by the upstream analysis steps.
# An exact filename must occur exactly once within the designated results
# directory; duplicate copies or missing exports must not be selected silently.
source_table <- function(source, filename) {
  matches <- list.files(source, pattern = "\\.csv$", recursive = TRUE,
                        full.names = TRUE)
  matches <- matches[basename(matches) == filename]
  if (length(matches) != 1L) {
    stop("Expected exactly one ", filename, " under ", source,
         "; found ", length(matches),
         ". Check the upstream results or specify the correct results directory.")
  }
  utils::read.csv(matches, check.names = FALSE)
}

# Literature-based target universe is defined independently of the present
# irradiation data. Keep gene identifiers and publication symbols separate.
# Version suffixes on Ensembl IDs are removed for joining to DESeq2 outputs;
# duplicate IDs after normalization are considered an input error.
read_targets <- function(path) {
  if (!file.exists(path)) stop("Target file does not exist: ", path)
  if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    # The curated target workbook stores the annotation on its first worksheet.
    df <- as.data.frame(readxl::read_excel(path, sheet = 1))
  } else {
    df <- utils::read.csv(path, check.names = FALSE)
  }
  needed <- c("Ensembl_ID", "DREAM")
  if (!all(needed %in% names(df))) stop("Target list lacks: ",
                                        paste(setdiff(needed, names(df)), collapse = ", "))
  df$ENSEMBL <- sub("\\..*$", "", as.character(df$Ensembl_ID))
  if (anyNA(df$ENSEMBL) || anyDuplicated(df$ENSEMBL)) {
    stop("Missing/duplicate Ensembl IDs in Engeland targets.")
  }
  if ("Gene_publication" %in% names(df)) {
    df$publication_symbol <- as.character(df$Gene_publication)
  } else if ("Gene" %in% names(df)) {
    df$publication_symbol <- as.character(df$Gene)
  } else {
    df$publication_symbol <- df$ENSEMBL
  }
  # Publication spreadsheet uses blank cells (NA), not explicit zeroes, for
  # targets without DREAM evidence. Preserve that meaning during filtering.
  flag <- trimws(as.character(df$DREAM))
  if (any(!is.na(df$DREAM) & !flag %in% c("", "0", "1"))) {
    stop("Unexpected values in Engeland DREAM column (expected 1/0/blank).")
  }
  df$DREAM <- ifelse(is.na(df$DREAM) | flag == "", 0L, as.integer(flag))
  df <- df[df$DREAM == 1L, c("ENSEMBL", "publication_symbol"), drop = FALSE]
  if (!nrow(df)) stop("No DREAM-positive genes found.")
  message("Publication DREAM-positive targets: ", nrow(df))
  df
}

# Shared import contract: verify required columns and unique Ensembl keys
# before merging or matching results. An invalid gene universe must fail early.
read_contrast <- function(source, filename, fields) {
  df <- source_table(source, filename)
  missing <- setdiff(fields, names(df))
  if (length(missing)) stop("Missing columns in ", filename, ": ",
                            paste(missing, collapse = ", "))
  df$ENSEMBL <- sub("\\..*$", "", as.character(df$ENSEMBL))
  if (anyNA(df$ENSEMBL) || anyDuplicated(df$ENSEMBL)) {
    stop("Missing or duplicate Ensembl IDs in ", filename)
  }
  df[, fields, drop = FALSE]
}

# Prespecified comparison family: 2 cell lines x 3 radiation qualities.
# The explicit ordering is used consistently in tables and plots.
contexts <- expand.grid(cell_line = c("HT1080", "MCF10A"),
                        radiation = c("7gy", "UVA", "UVB"),
                        stringsAsFactors = FALSE)
contexts <- contexts[order(match(contexts$cell_line, c("HT1080", "MCF10A")),
                           match(contexts$radiation, c("7gy", "UVA", "UVB"))), ]
contexts$contrast <- paste(contexts$cell_line, contexts$radiation, sep = "_")

# Hold interaction tables, within-genotype responses and signed GSEA rankings
# separately: each serves a different interpretative purpose downstream.
targets <- read_targets(targets_file)
all_interactions <- list()
all_responses <- list()
wald_ranks <- list()

for (i in seq_len(nrow(contexts))) {
  c_id <- contexts$contrast[i]
  cell <- contexts$cell_line[i]
  rad <- contexts$radiation[i]
  message("Loading ", c_id, " ...")

  # One interaction per cell line/radiation, imported twice:
  # original DESeq2 results provide stat, uncertainty and significance;
  # ashr provides shrunken log2FC for biological effect-size displays.
  # Shrinkage must not silently replace the original Wald ranking statistic.
  fn_raw <- paste0(c_id, "_p21_interaction_unshrunken.csv")
  fn_shr <- paste0(c_id, "_p21_interaction_ashr.csv")
  raw <- read_contrast(interaction_dir, fn_raw,
    c("ENSEMBL", "baseMean", "log2FoldChange", "lfcSE", "stat",
      "pvalue", "padj", "SYMBOL"))
  shr <- read_contrast(interaction_dir, fn_shr,
    c("ENSEMBL", "log2FoldChange", "pvalue", "padj"))
  names(shr)[2:4] <- c("interaction_ashr", "shrunken_p", "shrunken_padj")
  # Strictly require the same Ensembl universe before an inner join. Then
  # verify that the p values and BH-adjusted p values were not altered during
  # LFC shrinkage or mixed up with another contrast/version of the results.
  if (!setequal(raw$ENSEMBL, shr$ENSEMBL)) stop("ID universe mismatch: ", c_id)
  merged <- merge(raw, shr, by = "ENSEMBL", all = FALSE, sort = FALSE)
  ok_p <- isTRUE(all.equal(merged$pvalue, merged$shrunken_p,
                           check.attributes = FALSE, tolerance = 1e-10))
  ok_q <- isTRUE(all.equal(merged$padj, merged$shrunken_padj,
                           check.attributes = FALSE, tolerance = 1e-10))
  if (!ok_p || !ok_q) stop("Raw/shrunken inference mismatch in ", c_id)
  merged$contrast <- c_id
  all_interactions[[c_id]] <- merged
  # GSEA input universe = ALL genes with finite original Wald statistics.
  # Do not preselect significant genes or rank by ashr log2FC: that would
  # change the hypothesis addressed by this preranked enrichment analysis.
  # Positive rank => positive (KO minus WT) *interaction*, not necessarily
  # induction relative to the untreated control in either genotype.
  rankable <- merged[is.finite(merged$stat), c("ENSEMBL", "stat")]
  ranked <- rankable$stat
  names(ranked) <- rankable$ENSEMBL
  wald_ranks[[c_id]] <- sort(ranked, decreasing = TRUE)

  for (genotype in c("wt", "p21_neg")) {
  # Also import the two constituent irradiation-vs-control effects (WT/KO).
  # These 12 contrasts explain *how* each significant interaction arises;
  # a significant effect in only one genotype is not itself an interaction.
    fname <- paste0(cell, "_", genotype, "_", rad, "_vs_control_ashr.csv")
    response <- read_contrast(response_dir, fname,
      c("ENSEMBL", "log2FoldChange", "padj", "baseMean"))
    key <- paste(c_id, genotype, sep = "_")
    all_responses[[key]] <- response
  }
}

# ---- 02 | Original-Wald GSEA: one prespecified set, six contexts ------------

# A single prespecified DREAM set is tested against one signed genome-wide
# interaction ranking. Intersect targets with actually rankable genes: genes
# without valid statistics cannot participate in this context-specific test.
# fgsea returns normalized enrichment score (NES) and its leading edge.
run_one_gsea <- function(stats, dream_ids) {
  overlap <- intersect(dream_ids, names(stats))
  if (length(overlap) < 10L) stop("Fewer than 10 rankable DREAM genes.")
  # A two-sided scoreType permits enrichment at either tail of the ranking.
  # The sign of the NES is tied to the INTERACTION direction defined above.
  # Set-size limits are algorithm parameters, not post-hoc gene selection.
  fg <- fgsea::fgseaMultilevel(pathways = list(DREAM = overlap), stats = stats,
                              minSize = 10, maxSize = 1000, eps = 0,
                              scoreType = "std", nPermSimple = 10000)
  if (nrow(fg) != 1L) stop("Unexpected fgsea result: DREAM missing.")
  list(pvalue = as.numeric(fg$pval[1]), NES = as.numeric(fg$NES[1]),
       ES = as.numeric(fg$ES[1]),
       leading_edge = as.character(fg$leadingEdge[[1]]),
       set_size = as.integer(fg$size[1]))
}

# Primary analysis: finite Wald statistics, without a baseMean cutoff.
# Sensitivity analysis: repeat on genes with baseMean >= 10. This tests
# dependence on low-expression genes, NOT reproducibility in independent data.
results <- list()
for (i in seq_len(nrow(contexts))) {
  key <- contexts$contrast[i]
  ranks <- wald_ranks[[key]]
  primary <- run_one_gsea(ranks, targets$ENSEMBL)
  full <- all_interactions[[key]]
  expressed <- full[is.finite(full$stat) & is.finite(full$baseMean) &
                      full$baseMean >= 10, c("ENSEMBL", "stat")]
  stats_expressed <- expressed$stat
  names(stats_expressed) <- expressed$ENSEMBL
  sensitivity <- run_one_gsea(sort(stats_expressed, decreasing = TRUE), targets$ENSEMBL)
  results[[key]] <- list(primary = primary, sensitivity = sensitivity)
}

# One summary row per contrast; retain the ranked-universe size for auditing.
# Because the GSEA background differs by valid Wald statistics, NES values
# should be read with their context-specific set size and test statistics.
gsea <- do.call(rbind, lapply(seq_len(nrow(contexts)), function(i) {
  key <- contexts$contrast[i]
  a <- results[[key]]$primary; b <- results[[key]]$sensitivity
  data.frame(contrast = key, cell_line = contexts$cell_line[i],
             radiation = contexts$radiation[i],
             NES = a$NES, ES = a$ES, pvalue = a$pvalue,
             set_size = a$set_size, n_ranked = length(wald_ranks[[key]]),
             n_leading_edge = length(a$leading_edge),
             NES_baseMean10 = b$NES, p_baseMean10 = b$pvalue,
             n_ranked_baseMean10 = sum(is.finite(all_interactions[[key]]$stat) &
                                        all_interactions[[key]]$baseMean >= 10))
}))
# Primary multiple-testing family = six DREAM GSEA tests (1 set x 6 contexts).
# Apply BH separately to the six sensitivity p values. These FDR values are
# NOT the per-gene DESeq2 padj values imported above.
gsea$FDR_six <- stats::p.adjust(gsea$pvalue, method = "BH")
gsea$FDR_six_baseMean10 <- stats::p.adjust(gsea$p_baseMean10, method = "BH")
gsea$significant <- !is.na(gsea$FDR_six) & gsea$FDR_six < 0.05

dream_ids <- targets$ENSEMBL
if (nrow(gsea) != 6L) stop("Expected exactly six DREAM GSEA tests.")
message("Significant DREAM GSEAs (BH across six): ", sum(gsea$significant))

# ---- 03 | All DREAM gene responses (12 treatment effects and six deltas) ----

# Produce a full target-by-context table, including publication DREAM genes
# missing from a DESeq2 table. match() preserves target ordering and returns
# NA for unavailable features rather than dropping them from the export.
assemble_genes <- function(key) {
  raw <- all_interactions[[key]]
  sub <- raw[match(dream_ids, raw$ENSEMBL), , drop = FALSE]
  sub$ENSEMBL <- dream_ids
  sub$contrast <- key # Keep the context ID even for publication targets absent from DESeq2 results.
  sub$publication_symbol <- targets$publication_symbol
  # The interaction FDR comes from the imported gene-wise DESeq2 contrast.
  # Classify the sign using its original Wald log2FC; ashr shrinkage is used
  # later for reported magnitudes. An NA padj is not called significant.
  sub$interaction_significant <- !is.na(sub$padj) & sub$padj < 0.05
  sub$interaction_direction <- ifelse(!sub$interaction_significant, "not_significant",
                                      ifelse(sub$log2FoldChange > 0,
                                             "positive_KO_minus_WT", "negative_KO_minus_WT"))
  cmeta <- contexts[contexts$contrast == key, ]
  for (genotype in c("wt", "p21_neg")) {
    response <- all_responses[[paste(key, genotype, sep = "_")]]
    response <- response[match(dream_ids, response$ENSEMBL), , drop = FALSE]
    prefix <- if (genotype == "wt") "WT" else "KO"
    sub[[paste0(prefix, "_response_ashr")]] <- response$log2FoldChange
    sub[[paste0(prefix, "_response_padj")]] <- response$padj
    sub[[paste0(prefix, "_baseMean")]] <- response$baseMean
  }
  sub$cell_line <- cmeta$cell_line
  sub$radiation <- cmeta$radiation
  # Genotype-specific response calls require BOTH FDR < 0.05 and an effect
  # of the respective sign. No minimum absolute log2FC cutoff is imposed.
  # This answers whether a response is detectable *within* each genotype;
  # it does not establish a between-genotype difference by itself.
  sub$WT_repressed <- !is.na(sub$WT_response_padj) &
    sub$WT_response_padj < 0.05 & sub$WT_response_ashr < 0
  sub$KO_repressed <- !is.na(sub$KO_response_padj) &
    sub$KO_response_padj < 0.05 & sub$KO_response_ashr < 0
  sub$WT_induced <- !is.na(sub$WT_response_padj) &
    sub$WT_response_padj < 0.05 & sub$WT_response_ashr > 0
  sub$KO_induced <- !is.na(sub$KO_response_padj) &
    sub$KO_response_padj < 0.05 & sub$KO_response_ashr > 0
  # Descriptive repression patterns are useful for summarizing the target set.
  # In particular, "WT only" means KO was not FDR-repressed; it does NOT
  # prove the KO effect equals zero or differs significantly from the WT effect.
  sub$repression_pattern <- ifelse(sub$WT_repressed & sub$KO_repressed,
    "WT_and_KO_repressed", ifelse(sub$WT_repressed, "WT_repressed_only",
    ifelse(sub$KO_repressed, "KO_repressed_only", "neither_FDR_repressed")))
  # More precise interaction interpretation requires BOTH response and delta.
  # Regulatory labels add the formal interaction result to the response calls.
  # Labels are assigned in order: the broad "preserved" category is refined
  # when the interaction supports weaker/stronger KO repression. Remaining
  # significant interactions can have other signs or effect combinations.
  # These are phenotype patterns, not direct evidence of DREAM regulation.
  sub$regulatory_pattern <- "other_or_not_resolved"
  sub$regulatory_pattern[sub$WT_repressed & sub$KO_repressed] <- "repression_preserved"
  sub$regulatory_pattern[sub$WT_repressed & sub$KO_repressed &
                           sub$interaction_significant & sub$log2FoldChange > 0] <-
    "repression_preserved_but_attenuated_in_KO"
  sub$regulatory_pattern[sub$WT_repressed & sub$KO_repressed &
                           sub$interaction_significant & sub$log2FoldChange < 0] <-
    "repression_preserved_and_stronger_in_KO"
  sub$regulatory_pattern[sub$WT_repressed & !sub$KO_repressed &
                           sub$interaction_significant & sub$log2FoldChange > 0] <-
    "WT_repression_attenuated_KO_not_FDR_repressed"
  sub$regulatory_pattern[!sub$WT_repressed & sub$KO_repressed &
                           sub$interaction_significant & sub$log2FoldChange < 0] <-
    "KO_repression_stronger_WT_not_FDR_repressed"
  sub$regulatory_pattern[!sub$WT_repressed & !sub$KO_repressed &
                           sub$interaction_significant] <- "other_significant_interaction"
  # Leading-edge membership is defined algorithmically for EACH GSEA,
  # even when the overall GSEA fails BH significance; plotting below restricts
  # biological emphasis to significant contexts.
  sub$leading_edge <- sub$ENSEMBL %in% results[[key]]$primary$leading_edge
  sub$dream_target <- TRUE
  sub
}

dream_all <- do.call(rbind, lapply(contexts$contrast, assemble_genes))
rownames(dream_all) <- NULL

# Gene-level summary: patterns are descriptions, not evidence of DREAM activity.
patterns <- aggregate(rep(1L, nrow(dream_all)),
  by = list(contrast = dream_all$contrast,
            regulatory_pattern = dream_all$regulatory_pattern), FUN = sum)
names(patterns)[3] <- "n_genes"

# Check target coverage and potential low-count sensitivity by context.
# A low baseMean or missing padj describes the data; it is not an automatic
# justification for removing a gene from the primary prespecified gene set.
low_qc <- do.call(rbind, lapply(contexts$contrast, function(key) {
  a <- dream_all[dream_all$contrast == key, ]
  data.frame(contrast = key, targets_total = nrow(a),
             targets_with_wald = sum(is.finite(a$stat)),
             targets_baseMean_lt10 = sum(!is.na(a$baseMean) & a$baseMean < 10),
             targets_zero_baseMean = sum(a$baseMean == 0, na.rm = TRUE),
             targets_padj_NA = sum(is.na(a$padj)),
             targets_sig_interaction = sum(a$interaction_significant),
             WT_FDR_repressed = sum(a$WT_repressed, na.rm = TRUE),
             KO_FDR_repressed = sum(a$KO_repressed, na.rm = TRUE))
}))

# ---- 04 | GO-BP-anchored leading-edge themes (exploratory only) -------------

# GOALL includes ancestral Biological Process terms. Preserve the frequent
# replication/DDR overlap explicitly instead of assigning all double-annotated
# genes to replication. Display classes remain mutually exclusive; this is
# descriptive target-gene composition, NOT a separate GO enrichment analysis.
category_order <- c("Replication + DDR", "Replication (no DDR)",
                    "Mitosis / segregation", "DNA damage / checkpoint",
                    "Other cell cycle", "Chromatin / transcription",
                    "RNA processing", "Other / unassigned")
# As in the OTC Fig. 2C visualization: NPG/NRC colors and dark-blue outlines.
# Reorder NPG colors to group related biology without changing category order.
npg <- ggsci::pal_npg("nrc")(10)
category_color <- stats::setNames(npg[c(4, 6, 1, 5, 9, 3, 2, 10)], category_order)
curve_col <- npg[3]   # Teal enrichment curve
frame_col <- npg[4]   # Dark-blue frames (GSEA, GO bar, leading-edge genes)
# These are deliberately broad Biological Process anchors. GOALL supplies
# descendant-to-ancestor membership; assignment is annotation-based, not
# an over-representation/enrichment p-value. Priority of non-overlap classes
# follows the order in this list, and is therefore a reporting choice.
go_roots <- list(
  "DNA replication" = c("GO:0006260"),
  "Mitosis / segregation" = c("GO:0007059", "GO:0007067"),
  "DNA damage / checkpoint" = c("GO:0006974", "GO:0000075"),
  "Other cell cycle" = c("GO:0007049"),
  "Chromatin / transcription" = c("GO:0006325", "GO:0006355"),
  "RNA processing" = c("GO:0006396")
)

# Audit all original GO-root matches first. If a gene matches replication
# AND DNA-damage/checkpoint roots, give it a dedicated overlapping display
# class. Otherwise retain the original priority for the remaining classes.
# "Replication (no DDR)" means no match to THESE DDR GO roots, not that the
# gene has no biological DNA-damage function.
assign_go_category <- function(gene_ids, skip_go) {
  out <- data.frame(ENSEMBL = unique(gene_ids), go_category = "Other / unassigned")
  audit <- data.frame(ENSEMBL = out$ENSEMBL, assigned_category = out$go_category,
    n_BP_GOALL_terms = 0L, n_matched_categories = 0L, matched_categories = "", matched_root_GO_IDs = "",
    replication_GO = FALSE, DDR_GO = FALSE)
  if (skip_go) {
    audit$n_BP_GOALL_terms <- NA_integer_; audit$n_matched_categories <- NA_integer_
    audit$matched_categories <- "not assessed"; audit$matched_root_GO_IDs <- "not assessed"
    audit$replication_GO <- NA; audit$DDR_GO <- NA
    attr(out, "audit") <- audit; return(out)
  }
  suppressMessages({
    # Retrieve all inherited GO terms, then restrict to Biological Process.
    # Multiple GO assignments per gene are normal and must be audited before
    # collapsing to mutually exclusive categories for stacked-bar displays.
    go <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = out$ENSEMBL,
      keytype = "ENSEMBL", columns = c("GOALL", "ONTOLOGYALL"))
  })
  go <- go[!is.na(go$GOALL) & go$ONTOLOGYALL == "BP", ]
  go_by_gene <- split(go$GOALL, go$ENSEMBL)
  for (i in seq_len(nrow(out))) {
    # Genes with no matching BP annotations remain "Other / unassigned";
    # this is absence of a match to selected roots, not absence of function.
    ids <- unique(go_by_gene[[out$ENSEMBL[i]]]); if (!length(ids)) next
    hits <- names(go_roots)[vapply(go_roots, function(root) any(root %in% ids), logical(1))]
    roots <- unique(unlist(lapply(hits, function(cat) intersect(go_roots[[cat]], ids)), use.names = FALSE))
    rep_hit <- "DNA replication" %in% hits; ddr_hit <- "DNA damage / checkpoint" %in% hits
    audit$n_BP_GOALL_terms[i] <- length(ids)
    audit$n_matched_categories[i] <- length(hits)
    audit$matched_categories[i] <- paste(hits, collapse = "; ")
    audit$matched_root_GO_IDs[i] <- paste(roots, collapse = "; ")
    audit$replication_GO[i] <- rep_hit; audit$DDR_GO[i] <- ddr_hit
    # Explicitly preserve replication + DNA-damage/checkpoint overlaps.
    # For all other multi-hit genes, the first listed GO theme wins only
    # for visualization; all original matches remain in the audit table.
    if (rep_hit && ddr_hit) {
      out$go_category[i] <- "Replication + DDR"
    } else if (rep_hit) {
      out$go_category[i] <- "Replication (no DDR)"
    } else if (length(hits)) {
      out$go_category[i] <- hits[1L]
    }
  }
  audit$assigned_category <- out$go_category
  attr(out, "audit") <- audit
  out
}

# Export the GO audit independently of the figure so each displayed class
# can be traced back to its GO-root matches. With --skip-go all genes receive
# a plotting placeholder ("Other / unassigned"), NOT a biological assignment.
go_categories <- assign_go_category(dream_ids, opts$skip_go)
go_audit <- attr(go_categories, "audit")
go_audit$publication_symbol <- targets$publication_symbol[match(go_audit$ENSEMBL, targets$ENSEMBL)]
go_audit <- go_audit[, c("ENSEMBL", "publication_symbol", setdiff(names(go_audit), c("ENSEMBL", "publication_symbol")))]
stopifnot(nrow(go_audit) == length(dream_ids), identical(go_audit$assigned_category, go_categories$go_category))
if (!opts$skip_go) stopifnot(sum(go_audit$replication_GO & go_audit$DDR_GO) ==
  sum(go_categories$go_category == "Replication + DDR"))
utils::write.csv(go_audit, file.path(out_dir, "DREAM_GO_category_audit.csv"), row.names = FALSE, na = "")
if (opts$skip_go) message("GO-BP annotation skipped; category audit contains no biological assignments.") else
  message("GO-BP audit: ", sum(go_audit$replication_GO & go_audit$DDR_GO),
    " targets match replication AND DDR; ", sum(go_audit$n_matched_categories > 1L),
    " match multiple GO groups. See DREAM_GO_category_audit.csv.")
dream_all <- merge(dream_all, go_categories, by = "ENSEMBL", all.x = TRUE, sort = FALSE)
dream_all$go_category <- factor(dream_all$go_category, levels = category_order)

# Collect algorithmic leading edges for all six contexts, irrespective of
# overall GSEA significance. The output table is descriptive; the figure
# intentionally hides gene/composition panels for nonsignificant GSEAs.
leading <- dream_all[dream_all$leading_edge, c(
  "ENSEMBL", "SYMBOL", "publication_symbol", "contrast", "cell_line",
  "radiation", "stat", "interaction_ashr", "padj",
  "WT_response_ashr", "KO_response_ashr", "WT_response_padj",
  "KO_response_padj", "repression_pattern", "regulatory_pattern", "go_category")]
leading <- leading[order(match(leading$contrast, contexts$contrast), -abs(leading$stat)), ]

# Describe class composition relative to the *rankable* DREAM genes in each contrast.
# These proportions are descriptive; they are not independent enrichment tests.
# Denominator 1: leading-edge genes for that context.
# Denominator 2: all rankable DREAM genes in the corresponding interaction.
# The two proportions help recognize representation, but no significance
# test is computed for GO themes here; do not label as GO enrichment.
go_composition <- do.call(rbind, lapply(contexts$contrast, function(key) {
  background_ids <- intersect(dream_ids, names(wald_ranks[[key]]))
  background <- go_categories$go_category[match(background_ids, go_categories$ENSEMBL)]
  n_background <- as.integer(table(factor(background, levels = category_order)))
  hits <- leading$go_category[leading$contrast == key]
  n_leading <- as.integer(table(factor(hits, levels = category_order)))
  data.frame(contrast = key, category = category_order, n_leading = n_leading,
    fraction_leading = n_leading / length(hits), n_rankable_DREAM = n_background,
    fraction_rankable_DREAM = n_background / length(background_ids))
}))
utils::write.csv(go_composition, file.path(out_dir, "DREAM_GO_category_composition.csv"), row.names = FALSE)

# ---- 05 | Figures: NES overview and OTC-style GSEA / leading-edge panels ----

# Use the project's ekbSeq plot-export convention for consistent raster and
# vector results. Figures are outputs of existing statistical summaries.
save_dual <- function(plot, file, width, height, dpi = 350) {
  ekbSeq::export_plot_dual(file.path(fig_dir, file), plot_expr = plot,
    width = width, height = height, dpi = dpi, rasterize = FALSE)
}

# Heatmap is a compact overview of signed GSEA outcomes. Tile annotations
# display BH-significant results from the SIX-test GSEA family, not the count
# of individually significant DREAM genes.
## The same exact blue–white–red colors/breaks as module Figures 13/14.
plot_nes <- function(x) {
  x$radiation <- factor(x$radiation, levels = c("7gy", "UVA", "UVB"),
                        labels = c("7 Gy", "UVA", "UVB"))
  x$cell_line <- factor(x$cell_line, levels = c("MCF10A", "HT1080"))
  lim <- max(abs(x$NES), na.rm = TRUE)
  ggplot2::ggplot(x, ggplot2::aes(x = radiation, y = cell_line, fill = NES)) +
    ggplot2::geom_tile(color = "white", linewidth = 2.0) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.2f%s", NES,
                                  ifelse(significant, " *", ""))),
                       size = 5.0) +
    ggplot2::scale_fill_gradient2(low = "#385c98", mid = "#fafafa",
                                  high = "#e34a33", midpoint = 0,
                                  limits = c(-lim, lim), name = "NES") +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(title = "DREAM-target enrichment in p21 interaction rankings",
                  subtitle = "KO-minus-WT interactions | original Wald statistics",
                  x = NULL, y = NULL,
                  caption = "* BH FDR < 0.05 across six DREAM tests; sign denotes ranking direction, not expression direction") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"))
}

save_dual(plot_nes(gsea), "DREAM_interaction_NES_heatmap", 10, 4.1)

# Build comparable triptychs for ALL six KO-minus-WT contrasts, including
# MCF10A UVA (not significant after correction). Figure membership is fixed;
# GO composition summarizes target functions, not an additional enrichment test.
# Build the canonical GSEA trace from the original signed Wald ranking.
# The enriched tail and boxed leading-edge region depend on NES direction;
# NES < 0 does not imply that targets were transcriptionally repressed in KO.
plot_gsea_panel <- function(key, summary_row, ranks) {
  gene_set <- intersect(dream_ids, names(ranks))
  p <- fgsea::plotEnrichment(gene_set, ranks)

  # fgsea hardcodes its enrichment-curve color. Change only the line layer;
  # retain the original running scores, gene-set tick marks and ES statistics.
  for (j in seq_along(p$layers))
    if (inherits(p$layers[[j]]$geom, "GeomLine")) p$layers[[j]]$aes_params$colour <- curve_col

  # OTC-style box over the contributing leading-edge part of the ranking:
  # positive NES -> left/above zero; negative NES -> right/below zero.
  # Do not emphasize an algorithmic leading edge for nonsignificant GSEAs.
  if (isTRUE(summary_row$significant)) {
    le_pos <- match(results[[key]]$primary$leading_edge, names(ranks))
    le_pos <- le_pos[is.finite(le_pos)]
    if (!length(le_pos)) stop("No leading-edge positions for: ", key)
    positive <- summary_row$NES > 0
    p <- p + ggplot2::annotate("rect", xmin = if (positive) 1 else min(le_pos),
      xmax = if (positive) max(le_pos) else length(ranks),
      ymin = if (positive) 0 else -Inf, ymax = if (positive) Inf else 0,
      fill = NA, color = frame_col, linewidth = 0.85)
  }
  p + ggplot2::labs(title = sub("7gy", "7 Gy", sub("_", " | ", key), fixed = TRUE),
      subtitle = sprintf("NES %+.2f | BH FDR %.2g%s", summary_row$NES, summary_row$FDR_six,
        ifelse(summary_row$significant, "", " (n.s.)")),
      x = "Original KO-minus-WT interaction Wald ranking", y = "Enrichment score") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

# Compact, label-only stacked composition bars. Replication/DDR double-annotated
# genes are explicit; each gene still contributes exactly once per bar.
# Category fractions are of the leading edge, not all DREAM genes and not
# all tested genes. Each gene contributes to exactly one display category.
plot_composition <- function(x, legend_rows = 1L) {
  counts <- data.frame(category = factor(category_order, levels = category_order),
    n = as.integer(table(factor(x$go_category, levels = category_order))))
  counts$percent <- counts$n / sum(counts$n)
  ggplot2::ggplot(counts, ggplot2::aes(x = 1, y = percent, fill = category)) +
    ggplot2::geom_col(width = 0.62, color = "white", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(percent >= 0.08,
      paste0(round(percent * 100), "%"), "")), position = ggplot2::position_stack(vjust = 0.5),
      size = 3, fontface = "bold", color = "#172536") +
    ggplot2::annotate("rect", xmin = 0.69, xmax = 1.31, ymin = 0, ymax = 1,
      fill = NA, color = frame_col, linewidth = 0.85) +
    ggplot2::scale_fill_manual(values = category_color, breaks = category_order, drop = FALSE,
      name = "GO-BP category", guide = ggplot2::guide_legend(nrow = legend_rows, byrow = TRUE)) +
    ggplot2::scale_x_continuous(limits = c(0.58, 1.42), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1.00001), expand = c(0, 0)) +
    ggplot2::labs(title = "GO-BP composition", subtitle = paste0("n = ", nrow(x), " leading-edge genes")) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
      plot.subtitle = ggplot2::element_text(size = 8), plot.margin = ggplot2::margin(3, 3, 3, 3))
}

# Select the 25 strongest |Wald| contributors first, then group by the
# prespecified GO-BP priority; individual categories never drive selection.
# Top 25 by absolute original Wald statistic within a *selected* leading
# edge, then reorder visually by GO theme. Ordering by GO does not influence
# which genes are chosen and must not be read as a significance hierarchy.
prepare_top25 <- function(x) {
  x <- utils::head(x[order(-abs(x$stat)), , drop = FALSE], 25L)
  x <- x[order(match(as.character(x$go_category), category_order), -abs(x$stat)), , drop = FALSE]
  x$plot_label <- ifelse(!is.na(x$SYMBOL) & nzchar(x$SYMBOL), x$SYMBOL, x$publication_symbol)
  x$plot_pos <- seq_len(nrow(x))
  x$pos <- ((x$plot_pos - 1L) %% 13L) + 1L
  x$column <- ((x$plot_pos - 1L) %/% 13L) + 1L
  x
}

# OTC-style gene panel: framed list, two columns and solid colored squares.
# The title is inside the same dark-blue frame as the 25 gene labels.
plot_leading_genes <- function(x) {
  x$square_x <- ifelse(x$column == 1L, 0.95, 2.01)
  x$text_x <- x$square_x + 0.10
  ggplot2::ggplot(x) +
    ggplot2::geom_point(ggplot2::aes(x = square_x, y = -pos, color = go_category),
      shape = 15, size = 2.65, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(x = text_x, y = -pos, label = plot_label),
      hjust = 0, size = 3, show.legend = FALSE) +
    ggplot2::annotate("rect", xmin = 0.82, xmax = 3.08, ymin = -13.64, ymax = 0.65,
      fill = NA, color = frame_col, linewidth = 0.85) +
    ggplot2::annotate("text", x = 0.93, y = 0.12, label = "Leading-edge genes (top 25)",
      hjust = 0, fontface = "bold", size = 3.35) +
    ggplot2::scale_color_manual(values = category_color, breaks = category_order, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0.78, 3.12), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(-13.79, 0.82), expand = c(0, 0)) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(plot.margin = ggplot2::margin(3, 4, 3, 4))
}

# A nonsignificant GSEA retains an algorithmic leading edge, but we deliberately
# do not display or interpret its composition as an enriched gene set.
plot_no_leading <- function(title) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5,
      label = "No leading edge shown\n(not significantly enriched)",
      size = 3.15, color = "grey45", lineheight = 1.25) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(panel.border = ggplot2::element_rect(color = frame_col, fill = NA, linewidth = 0.85),
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      plot.margin = ggplot2::margin(3, 4, 3, 4))
}

# HT1080: 7 Gy, UVA, UVB (top); MCF10A: 7 Gy, UVA, UVB (bottom).
# Flat 2 x 9 patchwork grid collects GO fills into ONE shared bottom legend;
# the same triptychs are also exported separately for closer inspection.
# GSEA curves remain visible for all six predefined experiments.
# GO and gene panels are shown only if that context's GSEA passes BH FDR.
# This keeps exploratory leading-edge composition distinct from evidence
# for statistically significant enrichment.
significant_keys <- gsea$contrast[gsea$significant]
all_triptychs <- setNames(vector("list", nrow(contexts)), contexts$contrast)
six_panels <- list(); plotted_genes <- vector("list", nrow(contexts))
for (i in seq_len(nrow(contexts))) {
  key <- contexts$contrast[i]
  x <- leading[leading$contrast == key, , drop = FALSE]
  significant <- isTRUE(gsea$significant[gsea$contrast == key])
  a <- plot_gsea_panel(key, gsea[gsea$contrast == key, ], wald_ranks[[key]])
  if (significant) {
    top25 <- prepare_top25(x)
    plotted_genes[[i]] <- top25[, c("contrast", "ENSEMBL", "plot_label", "go_category", "stat", "plot_pos")]
    b <- plot_composition(x, legend_rows = 1L); c <- plot_leading_genes(top25)
    b_single <- plot_composition(x, legend_rows = 2L)
  } else {
    # Preserve the GSEA curve and its nonsignificant FDR, but omit the GO/gene lists.
    b <- plot_no_leading("GO-BP composition"); c <- plot_no_leading("Leading-edge genes")
    b_single <- b
  }
  # The standalone triptych collects the GO legend beneath ALL three panels;
  # this prevents it from encroaching on the bar or the gene list.
  all_triptychs[[key]] <- patchwork::wrap_plots(a, b_single, c, ncol = 3,
    widths = c(2.7, 1, 1.8), guides = "collect") &
    ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal",
      legend.title = ggplot2::element_text(size = 8), legend.text = ggplot2::element_text(size = 8),
      legend.key.size = grid::unit(3, "mm"), legend.spacing.x = grid::unit(1.5, "mm"))
  save_dual(all_triptychs[[key]], paste0("DREAM_GSEA_", key, "_leading_edge"), 17, 5.4)
  six_panels <- c(six_panels, list(a, b, c))
}
utils::write.csv(do.call(rbind, plotted_genes), file.path(out_dir, "DREAM_top25_plotted.csv"), row.names = FALSE)

# Two shorter rows of nine component plots = three triptychs per row. Unlike
# nested patchworks, the flat grid reliably collects identical fill legends.
p_six <- patchwork::wrap_plots(six_panels, ncol = 9, byrow = TRUE,
    widths = rep(c(2.7, 1, 1.8), 3), heights = c(1, 1), guides = "collect") +
  patchwork::plot_annotation(
    title = "DREAM-target enrichment in p21-dependent irradiation responses",
    subtitle = "HT1080 (top) and MCF10A (bottom) | 7 Gy, UVA, UVB | KO-minus-WT interactions",
    caption = "GSEA: original Wald ranking; BH FDR across six tests. GO-BP target functions: DDR = DNA damage response/checkpoint; 25 genes selected by |Wald|, then GO-grouped.")
p_six <- p_six & ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal",
  legend.text = ggplot2::element_text(size = 9), legend.title = ggplot2::element_text(size = 9),
  legend.key.size = grid::unit(3.5, "mm"))
save_dual(p_six, "DREAM_GSEA_six_contexts", width = 32, height = 11.7, dpi = 220)

# Membership matrix for direct comparison of the leading edges; includes the
# descriptive MCF10A UVA set despite its nonsignificant GSEA result.
# Presence matrix for overlap inspection: presence counts are not tests of
# agreement/replicability across irradiation qualities or cell lines.
# The matrix intentionally includes a nonsignificant context; use FDR_six
# when interpreting its context-specific leading-edge entries.
union_ids <- unique(leading$ENSEMBL)
le_presence <- data.frame(ENSEMBL = union_ids,
  SYMBOL = targets$publication_symbol[match(union_ids, targets$ENSEMBL)],
  go_category = as.character(go_categories$go_category[match(union_ids, go_categories$ENSEMBL)]))
for (key in contexts$contrast)
  le_presence[[key]] <- as.integer(union_ids %in% leading$ENSEMBL[leading$contrast == key])
le_presence$n_contexts <- rowSums(le_presence[, contexts$contrast, drop = FALSE])
le_presence <- le_presence[order(-le_presence$n_contexts,
  match(le_presence$go_category, category_order), le_presence$SYMBOL), ]
utils::write.csv(le_presence, file.path(out_dir, "DREAM_leading_edge_presence.csv"), row.names = FALSE)

# ---- 06 | Preservation / attenuation of irradiation responses --------------

# Per-gene x/y coordinates are shrunken WT/KO irradiation effects.
# The diagonal is equal response (KO = WT), while point COLOR is the formal
# KO-minus-WT interaction FDR test. This avoids inferring a difference merely
# because one genotype has a significant response and the other does not.
plot_response_scatter <- function(x) {
  x$radiation <- factor(x$radiation, levels = c("7gy", "UVA", "UVB"),
                        labels = c("7 Gy", "UVA", "UVB"))
  x$interaction_direction <- factor(x$interaction_direction,
    levels = c("not_significant", "positive_KO_minus_WT", "negative_KO_minus_WT"))
  ggplot2::ggplot(x, ggplot2::aes(x = WT_response_ashr, y = KO_response_ashr)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey85") +
    ggplot2::geom_vline(xintercept = 0, color = "grey85") +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey50",
                         linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(color = interaction_direction),
                        size = 1.75, alpha = .72, na.rm = TRUE) +
    ggplot2::facet_grid(cell_line ~ radiation) +
    ggplot2::coord_equal() +
    ggplot2::scale_color_manual(values = c("not_significant" = "#b8bec5",
      "positive_KO_minus_WT" = "#e34a33", "negative_KO_minus_WT" = "#385c98"),
      labels = c("Not significant", "Positive Δp21", "Negative Δp21")) +
    ggplot2::labs(title = "DREAM-target irradiation responses in WT and p21 KO",
      subtitle = "Dashed line: equal WT and KO response; colors: formal interaction FDR < 0.05",
      x = "WT irradiation log2FC (ashr)", y = "p21 KO irradiation log2FC (ashr)",
      color = "p21 interaction",
      caption = "Significant WT/KO irradiation responses and formal interactions are distinct tests.") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
}

save_dual(plot_response_scatter(dream_all), "DREAM_WT_KO_response_scatter", 12, 7)

# Heatmap focuses on genes from significant interactions, but never replaces
# full-population descriptive summaries or gene-level statistical results.
# Heatmap gene selection is intentionally biased toward strong/recurrent
# leading-edge signals from GSEA-significant contexts. Rank genes first by
# number of contexts, then by maximum absolute *shrunken interaction* effect;
# cap at 50 for legibility. The heatmap illustrates candidate patterns,
# not a representative genome-wide view or independent validation.
if (length(significant_keys)) {
  sig_le <- leading[leading$contrast %in% significant_keys, , drop = FALSE]
  counts <- table(sig_le$ENSEMBL)
  ids <- names(counts)
  importance <- vapply(ids, function(id) {
    max(abs(sig_le$interaction_ashr[sig_le$ENSEMBL == id]), na.rm = TRUE)
  }, numeric(1))
  ids <- ids[order(-as.integer(counts[ids]), -importance)]
  ids <- head(ids, 50L)
  hm <- dream_all[dream_all$ENSEMBL %in% ids,
                  c("ENSEMBL", "SYMBOL", "publication_symbol", "cell_line",
                    "radiation", "WT_response_ashr", "KO_response_ashr")]
  hm <- rbind(
    data.frame(hm[c("ENSEMBL", "SYMBOL", "publication_symbol", "cell_line", "radiation")],
               genotype = "WT", effect = hm$WT_response_ashr),
    data.frame(hm[c("ENSEMBL", "SYMBOL", "publication_symbol", "cell_line", "radiation")],
               genotype = "KO", effect = hm$KO_response_ashr))
  hm$label <- ifelse(!is.na(hm$SYMBOL) & nzchar(hm$SYMBOL),
                     hm$SYMBOL, hm$publication_symbol)
  lab <- hm[match(ids, hm$ENSEMBL), c("ENSEMBL", "label")]
  lab$label <- make.unique(lab$label)
  hm$label <- lab$label[match(hm$ENSEMBL, lab$ENSEMBL)]
  hm$label <- factor(hm$label, levels = rev(lab$label))
  hm$condition <- paste(hm$radiation, hm$genotype, sep = "_")
  hm$condition <- factor(hm$condition,
    levels = c("7gy_WT", "7gy_KO", "UVA_WT", "UVA_KO", "UVB_WT", "UVB_KO"))
  h_lim <- max(1, ceiling(max(abs(hm$effect), na.rm = TRUE) * 2) / 2)
  p <- ggplot2::ggplot(hm, ggplot2::aes(x = condition, y = label, fill = effect)) +
    ggplot2::geom_tile(color = "white", linewidth = .15) +
    ggplot2::facet_grid(. ~ cell_line, scales = "free_x") +
    ggplot2::scale_fill_gradient2(low = "#385c98", mid = "white",
      high = "#e34a33", limits = c(-h_lim, h_lim), name = "ashr log2FC") +
    ggplot2::labs(title = "Leading-edge DREAM genes | irradiation responses",
      subtitle = "Union of leading edges from significant DREAM GSEA contrasts; up to 50 genes",
      x = NULL, y = NULL, caption = "Selected genes are not representative of all DREAM targets.") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_dual(p, "DREAM_leading_edge_response_heatmap", 11.5,
            max(6, min(12, 2.4 + length(ids) * .17)))
}

# ---- 07 | Auditable Excel and complete CSV exports --------------------------

# Record paths to the source results directories (not temporary transport
# archives), alongside analysis settings and software versions. Individual
# CSV filenames are fixed in the input loop above.
meta <- data.frame(Item = c(
  "Analysis", "Study", "Input targets", "GSEA ranking", "Gene set",
  "Primary GSEA universe", "Primary GSEA test family", "GSEA threshold",
  "Sensitivity universe", "Response effects", "Interaction definition",
  "Interaction FDR", "GO categories", "Gene selection for heatmap",
  "Targets file", "Interaction results directory", "Response results directory",
  "Targets MD5",
  "GO annotation database", "GO category assignment priority",
  "Six-context GSEA figure", "Top-25 gene display", "R version", "ekbSeq version", "fgsea version", "Date"),
 Value = c(
  "09_DREAM_analysis.R", "p21 irradiation responses | HT1080 & MCF10A",
  "Engeland et al. 2022, PMID 35361964", "Original, unshrunken DESeq2 Wald statistic",
  paste0(nrow(targets), " DREAM=1 publication targets"),
  "All unique genes with finite original Wald statistic in each interaction",
  "One DREAM set x six contrasts, BH across all six p-values", "BH FDR < 0.05",
  "Genes with finite Wald statistic and baseMean >= 10 (secondary)",
  "Original ashr-shrunken irradiation vs control log2FC",
  "(KOirr-KOctrl)-(WTirr-WTctrl) (processing adjusted)",
  "Original DESeq2 BH padj from six individual interactions",
  if (opts$skip_go) "Disabled (--skip-go)" else
    "GOALL BP ancestry; replication/DDR overlaps shown explicitly, remaining classes prioritized; NOT enrichment",
  "Leading-edge union of BH-significant GSEAs only; top 50 by recurrence and |ashr|",
  normalizePath(targets_file, winslash = "/"),
  normalizePath(interaction_dir, winslash = "/"),
  normalizePath(response_dir, winslash = "/"),
  unname(tools::md5sum(targets_file)),
  if (opts$skip_go) "Not used (--skip-go)" else paste0("org.Hs.eg.db ", utils::packageVersion("org.Hs.eg.db")),
  paste0("Replication AND DDR -> Replication + DDR; replication without DDR -> Replication (no DDR); ",
    "other GO groups prioritized: ", paste(names(go_roots)[-1], collapse = " > "), ". DDR roots: GO:0006974, GO:0000075"),
  "All six GSEA curves; NPG-colored GO composition and top 25 leading genes shown only for BH-significant contexts",
  "In BH-significant contexts: select 25 leading genes by largest |Wald statistic|, then group by GO-BP priority",
  R.version.string, as.character(utils::packageVersion("ekbSeq")),
  as.character(utils::packageVersion("fgsea")), as.character(Sys.Date())))

# Export a flat, interpretable gene-level table (raw vs shrunken remain distinct).
cols <- c("cell_line", "radiation", "contrast", "ENSEMBL", "SYMBOL",
 "publication_symbol", "baseMean", "stat", "log2FoldChange", "lfcSE",
 "interaction_ashr", "pvalue", "padj", "interaction_direction",
 "WT_response_ashr", "KO_response_ashr", "WT_response_padj",
 "KO_response_padj", "WT_repressed", "KO_repressed", "WT_induced", "KO_induced",
 "repression_pattern", "regulatory_pattern", "leading_edge", "go_category")
dream_export <- dream_all[, cols]
names(dream_export)[names(dream_export) == "log2FoldChange"] <- "interaction_Wald_log2FC"
names(dream_export)[names(dream_export) == "lfcSE"] <- "interaction_Wald_SE"

# Complete machine-readable tables are the primary analysis deliverables.
# Publication-specific gene selections should reference these full tables,
# rather than silently redefining thresholds when assembling a manuscript.
utils::write.csv(dream_export, file.path(out_dir, "DREAM_gene_results_all.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(gsea, file.path(out_dir, "DREAM_GSEA_summary.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(leading, file.path(out_dir, "DREAM_leading_edge_genes.csv"),
                 row.names = FALSE, na = "")
utils::write.csv(patterns, file.path(out_dir, "DREAM_response_patterns.csv"),
                 row.names = FALSE, na = "")

# WORKFLOW NOTE: the existing Excel export is preserved for compatibility.
# Current Research Analysis Principles favor assembling manuscript/supplement
# workbooks later in a dedicated export script. Revising that workflow is
# separate from correcting the input paths in this script.
sheets <- list("GSEA_summary" = gsea,
               "Response_patterns" = patterns,
               "DREAM_gene_results" = dream_export,
               "Leading_edges" = leading,
               "Leading_edge_GO" = go_categories,
               "Low_count_QC" = low_qc,
               "Metadata" = meta)
ekbSeq::export_excel(sheets, file.path(out_dir, "DREAM_analysis_results.xlsx"))

message("\nDREAM analysis complete. Outputs: ", normalizePath(out_dir),
        "\nSignificant GSEA: ", paste(significant_keys, collapse = ", "))
print(gsea[, c("contrast", "NES", "pvalue", "FDR_six",
               "NES_baseMean10", "FDR_six_baseMean10")], row.names = FALSE)