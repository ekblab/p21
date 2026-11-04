## ************************************************************************************************
## p21 bulk RNA-seq
## 07 - IPA multi-observation export (24 h)
##
## Goal:
##   Create ONE Ingenuity Pathway Analysis (IPA) upload file with 18 observations:
##     - 12 irradiation-vs-control responses (2 cell lines x WT/p21-KO x 3 irradiations)
##     -  6 formal p21 interactions (KO response - WT response)
##
## Input: six internally consistent ashr-shrunken interaction CSVs from script 04.
##   These already contain the separately shrunken WT/KO responses from script 02,
##   their respective DESeq2 BH-adjusted p-values, AND the directly shrunken
##   p21 interaction with its own BH-adjusted p-value.
##
## No refitting, no additional tests, no LFC/significance filter, no extra outputs.
## The Wald tests were performed on unshrunken estimates in scripts 02/04;
## log2FC shrinkage does not change the adjusted p-values.
##
## Output: Results/07_IPA/p21_IPA_18_observations.txt (tab-delimited)
##   ID | <observation> Expr Log Ratio | <observation> Expr FDR | ...
##
## In IPA: map the Ensembl ID, log2FC columns to 'Expr Log Ratio', and adjusted
## p-values to 'Expr FDR'. Interactions are KO-minus-WT differences of irradiation
## effects; their signs are NOT absolute expression levels or activation scores.
##
## Run after scripts 02 and 04 whenever the upstream analysis changes.
## ************************************************************************************************

## ******
## Dependencies and paths
if (!requireNamespace("here", quietly = TRUE))
  stop("Missing R package 'here'. Install it with install.packages('here').")
here::i_am("R/07_IPA_export.R")

input_dir <- here::here("Results", "04_p21_interactions", "tables", "shrunken")
output_dir <- here::here("Results", "07_IPA")
output_file <- file.path(output_dir, "p21_IPA_18_observations.txt")

cell_lines <- c("HT1080", "MCF10A")
radiations <- c("7gy", "UVA", "UVB")
radiation_label <- function(x) if (x == "7gy") "7GY" else x

## ******
## Define all 18 observations in a fixed, biologically readable column order.
response_labels <- unlist(lapply(cell_lines, function(cell_line) {
  unlist(lapply(c("WT", "p21KO"), function(genotype) {
    vapply(radiations, function(rad)
      paste(cell_line, genotype, radiation_label(rad), sep = "_"), character(1))
  }), use.names = FALSE)
}), use.names = FALSE)

interaction_labels <- unlist(lapply(cell_lines, function(cell_line) {
  vapply(radiations, function(rad)
    paste(cell_line, "p21KOminusWT", radiation_label(rad), sep = "_"), character(1))
}), use.names = FALSE)

labels <- c(response_labels, interaction_labels)
stopifnot(length(labels) == 18L, !anyDuplicated(labels))

source_paths <- unlist(lapply(cell_lines, function(cell_line)
  vapply(radiations, function(rad)
    file.path(input_dir, paste0(cell_line, "_", rad, "_p21_interaction_ashr.csv")),
    character(1))), use.names = FALSE)

missing_files <- source_paths[!file.exists(source_paths)]
if (length(missing_files)) {
  stop("Missing ashr interaction results from script 04:\n  ",
       paste(missing_files, collapse = "\n  "),
       "\nRe-run R/04_p21_interactions.R before exporting to IPA.")
}

## ******
## Convert CSV values to valid IPA cells without rounding the original numbers.
## NAs become empty cells; finite numerical values retain their original text.
ipa_numeric <- function(x, context, fdr = FALSE) {
  x <- trimws(x)
  x[is.na(x) | x %in% c("NA", "NaN", "NULL", "")] <- ""
  present <- nzchar(x)
  if (any(present)) {
    parsed <- suppressWarnings(as.numeric(x[present]))
    if (any(!is.finite(parsed)))
      stop("Non-numeric or non-finite value in ", context, ".")
    if (fdr && any(parsed < 0 | parsed > 1))
      stop("Adjusted p-values outside [0, 1] in ", context, ".")
  }
  x
}

required <- c("ENSEMBL", "log2FoldChange", "padj", "lfc_method",
              "WT_response_log2FC", "WT_response_padj",
              "p21neg_response_log2FC", "p21neg_response_padj")

genes <- NULL
values <- NULL
filled <- rep(FALSE, length(labels))

## ************************************************************************************************
## Read the six interaction tables; each provides three observations.
for (cell_line in cell_lines) {
  for (rad in radiations) {
    path <- file.path(input_dir,
                      paste0(cell_line, "_", rad, "_p21_interaction_ashr.csv"))
    message("Reading: ", basename(path))
    tab <- utils::read.csv(path, colClasses = "character", check.names = FALSE,
                           stringsAsFactors = FALSE, na.strings = c("", "NA"))

    if (!all(required %in% names(tab)))
      stop("Unexpected columns in ", path, "; missing: ",
           paste(setdiff(required, names(tab)), collapse = ", "))
    if (!nrow(tab) || anyNA(tab$ENSEMBL) || any(!nzchar(tab$ENSEMBL)) ||
        anyDuplicated(tab$ENSEMBL))
      stop("Missing or duplicate Ensembl gene IDs in ", path)
    if (any(is.na(tab$lfc_method)) ||
        any(tab$lfc_method != "ashr_shrunken"))
      stop("Input is not consistently labeled 'ashr_shrunken': ", path)

    ## All six tables must contain precisely the same gene IDs; row order may differ.
    if (is.null(genes)) {
      genes <- sort(tab$ENSEMBL)
      values <- matrix("", nrow = length(genes), ncol = 2L * length(labels))
    } else if (nrow(tab) != length(genes) || !setequal(tab$ENSEMBL, genes)) {
      stop("Gene universe differs between script-04 exports: ", path,
           "\nRegenerate all six interaction tables from the same analysis.")
    }
    tab <- tab[match(genes, tab$ENSEMBL), , drop = FALSE]

    rad_label <- radiation_label(rad)
    fields <- list(
      list(label = paste(cell_line, "WT", rad_label, sep = "_"),
           lfc = "WT_response_log2FC", fdr = "WT_response_padj"),
      list(label = paste(cell_line, "p21KO", rad_label, sep = "_"),
           lfc = "p21neg_response_log2FC", fdr = "p21neg_response_padj"),
      list(label = paste(cell_line, "p21KOminusWT", rad_label, sep = "_"),
           lfc = "log2FoldChange", fdr = "padj")
    )

    for (spec in fields) {
      position <- match(spec$label, labels)
      if (is.na(position) || filled[position])
        stop("Duplicate or unknown observation: ", spec$label)
      fc <- ipa_numeric(tab[[spec$lfc]], paste(spec$label, spec$lfc))
      q  <- ipa_numeric(tab[[spec$fdr]], paste(spec$label, spec$fdr), fdr = TRUE)
      if (any(!nzchar(fc) & nzchar(q)))
        stop("Adjusted p-value without a log2FC in ", spec$label)

      values[, 2L * position - 1L] <- fc
      values[, 2L * position] <- q
      filled[position] <- TRUE
      message(sprintf("  %-29s  FDR < 0.05: %d",
                      spec$label, sum(nzchar(q) & as.numeric(replace(q, !nzchar(q), "1")) < 0.05)))
    }
  }
}

stopifnot(all(filled), length(genes) > 0L, ncol(values) == 36L,
          !anyDuplicated(genes))

## ************************************************************************************************
## Write the single tab-separated IPA file.
headers <- c("ID", as.vector(rbind(paste0(labels, " Expr Log Ratio"),
                                    paste0(labels, " Expr FDR"))))
export <- data.frame(ID = genes, values, check.names = FALSE,
                     stringsAsFactors = FALSE)
names(export) <- headers

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(export, file = output_file, sep = "\t", quote = FALSE,
                   row.names = FALSE, col.names = TRUE, na = "")
message(sprintf("\nIPA export complete: %s\n  %d genes | %d observations | %d columns",
                output_file, nrow(export), length(labels), ncol(export)))
message("Upload as Expr Log Ratio + Expr FDR (BH-adjusted p-values).")
