## ************************************************************************************************
## p21 bulk RNA-seq - Shared helpers
##
## Only reusable operations required across the planned scripts belong here.
## Primary and downstream results are read/written as CSV or XLSX, NOT as RDS.
## The only persistent RDS is Data/processed_data.rds.
## ************************************************************************************************

## ******
## validate the model matrix before fitting a DESeq2 model
check_design <- function(metadata, formula) {
  mm <- stats::model.matrix(formula, data = metadata)
  if (qr(mm)$rank < ncol(mm)) stop("Rank-deficient model: ", deparse(formula))
  message(deparse(formula), ": ", nrow(mm), " samples, ", ncol(mm), " coefficients (full rank).")
  invisible(mm)
}

## ******
## construct canonical metadata group names (24-h experiment)
group_name <- function(cell_line, genotype, radiation) {
  state <- ifelse(radiation == "control", "control", paste0("irradiated_", radiation))
  paste(cell_line, genotype, state, "24h", sep = "_")
}

## ******
## signed group weights for treatment vs matched control
response_weights <- function(cell_line, genotype, radiation) {
  if (radiation == "control") stop("A response requires an irradiated condition.")
  stats::setNames(c(1, -1), c(group_name(cell_line, genotype, radiation),
                              group_name(cell_line, genotype, "control")))
}

## p21-negative minus WT irradiation response (difference-in-differences)
interaction_weights <- function(cell_line, radiation) {
  c(response_weights(cell_line, "p21_neg", radiation), -response_weights(cell_line, "wt", radiation))
}

## ******
## convert explicitly named group weights into a DESeq2 numeric contrast
## Nuisance covariates are held at the same level in all artificial design rows.
## Hence a treatment contrast is adjusted for processing, not an observed-profile contrast.
group_contrast <- function(dds, weights) {
  if (!length(weights) || is.null(names(weights)) || anyNA(names(weights)))
    stop("Named group weights are required.")
  if (anyNA(weights) || any(!is.finite(weights))) stop("Invalid group weights.")

  meta <- as.data.frame(SummarizedExperiment::colData(dds))
  if (!all(names(weights) %in% levels(meta$group))) stop("Unknown group in contrast.")
  collapsed <- tapply(as.numeric(weights), names(weights), sum)
  weights <- stats::setNames(as.numeric(collapsed), dimnames(collapsed)[[1]])
  if (!isTRUE(all.equal(sum(weights), 0, tolerance = 1e-9))) stop("Group weights must sum to zero.")

  reference_rows <- meta[rep(1L, length(weights)), , drop = FALSE]
  reference_rows$group <- factor(names(weights), levels = levels(meta$group))
  reference_rows$processing <- factor(levels(meta$processing)[1], levels = levels(meta$processing))
  if ("sequencing_run" %in% all.vars(DESeq2::design(dds)))
    reference_rows$sequencing_run <- factor(levels(meta$sequencing_run)[1], levels = levels(meta$sequencing_run))

  design_matrix <- stats::model.matrix(DESeq2::design(dds), data = meta)
  reference_matrix <- stats::model.matrix(DESeq2::design(dds), data = reference_rows)
  if (!identical(colnames(design_matrix), colnames(reference_matrix))) stop("Design column mismatch.")
  if (qr(design_matrix)$rank < ncol(design_matrix)) stop("Rank-deficient DESeq2 design.")

  contrast <- as.numeric(crossprod(weights, reference_matrix))
  if (length(contrast) != length(DESeq2::resultsNames(dds))) stop("DESeq2 coefficient mismatch.")
  names(contrast) <- DESeq2::resultsNames(dds)
  if (any(abs(contrast[grepl("^(processing|sequencing_run)", names(contrast))]) > 1e-9))
    stop("Nuisance effects did not cancel in the contrast.")
  contrast
}

## ******
## merge a DESeq2 results object with the one-row-per-feature annotation
result_table <- function(res, annotation) {
  tab <- data.frame(ENSEMBL = rownames(res), as.data.frame(res), check.names = FALSE, row.names = NULL)
  if (anyDuplicated(annotation$ENSEMBL)) stop("Duplicate annotation feature IDs.")
  idx <- match(tab$ENSEMBL, annotation$ENSEMBL)
  if (anyNA(idx)) stop("DESeq2 results contain features without annotation rows.")
  cbind(tab, annotation[idx, setdiff(names(annotation), "ENSEMBL"), drop = FALSE])
}
