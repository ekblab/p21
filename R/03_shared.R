## p21 / 03 | Minimal shared project paths and plotting conventions.
## Always sourced by one of the three stages. No expensive calculations.
if (!requireNamespace("here", quietly = TRUE)) stop("Missing package: here")
here::i_am("R/03_response_architecture.R")
p21_result_root <- here::here("Results", "03_response_architecture")
p21_module_dir <- file.path(p21_result_root, "02_MODULE_L2")
p21_evaluation_dir <- file.path(p21_module_dir, "k_evaluation")
p21_pca_dir <- file.path(p21_result_root, "01_PCA")
p21_sensitivity_dir <- file.path(p21_result_root, "03_PROCESSING_MATCHED_SENSITIVITY")
p21_module_order <- c("M01", "M08", "M02", "M06", "M03", "M05", "M04", "M07")
p21_module_colors <- c(M01 = "#E64B35", M02 = "#4DBBD5", M03 = "#00A087",
  M04 = "#3C5488", M05 = "#F39B7F", M06 = "#8491B4", M07 = "#91D1C2", M08 = "#7E6148")
p21_check <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required R packages: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}
p21_export <- function(path, plot, width, height) {
  ekbSeq::export_plot_dual(path, plot, width = width, height = height, rasterize = FALSE)
}
