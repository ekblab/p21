# Post-run file and schema checks; run from the project root.
source("R/project_helpers.R")
tags <- expand.grid(line = c("HT1080", "MCF10A"), genotype = c("wt", "p21_neg"),
  radiation = c("7gy", "UVA", "UVB"), stringsAsFactors = FALSE)
expected <- character()
for (i in seq_len(nrow(tags))) {
  line <- tags$line[i]; genotype <- tags$genotype[i]; radiation <- tags$radiation[i]
  trt <- group_name(line, genotype, radiation); ctrl <- group_name(line, genotype, "control")
  expected <- c(expected, response_path(line, genotype, radiation),
    file.path("Results/01_primary_DE/tables", paste0(trt, "_vs_", ctrl, "_unshrunken.xlsx")),
    file.path("Results/01_primary_DE/tables", paste0(trt, "_vs_", ctrl, "_ashr.xlsx")))
}
for (line in c("HT1080", "MCF10A")) for (radiation in c("7gy", "UVA", "UVB"))
  expected <- c(expected, interaction_path(line, radiation),
    file.path("Results/02_p21_interactions/response_scatter", paste0(line, "_", radiation, ".png")))
for (line in c("HT1080", "MCF10A")) for (pair in c("UVB_vs_UVA", "IR_vs_UVA", "IR_vs_UVB"))
  expected <- c(expected, file.path("Results/04_radiation_quality/tables", paste0(line, "_", pair, "_p21_difference.xlsx")))
for (radiation in c("7gy", "UVA", "UVB"))
  expected <- c(expected, file.path("Results/07_cellline_context/tables", paste0(radiation, "_cellline_difference_in_p21.xlsx")))
expected <- c(expected, "Results/01_primary_DE/historical_audit.csv",
  "Results/03_response_architecture/response_matrix.rds",
  "Results/05_functional/GSEA/hallmark_NES_matrix.png",
  "R/p21_QC_report.html", "R/p21_analysis_report.html")
missing <- expected[!file.exists(expected)]
if (length(missing)) stop("Missing outputs:\n", paste(missing, collapse = "\n"))
res <- read_result(response_path("HT1080", "wt", "UVB"))
stopifnot(all(c("unshrunken", "shrunken") %in% names(res)),
  identical(res$unshrunken$ENSEMBL, res$shrunken$ENSEMBL),
  all(c("log2FoldChange", "lfcSE", "stat", "pvalue", "padj") %in% names(res$unshrunken)))
cat("PASS:", length(expected), "expected files and representative result schemas.\n")
