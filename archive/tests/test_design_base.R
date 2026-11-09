# Base-R-only smoke tests; run from the project root, without DESeq2.
source("R/project_helpers.R")
meta <- read.csv("Data/metadata_v0.1_20250514.csv", check.names = FALSE)
counts <- read.delim("Data/p21_counts_counting_summary_summary_unique.tsv", check.names = FALSE)
stopifnot(nrow(meta) == 96L, ncol(counts) - 2L == 96L, nrow(counts) == 61920L,
  setequal(names(counts)[-(1:2)], as.character(meta$ID)), !anyDuplicated(counts$gene_id))
exclude <- meta$processing == "UKD" & meta$cell_type == "MCF10A_wt" & meta$irradiation == "control"
stopifnot(sum(exclude) == 5L)
meta <- meta[!exclude, ]; meta$group <- factor(meta$sample)
meta$processing <- factor(meta$processing); meta$sequencing_run <- factor(meta$run)
x <- check_design(meta, ~ processing + group)
check_design(meta, ~ processing + sequencing_run + group)
stopifnot(nrow(meta) == 91L, nlevels(meta$group) == 16L)

for (line in c("HT1080", "MCF10A")) for (genotype in c("wt", "p21_neg"))
  for (radiation in c("7gy", "UVA", "UVB")) {
    w <- response_weights(line, genotype, radiation)
    stopifnot(identical(unname(w), c(1, -1)))
    z <- group_contrast_matrix(meta, ~ processing + group, w)
    stopifnot(length(z) == ncol(x), abs(z["processingUKD"]) < 1e-9)
  }
for (line in c("HT1080", "MCF10A")) for (radiation in c("7gy", "UVA", "UVB")) {
  w <- interaction_weights(line, radiation)
  stopifnot(identical(unname(w), c(1, -1, -1, 1)))
  z <- group_contrast_matrix(meta, ~ processing + group, w)
  stopifnot(length(z) == ncol(x), abs(z["processingUKD"]) < 1e-9)
  run_z <- group_contrast_matrix(meta, ~ processing + sequencing_run + group, w)
  stopifnot(all(abs(run_z[grepl("^(processing|sequencing_run)", names(run_z))]) < 1e-9))
}
for (line in c("HT1080", "MCF10A")) for (pair in list(c("UVB", "UVA"), c("7gy", "UVA"), c("7gy", "UVB"))) {
  z <- group_contrast_matrix(meta, ~ processing + group,
    c(interaction_weights(line, pair[1]), -interaction_weights(line, pair[2])))
  stopifnot(length(z) == ncol(x), abs(z["processingUKD"]) < 1e-9)
}
for (radiation in c("7gy", "UVA", "UVB")) {
  z <- group_contrast_matrix(meta, ~ processing + group,
    c(interaction_weights("HT1080", radiation), -interaction_weights("MCF10A", radiation)))
  stopifnot(length(z) == ncol(x), abs(z["processingUKD"]) < 1e-9)
}

# Synthetic fitted coefficients: KO UVB response 2, WT UVB response 1;
# a large processing effect must cancel from the p21 contrast.
beta <- setNames(rep(0, ncol(x)), colnames(x)); beta["processingUKD"] <- 7
beta[paste0("group", group_name("HT1080", "p21_neg", "UVB"))] <- 2
beta[paste0("group", group_name("HT1080", "wt", "UVB"))] <- 1
stopifnot(sum(beta * group_contrast_matrix(meta, ~ processing + group, interaction_weights("HT1080", "UVB"))) == 1)
legacy_processing <- function(line, genotype, radiation) {
  trt <- group_name(line, genotype, radiation); ctrl <- group_name(line, genotype, "control")
  a <- x[meta$group == trt, , drop = FALSE]; b <- x[meta$group == ctrl, , drop = FALSE]
  (colMeans(unique(a)) - colMeans(unique(b)))["processingUKD"]
}
stopifnot(legacy_processing("HT1080", "wt", "UVA") == 0.5,
  legacy_processing("MCF10A", "wt", "UVA") == 1)
cat("PASS: inputs, 96 -> 91, full-rank designs, 12 responses, 6 p21, 6 quality and 3 cell-line contrasts; synthetic sign and legacy processing weights.\n")
