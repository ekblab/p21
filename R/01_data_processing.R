## ************************************************************************************************
## p21 bulk RNA-seq
## 01 - Data import, annotation and differential expression model
##
## Current primary model (Model A):
##   ~ processing + group
## Alternative design for later sensitivity analysis (Model B):
##   ~ processing + sequencing_run + group
##   group = cell line x p21 genotype x irradiation (24 h)
##
## Decisions:
##   - Exclude the five UKD MCF10A WT controls, as prespecified.
##   - Remove synthetic ERCC/SIRV spike-ins BEFORE normalization and DESeq2.
##   - Assume shared additive technical effects (processing; optionally run).
##   - Only fit Model A here; evaluate Model B in the primary-DE script.
##   - Do not prefilter low counts beyond spike-ins; DESeq2 results() will apply
##     independent filtering to individual contrasts in the next script.
##   - Define DEGs later using padj < 0.05 (BH), without an LFC cutoff.
##
## Output:
##   Data/processed_data.rds - the ONE shared processed RDS file
##   Results/01_data_processing/tables/*.csv - transparent QC/audit tables
## ************************************************************************************************

rm(list = ls())
gc()
source("R/00_setup.R")

## ******
## input and output paths
counts_file <- "Data/p21_counts_counting_summary_summary_unique.tsv"
metadata_file <- "Data/metadata_v0.1_20250514.csv"
output_dir <- file.path("Results", "01_data_processing", "tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ************************************************************************************************
## Import count matrix
##
## Preserve the count-file ENSEMBL identifiers for matching to older results.
## The 'gene' column contains input symbols and is NOT part of the count matrix.
counts_input <- utils::read.delim(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("gene", "gene_id") %in% names(counts_input))) stop("Count file must contain 'gene' and 'gene_id'.")
if (anyNA(counts_input$gene_id) || any(!nzchar(as.character(counts_input$gene_id))) ||
    anyDuplicated(counts_input$gene_id)) stop("Invalid or duplicated feature IDs.")

sample_columns <- setdiff(names(counts_input), c("gene", "gene_id"))
counts_all <- as.matrix(counts_input[, sample_columns, drop = FALSE])
if (!is.numeric(counts_all) || anyNA(counts_all) || any(counts_all < 0) ||
    any(counts_all != floor(counts_all)) || any(counts_all > .Machine$integer.max))
  stop("Count matrix must contain non-negative integer counts without missing values.")
storage.mode(counts_all) <- "integer"
rownames(counts_all) <- as.character(counts_input$gene_id)
colnames(counts_all) <- ifelse(grepl("^X[0-9]+$", colnames(counts_all)),
                               colnames(counts_all), paste0("X", colnames(counts_all)))

feature_info <- data.frame(ENSEMBL = as.character(counts_input$gene_id),
                           input_symbol = as.character(counts_input$gene), stringsAsFactors = FALSE)
rm(counts_input)

## ************************************************************************************************
## Import metadata and verify sample correspondence
meta_all <- utils::read.csv(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
required_fields <- c("ID", "cell_type", "irradiation", "sample", "processing", "run")
if (!all(required_fields %in% names(meta_all))) stop("Missing required metadata fields.")
if (anyNA(meta_all$ID) || any(!nzchar(as.character(meta_all$ID)))) stop("Missing sample IDs in metadata.")
meta_all$ID <- ifelse(grepl("^X[0-9]+$", as.character(meta_all$ID)),
                      as.character(meta_all$ID), paste0("X", meta_all$ID))

if (anyNA(meta_all$ID) || anyDuplicated(meta_all$ID) || anyDuplicated(colnames(counts_all)))
  stop("Missing or duplicate sample IDs.")
if (!setequal(meta_all$ID, colnames(counts_all))) stop("Sample IDs in metadata and count matrix differ.")
counts_all <- counts_all[, meta_all$ID, drop = FALSE]
stopifnot(identical(colnames(counts_all), meta_all$ID), nrow(meta_all) == 96L)

## ******
## derive biologically interpretable variables and verify original group labels
meta_all$cell_line <- sub("_.*", "", meta_all$cell_type)
meta_all$genotype <- ifelse(grepl("p21_neg$", meta_all$cell_type), "p21_neg", "wt")
meta_all$radiation <- ifelse(meta_all$irradiation == "control", "control",
                             sub("^irradiated_", "", meta_all$irradiation))
meta_all$sequencing_run <- meta_all$run
meta_all$group <- as.character(meta_all$sample)

expected_groups <- group_name(meta_all$cell_line, meta_all$genotype, meta_all$radiation)
if (!identical(unname(expected_groups), meta_all$group)) stop("Original sample labels differ from expected 24-h groups.")
if (!setequal(meta_all$cell_line, c("HT1080", "MCF10A")) ||
    !setequal(meta_all$genotype, c("wt", "p21_neg")) ||
    !setequal(meta_all$radiation, c("control", "7gy", "UVA", "UVB")))
  stop("Unexpected experimental conditions in metadata.")

## ************************************************************************************************
## Exclude the previously identified MCF10A WT controls
excluded <- meta_all$processing == "UKD" & meta_all$cell_type == "MCF10A_wt" &
            meta_all$irradiation == "control"
if (sum(excluded) != 5L) stop("The expected five UKD MCF10A WT controls were not found.")
meta_all$excluded <- excluded

sample_decisions <- data.frame(ID = meta_all$ID, excluded = excluded,
                               reason = ifelse(excluded, "Prespecified UKD MCF10A WT control exclusion", ""))
utils::write.csv(sample_decisions, file.path(output_dir, "sample_decisions.csv"), row.names = FALSE)

## ************************************************************************************************
## Remove synthetic spike-ins BEFORE size-factor and dispersion estimation
##
## Match both identifiers and supplied symbols. Only ERCC-* synthetic IDs are
## removed: human DNA-repair genes such as ERCC1/ERCC2 must remain.
## SIRV2, SIRV5 and isoforms starting with SIRV1/SIRV2/... are also removed.
is_spike <- grepl("^ERCC[-_][0-9]+", feature_info$ENSEMBL, ignore.case = TRUE) |
            grepl("^ERCC[-_][0-9]+", feature_info$input_symbol, ignore.case = TRUE) |
            grepl("^SIRV([0-9]|[-_.]|$)", feature_info$ENSEMBL, ignore.case = TRUE) |
            grepl("^SIRV([0-9]|[-_.]|$)", feature_info$input_symbol, ignore.case = TRUE)

spike_table <- feature_info[is_spike, , drop = FALSE]
spike_table$total_counts <- rowSums(counts_all[is_spike, , drop = FALSE])
spike_table$n_detected_samples <- rowSums(counts_all[is_spike, , drop = FALSE] > 0)
utils::write.csv(spike_table, file.path(output_dir, "excluded_spike_ins.csv"), row.names = FALSE)

counts_endogenous <- counts_all[!is_spike, , drop = FALSE]
feature_info <- feature_info[!is_spike, , drop = FALSE]
if (!nrow(counts_endogenous)) stop("No features remain after spike-in filtering.")
stopifnot(identical(rownames(counts_endogenous), feature_info$ENSEMBL))
message("Excluded ", sum(is_spike), " ERCC/SIRV spike-in features; ", nrow(counts_endogenous), " features remain.")

## ******
## Record count-level metrics for all 96 samples, including excluded controls.
## Counts below are descriptive QC metrics and are not DESeq2-normalized.
sample_qc <- data.frame(ID = meta_all$ID, excluded = excluded,
                        raw_library_counts = colSums(counts_all),
                        endogenous_library_counts = colSums(counts_endogenous),
                        raw_detected_features = colSums(counts_all > 0),
                        endogenous_detected_features = colSums(counts_endogenous > 0))
utils::write.csv(sample_qc, file.path(output_dir, "sample_count_summary.csv"), row.names = FALSE)

## ************************************************************************************************
## Construct analysis dataset (91 samples)
meta <- meta_all[!excluded, , drop = FALSE]
rownames(meta) <- meta$ID
counts <- counts_endogenous[, meta$ID, drop = FALSE]
stopifnot(nrow(meta) == 91L, identical(colnames(counts), rownames(meta)))

meta$processing <- factor(meta$processing, levels = c("EKB", "UKD"))
meta$sequencing_run <- factor(meta$sequencing_run)
meta$cell_line <- factor(meta$cell_line, levels = c("HT1080", "MCF10A"))
meta$genotype <- factor(meta$genotype, levels = c("wt", "p21_neg"))
meta$radiation <- factor(meta$radiation, levels = c("control", "7gy", "UVA", "UVB"))
meta$group <- factor(meta$group)
if (anyNA(meta$processing) || anyNA(meta$group) || nlevels(meta$group) != 16L)
  stop("Invalid factors or missing experimental groups.")

## ******
## Export complete design count table, including zero-frequency combinations
study_design <- as.data.frame(table(meta$cell_line, meta$genotype, meta$radiation,
                                    meta$processing, meta$sequencing_run))
names(study_design) <- c("cell_line", "genotype", "radiation", "processing", "sequencing_run", "n")
utils::write.csv(study_design, file.path(output_dir, "design_counts.csv"), row.names = FALSE)

## ************************************************************************************************
## Feature annotation using the installed org.Hs.eg.db version
##
## Database keys are mapped using ENSEMBL IDs without version suffixes, but
## original feature identifiers are retained. No online biomaRt query is needed.
## For one-to-many mappings, prefer a mapped symbol and keep one record per gene.
annotation <- feature_info
annotation$ENSEMBL_base <- sub("\\.[0-9]+$", "", annotation$ENSEMBL)
annotation$SYMBOL <- NA_character_
annotation$ENTREZID <- NA_character_
annotation$GENENAME <- NA_character_
annotation$GENETYPE_AnnoDBI <- NA_character_

orgdb <- org.Hs.eg.db::org.Hs.eg.db
mapping_columns <- intersect(c("SYMBOL", "ENTREZID", "GENENAME", "GENETYPE"),
                             AnnotationDbi::columns(orgdb))
ensembl_keys <- unique(annotation$ENSEMBL_base[grepl("^ENSG[0-9]+$", annotation$ENSEMBL_base)])
if (length(ensembl_keys)) {
  mapping <- suppressMessages(AnnotationDbi::select(orgdb, keys = ensembl_keys,
                                                      keytype = "ENSEMBL", columns = mapping_columns))
  mapping <- mapping[order(mapping$ENSEMBL, is.na(mapping$SYMBOL), mapping$SYMBOL, na.last = TRUE), ]
  mapping <- mapping[!duplicated(mapping$ENSEMBL), , drop = FALSE]
  idx <- match(annotation$ENSEMBL_base, mapping$ENSEMBL)
  for (field in intersect(c("SYMBOL", "ENTREZID", "GENENAME"), names(mapping)))
    annotation[[field]] <- as.character(mapping[[field]][idx])
  if ("GENETYPE" %in% names(mapping))
    annotation$GENETYPE_AnnoDBI <- as.character(mapping$GENETYPE[idx])
}

missing_symbol <- is.na(annotation$SYMBOL) | annotation$SYMBOL == ""
annotation$SYMBOL[missing_symbol] <- annotation$input_symbol[missing_symbol]
stopifnot(!anyDuplicated(annotation$ENSEMBL), identical(annotation$ENSEMBL, rownames(counts)))
utils::write.csv(annotation, file.path(output_dir, "feature_annotation.csv"), row.names = FALSE)

## ************************************************************************************************
## Define and validate both statistical models
##
## Model A (current primary): One shared, additive processing effect.
## Model B (sensitivity): Additionally adjust for the sequencing run.
## Both nuisance effects are held constant in later treatment contrasts.
## Rank checks establish estimability, not that batch effects are truly additive.
primary_design <- ~ processing + group
alternative_design <- ~ processing + sequencing_run + group

primary_matrix <- check_design(meta, primary_design)
alternative_matrix <- check_design(meta, alternative_design)

model_audit <- data.frame(
  model = c("processing + group", "processing + sequencing_run + group"),
  formula = c(deparse(primary_design), deparse(alternative_design)),
  n_samples = c(nrow(primary_matrix), nrow(alternative_matrix)),
  n_coefficients = c(ncol(primary_matrix), ncol(alternative_matrix)),
  matrix_rank = c(qr(primary_matrix)$rank, qr(alternative_matrix)$rank))
utils::write.csv(model_audit, file.path(output_dir, "model_design_audit.csv"), row.names = FALSE)

## ************************************************************************************************
## Fit only Model A at the processing stage
##
## Keep the existing DESeq2 fit for QC and as the reference for downstream work.
## Model B will be RE-FITTED from the endogenous counts in 02_primary_DE.R;
## checking its design here does not constitute fitting or validating the model.
## No batch-corrected count matrix is used for statistical testing.
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = meta, design = primary_design)
dds <- DESeq2::DESeq(dds, quiet = TRUE, parallel = FALSE)

## ******
## Unadjusted VST for exploratory QC
## Both processing-only and processing+run display corrections are computed
## later in the QC report, never applied to the DESeq2 count matrix.
vsd <- DESeq2::vst(dds, blind = FALSE)

## ************************************************************************************************
## Save the single processed object used by all subsequent scripts
processed_data <- list(dds = dds, vsd = vsd, annotation = annotation, metadata_all = meta_all)
saveRDS(processed_data, "Data/processed_data.rds")

message("Processing complete: ", ncol(dds), " samples, ", nrow(dds),
        " endogenous features. Saved Data/processed_data.rds")
