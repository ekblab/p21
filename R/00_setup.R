## ************************************************************************************************
## p21 bulk RNA-seq
## 00 - Project setup
##
## Run scripts from the project root (the directory containing R/ and Data/).
## Check essential dependencies and record the software environment.
## This script does not load data, run DESeq2 or write intermediate RDS files.
## ************************************************************************************************

## ******
## project structure
if (!file.exists("R/project_helpers.R")) stop("Run this script from the p21 project root.")

input_files <- c("Data/p21_counts_counting_summary_summary_unique.tsv", "Data/metadata_v0.1_20250514.csv")
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) stop("Missing input files: ", paste(missing_files, collapse = ", "))

## ******
## packages required for processing and the planned primary DE analysis
## Packages needed only by later analyses will be checked in those scripts.
required_packages <- c("DESeq2", "SummarizedExperiment", "AnnotationDbi", "org.Hs.eg.db", "ashr", "ekbSeq")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
if (utils::packageVersion("ekbSeq") < "1.2.0") stop("ekbSeq >= 1.2.0 is required.")

source("R/project_helpers.R")

## ******
## output directories and reproducibility
setup_dir <- file.path("Results", "00_setup")
dir.create(setup_dir, recursive = TRUE, showWarnings = FALSE)

package_versions <- vapply(required_packages, function(pkg) as.character(utils::packageVersion(pkg)), character(1))
ekb_sha <- utils::packageDescription("ekbSeq")$GithubSHA1
if (is.null(ekb_sha) || is.na(ekb_sha) || !nzchar(ekb_sha)) ekb_sha <- "not recorded"

writeLines(c(capture.output(utils::sessionInfo()), "", "Required package versions:",
             paste(names(package_versions), package_versions, sep = " = "),
             paste("ekbSeq GitHub SHA:", ekb_sha)), file.path(setup_dir, "session_info.txt"))

message("p21 setup successful: essential packages and input files found.")
