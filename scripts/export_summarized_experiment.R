#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("usage: export_summarized_experiment.R COLLAPSED_COUNTS.tsv[.gz] OUTPUT.rds SAMPLE_NAME [MAPPING.tsv.gz]")
}

suppressPackageStartupMessages({
    library(Matrix)
    library(S4Vectors)
    library(SummarizedExperiment)
})

counts_path <- normalizePath(args[[1]], mustWork = TRUE)
output_path <- args[[2]]
sample_name <- args[[3]]
mapping_path <- if (length(args) >= 4) normalizePath(args[[4]], mustWork = TRUE) else NA_character_

tbl <- read.delim(counts_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("barcode", "sample_index", "umi", "read_count")
if (!all(required %in% colnames(tbl))) {
    stop("collapsed count table must contain: ", paste(required, collapse = ", "))
}

feature_id <- paste(tbl$barcode, tbl$sample_index, tbl$umi, sep = "_")
if (anyDuplicated(feature_id)) stop("feature identifiers are not unique")

count_matrix <- sparseMatrix(
    i = seq_len(nrow(tbl)),
    j = rep.int(1L, nrow(tbl)),
    x = as.numeric(tbl$read_count),
    dims = c(nrow(tbl), 1L),
    dimnames = list(feature_id, sample_name)
)

row_data <- DataFrame(
    sequence = feature_id,
    barcode = tbl$barcode,
    sampleIndex = tbl$sample_index,
    umi = tbl$umi,
    row.names = feature_id
)
col_data <- DataFrame(Name = sample_name, row.names = sample_name)

se <- SummarizedExperiment(
    assays = list(counts = count_matrix),
    rowData = row_data,
    colData = col_data,
    metadata = list(
        countType = "reads",
        collapseScope = "barcode_only",
        preservedFields = c("sampleIndex", "umi"),
        sourceCollapsedCounts = counts_path,
        barcodeMapping = mapping_path,
        processingInfo = paste("Exported by pymutscan 0.3.1 on", Sys.time())
    )
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(se, output_path)
cat("saved", output_path, "with", nrow(se), "rows and", sum(assay(se, "counts")), "reads\n")
