# Test script for the EBS -> Bioflow conversion
# Runs getBioflowRData() on the sample data and checks the generated object
# against the three interoperability use cases: MTA, F1 qa/qc, Pedigree qa/qc.

# Define paths to test data files (relative to project root)
phenotype_file    <- "test/bioflow_pheno_data.csv"
pedigree_file     <- "test/PedF1.csv"
genotype_vcf_file <- "test/bioflow_genotype_data_fix.vcf"
output_file       <- "bioflow_input_test"
output_path       <- "test"

# Verify test data files exist
cat("Checking for test data files...\n")
if (!file.exists(phenotype_file)) {
  stop(sprintf("Phenotype file not found: %s", phenotype_file), call. = FALSE)
}
if (!file.exists(pedigree_file)) {
  stop(sprintf("Pedigree file not found: %s", pedigree_file), call. = FALSE)
}
if (!file.exists(genotype_vcf_file)) {
  stop(sprintf("Genotype VCF file not found: %s", genotype_vcf_file), call. = FALSE)
}
cat("All test data files found\n")

# Source and run the function to collect all the codebase into a
# single script
source("scripts/bundle_scripts.R")
bundle_bioflow_scripts()

# Source only the function definitions from the main script
source("scripts/bundled_getBioflowRdata.R")
#source("https://raw.githubusercontent.com/Breeding-Analytics/EBS-hackathon-26/refs/heads/main/scripts/bundled_getBioflowRdata.R")

# Run the function
cat("\nExecuting getBioflowRData...\n")
traits <- c("Maize_Plant_Height")

out <- getBioflowRData(
  phenotypeFile = phenotype_file,
  pedigreeFile  = pedigree_file,
  genotypeFile  = genotype_vcf_file,
  traits        = traits,
  outputPath    = output_path,
  outputFile    = output_file,
  # PedF1.csv calls the cross-type column "entry" and the origin year "year".
  # Auto-detection resolves both, but being explicit documents the intent.
  pedigreeMapping = list(crossType = "entry", yearOfOrigin = "year")
)

result   <- out$result
out_file <- out$file

# ---------------------------------------------------------------------------
# Verify output
# ---------------------------------------------------------------------------
cat("\nVerifying output...\n")

failures <- character(0)
expect <- function(condition, what) {
  if (isTRUE(condition)) {
    cat("  ok   ", what, "\n", sep = "")
  } else {
    cat("  FAIL ", what, "\n", sep = "")
    failures <<- c(failures, what)
  }
}

# the file was written where we asked for it
expect(file.exists(out_file), paste("output written to", out_file))
expect(basename(out_file) == paste0(output_file, ".RData"),
       "outputFile argument is honoured")

# pedigree made it in
ped <- result$data$pedigree
meta <- result$metadata$pedigree
expect(!is.null(ped) && nrow(ped) > 0, "data$pedigree is populated")
expect(all(c("designation", "mother", "father") %in% colnames(ped)),
       "pedigree has designation, mother and father")
expect("crossType" %in% colnames(ped), "pedigree has crossType")
expect(any(ped$crossType == "F1"), "crossType contains 'F1'")
expect(!any(is.na(ped$crossType)), "crossType has no NA values")
expect("sample_id" %in% colnames(ped), "pedigree carries sample_id")

# real parentage, not the old all-NA placeholder
expect(any(!is.na(ped$mother)), "mother column has real values")
expect(any(!is.na(ped$father)), "father column has real values")

# every column of the source file survived
src_cols <- colnames(read.csv(pedigree_file, check.names = FALSE))
cat("    source pedigree columns: ", paste(src_cols, collapse = ", "), "\n", sep = "")
cat("    object pedigree columns: ", paste(colnames(ped), collapse = ", "), "\n", sep = "")
expect(ncol(ped) >= length(src_cols),
       "no pedigree column was dropped on the way in")

# the invariant qaPed depends on
expect(nrow(meta) == ncol(ped),
       "metadata$pedigree has one row per pedigree column")
expect(identical(meta$value[meta$value %in% colnames(ped)], colnames(ped)),
       "metadata$pedigree row order matches pedigree column order")

# the two regressions fixed in this change
qa_id <- as.character(result$status$analysisId[result$status$module == "qaGeno"])
expect(!is.null(result$modifications$geno_imp) &&
         qa_id %in% names(result$modifications$geno_imp),
       "modifications$geno_imp is keyed by the qaGeno analysisId")
if (is.null(result$modifications$geno_imp[[qa_id]])) {
  cat("    (this VCF had no missing calls, so the log is legitimately NULL)\n")
}
sta_rows <- result$status[result$status$module == "sta", , drop = FALSE]
expect(nrow(sta_rows) == 1, "exactly one 'sta' status row")
expect(nrow(sta_rows) > 0 && !is.na(sta_rows$analysisIdName[1]) &&
         nzchar(sta_rows$analysisIdName[1]),
       "the 'sta' status row is named")

# geno_imp keyed by the qaGeno stamp
qa_ids <- as.character(result$status$analysisId[result$status$module == "qaGeno"])
expect(length(qa_ids) > 0 && all(qa_ids %in% names(result$data$geno_imp)),
       "data$geno_imp is keyed by the qaGeno analysisId")

# canonical predictions
expect("effectType" %in% colnames(result$predictions),
       "predictions carry effectType")

cat("\n")
if (length(failures) > 0) {
  stop(sprintf("%d check(s) failed:\n  - %s",
               length(failures), paste(failures, collapse = "\n  - ")), call. = FALSE)
}
cat("Test completed successfully.\n")
cat("Note: the readiness report above states which use cases this object can\n")
cat("drive. F1 qa/qc and Pedigree qa/qc additionally require the pedigree IDs\n")
cat("to match the VCF sample names.\n")
