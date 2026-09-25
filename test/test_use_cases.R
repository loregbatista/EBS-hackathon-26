# End-to-end verification of the three EBS -> Bioflow use cases.
#
# This goes further than the readiness report: it feeds the generated object to
# the actual Bioflow pipeline function behind F1 qa/qc
# (cgiarPipeline::individualVerification) so a pass means Bioflow really accepts
# the file, not just that our own checks are satisfied.

source("scripts/bundle_scripts.R")
bundle_bioflow_scripts()
source("scripts/bundled_getBioflowRdata.R")
source("test/make_coherent_fixture.R")

cat("\n### Building a coherent EBS-style fixture ###\n")
fx <- make_coherent_fixture(outDir = "test")
cat("  phenotypes: ", fx$phenotypeFile, "\n", sep = "")
cat("  pedigree  : ", fx$pedigreeFile, "\n", sep = "")
cat("  genotypes : ", fx$genotypeFile, " (", fx$nMarkers, " markers)\n", sep = "")
cat("  mislabelled F1s (should be flagged by QA): ",
    paste(fx$mislabelled, collapse = ", "), "\n", sep = "")

cat("\n### Running getBioflowRData ###\n")
out <- getBioflowRData(
  phenotypeFile = fx$phenotypeFile,
  pedigreeFile  = fx$pedigreeFile,
  genotypeFile  = fx$genotypeFile,
  traits        = fx$trait,
  outputPath    = "test",
  outputFile    = "coherent_bioflow_input",
  pedigreeMapping = list(crossType = "entry", yearOfOrigin = "year")
)
result <- out$result

failures <- character(0)
expect <- function(condition, what) {
  if (isTRUE(condition)) cat("  ok   ", what, "\n", sep = "")
  else { cat("  FAIL ", what, "\n", sep = ""); failures <<- c(failures, what) }
}

# ---------------------------------------------------------------------------
cat("\n### Use-case readiness ###\n")
report <- validate_bioflow_object(result, verbose = FALSE)
for (nm in names(report)) {
  ready <- isTRUE(attr(report[[nm]], "ready"))
  expect(ready, paste0(nm, " is ready"))
  if (!ready) {
    bad <- report[[nm]][!report[[nm]]$ok, ]
    for (i in seq_len(nrow(bad))) {
      cat("         -> ", bad$requirement[i], ": ", bad$detail[i], "\n", sep = "")
    }
  }
}

# ---------------------------------------------------------------------------
cat("\n### Pedigree passthrough ###\n")
ped <- result$data$pedigree
srcPed <- read.csv(fx$pedigreeFile, check.names = FALSE)
cat("  source columns: ", paste(colnames(srcPed), collapse = ", "), "\n", sep = "")
cat("  object columns: ", paste(colnames(ped), collapse = ", "), "\n", sep = "")
expect(all(c("seed_source", "nursery_code") %in% colnames(ped)),
       "non-canonical columns were carried through verbatim")
expect(identical(as.character(ped$seed_source), as.character(srcPed$seed_source)),
       "carried-through column values are unchanged")
expect("other" %in% colnames(ped),
       "'plant_no' was recognised as Bioflow's 'other' parameter")
expect(identical(as.integer(ped$other), as.integer(srcPed$plant_no)),
       "'other' preserved the plant numbers")
expect(all(c("designation","sample_id","mother","father","crossType","yearOfOrigin")
           %in% colnames(ped)),
       "all canonical parameters resolved")
expect(identical(
         result$metadata$pedigree$value[result$metadata$pedigree$value %in% colnames(ped)],
         colnames(ped)),
       "metadata/column order invariant holds")
expect(sum(ped$crossType == "F1") == length(fx$f1) + 2,
       "F1 rows counted correctly (12 designations + 2 extra samples)")

# ---------------------------------------------------------------------------
cat("\n### Running the real Bioflow F1 qa/qc pipeline function ###\n")
qaGenoId <- as.character(result$status$analysisId[result$status$module == "qaGeno"])
verified <- try(
  cgiarPipeline::individualVerification(
    object = result,
    analysisIdForGenoModifications = qaGenoId,
    markersToBeUsed = adegenet::locNames(result$data$geno_imp[[qaGenoId]]),
    colsForExpecGeno = c("mother", "father"),
    ploidy = as.numeric(result$metadata$geno$value[result$metadata$geno$parameter == "ploidity"]),
    sc_filter = NULL,
    het = 10,
    matchThres = c(0.9, 0, 0.7)
  ),
  silent = TRUE
)

expect(!inherits(verified, "try-error"),
       "cgiarPipeline::individualVerification() ran on the generated object")
if (inherits(verified, "try-error")) {
  cat("         -> ", as.character(verified), "\n", sep = "")
} else {
  expect("gVerif" %in% verified$status$module, "a gVerif status row was written")
  gv <- verified$predictions[verified$predictions$module == "gVerif", ]
  pm <- gv[gv$trait == "probMatch", c("designation", "predictedValue")]
  pm <- pm[!is.na(pm$predictedValue), ]
  expect(nrow(pm) > 0, "probMatch was computed for at least one F1")

  correct <- pm[!pm$designation %in% fx$mislabelled, "predictedValue"]
  wrong   <- pm[pm$designation %in% fx$mislabelled, "predictedValue"]
  cat("  probMatch, correctly labelled F1s: mean ", round(mean(correct), 3),
      " (n=", length(correct), ")\n", sep = "")
  if (length(wrong) > 0) {
    cat("  probMatch, mislabelled F1s      : mean ", round(mean(wrong), 3),
        " (n=", length(wrong), ")\n", sep = "")
    expect(mean(correct) > mean(wrong),
           "correctly labelled F1s score higher than the mislabelled ones")
  }
}

# ---------------------------------------------------------------------------
cat("\n")
if (length(failures) > 0) {
  stop(sprintf("%d check(s) failed:\n  - %s", length(failures),
               paste(failures, collapse = "\n  - ")), call. = FALSE)
}
cat("All use-case checks passed.\n")
