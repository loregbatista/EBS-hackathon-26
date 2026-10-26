# Regression test for the orientation-aware F1 progeny key in
# validate_bioflow_object.R::validate_use_case_f1().
#
# Covers the three identifier-orientation cases that previously crashed F1
# qa/qc when a sample_id column existed but the VCF was keyed by designation:
#
#   1. VCF keyed by designation + pedigree WITHOUT sample_id
#   2. VCF keyed by sample_id    + pedigree WITH    sample_id
#   3. VCF keyed by designation  + pedigree WITH    sample_id  (the crash case)
#
# In every case the F1 progeny must be found in the marker matrix (nProgeny > 0)
# so cgiarPipeline::individualVerification can run.
#
# Run with:  Rscript test/test_progeny_key_orientation.R

source("scripts/validate_bioflow_object.R")

failures <- character(0)
expect <- function(condition, what) {
  if (isTRUE(condition)) cat("  ok   ", what, "\n", sep = "")
  else { cat("  FAIL ", what, "\n", sep = ""); failures <<- c(failures, what) }
}

# A tiny marker matrix whose rownames are the genotype individual names. This is
# all validate_use_case_f1's progeny-overlap block needs (vbo_geno_ind_names
# falls back to rownames(as.matrix(geno)) for a non-genlight object).
make_geno <- function(ind_names, n_markers = 10) {
  m <- matrix(0L, nrow = length(ind_names), ncol = n_markers)
  rownames(m) <- ind_names
  m
}

# Build a minimal Bioflow-style object. The pedigree is stored under canonical
# column names with an identity metadata mapping (matching how the script emits
# it), so vbo_canonical_pedigree returns it unchanged.
make_object <- function(ped, geno_names) {
  meta <- data.frame(parameter = colnames(ped), value = colnames(ped),
                     stringsAsFactors = FALSE)
  list(
    data = list(
      geno     = make_geno(geno_names),
      geno_imp = setNames(list(make_geno(geno_names)), "1"),
      pedigree = ped
    ),
    metadata = list(
      geno     = data.frame(parameter = "ploidity", value = "2",
                            stringsAsFactors = FALSE),
      pedigree = meta
    ),
    status = data.frame(module = "qaGeno", analysisId = "1",
                        stringsAsFactors = FALSE)
  )
}

# Pull the computed progeny overlap out of the readiness report. The relevant
# requirement string starts with "F1 progeny (".
progeny_row <- function(report_f1) {
  report_f1[grepl("^F1 progeny \\(", report_f1$requirement), , drop = FALSE]
}

# ----- Case 1: VCF designation, pedigree without sample_id -----------------
cat("\n### Case 1: VCF designation + pedigree without sample_id ###\n")
ped1 <- data.frame(
  designation = c("F1_A", "F1_B", "MOM1", "DAD1"),
  mother      = c("MOM1", "MOM1", NA, NA),
  father      = c("DAD1", "DAD1", NA, NA),
  crossType   = c("F1", "F1", "parent", "parent"),
  stringsAsFactors = FALSE
)
obj1 <- make_object(ped1, geno_names = c("F1_A", "F1_B", "MOM1", "DAD1"))
rep1 <- validate_use_case_f1(obj1)
pr1  <- progeny_row(rep1)
expect(nrow(pr1) == 1 && pr1$ok, "case 1: F1 progeny found via designation")

# ----- Case 2: VCF sample_id, pedigree with sample_id ----------------------
cat("\n### Case 2: VCF sample_id + pedigree with sample_id ###\n")
ped2 <- data.frame(
  designation = c("F1_A", "F1_B", "MOM1", "DAD1"),
  sample_id   = c("S001", "S002", "S100", "S200"),
  mother      = c("MOM1", "MOM1", NA, NA),
  father      = c("DAD1", "DAD1", NA, NA),
  crossType   = c("F1", "F1", "parent", "parent"),
  stringsAsFactors = FALSE
)
# Genotype matrix keyed by sample_id for the F1s (parents by designation, as
# individualVerification matches parents on designation).
obj2 <- make_object(ped2, geno_names = c("S001", "S002", "MOM1", "DAD1"))
rep2 <- validate_use_case_f1(obj2)
pr2  <- progeny_row(rep2)
expect(nrow(pr2) == 1 && pr2$ok, "case 2: F1 progeny found via sample_id")
expect(grepl("sample_id", pr2$requirement), "case 2: progeny key resolved to sample_id")

# ----- Case 3: VCF designation, pedigree WITH sample_id (crash case) --------
cat("\n### Case 3: VCF designation + pedigree with sample_id (crash case) ###\n")
ped3 <- ped2   # same pedigree, including the sample_id column
# But the VCF is keyed by DESIGNATION, not sample_id.
obj3 <- make_object(ped3, geno_names = c("F1_A", "F1_B", "MOM1", "DAD1"))
rep3 <- validate_use_case_f1(obj3)
pr3  <- progeny_row(rep3)
expect(nrow(pr3) == 1 && pr3$ok,
       "case 3: F1 progeny found via designation despite sample_id column")
expect(grepl("designation", pr3$requirement),
       "case 3: progeny key fell back to designation (no spurious sample_id switch)")

cat("\n")
if (length(failures) > 0) {
  stop(sprintf("%d check(s) failed:\n  - %s", length(failures),
               paste(failures, collapse = "\n  - ")), call. = FALSE)
}
cat("All progeny-key orientation checks passed.\n")
