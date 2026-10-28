# Tests for the useBP opt-in gating (spec tpp-review Task 5.1)
#
# Property 8: useBP gating is fail-fast and side-effect-free on failure
# Validates: Requirements 3.2, 3.4, 3.5, 3.6
#   - useBP = FALSE reads no files and populates nothing, regardless of bpid/msid.
#   - useBP = TRUE with each unmet precondition writes NO RData and records the
#     appropriate diagnostic.
#   - "Usable bpid" requires >= 1 non-NA bpid (an msid-only column is not usable).
#
# The gate is exercised primarily as the pure helpers tpp_has_usable_bpid() and
# tpp_gate_decision() (no heavy pipeline needed). A single end-to-end check uses
# getBioflowRData() with the MOA sample + a tiny temp outputPath to assert no
# RData is written on failure branches. All temp outputs are cleaned up.
#
# Run with R 4.4.2 (testthat is not installed under 4.6.1):
#   & "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" -e "testthat::test_dir('tests/testthat')"

library(testthat)

# ---- locate scripts + real sample data relative to the repo ------------------
find_repo_root <- function(start) {
  d <- normalizePath(start, mustWork = TRUE)
  for (i in seq_len(8L)) {
    if (file.exists(file.path(d, "scripts", "tpp_gating_helpers.R"))) {
      return(d)
    }
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Could not locate EBS-hackathon-26 repo root from: ", start)
}

repo_root      <- find_repo_root(getwd())
workspace_root <- normalizePath(file.path(repo_root, ".."), mustWork = TRUE)
scripts_dir    <- file.path(repo_root, "scripts")

source(file.path(scripts_dir, "tpp_id_helpers.R"))
source(file.path(scripts_dir, "tpp_population_helpers.R"))
source(file.path(scripts_dir, "tpp_gating_helpers.R"))

moa_pheno <- file.path(repo_root, "use_cases", "1_MOA",
                       "bioflow_phenotype_data_MOA.csv")
bp_file   <- file.path(workspace_root, "PAM_module", "TPP_examples",
                       "TPP_traits_for_Bioflow.csv")
tr_file   <- file.path(repo_root, "use_cases", "1_MOA", "variable_tr_mapping.csv")

# A tiny in-memory pheno fixture: bpid usable, msid empty (mirrors MOA).
pheno_usable_bpid <- function() {
  data.frame(
    bpid = c("DBP0003", "DBP0003"),
    msid = c("", ""),
    WH_PH_M_CM = c(110, 95),
    stringsAsFactors = FALSE
  )
}

# ---- tpp_has_usable_bpid (Requirement 3.3) ----------------------------------

test_that("usable bpid requires the column present AND >= 1 non-NA value", {
  # Column present with a real value -> usable.
  expect_true(tpp_has_usable_bpid(pheno_usable_bpid()))

  # Column present but every value NA / blank -> NOT usable.
  expect_false(tpp_has_usable_bpid(
    data.frame(bpid = c(NA, NA), msid = c("5", "5"), stringsAsFactors = FALSE)))
  expect_false(tpp_has_usable_bpid(
    data.frame(bpid = c("", "  "), stringsAsFactors = FALSE)))

  # Column absent -> NOT usable (even with an msid column present).
  expect_false(tpp_has_usable_bpid(
    data.frame(msid = c("MS00054", "MS00054"), stringsAsFactors = FALSE)))

  # NULL / non-data.frame -> NOT usable.
  expect_false(tpp_has_usable_bpid(NULL))
})

test_that("an msid-only column is not usable (msid alone is insufficient)", {
  # msid present and populated, but no bpid column -> not usable.
  msid_only <- data.frame(msid = c("MS00054", "MS00054"), stringsAsFactors = FALSE)
  expect_false(tpp_has_usable_bpid(msid_only))
})

# ---- tpp_gate_decision: useBP = FALSE (Requirement 3.2) ---------------------

test_that("useBP = FALSE proceeds without populating and reads no files, regardless of bpid/msid", {
  # Even with fully usable bpid+msid, useBP=FALSE must not populate.
  dec <- tpp_gate_decision(
    useBP = FALSE,
    data_pheno = data.frame(bpid = "DBP0003", msid = "MS00054",
                            stringsAsFactors = FALSE),
    phenotypeFile = moa_pheno,
    BP_file = bp_file,
    TR_mapping_file = "variable_tr_mapping.csv"
  )
  expect_true(dec$proceed)
  expect_false(dec$populate)
  expect_identical(dec$status, "no_optin")

  # And with no bpid at all, still a clean no-op proceed.
  dec2 <- tpp_gate_decision(
    useBP = FALSE,
    data_pheno = data.frame(x = 1), phenotypeFile = moa_pheno,
    BP_file = NULL, TR_mapping_file = NULL
  )
  expect_true(dec2$proceed)
  expect_false(dec2$populate)
  expect_identical(dec2$status, "no_optin")
})

# ---- tpp_gate_decision: useBP = TRUE STOP branches --------------------------

test_that("useBP = TRUE with no usable bpid -> STOP (Requirement 3.4)", {
  dec <- tpp_gate_decision(
    useBP = TRUE,
    data_pheno = data.frame(msid = "MS00054", stringsAsFactors = FALSE),
    phenotypeFile = moa_pheno,
    BP_file = bp_file,
    TR_mapping_file = "variable_tr_mapping.csv"
  )
  expect_false(dec$proceed)
  expect_false(dec$populate)
  expect_identical(dec$status, "stop_no_usable_bpid")
  expect_identical(dec$requirement, "3.4")
  expect_match(dec$reason, "usable bpid")
})

test_that("useBP = TRUE, usable bpid, BP_file NULL/unreadable -> STOP (Requirement 3.5)", {
  # NULL BP_file.
  dec_null <- tpp_gate_decision(
    useBP = TRUE, data_pheno = pheno_usable_bpid(),
    phenotypeFile = moa_pheno, BP_file = NULL,
    TR_mapping_file = "variable_tr_mapping.csv"
  )
  expect_false(dec_null$proceed)
  expect_identical(dec_null$status, "stop_bp_file")
  expect_identical(dec_null$requirement, "3.5")

  # Unreadable BP_file path.
  dec_missing <- tpp_gate_decision(
    useBP = TRUE, data_pheno = pheno_usable_bpid(),
    phenotypeFile = moa_pheno,
    BP_file = file.path(tempdir(), "no-such-bp-file.csv"),
    TR_mapping_file = "variable_tr_mapping.csv"
  )
  expect_false(dec_missing$proceed)
  expect_identical(dec_missing$status, "stop_bp_file")
  expect_identical(dec_missing$requirement, "3.5")
})

test_that("useBP = TRUE, usable bpid, TR_mapping_file NULL/unresolvable -> STOP (Requirement 3.6)", {
  # NULL TR_mapping_file (BP_file is readable).
  dec_null <- tpp_gate_decision(
    useBP = TRUE, data_pheno = pheno_usable_bpid(),
    phenotypeFile = moa_pheno, BP_file = bp_file, TR_mapping_file = NULL
  )
  expect_false(dec_null$proceed)
  expect_identical(dec_null$status, "stop_tr_mapping_file")
  expect_identical(dec_null$requirement, "3.6")

  # Unresolvable filename in the phenotype directory.
  dec_missing <- tpp_gate_decision(
    useBP = TRUE, data_pheno = pheno_usable_bpid(),
    phenotypeFile = moa_pheno, BP_file = bp_file,
    TR_mapping_file = "not-a-real-mapping.csv"
  )
  expect_false(dec_missing$proceed)
  expect_identical(dec_missing$status, "stop_tr_mapping_file")
  expect_identical(dec_missing$requirement, "3.6")
})

# ---- tpp_gate_decision: proceed-to-populate (Requirement 3.7) ---------------

test_that("useBP = TRUE with usable bpid and both files readable -> populate", {
  dec <- tpp_gate_decision(
    useBP = TRUE, data_pheno = pheno_usable_bpid(),
    phenotypeFile = moa_pheno, BP_file = bp_file,
    TR_mapping_file = "variable_tr_mapping.csv"
  )
  expect_true(dec$proceed)
  expect_true(dec$populate)
  expect_identical(dec$status, "proceed_populate")
  # TR mapping path is resolved against the phenotype directory.
  expect_true(file.exists(dec$tr_mapping_path))
  expect_identical(normalizePath(dec$tr_mapping_path), normalizePath(tr_file))
})

# ---- gating diagnostics JSON writer -----------------------------------------

test_that("write_gating_diagnostics_json records the gating failure and writes a file", {
  skip_if_not_installed("jsonlite")
  dec <- tpp_gate_decision(
    useBP = TRUE, data_pheno = pheno_usable_bpid(),
    phenotypeFile = moa_pheno, BP_file = NULL,
    TR_mapping_file = "variable_tr_mapping.csv"
  )
  out <- file.path(tempdir(), "gating-diag-test.json")
  on.exit(unlink(out), add = TRUE)

  write_gating_diagnostics_json(dec, out)
  expect_true(file.exists(out))

  payload <- jsonlite::fromJSON(out)
  expect_false(payload$tpp_gating$rdata_written)
  expect_false(payload$tpp_gating$populated)
  expect_identical(payload$tpp_gating$status, "stop_bp_file")
  expect_identical(payload$tpp_gating$requirement, "3.5")
})

# ---- end-to-end: no RData is written on failure branches --------------------

test_that("getBioflowRData(useBP=TRUE) on MOA sample with missing BP_file writes NO RData", {
  skip_if_not(file.exists(moa_pheno), "MOA sample phenotype file not available")

  # Source the full (bundled) script so getBioflowRData is defined. Use the
  # bundle so this exercises exactly what ships.
  bundle <- file.path(scripts_dir, "bundled_getBioflowRdata.R")
  skip_if_not(file.exists(bundle), "bundled script not available")

  # ensure_* calls inside getBioflowRData need network/packages; we only care
  # that the gate stops BEFORE any RData is written. Source in a child env.
  env <- new.env()
  suppressWarnings(suppressMessages(sys.source(bundle, envir = env)))
  skip_if_not(exists("getBioflowRData", envir = env),
              "getBioflowRData not defined by the bundle")

  out_dir <- file.path(tempdir(), paste0("tpp_gate_e2e_", as.integer(Sys.time())))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  # useBP = TRUE, MOA has usable bpid (DBP0003) but BP_file is missing -> STOP.
  err <- tryCatch(
    env$getBioflowRData(
      phenotypeFile = moa_pheno,
      traits = "WH_PH_M_CM",
      outputPath = out_dir,
      outputFile = "gate_fail",
      validate = FALSE,
      runSOA = FALSE,
      useBP = TRUE,
      BP_file = file.path(tempdir(), "no-such-bp-file.csv"),
      TR_mapping_file = "variable_tr_mapping.csv"
    ),
    error = function(e) e
  )

  # The function must have stopped (fail-fast).
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "no RData written")

  # Side-effect-free: NO .RData file exists in the output directory.
  rdata_files <- list.files(out_dir, pattern = "\\.RData$", full.names = TRUE)
  expect_length(rdata_files, 0L)

  # But a gating diagnostics JSON WAS recorded.
  diag_files <- list.files(out_dir, pattern = "\\.diagnostics\\.json$", full.names = TRUE)
  expect_true(length(diag_files) >= 1L)
})
