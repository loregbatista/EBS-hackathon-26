# Tests for TPP file ingestion + BP/MS matching helpers (spec tpp-review Task 4.1)
#
# Validates: Requirements 2.2, 3.1
#   - BP/MS lookup on the real sample data: DBP0003 -> BP00003 -> TPP set (TPP00697).
#   - TR_Id resolver maps WH_PH_M_CM -> TR-00388 (and the other MOA mappings).
#
# Run with R 4.4.2 (testthat is not installed under 4.6.1):
#   & "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" -e "testthat::test_dir('tests/testthat')"

library(testthat)

# ---- locate scripts + real sample data relative to the repo ------------------
# This test file lives at EBS-hackathon-26/tests/testthat/. The repo root is two
# levels up; the workspace root (which holds PAM_module/) is one level above that.
# Walk up from the current working directory until the repo marker (scripts/
# tpp_population_helpers.R) is found, so the test works whether invoked from the
# repo root, the tests/ dir, or the testthat/ dir.
find_repo_root <- function(start) {
  d <- normalizePath(start, mustWork = TRUE)
  for (i in seq_len(8L)) {
    if (file.exists(file.path(d, "scripts", "tpp_population_helpers.R"))) {
      return(d)
    }
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Could not locate EBS-hackathon-26 repo root from: ", start)
}

repo_root <- find_repo_root(getwd())
workspace_root <- normalizePath(file.path(repo_root, ".."), mustWork = TRUE)

scripts_dir <- file.path(repo_root, "scripts")
source(file.path(scripts_dir, "tpp_id_helpers.R"))
source(file.path(scripts_dir, "tpp_population_helpers.R"))

tpp_file <- file.path(workspace_root, "PAM_module", "TPP_examples",
                      "TPP_traits_for_Bioflow.csv")
tr_mapping_file <- file.path(repo_root, "use_cases", "1_MOA",
                             "variable_tr_mapping.csv")

test_that("sample data files are present", {
  expect_true(file.exists(tpp_file))
  expect_true(file.exists(tr_mapping_file))
})

# ---- ingestion ---------------------------------------------------------------
test_that("read_tpp_file reads the TPP file and validates key columns", {
  tpp <- read_tpp_file(tpp_file)
  expect_s3_class(tpp, "data.frame")
  expect_true(nrow(tpp) > 0)
  for (col in c("TraitId", "TraitName", "TraitRequirement", "ScaleOption",
                "DesiredScore", "BestCheck", "TPPId", "TPPCode", "BPId",
                "BPCode", "MSId", "MSShortName")) {
    expect_true(col %in% names(tpp), info = col)
  }
})

test_that("read_tpp_file fails fast on a missing file", {
  expect_error(read_tpp_file(file.path(tempdir(), "does-not-exist.csv")),
               "not found")
})

# ---- Requirement 3.1 / BP + MS lookup ---------------------------------------
test_that("BP lookup: DBP0003 -> BP00003 -> TPP00697 (Requirement 3.1)", {
  tpp <- read_tpp_file(tpp_file)

  # Pheno bpid from the MOA sample.
  bp_rows <- match_bp(tpp, "DBP0003")
  expect_true(nrow(bp_rows) > 0)

  # Every matched row normalizes to BP key 3 and carries BPCode BP00003.
  expect_true(all(normalize_id(bp_rows$BPId) == 3L))
  expect_true(all(bp_rows$BPCode == "BP00003"))

  # The BP resolves to the expected TPP set (TPP00697 is the BP00003 TPP).
  tpps <- distinct_tpps(bp_rows)
  expect_true("TPP00697" %in% tpps$TPPCode)
})

test_that("BP lookup reconciles equivalent ID forms", {
  tpp <- read_tpp_file(tpp_file)
  n_prefixed <- nrow(match_bp(tpp, "DBP0003"))
  n_coded    <- nrow(match_bp(tpp, "BP00003"))
  n_numeric  <- nrow(match_bp(tpp, 3))
  expect_equal(n_prefixed, n_coded)
  expect_equal(n_prefixed, n_numeric)
  expect_true(n_prefixed > 0)
})

test_that("MS lookup: MS00054 -> TPP00697 and MS<->TPP is 1:1", {
  tpp <- read_tpp_file(tpp_file)

  # BP00003's market segment is MS00054 (MSId 54).
  ms_rows <- match_ms(tpp, "MS00054")
  expect_true(nrow(ms_rows) > 0)
  expect_true(all(normalize_id(ms_rows$MSId) == 54L))

  # MS <-> TPP is 1:1: a market-segment match resolves to exactly one TPP.
  tpps <- distinct_tpps(ms_rows)
  expect_equal(nrow(tpps), 1L)
  expect_equal(tpps$TPPCode, "TPP00697")
})

test_that("unmatched / empty IDs return zero rows, not an error", {
  tpp <- read_tpp_file(tpp_file)
  expect_equal(nrow(match_bp(tpp, "")), 0L)
  expect_equal(nrow(match_bp(tpp, NA)), 0L)
  expect_equal(nrow(match_bp(tpp, "DBP99999")), 0L)
  expect_equal(nrow(match_ms(tpp, "")), 0L)
})

# ---- Requirement 2.2 / TR_Id resolver ---------------------------------------
test_that("TR_Id resolver maps pheno variables to TR_Ids (Requirement 2.2)", {
  tr <- read_tr_mapping(tr_mapping_file)

  # Forward: pheno variable abbreviation -> TR_Id.
  expect_equal(tr$tr_id_for_variable("WH_PH_M_CM"),  "TR-00388")
  expect_equal(tr$tr_id_for_variable("GYDW_CMP_THA"), "TR-00016")
  expect_equal(tr$tr_id_for_variable("MAT_CALC_DAY"), "TR-00019")

  # Padding/empty rows are dropped.
  expect_equal(nrow(tr$mapping), 3L)
})

test_that("TR_Id resolver maps TR_Ids back to pheno variables via normalized key", {
  tr <- read_tr_mapping(tr_mapping_file)

  # Reverse resolution, robust to formatting via normalize_id.
  expect_equal(tr$variable_for_tr_id("TR-00388"), "WH_PH_M_CM")
  expect_equal(tr$variable_for_tr_id("388"),      "WH_PH_M_CM")   # normalized form
  expect_equal(resolve_pheno_trait(tr, "TR-00016"), "GYDW_CMP_THA")

  # A TR_Id with no pheno mapping resolves to NA (Requirement 2.4 leaves it empty).
  expect_true(is.na(tr$variable_for_tr_id("TR-99999")))
})

test_that("resolve_pheno_trait links a TPP-file TraitId to its pheno trait", {
  tpp <- read_tpp_file(tpp_file)
  tr  <- read_tr_mapping(tr_mapping_file)

  # TR-00388 (Plant Height) exists in the TPP file and maps to WH_PH_M_CM.
  expect_true(any(tpp$TraitId == "TR-00388"))
  expect_equal(resolve_pheno_trait(tr, "TR-00388"), "WH_PH_M_CM")
})
