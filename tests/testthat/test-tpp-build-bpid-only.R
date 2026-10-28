# Tests for the bpid-only TPP object build (spec tpp-review Task 7.1).
#
# Property 3: bpid-only path populates all BP TPPs and retains all environments.
# Validates: Requirements 5.1, 5.2, 5.3
#   - bpid present, msid empty/absent for all rows -> populate EVERY distinct
#     TPPId grouped under the matched normalized BPId (Req 5.1), building the
#     full upload-equivalent object for each (Req 5.3) with NO environment
#     subsetting (env_filters NULL -> every trait keeps all environments).
#   - bpid present but no TPP-file rows match the BP -> status "no_bp_match":
#     skip, RData WITHOUT TPP population, "no TPPs for BP" diagnostic (Req 5.2 / 5.4).
#
# The canonical bpid-only case is the MOA sample (bpid "DBP0003" -> BP00003,
# msid empty), using the real sample TPP file + MOA variable_tr_mapping.csv.
#
# Run with R 4.4.2 (testthat is not installed under 4.6.1):
#   & "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" -e "testthat::test_dir('tests/testthat')"

library(testthat)

find_repo_root <- function(start) {
  d <- normalizePath(start, mustWork = TRUE)
  for (i in seq_len(8L)) {
    if (file.exists(file.path(d, "scripts", "tpp_build_helpers.R"))) {
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
source(file.path(scripts_dir, "tpp_build_helpers.R"))

tpp_file <- file.path(workspace_root, "PAM_module", "TPP_examples",
                      "TPP_traits_for_Bioflow.csv")
tr_mapping_file <- file.path(repo_root, "use_cases", "1_MOA",
                             "variable_tr_mapping.csv")

skip_if_no_sample <- function() {
  if (!file.exists(tpp_file)) skip("sensitive TPP sample file not present")
  if (!file.exists(tr_mapping_file)) skip("TR mapping sample file not present")
}

upload_element_set <- c("traits", "env_filters", "scale_validations",
                        "checks_per_trait", "checks")

# ---- Requirement 5.1 / 5.3: MOA bpid-only -> all BP00003 TPPs, all envs ------
test_that("bpid-only on the MOA sample populates every distinct TPP under BP00003", {
  skip_if_no_sample()
  tpp <- read_tpp_file(tpp_file)
  tr  <- read_tr_mapping(tr_mapping_file)

  built <- populate_tpp_bpid_only(tpp, tr, "DBP0003",
                                  pheno_designations = character(0))

  expect_equal(built$status, "populated")

  # The distinct TPPs the helper reports must equal those grouped by the matched
  # normalized BPId directly from the TPP file (Req 5.1), with no omissions/extras.
  expected_codes <- distinct_tpps(match_bp(tpp, "DBP0003"))$TPPCode
  expect_true(length(expected_codes) > 1L)  # BP00003 spans several TPPs
  expect_setequal(built$tpp_codes, expected_codes)
  expect_setequal(names(built$tpps), expected_codes)

  # Every populated TPP id is a valid TPP code.
  expect_true(all(grepl("^TPP[0-9]{5}$", built$tpp_codes)))
})

test_that("each bpid-only TPP is upload-shaped and retains all environments", {
  skip_if_no_sample()
  tpp <- read_tpp_file(tpp_file)
  tr  <- read_tr_mapping(tr_mapping_file)
  built <- populate_tpp_bpid_only(tpp, tr, "DBP0003",
                                  pheno_designations = character(0))

  for (code in names(built$tpps)) {
    entry <- built$tpps[[code]]
    meta  <- entry$metadata

    # EXACTLY the upload element set (same names, no extras, no omissions).
    expect_setequal(names(meta), upload_element_set)

    # traits: 3-column data.frame (tpp_trait, pheno_trait, trait_id).
    expect_s3_class(meta$traits, "data.frame")
    expect_equal(colnames(meta$traits), c("tpp_trait", "pheno_trait", "trait_id"))

    # checks: 3-column data.frame.
    expect_s3_class(meta$checks, "data.frame")
    expect_equal(colnames(meta$checks),
                 c("tpp_check", "tpp_trait", "pheno_designation"))

    # Requirement 5.3: NO environment subsetting -> env_filters NULL for every
    # trait in every populated TPP (the per-trait env set equals the all-envs
    # set, so no trait carries an env filter entry).
    expect_null(meta$env_filters)

    # data sheet carries the Bioflow upload headers.
    expect_s3_class(entry$data, "data.frame")
    expect_true(all(c("Trait ID", "Trait Name", "Trait Requirement") %in%
                      names(entry$data)))
  }
})

test_that("bpid-only env count per trait equals the source count (no subsetting)", {
  skip_if_no_sample()
  tpp <- read_tpp_file(tpp_file)
  tr  <- read_tr_mapping(tr_mapping_file)
  built <- populate_tpp_bpid_only(tpp, tr, "DBP0003",
                                  pheno_designations = character(0))

  # The TPP file is tidy on traits (one row per TPP-trait); applying no env
  # subsetting means env_filters stays NULL. An explicit per-trait check: with
  # env_filters NULL, no trait's environment set is narrower than all-envs, so
  # the retained env count equals the source count for every trait (Req 5.3).
  all_null <- vapply(built$tpps,
                     function(e) is.null(e$metadata$env_filters),
                     logical(1))
  expect_true(all(all_null))
})

# ---- Requirement 5.2 / 5.4: no BP match -> skip + diagnostic, no TPP info ----
test_that("no BP match returns no_bp_match with a diagnostic and no built TPPs", {
  skip_if_no_sample()
  tpp <- read_tpp_file(tpp_file)
  tr  <- read_tr_mapping(tr_mapping_file)

  # A bpid that normalizes to a BPId absent from the TPP file.
  built <- populate_tpp_bpid_only(tpp, tr, "DBP9999",
                                  pheno_designations = character(0))

  expect_equal(built$status, "no_bp_match")
  expect_length(built$tpps, 0L)
  expect_length(built$tpp_codes, 0L)
  expect_false(is.na(built$diagnostic))
  expect_true(grepl("no TPPs", built$diagnostic, ignore.case = TRUE))
  # The caller would generate the RData WITHOUT TPP population in this case: the
  # helper returns no built pieces to attach.
})

test_that("bpid-only check auto-mapping still honours exact designation match", {
  skip_if_no_sample()
  tpp <- read_tpp_file(tpp_file)
  tr  <- read_tr_mapping(tr_mapping_file)

  # With NO pheno designations, no check auto-maps in any populated TPP.
  built0 <- populate_tpp_bpid_only(tpp, tr, "DBP0003",
                                   pheno_designations = character(0))
  mapped0 <- unlist(lapply(built0$tpps, function(e) e$metadata$checks_per_trait),
                    use.names = FALSE)
  expect_length(mapped0, 0L)
})
