# Tests for the full-path TPP object build (spec tpp-review Task 6.1).
#
# Property 2: Full-path build produces exactly the upload element set and shapes.
# Validates: Requirements 4.2, 4.3, 4.4
#   - exactly-one MS match -> metadata$TPP[[id]] has exactly the element set
#     {traits, env_filters, scale_validations, checks_per_trait, checks} with the
#     upload-produced shapes (traits 3-col df; checks 3-col df; env_filters NULL
#     or a named list keyed by Trait ID).
#   - zero MS matches  -> skip + "no matching TPP" diagnostic.
#   - >1 distinct TPP  -> skip + "1:1 MS<->TPP violated" diagnostic.
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

# Expected element set an interactive upload produces for metadata$TPP[[id]].
upload_element_set <- c("traits", "env_filters", "scale_validations",
                        "checks_per_trait", "checks")

# ---- Requirement 4.4 / exactly one TPP: upload element set + shapes ----------
test_that("full-path build on TPP00697 yields the exact upload element set", {
  skip_if_no_sample()
  tpp <- read_tpp_file(tpp_file)
  tr  <- read_tr_mapping(tr_mapping_file)

  # MS00054 (BP00003) is the realistic single-TPP market segment -> TPP00697.
  built <- populate_tpp_full_path(tpp, tr, "MS00054", pheno_designations = character(0))

  expect_equal(built$status, "populated")
  expect_equal(built$tpp_id, "TPP00697")
  expect_true(grepl("^TPP[0-9]{5}$", built$tpp_id))

  meta <- built$metadata
  expect_true(is.list(meta))
  # EXACTLY the upload element set (same names, no extras, no omissions).
  expect_setequal(names(meta), upload_element_set)
})

test_that("full-path build produces upload-shaped traits and checks data.frames", {
  skip_if_no_sample()
  tpp <- read_tpp_file(tpp_file)
  tr  <- read_tr_mapping(tr_mapping_file)
  built <- populate_tpp_full_path(tpp, tr, "MS00054", pheno_designations = character(0))
  meta <- built$metadata

  # traits: 3-column data.frame (tpp_trait, pheno_trait, trait_id).
  expect_s3_class(meta$traits, "data.frame")
  expect_equal(colnames(meta$traits), c("tpp_trait", "pheno_trait", "trait_id"))

  # Only TR-matched traits are included (Req 4.5). The MOA TR mapping resolves
  # TR-00016 (GYDW_CMP_THA) and TR-00019 (MAT_CALC_DAY), both in TPP00697.
  expect_true(all(c("TR-00016", "TR-00019") %in% meta$traits$trait_id))
  expect_true(all(nzchar(meta$traits$pheno_trait)))
  expect_setequal(meta$traits$pheno_trait, c("GYDW_CMP_THA", "MAT_CALC_DAY"))

  # checks: 3-column data.frame (tpp_check, tpp_trait, pheno_designation).
  expect_s3_class(meta$checks, "data.frame")
  expect_equal(colnames(meta$checks),
               c("tpp_check", "tpp_trait", "pheno_designation"))

  # env_filters: NULL or a named list keyed by Trait ID. Full path keeps all
  # envs -> NULL.
  expect_true(is.null(meta$env_filters) ||
                (is.list(meta$env_filters) && !is.null(names(meta$env_filters))))
  expect_null(meta$env_filters)

  # scale_validations / checks_per_trait: named lists keyed by Trait Name.
  expect_true(is.list(meta$scale_validations))
  expect_true(is.list(meta$checks_per_trait))
})

test_that("data$TPP sheet carries the Bioflow upload headers, checks standardized", {
  skip_if_no_sample()
  tpp <- read_tpp_file(tpp_file)
  tr  <- read_tr_mapping(tr_mapping_file)
  built <- populate_tpp_full_path(tpp, tr, "MS00054", pheno_designations = character(0))

  sheet <- built$data
  expect_s3_class(sheet, "data.frame")
  expect_true(all(c("Trait ID", "Trait Name", "Trait Requirement",
                    "Desired Score", "Scale Option", "Best Checks") %in% names(sheet)))
  # One row per TPP00697 trait (21 in the sample).
  expect_equal(nrow(sheet), 21L)

  # The checks column is standardized: ";" is converted to ", " (no bare ";").
  best <- sheet[["Best Checks"]]
  best_nonNA <- best[!is.na(best)]
  if (length(best_nonNA) > 0) {
    expect_false(any(grepl(";", best_nonNA)))
  }
})

# ---- Requirement 4.7 / 4.8: exact-match check auto-mapping -------------------
test_that("checks auto-map only on an exact pheno-designation match", {
  skip_if_no_sample()
  tpp <- read_tpp_file(tpp_file)
  tr  <- read_tr_mapping(tr_mapping_file)

  # TPP00697's mapped traits reference the best-check string
  # "INIAF Okinawa & INIAF Tropical (Bolivia)" and "Cortazar, alondra, ...".
  # Feed one exact designation so it auto-maps; an unrelated one stays unmapped.
  pheno_desigs <- c("Cortazar", "SOME_OTHER_LINE")
  built <- populate_tpp_full_path(tpp, tr, "MS00054", pheno_designations = pheno_desigs)
  meta <- built$metadata

  mapped_desigs <- unlist(meta$checks_per_trait, use.names = FALSE)
  # Only exact matches appear; the non-present designation never maps.
  expect_false("SOME_OTHER_LINE" %in% mapped_desigs)
  # Any mapped designation must be one we actually supplied.
  if (length(mapped_desigs) > 0) {
    expect_true(all(mapped_desigs %in% pheno_desigs))
  }

  # With NO pheno designations, nothing auto-maps (checks_per_trait empty).
  built0 <- populate_tpp_full_path(tpp, tr, "MS00054", pheno_designations = character(0))
  expect_equal(length(built0$metadata$checks_per_trait), 0L)
})

# ---- Requirement 4.2: zero MS matches -> skip + diagnostic -------------------
test_that("zero MS matches skips and records a no-matching-TPP diagnostic", {
  skip_if_no_sample()
  tpp <- read_tpp_file(tpp_file)
  tr  <- read_tr_mapping(tr_mapping_file)

  built <- populate_tpp_full_path(tpp, tr, "MS99999", pheno_designations = character(0))
  expect_equal(built$status, "no_match")
  expect_null(built$data)
  expect_null(built$metadata)
  expect_true(grepl("no matching TPP", built$diagnostic, ignore.case = TRUE))
})

# ---- Requirement 4.3: >1 distinct TPP -> skip + 1:1-violated diagnostic ------
test_that(">1 distinct TPP for a single msid skips with a 1:1-violated diagnostic", {
  # Construct a synthetic TPP file where ONE market segment (MSId 54) maps to
  # TWO distinct TPPs, violating the 1:1 rule.
  synthetic <- data.frame(
    TraitId          = c("TR-00016", "TR-00019"),
    TraitName        = c("Yield", "Maturity"),
    TraitType        = c("Agronomic", "Phenological"),
    TraitRequirement = c("Essential: Improve", "Essential: Maintain"),
    ScaleOption      = c("tons/ha", "days"),
    DesiredScore     = c("Higher than 5", "Higher than 70"),
    BestCheck        = c(NA_character_, NA_character_),
    ValueProposition = c("x", "y"),
    TPPId            = c("697", "999"),
    TPPCode          = c("TPP00697", "TPP00999"),
    BPId             = c("3", "3"),
    BPCode           = c("BP00003", "BP00003"),
    MSId             = c("54", "54"),
    MSShortName      = c("Wheat Spring LAC MS00054", "Wheat Spring LAC MS00054"),
    stringsAsFactors = FALSE
  )
  tr <- read_tr_mapping(tr_mapping_file)

  built <- populate_tpp_full_path(synthetic, tr, "MS00054",
                                  pheno_designations = character(0))
  expect_equal(built$status, "not_1to1")
  expect_null(built$data)
  expect_null(built$metadata)
  expect_true(grepl("1:1", built$diagnostic, fixed = TRUE))
  expect_setequal(built$tpp_codes, c("TPP00697", "TPP00999"))
})

# ---- structural-equivalence sanity: an upload-equivalent trait row -----------
test_that("built traits rows carry trait_id verbatim from the TPP TraitId", {
  skip_if_no_sample()
  tpp <- read_tpp_file(tpp_file)
  tr  <- read_tr_mapping(tr_mapping_file)
  built <- populate_tpp_full_path(tpp, tr, "MS00054", pheno_designations = character(0))

  # trait_id equals the TPP TraitId (not a name), per Requirement 2.
  expect_true(all(grepl("^TR-[0-9]{5}$", built$metadata$traits$trait_id)))
})
