# Copyright (C) 2026 Enterprise Breeding System
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# -------------------------------------------------------------------------------------
# Name             : validate_bioflow_object
# Description      : Check whether a generated Bioflow object can actually drive
#                    the three EBS interoperability use cases:
#                      1. Multi Trial Analysis (MTA)
#                      2. F1 qa/qc            (module string "gVerif")
#                      3. Pedigree qa/qc      (module string "qaPed")
# Maintainer       : Bioflow core development team
# ---------------------------------------------------------------------------------------------
#
# The point of this file is to turn the Shiny modules' silent/late failures into
# an explicit checklist at export time. Each requirement below was taken from
# the module or pipeline code that enforces it, and the source is named in the
# message so a failure is actionable.
# ---------------------------------------------------------------------------------------------

#' Fetch the genotype individual names, whatever the geno slot holds
#'
#' @param object A Bioflow result object.
#' @return Character vector, possibly empty.
vbo_geno_ind_names <- function(object) {
  geno <- object$data$geno
  if (is.null(geno)) return(character(0))
  out <- try(
    if (inherits(geno, "genlight")) adegenet::indNames(geno) else rownames(as.matrix(geno)),
    silent = TRUE
  )
  if (inherits(out, "try-error") || is.null(out)) character(0) else as.character(out)
}

#' Rename a pedigree to canonical parameter names, the way Bioflow does
#'
#' @param object A Bioflow result object.
#' @return The pedigree data.frame with canonical column names, or NULL.
vbo_canonical_pedigree <- function(object) {
  ped  <- object$data$pedigree
  meta <- object$metadata$pedigree
  if (is.null(ped) || is.null(meta) || nrow(ped) == 0) return(NULL)
  cols <- colnames(ped)
  idx  <- match(cols, meta$value)
  newNames <- ifelse(is.na(idx), cols, meta$parameter[idx])
  colnames(ped) <- newNames
  ped
}

#' Record one check result
#' @noRd
vbo_check <- function(checks, ok, requirement, detail = "") {
  rbind(checks, data.frame(
    ok          = isTRUE(ok),
    requirement = requirement,
    detail      = detail,
    stringsAsFactors = FALSE
  ))
}

#' Validate the object for Multi Trial Analysis
#'
#' Requirements read off bioflow/R/mod_mtaLMMsolveApp.R and
#' cgiarPipeline::metLMMsolver.
#' @noRd
validate_use_case_mta <- function(object) {
  checks <- NULL
  status <- object$status

  # The STA step now stamps module "soa" (older objects used "sta").
  staModules <- c("soa", "sta")
  # MTA's real requirements: phenotypic data is present, and SOA (Single
  # Occurrence/Trial Analysis) ran successfully on it.
  hasPheno <- !is.null(object$data$pheno) && nrow(object$data$pheno) > 0
  checks <- vbo_check(checks, hasPheno,
    "phenotypic data is present",
    "MTA is a phenotype analysis; with no phenotype data there is nothing to analyze")

  hasSta <- !is.null(status) && any(staModules %in% status$module)
  checks <- vbo_check(checks, hasSta,
    "SOA (Single Trial Analysis) has been run",
    "mod_mtaLMMsolveApp.R blocks the run when there is no 'soa' (STA) stamp in status")

  staIds <- if (hasSta) unique(status$analysisId[status$module %in% staModules]) else numeric(0)

  checks <- vbo_check(checks, hasSta && is.numeric(status$analysisId),
    "status$analysisId is numeric",
    "the STA dropdown calls as.POSIXct(analysisId); a character id yields NA labels")

  # One status row per sta analysisId, otherwise names<- errors and the
  # dropdown silently stays empty.
  dupSta <- hasSta && (sum(status$module %in% staModules) != length(staIds))
  checks <- vbo_check(checks, !dupSta,
    "no duplicate 'sta' status rows",
    "names(traitsMta) <- paste(analysisIdName, ...) errors on length mismatch")

  preds <- object$predictions
  predOk <- !is.null(preds) && nrow(preds) > 0 && any(preds$analysisId %in% staIds)
  checks <- vbo_check(checks, predOk,
    "predictions exist for the STA analysisId",
    "metLMMsolver stops with 'Not enough data ... perform an STA' when < 2 rows")

  requiredPredCols <- c("analysisId", "trait", "environment", "designation",
                        "predictedValue", "stdError", "reliability", "entryType")
  missingPredCols <- if (is.null(preds)) requiredPredCols else setdiff(requiredPredCols, colnames(preds))
  checks <- vbo_check(checks, length(missingPredCols) == 0,
    "predictions carry the columns MTA reads",
    if (length(missingPredCols)) paste("missing:", paste(missingPredCols, collapse = ", ")) else "")

  # metLMMsolver does predictionsBind[, colnames(phenoDTfile$predictions)], so a
  # non-canonical extra column makes the final rbind fail.
  canonicalPredCols <- c("module", "analysisId", "pipeline", "trait", "gid",
                         "designation", "mother", "father", "effectType",
                         "entryType", "environment", "predictedValue",
                         "stdError", "reliability")
  extraPredCols <- if (is.null(preds)) character(0) else setdiff(colnames(preds), canonicalPredCols)
  checks <- vbo_check(checks, length(extraPredCols) == 0,
    "predictions has no non-canonical extra column",
    if (length(extraPredCols)) paste("extra:", paste(extraPredCols, collapse = ", "),
                                     "- metLMMsolver's final rbind would fail") else "")

  # designationEffectType drives a branch that errors on NA.
  modeling <- object$modeling
  detOk <- !is.null(modeling) && any(modeling$parameter == "designationEffectType" &
                                       modeling$analysisId %in% staIds)
  checks <- vbo_check(checks, detOk,
    "modeling has designationEffectType for the STA id",
    "metLMMsolver evaluates if(... == 'BLUP'); absence gives 'missing value where TRUE/FALSE needed'")

  # Environment filtering is driven by metrics.
  metrics <- object$metrics
  h2Params <- c("plotH2", "H2", "meanR2", "r2",
                paste0(rep(c("plotH2", "H2", "meanR2", "r2"), each = 3),
                       c("_designation", "_mother", "_father")))
  metricOk <- !is.null(metrics) && any(metrics$parameter %in% h2Params &
                                         metrics$analysisId %in% staIds)
  checks <- vbo_check(checks, metricOk,
    "metrics carry an H2/r2 parameter for the STA id",
    "metLMMsolver keeps only environments whose heritability falls in [heritLB, heritUB]")

  mappedPheno <- c("environment", "designation", "trait")
  metaPheno   <- object$metadata$pheno
  missingMap  <- if (is.null(metaPheno)) mappedPheno else setdiff(mappedPheno, metaPheno$parameter)
  checks <- vbo_check(checks, length(missingMap) == 0,
    "metadata$pheno maps environment, designation and trait",
    if (length(missingMap)) paste("missing:", paste(missingMap, collapse = ", "),
                                  "- the MTA input panels stay hidden") else "")

  envColOk <- !is.null(object$data$pheno) && "environment" %in% colnames(object$data$pheno)
  checks <- vbo_check(checks, envColOk,
    "data$pheno has a literal 'environment' column",
    "cgiarPipeline::summaryWeather() aggregates on it when data$weather is absent")

  # Environments with a single prediction row per trait are dropped.
  nEnv <- 0
  if (predOk) {
    sub <- preds[preds$analysisId %in% staIds, , drop = FALSE]
    tab <- table(sub$environment, sub$trait)
    nEnv <- sum(apply(tab, 1, function(r) any(r > 1)))
  }
  checks <- vbo_check(checks, nEnv >= 2,
    "at least 2 environments with >1 prediction row",
    paste0("found ", nEnv, "; mod_mtaLMMsolveApp.R zeroes any env x trait cell <= 1"))

  checks
}

#' Validate the object for F1 qa/qc ("gVerif")
#'
#' Requirements read off bioflow/R/mod_hybridityApp.R and
#' cgiarPipeline::individualVerification.
#' @noRd
validate_use_case_f1 <- function(object) {
  checks <- NULL

  isGenlight <- inherits(object$data$geno, "genlight")
  checks <- vbo_check(checks, isGenlight,
    "data$geno is a genlight object",
    "individualVerification only builds Markers when class(geno)[1] == 'genlight'")

  qaGenoIds <- if (!is.null(object$status)) {
    as.character(object$status$analysisId[object$status$module == "qaGeno"])
  } else character(0)
  checks <- vbo_check(checks, length(qaGenoIds) > 0,
    "status contains a 'qaGeno' row",
    "the QA-geno stamp dropdown is built from status[module == 'qaGeno']")

  impNames <- names(object$data$geno_imp)
  keyOk <- length(qaGenoIds) > 0 && !is.null(impNames) && any(qaGenoIds %in% impNames)
  checks <- vbo_check(checks, keyOk,
    "data$geno_imp is keyed by the qaGeno analysisId",
    if (!keyOk) paste0("geno_imp names: ", paste(impNames, collapse = ", "),
                       " vs qaGeno id(s): ", paste(qaGenoIds, collapse = ", ")) else "")

  metaGeno  <- object$metadata$geno
  ploidyOk  <- !is.null(metaGeno) && "ploidity" %in% metaGeno$parameter
  checks <- vbo_check(checks, ploidyOk,
    "metadata$geno has a 'ploidity' row",
    "mod_hybridityApp.R reads ploidy from metadata$geno[parameter == 'ploidity'] (note the spelling)")

  ped <- vbo_canonical_pedigree(object)
  checks <- vbo_check(checks, !is.null(ped),
    "data$pedigree is present",
    "mod_hybridityApp.R shows 'Pedigree information is required to run this module'")

  hasCrossType <- !is.null(ped) && "crossType" %in% colnames(ped)
  checks <- vbo_check(checks, hasCrossType,
    "pedigree has a crossType column",
    "individualVerification stops without it")

  hasF1 <- hasCrossType && any(ped$crossType == "F1", na.rm = TRUE)
  checks <- vbo_check(checks, hasF1,
    "crossType contains the literal 'F1'",
    "the comparison is == \"F1\", case-sensitive")

  naCross <- hasCrossType && any(is.na(ped$crossType))
  checks <- vbo_check(checks, !naCross,
    "crossType has no missing values",
    "ped[ped$crossType == 'F1', ] uses == so NA rows become all-NA phantom rows")

  # ID overlap with the marker matrix: the single most common real failure.
  indNames <- vbo_geno_ind_names(object)
  if (!is.null(ped) && hasF1 && length(indNames) > 0) {
    f1 <- ped[which(ped$crossType == "F1"), , drop = FALSE]
    # Orientation-aware progeny key: do not assume the genlight is keyed by
    # sample_id merely because a sample_id column exists. Pick whichever key's
    # values actually overlap the marker-matrix individual names, so a
    # designation-keyed VCF (with or without a sample_id column) still matches.
    progenyKey <- "designation"
    if ("sample_id" %in% colnames(ped)) {
      nBySid   <- length(intersect(unique(stats::na.omit(f1$sample_id)),   indNames))
      nByDesig <- length(intersect(unique(stats::na.omit(f1$designation)), indNames))
      if (nBySid >= nByDesig && nBySid > 0) progenyKey <- "sample_id"
    }
    nProgeny <- length(intersect(unique(f1[[progenyKey]]), indNames))
    nMother  <- length(intersect(unique(stats::na.omit(f1$mother)), indNames))
    nFather  <- length(intersect(unique(stats::na.omit(f1$father)), indNames))

    checks <- vbo_check(checks, nProgeny > 0,
      paste0("F1 progeny (", progenyKey, ") found in the marker matrix"),
      paste0(nProgeny, " of ", length(unique(f1[[progenyKey]])),
             " matched; individualVerification stops with 'None of your progeny genotypes have marker information'"))
    checks <- vbo_check(checks, nMother > 0,
      "F1 mothers found in the marker matrix",
      paste0(nMother, " matched; parents must be genotyped too"))
    checks <- vbo_check(checks, nFather > 0,
      "F1 fathers found in the marker matrix",
      paste0(nFather, " matched; parents must be genotyped too"))
  }

  nMarkers <- if (isGenlight) adegenet::nLoc(object$data$geno) else 0
  checks <- vbo_check(checks, nMarkers >= 8,
    "at least 8 markers",
    paste0("found ", nMarkers, "; crossVerification returns NA probMatch below min_markers = 8"))

  checks
}

#' Validate the object for Pedigree qa/qc ("qaPed")
#'
#' Requirements read off bioflow/R/mod_qaPedApp.R.
#' @noRd
validate_use_case_ped <- function(object) {
  checks <- NULL

  hasGeno <- !is.null(object$data$geno)
  checks <- vbo_check(checks, hasGeno,
    "data$geno is present",
    "qaPed needs the raw genlight for loc.all and the marker matrix")

  locAllOk <- hasGeno && !is.null(try(object$data$geno$loc.all, silent = TRUE)) &&
    length(object$data$geno$loc.all) > 0
  checks <- vbo_check(checks, locAllOk,
    "data$geno has allele strings in loc.all",
    "plot_impossible() splits glgeno$loc.all on '/' to build IUPAC triplets")

  qaGenoIds <- if (!is.null(object$status)) {
    as.character(object$status$analysisId[object$status$module == "qaGeno"])
  } else character(0)
  checks <- vbo_check(checks, length(qaGenoIds) > 0,
    "status contains a 'qaGeno' row",
    "req(input$version2qaPed) blocks the module until a qaGeno stamp exists")

  impNames <- names(object$data$geno_imp)
  keyOk <- length(qaGenoIds) > 0 && !is.null(impNames) && any(qaGenoIds %in% impNames)
  checks <- vbo_check(checks, keyOk,
    "data$geno_imp is keyed by the qaGeno analysisId",
    "qaPed resolves geno_imp by which(names(geno_imp) == version2qaPed)")

  meta <- object$metadata$pedigree
  ped  <- object$data$pedigree

  fatherVal <- if (!is.null(meta) && "father" %in% meta$parameter) {
    meta$value[meta$parameter == "father"]
  } else NULL
  fatherOk <- !is.null(fatherVal) && length(fatherVal) > 0 &&
    !is.na(fatherVal[1]) && nzchar(fatherVal[1])
  checks <- vbo_check(checks, fatherOk,
    "metadata$pedigree maps 'father' to a non-empty column",
    "is_valid_fatherCol() gates the module's green 'data is complete' message")

  requiredPed <- c("designation", "mother", "father")
  pedCanon <- vbo_canonical_pedigree(object)
  missingPed <- if (is.null(pedCanon)) requiredPed else setdiff(requiredPed, colnames(pedCanon))
  checks <- vbo_check(checks, length(missingPed) == 0,
    "pedigree resolves designation, mother and father",
    if (length(missingPed)) paste("missing:", paste(missingPed, collapse = ", ")) else "")

  # The positional-rename invariant. qaPed corrupts columns silently otherwise.
  alignOk <- FALSE
  alignDetail <- "pedigree or metadata missing"
  if (!is.null(ped) && !is.null(meta) && nrow(ped) > 0) {
    alignOk <- !inherits(
      try(assert_pedigree_alignment(ped, meta), silent = TRUE), "try-error"
    )
    if (!alignOk) {
      alignDetail <- "metadata$pedigree rows must match data$pedigree columns one-to-one, in order"
    } else {
      alignDetail <- ""
    }
  }
  checks <- vbo_check(checks, alignOk,
    "metadata$pedigree aligns positionally with data$pedigree",
    alignDetail)

  # Triplet overlap with the marker matrix.
  indNames <- vbo_geno_ind_names(object)
  if (!is.null(pedCanon) && length(missingPed) == 0 && length(indNames) > 0) {
    complete <- !is.na(pedCanon$designation) & !is.na(pedCanon$mother) & !is.na(pedCanon$father)
    genotyped <- complete &
      pedCanon$designation %in% indNames &
      pedCanon$mother %in% indNames &
      pedCanon$father %in% indNames
    nTriplets <- sum(genotyped)
    checks <- vbo_check(checks, nTriplets > 0,
      "at least one fully genotyped designation/mother/father triplet",
      paste0(nTriplets, " of ", nrow(pedCanon),
             " row(s) have all three individuals in the marker matrix; ",
             "without any, every metric is 'NO GENO DATA' and mclust::Mclust() errors"))

    # plot_GRM() runs mclust::Mclust(..., G = 2) on the matched triplets, which
    # needs at least 2 observations to fit a 2-component model.
    checks <- vbo_check(checks, nTriplets >= 2,
      "at least 2 fully genotyped triplets for clustering",
      paste0("found ", nTriplets,
             "; plot_GRM() calls mclust::Mclust(GLM_metric, G = 2), which errors with < 2 points"))

    # plot_impossible() / hist_wrongG only run the real comparison when the share
    # of rows with all three of designation/mother/father present clears the
    # module's gate: (count_complete * 100) / nrow(peddata) >= 0.8. (This is the
    # module's literal expression; mirror it so readiness matches module behaviour.)
    nComplete <- sum(complete)
    gateOk <- nrow(pedCanon) > 0 && ((nComplete * 100) / nrow(pedCanon)) >= 0.8
    checks <- vbo_check(checks, gateOk,
      "triplet-completeness gate for the comparison is met",
      paste0(nComplete, " of ", nrow(pedCanon),
             " row(s) have designation, mother and father all present; ",
             "plot_impossible() takes the 'Not enough matching information' branch otherwise"))
  }

  checks
}

# ---------------------------------------------------------------------------
# Human-friendly diagnostics (for EBS software engineers, not QG specialists)
# ---------------------------------------------------------------------------

#' One-line plain-language description of what each use case is for.
#' @noRd
VBO_USECASE_INFO <- list(
  "MTA (Multi Trial Analysis)" = list(
    id      = "mta",
    purpose = "Combine trials across environments to estimate genetic values (needs phenotypes and a Single Trial Analysis)."
  ),
  "F1 qa/qc (gVerif)" = list(
    id      = "f1_qaqc",
    purpose = "Check that each F1 hybrid is really the cross of its recorded mother and father, using marker data."
  ),
  "Pedigree qa/qc (qaPed)" = list(
    id      = "pedigree_qaqc",
    purpose = "Check recorded parentage of later-generation lines against marker data (designation/mother/father triplets)."
  )
)

#' Plain-language meaning + concrete fix action, keyed by the check's
#' `requirement` string. Target reader: an EBS software engineer preparing the
#' phenotype/pedigree/VCF files and running the Rscript. Keep this copy free of
#' QG/Bioflow internals; say what is wrong with the INPUT and what to change.
#'
#' Any check not listed here falls back to its technical `detail` string, so the
#' JSON is always complete even if copy is missing.
#' @noRd
VBO_FRIENDLY <- list(
  # ---- MTA -------------------------------------------------------------
  "phenotypic data is present" = list(
    meaning = "No phenotypic data was supplied, so there is nothing for MTA to analyze.",
    action  = "Provide a phenotype file; MTA is a phenotype analysis and cannot run without one."),
  "SOA (Single Trial Analysis) has been run" = list(
    meaning = "Single Trial Analysis (SOA) has not run on this object, so trial-combining (MTA) cannot start.",
    action  = "Provide phenotypic data so the script runs model-based QA and SOA automatically. If phenotypes are present but SOA still did not run, check the earlier log for the resMod/staLMM warning."),
  "status$analysisId is numeric" = list(
    meaning = "The analysis identifier is not a number, which the analysis stamp menu cannot read.",
    action  = "Re-generate the object with this script (it creates numeric analysis IDs); do not hand-edit the IDs."),
  "no duplicate 'sta' status rows" = list(
    meaning = "There is more than one Single Trial Analysis entry, which confuses the analysis menu.",
    action  = "Re-generate the object from a single run; avoid merging multiple exports."),
  "predictions exist for the STA analysisId" = list(
    meaning = "No trial results were produced, usually because phenotypes or traits were missing or too sparse.",
    action  = "Check that the phenotype file has the chosen traits with enough data, and run with runSOA = TRUE."),
  "predictions carry the columns MTA reads" = list(
    meaning = "The results table is missing columns the trial-combining analysis needs.",
    action  = "Re-generate with this script unmodified; if it persists, the phenotype file may be missing required fields."),
  "predictions has no non-canonical extra column" = list(
    meaning = "The results table has an unexpected extra column that breaks trial-combining.",
    action  = "Re-generate with this script unmodified; do not add columns to the exported object."),
  "modeling has designationEffectType for the STA id" = list(
    meaning = "The Single Trial Analysis did not record how genotype effects were modelled.",
    action  = "Re-run with runSOA = TRUE so Single Trial Analysis records this automatically."),
  "metrics carry an H2/r2 parameter for the STA id" = list(
    meaning = "No heritability/reliability metric was produced by Single Trial Analysis.",
    action  = "Re-run with runSOA = TRUE and ensure traits have enough replicated data."),
  "metadata$pheno maps environment, designation and trait" = list(
    meaning = "The phenotype column mapping is incomplete (environment, designation and trait must all be mapped).",
    action  = "Check the phenotype file has the expected EBS columns so the script can map environment, germplasmName and the traits."),
  "data$pheno has a literal 'environment' column" = list(
    meaning = "The phenotype data has no 'environment' column for grouping trials.",
    action  = "Provide the EBS phenotype fields the script builds 'environment' from (stage, year, season, site, experiment, occurrence)."),
  "at least 2 environments with >1 prediction row" = list(
    meaning = "Trial-combining needs results from at least two environments; fewer were found.",
    action  = "Provide phenotype data spanning at least two environments (occurrences)."),

  # ---- shared geno -----------------------------------------------------
  "data$geno is a genlight object" = list(
    meaning = "The genotype data was not loaded into the expected marker format.",
    action  = "Provide a valid VCF as the genotype file so the script can build the marker matrix."),
  "status contains a 'qaGeno' row" = list(
    meaning = "The genotype QA step is missing, so marker-based modules cannot start.",
    action  = "Provide a genotype (VCF) file; the script runs genotype QA automatically when one is supplied."),
  "data$geno_imp is keyed by the qaGeno analysisId" = list(
    meaning = "The imputed marker data is not linked to the genotype QA step.",
    action  = "Re-generate with this script unmodified (it links these automatically); provide a valid VCF."),
  "metadata$geno has a 'ploidity' row" = list(
    meaning = "The genotype ploidy was not recorded.",
    action  = "Re-generate with this script unmodified; provide a VCF with genotype calls it can read."),
  "at least 8 markers" = list(
    meaning = "There are too few markers for a reliable verification (at least 8 are needed).",
    action  = "Provide a VCF with at least 8 usable markers (more is better)."),

  # ---- F1 qa/qc --------------------------------------------------------
  "data$pedigree is present" = list(
    meaning = "No pedigree information was supplied.",
    action  = "Provide the pedigree CSV as the pedigree file."),
  "pedigree has a crossType column" = list(
    meaning = "The pedigree does not say which individuals are F1 hybrids.",
    action  = "Add a crossType column (or an 'entry'/type column) marking F1 hybrids; parents can be left as 'parent'."),
  "crossType contains the literal 'F1'" = list(
    meaning = "No rows are marked as F1 hybrids, so there is nothing to verify as a hybrid.",
    action  = "Label the hybrid rows as F1 (e.g. 'F1' or 'hybrid') in the crossType column. Note: labels like 'F6' are NOT treated as F1."),
  "crossType has no missing values" = list(
    meaning = "Some rows have a blank crossType, which corrupts the hybrid selection.",
    action  = "Fill every crossType cell (use 'parent' for non-hybrids); the script fills blanks with 'parent' automatically when it can."),
  "F1 progeny (sample_id) found in the marker matrix" = list(
    meaning = "The genotyped F1 hybrids could not be matched to the pedigree.",
    action  = "Make sure the F1 samples in the VCF match the pedigree (by sample_code when provided, otherwise by designation)."),
  "F1 progeny (designation) found in the marker matrix" = list(
    meaning = "The genotyped F1 hybrids could not be matched to the pedigree by designation.",
    action  = "Make sure the F1 sample names in the VCF match the pedigree designations (or provide a sample_code column that maps to them)."),
  "F1 mothers found in the marker matrix" = list(
    meaning = "The mothers of the F1 hybrids are not genotyped, so the cross cannot be verified.",
    action  = "Include the mother (female parent) genotypes in the VCF; their designations must match the pedigree 'mother' column."),
  "F1 fathers found in the marker matrix" = list(
    meaning = "The fathers of the F1 hybrids are not genotyped, so the cross cannot be verified.",
    action  = "Include the father (male parent) genotypes in the VCF; their designations must match the pedigree 'father' column."),

  # ---- Pedigree qa/qc --------------------------------------------------
  "data$geno is present" = list(
    meaning = "No genotype data was supplied for pedigree verification.",
    action  = "Provide the VCF genotype file."),
  "data$geno has allele strings in loc.all" = list(
    meaning = "The marker allele information needed for pedigree checks is missing.",
    action  = "Provide a standard VCF with REF/ALT alleles so the script can record allele strings."),
  "metadata$pedigree maps 'father' to a non-empty column" = list(
    meaning = "The pedigree does not map a father (male parent) column.",
    action  = "Add a father/male-parent column to the pedigree file."),
  "pedigree resolves designation, mother and father" = list(
    meaning = "The pedigree is missing one of the required columns: designation, mother or father.",
    action  = "Ensure the pedigree file has designation, mother and father columns (empty parents allowed)."),
  "metadata$pedigree aligns positionally with data$pedigree" = list(
    meaning = "The pedigree columns and their descriptions are out of step (an internal consistency problem).",
    action  = "Re-generate with this script unmodified; do not reorder or rename columns in the exported object."),
  "at least one fully genotyped designation/mother/father triplet" = list(
    meaning = "No line has its designation, mother AND father all genotyped, so no triplet can be checked.",
    action  = "Genotype the lines together with both their parents, and make sure all three names match the VCF (as designations)."),
  "at least 2 fully genotyped triplets for clustering" = list(
    meaning = "Fewer than two complete (line + both parents genotyped) triplets exist; the statistics need at least two.",
    action  = "Provide at least two lines that have both parents genotyped, with names matching the VCF (as designations)."),
  "triplet-completeness gate for the comparison is met" = list(
    meaning = "Too few pedigree rows list a designation, mother and father together for the comparison to run.",
    action  = "Fill in mother and father for the lines you want checked so most rows have all three recorded.")
)

#' Attach friendly meaning/action to a checks data.frame
#' @noRd
vbo_decorate_checks <- function(tbl) {
  n <- nrow(tbl)
  meaning <- character(n)
  action  <- character(n)
  for (i in seq_len(n)) {
    f <- VBO_FRIENDLY[[tbl$requirement[i]]]
    if (!is.null(f)) {
      meaning[i] <- f$meaning
      action[i]  <- f$action
    } else {
      # Fall back to the technical detail so nothing is ever blank.
      meaning[i] <- tbl$detail[i]
      action[i]  <- ""
    }
  }
  tbl$meaning <- meaning
  tbl$action  <- action
  tbl
}

#' Build the diagnostics payload (a plain list ready for jsonlite::toJSON)
#' @noRd
vbo_json_payload <- function(useCases, ready, inputs = NULL, execChecks = NULL) {
  overallReady <- as.logical(ready)
  readyNames   <- names(ready)[overallReady]

  useCaseList <- lapply(names(useCases), function(nm) {
    tbl  <- vbo_decorate_checks(useCases[[nm]])
    info <- VBO_USECASE_INFO[[nm]]
    isReady <- isTRUE(ready[[nm]])

    failed <- tbl[!tbl$ok, , drop = FALSE]
    checksOut <- lapply(seq_len(nrow(tbl)), function(i) {
      list(
        passed            = isTRUE(tbl$ok[i]),
        requirement       = tbl$requirement[i],
        what_it_means     = tbl$meaning[i],
        how_to_fix        = if (!tbl$ok[i]) tbl$action[i] else "",
        technical_detail  = tbl$detail[i]
      )
    })

    list(
      use_case        = nm,
      id              = if (!is.null(info)) info$id else NA_character_,
      what_it_checks  = if (!is.null(info)) info$purpose else "",
      ready           = isReady,
      status          = if (isReady) "READY" else "NOT READY",
      summary         = if (isReady) {
        paste0("This file can run ", nm, ".")
      } else {
        paste0("This file cannot run ", nm, " yet: ",
               nrow(failed), " requirement(s) not met. See 'how_to_fix' below.")
      },
      checks_passed   = sum(tbl$ok),
      checks_total    = nrow(tbl),
      blocking_issues = unname(failed$requirement),
      checks          = checksOut,
      execution_checks = {
        er <- if (!is.null(execChecks)) execChecks[[nm]] else NULL
        if (is.null(er) || nrow(er) == 0) {
          list()
        } else {
          lapply(seq_len(nrow(er)), function(i) {
            list(
              check  = er$label[i],
              id     = er$id[i],
              ran    = isTRUE(er$ran[i]),
              status = if (isTRUE(er$ran[i])) "RAN" else "FAILED",
              note   = er$note[i],
              detail = er$detail[i]
            )
          })
        }
      }
    )
  })

  list(
    schema_version = "1.0",
    generated_at   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    audience       = "EBS software engineers preparing Bioflow input files",
    inputs         = inputs,
    overall = list(
      ready_use_cases     = unname(readyNames),
      not_ready_use_cases = unname(names(ready)[!overallReady]),
      ready_count         = sum(overallReady),
      total_use_cases     = length(ready),
      verdict             = if (all(overallReady)) {
        "All use cases are ready."
      } else if (any(overallReady)) {
        paste0(sum(overallReady), " of ", length(ready),
               " use cases are ready. The rest list what to fix under 'use_cases'.")
      } else {
        "No use cases are ready yet. See each use case's 'how_to_fix' guidance."
      }
    ),
    use_cases = useCaseList
  )
}

#' Write the diagnostics JSON to disk
#' @noRd
write_diagnostics_json <- function(payload, path) {
  if (is.null(path) || !nzchar(path)) return(invisible(NULL))
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    cli::cli_warn("Package 'jsonlite' is not available; skipping diagnostics JSON.")
    return(invisible(NULL))
  }
  json <- jsonlite::toJSON(payload, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null")
  writeLines(json, path, useBytes = TRUE)
  cli::cli_inform(paste0("Diagnostics written to: ", path))
  invisible(path)
}

# ---------------------------------------------------------------------------
# Execution checks: actually attempt the target Bioflow/cgiarPipeline function
# for each use case and record whether it ran. These go beyond the static
# readiness checks by exercising the real analysis engines.
# ---------------------------------------------------------------------------

#' Record one execution-check result
#' @noRd
vbo_exec <- function(records, id, label, ran, detail = "", note = "") {
  rbind(records, data.frame(
    id = id, label = label, ran = isTRUE(ran),
    detail = detail, note = note, stringsAsFactors = FALSE
  ))
}

#' Short one-line error text from a try-error / condition.
#' @noRd
vbo_err1 <- function(x) {
  msg <- if (inherits(x, "try-error")) {
    attr(x, "condition")$message
  } else if (inherits(x, "condition")) {
    conditionMessage(x)
  } else as.character(x)
  msg <- paste(msg, collapse = " ")
  msg <- gsub("[[:space:]]+", " ", trimws(msg))
  if (nchar(msg) > 300) paste0(substr(msg, 1, 300), " ...") else msg
}

#' UC1 (MOA): try metLMMsolver (main effects + GBLUP) and metASREML.
#' Reports whether the heritability filter had to be bypassed for the model to
#' keep any environment.
#' @noRd
vbo_exec_mta <- function(object) {
  rec <- NULL
  status <- object$status
  soaId <- unique(status$analysisId[status$module %in% c("soa", "sta")])
  if (length(soaId) == 0) {
    return(vbo_exec(rec, "moa", "Run MOA", FALSE,
                    "No SOA (soa) stamp present; MOA has no single-trial results to combine."))
  }
  soaId <- soaId[length(soaId)]
  traitsSoa <- unique(object$predictions$trait[object$predictions$analysisId %in% soaId])
  traitsSoa <- traitsSoa[!is.na(traitsSoa)]
  if (length(traitsSoa) == 0) {
    return(vbo_exec(rec, "moa", "Run MOA", FALSE, "No SOA predictions with a trait to analyze."))
  }
  trait1   <- traitsSoa[1]
  genoId   <- unique(status$analysisId[status$module == "qaGeno"])
  genoId   <- if (length(genoId)) genoId[length(genoId)] else NULL
  haveGeno <- !is.null(object$data$geno)

  # helper: run metLMMsolver with a given filter and (optional) GBLUP kernel;
  # returns list(ok, npred, err)
  tryMet <- function(bypass, gblup) {
    args <- list(phenoDTfile = object, analysisId = soaId, trait = trait1,
                 fixedTerm = list("1"), randomTerm = list("designation"),
                 verbose = FALSE)
    if (bypass) { args$heritLB <- 0; args$heritUB <- 1 }
    if (gblup) {
      args$expCovariates  <- list("genoA")
      args$analysisIdGeno <- genoId
      # metLMMsolver's marker-kernel path reads nPC["geno"]; pass 0 to use the
      # full kernel (no PCA reduction), mirroring the Bioflow UI default.
      args$nPC <- c(geno = 0)
    } else {
      args$expCovariates <- list("none")
    }
    out <- tryCatch(do.call(cgiarPipeline::metLMMsolver, args),
                    error = function(e) e, warning = function(w) w)
    if (inherits(out, "condition")) return(list(ok = FALSE, npred = 0, err = vbo_err1(out)))
    np <- tryCatch(sum(out$status$module == "moaLmms"), error = function(e) 0)
    list(ok = TRUE, npred = np, err = "")
  }

  # ---- main-effects LMMsolve: default filter, then bypass if needed --------
  d  <- tryMet(bypass = FALSE, gblup = FALSE)
  if (d$ok) {
    rec <- vbo_exec(rec, "moa_lmms_main", "MOA main-effects (LMMsolve)", TRUE,
                    "", "ran with default heritability filter [0.15, 0.95]")
  } else {
    b <- tryMet(bypass = TRUE, gblup = FALSE)
    if (b$ok) {
      rec <- vbo_exec(rec, "moa_lmms_main", "MOA main-effects (LMMsolve)", TRUE,
                      paste0("Default filter failed/empty (", d$err, ")"),
                      "ran ONLY after bypassing the heritability filter (heritLB=0, heritUB=1)")
    } else {
      rec <- vbo_exec(rec, "moa_lmms_main", "MOA main-effects (LMMsolve)", FALSE,
                      paste0("default filter: ", d$err, " | bypassed filter: ", b$err))
    }
  }

  # ---- GBLUP LMMsolve ------------------------------------------------------
  if (haveGeno && !is.null(genoId)) {
    d <- tryMet(bypass = FALSE, gblup = TRUE)
    if (d$ok) {
      rec <- vbo_exec(rec, "moa_lmms_gblup", "MOA GBLUP (LMMsolve)", TRUE, "",
                      "ran with default heritability filter [0.15, 0.95]")
    } else {
      b <- tryMet(bypass = TRUE, gblup = TRUE)
      if (b$ok) {
        rec <- vbo_exec(rec, "moa_lmms_gblup", "MOA GBLUP (LMMsolve)", TRUE,
                        paste0("Default filter failed/empty (", d$err, ")"),
                        "ran ONLY after bypassing the heritability filter (heritLB=0, heritUB=1)")
      } else {
        rec <- vbo_exec(rec, "moa_lmms_gblup", "MOA GBLUP (LMMsolve)", FALSE,
                        paste0("default filter: ", d$err, " | bypassed filter: ", b$err))
      }
    }
  } else {
    rec <- vbo_exec(rec, "moa_lmms_gblup", "MOA GBLUP (LMMsolve)", FALSE,
                    "No genotype/qaGeno stamp available for a marker kernel.")
  }

  # ---- ASREML main effects (if asreml is installed) ------------------------
  if (requireNamespace("asreml", quietly = TRUE)) {
    # asreml must be attached: metASREML relies on asreml's family helpers
    # (e.g. asr_gaussian) being on the search path, not just the namespace.
    suppressWarnings(suppressMessages(try(attachNamespace("asreml"), silent = TRUE)))
    # metASREML needs covMod/addG lists aligned to randomTerm ("none" = no kernel).
    aOut <- tryCatch(
      cgiarPipeline::metASREML(phenoDTfile = object, analysisId = soaId, trait = trait1,
                               fixedTerm = list("1"), randomTerm = list("designation"),
                               covMod = list("none"), addG = list("none"),
                               heritLB = 0, heritUB = 1, verbose = FALSE),
      error = function(e) e, warning = function(w) w
    )
    if (inherits(aOut, "condition")) {
      rec <- vbo_exec(rec, "moa_asreml_main", "MOA main-effects (ASREML)", FALSE,
                      vbo_err1(aOut))
    } else {
      rec <- vbo_exec(rec, "moa_asreml_main", "MOA main-effects (ASREML)", TRUE, "",
                      "ran with heritability filter bypassed (heritLB=0, heritUB=1)")
    }
  } else {
    rec <- vbo_exec(rec, "moa_asreml_main", "MOA main-effects (ASREML)", FALSE,
                    "asreml package is not installed in this environment.")
  }
  rec
}

#' UC2 (F1 qa/qc): try cgiarPipeline::individualVerification end to end.
#' @noRd
vbo_exec_f1 <- function(object) {
  rec <- NULL
  qaGenoId <- as.character(object$status$analysisId[object$status$module == "qaGeno"])
  if (length(qaGenoId) == 0) {
    return(vbo_exec(rec, "f1_verif", "Run F1 qa/qc verification", FALSE,
                    "No qaGeno stamp; individualVerification needs the imputed markers."))
  }
  qaGenoId <- qaGenoId[length(qaGenoId)]
  ploidy <- tryCatch(as.numeric(object$metadata$geno$value[object$metadata$geno$parameter == "ploidity"]),
                     error = function(e) 2)
  if (length(ploidy) == 0 || is.na(ploidy[1])) ploidy <- 2
  markers <- tryCatch(adegenet::locNames(object$data$geno_imp[[qaGenoId]]), error = function(e) NULL)

  out <- tryCatch(
    cgiarPipeline::individualVerification(
      object = object,
      analysisIdForGenoModifications = qaGenoId,
      markersToBeUsed = markers,
      colsForExpecGeno = c("mother", "father"),
      ploidy = ploidy[1],
      sc_filter = NULL, het = 10, matchThres = c(0.9, 0, 0.7)
    ),
    error = function(e) e, warning = function(w) w
  )
  if (inherits(out, "condition")) {
    rec <- vbo_exec(rec, "f1_verif", "Run F1 qa/qc verification", FALSE, vbo_err1(out))
  } else {
    gv <- tryCatch(out$predictions[out$predictions$module == "gVerif", , drop = FALSE],
                   error = function(e) NULL)
    nPM <- if (!is.null(gv)) sum(gv$trait == "probMatch" & !is.na(gv$predictedValue)) else 0
    rec <- vbo_exec(rec, "f1_verif", "Run F1 qa/qc verification", TRUE, "",
                    paste0("individualVerification ran; probMatch computed for ", nPM, " F1(s)."))
  }
  rec
}

#' UC3 (Pedigree qa/qc): reproduce the qaPed computation path from
#' bioflow/R/mod_qaPedApp.R (plot_impossible + plot_GRM) and report the exact
#' point of failure if it cannot complete.
#' @noRd
vbo_exec_ped <- function(object) {
  rec <- NULL

  if (!requireNamespace("bioflow", quietly = TRUE)) {
    return(vbo_exec(rec, "ped_qaped", "Run Pedigree qa/qc", FALSE,
                    "bioflow package not installed (needed for InputInfo/CompareTriplets data)."))
  }
  # mclust::Mclust dispatches to mclustBIC, which must be on the search path.
  suppressWarnings(suppressMessages(try(attachNamespace("mclust"), silent = TRUE)))
  qaGenoId <- object$status$analysisId[object$status$module == "qaGeno"]
  if (length(qaGenoId) == 0) {
    return(vbo_exec(rec, "ped_qaped", "Run Pedigree qa/qc", FALSE, "No qaGeno stamp."))
  }
  qaGenoId <- qaGenoId[length(qaGenoId)]

  res <- tryCatch({
    geno_imp <- cgiarBase::resolveGenoStamp(object$data$geno_imp, qaGenoId)
    glgeno   <- object$data$geno
    gi       <- object$data$geno_imp[[geno_imp]]
    # mimic modifyGeno: subset raw geno to imputed loci/inds
    ln <- which(adegenet::locNames(glgeno) %in% adegenet::locNames(gi))
    inn <- which(adegenet::indNames(glgeno) %in% adegenet::indNames(gi))
    geno <- glgeno[inn, ln]
    Markers <- as.data.frame(gi)

    peddata <- object$data$pedigree
    metaped <- object$metadata$pedigree

    # --- Step 1: triplet comparison (plot_impossible) ---
    uno <- data.frame(value = names(peddata))
    names(peddata) <- metaped$parameter[which(metaped$value %in% uno$value == TRUE)]
    resDMF <- lapply(split(peddata[, c("designation", "mother", "father")], 1:dim(peddata)[1]), unlist)
    dimres <- unlist(lapply(resDMF, length))
    gate <- (length(which(dimres == 3)) * 100) / dim(peddata)[1]
    if (gate < 0.8) {
      stop(paste0("triplet-completeness gate not met (", round(gate, 3),
                  " < 0.8): too few designation/mother/father triplets to compare."))
    }
    geno_matrix <- as.matrix(glgeno)
    nTrip <- sum(vapply(resDMF[which(dimres == 3)],
                        function(idx) all(idx %in% rownames(geno_matrix)), logical(1)))
    if (nTrip == 0) {
      stop("no fully genotyped designation/mother/father triplet (every triplet would be 'NO GENO DATA').")
    }

    # --- Step 2: GRM + clustering (plot_GRM) ---
    useful <- which(dimres == 3)
    pedU <- peddata[useful, ]
    # geno is a genlight (S4); convert to a numeric matrix before arithmetic.
    genoMat <- as.matrix(geno)
    ploidyFactor <- max(genoMat, na.rm = TRUE) / 2
    G <- sommer::A.mat(genoMat - ploidyFactor)
    resM <- lapply(split(pedU[, c("designation", "mother")], 1:dim(pedU)[1]), unlist)
    resM <- lapply(resM, function(x) match(x, rownames(G)))
    GLM_metric <- unlist(lapply(resM, function(x) G[x[1], x[2]]))
    GLM_metric <- GLM_metric[!is.na(GLM_metric)]
    if (length(GLM_metric) < 2) {
      stop(paste0("only ", length(GLM_metric),
                  " line-mother relationship value(s); mclust::Mclust(G=2) needs at least 2."))
    }
    fit <- mclust::Mclust(GLM_metric, G = 2, verbose = FALSE)
    if (is.null(fit)) stop("mclust::Mclust(G = 2) failed to fit two components.")
    list(nTrip = nTrip, nGLM = length(GLM_metric))
  }, error = function(e) e, warning = function(w) w)

  if (inherits(res, "condition")) {
    rec <- vbo_exec(rec, "ped_qaped", "Run Pedigree qa/qc", FALSE, vbo_err1(res))
  } else {
    rec <- vbo_exec(rec, "ped_qaped", "Run Pedigree qa/qc", TRUE, "",
                    paste0("qaPed computation ran: ", res$nTrip,
                           " genotyped triplet(s), ", res$nGLM, " clustering value(s)."))
  }
  rec
}

#' Dispatch execution checks for all use cases.
#' @noRd
run_execution_checks <- function(object) {
  list(
    "MTA (Multi Trial Analysis)" = vbo_exec_mta(object),
    "F1 qa/qc (gVerif)"          = vbo_exec_f1(object),
    "Pedigree qa/qc (qaPed)"     = vbo_exec_ped(object)
  )
}

#' Report whether a generated object can run the three EBS use cases
#'
#' @param object A Bioflow result object (the list saved as `result`).
#' @param verbose Print a human-readable report. Default TRUE.
#' @param diagnosticsFile Optional path. When given, an engineer-friendly
#'   diagnostics JSON (readiness per use case plus plain-language meaning and
#'   fix actions) is written there.
#' @param inputs Optional named list of input file paths recorded in the JSON.
#' @return Invisibly, a named list of data.frames (one per use case), each with
#'   `ok`, `requirement` and `detail` columns, plus a `ready` attribute.
validate_bioflow_object <- function(object, verbose = TRUE, diagnosticsFile = NULL, inputs = NULL, runExecution = FALSE) {

  useCases <- list(
    "MTA (Multi Trial Analysis)" = validate_use_case_mta(object),
    "F1 qa/qc (gVerif)"          = validate_use_case_f1(object),
    "Pedigree qa/qc (qaPed)"     = validate_use_case_ped(object)
  )

  ready <- vapply(useCases, function(x) all(x$ok), logical(1))
  for (nm in names(useCases)) attr(useCases[[nm]], "ready") <- ready[[nm]]

  # Optionally attempt to actually run each use case's target function.
  execChecks <- NULL
  if (isTRUE(runExecution)) {
    execChecks <- run_execution_checks(object)
  }

  if (isTRUE(verbose)) {
    cat("\n=========================================================\n")
    cat(" Bioflow use-case readiness report\n")
    cat("=========================================================\n")
    for (nm in names(useCases)) {
      tbl <- useCases[[nm]]
      cat("\n", ifelse(ready[[nm]], "[READY]    ", "[NOT READY]"), " ", nm, "\n", sep = "")
      for (i in seq_len(nrow(tbl))) {
        mark <- if (tbl$ok[i]) "  ok   " else "  FAIL "
        cat(mark, tbl$requirement[i], "\n", sep = "")
        if (!tbl$ok[i] && nzchar(tbl$detail[i])) {
          cat("         -> ", tbl$detail[i], "\n", sep = "")
        }
      }
    }
    cat("\n---------------------------------------------------------\n")
    cat("Summary: ", sum(ready), " of ", length(ready), " use case(s) ready\n", sep = "")
    cat("---------------------------------------------------------\n\n")

    if (!is.null(execChecks)) {
      cat("=========================================================\n")
      cat(" Execution checks (did the target function actually run?)\n")
      cat("=========================================================\n")
      for (nm in names(execChecks)) {
        er <- execChecks[[nm]]
        cat("\n", nm, "\n", sep = "")
        if (is.null(er) || nrow(er) == 0) { cat("  (no execution check)\n"); next }
        for (i in seq_len(nrow(er))) {
          mark <- if (er$ran[i]) "  RAN   " else "  FAILED"
          cat(mark, " ", er$label[i], "\n", sep = "")
          if (nzchar(er$note[i]))   cat("          note: ", er$note[i], "\n", sep = "")
          if (!er$ran[i] && nzchar(er$detail[i])) cat("          why : ", er$detail[i], "\n", sep = "")
          if (er$ran[i] && nzchar(er$detail[i]))  cat("          info: ", er$detail[i], "\n", sep = "")
        }
      }
      cat("\n---------------------------------------------------------\n\n")
    }
  }

  # Write the machine-readable, engineer-friendly diagnostics JSON.
  if (!is.null(diagnosticsFile) && nzchar(diagnosticsFile)) {
    payload <- vbo_json_payload(useCases, ready, inputs = inputs, execChecks = execChecks)
    write_diagnostics_json(payload, diagnosticsFile)
  }

  invisible(useCases)
}
