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

  hasSta <- !is.null(status) && "sta" %in% status$module
  checks <- vbo_check(checks, hasSta,
    "status contains a 'sta' row",
    "mod_mtaLMMsolveApp.R blocks the run when sum(status$module %in% 'sta') == 0")

  staIds <- if (hasSta) unique(status$analysisId[status$module == "sta"]) else numeric(0)

  checks <- vbo_check(checks, hasSta && is.numeric(status$analysisId),
    "status$analysisId is numeric",
    "the STA dropdown calls as.POSIXct(analysisId); a character id yields NA labels")

  # One status row per sta analysisId, otherwise names<- errors and the
  # dropdown silently stays empty.
  dupSta <- hasSta && (sum(status$module == "sta") != length(staIds))
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
    progenyKey <- if ("sample_id" %in% colnames(ped)) "sample_id" else "designation"
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
  }

  checks
}

#' Report whether a generated object can run the three EBS use cases
#'
#' @param object A Bioflow result object (the list saved as `result`).
#' @param verbose Print a human-readable report. Default TRUE.
#' @return Invisibly, a named list of data.frames (one per use case), each with
#'   `ok`, `requirement` and `detail` columns, plus a `ready` attribute.
validate_bioflow_object <- function(object, verbose = TRUE) {

  useCases <- list(
    "MTA (Multi Trial Analysis)" = validate_use_case_mta(object),
    "F1 qa/qc (gVerif)"          = validate_use_case_f1(object),
    "Pedigree qa/qc (qaPed)"     = validate_use_case_ped(object)
  )

  ready <- vapply(useCases, function(x) all(x$ok), logical(1))
  for (nm in names(useCases)) attr(useCases[[nm]], "ready") <- ready[[nm]]

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
  }

  invisible(useCases)
}
