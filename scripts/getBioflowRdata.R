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
# Name             : getBioflowRData
# Description      : Create RData that will be used by Bioflow, ready for
#                    Single Trial Analysis (STA), Multi Trial Analysis (MTA),
#                    F1 qa/qc and Pedigree qa/qc.
# R Version        : 4.2.3
# Pkg Dependency   : vcfR, adegenet, cli, rlang, remotes, openssl, stringr,
#                    purrr, glue, Matrix, cgiarBase, cgiarPipeline
# ---------------------------------------------------------------------------------------------
# Author           : Alaine A. Gulles
# Author Email     : a.gulles@cgiar.org
# Date             : 2024.10.15
# Date Modified    : 2026.09.24
# Maintainer       : Bioflow core development team
# Maintainer Email :
# Script Version   : 2.0
# Command          : getBioflowRData(phenotypeFile = "C:/sta_input_phenotype_v3_rev.csv",
#                             pedigreeFile  = "C:/pedigree.csv",
#                             genotypeFile  = "C:/markers.vcf",
#                             traits = c("YLD_TON","FLW50"),
#                             outputPath = "C:/results",
#                             outputFile = "Request1",
#                             requestId = "1010")
# ---------------------------------------------------------------------------------------------
#' @name getBioflowRData
#' @aliases getBioflowRData
#' @title Create R object for Bioflow with phenotypic, pedigree and genotypic QA/QC
#'
#' @description Creates the `result` object Bioflow loads, covering three EBS
#'   interoperability use cases:
#'   \enumerate{
#'     \item Multi Trial Analysis (MTA) - needs phenotypes; STA is run here.
#'     \item F1 qa/qc - needs pedigree (designation, mother, father, crossType,
#'           optionally sample_id) plus genotypes for progeny AND parents.
#'     \item Pedigree qa/qc - needs pedigree (designation, mother, father) plus
#'           genotypes for the designation/mother/father triplets.
#'   }
#'   A readiness report is printed at the end stating which of the three the
#'   generated file can actually drive.
#'
#' @param phenotypeFile a string indicating the path and file name in csv format containing the QCed phenotypic data for analysis from EBS
#' @param pedigreeFile NULL or string indicating the path and file name in csv format containing the pedigree data associated with the phenotype file. Required for the F1 qa/qc and Pedigree qa/qc use cases. Every column of this file is carried into the R object; designation, mother and father are mandatory.
#' @param genotypeFile NULL or string indicating the path and file name of VCF file containing the genotypic data associated with the phenotype file. Required for the F1 qa/qc and Pedigree qa/qc use cases; optional for MTA.
#' @param traits a character vector indicating the trait that will be use for analysis
#' @param outputPath string indicating the path where the RData will be save
#' @param outputFile string indicating the RData file name that will be used. When NULL the file is named with the md5 hash of the phenotype analysisId (the EBS content-addressed convention).
#' @param requestId string indicating the analysis request id generated in EBS
#' @param pedigreeMapping optional named list mapping Bioflow pedigree parameters to columns of `pedigreeFile`, e.g. `list(crossType = "entry", yearOfOrigin = "year")`. Only needed when auto-detection cannot resolve a column.
#' @param deriveCrossType logical; when `pedigreeFile` has no crossType column, infer one from parent completeness so the F1 qa/qc module can run. Default TRUE.
#' @param runSOA logical or NULL; run model-based outlier detection (`cgiarPipeline::resMod`, written as a `qaMb` stamp) followed by Single Trial Analysis (`cgiarPipeline::staLMM`, written as a `soa` stamp) so the object arrives MTA-ready. Default NULL = auto: runs whenever phenotypic data is present. Set FALSE to export a pedigree/genotype-only object for the qa/qc use cases.
#' @param validate logical; print the use-case readiness report. Default TRUE.
#' @param diagnosticsFile optional path for the diagnostics JSON. When NULL it is written beside the RData as `<outputFile>.diagnostics.json`. The JSON reports, per use case (MTA, F1 qa/qc, Pedigree qa/qc), whether the file is READY plus, for each unmet requirement, a plain-language meaning and a concrete fix action aimed at the EBS engineer preparing the files.
#'
#' @param runExecution logical; when TRUE, after validation the target Bioflow/cgiarPipeline analysis for each use case is actually attempted (MOA via metLMMsolver + metASREML, F1 via individualVerification, Pedigree qa/qc via the qaPed computation) and the outcome is added to the console log and the diagnostics JSON under `execution_checks`. Default FALSE. These run the real analysis engines, so they are slower.
#'
#' @param useBP logical master opt-in switch for automatic TPP (Target Product Profile) population from the phenotype `bpid`/`msid` columns. Default FALSE. When FALSE the script generates the RData exactly as before and never reads `BP_file` or `TR_mapping_file`, never populating any TPP information, regardless of whether the phenotype data carries `bpid`/`msid`. When TRUE the script requires usable `bpid` information (the `bpid` column present with at least one non-NA value) plus two readable file arguments, and is fail-fast: if any precondition is unmet it stops WITHOUT writing any RData file and records a diagnostic in the diagnostics JSON.
#' @param BP_file NULL or string giving the FULL path (including filename) to the sensitive breeding-portal TPP file (`TPP_traits_for_Bioflow.csv`, held only by the EBS team under `PAM_module/TPP_examples/`). Used verbatim. Read only when `useBP = TRUE` and the phenotype data has usable `bpid` information; required (non-NULL and readable) in that case. Default NULL.
#' @param TR_mapping_file NULL or string giving the FILENAME ONLY of the trait-mapping file (e.g. `"variable_tr_mapping.csv"`, columns `ebs_variable_abbrev` -> `bp_tr_id`). Resolved by the script as `file.path(dirname(phenotypeFile), TR_mapping_file)`, so it must sit in the SAME directory as the phenotype/genotype/pedigree inputs. Read only when `useBP = TRUE` and the phenotype data has usable `bpid` information; required (resolvable and readable) in that case. Default NULL.
#'
#' @return Invisibly, a list with `result` (the object written), `file` (the
#'   path it was written to) and `diagnosticsFile` (the diagnostics JSON path).
# -------------------------------------------------------------------------------------

getBioflowRData <- function(phenotypeFile = NULL, pedigreeFile = NULL, genotypeFile = NULL,
                     traits = NULL, outputPath, outputFile = NULL, requestId = NULL,
                     pedigreeMapping = NULL, deriveCrossType = TRUE,
                     runSOA = NULL, validate = TRUE, diagnosticsFile = NULL,
                     runExecution = FALSE,
                     useBP = FALSE, BP_file = NULL, TR_mapping_file = NULL) {

  # Now we can use ensure_cran_packages since packages_verification.R has been sourced
  ensure_cran_packages(c("vcfR", "adegenet", "cli", "rlang", "remotes",
                         "openssl", "stringr", "purrr", "glue", "Matrix",
                         "jsonlite"))
  # Install if not cgiarPipelines and cgiarBase using remotes
  ensure_github_packages()

  # If null dont bind
  status <- NULL
  analysisIdPheno <- NULL
  analysisIdGeno <- NULL

  havePheno <- !is.null(phenotypeFile) && nzchar(phenotypeFile)
  haveGeno <- !is.null(genotypeFile) && nzchar(genotypeFile)

  # runSOA controls model-based QA + Single Trial Analysis. Default (NULL) is
  # "auto": run whenever phenotypic data is present. TRUE/FALSE force it.
  if (is.null(runSOA)) {
    runSOA <- havePheno
  }

  # --------------------------------------------------------------------------
  # Opt-in TPP population gate (Requirement 3, spec tpp-review Task 5).
  #
  # Evaluated BEFORE any RData is built or written so that a gating failure is
  # genuinely fail-fast and side-effect-free: on any STOP branch NO RData file
  # is written, a diagnostic is recorded in the diagnostics JSON, and the
  # function stops. When useBP = FALSE this is a zero-cost no-op (no files are
  # read) and the RData is generated exactly as before.
  #
  # `tpp_gate_decision()` / `tpp_has_usable_bpid()` /
  # `write_gating_diagnostics_json()` live in scripts/tpp_gating_helpers.R
  # (bundled before this file). The decision is retained in `tppGate` and
  # consumed at the population hand-off point further below.
  tppGate <- NULL
  if (isTRUE(useBP)) {
    # Lightweight raw read solely for the "usable bpid" check. The main pheno
    # read below is unchanged; this does not alter it. We only read when the
    # user has opted in, so the useBP = FALSE path reads nothing.
    gatePheno <- if (havePheno) {
      tryCatch(
        read.csv(phenotypeFile, encoding = "utf-8", check.names = FALSE),
        error = function(e) NULL
      )
    } else {
      NULL
    }

    tppGate <- tpp_gate_decision(
      useBP           = useBP,
      data_pheno      = gatePheno,
      phenotypeFile   = phenotypeFile,
      BP_file         = BP_file,
      TR_mapping_file = TR_mapping_file
    )

    if (!isTRUE(tppGate$proceed)) {
      # Fail-fast STOP: write NO RData. Record the gating failure in a
      # diagnostics JSON at the same path the normal diagnostics would use.
      gatingDiagFile <- diagnosticsFile
      if (is.null(gatingDiagFile) || !nzchar(gatingDiagFile)) {
        gatingBase <- if (!is.null(outputFile) && nzchar(outputFile)) {
          outputFile
        } else {
          "tpp_gating_failure"
        }
        gatingDiagFile <- file.path(outputPath, paste0(gatingBase, ".diagnostics.json"))
      }
      if (!dir.exists(outputPath)) {
        dir.create(outputPath, recursive = TRUE)
      }
      gatingInputs <- list(
        phenotypeFile = if (!is.null(phenotypeFile) && nzchar(phenotypeFile)) phenotypeFile else NA_character_,
        BP_file       = if (!is.null(BP_file) && nzchar(as.character(BP_file)[1L])) BP_file else NA_character_,
        TR_mapping_file = if (!is.null(TR_mapping_file) && nzchar(as.character(TR_mapping_file)[1L])) TR_mapping_file else NA_character_
      )
      write_gating_diagnostics_json(tppGate, gatingDiagFile, inputs = gatingInputs)
      stop(
        sprintf(
          "TPP population gate failed (Requirement %s); no RData written. %s See diagnostics: %s",
          tppGate$requirement, tppGate$reason, gatingDiagFile
        ),
        call. = FALSE
      )
    }
  }

  if(havePheno) {
    analysisIdPheno <- round(as.numeric(Sys.time()), 0)
    cat(paste0("Phenotypic QA/QC analysisId: ", analysisIdPheno, "\n"))

    # --- data -> pheno
    data_pheno <- read.csv(phenotypeFile, encoding = 'utf-8', check.names = F)

    missingTraits <- setdiff(traits, colnames(data_pheno))
    if (length(missingTraits) > 0) {
      stop(paste0("These traits are not columns of the phenotype file: ",
                  paste(missingTraits, collapse = ", ")), call. = FALSE)
    }

    # remove non-alphanumeric character in occurrenceName and revise how enviroment is created
    data_pheno$environment <- paste0("env",
                                     paste(data_pheno$breedingStage,
                                           data_pheno$year,
                                           data_pheno$season,
                                           gsub("[^a-zA-Z0-9]", "", data_pheno$site),
                                           gsub("[^a-zA-Z0-9]", "", data_pheno$experimentName),
                                           gsub("[^a-zA-Z0-9]", "", data_pheno$occurrenceName), sep = "_"))

    # Create dummy columns to mimic STA behiavor
    data_pheno[,paste0(traits, "-residual")] <- NA

    # check if variable are present
    if (any(is.na(data_pheno$design)) || any(data_pheno$design == "")) {
      stop("Design information are missing for some or all rows.")
    }

    # check design parameter base, rep in EBS pertains to the number of times that entry appear in the occurrence
    data_pheno[data_pheno$design == "Partially Replicated", "rep"] <- NA
    data_pheno[data_pheno$design == "Augmented", "rep"] <- NA
    data_pheno[data_pheno$design == "Augmented RCBD", "rep"] <- NA # assumes that blockNumber is not NA

    # --- metadata -> pheno
    # NOTE: rep() is bugged with python
    # metadata_pheno <- data.frame(
    #   parameter = as.vector(c("stage", "year", "season", "location", "trial", "study", "rep", "iBlock", "row", "col", "designation", "gid", "entryType", rep("trait", nTrait), "environment", "source")),
    #   value = as.vector(c("breedingStage", "year", "season", "site", "experiment", "occurrenceName", "rep", "blockNumber", "paY", "paX", "germplasmName", "germplasmDbId", "entryType", traits, "environment", "ebs-ba"))
    # )

    metadata_pheno1 <- data.frame(
      parameter = c("stage", "year", "season", "location", "trial", "study", "rep", "iBlock", "row", "col", "designation", "gid", "entryType"),
      value = c("breedingStage", "year", "season", "site", "experimentName", "occurrenceName", "rep", "blockNumber", "paY", "paX", "germplasmName", "germplasmDbId", "entryType")
    )
    metadata_pheno2 <- data.frame(
      parameter = "trait",
      value = traits
    )
    if (is.null(requestId)) {
      metadata_pheno3 <- data.frame(
        parameter = c("environment", "source", "sourceId"),
        value = c("environment", "ebs-ba", NA)
      )
    } else {
      metadata_pheno3 <- data.frame(
        parameter = c("environment", "source", "sourceId"),
        value = c("environment", "ebs-ba", requestId)
      )
    }

    metadata_pheno <- rbind(rbind(metadata_pheno1, metadata_pheno2), metadata_pheno3)

    metadata_pheno_parameter_size <- length(metadata_pheno$parameter)
    metadata_pheno_value_size <- length(metadata_pheno$value)
    if (metadata_pheno_parameter_size != metadata_pheno_value_size) {
      stop(paste0("metadata length mismatch: ",
                  "parameter=", metadata_pheno_parameter_size,
                  " vs. value=", metadata_pheno_value_size))
    }

    # --- modifications -> pheno
    modifications_pheno <- data.frame(
      module = "qaRaw",
      analysisId = analysisIdPheno,
      trait = traits,
      reason = "none",
      row = NA,
      value = NA
    )
    row.names(modifications_pheno) <- traits

    # --- status
    # One row per analysis actually performed. analysisIdName must be a real
    # string: MTA and the qa/qc modules label their dropdowns with it.
    status <- data.frame(
      module = "qaRaw",
      analysisId = analysisIdPheno,
      analysisIdName = "qa_ebs_pdm",
      stringsAsFactors = FALSE
    )

    # --- modeling
    modeling <- data.frame(
      module = "qaRaw",
      analysisId = analysisIdPheno,
      trait = traits,
      environment = NA,
      parameter = "outlierCoefOutqPheno",
      value = NA
    )
    row.names(modeling) <- traits

    # --- empty pheno QA/QC stamp (the analysisId handed to resMod/staLMM)
    # EBS performs the experimental-design, missingness and consistency checks
    # upstream, so here we only record the stamp (a status row, no row-level
    # modifications). resMod()/staLMM() filter modifications$pheno by this id;
    # passing the qaRaw placeholder (whose rows are NA) would crash resMod on
    # `mydata[as.numeric(row), ] <- NA`, so we use a dedicated empty stamp with
    # no modification rows. Module 'qaConsistency' is one of the stamps Bioflow's
    # model-based QA dropdown recognises.
    analysisIdPhenoQA <- round(as.numeric(Sys.time()), 0) + 2
    status <- rbind(status, data.frame(
      module = "qaConsistency",
      analysisId = analysisIdPhenoQA,
      analysisIdName = "qa_ebs_consist",
      stringsAsFactors = FALSE
    ))
  }
  # --- data -> pedigree
  # Read the pedigree file when one is supplied. Every column of the file is
  # carried through; designation, mother and father are mandatory. Falling back
  # to a placeholder keeps the schema valid but disables the two qa/qc use cases.
  if (!is.null(pedigreeFile) && nzchar(pedigreeFile)) {
    pedigreeBundle <- read_pedigree_file(
      pedigreeFile    = pedigreeFile,
      mapping         = pedigreeMapping,
      deriveCrossType = deriveCrossType,
      restrictTo      = if (havePheno) data_pheno$germplasmName else NULL
    )

  } else {
    pedigreeBundle <- synthesize_pedigree_from_pheno(data_pheno$germplasmName)
  }
  data_pedigree     <- pedigreeBundle$data
  metadata_pedigree <- pedigreeBundle$metadata

  # --- data -> geno
  # Genotypes are optional: MTA works without them, but both qa/qc use cases
  # need them (and F1 qa/qc needs the parents genotyped, not just the progeny).

  data_geno              <- NULL
  data_geno_imp          <- NULL
  metadata_geno          <- NULL
  modifications_geno     <- NULL
  modifications_geno_imp <- NULL
  analysisIdGeno         <- NULL
  genoPloidy             <- NULL

  if (haveGeno) {
    analysisIdGeno <- round(as.numeric(Sys.time()), 0) + 1
    cat(paste0("Genotype QA/QC analysisId: ", analysisIdGeno, "\n"))


    data_geno <- read_vcf(genotypeFile)

    # EBS VCFs key samples by sample_id (sample_code), but F1 qa/qc matches
    # parents by designation and Pedigree qa/qc matches designation/mother/father
    # triplets by designation. Translate parental sample_ids to their designation
    # (keeping F1 progeny under their sample_id) so both qa/qc modules can match
    # the parents. No-op when the pedigree carries no sample_id mapping or the
    # VCF is already designation-keyed.
    data_geno <- align_geno_to_pedigree(data_geno, data_pedigree)

    genoPloidy <- max(adegenet::ploidy(data_geno))

    filter_1 <- list("maf", ">=", 0)

    filtering_steps <- list(
      filter_1
    )

    filt_seq <- lapply(filtering_steps, function(x){
      setNames(as.list(x), c("param", "operator", "threshold"))
    })

    filt_gl <- apply_sequence_filtering(data_geno, filt_seq)

    # --- metadata -> geno
    metadata_geno <- data.frame(
      parameter = as.vector(c("input_format", "ploidity")),
      value = as.vector(c("vcf", genoPloidy))
    )

    # --- modifications -> geno
    filt_mods <- get_filter_log(data_geno, filt_gl$filt_log)

    if(dim(filt_mods)[1] > 0) {
      modifications_geno <- filt_mods
    } else {
      modifications_geno <- data.frame(reason = c(NA),
                                row = c(NA),
                                col = c(NA),
                                value = c(NA))
    }

    modifications_geno$analysisId <- analysisIdGeno
    modifications_geno$analysisIdName <- "qa_ebs_mda"
    modifications_geno$module <- "qaGeno"

    # --- data -> geno_imp
    # Keyed by the qaGeno analysisId as a character string: every consumer
    # resolves it with which(names(geno_imp) == <stamp>).
    imp_freq <- impute_gl(filt_gl$gl,
                          ploidity = genoPloidy,
                          method = 'frequency')

    data_geno_imp <- list()
    data_geno_imp[[as.character(analysisIdGeno)]] <- imp_freq$gl

    # --- modifications -> geno_imp
    # impute_gl() returns list(gl = ..., log = ...); reading
    # imp_freq$imputation_log$log silently yielded NULL and discarded the log.
    # Single-bracket assignment with list() keeps the analysisId key present even
    # when the log is NULL (nothing needed imputing); `[[<-` would delete it.
    modifications_geno_imp <- list()
    modifications_geno_imp[as.character(analysisIdGeno)] <- list(imp_freq$log)
    if (is.null(imp_freq$log)) {
      cli::cli_inform(paste0(
        "No missing genotype calls to impute; modifications$geno_imp[['",
        analysisIdGeno, "']] recorded as NULL."
      ))
    }

    if (is.null(status)) {
      status <- data.frame(
        module = "qaGeno",
        analysisId = analysisIdGeno,
        analysisIdName = "qa_ebs_mda",
        stringsAsFactors = FALSE)
    } else {
      status <- rbind(status, data.frame(
        module = "qaGeno",
        analysisId = analysisIdGeno,
        analysisIdName = "qa_ebs_mda",
        stringsAsFactors = FALSE
    ))
    }

  } else {
    cli::cli_warn(paste0(
      "No genotype file supplied. F1 qa/qc and Pedigree qa/qc both require ",
      "marker data and will not be available on this object."
    ))
  }


  # --- Create final R object
  result <- list(
    data = list(
      pheno = if (havePheno) data_pheno else NULL,  # data.frame
      pedigree = data_pedigree, # data.frame
      geno = if (haveGeno) data_geno else NULL,  # genlight object
      geno_imp = if (haveGeno) data_geno_imp else NULL  # named list of genlight objects, keyed by analysisId
    ),
    metadata = list(
      pheno = if (havePheno) metadata_pheno else NULL,  # data.frame
      pedigree = metadata_pedigree,  # data.frame
      geno = if (haveGeno) metadata_geno else NULL  # data.frame
    ),
    modifications = list(
      pheno = if (havePheno) modifications_pheno else NULL, # data.frame
      geno = if (haveGeno) modifications_geno else NULL,  # data.frame
      geno_imp =  if (haveGeno) modifications_geno_imp else NULL # named list, keyed by analysisId
    ),
    status = status,  # data.frame
    modeling = if (havePheno) modeling else NULL  # data.frame
  )

  # Drop the slots we could not populate so Bioflow's is.null() guards fire
  # correctly instead of seeing an empty-but-present element.
  result$data          <- result$data[!vapply(result$data, is.null, logical(1))]
  result$metadata      <- result$metadata[!vapply(result$metadata, is.null, logical(1))]
  result$modifications <- result$modifications[!vapply(result$modifications, is.null, logical(1))]

  # --- Model-based outlier detection (qaMb), then Single Trial Analysis
  # The updated cgiarPipeline::staLMM() returns NULL unless a 'qaMb' stamp
  # exists in status, so SOA now requires model-based outlier detection to run
  # first. This mimics Bioflow's qaSta "Model-Based Outlier Detection" step:
  # cgiarPipeline::resMod() fits, per trait x environment, the SOA linear model
  # and flags a record as an outlier when >= 2 of 3 criteria hold
  # (Cook's distance > 4/N, |studentized residual| > 3.5, leverage > 2p/N).
  # EBS performs the upstream pheno QA (experimental design / missingness /
  # consistency) itself, so those stamps are not generated here; the existing
  # qaRaw stamp (analysisIdPheno) is the pheno-QA analysisId passed downstream.
  if (isTRUE(runSOA)) {

    # 1) Compute residual diagnostics and outlier flags for the mapped traits.
    resModOut <- try(
      cgiarPipeline::resMod(
        phenoDTfile = result,
        analysisId  = analysisIdPhenoQA,
        trait       = traits,
        genoUnit    = "designation",
        addRes      = FALSE,
        verbose     = FALSE
      ),
      silent = TRUE
    )

    if (inherits(resModOut, "try-error") || is.null(resModOut) || nrow(resModOut) == 0) {
      cli::cli_warn(paste0(
        "Model-based outlier detection (cgiarPipeline::resMod) produced no ",
        "result; the qaMb stamp and SOA will be skipped. Detail: ",
        if (inherits(resModOut, "try-error")) as.character(resModOut) else "empty result"
      ))
    } else {
      # 2) Write the qaMb stamp exactly as Bioflow's qaSta module does:
      #    - outlier rows go to modifications$pheno (module 'qaMb', reason
      #      'outlierRes'); cookD/leverage thresholds go to modeling; a 'qaMb'
      #      status row records the stamp.
      analysisIdMb <- as.numeric(Sys.time())
      tmpOut <- resModOut[which(resModOut$outlierFinal == "outlier"), , drop = FALSE]

      if (nrow(tmpOut) > 0) {
        modMb <- data.frame(
          module     = "qaMb",
          analysisId = analysisIdMb,
          trait      = tmpOut[, "trait"],
          reason     = "outlierRes",
          row        = tmpOut[, "idRowR"],
          value      = tmpOut[, "Observed"],
          stringsAsFactors = FALSE
        )
        if (is.null(result$modifications$pheno)) {
          result$modifications$pheno <- modMb
        } else {
          result$modifications$pheno <- rbind(
            result$modifications$pheno, modMb[, colnames(result$modifications$pheno)]
          )
        }
      }

      # Per trait x environment thresholds (Cook's D = 4/N, leverage = 2p/N).
      dfU <- resModOut[!duplicated(resModOut[, c("trait", "environment")]), , drop = FALSE]
      modelingMb <- rbind(
        data.frame(module = "qaMb", analysisId = analysisIdMb, trait = dfU$trait,
                   environment = dfU$environment, parameter = "cookD",
                   value = round(4 / dfU$N, 4), stringsAsFactors = FALSE),
        data.frame(module = "qaMb", analysisId = analysisIdMb, trait = dfU$trait,
                   environment = dfU$environment, parameter = "leverage",
                   value = round((2 * dfU$p) / dfU$N, 4), stringsAsFactors = FALSE)
      )
      if (is.null(result$modeling)) {
        result$modeling <- modelingMb
      } else {
        result$modeling <- rbind(result$modeling, modelingMb[, colnames(result$modeling)])
      }

      statusMb <- data.frame(module = "qaMb", analysisId = analysisIdMb,
                             analysisIdName = "qa_ebs_mb", stringsAsFactors = FALSE)
      result$status <- rbind(result$status, statusMb[, colnames(result$status)])

      cli::cli_inform(paste0(
        "Model-based outlier detection tagged ", nrow(tmpOut),
        " record(s) as outliers across ", nrow(dfU),
        " trait x environment combination(s); qaMb stamp ", analysisIdMb, " added."
      ))
    }

    # 3) Single Trial Analysis (SOA). staLMM() requires a qaMb stamp and returns
    #    NULL otherwise, so only run it when the qaMb step above succeeded, and
    #    never overwrite a good object with a NULL/failed result.
    haveQaMb <- !is.null(result$status) && "qaMb" %in% result$status$module
    if (haveQaMb) {
      resSta <- try(
        cgiarPipeline::staLMM(phenoDTfile = result, analysisId=analysisIdPhenoQA,
                                    trait=traits,
                                    traitFamily = NULL,
                                    fixedTerm = NULL,
                                    returnFixedGeno = "TRUE",
                                    genoUnit = "designation",
                                    rowColRole = "spatial",
                                    verbose = "FALSE",
                                    maxit = 35),
        silent = TRUE
      )
      if (inherits(resSta, "try-error") || is.null(resSta)) {
        cli::cli_warn(paste0(
          "Single Trial Analysis (cgiarPipeline::staLMM) did not return a result; ",
          "the object is kept without an SOA stamp. Detail: ",
          if (inherits(resSta, "try-error")) as.character(resSta) else "staLMM returned NULL"
        ))
      } else {
        result <- resSta
      }
    } else {
      cli::cli_warn(paste0(
        "No qaMb stamp present (model-based outlier detection did not run), ",
        "so Single Trial Analysis was skipped."
      ))
    }

    # staLMM adds its own "sta" status row but leaves analysisIdName as NA, which
    # renders as "NA_<timestamp>" in the MTA stamp dropdown. Name it in place -
    # appending another row instead would create a duplicate "sta" entry and
    # break that dropdown (names(traitsMta) <- paste(...) length mismatch).
    if (!is.null(result$status) && any(c("sta", "soa") %in% result$status$module)) {
      if (!"analysisIdName" %in% colnames(result$status)) {
        result$status$analysisIdName <- NA_character_
      }
      result$status$analysisIdName <- as.character(result$status$analysisIdName)
      staRows <- which(result$status$module %in% c("sta", "soa"))
      needsName <- staRows[is.na(result$status$analysisIdName[staRows]) |
                             !nzchar(result$status$analysisIdName[staRows])]
      if (length(needsName) > 0) {
        result$status$analysisIdName[needsName] <- "ebs_soa_ph"
      }
    }

    # staLMM returns a 13-column predictions table without effectType. MTA has a
    # shim that adds it, but other consumers (and the canonical schema) expect
    # it, so add it here rather than relying on that shim.
    if (!is.null(result$predictions) && !"effectType" %in% colnames(result$predictions)) {
      result$predictions$effectType <- NA_character_
    }
  }

  # --------------------------------------------------------------------------
  # TPP population hand-off (Requirement 3.7 -> Requirements 4 & 5).
  #
  # Reached only when the opt-in gate passed with populate = TRUE (useBP = TRUE,
  # usable bpid, and both BP_file + TR_mapping_file readable). At this point the
  # `result` object is fully assembled, so this is where automatic TPP
  # population plugs in. Tasks 6 and 7 implement the actual build of
  # `result$data$TPP` / `result$metadata$TPP`:
  #   * Task 6 - full path  (a record carries BOTH bpid AND msid)
  #   * Task 7 - bpid-only path (msid empty/absent across rows) + skip path
  #
  # Everything the build step needs is prepared here and left in `tppPopInputs`:
  #   - the loaded TPP file (BP_file) via read_tpp_file()
  #   - the loaded TR mapping via read_tr_mapping() (resolved path from the gate)
  #   - the BP-matched TPP rows for this phenotype's bpid via match_bp()
  # These reuse the Task 4 ingestion helpers (scripts/tpp_population_helpers.R).
  tppPopInputs <- NULL
  if (!is.null(tppGate) && isTRUE(tppGate$populate)) {
    tpp_data    <- read_tpp_file(tppGate$bp_file)
    tr_mapping  <- read_tr_mapping(tppGate$tr_mapping_path)

    # The gate guaranteed usable bpid; take the first non-NA bpid as the key.
    pheno_bpid  <- if (havePheno && "bpid" %in% names(data_pheno)) data_pheno$bpid else NULL
    matched_bp  <- match_bp(tpp_data, pheno_bpid)

    tppPopInputs <- list(
      tpp_data   = tpp_data,    # full TPP file (character columns)
      tr_mapping = tr_mapping,  # resolver list from read_tr_mapping()
      matched_bp = matched_bp,  # TPP rows whose normalized BPId matches the pheno bpid
      bp_file    = tppGate$bp_file,
      tr_mapping_path = tppGate$tr_mapping_path
    )

    cli::cli_inform(paste0(
      "useBP gate passed: loaded TPP file (", nrow(tpp_data), " rows), ",
      "TR mapping (", nrow(tr_mapping$mapping), " entries), matched ",
      nrow(matched_bp), " TPP row(s) for the phenotype bpid."
    ))

    # >>> POPULATION HAND-OFF POINT <<<
    # Build result$data$TPP[[tpp_id]] + result$metadata$TPP[[tpp_id]] here,
    # BEFORE the RData is saved below, reusing the Task 6 build helpers
    # (scripts/tpp_build_helpers.R). Diagnostics from skip branches are
    # accumulated in `tppPopDiagnostics` and surfaced on the console; the
    # diagnostics JSON itself is written by validate_bioflow_object() below.
    tppPopDiagnostics <- list()

    # Pheno designations drive the exact-match check auto-mapping (Req 4.7/4.8).
    pheno_designations <- if (havePheno && "germplasmName" %in% names(data_pheno)) {
      unique(as.character(data_pheno$germplasmName))
    } else {
      character(0)
    }

    # Distinguish records carrying a non-empty msid (full path, Task 6 /
    # Requirement 4) from the msid-empty case (bpid-only path, Task 7 /
    # Requirement 5). A record is on the full path when it has a usable bpid AND
    # a non-empty msid. The TPP file is tidy on traits, so we drive population by
    # the DISTINCT non-empty msid values present in the pheno data.
    pheno_msid <- if (havePheno && "msid" %in% names(data_pheno)) {
      trimws(as.character(data_pheno$msid))
    } else {
      character(0)
    }
    full_path_msids <- unique(pheno_msid[!is.na(pheno_msid) & nzchar(pheno_msid)])

    if (length(full_path_msids) > 0L) {
      # ---- FULL PATH (Requirement 4): both bpid and msid present -------------
      for (ms in full_path_msids) {
        built <- populate_tpp_full_path(
          tpp_data           = tpp_data,
          tr_mapping         = tr_mapping,
          msid               = ms,
          pheno_designations = pheno_designations
        )

        if (identical(built$status, "populated")) {
          if (is.null(result$data))     result$data     <- list()
          if (is.null(result$metadata)) result$metadata <- list()
          if (is.null(result$data$TPP))     result$data$TPP     <- list()
          if (is.null(result$metadata$TPP)) result$metadata$TPP <- list()

          result$data$TPP[[built$tpp_id]]     <- built$data
          result$metadata$TPP[[built$tpp_id]] <- built$metadata

          cli::cli_inform(paste0(
            "TPP population (full path): msid '", ms, "' -> ", built$tpp_id,
            " (", nrow(built$metadata$traits), " mapped trait(s))."
          ))
        } else {
          # Requirement 4.2 / 4.3: skip, leave the object unchanged for this
          # record, and record the diagnostic. (status "no_match" or "not_1to1")
          tppPopDiagnostics[[length(tppPopDiagnostics) + 1L]] <- list(
            msid       = ms,
            status     = built$status,
            diagnostic = built$diagnostic
          )
          cli::cli_warn(built$diagnostic)
        }
      }
    } else {
      # ---- bpid-only path (Requirement 5): msid empty/absent across rows -----
      # No single market segment to locate one TPP by, so populate EVERY distinct
      # TPP under the matched BP, retaining ALL environments for every trait (no
      # subsetting -> build_tpp_metadata emits env_filters = NULL). The gate
      # guaranteed usable bpid; `pheno_bpid` is the pheno bpid column.
      built <- populate_tpp_bpid_only(
        tpp_data           = tpp_data,
        tr_mapping         = tr_mapping,
        bpid               = pheno_bpid,
        pheno_designations = pheno_designations
      )

      if (identical(built$status, "populated")) {
        if (is.null(result$data))     result$data     <- list()
        if (is.null(result$metadata)) result$metadata <- list()
        if (is.null(result$data$TPP))     result$data$TPP     <- list()
        if (is.null(result$metadata$TPP)) result$metadata$TPP <- list()

        for (tpp_code in names(built$tpps)) {
          entry <- built$tpps[[tpp_code]]
          result$data$TPP[[entry$tpp_id]]     <- entry$data
          result$metadata$TPP[[entry$tpp_id]] <- entry$metadata

          cli::cli_inform(paste0(
            "TPP population (bpid-only): ", entry$tpp_id, " (",
            nrow(entry$metadata$traits), " mapped trait(s), all environments retained)."
          ))
        }

        cli::cli_inform(paste0(
          "TPP population (bpid-only): populated ", length(built$tpp_codes),
          " distinct TPP(s) under the matched BP (", 
          paste(built$tpp_codes, collapse = ", "), ")."
        ))
      } else {
        # Requirement 5.2 / 5.4: no BP rows -> skip, generate the RData WITHOUT
        # TPP population, record the diagnostic. status "no_bp_match".
        tppPopDiagnostics[[length(tppPopDiagnostics) + 1L]] <- list(
          msid       = NA_character_,
          status     = built$status,
          diagnostic = built$diagnostic
        )
        cli::cli_warn(built$diagnostic)
      }
    }

    # Retain the accumulated diagnostics for downstream reporting / tests.
    attr(tppPopInputs, "tpp_population_diagnostics") <- tppPopDiagnostics
  }

  # check if folder exist or not
  if (!dir.exists(outputPath)) {
    dir.create(outputPath, recursive = TRUE)
  }

  # Honour the requested file name; fall back to the md5 convention when none
  # is given (previously the argument was accepted and then always overwritten).
  if (is.null(outputFile) || !nzchar(outputFile)) {
    if (!is.null(analysisIdPheno)){
      outputFile <- openssl::md5(as.character(analysisIdPheno))
    } else {
      outputFile <- openssl::md5(as.character(analysisIdGeno))
    }
  }
  outFilePath <- file.path(outputPath, paste0(outputFile, ".RData"))
  save(result, file = outFilePath)
  cat(paste0("Saved: ", outFilePath, "\n"))

  # Diagnostics JSON: default beside the .RData as <outputFile>.diagnostics.json.
  if (is.null(diagnosticsFile)) {
    diagnosticsFile <- file.path(outputPath, paste0(outputFile, ".diagnostics.json"))
  }
  diagInputs <- list(
    phenotypeFile = if (!is.null(phenotypeFile) && nzchar(phenotypeFile)) phenotypeFile else NA_character_,
    pedigreeFile  = if (!is.null(pedigreeFile)  && nzchar(pedigreeFile))  pedigreeFile  else NA_character_,
    genotypeFile  = if (!is.null(genotypeFile)  && nzchar(genotypeFile))  genotypeFile  else NA_character_,
    outputFile    = outFilePath
  )

  if (isTRUE(validate)) {
    validate_bioflow_object(result, diagnosticsFile = diagnosticsFile, inputs = diagInputs,
                            runExecution = runExecution)
  } else {
    # Even when the console report is suppressed, still emit the JSON so the
    # EBS engineer always gets a machine-readable diagnostics artifact.
    validate_bioflow_object(result, verbose = FALSE,
                            diagnosticsFile = diagnosticsFile, inputs = diagInputs,
                            runExecution = runExecution)
  }

  invisible(list(result = result, file = outFilePath, diagnosticsFile = diagnosticsFile))

} # end of staRData fxn
