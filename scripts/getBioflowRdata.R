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
#' @param runSTA logical; run `cgiarPipeline::staLMM()` so the object arrives MTA-ready. Default TRUE. Set FALSE to export a pedigree/genotype-only object for the qa/qc use cases.
#' @param validate logical; print the use-case readiness report. Default TRUE.
#'
#' @return Invisibly, a list with `result` (the object written) and `file` (the
#'   path it was written to).
# -------------------------------------------------------------------------------------

getBioflowRData <- function(phenotypeFile, pedigreeFile = NULL, genotypeFile = NULL,
                     traits, outputPath, outputFile = NULL, requestId = NULL,
                     pedigreeMapping = NULL, deriveCrossType = TRUE,
                     runSTA = TRUE, validate = TRUE) {

  # Now we can use ensure_cran_packages since packages_verification.R has been sourced
  ensure_cran_packages(c("vcfR", "adegenet", "cli", "rlang", "remotes",
                         "openssl", "stringr", "purrr", "glue", "Matrix"))
  # Install if not cgiarPipelines and cgiarBase using remotes
  ensure_github_packages()

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

  # --- data -> pedigree
  # Read the pedigree file when one is supplied. Every column of the file is
  # carried through; designation, mother and father are mandatory. Falling back
  # to a placeholder keeps the schema valid but disables the two qa/qc use cases.
  if (!is.null(pedigreeFile) && nzchar(pedigreeFile)) {
    pedigreeBundle <- read_pedigree_file(
      pedigreeFile    = pedigreeFile,
      mapping         = pedigreeMapping,
      deriveCrossType = deriveCrossType,
      restrictTo      = data_pheno$germplasmName
    )
  } else {
    pedigreeBundle <- synthesize_pedigree_from_pheno(data_pheno$germplasmName)
  }
  data_pedigree     <- pedigreeBundle$data
  metadata_pedigree <- pedigreeBundle$metadata

  # --- data -> geno
  # Genotypes are optional: MTA works without them, but both qa/qc use cases
  # need them (and F1 qa/qc needs the parents genotyped, not just the progeny).
  haveGeno <- !is.null(genotypeFile) && nzchar(genotypeFile)

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
  } else {
    cli::cli_warn(paste0(
      "No genotype file supplied. F1 qa/qc and Pedigree qa/qc both require ",
      "marker data and will not be available on this object."
    ))
  }

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
  if (haveGeno) {
    status <- rbind(status, data.frame(
      module = "qaGeno",
      analysisId = analysisIdGeno,
      analysisIdName = "qa_ebs_mda",
      stringsAsFactors = FALSE
    ))
  }

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

  # --- Create final R object
  result <- list(
    data = list(
      pheno = data_pheno,  # data.frame
      pedigree = data_pedigree, # data.frame
      geno = data_geno,  # genlight object
      geno_imp = data_geno_imp  # named list of genlight objects, keyed by analysisId
    ),
    metadata = list(
      pheno = metadata_pheno,  # data.frame
      pedigree = metadata_pedigree,  # data.frame
      geno = metadata_geno  # data.frame
    ),
    modifications = list(
      pheno = modifications_pheno, # data.frame
      geno = modifications_geno,  # data.frame
      geno_imp = modifications_geno_imp  # named list, keyed by analysisId
    ),
    status = status,  # data.frame
    modeling = modeling  # data.frame
  )

  # Drop the slots we could not populate so Bioflow's is.null() guards fire
  # correctly instead of seeing an empty-but-present element.
  result$data          <- result$data[!vapply(result$data, is.null, logical(1))]
  result$metadata      <- result$metadata[!vapply(result$metadata, is.null, logical(1))]
  result$modifications <- result$modifications[!vapply(result$modifications, is.null, logical(1))]

  # --- Single Trial Analysis, so the object arrives MTA-ready
  if (isTRUE(runSTA)) {
    result <- cgiarPipeline::staLMM(phenoDTfile = result, analysisId=analysisIdPheno,
                                    trait=traits,
                                    traitFamily = NULL,
                                    fixedTerm = NULL,
                                    returnFixedGeno = "TRUE",
                                    genoUnit = "designation",
                                    rowColRole = "spatial",
                                    verbose = "FALSE",
                                    maxit = 35)

    # staLMM adds its own "sta" status row but leaves analysisIdName as NA, which
    # renders as "NA_<timestamp>" in the MTA stamp dropdown. Name it in place -
    # appending another row instead would create a duplicate "sta" entry and
    # break that dropdown (names(traitsMta) <- paste(...) length mismatch).
    if (!is.null(result$status) && "sta" %in% result$status$module) {
      if (!"analysisIdName" %in% colnames(result$status)) {
        result$status$analysisIdName <- NA_character_
      }
      result$status$analysisIdName <- as.character(result$status$analysisIdName)
      staRows <- which(result$status$module == "sta")
      needsName <- staRows[is.na(result$status$analysisIdName[staRows]) |
                             !nzchar(result$status$analysisIdName[staRows])]
      if (length(needsName) > 0) {
        result$status$analysisIdName[needsName] <- "ebs_sta_ph"
      }
    }

    # staLMM returns a 13-column predictions table without effectType. MTA has a
    # shim that adds it, but other consumers (and the canonical schema) expect
    # it, so add it here rather than relying on that shim.
    if (!is.null(result$predictions) && !"effectType" %in% colnames(result$predictions)) {
      result$predictions$effectType <- NA_character_
    }
  }

  # check if folder exist or not
  if (!dir.exists(outputPath)) {
    dir.create(outputPath, recursive = TRUE)
  }

  # Honour the requested file name; fall back to the md5 convention when none
  # is given (previously the argument was accepted and then always overwritten).
  if (is.null(outputFile) || !nzchar(outputFile)) {
    outputFile <- openssl::md5(as.character(analysisIdPheno))
  }
  outFilePath <- file.path(outputPath, paste0(outputFile, ".RData"))
  save(result, file = outFilePath)
  cat(paste0("Saved: ", outFilePath, "\n"))

  if (isTRUE(validate)) {
    validate_bioflow_object(result)
  }

  invisible(list(result = result, file = outFilePath))

} # end of staRData fxn
