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
# Name             : read_pedigree_functions
# Description      : Read an EBS pedigree csv file and turn it into the
#                    data$pedigree / metadata$pedigree pair expected by Bioflow.
# Maintainer       : Bioflow core development team
# ---------------------------------------------------------------------------------------------
#
# WHY THIS FILE EXISTS
# --------------------
# Three Bioflow use cases need pedigree information:
#
#   1. MTA            - optional (only the SCA/GCA models need mother/father)
#   2. F1 qa/qc       - REQUIRED: designation, mother, father, crossType
#                       (+ sample_id when several genotyped samples share a designation)
#   3. Pedigree qa/qc - REQUIRED: designation, mother, father
#
# CRITICAL COMPATIBILITY CONSTRAINT
# ---------------------------------
# bioflow/R/mod_qaPedApp.R renames the pedigree columns *positionally*:
#
#     uno = data.frame(value = names(peddata))
#     names(peddata) = metaped$parameter[which(metaped$value %in% uno$value == T)]
#
# The right-hand side is evaluated in `metadata$pedigree` ROW order, and is then
# assigned to `names(peddata)` in `data$pedigree` COLUMN order. There is no
# name-wise matching. Consequently, for the Pedigree qa/qc module to read the
# right columns, the object MUST satisfy:
#
#   (a) every column of data$pedigree has exactly one row in metadata$pedigree
#       whose `value` equals that column name; and
#   (b) those metadata rows appear in the same order as the columns.
#
# If (a) or (b) is violated, qaPed does not error out - it silently mislabels
# columns and reports nonsense. Everything below is built to guarantee the
# invariant by construction, and `assert_pedigree_alignment()` re-checks it.
#
# The F1 qa/qc path (cgiarPipeline::individualVerification) is name-based via
# cgiarBase::replaceValues() and so is insensitive to order, but it is satisfied
# by the same structure.
# ---------------------------------------------------------------------------------------------

# Canonical Bioflow pedigree parameters, in the order we emit them.
# Mirrors the vocabulary offered by bioflow/R/mod_getDataPed.R.
PED_CANONICAL_PARAMS <- c(
  "designation",
  "sample_id",
  "mother",
  "father",
  "crossType",
  "yearOfOrigin",
  "other",
  "batch",
  "plate",
  "position"
)

# Parameters without which no pedigree-dependent module can run.
PED_REQUIRED_PARAMS <- c("designation", "mother", "father")

# Conservative header synonyms. Keys are canonical parameters, values are
# normalised (lowercased, non-alphanumeric stripped) header spellings.
# `entry` is deliberately NOT listed for crossType: it is far too generic to
# claim by name alone and is instead picked up by content sniffing below.
PED_HEADER_SYNONYMS <- list(
  designation  = c("designation", "germplasmname", "germplasm", "germplasmcode",
                   "line", "linename", "genotype", "genotypename", "entryname",
                   "accession", "accessionname"),
  sample_id    = c("sampleid", "sample", "samplename", "samplecode",
                   "dnasampleid", "sampledbid", "tissuesampleid", "sampleno"),
  mother       = c("mother", "mothername", "female", "femaleparent", "femalename",
                   "parent1", "parenta", "hybridparent1", "seedparent", "p1"),
  father       = c("father", "fathername", "male", "maleparent", "malename",
                   "parent2", "parentb", "hybridparent2", "pollenparent", "p2"),
  crossType    = c("crosstype", "crossclass", "crosscategory", "individualtype",
                   "germplasmtype", "materialtype", "pedigreetype"),
  yearOfOrigin = c("yearoforigin", "year", "originyear", "yob", "cycleyear"),
  other        = c("other", "plantnumber", "plantno", "plantnr", "othermetadata"),
  batch        = c("batch", "batchnumber", "batchno", "batchid"),
  plate        = c("plate", "platenumber", "plateno", "plateid"),
  position     = c("position", "well", "wellposition", "platepos", "positionid")
)

# Tokens recognised as "this individual is an F1 to be verified".
PED_F1_TOKENS     <- c("f1", "f_1", "hybrid", "cross", "progeny")
# Tokens recognised as "this individual is a parent / not an F1".
PED_PARENT_TOKENS <- c("parent", "p", "line", "inbred", "founder", "male", "female")

#' Normalise a header for synonym matching
#'
#' Lowercases and removes everything that is not a letter or digit so that
#' `Hybrid_Parent1`, `hybrid parent 1` and `HybridParent1` all collapse to the
#' same key.
#'
#' @param x Character vector of column names.
#' @return Character vector of normalised names.
ped_normalise_header <- function(x) {
  gsub("[^a-z0-9]", "", tolower(as.character(x)))
}

#' Blank-to-NA coercion for identifier columns
#'
#' EBS/CSV exports represent unknown parents in several ways. Bioflow expects a
#' real `NA`, and `cgiarPipeline::individualVerification` compares with `==`,
#' so empty strings and the literal text "NA" must not survive.
#'
#' @param x Vector to clean.
#' @return Character vector with blanks/placeholders as NA.
ped_blank_to_na <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "na", "N/A", "n/a", "NULL", "null", "-", ".", "?")] <- NA_character_
  x
}

#' Detect which canonical parameter a column corresponds to
#'
#' Resolution order: explicit user mapping, then header synonyms, then (for
#' crossType only) content sniffing.
#'
#' @param headers Character vector of the pedigree file's column names.
#' @param pedigreeData The pedigree data.frame, used for content sniffing.
#' @param mapping Optional named list/vector, canonical parameter -> column name,
#'   supplied by the caller to override or complete auto-detection.
#' @return Named character vector: names are canonical parameters, values are the
#'   matched column names. Only resolved parameters are present.
detect_pedigree_columns <- function(headers, pedigreeData = NULL, mapping = NULL) {
  resolved <- character(0)
  taken    <- character(0)

  # ---- 1. explicit caller-supplied mapping wins -----------------------------
  if (!is.null(mapping) && length(mapping) > 0) {
    mapping <- unlist(mapping)
    unknownParams <- setdiff(names(mapping), PED_CANONICAL_PARAMS)
    if (length(unknownParams) > 0) {
      cli::cli_abort(paste0(
        "`pedigreeMapping` contains parameters Bioflow does not recognise: ",
        paste(unknownParams, collapse = ", "), ". Valid values are: ",
        paste(PED_CANONICAL_PARAMS, collapse = ", "), "."
      ))
    }
    for (param in names(mapping)) {
      col <- mapping[[param]]
      if (is.na(col) || !nzchar(col)) next
      if (!col %in% headers) {
        cli::cli_abort(paste0(
          "`pedigreeMapping` maps ", param, " to column '", col,
          "' which is not present in the pedigree file. Available columns: ",
          paste(headers, collapse = ", "), "."
        ))
      }
      resolved[param] <- col
      taken <- c(taken, col)
    }
  }

  # ---- 2. header synonyms ---------------------------------------------------
  normHeaders <- ped_normalise_header(headers)
  for (param in PED_CANONICAL_PARAMS) {
    if (param %in% names(resolved)) next
    hits <- which(normHeaders %in% PED_HEADER_SYNONYMS[[param]] & !(headers %in% taken))
    if (length(hits) == 0) next
    if (length(hits) > 1) {
      # Prefer the synonym listed earliest (most canonical spelling).
      rank <- match(normHeaders[hits], PED_HEADER_SYNONYMS[[param]])
      hits <- hits[order(rank)]
      cli::cli_inform(paste0(
        "Several columns could be '", param, "' (",
        paste(headers[hits], collapse = ", "), "); using '", headers[hits[1]],
        "'. Pass `pedigreeMapping` to choose explicitly."
      ))
    }
    resolved[param] <- headers[hits[1]]
    taken <- c(taken, headers[hits[1]])
  }

  # ---- 3. content sniffing for crossType ------------------------------------
  # Needed because EBS exports often call this column `entry`, which is too
  # generic to map by name. A column qualifies only if it actually contains an
  # F1 token, which is exactly what the F1 module requires.
  if (!("crossType" %in% names(resolved)) && !is.null(pedigreeData)) {
    for (col in setdiff(headers, taken)) {
      vals <- ped_blank_to_na(pedigreeData[[col]])
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0) next
      normVals <- ped_normalise_header(vals)
      uniqueVals <- unique(normVals)
      # Must contain an F1 marker and be a small controlled vocabulary.
      if (any(uniqueVals %in% PED_F1_TOKENS) &&
          all(uniqueVals %in% c(PED_F1_TOKENS, PED_PARENT_TOKENS)) &&
          length(uniqueVals) <= 5) {
        resolved["crossType"] <- col
        taken <- c(taken, col)
        cli::cli_inform(paste0(
          "Using column '", col, "' as crossType (its values look like cross types: ",
          paste(utils::head(unique(vals), 4), collapse = ", "), ")."
        ))
        break
      }
    }
  }

  resolved
}

#' Normalise crossType values to what Bioflow compares against
#'
#' Both `cgiarPipeline::individualVerification` and `mod_hybridityApp.R` test
#' `crossType == "F1"` with that exact spelling, and they use `==` rather than
#' `which()`, so an NA produces phantom all-NA rows that later trip the
#' "Designation column in the pedigree file cannot have missing data" error.
#' Everything that is not an F1 therefore becomes a non-NA label.
#'
#' @param x Raw crossType column.
#' @return Character vector using the literal "F1" for progeny.
normalise_cross_type <- function(x) {
  raw  <- ped_blank_to_na(x)
  norm <- ped_normalise_header(raw)

  out <- raw
  out[norm %in% PED_F1_TOKENS] <- "F1"
  # Any non-F1 label is passed through untouched (Bioflow never compares it),
  # but missing values must be filled or the F1 subset breaks.
  out[is.na(out)] <- "parent"
  out
}

#' Derive a crossType column when the pedigree file has none
#'
#' An individual with both parents recorded is treated as an F1 to be verified;
#' everything else is a parent. This is a heuristic, so it is always logged.
#'
#' @param mother,father Cleaned parent columns.
#' @return Character vector of "F1"/"parent".
derive_cross_type <- function(mother, father) {
  isF1 <- !is.na(mother) & !is.na(father)
  ifelse(isF1, "F1", "parent")
}

#' Assert the metadata/column alignment invariant required by qaPed
#'
#' See the header comment: Pedigree qa/qc renames columns positionally, so a
#' mismatch corrupts results silently. Fail loudly here instead.
#'
#' @param pedigreeData data$pedigree.
#' @param pedigreeMetadata metadata$pedigree.
#' @return TRUE invisibly, or aborts.
assert_pedigree_alignment <- function(pedigreeData, pedigreeMetadata) {
  cols <- colnames(pedigreeData)
  # Reproduce mod_qaPedApp.R's expression exactly.
  renamed <- pedigreeMetadata$parameter[which(pedigreeMetadata$value %in% cols)]

  if (length(renamed) != length(cols)) {
    cli::cli_abort(paste0(
      "metadata$pedigree does not describe data$pedigree one-to-one: ",
      length(cols), " column(s) but ", length(renamed),
      " matching metadata row(s). Bioflow's Pedigree qa/qc module renames ",
      "pedigree columns positionally and would mislabel them."
    ))
  }
  mappedInOrder <- pedigreeMetadata$value[which(pedigreeMetadata$value %in% cols)]
  if (!identical(mappedInOrder, cols)) {
    cli::cli_abort(paste0(
      "metadata$pedigree row order does not match data$pedigree column order.\n",
      "  columns : ", paste(cols, collapse = ", "), "\n",
      "  metadata: ", paste(mappedInOrder, collapse = ", "), "\n",
      "Bioflow's Pedigree qa/qc module renames pedigree columns positionally, ",
      "so this would silently mislabel them."
    ))
  }
  invisible(TRUE)
}

#' Read an EBS pedigree csv into Bioflow's data/metadata pair
#'
#' Every column of the input file is carried into `data$pedigree`. Columns that
#' map onto a canonical Bioflow parameter are emitted first (in canonical
#' order); any remaining columns are appended and mapped to themselves so that
#' no information is lost and the positional-rename invariant still holds.
#'
#' @param pedigreeFile Path to the pedigree csv.
#' @param mapping Optional named list, canonical parameter -> column name, to
#'   override auto-detection (e.g. `list(crossType = "entry")`).
#' @param deriveCrossType Logical. When the file has no crossType column, infer
#'   one from parent completeness so the F1 qa/qc module can run. Default TRUE.
#' @param restrictTo Optional character vector of designations (typically the
#'   germplasm present in the phenotype file). When supplied, a `note` is
#'   returned describing the overlap; rows are NOT dropped, because parents are
#'   legitimately absent from the phenotype file.
#' @return A list with `data` (data.frame), `metadata` (parameter/value
#'   data.frame) and `notes` (character vector of diagnostics).
read_pedigree_file <- function(pedigreeFile,
                               mapping = NULL,
                               deriveCrossType = TRUE,
                               restrictTo = NULL) {

  if (is.null(pedigreeFile) || !nzchar(pedigreeFile)) {
    cli::cli_abort("`pedigreeFile` is NULL or empty; nothing to read.")
  }
  if (!file.exists(pedigreeFile)) {
    cli::cli_abort("`pedigreeFile` does not exist: {pedigreeFile}")
  }

  raw <- utils::read.csv(pedigreeFile, encoding = "utf-8", check.names = FALSE,
                         stringsAsFactors = FALSE, na.strings = c("", "NA"))

  if (nrow(raw) == 0) {
    cli::cli_abort("`pedigreeFile` has no data rows: {pedigreeFile}")
  }

  headers <- colnames(raw)
  notes   <- character(0)

  resolved <- detect_pedigree_columns(headers, pedigreeData = raw, mapping = mapping)

  missingRequired <- setdiff(PED_REQUIRED_PARAMS, names(resolved))
  if (length(missingRequired) > 0) {
    cli::cli_abort(paste0(
      "The pedigree file is missing required column(s): ",
      paste(missingRequired, collapse = ", "), ".\n",
      "Found columns: ", paste(headers, collapse = ", "), ".\n",
      "Bioflow needs at least designation, mother and father. Use ",
      "`pedigreeMapping = list(", missingRequired[1], " = \"<column>\")` to map explicitly."
    ))
  }

  # ---- assemble the canonical block, in canonical order ---------------------
  outData <- list()
  outMeta <- list()

  for (param in PED_CANONICAL_PARAMS) {
    if (!(param %in% names(resolved))) next
    srcCol <- resolved[[param]]
    values <- raw[[srcCol]]

    if (param %in% c("designation", "sample_id", "mother", "father")) {
      values <- ped_blank_to_na(values)
    }
    if (param == "crossType") {
      before <- ped_blank_to_na(values)
      values <- normalise_cross_type(values)
      nFilled <- sum(is.na(before))
      if (nFilled > 0) {
        notes <- c(notes, paste0(
          "crossType was missing on ", nFilled,
          " row(s); set to 'parent' so the F1 subset does not break."
        ))
      }
    }

    # Emit the column under its canonical name. Doing so makes the object work
    # on both the renamed and the raw code paths in mod_hybridityApp.R (which
    # reads result$data$pedigree$sample_id literally in places).
    outData[[param]] <- values
    outMeta[[length(outMeta) + 1L]] <- data.frame(parameter = param, value = param,
                                                 stringsAsFactors = FALSE)
    if (!identical(srcCol, param)) {
      notes <- c(notes, paste0("Mapped '", srcCol, "' -> ", param, "."))
    }
  }

  # ---- derive crossType if absent ------------------------------------------
  if (!("crossType" %in% names(outData)) && isTRUE(deriveCrossType)) {
    derived <- derive_cross_type(outData$mother, outData$father)
    nF1 <- sum(derived == "F1")
    outData$crossType <- derived
    outMeta[[length(outMeta) + 1L]] <- data.frame(parameter = "crossType", value = "crossType",
                                                 stringsAsFactors = FALSE)
    notes <- c(notes, paste0(
      "No crossType column found; derived one from parent completeness (",
      nF1, " row(s) flagged 'F1', ", nrow(raw) - nF1, " 'parent'). ",
      "Supply a real crossType column if this is not what you want."
    ))
  }

  # reportqaPed.Rmd subsets on crossType whenever sample_id or other is mapped,
  # so those must never appear without it.
  if (any(c("sample_id", "other") %in% names(outData)) && !("crossType" %in% names(outData))) {
    cli::cli_abort(paste0(
      "The pedigree maps sample_id/other but has no crossType column. ",
      "Bioflow's Pedigree qa/qc dashboard subsets on crossType in that case ",
      "and would fail. Provide a crossType column or set deriveCrossType = TRUE."
    ))
  }

  # ---- carry every remaining column verbatim -------------------------------
  extras <- setdiff(headers, unname(resolved))
  for (col in extras) {
    # Guard against an extra column colliding with a canonical name we emitted.
    target <- col
    if (target %in% names(outData)) {
      target <- paste0(col, "_orig")
      notes <- c(notes, paste0(
        "Extra column '", col, "' renamed to '", target, "' to avoid clashing ",
        "with a mapped Bioflow parameter."
      ))
    }
    outData[[target]] <- raw[[col]]
    # Identity mapping keeps the positional-rename invariant intact and makes
    # the column visible to Bioflow's column pickers.
    outMeta[[length(outMeta) + 1L]] <- data.frame(parameter = target, value = target,
                                                 stringsAsFactors = FALSE)
  }
  if (length(extras) > 0) {
    notes <- c(notes, paste0(
      "Carried ", length(extras), " additional pedigree column(s) through to the ",
      "R object: ", paste(extras, collapse = ", "), "."
    ))
  }

  pedigreeData     <- as.data.frame(outData, stringsAsFactors = FALSE, check.names = FALSE)
  pedigreeMetadata <- do.call(rbind, outMeta)
  rownames(pedigreeData)     <- NULL
  rownames(pedigreeMetadata) <- NULL

  # yearOfOrigin is declared integer by cgiarBase::create_getData_object().
  if ("yearOfOrigin" %in% colnames(pedigreeData)) {
    suppressWarnings(
      pedigreeData$yearOfOrigin <- as.integer(ped_blank_to_na(pedigreeData$yearOfOrigin))
    )
  }

  if (any(is.na(pedigreeData$designation))) {
    nBad <- sum(is.na(pedigreeData$designation))
    cli::cli_abort(paste0(
      "The designation column has ", nBad, " missing value(s). Bioflow's F1 ",
      "qa/qc module aborts on missing designations. Please correct the file."
    ))
  }

  assert_pedigree_alignment(pedigreeData, pedigreeMetadata)

  # ---- diagnostics ---------------------------------------------------------
  if (!is.null(restrictTo)) {
    shared <- length(intersect(unique(pedigreeData$designation), unique(restrictTo)))
    notes  <- c(notes, paste0(
      "Pedigree/phenotype designation overlap: ", shared, " of ",
      length(unique(pedigreeData$designation)), " pedigree designation(s) ",
      "also appear in the phenotype file."
    ))
  }

  nWithBoth <- sum(!is.na(pedigreeData$mother) & !is.na(pedigreeData$father))
  notes <- c(notes, paste0(
    "Read ", nrow(pedigreeData), " pedigree row(s); ", nWithBoth,
    " have both parents recorded."
  ))
  if ("crossType" %in% colnames(pedigreeData)) {
    notes <- c(notes, paste0(
      "crossType == 'F1' on ", sum(pedigreeData$crossType == "F1"), " row(s)."
    ))
  }

  for (n in notes) cli::cli_inform(n)

  list(data = pedigreeData, metadata = pedigreeMetadata, notes = notes)
}

#' Build the placeholder pedigree used when no pedigree file is supplied
#'
#' Keeps the schema of `cgiarBase::create_getData_object()` so nothing
#' downstream trips on a missing column, but contains no parentage. Callers are
#' expected to tell the user that the pedigree-dependent use cases will not run.
#'
#' @param designations Character vector of germplasm names from the phenotype file.
#' @return A list with `data`, `metadata` and `notes`, as `read_pedigree_file()`.
synthesize_pedigree_from_pheno <- function(designations) {
  designations <- unique(ped_blank_to_na(designations))
  designations <- designations[!is.na(designations)]

  pedigreeData <- data.frame(
    designation  = designations,
    mother       = NA_character_,
    father       = NA_character_,
    yearOfOrigin = NA_integer_,
    stringsAsFactors = FALSE
  )
  pedigreeMetadata <- data.frame(
    parameter = c("designation", "mother", "father", "yearOfOrigin"),
    value     = c("designation", "mother", "father", "yearOfOrigin"),
    stringsAsFactors = FALSE
  )

  assert_pedigree_alignment(pedigreeData, pedigreeMetadata)

  notes <- paste0(
    "No pedigree file supplied: built a placeholder pedigree from ",
    nrow(pedigreeData), " phenotype designation(s) with empty parentage. ",
    "F1 qa/qc and Pedigree qa/qc CANNOT run on this object."
  )
  cli::cli_warn(notes)

  list(data = pedigreeData, metadata = pedigreeMetadata, notes = notes)
}
