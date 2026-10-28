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
# Name             : tpp_population_helpers
# Description      : Reusable ingestion + BP/MS matching helpers for automatic
#                    TPP population in the Bioflow RData generation script.
#                    This file covers ONLY file ingestion and ID matching
#                    (spec tpp-review Task 4); the opt-in `useBP` gating and the
#                    `getBioflowRData()` argument wiring are added by Task 5, and
#                    the actual `data$TPP` / `metadata$TPP` object build by
#                    Tasks 6/7. These helpers are written to be called with
#                    explicit file paths so Task 5 can wire them to the formal
#                    `BP_file` / `TR_mapping_file` arguments without changing
#                    their bodies.
# Maintainer       : Bioflow core development team
# ---------------------------------------------------------------------------------------------
#
# DEPENDENCY
# ----------
# BP/MS matching relies on `normalize_id()` from scripts/tpp_id_helpers.R, which
# is bundled BEFORE this file in scripts/bundle_scripts.R. Do NOT reimplement the
# normalization rule here; always call `normalize_id()`.
#
# TPP FILE (BP_file): TPP_traits_for_Bioflow.csv
# ----------------------------------------------
# Long/tidy, one row per TPP-trait. Key columns consumed by the matching logic
# (verified against the real sample):
#
#   TraitId          e.g. "TR-00388"   -> trait matching (via variable_tr_mapping)
#   TraitName        e.g. "Plant Height"
#   TraitRequirement e.g. "Essential: Maintain"
#   ScaleOption      scale bounds source
#   DesiredScore     desired-score source
#   BestCheck        best-check designation(s)
#   TPPId            e.g. "697"        -> TPP grouping key (coded form TPPCode "TPP00697")
#   TPPCode          e.g. "TPP00697"
#   BPId             e.g. "3"          -> BP matching key   (coded form BPCode "BP00003")
#   BPCode           e.g. "BP00003"
#   MSId             e.g. "54"         -> MS matching key   (coded form MSShortName "... MS00054")
#   MSShortName      e.g. "Wheat Spring LAC MS00054"
#
# MS <-> TPP is 1:1.
#
# variable_tr_mapping.csv (TR_mapping_file):
# ------------------------------------------
# Columns `ebs_variable_abbrev` -> `bp_tr_id`, e.g. WH_PH_M_CM -> TR-00388.
# Maps pheno variable abbreviations (phenotype column names) to TR_Ids so that
# a TPP row's `TraitId` can be resolved back to the pheno trait column.

# ---- Required TPP-file columns ----------------------------------------------
# The columns the ingestion + matching logic depends on. Validation fails fast
# if any are missing so Task 5's gating can surface a clear diagnostic.
.TPP_REQUIRED_COLS <- c(
  "TraitId", "TraitName", "TraitRequirement", "ScaleOption", "DesiredScore",
  "BestCheck", "TPPId", "TPPCode", "BPId", "BPCode", "MSId", "MSShortName"
)

#' Read and validate the breeding-portal TPP file (BP_file)
#'
#' Reads `TPP_traits_for_Bioflow.csv` (long/tidy, one row per TPP-trait) from an
#' explicit path and validates that the key columns required for matching are
#' present. All columns are read as character to keep IDs (e.g. "3", "54")
#' comparison-stable; matching is always performed through `normalize_id()`.
#'
#' NOTE (Task 5 wiring): the caller passes the full path. When wired to the
#' formal argument this is `BP_file`.
#'
#' @param path Full path (including filename) to the TPP CSV.
#'
#' @return A data.frame of the TPP file with all columns as character, carrying
#'   attribute `"tpp_file_path"` for diagnostics.
#' @export
read_tpp_file <- function(path) {
  if (is.null(path) || !is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("read_tpp_file: `path` must be a single non-empty file path.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("read_tpp_file: BP_file not found or unreadable: %s", path), call. = FALSE)
  }

  tpp <- utils::read.csv(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    colClasses = "character",
    check.names = TRUE,
    na.strings = c("NA", "")
  )

  missing_cols <- setdiff(.TPP_REQUIRED_COLS, names(tpp))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "read_tpp_file: TPP file is missing required column(s): %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  attr(tpp, "tpp_file_path") <- path
  tpp
}

#' Read variable_tr_mapping.csv and build a TR_Id <-> pheno-trait resolver
#'
#' The mapping file has columns `ebs_variable_abbrev` -> `bp_tr_id` (e.g.
#' WH_PH_M_CM -> TR-00388). The real sample contains trailing all-empty rows and
#' padding columns; these are dropped. The returned resolver exposes lookups in
#' both directions plus a normalized TR_Id key (via `normalize_id`) so trait
#' matching against the TPP file's `TraitId` is robust to formatting.
#'
#' NOTE (Task 5 wiring): the caller passes the resolved path. When wired to the
#' formal argument this is `file.path(dirname(phenotypeFile), TR_mapping_file)`.
#'
#' @param path Full path to variable_tr_mapping.csv.
#'
#' @return A list with:
#'   \describe{
#'     \item{mapping}{data.frame with columns `ebs_variable_abbrev`, `bp_tr_id`,
#'       `tr_key` (normalized TR_Id), non-empty rows only.}
#'     \item{tr_id_for_variable}{function(variable) -> bp_tr_id (or NA).}
#'     \item{variable_for_tr_id}{function(tr_id) -> ebs_variable_abbrev, matched
#'       via normalized TR_Id key (or NA).}
#'   }
#' @export
read_tr_mapping <- function(path) {
  if (is.null(path) || !is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("read_tr_mapping: `path` must be a single non-empty file path.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("read_tr_mapping: TR_mapping_file not found or unreadable: %s", path), call. = FALSE)
  }

  raw <- utils::read.csv(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    colClasses = "character",
    check.names = TRUE,
    na.strings = c("NA", "")
  )

  required <- c("ebs_variable_abbrev", "bp_tr_id")
  missing_cols <- setdiff(required, names(raw))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "read_tr_mapping: mapping file is missing required column(s): %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  variable <- trimws(raw[["ebs_variable_abbrev"]])
  tr_id    <- trimws(raw[["bp_tr_id"]])

  # Drop the padding / all-empty rows the export leaves behind.
  keep <- !is.na(variable) & nzchar(variable) & !is.na(tr_id) & nzchar(tr_id)
  variable <- variable[keep]
  tr_id    <- tr_id[keep]

  mapping <- data.frame(
    ebs_variable_abbrev = variable,
    bp_tr_id            = tr_id,
    tr_key              = normalize_id(tr_id),
    stringsAsFactors    = FALSE
  )

  tr_id_for_variable <- function(v) {
    v <- as.character(v)
    out <- mapping$bp_tr_id[match(v, mapping$ebs_variable_abbrev)]
    out
  }

  variable_for_tr_id <- function(id) {
    key <- normalize_id(id)
    out <- mapping$ebs_variable_abbrev[match(key, mapping$tr_key)]
    out
  }

  list(
    mapping            = mapping,
    tr_id_for_variable = tr_id_for_variable,
    variable_for_tr_id = variable_for_tr_id
  )
}

#' Match a pheno `bpid` to TPP rows by normalized BPId
#'
#' Normalizes the pheno `bpid` (e.g. "DBP0003" -> 3) and the TPP file's `BPId`
#' column, then returns the subset of TPP rows whose normalized `BPId` equals it.
#'
#' @param tpp A data.frame from `read_tpp_file()`.
#' @param bpid The pheno `bpid` value(s); a single value is expected but a vector
#'   is tolerated (its first non-NA normalized key is used).
#'
#' @return A data.frame: the matching subset of `tpp` (zero rows if no match).
#' @export
match_bp <- function(tpp, bpid) {
  key <- normalize_id(bpid)
  key <- key[!is.na(key)]
  if (length(key) == 0L) {
    return(tpp[integer(0), , drop = FALSE])
  }
  key <- key[1L]
  tpp[which(normalize_id(tpp$BPId) == key), , drop = FALSE]
}

#' Match a pheno `msid` to TPP rows by normalized MSId
#'
#' Normalizes the pheno `msid` and the TPP file's `MSId` column, then returns the
#' subset of TPP rows whose normalized `MSId` equals it. MS <-> TPP is 1:1, so a
#' successful match should resolve to exactly one distinct TPPId (callers in
#' Task 6 enforce the 0 / 1 / >1 decision).
#'
#' @param tpp A data.frame from `read_tpp_file()`.
#' @param msid The pheno `msid` value(s).
#'
#' @return A data.frame: the matching subset of `tpp` (zero rows if no match).
#' @export
match_ms <- function(tpp, msid) {
  key <- normalize_id(msid)
  key <- key[!is.na(key)]
  if (length(key) == 0L) {
    return(tpp[integer(0), , drop = FALSE])
  }
  key <- key[1L]
  tpp[which(normalize_id(tpp$MSId) == key), , drop = FALSE]
}

#' Distinct TPPs (TPPId -> TPPCode) present in a set of TPP rows
#'
#' Helper used by both the full path (Task 6) and the bpid-only path (Task 7) to
#' enumerate the distinct TPPs a match resolved to.
#'
#' @param tpp_rows A data.frame subset of a TPP file.
#'
#' @return A data.frame with columns `TPPId`, `TPPCode`, `tpp_key` (normalized
#'   TPPId), one row per distinct TPP, ordered by `tpp_key`.
#' @export
distinct_tpps <- function(tpp_rows) {
  if (nrow(tpp_rows) == 0L) {
    return(data.frame(TPPId = character(0), TPPCode = character(0),
                      tpp_key = integer(0), stringsAsFactors = FALSE))
  }
  idx <- !duplicated(tpp_rows$TPPId)
  out <- data.frame(
    TPPId    = tpp_rows$TPPId[idx],
    TPPCode  = tpp_rows$TPPCode[idx],
    tpp_key  = normalize_id(tpp_rows$TPPId[idx]),
    stringsAsFactors = FALSE
  )
  out[order(out$tpp_key), , drop = FALSE]
}

#' Resolve a TPP row's TraitId to the matching pheno trait column
#'
#' Uses the resolver from `read_tr_mapping()` to turn a TPP-file `TraitId`
#' (e.g. "TR-00388") into the pheno variable abbreviation / column name
#' (e.g. "WH_PH_M_CM"), matching on the normalized TR_Id key. Returns `NA` when
#' no pheno trait maps to the TraitId (Requirement 2.4 leaves `pheno_trait`
#' empty in that case; the build tasks handle that).
#'
#' @param tr_mapping The list returned by `read_tr_mapping()`.
#' @param trait_id   One or more TPP-file `TraitId` values.
#'
#' @return A character vector of pheno variable abbreviations (NA where unmatched).
#' @export
resolve_pheno_trait <- function(tr_mapping, trait_id) {
  tr_mapping$variable_for_tr_id(trait_id)
}
