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
# Name             : tpp_gating_helpers
# Description      : Opt-in `useBP` gating for automatic TPP population in the
#                    Bioflow RData generation script (spec tpp-review Task 5,
#                    Requirement 3). These are PURE, side-effect-free decision
#                    helpers plus a small diagnostics-JSON writer used on the
#                    fail-fast STOP branches. `getBioflowRData()` consults
#                    `tpp_gate_decision()` BEFORE building or writing any RData;
#                    when the decision is a STOP it writes a gating-failure
#                    diagnostics JSON (via `write_gating_diagnostics_json()`) and
#                    stops WITHOUT writing an RData file.
# Maintainer       : Bioflow core development team
# ---------------------------------------------------------------------------------------------
#
# DESIGN (Requirement 3, design.md "RData script entry point and opt-in gating")
# ------------------------------------------------------------------------------
# Three new, backward-compatible arguments drive the gate:
#
#   useBP           logical, default FALSE  - master opt-in switch.
#   BP_file         character, default NULL  - FULL path (incl. filename) to the
#                                              sensitive breeding-portal TPP file.
#   TR_mapping_file character, default NULL  - FILENAME ONLY, resolved as
#                                              file.path(dirname(phenotypeFile),
#                                                        TR_mapping_file).
#
# "Usable bpid info": the `bpid` column is PRESENT AND at least one row has a
# non-NA `bpid` value. An `msid` column alone is NOT sufficient and `msid` is
# NOT required.
#
# Gating decision flow (evaluated before any TPP population begins):
#
#   useBP == FALSE                          -> proceed=TRUE,  populate=FALSE (RData as today)
#   useBP == TRUE & !usable-bpid            -> STOP (no RData), diagnostic 3.4
#   useBP == TRUE &  usable-bpid & BP_file NULL/unreadable
#                                           -> STOP (no RData), diagnostic 3.5
#   useBP == TRUE &  usable-bpid & TR_mapping_file NULL/unresolvable
#                                           -> STOP (no RData), diagnostic 3.6
#   else                                    -> proceed=TRUE,  populate=TRUE
#
# On any STOP branch the caller must be fail-fast and side-effect-free: write NO
# RData, record the diagnostic in the diagnostics JSON, and stop.

# ---- "Usable bpid" check (Requirement 3.3) ----------------------------------

#' Does the phenotype data carry usable breeding-pipeline (bpid) information?
#'
#' Usable means: the `bpid` column is PRESENT AND at least one row has a non-NA
#' `bpid` value. Blank/whitespace-only strings are treated as NA (the pheno CSV
#' is read with `check.names = FALSE`, so an empty cell arrives as "").
#' An `msid` column alone is NOT sufficient and `msid` is NOT required.
#'
#' @param data_pheno A phenotype data.frame (or NULL).
#'
#' @return `TRUE` when the `bpid` column exists and has >= 1 non-NA, non-blank
#'   value; `FALSE` otherwise (including when `data_pheno` is NULL or has no
#'   `bpid` column).
#' @export
tpp_has_usable_bpid <- function(data_pheno) {
  if (is.null(data_pheno) || !is.data.frame(data_pheno)) {
    return(FALSE)
  }
  if (!("bpid" %in% names(data_pheno))) {
    return(FALSE)
  }
  bpid <- data_pheno[["bpid"]]
  # Treat blank / whitespace-only strings as "no value" alongside NA.
  chr <- trimws(as.character(bpid))
  usable <- !is.na(bpid) & !is.na(chr) & nzchar(chr)
  any(usable)
}

# ---- Gating decision (Requirement 3.2, 3.4, 3.5, 3.6, 3.7) -------------------

#' Decide whether automatic TPP population may proceed (pure, side-effect-free)
#'
#' Implements the fail-fast opt-in gate. This function reads NO files and writes
#' NOTHING; it only inspects the already-loaded phenotype data and the file
#' arguments and returns a decision record. The caller performs the side effects
#' (writing the diagnostic, stopping, or proceeding).
#'
#' `BP_file` readability is checked with `file.exists()` on the verbatim full
#' path. `TR_mapping_file` is resolved as `file.path(dirname(phenotypeFile),
#' TR_mapping_file)` and checked the same way; the resolved path is returned in
#' the decision so the caller / population step can reuse it.
#'
#' @param useBP           logical master switch.
#' @param data_pheno      phenotype data.frame (already loaded) or NULL.
#' @param phenotypeFile   path to the phenotype file (used to derive the data
#'   directory via `dirname()` for resolving `TR_mapping_file`).
#' @param BP_file         full path (incl. filename) to the TPP file, or NULL.
#' @param TR_mapping_file bare filename for the TR mapping, or NULL.
#'
#' @return A list with:
#'   \describe{
#'     \item{proceed}{logical - TRUE when the RData should still be generated
#'       (useBP=FALSE, or useBP=TRUE and all preconditions met). FALSE means a
#'       fail-fast STOP: write NO RData.}
#'     \item{populate}{logical - TRUE only when TPP population should run.}
#'     \item{status}{one of "no_optin", "stop_no_usable_bpid",
#'       "stop_bp_file", "stop_tr_mapping_file", "proceed_populate".}
#'     \item{requirement}{the Requirement 3.x tag for a STOP, else NA.}
#'     \item{reason}{plain-language diagnostic message (NA when not a STOP).}
#'     \item{bp_file}{the BP_file path (echoed for the population step).}
#'     \item{tr_mapping_path}{the resolved TR_mapping_file path (or NA).}
#'   }
#' @export
tpp_gate_decision <- function(useBP, data_pheno, phenotypeFile,
                              BP_file = NULL, TR_mapping_file = NULL) {

  mk <- function(proceed, populate, status, requirement = NA_character_,
                 reason = NA_character_, tr_mapping_path = NA_character_) {
    list(
      proceed         = proceed,
      populate        = populate,
      status          = status,
      requirement     = requirement,
      reason          = reason,
      bp_file         = if (is.null(BP_file)) NA_character_ else BP_file,
      tr_mapping_path = tr_mapping_path
    )
  }

  # Requirement 3.2: useBP=FALSE -> generate RData as today, no file reads,
  # no TPP population, regardless of bpid/msid.
  if (!isTRUE(useBP)) {
    return(mk(proceed = TRUE, populate = FALSE, status = "no_optin"))
  }

  # useBP == TRUE from here on.

  # Requirement 3.4: no usable bpid -> fail-fast STOP.
  if (!tpp_has_usable_bpid(data_pheno)) {
    return(mk(
      proceed     = FALSE,
      populate    = FALSE,
      status      = "stop_no_usable_bpid",
      requirement = "3.4",
      reason      = paste0(
        "useBP = TRUE but the phenotype file lacks usable bpid information ",
        "(the 'bpid' column is absent, or every bpid value is NA/blank). ",
        "At least one occurrence with a non-NA bpid is required."
      )
    ))
  }

  # Requirement 3.5: BP_file NULL or not readable -> fail-fast STOP.
  bp_ok <- !is.null(BP_file) &&
    is.character(BP_file) && length(BP_file) == 1L && nzchar(BP_file) &&
    file.exists(BP_file)
  if (!bp_ok) {
    return(mk(
      proceed     = FALSE,
      populate    = FALSE,
      status      = "stop_bp_file",
      requirement = "3.5",
      reason      = paste0(
        "useBP = TRUE and the phenotype file has usable bpid information, but ",
        "BP_file is ",
        if (is.null(BP_file) || !nzchar(as.character(BP_file)[1L])) {
          "NULL/empty"
        } else {
          paste0("not readable at: ", BP_file)
        },
        ". BP_file must be a full path (including filename) to the ",
        "breeding-portal TPP file (TPP_traits_for_Bioflow.csv)."
      )
    ))
  }

  # Requirement 3.6: TR_mapping_file NULL or unresolvable in dirname(phenotypeFile).
  tr_null <- is.null(TR_mapping_file) ||
    !is.character(TR_mapping_file) || length(TR_mapping_file) != 1L ||
    !nzchar(TR_mapping_file)
  data_dir <- if (!is.null(phenotypeFile) && nzchar(phenotypeFile)) {
    dirname(phenotypeFile)
  } else {
    NA_character_
  }
  tr_resolved <- if (!tr_null && !is.na(data_dir)) {
    file.path(data_dir, TR_mapping_file)
  } else {
    NA_character_
  }
  tr_ok <- !tr_null && !is.na(tr_resolved) && file.exists(tr_resolved)
  if (!tr_ok) {
    return(mk(
      proceed     = FALSE,
      populate    = FALSE,
      status      = "stop_tr_mapping_file",
      requirement = "3.6",
      reason      = paste0(
        "useBP = TRUE and the phenotype file has usable bpid information, but ",
        "TR_mapping_file is ",
        if (tr_null) {
          "NULL/empty"
        } else {
          paste0("not resolvable in the phenotype directory (expected at: ",
                 tr_resolved, ")")
        },
        ". TR_mapping_file must be a filename (e.g. 'variable_tr_mapping.csv') ",
        "sitting in the same directory as the phenotype input."
      ),
      tr_mapping_path = tr_resolved
    ))
  }

  # Requirement 3.7: all preconditions met -> proceed to population.
  mk(
    proceed         = TRUE,
    populate        = TRUE,
    status          = "proceed_populate",
    tr_mapping_path = tr_resolved
  )
}

# ---- Gating-failure diagnostics JSON (Requirement 3.4 / 3.5 / 3.6) ----------

#' Write a gating-failure diagnostics JSON (fail-fast STOP branches)
#'
#' On a STOP branch no RData is written, so `validate_bioflow_object()` never
#' runs and never emits its diagnostics JSON. This writer produces a small,
#' engineer-friendly JSON at the same path the normal diagnostics would use, so
#' the EBS engineer always gets a machine-readable artifact explaining why the
#' run stopped. It reuses the same jsonlite settings as the main diagnostics
#' writer for consistency.
#'
#' @param decision        the STOP decision list from `tpp_gate_decision()`.
#' @param path            destination path for the diagnostics JSON.
#' @param inputs          optional named list of input file paths (echoed).
#'
#' @return Invisibly, `path` (or NULL when jsonlite is unavailable / path empty).
#' @export
write_gating_diagnostics_json <- function(decision, path, inputs = NULL) {
  if (is.null(path) || !nzchar(path)) return(invisible(NULL))
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    cli::cli_warn("Package 'jsonlite' is not available; skipping gating diagnostics JSON.")
    return(invisible(NULL))
  }

  payload <- list(
    schema_version = "1.0",
    generated_at   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    audience       = "EBS software engineers preparing Bioflow input files",
    inputs         = inputs,
    overall = list(
      ready_use_cases     = list(),
      not_ready_use_cases = list(),
      ready_count         = 0L,
      total_use_cases     = 0L,
      verdict             = paste0(
        "No RData was written: automatic TPP population (useBP = TRUE) was ",
        "requested but a precondition was not met. See 'tpp_gating' below."
      )
    ),
    tpp_gating = list(
      useBP        = TRUE,
      populated    = FALSE,
      rdata_written = FALSE,
      status       = decision$status,
      requirement  = decision$requirement,
      what_it_means = decision$reason,
      how_to_fix   = tpp_gating_fix_action(decision$status)
    ),
    use_cases = list()
  )

  json <- jsonlite::toJSON(payload, pretty = TRUE, auto_unbox = TRUE,
                           null = "null", na = "null")
  writeLines(json, path, useBytes = TRUE)
  cli::cli_inform(paste0("Gating diagnostics written to: ", path))
  invisible(path)
}

#' Plain-language fix action for a gating STOP status
#' @noRd
tpp_gating_fix_action <- function(status) {
  switch(
    status,
    stop_no_usable_bpid = paste0(
      "Provide a phenotype file whose 'bpid' column has at least one non-NA, ",
      "non-blank value, or set useBP = FALSE to generate the RData without TPP ",
      "population."
    ),
    stop_bp_file = paste0(
      "Pass BP_file as the full path (including filename) to the breeding-portal ",
      "TPP file (TPP_traits_for_Bioflow.csv) and ensure the file is readable."
    ),
    stop_tr_mapping_file = paste0(
      "Place the TR mapping file (e.g. variable_tr_mapping.csv) in the same ",
      "directory as the phenotype input and pass its bare filename as ",
      "TR_mapping_file."
    ),
    ""
  )
}
