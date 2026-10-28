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
# Name             : tpp_build_helpers
# Description      : Build the auto-populated TPP object (data$TPP[[id]] +
#                    metadata$TPP[[id]]) EXACTLY as a Bioflow TPP Excel upload
#                    produces it, from the sensitive breeding-portal TPP file.
#                    This file covers the FULL PATH (spec tpp-review Task 6,
#                    Requirement 4): a phenotype record carrying BOTH a usable
#                    bpid AND a non-empty msid. The TPP is located via
#                    msid -> MSId (1:1), and the object is built for that single
#                    TPP. The bpid-only path (Task 7, Requirement 5) will add a
#                    sibling builder that reuses the same primitives with no
#                    environment subsetting.
# Maintainer       : Bioflow core development team
# ---------------------------------------------------------------------------------------------
#
# DEPENDENCIES
# ------------
#   * normalize_id()          - scripts/tpp_id_helpers.R          (bundled first)
#   * match_ms(), distinct_tpps(), resolve_pheno_trait()
#                             - scripts/tpp_population_helpers.R  (bundled before this)
#
# STRUCTURAL EQUIVALENCE (the guiding principle)
# ----------------------------------------------
# A Bioflow TPP upload (mod_getDataTPP.R / fct_tpp_helpers.R) produces:
#
#   data$TPP[[tpp_id]]         raw standardized sheet data.frame with the Excel
#                              upload headers. Only the checks column is
#                              standardized (standardize_checks_column):
#                              gsub [;/|+] -> comma, split, trim, drop empties,
#                              rejoin with ", "; NA preserved.
#
#   metadata$TPP[[tpp_id]]     list(
#     $traits              data.frame(tpp_trait, pheno_trait, trait_id)
#     $env_filters         NULL or named list keyed by Trait ID (only traits
#                          whose kept-env set differs from the all-envs set)
#     $scale_validations   named by Trait Name:
#                            list(is_violation, pct, scale_min, scale_max)
#     $checks_per_trait    named by Trait Name: char vector(s) of designations
#     $checks              data.frame(tpp_check, tpp_trait, pheno_designation)
#   )
#
# The Excel upload column headers (authoritative, from the Bioflow upload
# template) are, in order:
#   "Trait ID", "Trait Type", "Trait Name", "Trait Requirement",
#   "Desired Score", "Scale Option", "Best Checks", "Value Proposition".
# The sensitive CSV (BP_file) is long/tidy with no-space column names; the map
# below reproduces the upload sheet from it so downstream consumers (and the
# Review TPP module) treat an auto-populated TPP identically to an uploaded one.
#
# ENVIRONMENTS ON THE FULL PATH (Requirement 4.6)
# -----------------------------------------------
# Env per mapped trait = every occurrence mapped to that msid under the located
# TPP. The TPP file is tidy on TRAITS, not on environments, so "every occurrence
# mapped to that msid" is the full environment set for the data. Because the kept
# env set therefore equals the all-envs set, no trait gets an env_filters entry
# and $env_filters is NULL -- exactly what an upload with no filtering stores.

# ---- Upload-sheet column map (CSV name -> Bioflow upload header) -------------
# Ordered to match the Bioflow TPP Excel upload template. The checks column
# ("Best Checks") is detected downstream by the "check" substring rule, exactly
# like the upload flow's get_default_checks_col().
.TPP_UPLOAD_COL_MAP <- c(
  "Trait ID"          = "TraitId",
  "Trait Type"        = "TraitType",
  "Trait Name"        = "TraitName",
  "Trait Requirement" = "TraitRequirement",
  "Desired Score"     = "DesiredScore",
  "Scale Option"      = "ScaleOption",
  "Best Checks"       = "BestCheck",
  "Value Proposition" = "ValueProposition"
)

#' Default checks column name (port of bioflow get_default_checks_col)
#'
#' Returns the first column name containing "check" (case-insensitive), or NULL
#' if none. Byte-for-byte behavioral copy of bioflow's `get_default_checks_col`
#' so the auto-populated sheet's checks column is detected identically to the
#' upload flow.
#'
#' @param col_names Character vector of column names.
#' @return The first matching column name, or NULL.
#' @export
get_default_checks_col <- function(col_names) {
  matches <- grepl("check", col_names, ignore.case = TRUE)
  if (any(matches)) {
    return(col_names[which(matches)[1]])
  }
  NULL
}

#' Parse a "Scale Option" text into numeric min/max (port of bioflow helper)
#'
#' Behavioral copy of bioflow's `tpp_parse_scale_option` (fct_tpp_helpers.R):
#' supports "X to Y" and "X - Y" / "X-Y" numeric ranges; returns
#' `list(min, max)` or NULL when unparseable. Used so auto-populated
#' `scale_validations` entries carry the same bounds the upload flow derives.
#'
#' @param scale_option Character scalar from the "Scale Option" column.
#' @return list(min, max) or NULL.
#' @export
tpp_parse_scale_option <- function(scale_option) {
  if (is.null(scale_option) || length(scale_option) != 1 || is.na(scale_option)) {
    return(NULL)
  }
  scale_option <- trimws(as.character(scale_option))
  if (nchar(scale_option) == 0) {
    return(NULL)
  }

  match_to <- regmatches(scale_option,
    regexec("^\\s*([0-9]*\\.?[0-9]+)\\s+to\\s+([0-9]*\\.?[0-9]+)\\s*$",
            scale_option, ignore.case = TRUE))[[1]]
  if (length(match_to) == 3) {
    min_val <- suppressWarnings(as.numeric(match_to[2]))
    max_val <- suppressWarnings(as.numeric(match_to[3]))
    if (!is.na(min_val) && !is.na(max_val)) {
      return(list(min = min(min_val, max_val), max = max(min_val, max_val)))
    }
  }

  match_dash <- regmatches(scale_option,
    regexec("^\\s*([0-9]*\\.?[0-9]+)\\s*-\\s*([0-9]*\\.?[0-9]+)\\s*$",
            scale_option))[[1]]
  if (length(match_dash) == 3) {
    min_val <- suppressWarnings(as.numeric(match_dash[2]))
    max_val <- suppressWarnings(as.numeric(match_dash[3]))
    if (!is.na(min_val) && !is.na(max_val)) {
      return(list(min = min(min_val, max_val), max = max(min_val, max_val)))
    }
  }

  NULL
}

#' Standardize a checks column vector (port of bioflow standardize_checks_column)
#'
#' Byte-for-byte behavioral copy of bioflow's
#' `R/mod_getDataTPP.R::standardize_checks_column()` so the auto-populated
#' `data$TPP[[id]]` checks column is standardized identically to an upload:
#' replace `[;/|+]` with a comma, split on comma, trim each entry, drop empties,
#' rejoin with ", "; preserve NA (and the literal string "NA").
#'
#' @param col_vector A character vector (coerced if not character).
#' @return A character vector the same length as `col_vector`.
#' @export
tpp_standardize_checks_column <- function(col_vector) {
  if (!is.character(col_vector)) {
    col_vector <- as.character(col_vector)
  }
  vapply(col_vector, function(cell) {
    if (is.na(cell) || cell == "NA") {
      return(NA_character_)
    }
    standardized <- gsub("[;/|\\+]", ",", cell)
    entries <- strsplit(standardized, ",", fixed = TRUE)[[1]]
    entries <- trimws(entries)
    entries <- entries[nchar(entries) > 0]
    if (length(entries) == 0) {
      return(NA_character_)
    }
    paste(entries, collapse = ", ")
  }, character(1), USE.NAMES = FALSE)
}

#' Build the raw standardized upload sheet for a single TPP's rows
#'
#' Reproduces `data$TPP[[tpp_id]]` as the Bioflow upload flow stores it: the
#' Excel upload sheet (upload header names, in template order) restricted to the
#' located TPP's trait rows, with the "Best Checks" column standardized exactly
#' like `standardize_checks_column`. Any upload column whose source CSV column is
#' absent is emitted as all-NA so the sheet shape is stable.
#'
#' @param tpp_rows data.frame subset of the TPP file for a single TPP (one row
#'   per TPP-trait), as returned by subsetting `read_tpp_file()` output.
#' @return A data.frame with the Bioflow upload headers; the checks column
#'   standardized.
#' @export
build_tpp_data_sheet <- function(tpp_rows) {
  stopifnot(is.data.frame(tpp_rows))
  upload_headers <- names(.TPP_UPLOAD_COL_MAP)
  n <- nrow(tpp_rows)

  cols <- lapply(upload_headers, function(h) {
    src <- .TPP_UPLOAD_COL_MAP[[h]]
    if (!is.null(src) && src %in% names(tpp_rows)) {
      as.character(tpp_rows[[src]])
    } else {
      rep(NA_character_, n)
    }
  })
  names(cols) <- upload_headers

  sheet <- as.data.frame(cols, stringsAsFactors = FALSE, check.names = FALSE)

  # Standardize the auto-detected checks column, mirroring the upload flow.
  checks_col <- get_default_checks_col(names(sheet))
  if (!is.null(checks_col)) {
    sheet[[checks_col]] <- tpp_standardize_checks_column(sheet[[checks_col]])
  }
  sheet
}

#' Build metadata$TPP[[id]]$traits for an auto-populated TPP (Requirement 4.5)
#'
#' Produces the 3-column traits data.frame (`tpp_trait`, `pheno_trait`,
#' `trait_id`) exactly as the upload flow stores it, INCLUDING ONLY TPP traits
#' whose normalized `TraitId` resolves to a pheno trait via the TR mapping
#' (Requirement 4.5: include only traits whose normalized TraitId equals a
#' normalized pheno TR_Id). `pheno_trait` is the resolved pheno column;
#' `trait_id` is the TPP TraitId verbatim.
#'
#' @param tpp_rows   data.frame subset for a single TPP.
#' @param tr_mapping resolver list from `read_tr_mapping()`.
#' @return A data.frame(tpp_trait, pheno_trait, trait_id); zero rows when nothing
#'   maps.
#' @export
build_tpp_traits_df <- function(tpp_rows, tr_mapping) {
  stopifnot(is.data.frame(tpp_rows))

  trait_names <- as.character(tpp_rows[["TraitName"]])
  trait_ids   <- as.character(tpp_rows[["TraitId"]])
  pheno       <- resolve_pheno_trait(tr_mapping, trait_ids)

  keep <- !is.na(pheno) & nzchar(pheno)

  data.frame(
    tpp_trait   = trait_names[keep],
    pheno_trait = pheno[keep],
    trait_id    = trait_ids[keep],
    stringsAsFactors = FALSE
  )
}

#' Build metadata$TPP[[id]]$scale_validations for an auto-populated TPP
#'
#' Keyed by Trait Name, each value is `list(is_violation, pct, scale_min,
#' scale_max)`, mirroring the upload flow. Because the RData script has no live
#' per-plot phenotype values to validate the scale against at build time, a
#' freshly auto-populated TPP carries no realized violations: entries are emitted
#' only for traits whose `ScaleOption` parses to a numeric range, with
#' `is_violation = FALSE` and `pct = 0` (the same shape the upload stores before
#' any pheno values flag a violation). This keeps the element present and
#' upload-shaped without inventing data.
#'
#' @param traits_df   the mapped traits data.frame from `build_tpp_traits_df()`.
#' @param tpp_rows    data.frame subset for a single TPP (source of ScaleOption).
#' @return A named list keyed by Trait Name (possibly empty `list()`).
#' @export
build_tpp_scale_validations <- function(traits_df, tpp_rows) {
  out <- list()
  if (is.null(traits_df) || nrow(traits_df) == 0L) return(out)
  if (!("ScaleOption" %in% names(tpp_rows))) return(out)

  # Index ScaleOption by Trait Name for lookup.
  scale_by_name <- tpp_rows[["ScaleOption"]]
  names(scale_by_name) <- as.character(tpp_rows[["TraitName"]])

  for (nm in traits_df$tpp_trait) {
    scale_text <- scale_by_name[[nm]]
    parsed <- tpp_parse_scale_option(scale_text)
    if (!is.null(parsed)) {
      out[[nm]] <- list(
        is_violation = FALSE,
        pct          = 0,
        scale_min    = parsed$min,
        scale_max    = parsed$max
      )
    }
  }
  out
}

#' Auto-map TPP best-check designations against the pheno designations
#'
#' Implements Requirement 4.7 / 4.8: a TPP best-check designation is auto-mapped
#' when it EXACTLY equals a designation string present in the pheno data;
#' unmatched check designations are left unmapped and processing continues.
#' Returns both upload-shaped structures:
#'   * `checks_per_trait` - named by Trait Name, the exactly-matched designations
#'   * `checks`           - data.frame(tpp_check, tpp_trait, pheno_designation),
#'                          one row per (check name, mapped designation); traits
#'                          with no mapped designation emit an unmapped row per
#'                          check name (pheno_designation = NA), mirroring the
#'                          upload flow.
#'
#' Only traits present in `traits_df` (i.e. mapped to a pheno trait) are
#' considered, matching the upload flow which only offers checks for mapped
#' traits whose standardized checks column is non-empty.
#'
#' @param traits_df         the mapped traits data.frame.
#' @param tpp_rows          data.frame subset for a single TPP.
#' @param pheno_designations character vector of designations present in the
#'   pheno data (used for the exact-match auto-mapping).
#' @return list(checks_per_trait = <named list>, checks = <data.frame>).
#' @export
build_tpp_checks <- function(traits_df, tpp_rows, pheno_designations) {
  checks_per_trait <- list()
  checks_rows <- list()

  if (is.null(pheno_designations)) pheno_designations <- character(0)
  pheno_designations <- as.character(pheno_designations)

  if (!is.null(traits_df) && nrow(traits_df) > 0L) {
    # Standardized "Best Checks" column keyed by Trait Name.
    checks_col <- get_default_checks_col(names(tpp_rows))
    best_by_name <- if (!is.null(checks_col)) {
      v <- tpp_standardize_checks_column(tpp_rows[[checks_col]])
      stats::setNames(v, as.character(tpp_rows[["TraitName"]]))
    } else {
      stats::setNames(rep(NA_character_, nrow(tpp_rows)),
                      as.character(tpp_rows[["TraitName"]]))
    }

    for (nm in traits_df$tpp_trait) {
      best <- best_by_name[[nm]]
      # Associated check names = the standardized best-check entries for the trait.
      associated <- character(0)
      if (!is.null(best) && !is.na(best) && nzchar(best)) {
        associated <- trimws(strsplit(best, ",", fixed = TRUE)[[1]])
        associated <- associated[nchar(associated) > 0 & associated != "NA"]
        associated <- unique(associated)
      }
      if (length(associated) == 0L) next  # no checks column -> trait skipped

      # Auto-map: a best-check designation that EXACTLY equals a pheno
      # designation is mapped; others are left unmapped (Req 4.7 / 4.8).
      matched <- associated[associated %in% pheno_designations]

      if (length(matched) > 0L) {
        checks_per_trait[[nm]] <- matched
        for (check_name in associated) {
          for (desig in matched) {
            checks_rows <- c(checks_rows, list(data.frame(
              tpp_check         = check_name,
              tpp_trait         = nm,
              pheno_designation = desig,
              stringsAsFactors  = FALSE
            )))
          }
        }
      } else {
        # Unmapped: one row per check name, pheno_designation = NA.
        for (check_name in associated) {
          checks_rows <- c(checks_rows, list(data.frame(
            tpp_check         = check_name,
            tpp_trait         = nm,
            pheno_designation = NA_character_,
            stringsAsFactors  = FALSE
          )))
        }
      }
    }
  }

  checks_df <- if (length(checks_rows) > 0L) {
    do.call(rbind, checks_rows)
  } else {
    data.frame(
      tpp_check         = character(0),
      tpp_trait         = character(0),
      pheno_designation = character(0),
      stringsAsFactors  = FALSE
    )
  }

  list(checks_per_trait = checks_per_trait, checks = checks_df)
}

#' Build the full metadata$TPP[[id]] list for an auto-populated TPP
#'
#' Assembles the complete upload-shaped metadata list with exactly the five
#' elements `{traits, env_filters, scale_validations, checks_per_trait, checks}`.
#'
#' `env_filters` is `NULL` on the full path: env per trait is every occurrence
#' mapped to the msid (the all-envs set), so no trait's kept-env set differs from
#' all-envs and no entry is stored -- identical to an upload with no filtering.
#'
#' @param tpp_rows           data.frame subset for a single TPP.
#' @param tr_mapping         resolver list from `read_tr_mapping()`.
#' @param pheno_designations designations present in the pheno data.
#' @return A named list(traits, env_filters, scale_validations,
#'   checks_per_trait, checks).
#' @export
build_tpp_metadata <- function(tpp_rows, tr_mapping, pheno_designations = NULL) {
  traits_df <- build_tpp_traits_df(tpp_rows, tr_mapping)
  checks    <- build_tpp_checks(traits_df, tpp_rows, pheno_designations)

  list(
    traits            = traits_df,
    env_filters       = NULL,  # full path keeps all envs -> no filter entry
    scale_validations = build_tpp_scale_validations(traits_df, tpp_rows),
    checks_per_trait  = checks$checks_per_trait,
    checks            = checks$checks
  )
}

#' Full-path TPP population for a single phenotype record (Requirement 4)
#'
#' Given the loaded TPP file and a phenotype record carrying BOTH a usable bpid
#' and a non-empty msid, locate the TPP via msid -> MSId, enforce the 1:1 rule,
#' and -- on exactly one match -- build the upload-equivalent `data$TPP[[id]]` +
#' `metadata$TPP[[id]]`. The function is PURE: it returns a decision record and
#' the built pieces; the caller attaches them to `result` and records diagnostics.
#'
#' Decision outcomes (Requirement 4.2 / 4.3 / 4.4):
#'   * zero MS matches  -> status "no_match"   (skip, leave object unchanged)
#'   * >1 distinct TPP  -> status "not_1to1"   (skip, 1:1 violated)
#'   * exactly one TPP  -> status "populated"  (tpp_id + data + metadata returned)
#'
#' @param tpp_data           the full TPP file data.frame (`read_tpp_file()`).
#' @param tr_mapping         resolver list from `read_tr_mapping()`.
#' @param msid               the phenotype record's msid (non-empty on this path).
#' @param pheno_designations designations present in the pheno data (for checks).
#' @return A list with:
#'   \describe{
#'     \item{status}{"no_match", "not_1to1", or "populated".}
#'     \item{tpp_id}{the TPPCode (e.g. "TPP00697") when populated, else NA.}
#'     \item{data}{the upload sheet data.frame when populated, else NULL.}
#'     \item{metadata}{the metadata list when populated, else NULL.}
#'     \item{diagnostic}{plain-language message for no_match / not_1to1, else NA.}
#'     \item{tpp_codes}{the distinct TPPCodes the msid matched (for diagnostics).}
#'   }
#' @export
populate_tpp_full_path <- function(tpp_data, tr_mapping, msid,
                                    pheno_designations = NULL) {
  ms_rows <- match_ms(tpp_data, msid)
  tpps    <- distinct_tpps(ms_rows)

  # Requirement 4.2: zero MS matches -> skip, leave object unchanged.
  if (nrow(tpps) == 0L) {
    return(list(
      status     = "no_match",
      tpp_id     = NA_character_,
      data       = NULL,
      metadata   = NULL,
      diagnostic = sprintf(
        "No matching TPP was found for msid '%s' (no TPP file row has a market segment whose normalized MSId equals it); TPP population skipped and data$TPP / metadata$TPP left unchanged.",
        as.character(msid)[1L]
      ),
      tpp_codes  = character(0)
    ))
  }

  # Requirement 4.3: more than one distinct TPP -> 1:1 violated, skip.
  if (nrow(tpps) > 1L) {
    return(list(
      status     = "not_1to1",
      tpp_id     = NA_character_,
      data       = NULL,
      metadata   = NULL,
      diagnostic = sprintf(
        "The 1:1 Market Segment <-> TPP relationship was violated: msid '%s' matched %d distinct TPPs (%s); TPP population skipped.",
        as.character(msid)[1L], nrow(tpps), paste(tpps$TPPCode, collapse = ", ")
      ),
      tpp_codes  = tpps$TPPCode
    ))
  }

  # Requirement 4.4: exactly one -> build the upload-equivalent object.
  tpp_code <- tpps$TPPCode[1L]
  tpp_rows <- ms_rows[ms_rows$TPPCode == tpp_code, , drop = FALSE]

  list(
    status     = "populated",
    tpp_id     = tpp_code,
    data       = build_tpp_data_sheet(tpp_rows),
    metadata   = build_tpp_metadata(tpp_rows, tr_mapping, pheno_designations),
    diagnostic = NA_character_,
    tpp_codes  = tpp_code
  )
}

# ---------------------------------------------------------------------------------------------
# BPID-ONLY PATH (spec tpp-review Task 7, Requirement 5)
# ---------------------------------------------------------------------------------------------
# When the phenotype data carries a usable bpid but the msid column is blank or
# absent for ALL rows, there is no single market segment to locate a single TPP
# by. Instead the script populates EVERY distinct TPP grouped under the matched
# BP, and for every trait in every such TPP it retains ALL environments with NO
# subsetting (Requirement 5.3): the per-trait environment count equals the count
# present in the source data.
#
# This reuses the exact same Task 6 build primitives (build_tpp_data_sheet,
# build_tpp_metadata) so an auto-populated TPP from the bpid-only path is
# structurally identical to one from an upload (structural equivalence). The
# "no environment subsetting" requirement is satisfied for free: build_tpp_metadata
# emits env_filters = NULL (no per-trait env restriction), which is precisely the
# upload shape for a TPP with no filtering -> every trait keeps all environments.
#
# The function is PURE: it returns a decision record plus the per-TPP built
# pieces; the caller attaches them to `result` and records diagnostics.

#' bpid-only TPP population across all TPPs under a matched BP (Requirement 5)
#'
#' Given the loaded TPP file, the TR mapping resolver, and a phenotype `bpid`
#' (msid empty/absent across all rows), identify every distinct TPP grouped under
#' the matched normalized `BPId` (Requirement 5.1) and build the full
#' upload-equivalent object for each, retaining ALL environments for every trait
#' (no subsetting; Requirement 5.3).
#'
#' Decision outcomes:
#'   * no BP rows match the normalized bpid -> status "no_bp_match": skip TPP
#'     population; the caller generates the RData WITHOUT TPP info and records a
#'     diagnostic (Requirement 5.2 / 5.4).
#'   * >=1 distinct TPP under the BP          -> status "populated": `tpps` holds
#'     one built entry per distinct TPP (tpp_id + data + metadata).
#'
#' @param tpp_data           the full TPP file data.frame (`read_tpp_file()`).
#' @param tr_mapping         resolver list from `read_tr_mapping()`.
#' @param bpid               the phenotype record's bpid (usable on this path).
#' @param pheno_designations designations present in the pheno data (for checks).
#' @return A list with:
#'   \describe{
#'     \item{status}{"no_bp_match" or "populated".}
#'     \item{tpps}{a named list (keyed by TPPCode) of entries
#'       `list(tpp_id, data, metadata)`; empty `list()` when no_bp_match.}
#'     \item{tpp_codes}{the distinct TPPCodes populated (character(0) when none).}
#'     \item{diagnostic}{plain-language message for no_bp_match, else NA.}
#'   }
#' @export
populate_tpp_bpid_only <- function(tpp_data, tr_mapping, bpid,
                                   pheno_designations = NULL) {
  bp_rows <- match_bp(tpp_data, bpid)

  # Requirement 5.2 / 5.4: no TPP file rows match the BP -> skip, RData without
  # TPP population, diagnostic.
  if (nrow(bp_rows) == 0L) {
    return(list(
      status     = "no_bp_match",
      tpps       = list(),
      tpp_codes  = character(0),
      diagnostic = sprintf(
        "No TPPs were found for bpid '%s' (no TPP file row has a breeding pipeline whose normalized BPId equals it); TPP population skipped and the RData is generated without TPP information.",
        as.character(bpid[!is.na(bpid)])[1L]
      )
    ))
  }

  # Requirement 5.1: group by normalized BPId and collect every distinct TPPId
  # under the matched BP.
  tpps <- distinct_tpps(bp_rows)

  built_list <- list()
  for (i in seq_len(nrow(tpps))) {
    tpp_code <- tpps$TPPCode[i]
    tpp_rows <- bp_rows[bp_rows$TPPCode == tpp_code, , drop = FALSE]

    # Requirement 5.3: build the full object with NO environment subsetting.
    # build_tpp_metadata emits env_filters = NULL -> every trait keeps all envs.
    built_list[[tpp_code]] <- list(
      tpp_id   = tpp_code,
      data     = build_tpp_data_sheet(tpp_rows),
      metadata = build_tpp_metadata(tpp_rows, tr_mapping, pheno_designations)
    )
  }

  list(
    status     = "populated",
    tpps       = built_list,
    tpp_codes  = tpps$TPPCode,
    diagnostic = NA_character_
  )
}
