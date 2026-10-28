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
# Name             : tpp_id_helpers
# Description      : Shared ID-normalization helper used for BP (Breeding
#                    Pipeline) and MS (Market Segment) matching between the
#                    phenotype data and the breeding-portal TPP file.
# Maintainer       : Bioflow core development team
# ---------------------------------------------------------------------------------------------
#
# WHY THIS FILE EXISTS
# --------------------
# BP and MS identifiers appear in several equivalent forms across the data:
#
#   * prefixed  : "DBP0003" (pheno bpid)
#   * coded     : "BP00003" / "MS00054" (TPP-file BPCode / MSShortName)
#   * numeric   : 3 / "3"  (TPP-file BPId / MSId)
#
# All of these denote the same underlying ID. Matching pheno `bpid`/`msid`
# against the TPP file's `BPId`/`MSId` therefore requires reconciling them to a
# single canonical key. `normalize_id` is that canonical reconciliation: it
# strips every non-digit character and then drops leading zeros, yielding an
# integer-like key.
#
# IMPORTANT - TWO-REPO CONSISTENCY
# --------------------------------
# An identical copy of `normalize_id` lives on the bioflow side at
#   bioflow/R/fct_tpp_id_helpers.R
# so that bioflow-side BP/MS matching reconciles IDs exactly the way this
# EBS-side copy does. The two copies MUST be kept byte-for-byte equivalent in
# behavior. If you change the normalization rule here, change it there too.

#' Normalize a BP / MS / TR identifier to a canonical integer-like key
#'
#' Reconciles coded ("BP00003"), prefixed ("DBP0003") and numeric (3, "3")
#' forms of the same identifier to a single canonical key by:
#'   1. coercing the input to character,
#'   2. removing every non-digit character,
#'   3. dropping leading zeros,
#'   4. returning the resulting integer key.
#'
#' The function is vectorized and type-stable: it always returns an integer
#' vector the same length as the input. Any element that is `NA`, empty, or
#' contains no digits normalizes to `NA_integer_` (treated downstream as
#' "no id"). Because the result is an integer, the function is idempotent:
#' `normalize_id(normalize_id(x))` equals `normalize_id(x)`.
#'
#' @param x A vector (character, numeric, integer, factor, or logical) of
#'   identifiers to normalize. Factors are coerced via their labels.
#'
#' @return An integer vector of the same length as `x`. `NA`, empty, and
#'   non-digit-only inputs become `NA_integer_`.
#'
#' @examples
#' normalize_id("DBP0003")            # 3
#' normalize_id(c("BP00003", "3", 3)) # 3 3 3
#' normalize_id(c(NA, "", "abc"))     # NA NA NA
#' @export
normalize_id <- function(x) {
  # Factors carry their labels, not the integer codes we want to read.
  if (is.factor(x)) {
    x <- as.character(x)
  }

  # Coerce to character so a numeric 3 and a character "3" travel the same path.
  chr <- as.character(x)

  # Strip every non-digit character: "DBP0003" -> "0003", "BP00003" -> "00003".
  digits_only <- gsub("[^0-9]", "", chr)

  # Drop leading zeros: "0003" -> "3", "00000" -> "". A value with no digits at
  # all (originally NA, empty, or non-digit-only) also lands here as "".
  trimmed <- sub("^0+", "", digits_only)

  # Anything that reduced to an empty string carries no usable id -> NA. We also
  # preserve NA from the original input (as.character(NA) -> NA -> ""), so every
  # no-id case funnels to NA_integer_.
  out <- suppressWarnings(as.integer(trimmed))
  out[is.na(chr) | trimmed == ""] <- NA_integer_

  out
}
