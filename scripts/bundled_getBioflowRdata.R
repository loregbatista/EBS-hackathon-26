# Auto-generated file. Do not edit directly.
# Source scripts are maintained in modular files under scripts/.
# Generated on: 2026-10-27 13:47:19

# ---- BEGIN: packages_verification.R ----
# CRAN packages the converter needs. vcfR/adegenet read and hold the markers;
# stringr/purrr/glue/Matrix are used by read_geno_functions.R; openssl generates
# the md5 output file name; cli/rlang provide the messaging and input checks.
BIOFLOW_CRAN_DEPENDENCIES <- c(
  "vcfR", "adegenet", "cli", "rlang", "remotes",
  "openssl", "stringr", "purrr", "glue", "Matrix"
)

ensure_cran_packages <- function(packages = BIOFLOW_CRAN_DEPENDENCIES,
                                repos = "https://cloud.r-project.org") {
  missing <- packages[!vapply(
    packages,
    function(pkg) requireNamespace(pkg, quietly = TRUE),
    FUN.VALUE = logical(1)
  )]

  if (length(missing) > 0) {
    tryCatch(
      install.packages(missing, repos = repos),
      error = function(e) {
        stop(
          sprintf("Failed to install required CRAN packages: %s", paste(missing, collapse = ", ")),
          call. = FALSE
        )
      }
    )
  }

  still_missing <- packages[!vapply(
    packages,
    function(pkg) requireNamespace(pkg, quietly = TRUE),
    FUN.VALUE = logical(1)
  )]

  if (length(still_missing) > 0) {
    stop(
      sprintf("Required CRAN packages are still unavailable: %s", paste(still_missing, collapse = ", ")),
      call. = FALSE
    )
  }
}

ensure_github_packages <- function(){
  # requireNamespace avoids attaching the packages (and the noisy warning that
  # require() emits when a package is absent) while still triggering install.
  if (!requireNamespace("cgiarPipeline", quietly = TRUE)) {
    remotes::install_github("Breeding-Analytics/cgiarPipeline")
  }
  if (!requireNamespace("cgiarBase", quietly = TRUE)) {
    remotes::install_github("Breeding-Analytics/cgiarBase")
  }

  still_missing <- c("cgiarPipeline", "cgiarBase")[
    !vapply(c("cgiarPipeline", "cgiarBase"),
            function(pkg) requireNamespace(pkg, quietly = TRUE),
            FUN.VALUE = logical(1))
  ]
  if (length(still_missing) > 0) {
    stop(
      sprintf("Required Breeding-Analytics packages are unavailable: %s",
              paste(still_missing, collapse = ", ")),
      call. = FALSE
    )
  }
}
# ---- END: packages_verification.R ----

# ---- BEGIN: read_geno_functions.R ----
#' @param geno_metadata Dataframe with the expected columns id, chrom, pos, ref, alt
#'
#' @return Same input dataframe but added a `filter` column if true the marker don't meet
#' the minimal requirements: no duplicated markers, physical position, bi-allelic, 
#'  allelic information, SNP.
#' @export
#'
#' @examples
#' geno_metadata(meta_daya)
process_metadata <- function(geno_metadata){
  # Expected fields id Chrom pos ref alt
  meta_columns <- c('id', 'chrom', 'pos', 'ref', 'alt')
  
  # Check if the columns in the input file match the expected column names
  if (length(intersect(meta_columns, colnames(geno_metadata))) != 5) {
    cli::cli_abort("metadata file doesn't have the expected column names,
                   have: {colnames(geno_metadata)}")
  }
  
  geno_metadata['filter'] <- FALSE
  
  # Flag markers with physical location
  geno_usable_idx <- which(!complete.cases(geno_metadata[,c('chrom','pos')]),)
  #geno_metadata[geno_usable_idx,'filter'] <- TRUE
  no_pos <- length(geno_usable_idx)
  
  if (no_pos > 0){
    cli::cli_inform("Were found {no_pos} markers
                   without physicall location, put into unk chrom and consecutive position")
    
    geno_metadata[geno_usable_idx, 'chrom'] <- "unk"
    geno_metadata[geno_usable_idx, 'pos'] <- seq(1, no_pos)
  }
  
  # Flag colocalized markers
  
  dup_pos_idx <- which(duplicated(geno_metadata[,c('chrom','pos')]))
  no_dup_pos <- length(dup_pos_idx)
  if (no_dup_pos > 0 ){
    cli::cli_inform("Were found {no_dup_pos} markers with duplicated position")
    geno_metadata[dup_pos_idx,'filter'] <- TRUE
  }
  
  # Check duplicated ids
  id_dup <- duplicated(geno_metadata$id)
  
  if(sum(id_dup) > 0){
    cli::cli_inform("Were found {sum(id_dup)} markers with duplicated id. \n
                  Renamed with chrom_pos nomenclature")  
  }
  
  # Count commas in alt column to get allele count
  allele_count <- nchar(geno_metadata$alt) - nchar(gsub(",", "", geno_metadata$alt)) + 1
  
  # Check if ref matches single ACGT pattern
  ref_is_acgt <- grepl("^[ACGT]$", geno_metadata$ref)
  ref_len <- ifelse(ref_is_acgt, nchar(geno_metadata$ref), NA)
  
  # Check if alt matches single ACGT pattern
  alt_is_acgt <- grepl("^[ACGT]$", geno_metadata$alt)
  alt_len <- ifelse(alt_is_acgt, nchar(gsub(",", "", geno_metadata$alt)), NA)
  
  # Update id for duplicates
  row_num <- seq_len(nrow(geno_metadata))
  dup_row_idx <- which(id_dup)
  geno_metadata$id[dup_row_idx] <- paste0(geno_metadata$chrom[dup_row_idx], "_", geno_metadata$pos[dup_row_idx])
  
  # Add new columns
  geno_metadata$allele_count <- allele_count
  geno_metadata$ref_len <- ref_len
  geno_metadata$alt_len <- alt_len
  
  # Flag no reference alleles
  ref_alt_na_idx <- which(is.na(geno_metadata$ref_len) | is.na(geno_metadata$alt_len))
  no_ref_alt <- length(ref_alt_na_idx)
  
  if(no_ref_alt > 0){
    cli::cli_inform("Were found {no_ref_alt} without ref and alt allele data")
    geno_metadata[ref_alt_na_idx,'filter'] <- TRUE
  }
  
  # flag multiallelic 
  multi_allelic_idx <- which(geno_metadata$allele_count > 1)
  multi_allelic <- length(multi_allelic_idx)
  
  if(multi_allelic > 0){
    cli::cli_warn("Were found {multi_allelic} multi-allelic markers")
    geno_metadata[multi_allelic_idx,'filter'] <- TRUE
  }

  # flag indels
  indel_idx <- which(geno_metadata$ref_len > 1)
  indel_idx <- unique(c(indel_idx, which(geno_metadata$alt_len > 1)))
  
  indel <- length(indel_idx)
  
  if(indel > 0){
    cli::cli_inform("Were found {indel} indel markers")
    geno_metadata[indel_idx,'filter'] <- TRUE
  }
  
  return(geno_metadata)
}

#' From locus genotype call data, get the allelic dosage given the alleles and ploidity
#'
#' This function takes a list of genotype calls, a named vector of allele counts,
#' and the ploidity level as input, and returns a list of allelic dosages for the
#' genotype calls. The allelic dosage is the count of the alternative allele in
#' the genotype call.
#'
#' The function first generates all possible genotype calls for the given alleles
#' and ploidity level using the `get_all_gt_calls` function. It then calculates
#' the allelic dosages for these possible genotype calls using the `convert_gt_to_dosage`
#' function, treating the second allele as the alternative allele.
#'
#' Finally, the function replaces the genotype calls in the input list with their
#' corresponding allelic dosages using the `replace_strings_with_integers` function.
#'
#' If a genotype call in the input list is not found in the set of possible genotype
#' calls, its allelic dosage will be set to NA.
#'
#' @param l A list of genotype calls, e.g., c("AG", "GG", "AA").
#' @param alleles ref and alternative allele, e.g., c("A", "C").
#' @param ploidity Integer. The ploidity level of the organism.
#'
#' @return A list of integers, representing the allelic dosages for the input
#'         genotype calls.
#' @export
#'
#' @examples
#' genotypes <- c("AG", "GG", "AA")
#' allele_counts <- c(A = 10, G = 20)
#' get_allelic_dosage(genotypes, allele_counts, 2)  # Returns list(1, 2, 0)
#'
#' # Example with missing genotype call
#' genotypes <- c("AG", "XX", "AA")
#' allele_counts <- c(A = 10, G = 20)
#' get_allelic_dosage(genotypes, allele_counts, 2)  # Returns list(1, NA, 0)
get_allelic_dosage <- function(l, alleles, ploidity, sep = "") {
  alleles_c <- unlist(stringr::str_split(alleles, "/"))
  # All possible genotype calls
  possible_gt_calls <- get_all_gt_calls(alleles_c, ploidity, sep)
  # Get dosage given the alternative allele
  possible_dosage <- convert_gt_to_dosage(possible_gt_calls, alleles_c[2], ploidity,sep)
  dosages <- replace_strings_with_integers(possible_dosage, l)
  return(dosages)
}

#' Get all possible genotype calls given a unique set of alleles
#'
#' This function generates all possible genotype calls for a given set of alleles
#' and ploidity level. The genotype calls are represented as strings of characters,
#' with each allele being a single character.
#'
#' The function uses a recursive approach to generate all possible combinations of
#' alleles for the specified ploidity level. For example, with two alleles "A" and "B",
#' and a ploidity of 2 (diploid), the function would generate the following genotype
#' calls: "AA", "AB", "BA", "BB".
#'
#' @param alleles List[String]. A list of unique alleles, e.g., c("A", "B", "C").
#' @param ploidity Integer. The ploidity level of the organism.
#'
#' @return List[String]. A list of all possible genotype calls for the given alleles
#'         and ploidity level.
#' @export
#'
#' @examples
#' get_all_gt_calls(c("A", "B"), 2)  # Returns c("AA", "AB", "BA", "BB")
#' get_all_gt_calls(c("A", "B", "C"), 3)  # Returns all 27 possible triploid calls
#' get_all_gt_calls(c("A"), 1)  # Returns c("A")
get_all_gt_calls <- function(alleles, ploidity, sep = "") {
  generate_calls <- function(prefix, ploidity, del = sep) {
    if (ploidity == 0) {
      if(nchar(sep) > 0){
        out <- substr(prefix, 1, nchar(prefix)-1)  
      } else {
        out <- substr(prefix, 1, nchar(prefix))  
      }
      
      return(out)
    }
    calls <- c()
    for (allele in alleles) {
      call <- paste0(prefix, allele,del)
      calls <- c(calls, generate_calls(call, ploidity - 1))
    }
    return(calls)
  }
  generate_calls("", ploidity)
}

#' Given a list of genotype calls, get the dosage of each one
#'
#' This function takes a list of genotype calls, an alternative allele, and the ploidity level
#' of the organism as input, and returns a list of allelic dosages corresponding to each
#' genotype call in the input list.
#'
#' @param locus List. A list of genotype calls, e.g., c("AG", "GG", "AA").
#' @param alt_allele String. The alternative allele, e.g., "A", "G".
#' @param ploidity Integer. The ploidity level of the organism, default is 2 (diploid).
#'
#' @return A list of integers, representing the allelic dosages of the alternative allele
#'         for each genotype call in the input list.
#' @export
#'
#' @examples
#' convert_gt_to_dosage(c("AG", "GG", "AA"), "A")  # Returns list(1, 0, 2)
#' convert_gt_to_dosage(c("AG", "GG", "AA"), "G")  # Returns list(1, 2, 0)
#' convert_gt_to_dosage(c("AAA", "GGG"), "A", 3)  # Returns list(3, 0) (triploid)
#' convert_gt_to_dosage(c(NA, "AG"), "A")  # Returns list(NA, 1)
convert_gt_to_dosage <- function(locus, alt_allele, ploidity = 2,sep="") {
  l <- sapply(locus,
              genocall_to_allelic_dosage,
              alt_allele = alt_allele,
              ploidity = ploidity,
              sep=sep)
  return(l)
}


#' Genotype call to allelic dosage of alternative allele
#'
#' This function takes a genotype call, an alternative allele, and the ploidity level
#' of the organism as input, and returns the allelic dosage of the alternative allele
#' in the genotype call.
#'
#' The genotype call is expected to be a string of characters representing the alleles,
#' with each allele being a single character. For example, "AG" represents a diploid
#' genotype with one allele being "A" and the other being "G".
#'
#' The allelic dosage is the count of the alternative allele in the genotype call.
#' For example, if the genotype call is "AG" and the alternative allele is "A", the
#' allelic dosage would be 1.
#'
#' If the genotype call is missing (represented as NA or an empty string), the
#' function returns NA.
#'
#' @param genotype_call String. Genotype call, e.g., "AG", "AAA" (for triploid).
#' @param alt_allele String. Alternative allele, e.g., "A", "G".
#' @param ploidity Integer. Ploidity level of the organism, default is 2 (diploid).
#'
#' @return Integer. The allelic dosage of the given genotype call for the alternative allele.
#' @export
#'
#' @examples
#' genocall_to_allelic_dosage("AG", "A")  # Returns 1
#' genocall_to_allelic_dosage("GG", "A")  # Returns 0
#' genocall_to_allelic_dosage("AAA", "A", 3)  # Returns 3 (triploid)
#' genocall_to_allelic_dosage(NA, "A")  # Returns NA
genocall_to_allelic_dosage <- function(genotype_call, alt_allele, ploidity = 2,sep="") {
  if (!is.na(nchar(genotype_call))) {
    # remove separators (and normalize phasing if present)
    genotype_call <- gsub("\\|", sep, genotype_call)
    if(sep != ""){
      genotype_call <- gsub(sep, "", genotype_call, fixed = TRUE)
    }
    # Genotype call successfully genotyped
    allele_length <- nchar(genotype_call) / ploidity
    
    # List with each allele as element
    split_genotype <- substring(genotype_call, seq(1, nchar(genotype_call), allele_length),
                                seq(allele_length, nchar(genotype_call), allele_length))
    
    # Matches of alt allele are the dosage
    dosage <- length(which(split_genotype == alt_allele))
    return(dosage)
  } else {
    # Genotype call missed
    return(NA)
  }
}

#' Replace a list of strings with their corresponding integer values
#'
#' This function takes two inputs: a named list or vector with string keys and integer values,
#' and a list of strings to be replaced. It replaces each string in the second list with the
#' corresponding integer value from the first list, based on the string keys.
#'
#' If a string in the second list does not have a corresponding key in the first list,
#' it will be replaced with NA.
#'
#' @param lookup_table A named list or vector with string keys and integer values.
#' @param strings_to_replace A list of strings to be replaced with their corresponding integer values.
#'
#' @return A list of integers, where each string in the input list has been replaced with its
#'         corresponding integer value from the lookup table, or NA if no match was found.
#' @export
#'
#' @examples
#' lookup <- c(A = 1, B = 2, C = 3)
#' strings <- c("B", "A", "D", "C")
#' replace_strings_with_integers(lookup, strings)  # Returns list(2, 1, NA, 3)
replace_strings_with_integers <- function(lookup_table, strings_to_replace) {
  # Use match to find the indices of the strings in the lookup table
  indices <- match(strings_to_replace, names(lookup_table))
  
  # Replace the strings with the corresponding integer values
  # or NA if no match was found
  integer_values <- lookup_table[indices]
  
  return(integer_values)
}

get_loc_missing <- function(gl) {
  # Get the number of occurrences of NAs (missing data) for each marker
  NA_counts <- adegenet::glNA(gl)
  
  # Divide the NA counts by the total number of samples to get the missing rate
  NA_counts <- NA_counts / adegenet::nInd(gl)/max(adegenet::ploidy(gl))
  
  return(NA_counts)
}

#' Get Locus Missingness
#'
#' Compute the missing rate for each locus (marker) in a genlight object.
#' The missing rate is a value between 0 and 1, where 0 indicates no missing data
#' for that locus, and 1 indicates that all samples have missing data for that locus.
#'
#' @param gl A genlight object.
#'
#' @return A numeric vector of length equal to the number of loci (markers),
#'   containing the missing rate for each locus.
#'
#' @export
#'
#' @examples
#' data(example_genlight)
#' loc_miss <- get_loc_missing(example_genlight)
#' head(loc_miss)
get_ind_missing <- function(gl) {
  # Convert the genlight object to a matrix and identify missing genotype calls
  mt <- is.na(as.matrix(gl))
  
  # Calculate the proportion of missing data for each individual (row)
  ind_miss <- Matrix::rowSums(mt) / adegenet::nLoc(gl)
  
  return(ind_miss)
}

#' Get Overall Missingness
#'
#' Compute the overall missing rate for a genlight object.
#' The overall missing rate is the proportion of missing genotype calls
#' across all individuals and loci in the dataset.
#'
#' @param gl A genlight object.
#'
#' @return A single numeric value representing the overall missing rate.
#'
#' @export
#'
#' @examples
#' data(example_genlight)
#' overall_miss <- get_overall_missingness(example_genlight)
#' print(overall_miss)
get_overall_missingness <- function(gl) {
  # Convert the genlight object to a matrix
  mt <- as.matrix(gl)
  
  # Identify missing genotype calls
  mt <- is.na(mt)
  
  # Calculate the overall missing rate
  overall_miss <- sum(mt) / (nrow(mt) * ncol(mt))
  
  return(overall_miss)
}

get_heterozygosity_metrics <- function(gl, ploidy = 2){
  # Boolean matrix of genotype calls where 0 > dosage < ploidy
  mt <- as.matrix(gl)
  het_ind_loc <- mt > 0 & mt < ploidy
  het_loc <- Matrix::colSums(het_ind_loc , na.rm = T)/adegenet::nInd(gl)
  het_ind <- Matrix::rowSums(het_ind_loc , na.rm = T)/adegenet::nLoc(gl)
  return(list(het_ind = het_ind, het_loc = het_loc))
}

get_loc_heterozygosity <- function(gl, ploidy = 2){
  mt <- as.matrix(gl)
  het_ind_loc <- mt > 0 & mt < ploidy
  het_loc <- Matrix::colSums(het_ind_loc , na.rm = T)/adegenet::nInd(gl)
  return(het_loc)
}

get_ind_heterozygosity <- function(gl, ploidy = 2){
  mt <- as.matrix(gl)
  het_ind_loc <- mt > 0 & mt < ploidy
  het_ind <- Matrix::rowSums(het_ind_loc , na.rm = T)/adegenet::nLoc(gl)
  return(het_ind)
}

get_maf <- function(gl, ploidy = 2){
  alf <- adegenet::glMean(gl)*(1/ploidy)
  maf <- ifelse(alf > 0.5, 1 - alf, alf)
  return(maf)
}

get_inbreeding <- function(gl){
  p <- 1 - gl@other$loc.metrics$maf
  he <- 2 * gl@other$loc.metrics$maf * p
  Fis <- 1 - (gl@other$loc.metrics$loc_het/he)
  return(Fis)
}

recalc_metrics <- function(gl){

  gl@other$loc.metrics <- data.frame(
    maf = get_maf(gl, max(adegenet::ploidy(gl))),
    loc_miss = get_loc_missing(gl),
    loc_het = get_loc_heterozygosity(gl, max(adegenet::ploidy(gl)))
  )
  
  # Use already calculated loc stats
  gl@other$loc.metrics$loc_Fis <- get_inbreeding(gl)
  
  gl@other$ind.metrics <- data.frame(
                                     ind_miss = get_ind_missing(gl),
                                     ind_het = get_ind_heterozygosity(gl, max(adegenet::ploidy(gl)))
                                     )
  return(gl)
}


get_overall_summary <- function(gl){
  ninds <- adegenet::nInd(gl)
  nlocs <- adegenet::nLoc(gl)
  overall_missiness <- mean(gl@other$ind.metrics$ind_miss)
  overall_heterozygosity <- mean(gl@other$ind.metrics$ind_het)
  overall_maf <- mean(gl@other$loc.metrics$maf, na.rm = T)
  
  out <- list(
    nind = ninds,
    nloc = nlocs,
    ov_miss = overall_missiness,
    ov_het = overall_heterozygosity,
    ov_maf = overall_maf
  )
  return(out)
}

# Function to infer ploidy from a single genotype call
get_ploidy_from_gt <- function(gt_string) {
  if (is.na(gt_string)) return(NA)
  # Count separators (/ or |) and add 1
  separators <- nchar(gt_string) - nchar(gsub("[/|]", "", gt_string))
  return(separators + 1)
}

#' read_vcf
#'
#' This function reads a VCF file (compressed or uncompressed) and converts it into a genlight object.
#'
#' @param path String. Path to the VCF file. It could be compressed.
#' @param ploidity Integer. Ploidity level of the organism. (Default = 2)
#' @param na_reps Vector. A vector containing the NA representations of genotype calls (default: empty).
#'
#' @return A genlight object.
#' @export
#'
#' @examples
#' fl = "https://github.com/Breeding-Analytics/cgiarGenomics/raw/main/tests/vcf_fmt/diploid.vcf.gz"
#' tempfl <- tempfile(pattern = 'diploid', fileext = '.vcf.gz')
#' download.file(fl, destfile = tempfl)
#' dat.dose.vcf = read_vcf(tempfl, ploidity = 2)
#' print(dat.dose.vcf)
#' plot(dat.dose.vcf)
read_vcf <- function(path, na_reps = c("-", "./."), sep="/") {
  
  if (!file.exists(path)){
    cli::cli_abort("`path` don't exist. Verify if is writed properly {path}")
  }
  # Read the VCF file
  vcf <- vcfR::read.vcfR(path)
  
  gt_matrix <- vcfR::extract.gt(vcf, return.alleles = TRUE)
  ploidy_matrix <- apply(gt_matrix, c(1, 2), get_ploidy_from_gt)
  # Get the most common ploidy level
  ploidity <- as.numeric(names(sort(table(ploidy_matrix), decreasing = TRUE))[1])

  # Get the metadata from the VCF file
  meta_vcf <- as.data.frame(vcfR::getFIX(vcf))
  
  # Rename columns and select only the needed ones
  meta <- data.frame(
    id = meta_vcf$ID,
    chrom = meta_vcf$CHROM,
    pos = meta_vcf$POS,
    ref = meta_vcf$REF,
    alt = meta_vcf$ALT,
    stringsAsFactors = FALSE
  )
  
  meta <- process_metadata(meta)
  mt <- t(vcfR::extract.gt(vcf, return.alleles = TRUE)[!meta$filter,])
  if (length(na_reps) > 0) {
    idx <- which(mt %in% na_reps)
    mt[idx] <- NA
  }
  
  # check ploidity
  mt_gt_str <- matrix(gsub(sep, "", mt), nrow = dim(mt)[1], ncol = dim(mt)[2])
  gc_len <- apply(mt_gt_str, 2, function(x) max(nchar(x), na.rm = TRUE))
  max_dosage <- max(gc_len, na.rm = TRUE)
  
  if(max_dosage < ploidity){
    cli::cli_warn("Max dosage ({max_dosage}) lower than ploidy lvl ({ploidity})")
  }
  
  if (max_dosage > ploidity){
    cli::cli_abort("Max dosage ({max_dosage}) higher than ploidy lvl ({ploidity})")
  }
  
  individuals <- rownames(mt)
  
  allele_set <- paste(meta$ref[!meta$filter], meta$alt[!meta$filter], sep=sep)

  gt <- mapply(function(col, arg, ploidity, sep) get_allelic_dosage(mt[,col], arg, ploidity, sep),
               col = seq(1, dim(mt)[2]), 
               arg = allele_set,
               ploidity = ploidity,
               sep = sep)
  
  gl <- new("genlight",
            gt,
            ploidy = ploidity,
            loc.names = meta$id[!meta$filter],
            ind.names = individuals,
            chromosome = meta$chrom[!meta$filter],
            position = meta$pos[!meta$filter])
  adegenet::alleles(gl) <- allele_set
  gl <- recalc_metrics(gl)
  return(gl)
}

#' Filter function
#' 
#' This function generalizes the filtering functions using the parameter name
#' and comparing in locus or individuals using the provided comparision operator.
#' A list indicating the used thresold, if the filter was performed over individuals
#' or locus and the indices of elements that meet the comparision.
#'
#' @param gl 
#' @param parameter 
#' @param threshold 
#' @param comparison_operator 
#'
#' @return
#' @export
#'
#' @examples
filter_gl <- function(gl, parameter, threshold, comparison_operator){
  
  # Verify if parameter exist on the gl and get the margin (ind, loc)
  filter_margin <- get_parameter_margin(gl, parameter)
  comparison_operator <- match.arg(comparison_operator, choices = c(">", ">=", "<", "<="))
  comparison_func <- match.fun(comparison_operator)
  # Verify if threshold is a float value
  if(threshold > 1 ){
    cli::cli_abort("`threshold`: {threshold} is greather of equal to 1, correct it")
  }
  
  if(filter_margin == 'loc'){
    index <- which(comparison_func(gl@other$loc.metrics[parameter], threshold))
    filter_out <- gl@loc.names[-c(index)]
    
  } else {
    index <- which(comparison_func(gl@other$ind.metrics[parameter], threshold))
    filter_out <- gl@ind.names[-c(index)]
    
  }
  
  out <- list(param = parameter,
              operator = comparison_operator,
              threshold = threshold,
              filter_margin = filter_margin,
              index = index,
              filter_out = filter_out)
  
  return(out)
}

get_parameter_margin <- function(gl, param_name){
  # Get the expected parameters
  loc_metric_names <- colnames(gl@other$loc.metrics)
  ind_metric_names <- colnames(gl@other$ind.metrics)
  
  param_name = match.arg(param_name, choices = c(loc_metric_names, ind_metric_names))
  
  if(param_name %in% loc_metric_names){
    by = 'loc'
  } else {
    by = 'ind'
  }
  return(by)
}

#' Apply a sequence of filterings over a gl object
#'
#' filt_sequence named list (param = param_name, threshold: t, operator: op)
#' @param gl 
#' @param filt_sequence 
#'
#' @return
#' @export
#'
#' @examples
apply_sequence_filtering <- function(gl, filt_sequence){
  if(!rlang::is_bare_list(filt_sequence)){
    cli::cli_abort("Provide a list of filter operations in `filt_sequence`")
  }
  if(!inherits(gl, "genlight")){
    cli::cli_abort("`gl` is not a genlight class")
  }
  
  # Allways add at the end a filter step to remove all NA loci and ind
  
  locNA_filt <- list("loc_miss", "<", 1)
  indNA_filt <- list("ind_miss", "<", 1)
  
  
  allNA_steps <- list(
    locNA_filt,
    indNA_filt
  )
  
  
  filt_NA_seq <- lapply(allNA_steps, function(x){
    setNames(as.list(x), c("param", "operator", "threshold"))
  })
  filt_sequence <- append(filt_sequence, filt_NA_seq)
  # Duplicate the gl object to perform the filtering
  working_gl <- gl
  filtering_log <- list()
  previous_margin <- ""
  
  for (i_step in 1:length(filt_sequence)){
    filt_step <- filt_sequence[[i_step]]
    param <- filt_step[['param']]
    threshold <- filt_step[['threshold']]
    operator <- filt_step[['operator']]
  
    i_filt_out <- filter_gl(working_gl,
                            parameter = param,
                            threshold = threshold,
                            comparison_operator = operator)
    
    i_margin <- get_parameter_margin(gl, param)
    
    if(length(i_filt_out$index) > 0){
      if(i_margin == "loc"){
        working_gl <- working_gl[,i_filt_out$index]
      } else {
        working_gl <- working_gl[i_filt_out$index,]
      }
    }
    working_gl <- recalc_metrics(working_gl)
    filtering_log[[glue::glue("{param}_{i_step}")]] <- i_filt_out
  }
  
  return(list(gl = working_gl, filt_log = filtering_log))
}

get_filter_log <- function(filter_step_log, geno_data){
  print("filter_processing...")
  base_loc_names <- adegenet::locNames(geno_data)
  base_ind_names <- adegenet::indNames(geno_data)
  out <- purrr::map_df(filter_step_log, function(filter_step){
    
    if(length(filter_step$filter_out) > 0){
      
      reason <- paste(filter_step$filter_margin,
                      filter_step$param,
                      filter_step$operator,
                      filter_step$threshold,
                      sep = '_')
      
      if(filter_step$filter_margin == 'loc'){
        loc_idx <- which(filter_step$filter_out %in% base_loc_names)
        col_data <- loc_idx
        row_data <- rep(NA, length(loc_idx))
        
      } else {
        ind_idx <- which(filter_step$filter_out %in% base_ind_names)
        col_data <- rep(NA, length(ind_idx))
        row_data <- ind_idx
      }
      
      filt_step_log <- data.frame(
        reason = rep(reason, length(filter_step$filter_out)),
        row = row_data,
        col = col_data,
        value = rep(NA, length(filter_step$filter_out))
      )
      return(filt_step_log)
    }
  })
  return(out)
}

#' Imputation with allele frequency
#' 
#' Assuming a bi-allelic marker, using the observed allelic frequency for one allele
#' is sampled the genotype call for any ploidity level
#'
#' @param q_frq 
#' @param ploidity 
#'
#' @return
#' @export
#'
#' @examples
i_freq_impute <- function(q_frq, ploidity = 2){
  if(!rlang::is_integerish(ploidity)){
    cli::cli_abort("`ploidity` is not an integer: {ploidity}")  
  }
  dosage <- 0
  for(i_chromatid in seq(ploidity)){
    i_dosage <- sample(c(1,0), size = 1, 
                       prob = c(q_frq, 1 - q_frq), replace = T)
    dosage <- dosage + i_dosage
  }
  return(dosage)
}

freq_impute <- function(gl, mt, ploidity){
  mt <- as.matrix(gl)
  # Get the allelic frequencies
  q_allele <- adegenet::glMean(gl)
  # Linear index of nas
  idx_na <- which(is.na(mt))
  na_loc_idx <- sapply(idx_na, function(x){
    loc_idx <- ceiling(x/nrow(mt))
    return(loc_idx)
  })
  
  if(length(na_loc_idx) > 0){
    imp <- unname(unlist(lapply(q_allele[na_loc_idx],
                                function(x) {
                                  return(as.numeric(i_freq_impute(q_frq = x, ploidity)))})))
    return(split(idx_na, imp))
  } else {
    return(NULL)
  }
}

apply_imputation <- function(mt, imp_dict){
  if (!is.null(imp_dict)){
    for (dosage in names(imp_dict)) {
      # Convert the list name to a numeric value
      num_dosage <- as.numeric(dosage)
      
      # Get the linear indices associated with this value
      idx <- imp_dict[[dosage]]
      
      # Assign the value to these positions in the matrix
      mt[idx] <- num_dosage
    }
  }
  return(mt)
}
#' Impute a gl object
#'  
#' Impute a genlight object using frequency or random forest method.
#' Returns a list with imputed genlight object and imputation log.
#'
#' @param gl genlight object
#' @param ploidity ploidy level (default=2)
#' @param method Imputation method: 'frequency' or 'random_forest' (default='frequency')
#' @param nflank Number of flanking markers for RF method (default=100)
#' @param ntree Number of trees for RF method (default=100)
#' @param seed Optional seed for reproducibility (RF method)
#'
#' @return List with elements:
#'   - gl: imputed genlight object
#'   - log: imputation dictionary (imputed positions grouped by dosage)
#' @export
#'
#' @examples
impute_gl <- function(gl, ploidity = 2, method = 'frequency', nflank = 100, ntree = 100, seed = NULL){
  
  loci_all_nas <- adegenet::glNA(gl)/ploidity == adegenet::nInd(gl)
  
  if(sum(loci_all_nas) > 0){
    cli::cli_warn("There are {sum(loci_all_nas)} loci with all missing data")
    # Filter out the all na loci
    all_notna_idxs <- which(!loci_all_nas)
    gl <- gl[,all_notna_idxs]
    mt <- as.matrix(gl)
  }
  
  
  nas_number <- sum(adegenet::glNA(gl))/ploidity
  number_imputations <- nas_number - (sum(loci_all_nas) * adegenet::nInd(gl))
  
  mt <- as.matrix(gl)
  
  
  cli::cli_inform("Missing genotype calls {number_imputations}")
  
  if(method == 'frequency'){
    imp_dict <- freq_impute(gl, mt, ploidity)
  } else if(method == 'random_forest'){
    cli::cli_inform("Imputing with Random Forest (nflank={nflank}, ntree={ntree})")
    imp_dict <- rf_impute(gl, nflank = nflank, ntree = ntree, seed = seed)
  } else {
    cli::cli_abort("Unknown imputation method: {method}. Use 'frequency' or 'random_forest'")
  }
  
  # apply the imputation creating a new gl instance
  imp_mt <- apply_imputation(mt, imp_dict)
  
  
  imp_gl <- new("genlight",
            imp_mt,
            ploidy = ploidity,
            loc.names = gl@loc.names,
            ind.names = gl@ind.names,
            chromosome = gl@chromosome,
            position = gl@position)
  
  adegenet::alleles(imp_gl) <- adegenet::alleles(gl)
  imp_gl <- recalc_metrics(imp_gl)
  
  return(list(gl = imp_gl, log = imp_dict))
}

get_filter_log <- function(ogl,filter_step_log){
  print("filter_processing...")
  base_loc_names <- adegenet::locNames(ogl)
  base_ind_names <- adegenet::indNames(ogl)
  out <- purrr::map_df(filter_step_log, function(filter_step){
    
    if(length(filter_step$filter_out) > 0){
      
      reason <- paste(filter_step$filter_margin,
                      filter_step$param,
                      filter_step$operator,
                      filter_step$threshold,
                      sep = '_')
      
      if(filter_step$filter_margin == 'loc'){
        loc_idx <- which(filter_step$filter_out %in% base_loc_names)
        col_data <- loc_idx
        row_data <- rep(NA, length(loc_idx))
        
      } else {
        ind_idx <- which(filter_step$filter_out %in% base_ind_names)
        col_data <- rep(NA, length(ind_idx))
        row_data <- ind_idx
      }
      
      filt_step_log <- data.frame(
        reason = rep(reason, length(filter_step$filter_out)),
        row = row_data,
        col = col_data,
        value = rep(NA, length(filter_step$filter_out))
      )
      return(filt_step_log)
    }
  })
  return(out)
}
# ---- END: read_geno_functions.R ----

# ---- BEGIN: read_pedigree_functions.R ----
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

#' Align genotype individual names to the pedigree identifier space
#'
#' EBS genotype exports key the VCF samples by sample_id (sample_code), but
#' Bioflow's F1 qa/qc matches parents by `designation` and Pedigree qa/qc
#' matches designation/mother/father triplets by `designation`. When the
#' pedigree provides a sample_id -> designation mapping, the parental samples in
#' the genlight must therefore be renamed to their designation, while the F1
#' progeny are kept under their sample_id (F1 qa/qc looks progeny up by
#' sample_id, and sibling F1s can share a designation).
#'
#' This mirrors bioflow/R/mod_getDataGeno.R::translate_ind_names() so the object
#' this script writes behaves the same as one prepared through the Bioflow UI.
#' Orientation-aware and idempotent: names already in the designation space are
#' left untouched, so a designation-keyed VCF, a sample_id-keyed VCF and a mixed
#' VCF all come out correct.
#'
#' @param gl A genlight object (the loaded genotypes).
#' @param pedigreeData data$pedigree, canonical column names.
#' @return `gl` with parental individual names translated to designations.
align_geno_to_pedigree <- function(gl, pedigreeData) {
  if (!inherits(gl, "genlight")) return(gl)
  if (is.null(pedigreeData) || !("sample_id" %in% colnames(pedigreeData)) ||
      !("designation" %in% colnames(pedigreeData))) {
    return(gl)
  }

  inds <- as.character(adegenet::indNames(gl))

  sid   <- as.character(pedigreeData$sample_id)
  desig <- as.character(pedigreeData$designation)

  # crossType, if present, tells us which rows are F1 progeny (kept as sample_id).
  isF1row <- if ("crossType" %in% colnames(pedigreeData)) {
    toupper(as.character(pedigreeData$crossType)) == "F1"
  } else {
    rep(FALSE, nrow(pedigreeData))
  }

  # Map each genlight name (assumed sample_id) to its pedigree row.
  row_idx <- match(inds, sid)

  final <- inds
  for (i in seq_along(inds)) {
    r <- row_idx[i]
    if (is.na(r)) next                                   # not a sample_id; leave as is
    if (isTRUE(isF1row[r])) next                         # F1 progeny: keep sample_id
    d <- desig[r]
    if (!is.na(d) && nzchar(d) && !identical(inds[i], d)) {
      final[i] <- d                                      # parent: sample_id -> designation
    }
  }

  if (!identical(final, inds)) {
    nTranslated <- sum(final != inds)
    cli::cli_inform(paste0(
      "Aligned ", nTranslated, " parental genotype sample(s) to their ",
      "designation so F1/Pedigree qa/qc can match parents; F1 progeny kept ",
      "under their sample_id."
    ))
    adegenet::indNames(gl) <- final
  }
  gl
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
# ---- END: read_pedigree_functions.R ----

# ---- BEGIN: validate_bioflow_object.R ----
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
    # helper: run metASREML with a given heritability filter; returns
    # list(ok, npred, err). Mirrors tryMet() for LMMsolve so the two engines
    # report on equal footing (default filter first, bypass only if needed).
    # metASREML needs covMod/addG lists aligned to randomTerm ("none" = no kernel).
    tryAsr <- function(bypass) {
      args <- list(phenoDTfile = object, analysisId = soaId, trait = trait1,
                   fixedTerm = list("1"), randomTerm = list("designation"),
                   covMod = list("none"), addG = list("none"), verbose = FALSE)
      if (bypass) { args$heritLB <- 0; args$heritUB <- 1 }
      out <- tryCatch(do.call(cgiarPipeline::metASREML, args),
                      error = function(e) e, warning = function(w) w)
      if (inherits(out, "condition")) return(list(ok = FALSE, npred = 0, err = vbo_err1(out)))
      np <- tryCatch(sum(out$status$module == "moaAsr"), error = function(e) 0)
      list(ok = TRUE, npred = np, err = "")
    }

    # Default filter first, then bypass if it failed or kept no environment.
    d <- tryAsr(bypass = FALSE)
    if (d$ok && d$npred > 0) {
      rec <- vbo_exec(rec, "moa_asreml_main", "MOA main-effects (ASREML)", TRUE, "",
                      "ran with default heritability filter [0.15, 0.95]")
    } else {
      b <- tryAsr(bypass = TRUE)
      if (b$ok) {
        dmsg <- if (d$ok) "default filter kept no environment" else d$err
        rec <- vbo_exec(rec, "moa_asreml_main", "MOA main-effects (ASREML)", TRUE,
                        paste0("Default filter failed/empty (", dmsg, ")"),
                        "ran ONLY after bypassing the heritability filter (heritLB=0, heritUB=1)")
      } else {
        dmsg <- if (d$ok) "kept no environment" else d$err
        rec <- vbo_exec(rec, "moa_asreml_main", "MOA main-effects (ASREML)", FALSE,
                        paste0("default filter: ", dmsg, " | bypassed filter: ", b$err))
      }
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
# ---- END: validate_bioflow_object.R ----

# ---- BEGIN: tpp_id_helpers.R ----
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
# ---- END: tpp_id_helpers.R ----

# ---- BEGIN: tpp_population_helpers.R ----
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
# ---- END: tpp_population_helpers.R ----

# ---- BEGIN: tpp_build_helpers.R ----
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
# ---- END: tpp_build_helpers.R ----

# ---- BEGIN: tpp_gating_helpers.R ----
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
# ---- END: tpp_gating_helpers.R ----

# ---- BEGIN: getBioflowRdata.R ----
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
# ---- END: getBioflowRdata.R ----

