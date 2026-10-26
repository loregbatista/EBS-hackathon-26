# Auto-generated file. Do not edit directly.
# Source scripts are maintained in modular files under scripts/.
# Generated on: 2026-10-26 17:06:20

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
vbo_json_payload <- function(useCases, ready, inputs = NULL) {
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
      checks          = checksOut
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
validate_bioflow_object <- function(object, verbose = TRUE, diagnosticsFile = NULL, inputs = NULL) {

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

  # Write the machine-readable, engineer-friendly diagnostics JSON.
  if (!is.null(diagnosticsFile) && nzchar(diagnosticsFile)) {
    payload <- vbo_json_payload(useCases, ready, inputs = inputs)
    write_diagnostics_json(payload, diagnosticsFile)
  }

  invisible(useCases)
}
# ---- END: validate_bioflow_object.R ----

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
#' @return Invisibly, a list with `result` (the object written), `file` (the
#'   path it was written to) and `diagnosticsFile` (the diagnostics JSON path).
# -------------------------------------------------------------------------------------

getBioflowRData <- function(phenotypeFile = NULL, pedigreeFile = NULL, genotypeFile = NULL,
                     traits = NULL, outputPath, outputFile = NULL, requestId = NULL,
                     pedigreeMapping = NULL, deriveCrossType = TRUE,
                     runSOA = NULL, validate = TRUE, diagnosticsFile = NULL) {

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
    validate_bioflow_object(result, diagnosticsFile = diagnosticsFile, inputs = diagInputs)
  } else {
    # Even when the console report is suppressed, still emit the JSON so the
    # EBS engineer always gets a machine-readable diagnostics artifact.
    validate_bioflow_object(result, verbose = FALSE,
                            diagnosticsFile = diagnosticsFile, inputs = diagInputs)
  }

  invisible(list(result = result, file = outFilePath, diagnosticsFile = diagnosticsFile))

} # end of staRData fxn
# ---- END: getBioflowRdata.R ----

