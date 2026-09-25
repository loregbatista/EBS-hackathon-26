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