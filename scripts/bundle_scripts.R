bundle_bioflow_scripts <- function(
  input_files = c(
    "packages_verification.R",
    "read_geno_functions.R",
    "read_pedigree_functions.R",
    "validate_bioflow_object.R",
    "getBioflowRdata.R"
  ),
  scripts_dir = "scripts",
  output_file = "bundled_getBioflowRdata.R"
) {
  scripts_dir <- normalizePath(scripts_dir, mustWork = TRUE)

  input_paths <- file.path(scripts_dir, input_files)
  missing_files <- input_paths[!file.exists(input_paths)]
  if (length(missing_files) > 0) {
    stop(
      sprintf(
        "Missing input script(s): %s",
        paste(basename(missing_files), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  output_path <- file.path(scripts_dir, output_file)

  header <- c(
    "# Auto-generated file. Do not edit directly.",
    "# Source scripts are maintained in modular files under scripts/.",
    sprintf("# Generated on: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    ""
  )

  blocks <- lapply(seq_along(input_paths), function(i) {
    block_title <- sprintf("# ---- BEGIN: %s ----", basename(input_paths[i]))
    block_footer <- sprintf("# ---- END: %s ----", basename(input_paths[i]))
    content <- readLines(input_paths[i], warn = FALSE, encoding = "UTF-8")
    c(block_title, content, block_footer, "")
  })

  bundled_content <- c(header, unlist(blocks, use.names = FALSE))
  writeLines(bundled_content, output_path, useBytes = TRUE)

  message(sprintf("Bundled script generated: %s", output_path))
  invisible(output_path)
}

if (identical(environment(), globalenv()) && !interactive()) {
  # Wrapped in local() so that sourcing this file does not leak `scripts_dir` /
  # `output_file` into the caller's global environment - doing so silently
  # overwrote same-named variables in scripts that source() this one.
  local({
    args <- commandArgs(trailingOnly = TRUE)

    scripts_dir <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "scripts"
    output_file <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "bundled_getBioflowRdata.R"

    # Only bundle when this file is the script being run, not when it is merely
    # sourced for its function definitions.
    invoked <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    if (length(invoked) == 0 || identical(basename(invoked[1]), "bundle_scripts.R")) {
      bundle_bioflow_scripts(
        scripts_dir = scripts_dir,
        output_file = output_file
      )
    }
  })
}
