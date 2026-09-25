# Verifies that the generated object can drive a Multi Trial Analysis in
# Bioflow's LMMsolve module using the "Main effects model + GEBV" combination.
#
# In mod_mtaLMMsolveApp.R that combination is radio = "mn_model" (Main effect)
# and radioModel = "geno_model" (GEBV as the surrogate of merit). This script
# reproduces the arguments that module builds for those settings and then runs
# the same three-function pipeline chain the Run button executes:
#
#   cgiarPipeline::premetLMMsolver()  -> gxe model bookkeeping
#   cgiarPipeline::metLMMsolver()     -> the mixed model
#   cgiarPipeline::postmetLMMsolver() -> derived predictions/metrics
#
# Defaults taken from the module:
#   nTermsFixed  = 1, leftSidesFixed1  = "environment"   (mn_model)
#   nTermsRandom = 1, leftSidesRandom1 = "designation"   (mn_model)
#                     rightSidesRandom1 = "genoA"        (geno_model / GEBV)
#   nPC = 0 for every covariate, heritLB/UB = 0.1/0.95, meanLB/UB = 0/1e6,
#   subsetGeno = -1, subsetPed = -1, maxIters = 35, useWeights = TRUE

source("scripts/bundled_getBioflowRdata.R")
source("test/make_coherent_fixture.R")

failures <- character(0)
expect <- function(condition, what) {
  if (isTRUE(condition)) cat("  ok   ", what, "\n", sep = "")
  else { cat("  FAIL ", what, "\n", sep = ""); failures <<- c(failures, what) }
}

cat("\n### Building the object ###\n")
fx <- make_coherent_fixture(outDir = "test")
out <- getBioflowRData(
  phenotypeFile = fx$phenotypeFile,
  pedigreeFile  = fx$pedigreeFile,
  genotypeFile  = fx$genotypeFile,
  traits        = fx$trait,
  outputPath    = "test",
  outputFile    = "coherent_bioflow_input",
  pedigreeMapping = list(crossType = "entry", yearOfOrigin = "year"),
  validate      = FALSE
)
dtMta <- out$result

staId    <- dtMta$status$analysisId[dtMta$status$module == "sta"]
qaGenoId <- as.character(dtMta$status$analysisId[dtMta$status$module == "qaGeno"])
trait    <- fx$trait

cat("  sta stamp    : ", staId, "\n", sep = "")
cat("  qaGeno stamp : ", qaGenoId, "\n", sep = "")
cat("  trait        : ", trait, "\n", sep = "")

# ---------------------------------------------------------------------------
# Preconditions the module itself enforces before enabling the Run button
# ---------------------------------------------------------------------------
cat("\n### Module preconditions ###\n")
expect("sta" %in% dtMta$status$module, "an 'sta' stamp exists (Run button guard)")
mapped <- length(which(c("environment","designation","trait") %in% dtMta$metadata$pheno$parameter))
expect(mapped == 3, "environment/designation/trait are mapped (input panels visible)")
# The GEBV covariate menu only offers genoA when data$geno exists.
covarChoices <- setdiff(names(dtMta$data), c("qtl","genodir","pheno"))
expect("geno" %in% covarChoices, "data$geno present, so 'genoA' is offered as a covariate")
expect(qaGenoId %in% names(dtMta$data$geno_imp),
       "the qaGeno stamp resolves in data$geno_imp (metLMMsolver looks it up by name)")

# ---------------------------------------------------------------------------
# Reproduce the module's argument construction
# ---------------------------------------------------------------------------
fixedTerm     <- list("environment")   # leftSidesFixed1
randomTerm    <- list("designation")   # leftSidesRandom1
expCovariates <- list("genoA")         # rightSidesRandom1, i.e. GEBV

# envsToInclude: the env x trait table, cells <= 1 zeroed, then the same
# apply() the module applies before passing it on.
dtProv <- dtMta$predictions[which(dtMta$predictions$analysisId %in% staId), ]
dtProvTable <- as.data.frame(do.call(rbind, list(with(dtProv, table(environment, trait)))))
bad <- which(dtProvTable <= 1, arr.ind = TRUE)
if (nrow(bad) > 0) dtProvTable[bad] <- 0
dtProvTable[which(dtProvTable > 1, arr.ind = TRUE)] <- 1
myEnvsTI <- if (nrow(dtProvTable) > 1) apply(dtProvTable, 2, function(z) z) else dtProvTable

cat("\n  environments offered to MTA:\n")
print(myEnvsTI)
expect(sum(as.matrix(myEnvsTI)) >= 2, "at least 2 environments are retained")

# traitFamily: the distribution grid defaults to all zeros, so the module falls
# back to the quasi identity family for every trait.
myFamily <- stats::setNames("quasi(link = 'identity',variance='constant')", trait)

# nPC: one entry per covariate choice, all zero by default.
nPCchoices <- c(setdiff(names(dtMta$data), c("qtl","genodir","pheno","geno_imp")),
                unique(dtProv$trait))
nPC <- stats::setNames(rep(0, length(nPCchoices)), nPCchoices)
cat("\n  nPC vector: ", paste(names(nPC), nPC, sep = "=", collapse = ", "), "\n", sep = "")

# ---------------------------------------------------------------------------
cat("\n### 1/3 premetLMMsolver ###\n")
result1 <- try(
  cgiarPipeline::premetLMMsolver(
    phenoDTfile = dtMta, fixedTerm = fixedTerm, randomTerm = randomTerm
  ), silent = TRUE)
expect(!inherits(result1, "try-error"), "premetLMMsolver() ran")
if (inherits(result1, "try-error")) {
  cat("         -> ", as.character(result1), "\n", sep = "")
} else {
  cat("  gxeModelNum: ", result1$gxeModelNum,
      " | gxeTerms: ", paste(unlist(result1$gxeTerms), collapse = ", "), "\n", sep = "")
}

cat("\n### 2/3 metLMMsolver (main effects + GEBV) ###\n")
result <- try(
  cgiarPipeline::metLMMsolver(
    phenoDTfile   = dtMta,
    analysisId    = staId,
    analysisIdGeno = qaGenoId,      # supplied because genoA is in the covariates
    fixedTerm     = fixedTerm,
    randomTerm    = randomTerm,
    expCovariates = expCovariates,
    envsToInclude = myEnvsTI,
    trait         = trait,
    traitFamily   = myFamily,
    useWeights    = TRUE,
    estHybrids    = TRUE,
    calculateSE   = TRUE,
    heritLB       = 0.1,
    heritUB       = 0.95,
    meanLB        = 0,
    meanUB        = 1000000,
    nPC           = nPC,
    subsetGeno    = -1,
    subsetPed     = -1,
    maxIters      = 35,
    verbose       = FALSE
  ), silent = TRUE)

expect(!inherits(result, "try-error"), "metLMMsolver() ran on the generated object")
if (inherits(result, "try-error")) {
  cat("         -> ", as.character(result), "\n", sep = "")
} else {
  mtaId <- result$status$analysisId[nrow(result$status)]
  expect("mtaLmms" %in% result$status$module, "an 'mtaLmms' status row was written")

  preds <- result$predictions[result$predictions$module == "mtaLmms" &
                                result$predictions$analysisId == mtaId, ]
  expect(nrow(preds) > 0, "MTA predictions were produced")
  cat("  predictions: ", nrow(preds), " rows\n", sep = "")
  cat("  effectTypes: ", paste(sort(unique(preds$effectType)), collapse = ", "), "\n", sep = "")

  # The GEBV check: a genomic main-effect model must yield designation
  # predictions for the genotyped material.
  desPreds <- preds[preds$effectType == "designation" & !is.na(preds$predictedValue), ]
  expect(nrow(desPreds) > 0, "GEBVs (designation effectType) were estimated")
  cat("  designations with a GEBV: ", length(unique(desPreds$designation)), "\n", sep = "")
  cat("  predictedValue range    : ",
      paste(round(range(desPreds$predictedValue), 2), collapse = " .. "), "\n", sep = "")

  # Confirm the genomic relationship actually entered the model rather than the
  # run silently degrading to a plain designation effect (TGV). LMMsolver
  # encodes the relationship inside grp(designation), so the formula itself says
  # nothing; metLMMsolver records the kernel used in a `kernels` modeling row.
  mdl <- result$modeling[result$modeling$analysisId == mtaId, ]
  mtr <- result$metrics[result$metrics$analysisId == mtaId, ]
  cat("  fixedFormula : ", paste(unique(mdl$value[mdl$parameter == "fixedFormula"]), collapse = " | "), "\n", sep = "")
  cat("  randomFormula: ", paste(unique(mdl$value[mdl$parameter == "randomFormula"]), collapse = " | "), "\n", sep = "")

  kernels <- mdl$value[mdl$parameter == "kernels"]
  cat("  kernels      : ", paste(kernels, collapse = ", "), "\n", sep = "")
  expect(any(kernels == "genoA"),
         "modeling records kernels = genoA, i.e. GEBV was the surrogate of merit")

  # How many environments the mixed model actually kept after the H2 filter.
  inc <- mdl[mdl$parameter == "includedInMta", c("environment","value")]
  cat("  environments included:\n")
  print(inc, row.names = FALSE)
  nEnvUsed <- sum(inc$value == "TRUE")
  expect(nEnvUsed >= 2,
         paste0("at least 2 environments entered the model (got ", nEnvUsed, ")"))

  varParams <- mtr$parameter[grepl("^Var", mtr$parameter)]
  cat("  variance components: ", paste(unique(varParams), collapse = ", "), "\n", sep = "")
  r2 <- mtr[grepl("^r2", mtr$parameter), c("parameter","value")]
  if (nrow(r2) > 0) {
    cat("  reliability/r2: ", paste(r2$parameter, round(r2$value, 3),
                                    sep = "=", collapse = ", "), "\n", sep = "")
    expect(any(r2$value > 0),
           "r2 > 0, so the model fitted rather than falling back to means")
  }

  # ---- GEBV vs TGV contrast ------------------------------------------------
  # Re-run with no covariance on the designation effect. If the genomic
  # relationship is genuinely being used, the two sets of predictions must
  # differ; identical values would mean genoA was silently ignored.
  cat("\n  --- contrast: same model with no genomic relationship (TGV) ---\n")
  tgv <- try(
    cgiarPipeline::metLMMsolver(
      phenoDTfile = dtMta, analysisId = staId, analysisIdGeno = NULL,
      fixedTerm = fixedTerm, randomTerm = randomTerm,
      expCovariates = list("none"),
      envsToInclude = myEnvsTI, trait = trait, traitFamily = myFamily,
      useWeights = TRUE, estHybrids = TRUE, calculateSE = TRUE,
      heritLB = 0.1, heritUB = 0.95, meanLB = 0, meanUB = 1000000,
      nPC = nPC, subsetGeno = -1, subsetPed = -1, maxIters = 35, verbose = FALSE
    ), silent = TRUE)

  if (inherits(tgv, "try-error")) {
    cat("  (TGV contrast run failed, skipping comparison)\n")
  } else {
    tgvId <- max(tgv$status$analysisId[tgv$status$module == "mtaLmms"])
    tp <- tgv$predictions[tgv$predictions$analysisId == tgvId &
                            tgv$predictions$effectType == "designation",
                          c("designation","predictedValue")]
    gp <- desPreds[, c("designation","predictedValue")]
    cmp <- merge(gp, tp, by = "designation", suffixes = c("_gebv","_tgv"))
    cmp <- cmp[stats::complete.cases(cmp), ]
    cat("  designations compared: ", nrow(cmp), "\n", sep = "")
    if (nrow(cmp) > 2) {
      rho <- suppressWarnings(stats::cor(cmp$predictedValue_gebv, cmp$predictedValue_tgv))
      maxAbsDiff <- max(abs(cmp$predictedValue_gebv - cmp$predictedValue_tgv))
      cat("  correlation GEBV vs TGV : ", round(rho, 4), "\n", sep = "")
      cat("  max |GEBV - TGV|        : ", round(maxAbsDiff, 4), "\n", sep = "")
      expect(maxAbsDiff > 1e-6,
             "GEBVs differ from the no-relationship run, so genoA shaped the estimates")
    }
  }

  cat("\n### 3/3 postmetLMMsolver ###\n")
  result2 <- try(
    cgiarPipeline::postmetLMMsolver(
      phenoDTfile = result, analysisId = mtaId,
      gxeModelNum = result1$gxeModelNum, gxeTerms = result1$gxeTerms
    ), silent = TRUE)
  expect(!inherits(result2, "try-error"), "postmetLMMsolver() ran")
  if (inherits(result2, "try-error")) {
    cat("         -> ", as.character(result2), "\n", sep = "")
  } else {
    cat("  predictions after postmet: ", nrow(result2$predictions), " rows\n", sep = "")
    cat("  effectTypes: ",
        paste(sort(unique(result2$predictions$effectType[
          result2$predictions$analysisId == mtaId])), collapse = ", "), "\n", sep = "")
    save(result2, file = "test/mta_gebv_output.RData")
    cat("  saved test/mta_gebv_output.RData\n")
  }
}

cat("\n")
if (length(failures) > 0) {
  stop(sprintf("%d check(s) failed:\n  - %s", length(failures),
               paste(failures, collapse = "\n  - ")), call. = FALSE)
}
cat("Main effects model with GEBV ran successfully on the generated object.\n")
