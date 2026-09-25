# Generates a small, internally consistent EBS-style dataset for verifying the
# three Bioflow interoperability use cases.
#
# The shipped fixtures cannot exercise them: PedF1.csv belongs to a different
# dataset than bioflow_pheno_data.csv (zero shared identifiers), the VCF holds
# only 4 markers (below crossVerification's 8-marker floor), and the phenotype
# file has a single occurrence (MTA needs >= 2 environments).
#
# What this builds:
#   * 7 homozygous parents (3 mothers x 4 fathers)
#   * 12 F1 designations, each genotyped under its own name so Pedigree qa/qc
#     can form designation/mother/father triplets
#   * 2 extra samples for one F1 designation, so the sample_id path used by
#     F1 qa/qc is exercised too
#   * 3 deliberately mislabelled F1s, so the QA metrics have two clusters and
#     mclust::Mclust(G = 2) has something real to separate
#   * 30 biallelic diploid SNPs
#   * phenotypes over 2 occurrences (= 2 Bioflow environments) with 2 reps

make_coherent_fixture <- function(outDir = "test", seed = 20260924) {
  set.seed(seed)

  mothers <- paste0("MP", 1:5)
  fathers <- paste0("FP", 1:5)
  parents <- c(mothers, fathers)

  crosses <- expand.grid(mother = mothers, father = fathers, stringsAsFactors = FALSE)
  crosses$designation <- paste0("F1_", sprintf("%02d", seq_len(nrow(crosses))))

  # Comfortably more markers than individuals, so the genomic relationship
  # matrix is well conditioned for the GEBV model.
  nMarkers <- 120
  markerIds <- paste0("M", sprintf("%03d", seq_len(nMarkers)))

  # ---- reference/alternate alleles -----------------------------------------
  bases <- c("A", "C", "G", "T")
  refAl <- character(nMarkers)
  altAl <- character(nMarkers)
  for (i in seq_len(nMarkers)) {
    pair <- sample(bases, 2)
    refAl[i] <- pair[1]
    altAl[i] <- pair[2]
  }

  # ---- parent genotypes: homozygous, so F1 dosages are deterministic -------
  # Homozygous parents also keep parental heterozygosity at 0%, which clears the
  # default parentHetThreshold in the F1 module.
  parentDose <- matrix(
    sample(c(0L, 2L), length(parents) * nMarkers, replace = TRUE),
    nrow = length(parents), dimnames = list(parents, markerIds)
  )

  # Guarantee every marker is polymorphic across the parent panel, otherwise it
  # carries no information and inflates the monomorphic count.
  for (j in seq_len(nMarkers)) {
    if (length(unique(parentDose[, j])) == 1L) {
      parentDose[sample(seq_along(parents), 1), j] <- 2L - parentDose[1, j]
    }
  }

  # ---- F1 genotypes from the TRUE parents ---------------------------------
  f1Dose <- matrix(0L, nrow = nrow(crosses), ncol = nMarkers,
                   dimnames = list(crosses$designation, markerIds))
  for (i in seq_len(nrow(crosses))) {
    f1Dose[i, ] <- as.integer(
      (parentDose[crosses$mother[i], ] + parentDose[crosses$father[i], ]) / 2L
    )
  }

  # ---- extra samples for one designation (multi-sample F1 path) ------------
  dupDesignation <- crosses$designation[1]
  dupSamples <- paste0(dupDesignation, c("_s2", "_s3"))
  dupDose <- rbind(f1Dose[dupDesignation, ], f1Dose[dupDesignation, ])
  rownames(dupDose) <- dupSamples

  genoDose <- rbind(parentDose, f1Dose, dupDose)

  # ---- pedigree ------------------------------------------------------------
  # Parents: no parentage of their own. sample_id equals the genotyped name.
  # `plant_no` maps onto Bioflow's `other` parameter; `seed_source` and
  # `nursery_code` have no canonical counterpart and must survive verbatim.
  pedParents <- data.frame(
    designation  = parents,
    sample_id    = parents,
    mother       = NA_character_,
    father       = NA_character_,
    entry        = "parent",
    year         = 2024L,
    plant_no     = NA_integer_,
    seed_source  = "genebank",
    nursery_code = paste0("N", sprintf("%03d", seq_along(parents))),
    stringsAsFactors = FALSE
  )

  pedF1 <- data.frame(
    designation  = crosses$designation,
    sample_id    = crosses$designation,
    mother       = crosses$mother,
    father       = crosses$father,
    entry        = "F1",
    year         = 2025L,
    plant_no     = seq_len(nrow(crosses)),
    seed_source  = "crossing_block",
    nursery_code = paste0("N", sprintf("%03d", 100 + seq_len(nrow(crosses)))),
    stringsAsFactors = FALSE
  )

  pedDup <- data.frame(
    designation  = dupDesignation,
    sample_id    = dupSamples,
    mother       = crosses$mother[1],
    father       = crosses$father[1],
    entry        = "F1",
    year         = 2025L,
    plant_no     = c(101L, 102L),
    seed_source  = "crossing_block",
    nursery_code = c("N201", "N202"),
    stringsAsFactors = FALSE
  )

  pedigree <- rbind(pedParents, pedF1, pedDup)

  # Deliberately mislabel four F1s by swapping in the wrong mother. Their
  # markers still come from the true mother, so QA should flag them.
  mislabelled <- crosses$designation[c(2, 7, 13, 21)]
  for (d in mislabelled) {
    trueMother <- crosses$mother[crosses$designation == d]
    wrongMother <- setdiff(mothers, trueMother)[1]
    pedigree$mother[pedigree$designation == d & pedigree$sample_id == d] <- wrongMother
  }

  pedigreePath <- file.path(outDir, "coherent_pedigree.csv")
  utils::write.csv(pedigree, pedigreePath, row.names = FALSE, na = "")

  # ---- VCF -----------------------------------------------------------------
  samples <- rownames(genoDose)
  gtOf <- function(d) c("0/0", "0/1", "1/1")[d + 1L]

  vcfLines <- c(
    "##fileformat=VCFv4.2",
    "##source=make_coherent_fixture.R",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
            samples), collapse = "\t")
  )
  for (j in seq_len(nMarkers)) {
    vcfLines <- c(vcfLines, paste(c(
      "1", j * 1000L, markerIds[j], refAl[j], altAl[j], ".", "PASS", ".", "GT",
      gtOf(genoDose[, j])
    ), collapse = "\t"))
  }
  vcfPath <- file.path(outDir, "coherent_geno.vcf")
  writeLines(vcfLines, vcfPath)

  # ---- phenotypes ----------------------------------------------------------
  # Designations are phenotyped (parents and F1s); genotypes live on samples.
  designations <- c(parents, crosses$designation)
  occurrences  <- c("occ_A", "occ_B")
  reps         <- 1:3

  pheno <- expand.grid(
    germplasmName  = designations,
    rep            = reps,
    occurrenceName = occurrences,
    stringsAsFactors = FALSE
  )
  n <- nrow(pheno)

  pheno$germplasmDbId <- match(pheno$germplasmName, designations)
  pheno$breedingStage <- "S1"
  pheno$year          <- 2025L
  pheno$season        <- "main"
  pheno$site          <- ifelse(pheno$occurrenceName == "occ_A", "Nairobi", "Kiboko")
  pheno$experimentName <- "exp_1"
  pheno$design        <- "RCBD"
  pheno$blockNumber   <- pheno$rep
  pheno$entryType     <- ifelse(pheno$germplasmName %in% parents, "parent", "test")
  pheno$plotNumber    <- seq_len(n)

  # a simple 2-D layout per occurrence so the spatial term has coordinates
  perOcc <- length(designations) * length(reps)
  side   <- ceiling(sqrt(perOcc))
  idxInOcc <- ave(seq_len(n), pheno$occurrenceName, FUN = seq_along)
  pheno$paX <- ((idxInOcc - 1) %% side) + 1L
  pheno$paY <- ((idxInOcc - 1) %/% side) + 1L

  # Trait with a genuinely polygenic signal: a weighted sum of marker dosages,
  # standardised, so genetic variance is large relative to the residual and
  # heritability is estimable in BOTH environments. A marker-driven trait is
  # also what makes a GEBV model meaningful - a trait unrelated to the markers
  # would give a genomic relationship nothing to explain.
  # The residual sd is chosen so plot-level heritability lands near 0.7, i.e.
  # comfortably inside MTA's default heritLB/heritUB of 0.1 / 0.95. Too clean a
  # trait pushes H2 above 0.95 and metLMMsolver then excludes the environment.
  markerEffects <- stats::rnorm(nMarkers)
  trueBV <- as.vector(scale(genoDose %*% markerEffects))
  names(trueBV) <- rownames(genoDose)

  bv <- trueBV[match(pheno$germplasmName, names(trueBV))]
  bv[is.na(bv)] <- 0
  envEffect <- ifelse(pheno$occurrenceName == "occ_A", 0, 12)
  pheno$Plant_Height_cm <- round(
    120 + 8 * bv + envEffect + stats::rnorm(n, 0, 5), 2
  )

  phenoPath <- file.path(outDir, "coherent_pheno.csv")
  utils::write.csv(pheno, phenoPath, row.names = FALSE)

  invisible(list(
    phenotypeFile = phenoPath,
    pedigreeFile  = pedigreePath,
    genotypeFile  = vcfPath,
    trait         = "Plant_Height_cm",
    parents       = parents,
    f1            = crosses$designation,
    mislabelled   = mislabelled,
    nMarkers      = nMarkers
  ))
}
