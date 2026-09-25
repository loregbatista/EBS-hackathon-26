# EBS-hackathon-26

Here we will store the code produced in the EBS-hackathon-26.

## Bioflow input conversion script

`scripts/bundled_getBioflowRdata.R` converts EBS phenotype, pedigree and
genotype (VCF) exports into the `result` object that Bioflow loads, and writes
it as an `.RData` file.

### Supported Bioflow use cases

| Use case | Needs phenotypes | Needs pedigree | Needs genotypes |
| --- | --- | --- | --- |
| 1. Multi Trial Analysis (MTA) | yes | only for SCA/GCA models | only for genomic models |
| 2. F1 qa/qc | no | yes | yes, progeny **and** parents |
| 3. Pedigree qa/qc | no | yes | yes, designation/mother/father triplets |

Single Trial Analysis is run during conversion (`cgiarPipeline::staLMM`), so the
object arrives MTA-ready. Pass `runSTA = FALSE` to skip it and export a
pedigree/genotype-only object for use cases 2 and 3.

After writing the file the script prints a **readiness report** stating, per use
case, which requirements are met and what is missing. Each check names the
Bioflow module or pipeline function that enforces it, so a failure is
actionable.

### What it does

- Checks for required CRAN packages and installs missing ones (`vcfR`,
  `adegenet`, `cli`, `rlang`, `remotes`, `openssl`, `stringr`, `purrr`, `glue`,
  `Matrix`), plus `cgiarBase` / `cgiarPipeline` from GitHub
- Reads the phenotype csv, derives the Bioflow `environment` column and applies
  the EBS design rules for `rep`
- Reads the pedigree csv: maps it onto Bioflow's pedigree vocabulary and carries
  every remaining column through unchanged
- Reads genotype data from VCF, compressed (`.gz`) or plain text, auto-detecting
  the ploidy level; filters and imputes markers, logging both
- Runs Single Trial Analysis and writes the `result` object as `.RData`

### Usage

```r
# NOTE: source(), not library() - this is a script, not a package
source("https://raw.githubusercontent.com/Breeding-Analytics/EBS-hackathon-26/refs/heads/main/scripts/bundled_getBioflowRdata.R")

traits <- c("Pollen_DAP_days","Silk_DAP_days","Plant_Height_cm","Ear_Height_cm",
            "Root_Lodging_plants","Stalk_Lodging_plants",
            "Yield_Mg_ha","Grain_Moisture","Twt_kg_m3")

out <- getBioflowRData(
  phenotypeFile = "test/pheno.csv",
  pedigreeFile  = "test/PedF1.csv",
  genotypeFile  = "test/inputF1.vcf",
  traits        = traits,
  outputPath    = "test",                 # directory
  outputFile    = "bioflow_input_test"    # file name, without the .RData suffix
)

out$file     # path of the file that was written
out$result   # the object itself
```

`outputPath` is the directory and `outputFile` is the file name. When
`outputFile` is `NULL` the file is named with the md5 hash of the phenotype
`analysisId`, which is the EBS content-addressed convention.

### Pedigree file

At minimum the pedigree csv needs **designation, mother and father**. Unknown
parents should be empty or `NA`.

Column names are auto-detected, case- and separator-insensitively, onto
Bioflow's pedigree vocabulary:

| Bioflow parameter | Recognised spellings (examples) | Required for |
| --- | --- | --- |
| `designation` | designation, germplasmName, line, genotype | all pedigree use cases |
| `mother` | mother, female, parent1, hybrid_parent1 | all pedigree use cases |
| `father` | father, male, parent2, hybrid_parent2 | all pedigree use cases |
| `crossType` | crossType, crossClass, materialType | F1 qa/qc |
| `sample_id` | sample_id, sampleName, dnaSampleId | F1 qa/qc with several samples per designation |
| `yearOfOrigin` | yearOfOrigin, year | realized genetic gain |
| `other` | other, plant_no, plantNumber | F1 qa/qc dashboard |
| `batch`, `plate`, `position` | batch, plate, well | Pedigree qa/qc plate layout |

**Every other column is carried into the R object unchanged**, so nothing in the
pedigree export is lost.

Override auto-detection when a header is ambiguous:

```r
out <- getBioflowRData(
  ...,
  pedigreeMapping = list(crossType = "entry", yearOfOrigin = "year")
)
```

Notes on `crossType`, which drives the F1 qa/qc module:

- Bioflow compares `crossType == "F1"`, exactly and case-sensitively. Recognised
  spellings (`f1`, `F1`, `hybrid`, `cross`, ...) are normalised to `F1`.
- Any other label marks a non-progeny row; `parent` is the usual choice.
- Missing values are filled with `parent`, because Bioflow's `==` comparison
  turns an `NA` into a phantom row that later aborts the run.
- If the file has no crossType column, one is derived from parent completeness
  (both parents recorded means F1). Pass `deriveCrossType = FALSE` to opt out.

### Identifier matching

This is the most common reason a file loads but a module then refuses to run.
The pedigree identifiers must match the **VCF sample names**:

- Pedigree qa/qc compares `designation`, `mother` and `father` against the
  marker matrix, so all three must be genotyped for a row to be evaluated.
- F1 qa/qc uses `sample_id` for the progeny when that column exists, and
  `designation` otherwise; `mother` and `father` must be genotyped either way.

The readiness report reports these overlaps as counts, so a mismatch shows up at
export time rather than as an empty dashboard.

### Column mapping and bundling

The scripts are maintained as modular files under `scripts/` and bundled into a
single sourceable file so the GitHub raw URL above works without installing
anything:

```bash
Rscript scripts/bundle_scripts.R
```

Bundle order is `packages_verification.R`, `read_geno_functions.R`,
`read_pedigree_functions.R`, `validate_bioflow_object.R`, `getBioflowRdata.R`.
Do not edit `scripts/bundled_getBioflowRdata.R` by hand; it is generated.

### Tests

```bash
# ingestion against the shipped sample files
Rscript test/test_create_bioflow_input.R

# end-to-end check of all three use cases on a coherent generated fixture,
# including a real cgiarPipeline::individualVerification() run
Rscript test/test_use_cases.R
```

The shipped sample files (`bioflow_pheno_data.csv`, `PedF1.csv`,
`bioflow_genotype_data_fix.vcf`) exercise the conversion but cannot drive the
use cases end to end: they come from different datasets and share no
identifiers, the VCF has 4 markers (below the 8-marker floor in
`cgiarBase::crossVerification`), and the phenotype file has a single occurrence
(MTA needs at least two environments). `test/make_coherent_fixture.R` generates
a consistent dataset for that purpose.
