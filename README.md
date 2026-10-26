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

When phenotypic data is supplied, the script runs **model-based outlier
detection** (`cgiarPipeline::resMod`, recorded as a `qaMb` stamp) followed by
**Single Occurrence / Trial Analysis** (`cgiarPipeline::staLMM`, recorded as a
`soa` stamp), so the object arrives MTA/MOA-ready. This mirrors Bioflow's qaSta
"Model-Based Outlier Detection" + SOA step: the EBS pheno QA/QC steps
(experimental design, missingness, consistency) are performed upstream in EBS,
and the only analytics replicated here is the model-based outlier detection that
EBS does not have.

`runSOA` controls this step. The default (`NULL`) is **auto**: SOA runs whenever
a phenotype file is present. Pass `runSOA = FALSE` to force a
pedigree/genotype-only object for use cases 2 and 3, or `runSOA = TRUE` to force
it on.

> **Requires an up-to-date `cgiarPipeline`.** The current `staLMM` returns
> `NULL` unless a `qaMb` stamp exists and records its output under module
> `soa` (previously `sta`). Install `cgiarPipeline` from a source new enough to
> export `resMod`; the script's dependency bootstrap handles this when it
> installs from GitHub.

After writing the file the script prints a **readiness report** stating, per use
case, which requirements are met and what is missing. Each check names the
Bioflow module or pipeline function that enforces it, so a failure is
actionable.

It also writes a **diagnostics JSON** next to the `.RData` (named
`<outputFile>.diagnostics.json`, or pass `diagnosticsFile` to choose the path).
The JSON is written for the EBS software engineer preparing the input files: for
each use case it gives a `READY` / `NOT READY` status and, for every unmet
requirement, a plain-language `what_it_means` and a concrete `how_to_fix`
action, alongside the original `technical_detail`. A top-level `overall.verdict`
summarises which use cases the file can drive. The path is returned as
`out$diagnosticsFile`.

### Production (Bioflow `main`) compatibility

The EBS team loads these objects in the production Bioflow instance and runs
**MOA** there. Because the model-based QA + SOA are performed here (outside
Bioflow) and recorded as `qaMb` + `soa` stamps, the object is loaded ready for
MOA: Bioflow `main`'s MOA module selects SOA stamps by the `soa` module, so a
script-generated object is both selectable and runnable for MOA without running
Single Trial Analysis inside the production UI.

### What it does

- Checks for required CRAN packages and installs missing ones (`vcfR`,
  `adegenet`, `cli`, `rlang`, `remotes`, `openssl`, `stringr`, `purrr`, `glue`,
  `Matrix`, `jsonlite`), plus `cgiarBase` / `cgiarPipeline` from GitHub
- Reads the phenotype csv, derives the Bioflow `environment` column and applies
  the EBS design rules for `rep`
- Reads the pedigree csv: maps it onto Bioflow's pedigree vocabulary and carries
  every remaining column through unchanged
- Reads genotype data from VCF, compressed (`.gz`) or plain text, auto-detecting
  the ploidy level; filters and imputes markers, logging both
- Aligns genotype individual names to the pedigree identifier space: when a
  `sample_id` column is mapped, parental samples keyed by `sample_id` in the VCF
  are translated to their `designation` (F1 progeny keep their `sample_id`), so
  the qa/qc modules can match parents. This is a no-op when the VCF is already
  designation-keyed.
- When phenotypes are present, runs model-based outlier detection
  (`cgiarPipeline::resMod` -> `qaMb` stamp) and Single Occurrence/Trial Analysis
  (`cgiarPipeline::staLMM` -> `soa` stamp)
- Writes the `result` object as `.RData` and an engineer-friendly diagnostics
  JSON alongside it

### Usage to run Single Occurence Analysis

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
  outputFile    = "bioflow_input_test",    # file name, without the .RData suffix if empty md5 generated
  runSOA        = TRUE                      # NULL (default) also runs SOA when phenotypes are present
)

out$file            # path of the file that was written
out$result          # the object itself
out$diagnosticsFile # path of the diagnostics JSON
```

### Usage to run F1 QA/QC or pedigree verification

```r
# NOTE: source(), not library() - this is a script, not a package
source("https://raw.githubusercontent.com/Breeding-Analytics/EBS-hackathon-26/refs/heads/main/scripts/bundled_getBioflowRdata.R")

out <- getBioflowRData(
  phenotypeFile = NULL,
  pedigreeFile  = "test/PedF1.csv",
  genotypeFile  = "test/inputF1.vcf",
  traits        = NULL,
  outputPath    = "test",                 # directory
  outputFile    = "bioflow_input_test",    # file name, without the .RData suffix if empty md5 generated
  runSOA        = FALSE                     # force pedigree/genotype-only (skip model-based QA + SOA)
)

out$file            # path of the file that was written
out$result          # the object itself
out$diagnosticsFile # path of the diagnostics JSON
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

The crossType value also steers how the genotype names are aligned, which in
turn decides which qa/qc module the object can drive:

- Rows labelled `F1` are treated as F1 progeny: their genotypes keep the
  `sample_id` key (F1 qa/qc looks progeny up by `sample_id`), while parents are
  translated to `designation`. Use this for the **F1 qa/qc** use case.
- Rows with any **non-F1** label (e.g. `F6` or a later-generation code) are not
  treated as progeny, so their genotypes are translated to `designation` too.
  This lets their `designation`/`mother`/`father` form genotyped triplets for
  the **Pedigree qa/qc** use case. The same input therefore targets F1 qa/qc or
  Pedigree qa/qc depending only on the crossType label.

### Identifier matching

This is the most common reason a file loads but a module then refuses to run.
Bioflow's qa/qc modules match pedigree identifiers against the **VCF sample
names**:

- Pedigree qa/qc compares `designation`, `mother` and `father` against the
  marker matrix, so all three must be genotyped (as designations) for a row to
  be evaluated.
- F1 qa/qc uses `sample_id` for the progeny when that column exists, and
  `designation` otherwise; `mother` and `father` must be genotyped either way.

EBS genotype exports usually key the VCF by `sample_id` (sample code), while the
pedigree references parents by `designation`. The script bridges this gap: when
a `sample_id -> designation` mapping is present it translates parental genotype
names to `designation` so the parents match, keeping F1 progeny under their
`sample_id`. When the VCF is already keyed by `designation`, no translation
happens. A consequence is that the F1 **designations** themselves must be in the
VCF for Pedigree qa/qc to form triplets (cross-string designations that are not
genotyped will not).

The readiness report (and the diagnostics JSON) report these overlaps as counts,
so a mismatch shows up at export time rather than as an empty dashboard.

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

### Use-case examples

`use_cases/` holds one worked example per use case, each with its input files
and the generated `.RData`, diagnostics JSON and a `*_checks_log.txt` capturing
the readiness report:

- `use_cases/1_MOA/` - a coherent phenotype + pedigree + genotype triplet;
  runs model-based QA + SOA and is ready for MTA/MOA (and, being fully coherent,
  the two qa/qc use cases as well).
- `use_cases/2_F1_qaqc/` - a real EBS F1 export (sample_id-keyed VCF); ready for
  F1 qa/qc.
- `use_cases/3_Pedigree_qaqc/` - the same F1 export with `crossType` relabelled
  from `F1` to a later-generation code, so it is ready for Pedigree qa/qc
  instead.

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
