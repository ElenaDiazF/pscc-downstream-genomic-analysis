# PSCC Downstream Genomic Analysis Repository

This repository contains selected scripts developed during the master thesis project for the analysis of penile squamous cell carcinoma (PSCC) whole-exome sequencing (WES) data.

The repository is a curated and documented version of the work performed during the project. It is intended to describe the analysis workflow and facilitate reproducibility without including raw sequencing data, patient-level clinical data, intermediate files, large outputs, credentials or private server-specific files.

## Repository organization

The repository is organized into two main sections:

``` text
pscc-downstream-genomic-analysis/
|
|-- README.md
|
|-- genomic_and_survival_analysis/
|   |-- README.md
|   |-- 1upload_and_initial_sample_checks.Rmd
|   |-- 2tp53_ihc_variants.Rmd
|   |-- 3maf_obj_creation.Rmd
|   |-- oncoplot_Generation.Rmd
|   |-- metrics_maf.Rmd
|   |-- genes_into_pathways.Rmd
|   |-- lollipop.Rmd
|   |-- boxplot_cyclind1.Rmd
|   |-- cohort_table_generation.Rmd
|   |-- detect_outlier_sample.R
|   |-- survival_recurrence.Rmd
|
|-- filtering_strategy_refinement/
|   |-- README.md
|   |-- filter_variants.R
|   |-- reports_list12.R
|   |-- submit_Reports_list12.sh 
|   |-- Changes_Server.Rmd
```

## 1. Genomic and survival analysis

This section contains the R/R Markdown workflow used after obtaining the final filtered somatic variant dataset. The workflow starts from filtered variant tables and metadata, performs initial sample checks, integrates TP53 mutational status with p53 IHC information, creates `maftools` MAF objects, and uses those objects for downstream analyses and visualizations.

Main analyses covered in this section include:

-   sample retention and metadata checks,
-   TP53 sequencing status and p53 IHC integration,
-   MAF object creation,
-   mutational landscape analysis,
-   HPV and TP53 subgroup comparisons,
-   pathway-level exploration,
-   Cyclin D1 IHC analysis,
-   cohort table generation,
-   TMB and outlier sample evaluation,
-   disease-specific survival and recurrence-free survival analysis.

See:

``` text
genomic_and_survival_analysis/README.md
```

## 2. Filtering strategy refinement

This section documents selected scripts and server-side changes related to refinement of the somatic variant filtering strategy. These scripts were used to improve filtering flexibility, evaluate variant list behaviour, export QC metrics, and document key decisions applied to FFPE-derived WES samples.

Main topics covered in this section include:

-   refinement of filtering parameters,
-   caller-support logic,
-   impact-specific filtering rules,
-   SNV and INDEL filtering decisions,
-   COSMIC/ClinVar-related rescue or prioritization options,
-   coverage and contamination export,
-   server-side pipeline/QC/reporting changes,
-   reproducibility improvements.

See:

``` text
filtering_strategy_refinement/README.md
```

## Data availability and privacy

This repository does not include raw sequencing files, BAM/VCF files, patient-level clinical information, REDCap credentials, tokens, private project outputs or large intermediate files.

The scripts may contain references to expected input objects or local paths used during development. These paths should be adapted when reproducing the analysis in a different environment.

## Suggested citation in the thesis

The code developed for the PSCC downstream genomic analysis and filtering strategy refinement workflow is available in this repository. The repository contains selected, documented, and sanitized scripts generated during the master thesis project.
