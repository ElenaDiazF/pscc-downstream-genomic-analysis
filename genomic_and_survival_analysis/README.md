# PSCC Genomic and Survival Data Analysis Workflow

This folder contains R/R Markdown scripts used for the PSCC downstream genomic analysis after creation of the filtered somatic variant dataset. The workflow starts from filtered variant tables and metadata, performs initial sample checks, integrates TP53 genetic status with p53 IHC information, creates `maftools` MAF objects and uses those objects for downstream analyses and visualizations.

Most downstream scripts expect objects created by earlier notebooks to already exist in the same R environment. Therefore, the first three scripts must be run before running the remaining downstream analysis scripts.

## Required execution order

Once the REDCap metadata and the filtered variant list from the server-side bioinformatic pipeline are available, run these three scripts first, in this exact order:

1.  `1upload_and_initial_sample_checks.Rmd`\
    Data upload and initial sample checks.
2.  `2tp53_ihc_variants.Rmd`\
    TP53 sequencing status and p53 IHC integration.
3.  `3maf_obj_creation.Rmd`\
    MAF-compatible table generation and MAF object creation.

Use the same R session or load a saved workspace containing the objects created by these three scripts.

## Core workflow scripts

### `1upload_and_initial_sample_checks.Rmd`

**Purpose:** Loads the filtered variant table, REDCap clinical data, sample metadata, and QC tables. It performs initial checks related to HPV status, metadata consistency, coverage, contamination, excluded samples and passing samples without filtered variants.

**Main inputs:**

-   `filtered_variants.fst`,
-   REDCap connection credentials from `tokens.R`,
-   `metadata-pscc.csv`,
-   `low_coverage_samples_mutational_load.tsv`,
-   `sample_contamination_stats.tsv`.

**Main outputs/objects:**

-   `list_analyzed`,
-   `oncomine_followup`,
-   `patient_hpv`,
-   `metadata_pscc`,
-   `low_coverage`,
-   `contamination`,
-   `metadata_passing`,
-   `metadata_not_passing`,
-   `tumor_not_passing_samples`,
-   `premalignant_not_passing_samples`,
-   `passing_without_variants`,
-   `metadata_summary`,
-   `variant_summary`,
-   helper functions such as `save_flextable_docx()`, `normalize_name()`, and `find_column()`.

**Dependencies:** Does not depend on previous scripts.

### `2tp53_ihc_variants.Rmd`

**Purpose:** Compares tumor TP53 variant status from sequencing with p53 IHC status from clinical follow-up data. It classifies tumor samples as TP53-mutated or TP53-wild-type and creates combined TP53/p53/HPV summary tables.

**Main inputs:** Objects from `1upload_and_initial_sample_checks.Rmd`, especially:

-   `list_analyzed`,
-   `oncomine_followup`,
-   `metadata_passing`,
-   `tissue_col`,
-   `tumor_not_passing_samples`,
-   `patient_hpv`,
-   `save_flextable_docx()`.

**Main outputs/objects:**

-   `list_analyzed_tumor`,
-   `ihc_table`,
-   `tp53mut`,
-   `samples_tp53mut`,
-   `samples_tp53wt`,
-   `metadata_passing_tumor_id`,
-   `ihc_genetics`,
-   Fisher test contingency matrix and results,
-   `ihc_genetics_complete`.

**Dependencies:** Run after `1upload_and_initial_sample_checks.Rmd`.

### `3maf_obj_creation.Rmd`

**Purpose:** Converts the filtered variant table into MAF-compatible tables, prepares clinical annotations, recodes clinical variables, and creates MAF objects for tumor, premalignant, HPV-associated, HPV-independent and TP53-status subgroups.

**Main inputs:** Objects from the first two scripts, especially:

-   `list_analyzed`,
-   `metadata_pscc`,
-   `metadata_not_passing`,
-   `oncomine_followup`,
-   `ihc_genetics`.

**Main outputs/objects:**

-   `maf_table_all`,
-   `filtered_variants`,
-   `metadata_filtered`,
-   `clinical_for_maf_all`,
-   `clinical_for_maf`,
-   `maf_table_tumor`,
-   `maf_tumor_premalignant`,
-   `maf_tumor`,
-   `maf_hpva`,
-   `maf_hpvi`,
-   `maf_hpvi_tp53wt`,
-   `maf_hpvi_tp53mut`,
-   `maf_objects`,
-   `maf_build_summary`,
-   `samples_hpva`,
-   `samples_hpvi`,
-   `samples_hpvi_tp53wt`,
-   `samples_hpvi_tp53mut`,
-   `ordered_samples`,
-   `group_order`,
-   `clinical_three_groups`,
-   `maf_three_groups`,
-   `top_genes_global`,
-   `oncoplot_ann_colors`,
-   helper functions including `prepare_maf_table()`, `add_dummy_samples()`, `top_genes()`, `build_maf_from_clinical()`, and `get_ordered_genes()`.

**Dependencies:** Run after `1upload_and_initial_sample_checks.Rmd` and `2tp53_ihc_variants.Rmd`.

## Downstream analysis scripts

The following scripts should be run only after the three core scripts above have been executed in order.

### `oncoplot_Generation.Rmd`

**Purpose:** Generates MAF summaries, tumor/premalignant and tumor-only oncoplots, HPV subgroup oncoplots, MAF comparisons, grouped oncoplots, and exports the tumor MAF object and source tables.

**Main inputs:** `maf_tumor_premalignant`, `maf_tumor`, `maf_hpva`, `maf_hpvi`, `maf_hpvi_tp53wt`, `maf_hpvi_tp53mut`, `maf_three_groups`, `ordered_samples`, `oncoplot_ann_colors`, `top_genes()`, `get_ordered_genes()` and `save_flextable_docx()`.

**Main outputs/objects:** MAF summary objects, oncoplot objects, HPV comparison results, grouped oncoplots, and files exported under `maf_tumor_export/`.

**Dependencies:** Requires `3maf_obj_creation.Rmd`.

### `metrics_maf.Rmd`

**Purpose:** Calculates tumor mutation burden, mutation count summaries, recurrent variant metrics, variant classification/type counts, SNV substitution classes and top mutated genes.

**Main inputs:** `maf_tumor`, `clinical_for_maf`, `list_analyzed`, and `save_flextable_docx()`.

**Main outputs/objects:** `tmb_table`, `median_tmb`, `tmb_with_clinical`, `median_per_group`, `test_global`, `test_results_hpv`, `sample_summary`, `variant_classification_counts`, `variant_type_counts`, `snv_class_counts`and `top_genes_table`.

**Dependencies:** Requires `3maf_obj_creation.Rmd`.

### `genes_into_pathways.Rmd`

**Purpose:** Assigns recurrently mutated tumor genes to functional pathways, creates pathway oncoplots, calculates pathway-level mutation counts and performs global and pairwise pathway comparisons across HPV/TP53 groups.

**Main inputs:** `maf_tumor`, `clinical_for_maf`, `ordered_samples`, `oncoplot_ann_colors`, `group_order`, and the pathway annotation file `gene_annotations_recurrence_3.xlsx`.

**Main outputs/objects:** `gene_summary`, `recurrent_genes`, `Robert_pathways`, pathway oncoplot objects, `sigpw_pathways`, `sample_groups`, `pathway_counts`, `stats_table` pairwise pathway plot lists, and pairwise statistics tables.

**Dependencies:** Requires `3maf_obj_creation.Rmd`.

### `lollipop.Rmd`

**Purpose:** Creates protein-domain lollipop plots for recurrently mutated genes in the tumor MAF object and compares HPV-independent versus HPV-associated tumors.

**Main inputs:** `maf_tumor`, `maf_hpvi`, and `maf_hpva`. These MAF objects must contain an amino-acid change column such as `hgvs_p`, `HGVSp_Short`, `Protein_Change`, or `AAChange`.

**Main outputs/objects:** `aa_column`, `top15_tumor_genes`, `genes_tumor`, `top15_hpv_genes`, `genes_hpv` and exported PNG lollipop plots.

**Dependencies:** Requires `3maf_obj_creation.Rmd`.

### `boxplot_cyclind1.Rmd`

**Purpose:** Summarizes and plots Cyclin D1 IHC percentages by HPV/TP53 mutational group and runs global and pairwise non-parametric tests.

**Main inputs:** `clinical_for_maf`, specifically `HPV_TP53_mutational_group` and `cyclind1_ihc_perc`.

**Main outputs/objects:** `ccd1_ihc_by_hpvtp53groups`, `ccdn1_kruskal`, and `ccdn1_wilcox`.

**Dependencies:** Requires `3maf_obj_creation.Rmd`.

### `cohort_table_generation.Rmd`

**Purpose:** Builds a grouped clinical cohort table using REDCap follow-up variables and HPV/TP53 group annotations, then exports it as a Word table.

**Main inputs:** `oncomine_followup`, `clinical_for_maf`, and `save_flextable_docx()`.

**Main outputs/objects:** `cohort_table_all`, `cohort`, a Word export named `cohort.docx`, and a saved R workspace at the path specified in the script.

**Dependencies:** Requires `1upload_and_initial_sample_checks.Rmd` for REDCap data and `3maf_obj_creation.Rmd` for `clinical_for_maf`.

### `detect_outlier_sample.R`

**Purpose:** Calculates TMB, applies Tukey's upper outlier cutoff (`Q3 + 2.5 * IQR`), identifies hypermutated samples, and plots TMB by sample.

**Main inputs:** `maf_tumor`, `clinical_for_maf_all`, and TMB-related objects from the MAF workflow. The plotting section also refers to `tmb_with_clinical`.

**Main outputs/objects:** `tmb_table_all`, `tmb_with_clinical_all`, `tmb_values`, `Q1`, `Q3`, `IQR_value`, `upper_cutoff`, and `hypermutated_samples`.

**Dependencies:** Requires `3maf_obj_creation.Rmd`; the plot may also require `metrics_maf.Rmd` if `tmb_with_clinical` is not already present.

### `survival_recurrence.Rmd`

**Purpose:** Performs recurrence and disease-specific survival analyses using TP53/sample metadata, REDCap clinical data, Kaplan-Meier plots and Cox regression models.

**Main inputs:** A TP53/sample metadata CSV file, REDCap clinical data from `tokens.R`, and survival-related R packages.

**Main outputs/objects:** `tp53`, `clinical`, `metadata_pscc_filtered`, `metadata_pscc_clinical`, `metadata_pscc_clinical_60m`, `group_genomic`, `group_ihc`, Kaplan-Meier fit and plot objects, `metadata_prepared`, `cox_uni_tbl_dds`, `cox_tp53_stage_dds` and `cox_uni_tbl_rec`.

**Dependencies:** Partly independent. It loads its own TP53/sample CSV and clinical data, but uses concepts and variables produced by the filtering/MAF workflow. Run the core workflow first when reproducing the full analysis.

## Practical reproduction notes

-   Run the first three scripts in the same R session, in the required order.
-   After `3maf_obj_creation.Rmd`, downstream scripts can be run as needed depending on the figures, tables, or analyses being reproduced.
-   Several scripts use absolute paths under a local project directory. These paths should be adapted when reproducing the workflow in another environment.
-   Some scripts write outputs, including Word tables, PNG plots, exported MAF files, and saved `.RData` workspaces.
-   This folder does not include raw data, patient-level clinical files, REDCap credentials, or large intermediate objects.
