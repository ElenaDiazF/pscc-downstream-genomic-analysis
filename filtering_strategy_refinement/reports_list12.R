cat("\n================================================\n")
cat("  PSCC WES Analysis Report Generation\n")
cat("================================================\n\n")

setwd("/PROJECTES/SQUAMOLAB/")
cat("[1/10] Loading libraries...\n")
library(data.table)

# Configure data.table threading
# With 200G memory and memory-optimized AF calculation (no large matrix),
# thread-local buffers (~16MB/thread for radix sort) are negligible
n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
max_threads <- max(n_threads - 1, 1)
setDTthreads(max_threads)
cat(sprintf(
    "✓ data.table configured to use %d threads\n",
    getDTthreads()
))

cat("[2/10] Loading squamr functions...\n")
source("code/squamr/R/utils.R")
source("code/squamr/R/contamination.R")
source("code/squamr/R/coverage.R")
source("code/squamr/R/merge_variants.R")
source("code/squamr/R/filter_variants.R")
source("code/squamr/R/ffperase.R")
source("code/squamr/R/misc.R")
source("code/squamr/R/graphics_utils.R")
source("code/squamr/R/plotting.R")
source("code/squamr/R/read_metadata.R")
cat("✓ All squamr functions loaded successfully\n\n")

sample_col <- "Sample_ID"
metadata_paths <- c("projects/pscc/wes/config/metadata-pscc.csv")
dirs_coverage <- c("projects/pscc/wes/analysis/postprocessing/mosdepth")
dirs_contamination <- c("projects/pscc/wes/analysis/calling/variant_calling")
dirs_variants <- c(
    "projects/pscc/wes/data/processed/vcf/FFPE/tumor_only",
    "projects/pscc/wes/data/processed/vcf/FFPE/matched",
    "projects/pscc/wes/data/processed/vcf/fresh/tumor_only",
    "projects/pscc/wes/data/processed/vcf/fresh/matched"
)
dirs_ffperase <- c("projects/pscc/wes/analysis/ffpe_denoise")

# ============================================================================
# Load metadata once for efficiency
# ============================================================================
cat("[3/10] Reading metadata...\n")
metadata <- read_metadata(
    metadata_paths = metadata_paths,
    sample_col = sample_col,
    normalize_names = TRUE
)
# After normalization, sample_col becomes lowercase
sample_col_norm <- tolower(sample_col)
cat(sprintf("✓ Loaded metadata for %d samples\n\n", nrow(metadata)))

# Ensure output directories exist
cat("[4/10] Setting up output directories...\n")
output_dirs <- c(
    "projects/pscc/wes/output/plots/qc",
    "projects/pscc/wes/output/plots/variants",
    "projects/pscc/wes/output/data"
)
for (d in output_dirs) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
cat("✓ Output directories ready\n\n")

# Coverage statistics
cat("[5/10] Processing QC metrics...\n")
cat("  → Reading coverage statistics...\n")
coverage_stats <- data.table::setDT(read_mosdepth_summaries(
    metadata = metadata,
    dirs = dirs_coverage,
    sample_col = sample_col_norm
))
plot_metric_by_sample(
    data = coverage_stats,
    file = "projects/pscc/wes/output/plots/qc/coverage_plot.png",
    metric_col = "mean",
    y_label = "Mean Coverage",
    color_by = "facility",
    threshold = 22,
    threshold_label = "22x",
    width = 10,
    height = 6
)

# Contamination
cat("  → Reading contamination statistics...\n")
contaminations_stats <- data.table::setDT(read_contamination_tables(
    metadata = metadata,
    dirs = dirs_contamination,
    sample_col = sample_col_norm
))
plot_metric_by_sample(
    data = contaminations_stats,
    file = "projects/pscc/wes/output/plots/qc/contamination_plot.png",
    metric_col = "contamination",
    y_label = "Contamination",
    color_by = "facility",
    width = 10,
    height = 6
)
cat("✓ QC metrics processed successfully\n\n")

# ============================================================================
# VARIANT PROCESSING WORKFLOW: Annotation → Counting → Filtering
# ============================================================================
#
# This pipeline processes variants through three main stages:
#
# 1. ANNOTATION (upstream: dbnsfp_annotate.sh)
#    - Extracts fields from dbNSFP v5.1a for 14 tools (configured in
#      config/params.yml, select_columns):
#      * Categorical _pred columns (majority-voted): SIFT, SIFT4G,
#        PolyPhen2-HDIV, PolyPhen2-HVAR, MutationTaster, MutationAssessor,
#        PROVEAN, MetaSVM, ClinPred, AlphaMissense
#      * Numeric _score columns (categorised to D/T/U): REVEL, MutScore,
#        DANN, CADD (column name: CADD_phred)
#    - Plus population frequencies: RegeneronME (ALL/EUR/SPA), gnomAD 4.1
#      (POPMAX and NFE)
#    - COVERAGE: Primarily MISSENSE SNVs (dbNSFP limitation)
#      * SNVs are NOT all missense - they include: missense, nonsense/stop-gain,
#        synonymous, start-loss, splice site variants
#      * Only MISSENSE SNVs get full functional annotation
#      * Stop-gain, start-loss, splice site SNVs: NO dbNSFP data (has_dbnsfp = FALSE)
#    - INDELs and unmatched variants: recovered by dbnsfp_annotate.sh and
#      appended with empty dbNSFP columns (".")
#
# 2. MERGING & COUNTING (merge_variants.R)
#    a) transform_predictions_scores():
#       - Converts 4 numeric score columns → categorical T/D/U:
#           REVEL    (tolerated <= 0.30, damaging >= 0.75)
#           MutScore (tolerated <= 0.30, damaging >= 0.75)
#           DANN     (tolerated <= 0.30, damaging >= 0.96)
#           CADD     (tolerated <= 10,   damaging >= 20)   <- column: cadd_phred
#       - Applies majority vote to 10 categorical _pred columns
#       - Sets has_dbnsfp flag (TRUE if >=1 real prediction/score exists)
#    b) add_num_damaging():
#       - Uses an explicit hardcoded list of 13 columns (deterministic):
#           _pred (9): SIFT, SIFT4G, PolyPhen2-HDIV, PolyPhen2-HVAR,
#                      MutationTaster, MutationAssessor, PROVEAN, MetaSVM,
#                      ClinPred
#           _score (3): REVEL, MutScore, DANN (after categorisation above)
#           _phred (1): CADD (cadd_phred, after categorisation above)
#       - AlphaMissense intentionally excluded (long-form string values)
#       - Damaging values: "D", "A", "M", "H", "LP", "P"
#       - Creates num_damaging column (0-13 range)
#
# 3. FILTERING (filter_variants.R)
#    12 sequential steps. Each step removes variants or entire samples.
#
#    Step 1 — Coverage (sample-level):
#      FFPE matched: 30x, FFPE tumor-only: 40x
#      Fresh matched: 30x, Fresh tumor-only: 40x
#      Entire sample dropped if mean mosdepth coverage < threshold.
#
#    Step 2 — Contamination (sample-level):
#      Entire sample dropped if VerifyBamID2 contamination > 5%.
#
#    Step 3 — Sample exclusion (sample-level):
#      Manually specified sample IDs removed entirely.
#
#    Step 4 — Low coverage samples (variant-level):
#      For samples flagged as borderline-coverage, a variant is REMOVED
#      only if it fails BOTH thresholds: num_callers < 4 AND num_damaging < 6.
#      OR logic: kept if num_callers >= 4 OR num_damaging >= 6.
#
#    Step 5 — VAF (variant-level):
#      SNV:   AF >= 2%  (AF = NA: kept)
#      INDEL: AF >= 4%  (AF = NA: kept)
#
#    Step 6 — Germline frequency (variant-level):
#      Remove if gnomAD v4.1 AF or 1000 Genomes AF >= 0.01% in any
#      population column (gnomAD popmax, NFE; 1000G all sub-populations).
#
#    Step 6.5 — Unique reads (variant-level; applied when UD column present):
#      UD (bam-readcount unique depth) >= 4.
#      Guards against PCR duplicate artefacts; NA UD passes through.
#
#    TIER 1 (Step 7): Baseline caller consensus for ALL variants
#      - Matched samples:    >=3 callers
#      - Tumor-only samples: >=4 callers
#      Impact-specific thresholds are NOT applied at this step;
#      that happens exclusively in Step 9.
#
#    (Step 8: specific_callers filter — disabled here; specific_callers = NULL
#     passes all variants through unchanged)
#
#    TIER 2 (Step 9): Impact-specific requirements — two paths:
#
#      INDEL PATH (treat_coding_indels_as_high = TRUE, default):
#        All coding INDELs (HIGH, MODERATE, LOW) use the HIGH caller threshold.
#        INDELs never have dbNSFP data; caller agreement is the only quality signal.
#        In-frame and splice-region INDELs are as clinically relevant as frameshifts
#        (e.g., EGFR exon 19 del, ERBB2 exon 20 ins). MODIFIER INDELs excluded.
#        - Matched:    requires 4+ callers
#        - Tumor-only: requires 5+ callers
#
#      SNV PATH:
#      a) HIGH impact SNVs (stop-gain/loss, start-loss):
#         - No functional scores available (has_dbnsfp = FALSE)
#         - Requires HIGHER caller count (4-5 callers)
#
#      b) MODERATE impact SNVs + has_dbnsfp = TRUE (missense):
#         - Requires num_damaging >=3 (out of 13 tools)
#         - Already passed baseline callers in Tier 1
#
#      c) MODERATE SNVs + has_dbnsfp = FALSE:
#         - Falls back to HIGHER caller count (same as HIGH SNVs)
#
#      Note on has_dbnsfp: computed from _pred and _score columns ONLY.
#      cadd_phred ends in _phred so it is NOT included in this check, even
#      though it contributes to num_damaging. A variant with only CADD data
#      and all other tools returning "." gets has_dbnsfp = FALSE and falls
#      through to the caller-count path.
#
#    Step 10 — Recurrence blacklist (variant-level):
#      Runs AFTER Steps 1-9 quality + impact filtering (not before).
#      Flags chr:pos:ref:alt combinations present in >30% of samples
#      → likely mapping artefacts / pseudogene reads.
#      Blacklist exported to: data/recurrence_blacklisted_variants.tsv
#
#    Step 11 — Cancer gene filter (variant-level):
#      Keeps only variants where is_oncogene, is_tumor_suppressor_gene,
#      or is_driver = TRUE/Yes (OncoKB + Compendium Cancer Genes).
#
#    Step 12 — Artifact removal (variant-level):
#      SOBDetector: removes sob_prediction = 'artifact'
#      FFPErase:    removes ffperase_prediction = 'True'
#
# EXAMPLE OUTCOMES:
#   v Coding INDEL:  4 callers → PASS (INDEL path, HIGH threshold)
#   v Missense SNV: 4 callers + 8/13 damaging → PASS (functional evidence)
#   v Frameshift SNV: 5 callers + 0 damaging → PASS (high caller support)
#   x Coding INDEL:  3 callers → FAIL (below HIGH threshold)
#   x Missense SNV: 4 callers + 2/13 damaging → FAIL (weak evidence)
#   x Frameshift SNV: 3 callers + 0 damaging → FAIL (low caller support)
#
# See filter_variants_workflow.md for full documentation of all 12 steps
# ============================================================================

# Merge variants
cat("[6/10] Processing variants data...\n")
# Step 1: Read and combine VCFs
# File-based caching allows resuming from checkpoints
if (!exists("combined_vcf")) {
    if (file.exists("projects/pscc/wes/output/data/combined_vcf.fst")) {
        cat("  → Loading combined VCF from checkpoint...\n")
        message("Resuming: Loading combined_vcf from checkpoint...")
        combined_vcf <- data.table::setDT(fst::read.fst(
            "projects/pscc/wes/output/data/combined_vcf.fst"
        ))
    } else {
        cat("  → Reading and combining VCF files...\n")
        message("Reading and combining VCF files...")
        combined_vcf <- read_and_combine_vcfs(
            sample_dirs = dirs_variants,
            vcf_extension = ".merged.annot.tsv",
            caller_names = c(),
            sample_id_col = sample_col_norm
        )
        fst::write.fst(
            combined_vcf,
            "projects/pscc/wes/output/data/combined_vcf.fst"
        )
    }
    cat(sprintf(
        "combined_vcf: %d rows, %.1f MB\n",
        nrow(combined_vcf),
        as.numeric(object.size(combined_vcf)) / 1024^2
    ))
}

# Step 2: Parse INFO and calculate AF
if (!exists("parsed_vcf")) {
    if (file.exists("projects/pscc/wes/output/data/parsed_vcf.fst")) {
        cat("  → Loading parsed VCF from checkpoint...\n")
        message("Resuming: Loading parsed_vcf from checkpoint...")
        parsed_vcf <- data.table::setDT(fst::read.fst(
            "projects/pscc/wes/output/data/parsed_vcf.fst"
        ))
    } else {
        if (!exists("combined_vcf")) {
            combined_vcf <- data.table::setDT(fst::read.fst(
                "projects/pscc/wes/output/data/combined_vcf.fst"
            ))
        }
        cat("  → Parsing VCF INFO fields...\n")
        message("Parsing VCF INFO fields...")
        parsed_vcf <- parse_vcf_info(
            combined_vcf,
            sample_id_col = sample_col_norm
        )
        fst::write.fst(
            parsed_vcf,
            "projects/pscc/wes/output/data/parsed_vcf.fst"
        )
    }
    # Drop leftover genotype/sample columns from old checkpoints (if any)
    # These are per-sample per-caller FORMAT columns only needed for AF calculation
    keep_cols <- !grepl("^(\\d+:)?[A-Z]{2,}\\d{4,}", names(parsed_vcf))
    drop_cols <- names(parsed_vcf)[!keep_cols]
    extra_drops <- intersect(c("INFO", "info", "FORMAT"), names(parsed_vcf))
    drop_cols <- c(drop_cols, extra_drops)
    if (length(drop_cols) > 0) {
        cat(sprintf(
            "  → Dropping %d leftover genotype columns from checkpoint\n",
            length(drop_cols)
        ))
        parsed_vcf[, (drop_cols) := NULL]
        gc()
    }
    cat(sprintf(
        "parsed_vcf: %d rows, %.1f MB\n",
        nrow(parsed_vcf),
        as.numeric(object.size(parsed_vcf)) / 1024^2
    ))
}
# Free combined_vcf now that parsing is complete
if (exists("combined_vcf")) {
    rm(combined_vcf)
    gc()
}

# Step 3: Merge with metadata
if (!exists("merged_data")) {
    if (file.exists("projects/pscc/wes/output/data/merged_with_metadata.fst")) {
        cat("  → Loading merged data from checkpoint...\n")
        message("Resuming: Loading merged_with_metadata from checkpoint...")
        merged_data <- data.table::setDT(fst::read.fst(
            "projects/pscc/wes/output/data/merged_with_metadata.fst"
        ))
    } else {
        if (!exists("parsed_vcf")) {
            parsed_vcf <- data.table::setDT(fst::read.fst(
                "projects/pscc/wes/output/data/parsed_vcf.fst"
            ))
        }
        cat("  → Merging variants with metadata...\n")
        message("Merging with metadata...")
        # merge_with_metadata modifies in-place (no copy) to save memory
        merged_data <- merge_with_metadata(
            variant_data = parsed_vcf,
            metadata = metadata,
            sample_id_col = sample_col_norm
        )
        # parsed_vcf and merged_data now point to the same data.table
        rm(parsed_vcf)
        gc()
        print(gc()) # Log memory state after merge
        fst::write.fst(
            merged_data,
            "projects/pscc/wes/output/data/merged_with_metadata.fst"
        )
    }
    cat(sprintf(
        "merged_data: %d rows, %.1f MB\n",
        nrow(merged_data),
        as.numeric(object.size(merged_data)) / 1024^2
    ))
} else {
    # Free parsed_vcf if merged_data already existed
    if (exists("parsed_vcf")) {
        rm(parsed_vcf)
        gc()
    }
}

# Step 4: Annotate with cancer genes
if (!exists("annotated_data")) {
    if (file.exists("projects/pscc/wes/output/data/annotated_data.fst")) {
        cat("  → Loading annotated data from checkpoint...\n")
        message("Resuming: Loading annotated_data from checkpoint...")
        annotated_data <- data.table::setDT(fst::read.fst(
            "projects/pscc/wes/output/data/annotated_data.fst"
        ))
    } else {
        if (!exists("merged_data")) {
            merged_data <- data.table::setDT(fst::read.fst(
                "projects/pscc/wes/output/data/merged_with_metadata.fst"
            ))
        }
        cat("  → Annotating with cancer genes...\n")
        message("Annotating with cancer genes...")
        # annotate_cancer_genes modifies in-place (no copy) to save memory
        annotated_data <- annotate_cancer_genes(
            variant_data = merged_data,
            oncokb_path = "references/oncoKB/cancerGeneList.tsv",
            compendium_path = "references/bbg_lab/Compendium_Cancer_Genes.tsv"
        )
        # merged_data and annotated_data now point to the same data.table
        rm(merged_data)
        gc()
        fst::write.fst(
            annotated_data,
            "projects/pscc/wes/output/data/annotated_data.fst"
        )
    }
    cat(sprintf(
        "annotated_data: %d rows, %.1f MB\n",
        nrow(annotated_data),
        as.numeric(object.size(annotated_data)) / 1024^2
    ))
} else {
    # Free merged_data if annotated_data already existed
    if (exists("merged_data")) {
        rm(merged_data)
        gc()
    }
}

# Step 5: Transform predictions and scores
if (!exists("merged_variants")) {
    if (file.exists("projects/pscc/wes/output/data/merged_variants.fst")) {
        cat("  → Loading merged variants from checkpoint...\n")
        message("Resuming: Loading merged_variants from checkpoint...")
        merged_variants <- data.table::setDT(fst::read.fst(
            "projects/pscc/wes/output/data/merged_variants.fst"
        ))
    } else {
        if (!exists("annotated_data")) {
            annotated_data <- data.table::setDT(fst::read.fst(
                "projects/pscc/wes/output/data/annotated_data.fst"
            ))
        }
        cat("  → Transforming predictions and scores...\n")
        message("Transforming predictions and scores...")
        merged_variants <- transform_predictions_scores(
            variant_data = annotated_data,
            transform_scores = TRUE,
            score_thresholds = list(
                revel_score = list(tolerated = 0.3, damaging = 0.75),
                mut_score_score = list(tolerated = 0.3, damaging = 0.75),
                dann_score = list(tolerated = 0.3, damaging = 0.96),
                cadd_phred = list(tolerated = 10, damaging = 20)
            ),
            transform_predictions = TRUE
        )
        fst::write.fst(
            merged_variants,
            "projects/pscc/wes/output/data/merged_variants.fst"
        )
    }
    cat(sprintf(
        "merged_variants: %d rows, %.1f MB\n",
        nrow(merged_variants),
        as.numeric(object.size(merged_variants)) / 1024^2
    ))
}
# Free annotated_data now that transformation is complete
if (exists("annotated_data")) {
    rm(annotated_data)
    gc()
}

# Read and merge FFPErase predictions (if available)
# if (!("ffperase_prediction" %in% names(merged_variants))) {
#     cat("  → Merging FFPErase predictions...\n")
#     message("Merging FFPErase predictions...")
#     ffperase_results <- data.table::setDT(read_ffperase_predictions(
#         metadata = metadata,
#         dirs = dirs_ffperase,
#         sample_col = sample_col_norm
#     ))
#     merged_variants <- data.table::setDT(merge_ffperase_predictions(
#         variant_data = merged_variants,
#         ffperase_data = ffperase_results,
#         sample_col = sample_col_norm
#     ))
#     # Re-save merged_variants with ffperase predictions
#     fst::write.fst(
#         merged_variants,
#         "projects/pscc/wes/output/data/merged_variants.fst"
#     )
#     # Free ffperase_results
#     rm(ffperase_results)
#     gc()
# } else {
#     cat("  → FFPErase predictions already present\n")
#     message("FFPErase predictions already present in merged_variants")
# }
cat("✓ Variants data processing complete\n\n")

# Artifact detection summary
cat("[7/10] Analyzing artifact predictions...\n")
artifact_counts <- list(total_variants = nrow(merged_variants))
if ("sob_prediction" %in% names(merged_variants)) {
    artifact_counts$sob_artifacts <- sum(
        !is.na(merged_variants$sob_prediction) &
            merged_variants$sob_prediction == "artifact",
        na.rm = TRUE
    )
} else {
    artifact_counts$sob_artifacts <- NA_integer_
}
if ("ffperase_prediction" %in% names(merged_variants)) {
    artifact_counts$ffperase_artifacts <- sum(
        !is.na(merged_variants$ffperase_prediction) &
            merged_variants$ffperase_prediction == "True",
        na.rm = TRUE
    )
} else {
    artifact_counts$ffperase_artifacts <- NA_integer_
}
if (
    "sob_prediction" %in%
        names(merged_variants) &&
        "ffperase_prediction" %in% names(merged_variants)
) {
    artifact_counts$both_artifacts <- sum(
        !is.na(merged_variants$sob_prediction) &
            !is.na(merged_variants$ffperase_prediction) &
            merged_variants$sob_prediction == "artifact" &
            merged_variants$ffperase_prediction == "True",
        na.rm = TRUE
    )
} else {
    artifact_counts$both_artifacts <- NA_integer_
}
artifact_summary <- data.table::as.data.table(artifact_counts)
cat(sprintf("  → Total variants: %d\n", artifact_counts$total_variants))
cat(sprintf(
    "  → SOB artifacts: %s\n",
    ifelse(
        is.na(artifact_counts$sob_artifacts),
        "N/A",
        artifact_counts$sob_artifacts
    )
))
cat(sprintf(
    "  → FFPErase artifacts: %s\n",
    ifelse(
        is.na(artifact_counts$ffperase_artifacts),
        "N/A",
        artifact_counts$ffperase_artifacts
    )
))
cat(sprintf(
    "  → Both tools: %s\n",
    ifelse(
        is.na(artifact_counts$both_artifacts),
        "N/A",
        artifact_counts$both_artifacts
    )
))
cat("✓ Artifact analysis complete\n\n")

# Filter variants before oncoplot
cat("[8/10] Filtering variants for list 1...\n")
if (!exists("filtered_variants_list1")) {
    if (file.exists("projects/pscc/wes/output/data/filtered_variants_list1.fst")) {
        cat("  → Loading filtered variants from checkpoint...\n")
        message("Resuming: Loading filtered_variants from checkpoint...")
        filtered_variants1 <- data.table::setDT(fst::read.fst(
            "projects/pscc/wes/output/data/filtered_variants_list1.fst"
        ))
    } else {
        # Ensure merged_variants is available (for jumping to this section)
        if (!exists("merged_variants")) {
            merged_variants <- data.table::setDT(fst::read.fst(
                "projects/pscc/wes/output/data/merged_variants.fst"
            ))
        }
        cat("  → Applying variant filters...\n")
        message("Filtering variants for list 1...")
        # FILTERING STRATEGY (see header comments for full workflow):
        # - Two-tier caller filtering: baseline (3-4 callers) + impact-specific (4-5 callers)
        # - INDEL PATH (treat_coding_indels_as_high = TRUE, default):
        #     All coding INDELs (HIGH/MODERATE/LOW) → HIGH caller threshold (4+ matched, 5+ tumor-only)
        #     MODIFIER INDELs (intronic/intergenic/UTR) excluded
        # - SNV PATH:
        #     HIGH SNVs (stop-gain): filtered by caller count (4-5 callers)
        #     MODERATE SNVs with dbNSFP (missense): filtered by damaging predictions (>=3/13 tools)
        #     MODERATE SNVs without dbNSFP: filtered by caller count
        # - Recurrence blacklist runs AFTER impact filtering (Step 10) to identify
        #   systematic artifacts among high-quality variants
        filtered_variants1 <- filter_variants(
            variant_data = merged_variants,
            sample_col = sample_col_norm,
            metadata = metadata,
            coverage_dirs = dirs_coverage,
            contamination_dirs = dirs_contamination,
            export_blacklisted_path = "projects/pscc/wes/output/data/recurrence_blacklisted_variants.tsv",
            export_low_coverage_path = "projects/pscc/wes/output/data/low_coverage_samples_mutational_load.tsv",
            export_passed_coverage_path = "projects/pscc/wes/output/data/passed_coverage_samples.tsv",
            export_coverage_path = "projects/pscc/wes/output/data/sample_coverage_stats.tsv",
            export_contamination_path = "projects/pscc/wes/output/data/sample_contamination_stats.tsv",
           filters = list(
                min_callers = 3,
                min_callers_high_impact = 3,
                min_callers_tumor_only = 3,
                min_callers_high_impact_tumor_only = 3,
                min_damaging = 4, 
                allowed_impacts = c("HIGH", "MODERATE"),
                require_cancer_gene = TRUE,
                remove_artifacts = TRUE,
                min_dp = 100, #variants over that value pass
                min_vaf_snv = 0.025, # 2% VAF threshold for SNVs
                min_vaf_indel = 0.04, # 4% VAF threshold for indels
                max_gnomad_af = 0.0001, # 0.01% MAF threshold for gnomAD POPMAX
                max_1000g_af = 0.0001, # 0.01% MAF threshold for 1000G EUR
                min_unique_reads = 4, # Minimum unique (non-duplicate) reads - protects against PCR artifacts
                max_recurrence_fraction = 0.3, # Blacklist variants in >30% of samples (mapping artifacts)
                exclude_samples = c("AS9750"), # Add sample IDs to exclude if needed
                low_cov_samples = NULL, # Add low coverage sample IDs if needed
                low_cov_min_callers = 3,
                low_cov_min_damaging = 4,
                # Sample-type-specific coverage filtering
                min_coverage_ffpe_matched = 22, # FFPE matched samples
                min_coverage_ffpe_tumor_only = 40, # FFPE tumor-only
                min_coverage_fresh_matched = 30, # Fresh matched
                min_coverage_fresh_tumor_only = 40, # Fresh tumor-only
                # Contamination filtering (optional)
                max_contamination = 0.05# e.g., 0.05 to exclude samples with >5% contamination
            ),
            verbose = TRUE
        )
        fst::write.fst(
            filtered_variants1,
            "projects/pscc/wes/output/data/filtered_variants_list1.fst"
        )
    }
    cat(sprintf(
        "filtered_variants: %d rows, %.1f MB\n",
        nrow(filtered_variants1),
        as.numeric(object.size(filtered_variants1)) / 1024^2
    ))
}
cat("✓ Variant filtering for list 1 complete\n\n")


cat("[8.2/10] Filtering variants for list 2...\n")
if (!exists("filtered_variants_list2")) {
    if (file.exists("projects/pscc/wes/output/data/filtered_variants_list2.fst")) {
        cat("  → Loading filtered variants from checkpoint...\n")
        message("Resuming: Loading filtered_variants from checkpoint...")
        filtered_variants2 <- data.table::setDT(fst::read.fst(
            "projects/pscc/wes/output/data/filtered_variants_list2.fst"
        ))
    } else {
        # Ensure merged_variants is available (for jumping to this section)
        if (!exists("merged_variants")) {
            merged_variants <- data.table::setDT(fst::read.fst(
                "projects/pscc/wes/output/data/merged_variants.fst"
            ))
        }
        cat("  → Applying variant filters...\n")
        message("Filtering variants for list 2...")
        # FILTERING STRATEGY (see header comments for full workflow):
        # - Two-tier caller filtering: baseline (3-4 callers) + impact-specific (4-5 callers)
        # - INDEL PATH (treat_coding_indels_as_high = TRUE, default):
        #     All coding INDELs (HIGH/MODERATE/LOW) → HIGH caller threshold (4+ matched, 5+ tumor-only)
        #     MODIFIER INDELs (intronic/intergenic/UTR) excluded
        # - SNV PATH:
        #     HIGH SNVs (stop-gain): filtered by caller count (4-5 callers)
        #     MODERATE SNVs with dbNSFP (missense): filtered by damaging predictions (>=3/13 tools)
        #     MODERATE SNVs without dbNSFP: filtered by caller count
        # - Recurrence blacklist runs AFTER impact filtering (Step 10) to identify
        #   systematic artifacts among high-quality variants
        filtered_variants2 <- filter_variants(
            variant_data = merged_variants,
            sample_col = sample_col_norm,
            metadata = metadata,
            coverage_dirs = dirs_coverage,
            contamination_dirs = dirs_contamination,
            export_blacklisted_path = "projects/pscc/wes/output/data/recurrence_blacklisted_variants.tsv",
            export_low_coverage_path = "projects/pscc/wes/output/data/low_coverage_samples_mutational_load.tsv",
            export_passed_coverage_path = "projects/pscc/wes/output/data/passed_coverage_samples.tsv",
            export_contamination_path = "projects/pscc/wes/output/data/sample_contamination_stats.tsv",
            filters = list(
                min_callers = 2,
                min_callers_high_impact = 2,
                min_callers_tumor_only = 2,
                min_callers_high_impact_tumor_only = 2,
                min_damaging = 7, 
                allowed_impacts = c("HIGH", "MODERATE"),
                required_caller_any = c("freebayes", "lofreq", "mutect2", "deepsomatic", "strelka"), 
                require_cosmic_id = TRUE,
                require_cancer_gene = TRUE,
                remove_artifacts = TRUE,
                min_dp = 100, #variants over that value pass
                min_vaf_snv = 0.025, # 2% VAF threshold for SNVs
                min_vaf_indel = 0.04, # 4% VAF threshold for indels
                max_gnomad_af = 0.0001, # 0.01% MAF threshold for gnomAD POPMAX
                max_1000g_af = 0.0001, # 0.01% MAF threshold for 1000G EUR
                min_unique_reads = 4, # Minimum unique (non-duplicate) reads - protects against PCR artifacts
                max_recurrence_fraction = 0.3, # Blacklist variants in >30% of samples (mapping artifacts)
                exclude_samples = c("AS9750"), # Add sample IDs to exclude if needed
                low_cov_samples = NULL, # Add low coverage sample IDs if needed
                low_cov_min_callers = 2,
                low_cov_min_damaging = 7,
                # Sample-type-specific coverage filtering
                min_coverage_ffpe_matched = 22, # FFPE matched samples
                min_coverage_ffpe_tumor_only = 40, # FFPE tumor-only
                min_coverage_fresh_matched = 30, # Fresh matched
                min_coverage_fresh_tumor_only = 40, # Fresh tumor-only
                # Contamination filtering (optional)
                max_contamination = 0.05# e.g., 0.05 to exclude samples with >5% contamination
            ),
            verbose = TRUE
        )
        fst::write.fst(
            filtered_variants2,
            "projects/pscc/wes/output/data/filtered_variants_list2.fst"
        )
    }
    cat(sprintf(
        "filtered_variants: %d rows, %.1f MB\n",
        nrow(filtered_variants2),
        as.numeric(object.size(filtered_variants2)) / 1024^2
    ))
}
cat("✓ Variant filtering for list 2 complete\n\n")

# Now we check both outputs filtered_variants1 and filtered_variants2 to see
# how many variants are in each list and how many overlap.
# This will help us understand the impact of the different filtering strategies
# on our final variant list.
cat("[8.3/10] Comparing filtered variant lists...\n")
if (exists("filtered_variants1") && exists("filtered_variants2")) {
    # Avoid duplicated variants when binding the two filtered lists.
    duplicate_cols <- intersect(names(filtered_variants2), names(filtered_variants1))
    if (length(duplicate_cols) == 0) {
        stop("No shared columns found between filtered_variants1 and filtered_variants2.")
    }

    filtered_variants_common <- dplyr::semi_join(
        filtered_variants1,
        filtered_variants2,
        by = duplicate_cols
    )

    filtered_variants1_only <- dplyr::anti_join(
        filtered_variants1,
        filtered_variants2,
        by = duplicate_cols
    )

    filtered_variants2_only <- dplyr::anti_join(
        filtered_variants2,
        filtered_variants1,
        by = duplicate_cols
    )

    filtered_variants <- dplyr::bind_rows(
        filtered_variants1,
        filtered_variants2_only
    )
    n_combined_before_dedup <- nrow(filtered_variants)
    setDT(filtered_variants)
    filtered_variants <- unique(filtered_variants, by = duplicate_cols)

    cat(sprintf("  -> filtered_variants1 rows: %d\n", nrow(filtered_variants1)))
    cat(sprintf("  -> filtered_variants2 rows: %d\n", nrow(filtered_variants2)))
    cat(sprintf("  -> rows in common: %d\n", nrow(filtered_variants_common)))
    cat(sprintf("  -> rows only in filtered_variants1: %d\n", nrow(filtered_variants1_only)))
    cat(sprintf("  -> rows only in filtered_variants2: %d\n", nrow(filtered_variants2_only)))
    cat(sprintf(
        "  -> duplicate rows removed from combined list: %d\n",
        n_combined_before_dedup - nrow(filtered_variants)
    ))

    # Remove recurrent artefacts observed in 5 or more samples.
    cat("  -> Removing recurrent artefacts observed in 5 or more samples...\n")

    list_analyzed <- filtered_variants

    # Variant identifier columns.
    id_cols <- c("number_chrom", "pos", "ref", "alt")
    recurrence_required_cols <- c(id_cols, "sample_id")
    missing_recurrence_cols <- setdiff(recurrence_required_cols, names(list_analyzed))
    if (length(missing_recurrence_cols) > 0) {
        stop(sprintf(
            "Missing columns required for recurrent artefact filtering: %s",
            paste(missing_recurrence_cols, collapse = ", ")
        ))
    }

    # Count unique samples per variant.
    recurrence_counts <- list_analyzed[
        ,
        .(n_samples = uniqueN(sample_id)),
        by = id_cols
    ]

    # Create a blacklist of recurrently observed variants.
    blacklisted <- recurrence_counts[
        n_samples >= 5
    ]

    # Export blacklisted variants.
    recurrently_observed_variants <- list_analyzed[
        blacklisted,
        on = id_cols,
        nomatch = 0
    ]

    # Remove blacklisted variants.
    list_analyzed_filtered <- list_analyzed[
        !blacklisted,
        on = id_cols
    ]

    cat(sprintf("  -> recurrent variant sites blacklisted: %d\n", nrow(blacklisted)))
    cat(sprintf(
        "  -> rows removed as recurrent artefacts: %d\n",
        nrow(list_analyzed) - nrow(list_analyzed_filtered)
    ))

    filtered_variants <- list_analyzed_filtered

    cat(sprintf(
        "filtered_variants: %d rows, %.1f MB\n",
        nrow(filtered_variants),
        as.numeric(object.size(filtered_variants)) / 1024^2
    ))

    fst::write.fst(
        filtered_variants,
        "projects/pscc/wes/output/data/filtered_variants.fst"
    )
} else {
    stop("filtered_variants1 and filtered_variants2 must exist before combining.")
}
cat("Filtered variant lists combined\n\n")

# Optional INDEL diagnostic analysis
# Tracks SNV/INDEL counts through each filter step to identify where INDELs
# are lost. Requires ~100GB RAM. Set to TRUE to run.
run_indel_diagnostics <- FALSE
if (run_indel_diagnostics) {
    cat("[8b/10] Running INDEL diagnostic analysis...\n")
    source("projects/pscc/wes/output/diagnose_indels.R")
    cat("✓ INDEL diagnostic analysis complete\n\n")
}

# Ensure filtered_variants is available for downstream plots (if jumping to this section)
cat("[9/10] Generating visualization plots...\n")
if (!exists("filtered_variants")) {
    message("Loading filtered_variants for QC plots...")
    filtered_variants <- data.table::setDT(fst::read.fst(
        "projects/pscc/wes/output/data/filtered_variants.fst"
    ))
}

# Variants classification (using combined filtered variants for QC)
cat("  → Creating SnpEff annotation plot...\n")
plot_snpeff_annotations(
    variant_data = filtered_variants,
    file = "projects/pscc/wes/output/plots/variants/snpeff_annotation_plot.png",
    sample_col = sample_col_norm,
    annotation_col = "annotation",
    top_n = 10,
    width = 10,
    height = 6
)

# Variant callers overlap (using combined filtered variants for QC)
# Dynamically detect caller columns (those ending with _called)
caller_cols <- grep("_called$", names(filtered_variants), value = TRUE)
callers_detected <- gsub("_called$", "", caller_cols)
if (length(callers_detected) > 0) {
    cat(sprintf(
        "  → Detected %d callers: %s\n",
        length(callers_detected),
        paste(callers_detected, collapse = ", ")
    ))
    message(sprintf(
        "Detected %d callers: %s",
        length(callers_detected),
        paste(callers_detected, collapse = ", ")
    ))

    # ==============================================================================
    # Plot 1: Variant Distribution by Number of Callers (Bar Plot)
    # ==============================================================================
    if ("num_callers" %in% names(filtered_variants)) {
        # Prepare data for caller distribution plot
        caller_dist <- as.data.frame(table(filtered_variants$num_callers))
        colnames(caller_dist) <- c("num_callers", "count")
        caller_dist$num_callers <- as.integer(as.character(
            caller_dist$num_callers
        ))
        caller_dist <- caller_dist[order(caller_dist$num_callers), ]

        # Create plot with tinyplot
        png(
            filename = "projects/pscc/wes/output/plots/variants/caller_distribution_plot.png",
            width = 10,
            height = 6,
            units = "in",
            res = 300
        )

        # Set up plotting parameters
        cat("  → Creating caller distribution plot...\n")
        # Set up plotting parameters
        par(
            mar = c(5, 5, 3, 2),
            family = "sans"
        )

        # Create barplot using base R (tinyplot alternative)
        bp <- barplot(
            caller_dist$count,
            names.arg = caller_dist$num_callers,
            xlab = "Number of Callers",
            ylab = "Number of Variants",
            main = "Variant Distribution by Number of Callers",
            col = "#5B8FA8",
            border = "white",
            las = 1,
            cex.lab = 1.2,
            cex.axis = 1.1,
            cex.main = 1.3
        )

        # Add count labels on top of bars
        text(
            x = bp,
            y = caller_dist$count,
            labels = format(caller_dist$count, big.mark = ","),
            pos = 3,
            cex = 0.9,
            col = "black"
        )

        # Add grid for readability
        grid(nx = NA, ny = NULL, col = "gray90", lty = 1)

        # Redraw bars on top of grid
        barplot(
            caller_dist$count,
            names.arg = caller_dist$num_callers,
            col = "#5B8FA8",
            border = "white",
            las = 1,
            add = TRUE
        )

        dev.off()

        message(sprintf(
            "Created caller distribution plot with %d unique caller counts",
            nrow(caller_dist)
        ))
    }

    # ==============================================================================
    # Plot 2: Caller Overlap Upset Plot
    # ==============================================================================
    cat("  → Creating caller overlap upset plot...\n")
    overlap_callers_plot <- plot_overlap_callers(
        merged_variants = filtered_variants,
        caller_names = callers_detected,
        title = NULL,
        fill_alpha = 0.5,
        stroke_size = 1,
        label_size = 4,
        set_name_size = 5
    )
} else {
    warning("No caller columns detected. Skipping overlap plot.")
    overlap_callers_plot <- NULL
}
if (!is.null(overlap_callers_plot)) {
    png(
        filename = "projects/pscc/wes/output/plots/variants/overlap_callers_plot.png",
        width = 10,
        height = 6,
        units = "in",
        res = 300
    )
    print(overlap_callers_plot)
    dev.off()
}

cat("✓ All visualization plots generated\n\n")

# Free merged_variants before loading filtered_variants to avoid double memory usage
if (exists("merged_variants")) {
    rm(merged_variants)
    gc()
}

# Ensure filtered_variants is available for MAF analysis (if jumping to this section)
cat("[10/10] Generating MAF analysis and oncoplots...\n")
if (!exists("filtered_variants")) {
    cat("  → Loading filtered variants for MAF analysis...\n")
    message("Loading filtered_variants for MAF analysis...")
    filtered_variants <- data.table::setDT(fst::read.fst(
        "projects/pscc/wes/output/data/filtered_variants.fst"
    ))
}

# MAF conversion for maftools using convert_to_maf helper
# First, deduplicate variants to avoid counting the same variant multiple times
# (e.g., when multiple transcripts/annotations exist for the same variant)
cat("  → Deduplicating variants for MAF conversion...\n")
dedup_cols <- c(sample_col_norm, "number_chrom", "pos", "ref", "alt", "gene")
filtered_variants_unique <- unique(
    filtered_variants,
    by = dedup_cols
)

message(sprintf(
    "Deduplicated variants: %d → %d (removed %d duplicate annotations)",
    nrow(filtered_variants),
    nrow(filtered_variants_unique),
    nrow(filtered_variants) - nrow(filtered_variants_unique)
))

cat("  → Preparing clinical data...\n")
clinical_cols <- intersect(
    c(
        sample_col_norm,
        "tissue_type",
        "preservation",
        "facility",
        "patient",
        "matched"
    ),
    names(filtered_variants_unique)
)
clinical_data <- unique(filtered_variants_unique[, ..clinical_cols])

cat("  → Converting to MAF format...\n")
maf_maftools <- convert_to_maf(
    variant_data = filtered_variants_unique,
    clinical_data = clinical_data,
    sample_col = sample_col_norm,
    annotation_col = "annotation",
    create_maf_object = TRUE
)

# Free large intermediate objects no longer needed
# Keep clinical_data — it's small and needed later for oncoplot annotations
rm(filtered_variants_unique)
gc()

# ============================================================================
# MAF Summary Visualizations
# ============================================================================

cat("  → Generating MAF summary plot...\n")
# MAF Summary Plot - dashboard view of mutation landscape
png(
    filename = "projects/pscc/wes/output/plots/variants/maf_summary.png",
    width = 10,
    height = 8,
    units = "in",
    res = 300
)
maftools::plotmafSummary(
    maf = maf_maftools,
    rmOutlier = TRUE,
    addStat = "median",
    dashboard = TRUE,
    titvRaw = FALSE
)
dev.off()

cat("  → Calculating transition/transversion ratios...\n")
# Transition and Transversions
titv_summary <- maftools::titv(
    maf = maf_maftools,
    plot = FALSE,
    useSyn = TRUE
)

cat("  → Plotting TiTv ratios...\n")
# Plot TiTv by sample
png(
    filename = "projects/pscc/wes/output/plots/variants/titv_plot.png",
    width = 10,
    height = 6,
    units = "in",
    res = 300
)
maftools::plotTiTv(res = titv_summary)
dev.off()

cat("  → Creating MAF barplot...\n")
# MAF barplot - variant classification per sample
png(
    filename = "projects/pscc/wes/output/plots/variants/mafbarplot.png",
    width = 12,
    height = 6,
    units = "in",
    res = 300
)
maftools::plotmafSummary(
    maf = maf_maftools,
    rmOutlier = TRUE,
    addStat = "median",
    dashboard = FALSE,
    titvRaw = FALSE
)
dev.off()

# ============================================================================
# Enhanced Oncoplots with VAF and Clinical Features
# ============================================================================

cat("  → Preparing enhanced oncoplot data...\n")
# Define number of genes to display in oncoplot
n_genes <- 15

# Get gene summary to understand ranking
gene_summary <- maftools::getGeneSummary(maf_maftools)
cat("\n=== Top 20 Genes by MAF Summary ===\n")
print(gene_summary[1:20, .(Hugo_Symbol, AlteredSamples, MutatedSamples, total)])

# Calculate mean VAF per gene for left bar plot
top_genes <- gene_summary[1:n_genes, Hugo_Symbol]

# Extract VAF data for the bar plot (following vignette section 0.4)
# VAF is stored as i_VAF (percentage) in our MAF object
vaf_col <- "i_VAF"
if (vaf_col %in% maftools::getFields(maf_maftools)) {
    gene_vaf <- maftools::subsetMaf(
        maf = maf_maftools,
        genes = top_genes,
        fields = vaf_col,
        mafObj = FALSE
    )[, .(mean(i_VAF, na.rm = TRUE)), by = Hugo_Symbol]
    colnames(gene_vaf) <- c("gene", "VAF")
    message(sprintf("Calculated mean VAF for %d genes", nrow(gene_vaf)))
} else {
    gene_vaf <- NULL
    warning(
        "VAF column (i_VAF) not found in MAF data. VAF bar plot will be skipped."
    )
}

# Prepare clinical features for annotation
available_features <- intersect(
    names(clinical_data),
    c("tissue_type", "preservation", "facility", "matched")
)

# Define muted color palettes for clinical annotations
muted_palettes <- list(
    palette1 = c("#5B8FA8", "#A85B5B", "#7A9B76", "#8B7355", "#6B8E9B"),
    palette2 = c("#8B7355", "#6B8E9B", "#9B8E6B", "#7B6B8E", "#8E7B6B"),
    palette3 = c("#7B6B8E", "#8E7B6B", "#6B8E7B", "#5D7A5D", "#7A5D5D"),
    palette4 = c("#5D7A5D", "#7A5D5D", "#6B7A8E", "#8E6B7A", "#7A8E6B")
)

# Dynamically generate colors based on actual values in clinical_data
annotation_colors <- list()
palette_idx <- 1
for (feat in available_features) {
    if (feat %in% names(clinical_data)) {
        unique_vals <- unique(as.character(clinical_data[[feat]]))
        unique_vals <- unique_vals[!is.na(unique_vals)]
        if (length(unique_vals) > 0) {
            # Get palette and cycle if needed
            pal <- muted_palettes[[palette_idx]]
            n_vals <- length(unique_vals)
            colors <- rep(pal, ceiling(n_vals / length(pal)))[1:n_vals]
            names(colors) <- unique_vals
            annotation_colors[[feat]] <- colors
            palette_idx <- (palette_idx %% length(muted_palettes)) + 1
        }
    }
}

cat("  → Generating enhanced oncoplot...\n")
# Enhanced oncoplot with TiTv, VAF bars, and multiple clinical features
png(
    filename = "projects/pscc/wes/output/plots/variants/oncoplot_enhanced.png",
    width = 24,
    height = 16,
    units = "in",
    res = 300
)
maftools::oncoplot(
    maf = maf_maftools,
    top = n_genes,
    draw_titv = FALSE,
    clinicalFeatures = available_features,
    annotationColor = annotation_colors,
    sortByAnnotation = FALSE, # Sort by mutation frequency, not by annotation
    leftBarData = if (!is.null(gene_vaf)) gene_vaf else NULL,
    leftBarLims = if (!is.null(gene_vaf)) c(0, 100) else NULL,
    showTumorSampleBarcodes = TRUE,
    removeNonMutated = TRUE,
    drawRowBar = TRUE,
    drawColBar = TRUE,
    fontSize = 0.6,
    SampleNamefontSize = 0.4,
    legendFontSize = 1.2,
    annotationFontSize = 1.2,
    gene_mar = 5,
    barcode_mar = 4
)
dev.off()

# ============================================================================
# EXOME-WIDE ONCOPLOT (all genes, no cancer gene filter)
# Same quality filters as above but with require_cancer_gene = FALSE (Step 11
# disabled), so all exome-captured genes are considered.
# ============================================================================

cat(
    "  → Filtering variants with list 1 filters for exome-wide oncoplot (no cancer gene restriction)...\n"
)
if (!exists("filtered_variants_exome1")) {
    if (
        file.exists("projects/pscc/wes/output/data/filtered_variants_exome1.fst")
    ) {
        cat("  → Loading exome-filtered variants from checkpoint...\n")
        message("Resuming: Loading filtered_variants_exome from checkpoint...")
        filtered_variants_exome1 <- data.table::setDT(fst::read.fst(
            "projects/pscc/wes/output/data/filtered_variants_exome1.fst"
        ))
    } else {
        if (!exists("merged_variants")) {
            merged_variants <- data.table::setDT(fst::read.fst(
                "projects/pscc/wes/output/data/merged_variants.fst"
            ))
        }
        cat(
            "  → Applying exome-wide variant filters (require_cancer_gene = FALSE)...\n"
        )
        message("Filtering variants for exome-wide oncoplot...")
        filtered_variants_exome1 <- filter_variants(
            variant_data = merged_variants,
            sample_col = sample_col_norm,
            metadata = metadata,
            coverage_dirs = dirs_coverage,
            contamination_dirs = dirs_contamination,
            export_blacklisted_path = NULL,
            export_low_coverage_path = NULL,
            export_coverage_path = NULL,
            filters = list(
                min_callers = 3,
                min_callers_high_impact = 3,
                min_callers_tumor_only = 3,
                min_callers_high_impact_tumor_only = 3,
                min_damaging = 4, # Missense requires >= 3 out of 13 damaging predictions
                allowed_impacts = c("HIGH", "MODERATE"),
                require_cancer_gene = FALSE,
                remove_artifacts = TRUE,
                min_vaf_snv = 0.025,
                min_vaf_indel = 0.04,
                min_dp = 100,
                max_gnomad_af = 0.0001,
                max_1000g_af = 0.0001,
                min_unique_reads = 4,
                max_recurrence_fraction = 0.3,
                exclude_samples = c("AS9750"),
                low_cov_samples = NULL,
                low_cov_min_callers = 3,
                low_cov_min_damaging = 4,
                min_coverage_ffpe_matched = 22,
                min_coverage_ffpe_tumor_only = 40,
                min_coverage_fresh_matched = 30,
                min_coverage_fresh_tumor_only = 40,
                max_contamination = 0.05
            ),
            verbose = TRUE
        )
        fst::write.fst(
            filtered_variants_exome1,
            "projects/pscc/wes/output/data/filtered_variants_exome1.fst"
        )
    }
    cat(sprintf(
        "filtered_variants_exome1: %d rows, %.1f MB\n",
        nrow(filtered_variants_exome1),
        as.numeric(object.size(filtered_variants_exome1)) / 1024^2
    ))
}


cat(
  "  → Filtering variants2 for exome-wide oncoplot (no cancer gene restriction)...\n"
)
if (!exists("filtered_variants_exome2")) {
  if (
    file.exists("projects/pscc/wes/output/data/filtered_variants_exome2.fst")
  ) {
    cat("  → Loading exome-filtered variants from checkpoint...\n")
    message("Resuming: Loading filtered_variants_exome from checkpoint...")
    filtered_variants_exome2 <- data.table::setDT(fst::read.fst(
      "projects/pscc/wes/output/data/filtered_variants_exome2.fst"
    ))
  } else {
    if (!exists("merged_variants")) {
      merged_variants <- data.table::setDT(fst::read.fst(
        "projects/pscc/wes/output/data/merged_variants.fst"
      ))
    }
    cat(
      "  → Applying exome-wide variant filters (require_cancer_gene = FALSE)...\n"
    )
    message("Filtering variants for exome-wide oncoplot...")
    filtered_variants_exome2 <- filter_variants(
      variant_data = merged_variants,
      sample_col = sample_col_norm,
      metadata = metadata,
      coverage_dirs = dirs_coverage,
      contamination_dirs = dirs_contamination,
      export_blacklisted_path = NULL,
      export_low_coverage_path = NULL,
      export_coverage_path = NULL,
      filters = list(
        min_callers = 2,
        min_callers_high_impact = 2,
        min_callers_tumor_only = 2,
        min_callers_high_impact_tumor_only = 2,
        min_damaging = 7, 
        allowed_impacts = c("HIGH", "MODERATE"),
        required_caller_any = c("freebayes", "lofreq", "mutect2", "deepsomatic", "strelka"), 
        require_cosmic_id = TRUE,
        require_cancer_gene = FALSE,
        remove_artifacts = TRUE,
        min_dp = 100, #variants over that value pass
        min_vaf_snv = 0.025, # 2% VAF threshold for SNVs
        min_vaf_indel = 0.04, # 4% VAF threshold for indels
        max_gnomad_af = 0.0001, # 0.01% MAF threshold for gnomAD POPMAX
        max_1000g_af = 0.0001, # 0.01% MAF threshold for 1000G EUR
        min_unique_reads = 4, # Minimum unique (non-duplicate) reads - protects against PCR artifacts
        max_recurrence_fraction = 0.3, # Blacklist variants in >30% of samples (mapping artifacts)
        exclude_samples = c("AS9750"), # Add sample IDs to exclude if needed
        low_cov_samples = NULL, # Add low coverage sample IDs if needed
        low_cov_min_callers = 2,
        low_cov_min_damaging = 7,
        # Sample-type-specific coverage filtering
        min_coverage_ffpe_matched = 22, # FFPE matched samples
        min_coverage_ffpe_tumor_only = 40, # FFPE tumor-only
        min_coverage_fresh_matched = 30, # Fresh matched
        min_coverage_fresh_tumor_only = 40, # Fresh tumor-only
        # Contamination filtering (optional)
        max_contamination = 0.05# e.g., 0.05 to exclude samples with >5% contamination
      ),
      verbose = TRUE
    )
    fst::write.fst(
      filtered_variants_exome2,
      "projects/pscc/wes/output/data/filtered_variants_exome2.fst"
    )
  }
  cat(sprintf(
    "filtered_variants_exome2: %d rows, %.1f MB\n",
    nrow(filtered_variants_exome2),
    as.numeric(object.size(filtered_variants_exome2)) / 1024^2
  ))
}

# Now we merge filtered_variants_exome1 and filtered_variants_exome2 to see
# how many variants are in each list and how many overlap.
# This will help us understand the impact of the different filtering strategies
# on our final variant list.
cat("[8.3/10] Merging filtered variant exome lists...\n")
if (exists("filtered_variants_exome1") && exists("filtered_variants_exome2")) {
  # Avoid duplicated variants when binding the two filtered lists.
  duplicate_cols_exome <- intersect(names(filtered_variants_exome2), names(filtered_variants_exome1))
  if (length(duplicate_cols_exome) == 0) {
    stop("No shared columns found between filtered_variants_exome1 and filtered_variants_exome2.")
  }
  
  filtered_variants_exome2_only <- dplyr::anti_join(
    filtered_variants_exome2,
    filtered_variants_exome1,
    by = duplicate_cols_exome
  )
  
  filtered_variants_exome <- dplyr::bind_rows(
    filtered_variants_exome1,
    filtered_variants_exome2_only
  )
  
  setDT(filtered_variants_exome)
  filtered_variants_exome <- unique(filtered_variants_exome, by = duplicate_cols_exome)

  cat(sprintf("filtered_variants: %d rows\n", nrow(filtered_variants)))
  cat(sprintf("filtered_variants_exome: %d rows\n", nrow(filtered_variants_exome)))
  
  fst::write.fst(
    filtered_variants_exome,
    "projects/pscc/wes/output/data/filtered_variants_exome.fst"
  )
} else {
  stop("filtered_variants_exome1 and filtered_variants_exome2 must exist before combining.")
}
cat("Filtered exome variant lists combined\n\n")









cat("  → Converting exome variants to MAF format...\n")
dedup_cols_exome <- c(
    sample_col_norm,
    "number_chrom",
    "pos",
    "ref",
    "alt",
    "gene"
)
filtered_variants_exome_unique <- unique(
    filtered_variants_exome,
    by = dedup_cols_exome
)
message(sprintf(
    "Exome deduplicated variants: %d → %d (removed %d duplicate annotations)",
    nrow(filtered_variants_exome),
    nrow(filtered_variants_exome_unique),
    nrow(filtered_variants_exome) - nrow(filtered_variants_exome_unique)
))

clinical_cols_exome <- intersect(
    c(
        sample_col_norm,
        "tissue_type",
        "preservation",
        "facility",
        "patient",
        "matched"
    ),
    names(filtered_variants_exome_unique)
)
clinical_data_exome <- unique(filtered_variants_exome_unique[,
    ..clinical_cols_exome
])

maf_exome <- convert_to_maf(
    variant_data = filtered_variants_exome_unique,
    clinical_data = clinical_data_exome,
    sample_col = sample_col_norm,
    annotation_col = "annotation",
    create_maf_object = TRUE
)

rm(filtered_variants_exome_unique)
gc()

cat("  → Generating exome-wide oncoplot...\n")
n_genes_exome <- 30

gene_summary_exome <- maftools::getGeneSummary(maf_exome)
cat("\n=== Top 30 Exome Genes by MAF Summary ===\n")
print(gene_summary_exome[
    1:min(30, nrow(gene_summary_exome)),
    .(Hugo_Symbol, AlteredSamples, MutatedSamples, total)
])

top_genes_exome <- gene_summary_exome[1:n_genes_exome, Hugo_Symbol]

# VAF per gene for left bar
vaf_col_exome <- "i_VAF"
if (vaf_col_exome %in% maftools::getFields(maf_exome)) {
    gene_vaf_exome <- maftools::subsetMaf(
        maf = maf_exome,
        genes = top_genes_exome,
        fields = vaf_col_exome,
        mafObj = FALSE
    )[, .(mean(i_VAF, na.rm = TRUE)), by = Hugo_Symbol]
    colnames(gene_vaf_exome) <- c("gene", "VAF")
} else {
    gene_vaf_exome <- NULL
}

# Reuse annotation colors logic for exome clinical features
available_features_exome <- intersect(
    names(clinical_data_exome),
    c("tissue_type", "preservation", "facility", "matched")
)
annotation_colors_exome <- list()
palette_idx_exome <- 1
for (feat in available_features_exome) {
    unique_vals <- unique(as.character(clinical_data_exome[[feat]]))
    unique_vals <- unique_vals[!is.na(unique_vals)]
    if (length(unique_vals) > 0) {
        pal <- muted_palettes[[palette_idx_exome]]
        n_vals <- length(unique_vals)
        colors <- rep(pal, ceiling(n_vals / length(pal)))[1:n_vals]
        names(colors) <- unique_vals
        annotation_colors_exome[[feat]] <- colors
        palette_idx_exome <- (palette_idx_exome %% length(muted_palettes)) + 1
    }
}

png(
    filename = "projects/pscc/wes/output/plots/variants/oncoplot_exome.png",
    width = 24,
    height = 18,
    units = "in",
    res = 300
)
maftools::oncoplot(
    maf = maf_exome,
    top = n_genes_exome,
    draw_titv = FALSE,
    clinicalFeatures = available_features_exome,
    annotationColor = annotation_colors_exome,
    sortByAnnotation = FALSE,
    leftBarData = if (!is.null(gene_vaf_exome)) gene_vaf_exome else NULL,
    leftBarLims = if (!is.null(gene_vaf_exome)) c(0, 100) else NULL,
    showTumorSampleBarcodes = TRUE,
    removeNonMutated = TRUE,
    drawRowBar = TRUE,
    drawColBar = TRUE,
    fontSize = 0.6,
    SampleNamefontSize = 0.4,
    legendFontSize = 1.2,
    annotationFontSize = 1.2,
    gene_mar = 5,
    barcode_mar = 4
)
dev.off()

cat("✓ Exome-wide oncoplot generated\n\n")

cat("\n================================================\n")
cat("  ✓ Report generation complete!\n")
cat("================================================\n")
cat(sprintf("  Output directory: projects/pscc/wes/output/\n"))
cat(sprintf(
    "  Plots: %d generated\n",
    length(list.files("projects/pscc/wes/output/plots", pattern = "\\.png$"))
))
cat(sprintf(
    "  Data files: %d saved\n",
    length(list.files("projects/pscc/wes/output/data", pattern = "\\.fst$"))
))
cat("================================================\n\n")
