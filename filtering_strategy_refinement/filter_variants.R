#' Filter Variants Based on Multiple Criteria
#'
#' Filters somatic variants based on customizable thresholds including number
#' of callers, damaging predictions, variant impact, presence in cancer gene
#' databases, artifact removal predictions, coverage, and contamination. This
#' function is designed to work with merged variant data from multiple callers
#' where not all callers may be available for all samples (e.g., tumor-only vs
#' matched samples).
#'
#' @param variant_data A tibble or data frame containing merged variant data.
#'   Expected to have columns like: num_callers, num_damaging, impact,
#'   is_oncogene, is_tumor_suppressor_gene, is_driver, and caller-specific
#'   columns (e.g., freebayes_called, mutect2_called).
#' @param sample_col Name of the column containing sample IDs (default:
#'   "sample_id")
#' @param metadata Optional data frame containing metadata. If provided,
#'   `metadata_paths` is ignored. Column names should be lowercase.
#' @param metadata_paths Vector of paths to metadata files for
#'   coverage/contamination filtering (default: NULL). Required if using
#'   min_coverage or max_contamination filters and metadata is not provided.
#' @param coverage_dirs Vector of directories containing mosdepth summary files
#'   (default: NULL). Required if using min_coverage filter.
#' @param contamination_dirs Vector of directories containing contamination
#'   table files (default: NULL). Required if using max_contamination filter.
#' @param filters A named list of filtering criteria. Available filters:
#'   \itemize{
#'     \item \code{min_callers}: Minimum number of callers that must call the
#'       variant (default: 2)
#'     \item \code{min_callers_high_impact}: Minimum number of callers for HIGH
#'       impact variants (stop-gain, start-loss) and for all coding INDELs
#'       (HIGH/MODERATE/LOW) when \code{treat_coding_indels_as_high = TRUE}
#'       (default: 3)
#'     \item \code{min_callers_tumor_only}: Minimum number of callers for tumor-only
#'       samples (default: NULL, uses min_callers if not specified)
#'     \item \code{min_callers_high_impact_tumor_only}: Minimum number of callers
#'       for HIGH impact variants and coding INDELs in tumor-only samples
#'       (default: NULL, uses min_callers_high_impact if not specified)
#'     \item \code{min_damaging}: Minimum number of damaging predictions for
#'       MODERATE impact variants (default: 3)
#'     \item \code{allowed_impacts}: Character vector of allowed variant impacts
#'       (default: c("HIGH", "MODERATE")). LOW and MODIFIER SNVs are excluded
#'       by default. INDELs bypass this check entirely — coding INDELs are
#'       handled by the INDEL PATH; MODIFIER INDELs are dropped within that path.
#'     \item \code{require_cancer_gene}: Logical, if TRUE requires variant to be
#'       in a known cancer gene (default: TRUE). Uses: variant_data[cgc_gene == TRUE], 
#'       a column taken from COSMIC cancer gene census. Column created by SnpEff.
#'     \item \code{remove_artifacts}: Logical, if TRUE removes variants predicted
#'       as artifacts by SeqObjectiveBias (default: TRUE)
#'     \item \code{specific_callers}: Character vector of specific callers to
#'       require (default: NULL, meaning any callers)
#'     \item \code{required_caller_any}: Character vector of callers; variant
#'       must be called by at least one of them (OR logic). Default: NULL (no
#'       requirement). Example: c("freebayes", "lofreq", "mutect2",
#'       "deepsomatic", "strelka")
#'     \item \code{require_clinvar_pathogenic}: Logical; if TRUE, keeps only
#'       variants with clinvar_clnsig matching Pathogenic or Likely_pathogenic
#'       (case-insensitive, excludes "conflicting"). Default: FALSE
#'#'     \item \code{require_cosmic_id}: Logical; if TRUE, keeps only variants
#'       where the \code{id} column is not NA and not ".". Default: FALSE 
#'     \item \code{require_cosmic_cgc}: Logical; if TRUE, keeps only variants
#'       where cosmic_cgc_v99 == "Yes". Gene Database of Cosmic Cancer Gene Census. Default: FALSE
#'     \item \code{exclude_samples}: Character vector of sample IDs to exclude
#'       (default: NULL)
#'     \item \code{min_vaf_snv}: Minimum variant allele frequency for SNVs
#'       (default: 0.02, i.e., 2%)
#'     \item \code{min_vaf_indel}: Minimum variant allele frequency for indels
#'       (default: 0.04, i.e., 4%)
#'     \item \code{max_gnomad_af}: Maximum gnomAD allele frequency for germline
#'       filtering (default: 0.0001, i.e., 0.01%)
#'     \item \code{max_1000g_af}: Maximum 1000 Genomes allele frequency for
#'       germline filtering (default: 0.0001, i.e., 0.01%)
#'     \item \code{min_dp}: Minimum total read depth (DP) at the variant
#'       position. Variants with \code{dp < min_dp} are removed regardless of
#'       any other criteria. (default: NULL, disabled)
#'     \item \code{min_unique_reads}: Minimum unique (non-duplicate) read depth
#'       at variant position (default: NULL). Requires UD field from bam-readcount
#'       annotation in variant_calling_par.sh. Set to e.g. 4 to filter low-support
#'       variants that may be artifacts from PCR duplicates.
#'     \item \code{treat_coding_indels_as_high}: Logical, if TRUE all coding INDELs
#'       (impact HIGH, MODERATE, or LOW) are filtered by the HIGH-impact caller
#'       threshold (\code{min_callers_high_impact}), regardless of SnpEff impact
#'       level. MODIFIER INDELs (intronic, intergenic, UTR) are still excluded.
#'       Requires \code{is_indel} column from \code{merge_variants()}. Follows
#'       best-practice cancer WES pipelines. Default: TRUE.
#'     \item \code{min_coverage}: Global minimum mean coverage threshold (default: NULL)
#'     \item \code{min_coverage_ffpe_matched}: Minimum coverage for FFPE matched samples
#'       (default: NULL, uses min_coverage if not specified)
#'     \item \code{min_coverage_ffpe_tumor_only}: Minimum coverage for FFPE tumor-only
#'       samples (default: NULL, uses min_coverage if not specified)
#'     \item \code{min_coverage_fresh_matched}: Minimum coverage for fresh matched samples
#'       (default: NULL, uses min_coverage if not specified)
#'     \item \code{min_coverage_fresh_tumor_only}: Minimum coverage for fresh tumor-only
#'       samples (default: NULL, uses min_coverage if not specified)
#'     \item \code{max_contamination}: Maximum contamination threshold (default:
#'       NULL)
#'     \item \code{low_cov_samples}: Character vector of low coverage sample IDs
#'       that require stricter filtering (default: NULL)
#'     \item \code{low_cov_min_callers}: Minimum callers for low coverage samples
#'       (default: 4)
#'     \item \code{low_cov_min_damaging}: Minimum damaging predictions for low
#'       coverage samples (default: 6)
#'   }
#' @param export_blacklisted_path Optional file path to export blacklisted
#'   recurrent variants before removal (default: NULL). When provided alongside
#'   \code{max_recurrence_fraction}, writes a TSV with per-position summaries
#'   including gene, annotation, sample count, and per-caller call percentages.
#'   Useful for documenting which variants were removed and investigating
#'   caller-specific patterns in systematic artifacts.
#' @param export_low_coverage_path Optional file path to export mutational load
#'   statistics for samples filtered due to insufficient coverage (default: NULL).
#'   When provided, writes a TSV with sample ID, coverage, variant count, and
#'   sample type. Useful for identifying samples with high mutational load that
#'   may benefit from resequencing at higher coverage.
#' @param export_passed_coverage_path Optional file path to export coverage
#'   statistics for samples that pass the coverage filter (default: NULL).
#'   When provided, writes a TSV with sample ID, total variant count,
#'   cancer gene variant count, mean coverage, threshold, and sample type.
#'   For matched samples, also includes the paired normal sample ID and its
#'   mean coverage. Useful for QC review of samples included in the analysis.
#' @param export_coverage_path Optional file path to export coverage statistics
#'  for all samples (default: NULL). When provided, writes a TSV with sample ID, 
#' mean coverage, threshold, and sample type for all samples in the coverage data. 
#' Useful for comprehensive documentation of coverage filtering and for identifying 
#' patterns across all samples.
#' @param export_contamination_path Optional file path to export contamination statistics
#'  for all samples (default: NULL). When provided, writes a TSV with sample ID, contamination
#'  level, and sample type for all samples in the contamination data. Useful for comprehensive 
#' documentation of contamination filtering and for identifying patterns across all samples.
#' @param verbose Logical indicating whether to print filtering statistics
#'   (default: TRUE)
#'
#' @return A filtered tibble/data frame containing only variants that pass all
#'   specified criteria.
#'
#' @details
#' The function applies filters in the following order:
#' \enumerate{
#'   \item Filters by coverage if \code{min_coverage} is specified
#'   \item Filters by contamination if \code{max_contamination} is specified
#'   \item Excludes specified samples if \code{exclude_samples} is provided
#'   \item Applies stricter filters to low coverage samples
#'   \item Filters by variant allele frequency (VAF)
#'   \item Filters by germline polymorphism frequency (gnomAD/1000G)
#'   \item Filters by total mininum depth \code{min_dp} if specified    
#'   \item Filters by unique read depth (if specified)
#'   \item Filters by minimum number of callers
#'   \item Filters by specific callers (if specified)
#'   \item Filters requiring at least one caller from required_callers_any (if specified)
#'   \item Filters requiring ClinVar Pathogenic or Likely_Pathogenic (if specified)
#'   \item Filters requiring COSMIC_id if specified
#'   \item Filters requiring COSMIC CGC membership \code{cosmic_cgc_v99} (if specified)
#'   \item Filters by variant impact and damaging predictions
#'   \item Blacklists recurrent variants exceeding \code{max_recurrence_fraction}
#'   \item Filters by cancer gene presence (if enabled)
#'   \item Removes predicted artifacts (if enabled)
#' }
#'
#' For HIGH impact variants (frameshifts, stop-gain/loss, INDELs), the function
#' uses caller-based filtering with \code{min_callers_high_impact}.
#' For MODERATE impact variants (missense), the function requires at least
#' \code{min_damaging} predictions if dbNSFP data is available, otherwise falls
#' back to caller-based filtering.
#'
#' Coverage and contamination filtering require the helper functions
#' \code{read_mosdepth_summaries} and \code{read_contamination_tables} to be
#' available in the environment (e.g., sourced from coverage.R and
#' contamination.R).
#'
#' @examples
#' \dontrun{
#' # Basic filtering with defaults
#' filtered_vars <- filter_variants(merged_variants)
#'
#' # Custom filtering with coverage and contamination
#' filtered_vars <- filter_variants(
#'   variant_data = merged_variants,
#'   sample_col = "Sample_ID",
#'   metadata_paths = c("config/metadata.csv"),
#'   coverage_dirs = c("analysis/mosdepth"),
#'   contamination_dirs = c("analysis/contamination"),
#'   filters = list(
#'     min_callers = 3,
#'     min_damaging = 5,
#'     allowed_impacts = c("HIGH", "MODERATE"),
#'     require_cancer_gene = TRUE,
#'     min_coverage = 60,
#'     max_contamination = 0.05,
#'     exclude_samples = c("SAMPLE1", "SAMPLE2"),
#'     low_cov_samples = c("SAMPLE3", "SAMPLE4"),
#'     low_cov_min_callers = 4
#'   )
#' )
#' }
#'
#' @export
filter_variants <- function(
    variant_data,
    sample_col = "sample_id",
    metadata = NULL,
    metadata_paths = NULL,
    coverage_dirs = NULL,
    contamination_dirs = NULL,
    filters = list(),
    export_blacklisted_path = NULL,
    export_low_coverage_path = NULL,
    export_passed_coverage_path = NULL,
    export_coverage_path = NULL, 
    export_contamination_path = NULL,
    verbose = TRUE
) {
    # Convert to data.table if not already (modifies by reference)
    if (!data.table::is.data.table(variant_data)) {
        data.table::setDT(variant_data)
    }

    # Store initial count
    n_initial <- nrow(variant_data)
    if (verbose) {
        message(sprintf("Starting with %d variants", n_initial))
    }

    # Set default filter values
    default_filters <- list(
        min_callers = 2,
        min_callers_high_impact = 3,
        min_callers_tumor_only = NULL, # Add +1 to min_callers for tumor-only if specified
        min_callers_high_impact_tumor_only = NULL, # Add +1 for high impact tumor-only
        min_damaging = 3,
        allowed_impacts = c("HIGH", "MODERATE"),
        require_cancer_gene = TRUE,
        remove_artifacts = TRUE,
        specific_callers = NULL,
        required_caller_any = NULL,
        require_clinvar_pathogenic = FALSE,
        require_cosmic_id = FALSE,
        require_cosmic_cgc = FALSE,
        exclude_samples = NULL,
        min_coverage = NULL,
        min_coverage_ffpe_matched = NULL,
        min_coverage_ffpe_tumor_only = NULL,
        min_coverage_fresh_matched = NULL,
        min_coverage_fresh_tumor_only = NULL,
        max_contamination = NULL,
        min_vaf_snv = 0.02,
        min_vaf_indel = 0.04,
        max_gnomad_af = NULL,
        max_1000g_af = NULL,
        min_dp = NULL,
        min_unique_reads = NULL,
        treat_coding_indels_as_high = TRUE,
        low_cov_samples = NULL,
        low_cov_min_callers = 4,
        low_cov_min_damaging = 6,
        max_recurrence_fraction = NULL
    )

    # Merge user filters with defaults
    filters <- modifyList(default_filters, filters)

    # Normalize sample_col to lowercase (metadata uses Sample_ID but after merge_with_metadata all columns are lowercase)
    variant_sample_col <- tolower(sample_col)

    # Validate sample column exists
    if (!variant_sample_col %in% names(variant_data)) {
        stop(sprintf(
            "Sample column '%s' not found in variant_data. Available columns: %s",
            variant_sample_col,
            paste(head(names(variant_data), 20), collapse = ", ")
        ))
    }

    # Diagnostic: Check for variants with matched==NA (shouldn't happen for tumor variants)
    if ("matched" %in% names(variant_data)) {
        n_na_matched <- sum(is.na(variant_data[["matched"]]))
        if (n_na_matched > 0 && verbose) {
            warning(sprintf(
                "%d variants have matched==NA. These represent normal samples and will use default thresholds. Consider excluding normal samples from variant data.",
                n_na_matched
            ))
        }
    }

    # Step 1: Filter by coverage if specified
    # Support both global min_coverage and sample-type-specific thresholds
    has_coverage_filter <- !is.null(filters$min_coverage) ||
        !is.null(filters$min_coverage_ffpe_matched) ||
        !is.null(filters$min_coverage_ffpe_tumor_only) ||
        !is.null(filters$min_coverage_fresh_matched) ||
        !is.null(filters$min_coverage_fresh_tumor_only)

    if (has_coverage_filter) {
        if (
            (is.null(metadata) && is.null(metadata_paths)) ||
                is.null(coverage_dirs)
        ) {
            warning(
                "min_coverage specified but metadata/metadata_paths or coverage_dirs missing. Skipping coverage filter."
            )
        } else {
            n_before <- nrow(variant_data)

            # Load coverage statistics (already merged with metadata by read_mosdepth_summaries)
            tryCatch(
                {
                    coverage_stats <- read_mosdepth_summaries(
                        metadata = metadata,
                        metadata_paths = metadata_paths,
                        dirs = coverage_dirs,
                        sample_col = tolower(sample_col)
                    )
                    data.table::setDT(coverage_stats)

                    # Validate required columns exist for sample-type-specific filtering
                    if (
                        !is.null(filters$min_coverage_ffpe_matched) ||
                            !is.null(filters$min_coverage_ffpe_tumor_only) ||
                            !is.null(filters$min_coverage_fresh_matched) ||
                            !is.null(filters$min_coverage_fresh_tumor_only)
                    ) {
                        required_cols <- c("preservation", "matched")
                        missing_cols <- setdiff(
                            required_cols,
                            names(coverage_stats)
                        )
                        if (length(missing_cols) > 0) {
                            stop(sprintf(
                                "Required metadata columns missing for sample-type-specific coverage filtering: %s. Available columns: %s",
                                paste(missing_cols, collapse = ", "),
                                paste(names(coverage_stats), collapse = ", ")
                            ))
                        }
                    }

                    # Determine threshold for each sample based on type
                    # Initialize threshold column if using global min_coverage
                    if (!is.null(filters$min_coverage)) {
                        coverage_stats[, threshold := filters$min_coverage]
                    } else {
                        coverage_stats[, threshold := NA_real_]
                    }

                    # Apply sample-type-specific thresholds if provided
                    # Use case-insensitive comparisons and handle NA values
                    if (!is.null(filters$min_coverage_ffpe_matched)) {
                        coverage_stats[
                            !is.na(preservation) &
                                !is.na(matched) &
                                tolower(preservation) == "ffpe" &
                                tolower(matched) == "yes",
                            threshold := filters$min_coverage_ffpe_matched
                        ]
                    }
                    if (!is.null(filters$min_coverage_ffpe_tumor_only)) {
                        coverage_stats[
                            !is.na(preservation) &
                                !is.na(matched) &
                                tolower(preservation) == "ffpe" &
                                tolower(matched) == "no",
                            threshold := filters$min_coverage_ffpe_tumor_only
                        ]
                    }
                    if (!is.null(filters$min_coverage_fresh_matched)) {
                        coverage_stats[
                            !is.na(preservation) &
                                !is.na(matched) &
                                tolower(preservation) == "fresh" &
                                tolower(matched) == "yes",
                            threshold := filters$min_coverage_fresh_matched
                        ]
                    }
                    if (!is.null(filters$min_coverage_fresh_tumor_only)) {
                        coverage_stats[
                            !is.na(preservation) &
                                !is.na(matched) &
                                tolower(preservation) == "fresh" &
                                tolower(matched) == "no",
                            threshold := filters$min_coverage_fresh_tumor_only
                        ]
                    }

                    # Identify samples below their respective threshold
                    low_cov_samples_detected <- coverage_stats[
                        !is.na(coverage_stats$threshold) &
                            coverage_stats$mean < coverage_stats$threshold,
                        sample
                    ]

                    # Count normal samples (matched==NA) for informational message
                    n_normals <- sum(
                        is.na(coverage_stats$matched),
                        na.rm = TRUE
                    )
                    if (verbose && n_normals > 0) {
                        message(sprintf(
                            "  Note: %d normal samples (matched==NA) in coverage data are not filtered (they don't produce variants)",
                            n_normals
                        ))
                    }

                    # Warn about samples without threshold assignment (missing metadata)
                    # Exclude normal samples (matched==NA) from this warning since they don't have variants
                    samples_without_threshold <- coverage_stats[
                        is.na(coverage_stats$threshold) &
                            !is.na(coverage_stats$sample) &
                            !is.na(coverage_stats$matched),
                        sample
                    ]
                    if (length(samples_without_threshold) > 0 && verbose) {
                        warning(sprintf(
                            "Coverage filtering: %d tumor samples have missing preservation/matched metadata and were not filtered (kept regardless of coverage). Set min_coverage to apply a global threshold.",
                            length(samples_without_threshold)
                        ))
                    }

                    if (length(low_cov_samples_detected) > 0) {
                        # Calculate mutational load for low coverage samples before removal
                        # This helps identify samples that might benefit from resequencing
                        low_cov_variant_counts <- variant_data[
                            variant_data[[variant_sample_col]] %in%
                                low_cov_samples_detected,
                            .(variant_count = .N),
                            by = variant_sample_col
                        ]

                        # Merge with coverage stats to get coverage and metadata
                        low_cov_summary <- coverage_stats[
                            sample %in% low_cov_samples_detected
                        ]
                        data.table::setnames(
                            low_cov_variant_counts,
                            variant_sample_col,
                            "sample"
                        )
                        low_cov_summary <- low_cov_summary[
                            low_cov_variant_counts,
                            on = "sample"
                        ]

                        # Sort by variant count (descending) to prioritize resequencing
                        data.table::setorder(low_cov_summary, -variant_count)

                        # Export if path provided
                        if (!is.null(export_low_coverage_path)) {
                            # Select relevant columns for export
                            export_cols <- intersect(
                                c(
                                    "sample",
                                    "variant_count",
                                    "mean",
                                    "threshold",
                                    "preservation",
                                    "matched",
                                    "tissue_type",
                                    "facility"
                                ),
                                names(low_cov_summary)
                            )
                            low_cov_export <- low_cov_summary[, ..export_cols]

                            # Ensure output directory exists
                            export_dir <- dirname(export_low_coverage_path)
                            if (!dir.exists(export_dir)) {
                                dir.create(
                                    export_dir,
                                    recursive = TRUE,
                                    showWarnings = FALSE
                                )
                            }

                            data.table::fwrite(
                                low_cov_export,
                                file = export_low_coverage_path,
                                sep = "\t",
                                quote = FALSE
                            )

                            if (verbose) {
                                message(sprintf(
                                    "  \u2713 Exported mutational load for %d low-coverage samples to: %s",
                                    nrow(low_cov_export),
                                    export_low_coverage_path
                                ))
                                # Show top candidates for resequencing
                                top_n <- min(5, nrow(low_cov_summary))
                                message(sprintf(
                                    "  Top %d candidates for resequencing (by variant count):",
                                    top_n
                                ))
                                for (i in seq_len(top_n)) {
                                    message(sprintf(
                                        "    %s: %d variants, %.1fx coverage (threshold: %.1fx)",
                                        low_cov_summary$sample[i],
                                        low_cov_summary$variant_count[i],
                                        low_cov_summary$mean[i],
                                        low_cov_summary$threshold[i]
                                    ))
                                }
                            }
                        }

                        # Filter by reference - more efficient
                        variant_data <- variant_data[
                            !(variant_data[[variant_sample_col]] %in%
                                low_cov_samples_detected)
                        ]

                        n_after <- nrow(variant_data)
                        if (verbose) {
                            message(sprintf(
                                "Excluded %d samples with insufficient coverage: removed %s variants from those samples (%s remaining). Samples removed: %s",
                                length(low_cov_samples_detected),
                                format(n_before - n_after, big.mark = ","),
                                format(n_after, big.mark = ","),
                                ifelse(length(low_cov_samples_detected) > 0,
                                        paste(low_cov_samples_detected, collapse = ", "),
                                        "None")
                            ))
                            # Show breakdown by sample type if verbose
                            removed_summary <- coverage_stats[
                                sample %in% low_cov_samples_detected,
                                .(
                                    n_samples = .N,
                                    mean_cov = mean(mean, na.rm = TRUE),
                                    threshold = first(threshold)
                                ),
                                by = .(preservation, matched)
                            ]
                            if (nrow(removed_summary) > 0) {
                                message("  Breakdown by sample type:")
                                for (i in seq_len(nrow(removed_summary))) {
                                    message(sprintf(
                                        "    %s/%s: %d samples (mean=%.1fx, threshold=%.1fx)",
                                        removed_summary$preservation[i],
                                        removed_summary$matched[i],
                                        removed_summary$n_samples[i],
                                        removed_summary$mean_cov[i],
                                        removed_summary$threshold[i]
                                    ))
                                }
                            }
                        }
                    } else {
                        if (verbose) {
                            message(
                                "All samples meet minimum coverage thresholds"
                            )
                        }
                    }

                    # Export all coverage data if requested
                    if (!is.null(export_coverage_path)) {
                        export_dir_coverage <- dirname(export_coverage_path)
                        if (!dir.exists(export_dir_coverage)) {
                            dir.create(
                                export_dir_coverage,
                                recursive = TRUE,
                                showWarnings = FALSE
                            )
                        }

                        data.table::fwrite(
                            coverage_stats,
                            file = export_coverage_path,
                            sep = "\t",
                            quote = FALSE
                        )

                        if (verbose) {
                            message(sprintf(
                                "  \u2713 Exported coverage data for %d samples to: %s",
                                nrow(coverage_stats),
                                export_coverage_path
                            ))
                        }
                    }

                    # Export passed coverage samples if requested
                    if (!is.null(export_passed_coverage_path)) {
                        # Tumor samples that passed coverage (have a threshold and are above it)
                        passed_tumor_samples <- coverage_stats[
                            !is.na(coverage_stats$threshold) &
                                !is.na(coverage_stats$matched) &
                                !(coverage_stats$sample %in% low_cov_samples_detected),
                            sample
                        ]

                        if (length(passed_tumor_samples) > 0) {
                            # Total variant counts per passing sample
                            passed_variant_counts_total <- variant_data[
                                variant_data[[variant_sample_col]] %in%
                                    passed_tumor_samples,
                                .(variant_count_total = .N),
                                by = variant_sample_col
                            ]
                            data.table::setnames(
                                passed_variant_counts_total,
                                variant_sample_col,
                                "sample"
                            )

                            # Cancer gene variant counts (is_oncogene, is_tumor_suppressor_gene, or is_driver)
                            cancer_cols <- intersect(
                                c("is_oncogene", "is_tumor_suppressor_gene", "is_driver"),
                                names(variant_data)
                            )
                            if (length(cancer_cols) > 0) {
                                cancer_mask <- Reduce(
                                    `|`,
                                    lapply(cancer_cols, function(col) {
                                        v <- variant_data[[col]]
                                        if (is.logical(v)) {
                                            !is.na(v) & v
                                        } else {
                                            !is.na(v) & tolower(as.character(v)) %in% c("true", "yes", "1")
                                        }
                                    })
                                )
                                passed_variant_counts_cancer <- variant_data[
                                    variant_data[[variant_sample_col]] %in%
                                        passed_tumor_samples & cancer_mask,
                                    .(variant_count_cancer = .N),
                                    by = variant_sample_col
                                ]
                                data.table::setnames(
                                    passed_variant_counts_cancer,
                                    variant_sample_col,
                                    "sample"
                                )
                            } else {
                                passed_variant_counts_cancer <- data.table::data.table(
                                    sample = character(0),
                                    variant_count_cancer = integer(0)
                                )
                            }

                            # Build passed summary from coverage_stats (tumor samples only)
                            passed_summary <- coverage_stats[
                                sample %in% passed_tumor_samples
                            ]

                            # Merge variant counts
                            passed_summary <- merge(
                                passed_summary,
                                passed_variant_counts_total,
                                by = "sample",
                                all.x = TRUE
                            )
                            passed_summary <- merge(
                                passed_summary,
                                passed_variant_counts_cancer,
                                by = "sample",
                                all.x = TRUE
                            )
                            # Fill NAs with 0 for samples with no variants
                            passed_summary[
                                is.na(variant_count_total),
                                variant_count_total := 0L
                            ]
                            passed_summary[
                                is.na(variant_count_cancer),
                                variant_count_cancer := 0L
                            ]

                            # Add normal sample coverage for matched samples
                            # Normal samples have matched == NA in coverage_stats
                            if ("patient" %in% names(coverage_stats)) {
                                normal_cov <- coverage_stats[
                                    is.na(coverage_stats$matched),
                                    .(normal_sample = sample, normal_mean_coverage = mean, patient)
                                ]
                                if (nrow(normal_cov) > 0 && "patient" %in% names(passed_summary)) {
                                    passed_summary <- merge(
                                        passed_summary,
                                        normal_cov,
                                        by = "patient",
                                        all.x = TRUE
                                    )
                                }
                            }

                            # Sort by mean coverage
                            data.table::setorder(passed_summary, mean)

                            # Select relevant columns for export
                            export_cols_passed <- intersect(
                                c(
                                    "sample",
                                    "variant_count_total",
                                    "variant_count_cancer",
                                    "mean",
                                    "threshold",
                                    "preservation",
                                    "matched",
                                    "tissue_type",
                                    "facility",
                                    "patient",
                                    "normal_sample",
                                    "normal_mean_coverage"
                                ),
                                names(passed_summary)
                            )
                            passed_export <- passed_summary[, ..export_cols_passed]

                            export_dir_passed <- dirname(export_passed_coverage_path)
                            if (!dir.exists(export_dir_passed)) {
                                dir.create(
                                    export_dir_passed,
                                    recursive = TRUE,
                                    showWarnings = FALSE
                                )
                            }

                            data.table::fwrite(
                                passed_export,
                                file = export_passed_coverage_path,
                                sep = "\t",
                                quote = FALSE
                            )

                            if (verbose) {
                                message(sprintf(
                                    "  \u2713 Exported coverage stats for %d passing samples to: %s",
                                    nrow(passed_export),
                                    export_passed_coverage_path
                                ))
                            }
                        }
                    }
                },
                error = function(e) {
                    warning(sprintf(
                        "Error reading coverage data: %s. Skipping coverage filter.",
                        e$message
                    ))
                }
            )
        }
    }

    # Step 2: Filter by contamination if specified
    if (!is.null(filters$max_contamination)) {
        if (
            (is.null(metadata) && is.null(metadata_paths)) ||
                is.null(contamination_dirs)
        ) {
            warning(
                "max_contamination specified but metadata/metadata_paths or contamination_dirs missing. Skipping contamination filter."
            )
        } else {
            n_before <- nrow(variant_data)

            if (verbose) {
                message(sprintf(
                    "Filtering by maximum contamination %.2f%%...",
                    filters$max_contamination * 100
                ))
            }

            # Load contamination statistics using existing function
            tryCatch(
                {
                    contamination_stats <- read_contamination_tables(
                        metadata = metadata,
                        metadata_paths = metadata_paths,
                        dirs = contamination_dirs,
                        sample_col = tolower(sample_col)
                    )
                    data.table::setDT(contamination_stats)

                    # Identify samples above contamination threshold
                    high_contam_samples <- contamination_stats[
                        contamination_stats$contamination >
                            filters$max_contamination,
                        sample
                    ]

                    if (length(high_contam_samples) > 0) {
                        # Filter by reference - more efficient
                        variant_data <- variant_data[
                            !(variant_data[[variant_sample_col]] %in%
                                high_contam_samples)
                        ]

                        n_after <- nrow(variant_data)
                        if (verbose) {
                            message(sprintf(
                                "Excluded %d samples with >%.2f%% contamination: removed %s variants from those samples (%s remaining). Samples: %s",
                                length(high_contam_samples),
                                filters$max_contamination * 100,
                                format(n_before - n_after, big.mark = ","),
                                format(n_after, big.mark = ","),
                                ifelse(length(high_contam_samples) > 0,
                                        paste(high_contam_samples, collapse = ", "),
                                        "None")
                            ))
                        }
                    } else {
                        if (verbose) {
                            message(
                                "All samples meet maximum contamination threshold"
                            )
                        }
                    }
                
                # Export all contamination data if requested
                    if (!is.null(export_contamination_path)) {
                        export_dir_contamination <- dirname(export_contamination_path)
                        if (!dir.exists(export_dir_contamination)) {
                            dir.create(
                                export_dir_contamination,
                                recursive = TRUE,
                                showWarnings = FALSE
                            )
                        }

                        data.table::fwrite(
                            contamination_stats,
                            file = export_contamination_path,
                            sep = "\t",
                            quote = FALSE
                        )

                        if (verbose) {
                            message(sprintf(
                                "  \u2713 Exported contamination data for %d samples to: %s",
                                nrow(contamination_stats),
                                export_contamination_path
                            ))
                        }
                    }
                },
                error = function(e) {
                    warning(sprintf(
                        "Error reading contamination data: %s. Skipping contamination filter.",
                        e$message
                    ))
                }
            )
        }
    }

    # Step 3: Exclude specified samples
    if (
        !is.null(filters$exclude_samples) && length(filters$exclude_samples) > 0
    ) {
        n_before <- nrow(variant_data)

        variant_data <- variant_data[
            !(variant_data[[variant_sample_col]] %in% filters$exclude_samples)
        ]

        n_after <- nrow(variant_data)
        if (verbose) {
            message(sprintf(
                "Excluded %d samples: removed %s variants from those samples (%s remaining). Samples: %s",
                length(filters$exclude_samples),
                format(n_before - n_after, big.mark = ","),
                format(n_after, big.mark = ","),
                paste(filters$exclude_samples, collapse = ", "),
                ifelse(length(filters$exclude_samples) > 0,
                        paste(filters$exclude_samples, collapse = ", "),
                        "None")
            ))
        }
    }

    # Step 4: Apply stricter filters for low coverage samples
    if (
        !is.null(filters$low_cov_samples) && length(filters$low_cov_samples) > 0
    ) {
        n_before <- nrow(variant_data)

        # Check if required columns exist
        if (
            !"num_callers" %in% names(variant_data) ||
                !"num_damaging" %in% names(variant_data)
        ) {
            warning(
                "Columns 'num_callers' or 'num_damaging' not found. Skipping low coverage sample filter."
            )
        } else {
            # Create mask for low coverage samples with strict criteria
            mask_low_cov <- variant_data[[variant_sample_col]] %in%
                filters$low_cov_samples
            mask_strict <- mask_low_cov &
                (variant_data$num_callers < filters$low_cov_min_callers) &
                (variant_data$num_damaging < filters$low_cov_min_damaging)

            # Remove variants that don't meet strict criteria
            variant_data <- variant_data[!mask_strict, ]

            n_after <- nrow(variant_data)
            if (verbose) {
                message(sprintf(
                    "Applied strict filters to %d low coverage samples: removed %d variants (%d remaining)",
                    length(filters$low_cov_samples),
                    n_before - n_after,
                    n_after
                ))
            }
        }
    }

    # Step 5: Filter by variant allele frequency (VAF)
    # SNVs: >= 2%, Indels: >= 4%
    if (!is.null(filters$min_vaf_snv) || !is.null(filters$min_vaf_indel)) {
        n_before <- nrow(variant_data)

        if ("af" %in% names(variant_data)) {
            # Determine variant type (SNV vs indel)
            # SNV: REF and ALT are both single nucleotide
            # Indel: REF or ALT length != 1
            if (
                "ref" %in% names(variant_data) && "alt" %in% names(variant_data)
            ) {
                # Direct filtering without intermediate columns - more memory efficient
                is_snv <- nchar(variant_data$ref) == 1 &
                    nchar(variant_data$alt) == 1
                is_indel <- !is_snv

                # Build combined filter condition
                # Explicitly handle NA AF: keep variants with NA AF
                # (they should not be silently removed by the comparison)
                af_is_na <- is.na(variant_data$af)
                if (
                    !is.null(filters$min_vaf_snv) &&
                        !is.null(filters$min_vaf_indel)
                ) {
                    keep_rows <- af_is_na |
                        (is_snv & variant_data$af >= filters$min_vaf_snv) |
                        (is_indel & variant_data$af >= filters$min_vaf_indel)
                    variant_data <- variant_data[keep_rows, ]
                } else if (!is.null(filters$min_vaf_snv)) {
                    keep_rows <- af_is_na |
                        !is_snv |
                        variant_data$af >= filters$min_vaf_snv
                    variant_data <- variant_data[keep_rows, ]
                } else if (!is.null(filters$min_vaf_indel)) {
                    keep_rows <- af_is_na |
                        !is_indel |
                        variant_data$af >= filters$min_vaf_indel
                    variant_data <- variant_data[keep_rows, ]
                }

                n_after <- nrow(variant_data)
                if (verbose) {
                    message(sprintf(
                        "Filtered by VAF (SNV >= %.1f%%, Indel >= %.1f%%): removed %d variants (%d remaining)",
                        filters$min_vaf_snv * 100,
                        filters$min_vaf_indel * 100,
                        n_before - n_after,
                        n_after
                    ))
                }
            } else {
                warning("REF or ALT columns not found, skipping VAF filter")
            }
        } else {
            warning("AF column not found, skipping VAF filter")
        }
    }

    # Step 6: Filter by germline polymorphism frequency
    # Remove variants with MAF >= 0.01% in 1000G EUR or gnomAD POPMAX
    if (!is.null(filters$max_gnomad_af) || !is.null(filters$max_1000g_af)) {
        n_before <- nrow(variant_data)
        n_gnomad_removed <- 0
        n_1000g_removed <- 0

        # Check for gnomAD population frequency columns
        # Common column names from dbNSFP: gnomad_exomes_af, gnomad_genomes_af,
        # gnomad_exomes_af_popmax, gnomad_genomes_af_popmax
        # Also match snake_case variants: gnom_ad (from camelCase cleaning of gnomAD)
        gnomad_cols <- grep(
            "gnomad.*af|gnomad.*popmax|gnom_ad.*af|gnom_ad.*popmax",
            names(variant_data),
            value = TRUE,
            ignore.case = TRUE
        )

        if (length(gnomad_cols) > 0 && !is.null(filters$max_gnomad_af)) {
            # Convert all columns to numeric at once if needed
            for (col in gnomad_cols) {
                if (!is.numeric(variant_data[[col]])) {
                    data.table::set(
                        variant_data,
                        j = col,
                        value = suppressWarnings(as.numeric(variant_data[[
                            col
                        ]]))
                    )
                }
            }

            # Filter once using all columns - more efficient
            n_before_col <- nrow(variant_data)
            keep_rows <- rep(TRUE, nrow(variant_data))
            for (col in gnomad_cols) {
                keep_rows <- keep_rows &
                    (is.na(variant_data[[col]]) |
                        variant_data[[col]] < filters$max_gnomad_af)
            }
            variant_data <- variant_data[keep_rows, ]
            n_gnomad_removed <- n_before_col - nrow(variant_data)

            if (verbose && n_gnomad_removed > 0) {
                message(sprintf(
                    "gnomAD filter: removed %d variants with MAF >= %.4f%%",
                    n_gnomad_removed,
                    filters$max_gnomad_af * 100
                ))
            }
        }

        # Check for 1000 Genomes population frequency columns
        # Common column names: 1000gp3_af, 1000gp3_eur_af, af_1kg, etc.
        # Exclude gnomAD and RegeneronME columns that also contain "eur" and "af"
        kg_cols <- grep(
            "1000g.*af|1000gp.*af|af.*1kg",
            names(variant_data),
            value = TRUE,
            ignore.case = TRUE
        )
        # Remove any gnomAD or RegeneronME columns that may have matched
        kg_cols <- kg_cols[
            !grepl("gnomad|gnom_ad|regeneron", kg_cols, ignore.case = TRUE)
        ]

        if (length(kg_cols) > 0 && !is.null(filters$max_1000g_af)) {
            # Convert all columns to numeric at once if needed
            for (col in kg_cols) {
                if (!is.numeric(variant_data[[col]])) {
                    data.table::set(
                        variant_data,
                        j = col,
                        value = suppressWarnings(as.numeric(variant_data[[
                            col
                        ]]))
                    )
                }
            }

            # Filter once using all columns - more efficient
            n_before_col <- nrow(variant_data)
            keep_rows <- rep(TRUE, nrow(variant_data))
            for (col in kg_cols) {
                keep_rows <- keep_rows &
                    (is.na(variant_data[[col]]) |
                        variant_data[[col]] < filters$max_1000g_af)
            }
            variant_data <- variant_data[keep_rows, ]
            n_1000g_removed <- n_before_col - nrow(variant_data)

            if (verbose && n_1000g_removed > 0) {
                message(sprintf(
                    "1000 Genomes filter: removed %d variants with MAF >= %.4f%%",
                    n_1000g_removed,
                    filters$max_1000g_af * 100
                ))
            }
        }

        n_after <- nrow(variant_data)
        n_total_removed <- n_before - n_after

        if (verbose && n_total_removed > 0) {
            message(sprintf(
                "Total germline polymorphisms removed: %d variants (%d remaining)",
                n_total_removed,
                n_after
            ))
        } else if (verbose && n_total_removed == 0) {
            if (length(gnomad_cols) == 0 && length(kg_cols) == 0) {
                message(
                    "No population frequency columns found, skipping germline filter"
                )
            } else {
                message("No common germline polymorphisms detected")
            }
        }
    }

    # Step 6.4: Filter by total read depth (DP)
    if (!is.null(filters$min_dp)) {
        n_before <- nrow(variant_data)

        if ("dp" %in% names(variant_data)) {
            variant_data <- variant_data[
                !is.na(variant_data$dp) &
                    variant_data$dp >= filters$min_dp,
            ]
            n_after <- nrow(variant_data)
            if (verbose) {
                message(sprintf(
                    "Filtered by minimum DP >= %d: removed %d variants (%d remaining)",
                    filters$min_dp,
                    n_before - n_after,
                    n_after
                ))
            }
        } else {
            warning("Column 'dp' not found, skipping DP filter")
        }
    }

    # Step 6.5: Filter by unique read depth (UD) if specified
    # UD is annotated by bam-readcount and excludes duplicate reads
    if (!is.null(filters$min_unique_reads)) {
        n_before <- nrow(variant_data)

        if ("ud" %in% names(variant_data)) {
            # Keep variants with sufficient unique reads OR where UD is not available
            # (to avoid losing variants from samples processed before UD annotation)
            variant_data <- variant_data[
                is.na(variant_data$ud) |
                    variant_data$ud >= filters$min_unique_reads,
            ]

            n_after <- nrow(variant_data)
            ud_available <- sum(!is.na(variant_data$ud))

            if (verbose) {
                message(sprintf(
                    "Filtered by minimum %d unique reads (UD): removed %d variants (%d remaining, %d with UD annotation)",
                    filters$min_unique_reads,
                    n_before - n_after,
                    n_after,
                    ud_available
                ))
            }
        } else {
            warning(
                "UD (unique depth) column not found. ",
                "Ensure VCFs were annotated with bam-readcount in variant_calling_par.sh. ",
                "Skipping unique reads filter."
            )
        }
    }

    # Step 7: Filter by number of callers (with sample-type-specific thresholds)
    if (!is.null(filters$min_callers)) {
        n_before <- nrow(variant_data)

        if ("num_callers" %in% names(variant_data)) {
            # Check if we need to apply sample-type-specific thresholds
            if (!is.null(filters$min_callers_tumor_only)) {
                # Check if 'matched' column exists (should be merged from metadata)
                if ("matched" %in% names(variant_data)) {
                    # Use existing matched column with case-insensitive comparison
                    tumor_only_mask <- !is.na(variant_data[["matched"]]) &
                        tolower(variant_data[["matched"]]) == "no"
                    matched_mask <- !is.na(variant_data[["matched"]]) &
                        tolower(variant_data[["matched"]]) == "yes"
                    unknown_mask <- is.na(variant_data[["matched"]])

                    # Filter by different thresholds
                    # Variants with unknown sample type use global threshold
                    variant_data <- variant_data[
                        !(tumor_only_mask &
                            variant_data$num_callers <
                                filters$min_callers_tumor_only) &
                            !(matched_mask &
                                variant_data$num_callers <
                                    filters$min_callers) &
                            !(unknown_mask &
                                variant_data$num_callers < filters$min_callers),
                    ]

                    n_after <- nrow(variant_data)
                    if (verbose) {
                        message(sprintf(
                            "Filtered by minimum callers (matched: %d, tumor-only: %d): removed %d variants (%d remaining)",
                            filters$min_callers,
                            filters$min_callers_tumor_only,
                            n_before - n_after,
                            n_after
                        ))
                    }
                } else {
                    warning(
                        "Sample-type-specific caller thresholds requested but 'matched' column not found in variant_data. ",
                        "Using global threshold instead."
                    )
                    # Fallback to global threshold
                    variant_data <- variant_data[
                        variant_data$num_callers >= filters$min_callers,
                    ]
                    n_after <- nrow(variant_data)
                    if (verbose) {
                        message(sprintf(
                            "Filtered by minimum %d callers: removed %d variants (%d remaining)",
                            filters$min_callers,
                            n_before - n_after,
                            n_after
                        ))
                    }
                }
            } else {
                # Apply global threshold
                variant_data <- variant_data[
                    variant_data$num_callers >= filters$min_callers,
                ]

                n_after <- nrow(variant_data)
                if (verbose) {
                    message(sprintf(
                        "Filtered by minimum %d callers: removed %d variants (%d remaining)",
                        filters$min_callers,
                        n_before - n_after,
                        n_after
                    ))
                }
            }
        } else {
            warning("Column 'num_callers' not found, skipping caller filter")
        }
    }

    # Step 8: Filter by specific callers if provided
    if (
        !is.null(filters$specific_callers) &&
            length(filters$specific_callers) > 0
    ) {
        n_before <- nrow(variant_data)

        # Build condition more efficiently using Reduce
        caller_cols <- paste0(filters$specific_callers, "_called")
        caller_cols <- caller_cols[caller_cols %in% names(variant_data)]

        if (length(caller_cols) > 0) {
            caller_mask <- Reduce(
                `|`,
                lapply(caller_cols, function(col) {
                    !is.na(variant_data[[col]]) &
                        variant_data[[col]] != 0 &
                        variant_data[[col]] != FALSE
                })
            )
            variant_data <- variant_data[caller_mask, ]
        }

        n_after <- nrow(variant_data)
        if (verbose) {
            message(sprintf(
                "Filtered by specific callers (%s): removed %d variants (%d remaining)",
                paste(filters$specific_callers, collapse = ", "),
                n_before - n_after,
                n_after
            ))
        }
    }

    # Step 8.5: Filter requiring at least one caller from required_caller_any
    if (
        !is.null(filters$required_caller_any) &&
            length(filters$required_caller_any) > 0
    ) {
        n_before <- nrow(variant_data)

        caller_cols <- paste0(filters$required_caller_any, "_called")
        caller_cols <- caller_cols[caller_cols %in% names(variant_data)]

        if (length(caller_cols) > 0) {
            any_caller_mask <- Reduce(
                `|`,
                lapply(caller_cols, function(col) {
                    !is.na(variant_data[[col]]) &
                        variant_data[[col]] != 0 &
                        variant_data[[col]] != FALSE
                })
            )
            variant_data <- variant_data[any_caller_mask, ]
        } else {
            warning("None of the required_caller_any columns found in data")
        }

        n_after <- nrow(variant_data)
        if (verbose) {
            message(sprintf(
                "Required any-caller filter (%s): removed %d variants (%d remaining)",
                paste(filters$required_caller_any, collapse = ", "),
                n_before - n_after,
                n_after
            ))
        }
    }

    # Step 8.6: Filter requiring ClinVar Pathogenic or Likely_pathogenic
    if (isTRUE(filters$require_clinvar_pathogenic)) {
        n_before <- nrow(variant_data)

        if ("clinvar_clnsig" %in% names(variant_data)) {
            clnsig <- variant_data$clinvar_clnsig
            clinvar_path_mask <- !is.na(clnsig) &
                grepl("pathogenic", clnsig, ignore.case = TRUE) &
                !grepl("conflicting", clnsig, ignore.case = TRUE)
            variant_data <- variant_data[clinvar_path_mask, ]
        } else {
            warning("Column 'clinvar_clnsig' not found, skipping require_clinvar_pathogenic filter")
        }

        n_after <- nrow(variant_data)
        if (verbose) {
            message(sprintf(
                "ClinVar pathogenic filter: removed %d variants (%d remaining)",
                n_before - n_after,
                n_after
            ))
        }
    }
    # Step 8.7: Filter requiring a COSMIC ID (id column not NA and not ".")
    if (isTRUE(filters$require_cosmic_id)) {
        n_before <- nrow(variant_data)

        if ("id" %in% names(variant_data)) {
            variant_data <- variant_data[
                !is.na(variant_data$id) &
                    variant_data$id != ".",
            ]
        } else {
            warning("Column 'id' not found, skipping require_cosmic_id filter")
        }

        n_after <- nrow(variant_data)
        if (verbose) {
            message(sprintf(
                "COSMIC ID filter (id not NA/\".\"): removed %d variants (%d remaining)",
                n_before - n_after,
                n_after
            ))
        }
    }
    # Step 8.8: Filter requiring COSMIC CGC membership
    if (isTRUE(filters$require_cosmic_cgc)) {
        n_before <- nrow(variant_data)

        if ("cosmic_cgc_v99" %in% names(variant_data)) {
            variant_data <- variant_data[
                !is.na(variant_data$cosmic_cgc_v99) &
                    variant_data$cosmic_cgc_v99 == "Yes",
            ]
        } else {
            warning("Column 'cosmic_cgc_v99' not found, skipping require_cosmic_cgc filter")
        }

        n_after <- nrow(variant_data)
        if (verbose) {
            message(sprintf(
                "COSMIC CGC filter (cosmic_cgc_v99 == 'Yes'): removed %d variants (%d remaining)",
                n_before - n_after,
                n_after
            ))
        }
    }

    # Step 9: Filter by variant impact and damaging predictions
    # This step applies ADDITIONAL caller/prediction requirements beyond the
    # baseline min_callers filter (Step 7). Different criteria is followed for indel or snv
    #
    # INDEL PATH (when treat_coding_indels_as_high = TRUE and is_indel column present):
    #   Only HIGH coding INDELs use the HIGH caller threshold.
    #   Rationale: INDELs never have dbNSFP data; caller agreement is the only quality
    #   signal. In-frame and splice-region INDELs are as clinically relevant as
    #   frameshift INDELs (e.g., EGFR exon 19 del, ERBB2 exon 20 ins).
    #   - impact HIGH + is_indel → min_callers_high_impact
    #   - impact/MODERATE/LOW/MODIFIER + is_indel → excluded (intronic, intergenic, UTR)
    #
    # SNV PATH (is_indel = FALSE):
    #   HIGH impact (frameshifts, stop-gain/loss):
    #     Filtered by HIGHER caller count: min_callers_high_impact
    #     If being tumor only: min_callers_high_impact_tumor_only
    #   MODERATE impact (missense):
    #     Requires dbNSFP
    #     requires num_damaging >= min_damaging (>=x out of 13 tools)
    #       13 tools: 9 _pred (SIFT, SIFT4G, PolyPhen2-HDIV, PolyPhen2-HVAR,
    #       MutationTaster, MutationAssessor, PROVEAN, MetaSVM, ClinPred)
    #       + 3 _score (REVEL, MutScore, DANN) + 1 _phred (CADD); AlphaMissense excluded
    #    For matched/unknown uses min_callers. For tumor_only it uses min_callers_tumor_only
    #   LOW/MODIFIER SNV: excluded (not in allowed_impacts)
    if (!is.null(filters$allowed_impacts) || !is.null(filters$min_damaging)) {
        n_before <- nrow(variant_data)

        # Build filter condition
        keep_mask <- rep(TRUE, nrow(variant_data))

        # Check if impact column exists
        impact_col <- NULL
        if ("impact" %in% names(variant_data)) {
            impact_col <- "impact"
        } else if ("IMPACT" %in% names(variant_data)) {
            impact_col <- "IMPACT"
        }

        if (!is.null(impact_col)) {
            # Detect INDELs: use is_indel flag from merge_variants() if available.
            # Fallback: compute from ref/alt lengths if is_indel column is absent
            # (supports checkpoints created before is_indel was added to parse_vcf_info).
            if (
                !"is_indel" %in% names(variant_data) &&
                    "ref" %in% names(variant_data) &&
                    "alt" %in% names(variant_data)
            ) {
                variant_data[,
                    is_indel := nchar(ref) != 1L | nchar(alt) != 1L
                ]
                if (verbose) {
                    message(
                        "  Note: is_indel column computed from ref/alt lengths (checkpoint predates merge_variants update)"
                    )
                }
            }
            indel_path_active <- isTRUE(filters$treat_coding_indels_as_high) &&
                "is_indel" %in% names(variant_data)
            is_indel_row <- if (indel_path_active) {
                !is.na(variant_data[["is_indel"]]) & variant_data[["is_indel"]]
            } else {
                rep(FALSE, nrow(variant_data))
            }

            # Apply different requirements based on impact level
            if (
                !is.null(filters$min_damaging) &&
                    "num_damaging" %in% names(variant_data) &&
                    "num_callers" %in% names(variant_data)
            ) {
                # Shared setup: impact values and dbNSFP flag
                impact_val <- variant_data[[impact_col]]
                impact_is_na <- is.na(impact_val)
                # Check for has_dbnsfp flag (added by transform_predictions_scores)
                # If absent, assume all variants have dbNSFP data (backward compat)
                has_dbnsfp_val <- if ("has_dbnsfp" %in% names(variant_data)) {
                    variant_data[["has_dbnsfp"]]
                } else {
                    rep(TRUE, nrow(variant_data))
                }

                # Check if we need sample-type-specific thresholds for HIGH impact
                if (!is.null(filters$min_callers_high_impact_tumor_only)) {
                    # Check if 'matched' column exists
                    if ("matched" %in% names(variant_data)) {
                        # Use existing matched column with case-insensitive comparison
                        is_tumor_only <- !is.na(variant_data[["matched"]]) &
                            tolower(variant_data[["matched"]]) == "no"
                        is_matched <- !is.na(variant_data[["matched"]]) &
                            tolower(variant_data[["matched"]]) == "yes"
                        is_unknown <- is.na(variant_data[["matched"]])

                        keep_mask <- (
                            # ==================================================
                            # INDEL PATH: all coding INDELs use HIGH threshold
                            # Coding = impact is HIGH, MODERATE, or LOW (not MODIFIER)
                            # ==================================================
                            # Coding INDEL HIGH - matched or unknown samples
                            (indel_path_active &
                                is_indel_row &
                                !impact_is_na &
                                (impact_val == "HIGH") &
                                (is_matched | is_unknown) &
                                (variant_data$num_callers >=
                                    filters$min_callers_high_impact)) |
                            # Coding INDEL HIGH - tumor-only samples
                            (indel_path_active &
                                is_indel_row &
                                !impact_is_na &
                                (impact_val == "HIGH") &
                                is_tumor_only &
                                (variant_data$num_callers >=
                                    filters$min_callers_high_impact_tumor_only)) |
                            # Coding INDEL MODERATE/LOW/MODIFIER as not having has_dbnsfp data, are dismissed
                            # ==================================================
                            # SNV PATH
                            # HIGH: caller count (dbNSFP not needed)
                            # MODERATE: num_damaging AND caller count
                            # ==================================================
                            # HIGH SNV - matched/unknown
                            (!impact_is_na &
                                (impact_val == "HIGH") &
                                !is_indel_row &
                                (is_matched | is_unknown) &
                                (variant_data$num_callers >=
                                    filters$min_callers_high_impact)) |
                            # HIGH SNV - tumor-only
                            (!impact_is_na &
                                (impact_val == "HIGH") &
                                !is_indel_row &
                                is_tumor_only &
                                (variant_data$num_callers >=
                                    filters$min_callers_high_impact_tumor_only)) |
                            # MODERATE SNV - matched/unknown: has_dbnsfp AND num_damaging AND caller count
                            (!impact_is_na &
                                (impact_val == "MODERATE") &
                                !is_indel_row &
                                has_dbnsfp_val &
                                (is_matched | is_unknown) &
                                (variant_data$num_damaging >= filters$min_damaging) &
                                (variant_data$num_callers >= filters$min_callers)) |
                            # MODERATE SNV - tumor-only: has_dbnsfp AND num_damaging AND caller count
                            (!impact_is_na &
                                (impact_val == "MODERATE") &
                                !is_indel_row &
                                has_dbnsfp_val &
                                is_tumor_only &
                                (variant_data$num_damaging >= filters$min_damaging) &
                                (variant_data$num_callers >= filters$min_callers_tumor_only))
                            # LOW and MODIFIER SNVs are excluded via allowed_impacts hard filter below
                        )
                    } else {
                        warning(
                            "Sample-type-specific HIGH impact thresholds requested but 'matched' column not found. ",
                            "Using global threshold instead."
                        )
                        keep_mask <- (
                            # Coding INDEL HIGH - global threshold
                            (indel_path_active &
                                is_indel_row &
                                !impact_is_na &
                                (impact_val == "HIGH") &
                                (variant_data$num_callers >=
                                    filters$min_callers_high_impact)) |
                            # HIGH SNV: caller count (dbNSFP irrelevant)
                            (!impact_is_na &
                                (impact_val == "HIGH") &
                                !is_indel_row &
                                (variant_data$num_callers >=
                                    filters$min_callers_high_impact)) |
                            # MODERATE SNV: has_dbnsfp AND num_damaging AND caller count
                            (!impact_is_na &
                                (impact_val == "MODERATE") &
                                !is_indel_row &
                                has_dbnsfp_val &
                                (variant_data$num_damaging >= filters$min_damaging) &
                                (variant_data$num_callers >= filters$min_callers))
                            # LOW and MODIFIER SNVs are excluded via allowed_impacts hard filter below
                        )
                    }
                } else {
                    # Global threshold (no sample-type-specific thresholds)
                    keep_mask <- (
                        # Coding INDEL HIGH - global threshold
                        (indel_path_active &
                            is_indel_row &
                            !impact_is_na &
                            (impact_val == "HIGH") &
                            (variant_data$num_callers >=
                                filters$min_callers_high_impact)) |
                        # Coding INDEL MODERATE - global threshold (requires has_dbnsfp)
                        (indel_path_active &
                            is_indel_row &
                            !impact_is_na &
                            (impact_val == "MODERATE") &
                            has_dbnsfp_val &
                            (variant_data$num_callers >=
                                filters$min_callers)) |
                        # HIGH SNV: caller count (dbNSFP irrelevant)
                        (!impact_is_na &
                            (impact_val == "HIGH") &
                            !is_indel_row &
                            (variant_data$num_callers >=
                                filters$min_callers_high_impact)) |
                        # MODERATE SNV: has_dbnsfp AND num_damaging AND caller count
                        (!impact_is_na &
                            (impact_val == "MODERATE") &
                            !is_indel_row &
                            has_dbnsfp_val &
                            (variant_data$num_damaging >= filters$min_damaging) &
                            (variant_data$num_callers >= filters$min_callers))
                        # LOW and MODIFIER SNVs are excluded via allowed_impacts hard filter below
                    )
                }
            }

            # Hard impact filter: enforce allowed_impacts as an absolute constraint
            # for ALL variants including INDELs. This removes LOW and MODIFIER
            # regardless of any path above.
            if (!is.null(filters$allowed_impacts)) {
                keep_mask <- keep_mask &
                    (variant_data[[impact_col]] %in% filters$allowed_impacts)
            }

            variant_data <- variant_data[keep_mask, ]

            n_after <- nrow(variant_data)
            if (verbose) {
                # Report MODERATE SNVs rescued by caller-count fallback (no dbNSFP)
                if ("has_dbnsfp" %in% names(variant_data)) {
                    n_mod_no_dbnsfp <- sum(
                        !is.na(variant_data[[impact_col]]) &
                            variant_data[[impact_col]] == "MODERATE" &
                            !is_indel_row[keep_mask] &
                            !variant_data[["has_dbnsfp"]]
                    )
                    if (n_mod_no_dbnsfp > 0) {
                        message(sprintf(
                            "  Note: %d MODERATE SNVs without dbNSFP data kept via caller-count fallback",
                            n_mod_no_dbnsfp
                        ))
                    }
                }
                # Report INDEL counts if indel_path_active
                if (indel_path_active && "is_indel" %in% names(variant_data)) {
                    n_indels_kept <- sum(
                        !is.na(variant_data[["is_indel"]]) &
                            variant_data[["is_indel"]]
                    )
                    if (n_indels_kept > 0) {
                        message(sprintf(
                            "  Note: %d coding INDELs kept (HIGH, min_callers_high_impact threshold)",
                            n_indels_kept
                        ))
                    }
                }
                message(sprintf(
                    "Filtered by impact and damaging predictions: removed %d variants (%d remaining)",
                    n_before - n_after,
                    n_after
                ))
            }
        } else {
            warning("Impact column not found, skipping impact filter")
        }
    }

    # Step 10: Blacklist recurrent variants (systematic artifacts)
    # Now that variants have passed quality and impact filters, identify those
    # where the exact same chr:pos:ref:alt appears in >max_recurrence_fraction
    # of samples - likely mapping artifacts (e.g., pseudogene misalignment,
    # repetitive regions). Calculating recurrence post-quality-filtering makes
    # the blacklist more interpretable (only high-quality recurrent variants).
    if (!is.null(filters$max_recurrence_fraction)) {
        n_before <- nrow(variant_data)

        # Validate threshold
        if (
            filters$max_recurrence_fraction <= 0 ||
                filters$max_recurrence_fraction >= 1
        ) {
            warning(
                "max_recurrence_fraction must be between 0 and 1 (exclusive). Skipping recurrence filter."
            )
        } else {
            # Count total unique samples in the dataset
            n_total_samples <- data.table::uniqueN(
                variant_data[[variant_sample_col]]
            )
            max_samples <- floor(
                n_total_samples * filters$max_recurrence_fraction
            )

            if (verbose) {
                message(sprintf(
                    "Recurrence blacklist: flagging variants in >%.0f%% of samples (>%d/%d samples)...",
                    filters$max_recurrence_fraction * 100,
                    max_samples,
                    n_total_samples
                ))
            }

            # Required columns for variant identity
            id_cols <- c("number_chrom", "pos", "ref", "alt")
            missing_id_cols <- setdiff(id_cols, names(variant_data))

            if (length(missing_id_cols) > 0) {
                warning(sprintf(
                    "Columns required for recurrence filter missing: %s. Skipping.",
                    paste(missing_id_cols, collapse = ", ")
                ))
            } else {
                # Count samples per unique variant position
                recurrence_counts <- variant_data[,
                    .(
                        n_samples = data.table::uniqueN(
                            get(variant_sample_col)
                        )
                    ),
                    by = id_cols
                ]

                # Identify blacklisted positions
                blacklisted <- recurrence_counts[
                    n_samples > max_samples
                ]

                if (nrow(blacklisted) > 0) {
                    # Log blacklisted variants
                    if (verbose) {
                        message(sprintf(
                            "  Found %d variant positions in >%d samples:",
                            nrow(blacklisted),
                            max_samples
                        ))
                        # Show top offenders
                        top_bl <- blacklisted[order(-n_samples)][
                            1:min(5, nrow(blacklisted))
                        ]
                        for (i in seq_len(nrow(top_bl))) {
                            message(sprintf(
                                "    %s:%d %s>%s (%d samples)",
                                top_bl$number_chrom[i],
                                top_bl$pos[i],
                                top_bl$ref[i],
                                top_bl$alt[i],
                                top_bl$n_samples[i]
                            ))
                        }
                    }

                    # Export blacklisted variants before removal
                    # IMPORTANT: Only variants that SURVIVED previous filtering steps
                    # and are being removed SPECIFICALLY by the blacklisting step are exported
                    if (!is.null(export_blacklisted_path)) {
                        # Get all occurrences of blacklisted positions that are
                        # in the CURRENT variant_data (i.e., survived previous filters)
                        blacklisted_variants <- variant_data[
                            blacklisted,
                            on = id_cols
                        ]

                        # Verify that we're only exporting variants being removed by THIS step
                        if (verbose) {
                            message(sprintf(
                                "  Exporting blacklist: %d positions affecting %d variants (these survived previous filters)",
                                nrow(blacklisted),
                                nrow(blacklisted_variants)
                            ))
                        }

                        caller_cols_bl <- grep(
                            "_called$",
                            names(blacklisted_variants),
                            value = TRUE
                        )
                        callers_bl <- gsub(
                            "_called$",
                            "",
                            caller_cols_bl
                        )

                        # Create summary: one row per blacklisted position
                        # with aggregated stats from ALL occurrences that survived previous filters
                        blacklist_summary <- blacklisted_variants[,
                            c(
                                list(
                                    gene = paste(
                                        unique(gene),
                                        collapse = "; "
                                    ),
                                    annotation = paste(
                                        unique(annotation),
                                        collapse = "; "
                                    ),
                                    n_samples = data.table::uniqueN(
                                        get(variant_sample_col)
                                    ),
                                    pct_samples = round(
                                        data.table::uniqueN(
                                            get(variant_sample_col)
                                        ) /
                                            n_total_samples *
                                            100,
                                        1
                                    ),
                                    mean_vaf_pct = round(
                                        mean(af * 100, na.rm = TRUE),
                                        1
                                    ),
                                    mean_callers = round(
                                        mean(num_callers, na.rm = TRUE),
                                        1
                                    )
                                ),
                                lapply(
                                    setNames(caller_cols_bl, callers_bl),
                                    function(col) {
                                        round(
                                            sum(get(col), na.rm = TRUE) /
                                                .N *
                                                100,
                                            1
                                        )
                                    }
                                )
                            ),
                            by = id_cols
                        ]

                        data.table::setnames(
                            blacklist_summary,
                            callers_bl,
                            paste0("pct_", callers_bl)
                        )
                        data.table::setorder(
                            blacklist_summary,
                            -n_samples,
                            -pct_samples
                        )

                        # Ensure output directory exists
                        export_dir <- dirname(export_blacklisted_path)
                        if (!dir.exists(export_dir)) {
                            dir.create(
                                export_dir,
                                recursive = TRUE,
                                showWarnings = FALSE
                            )
                        }

                        data.table::fwrite(
                            blacklist_summary,
                            file = export_blacklisted_path,
                            sep = "\t",
                            quote = FALSE
                        )
                        if (verbose) {
                            n_vars_removed <- nrow(blacklisted_variants)
                            message(sprintf(
                                "  \u2713 Exported %d blacklisted positions (%d variants being removed by blacklisting) to: %s",
                                nrow(blacklist_summary),
                                n_vars_removed,
                                export_blacklisted_path
                            ))
                            message(sprintf(
                                "  Note: Recurrence calculated from %d samples (post-quality and impact filtering)",
                                n_total_samples
                            ))
                        }
                    }

                    # Remove blacklisted variants via anti-join
                    variant_data <- variant_data[
                        !blacklisted,
                        on = id_cols
                    ]

                    n_after <- nrow(variant_data)
                    if (verbose) {
                        message(sprintf(
                            "  Recurrence blacklist: removed %s variants from %d positions (%s remaining)",
                            format(n_before - n_after, big.mark = ","),
                            nrow(blacklisted),
                            format(n_after, big.mark = ",")
                        ))
                    }
                } else {
                    if (verbose) {
                        message(
                            "  No recurrent artifact positions detected"
                        )
                    }
                }
            }
        }
    }

    # Step 11: Filter by cancer gene presence
    if (filters$require_cancer_gene) {
        n_before <- nrow(variant_data)

        if ("cgc_gene" %in% names(variant_data)) {
            variant_data <- variant_data[cgc_gene == TRUE]
            n_after <- nrow(variant_data)
            if (verbose) {
                message(sprintf(
                    "Filtered by cancer gene presence (cgc_gene == TRUE): removed %d variants (%d remaining)",
                    n_before - n_after,
                    n_after
                ))
            }
        } else {
            warning("Column 'cgc_gene' not found, skipping cancer gene filter")
        }
    }

    # Step 12: Remove artifacts
    if (filters$remove_artifacts) {
        n_before <- nrow(variant_data)
        n_sob_removed <- 0
        n_ffperase_removed <- 0

        # Check for SOB (SeqObjectiveBias) prediction column
        sob_col <- NULL
        if ("sob_prediction" %in% names(variant_data)) {
            sob_col <- "sob_prediction"
        } else if ("SOB_PREDICTION" %in% names(variant_data)) {
            sob_col <- "SOB_PREDICTION"
        }

        if (!is.null(sob_col)) {
            n_sob_before <- nrow(variant_data)
            keep_rows <- is.na(variant_data[[sob_col]]) |
                variant_data[[sob_col]] != "artifact"
            variant_data <- variant_data[keep_rows, ]
            n_sob_removed <- n_sob_before - nrow(variant_data)

            if (verbose && n_sob_removed > 0) {
                message(sprintf(
                    "SOBDetector: removed %d artifact variants",
                    n_sob_removed
                ))
            }
        }

        # Check for FFPErase prediction column
        ffperase_col <- NULL
        if ("ffperase_prediction" %in% names(variant_data)) {
            ffperase_col <- "ffperase_prediction"
        } else if ("ffperase_artifact" %in% names(variant_data)) {
            ffperase_col <- "ffperase_artifact"
        } else if ("FFPERASE_PREDICTION" %in% names(variant_data)) {
            ffperase_col <- "FFPERASE_PREDICTION"
        }

        if (!is.null(ffperase_col)) {
            n_ffperase_before <- nrow(variant_data)
            keep_rows <- is.na(variant_data[[ffperase_col]]) |
                variant_data[[ffperase_col]] != "True"
            variant_data <- variant_data[keep_rows, ]
            n_ffperase_removed <- n_ffperase_before - nrow(variant_data)

            if (verbose && n_ffperase_removed > 0) {
                message(sprintf(
                    "FFPErase: removed %d artifact variants",
                    n_ffperase_removed
                ))
            }
        }

        n_after <- nrow(variant_data)
        n_total_removed <- n_before - n_after

        if (verbose && n_total_removed > 0) {
            message(sprintf(
                "Total artifacts removed: %d variants (%d remaining)",
                n_total_removed,
                n_after
            ))
        } else if (verbose && n_total_removed == 0) {
            message(
                "No artifact prediction columns found or no artifacts detected"
            )
        }
    }

    # Final summary
    n_final <- nrow(variant_data)
    n_removed <- n_initial - n_final
    pct_removed <- round(100 * n_removed / n_initial, 2)

    if (verbose) {
        message(sprintf(
            "\nFiltering complete: %d/%d variants retained (%.2f%% removed)",
            n_final,
            n_initial,
            pct_removed
        ))
    }

    return(variant_data)
}

#' Filter and Annotate Gene Variants with Domain Information
#'
#' Extracts variants for a specific gene from merged WES variant data and
#' annotates them with functional domain information. Works with output from
#' \code{merge_variants()} and integrates seamlessly with existing squamr
#' filtering workflows. This is a general-purpose function that can be used
#' for any gene with defined protein domains.
#'
#' @param variant_data A data frame or tibble from \code{merge_variants()},
#'   containing at minimum: gene, hgvs_p (protein change notation). If a
#'   \code{hgvs_c} column is present, it is used as a fallback to infer
#'   \code{aa_pos} for splice variants where \code{hgvs_p} is NA.
#' @param gene_name Character string, name of the gene to filter (e.g., "TP53",
#'   "PIK3CA"). Must match the gene column in variant_data.
#' @param domains Named list of domain boundaries (amino acid positions).
#'   Each element should be a numeric vector of length 2: c(start, end).
#'   Must include at least one domain; if a "primary" domain is specified,
#'   it will be used for is_primary_domain classification.
#' @param primary_domain Character string, name of the primary domain of
#'   interest for creating binary classification (default: NULL). If specified,
#'   must exist in the domains list. Creates is_primary_domain and is_non_primary
#'   logical columns.
#'
#' @return A data frame with gene variants annotated with additional columns:
#'   \itemize{
#'     \item \code{aa_pos}: Amino acid position extracted from hgvs_p, or
#'       inferred from hgvs_c for splice variants (ceiling of coding position / 3)
#'     \item \code{domain}: Functional domain name (or "Other"/"Unknown")
#'     \item \code{is_primary_domain}: Logical, TRUE if variant is in primary domain
#'       (only if primary_domain specified)
#'     \item \code{is_non_primary}: Logical, TRUE if variant is in other functional
#'       domains but not primary (only if primary_domain specified)
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # TP53 with DNA-binding domain as primary
#' tp53_domains <- list(
#'   transactivation = c(1, 61),
#'   proline_rich = c(62, 101),
#'   dbd = c(102, 292),
#'   tetramerization = c(323, 356),
#'   c_terminal = c(357, 393)
#' )
#' tp53_vars <- filter_gene_variants(
#'   merged_variants,
#'   gene_name = "TP53",
#'   domains = tp53_domains,
#'   primary_domain = "dbd"
#' )
#'
#' # PIK3CA with helical domain as primary
#' pik3ca_domains <- list(
#'   p85_binding = c(1, 108),
#'   c2 = c(330, 487),
#'   helical = c(533, 693),
#'   kinase = c(697, 1068)
#' )
#' pik3ca_vars <- filter_gene_variants(
#'   merged_variants,
#'   gene_name = "PIK3CA",
#'   domains = pik3ca_domains,
#'   primary_domain = "helical"
#' )
#' }
filter_gene_variants <- function(
    variant_data,
    gene_name,
    domains,
    primary_domain = NULL
) {
    # Validate inputs
    if (!is.data.frame(variant_data)) {
        stop("'variant_data' must be a data frame.")
    }
    if (!is.character(gene_name) || length(gene_name) != 1) {
        stop("'gene_name' must be a single character string.")
    }
    if (!is.list(domains) || length(domains) == 0) {
        stop("'domains' must be a named list with at least one domain.")
    }
    if (!is.null(primary_domain) && !primary_domain %in% names(domains)) {
        stop(sprintf(
            "'primary_domain' (%s) not found in domains list.",
            primary_domain
        ))
    }

    required_cols <- c("gene", "hgvs_p")
    missing_cols <- setdiff(required_cols, colnames(variant_data))
    if (length(missing_cols) > 0) {
        stop(
            "Required column(s) missing: ",
            paste(missing_cols, collapse = ", ")
        )
    }

    # Convert to data.table if not already
    if (!data.table::is.data.table(variant_data)) {
        data.table::setDT(variant_data)
    } else {
        variant_data <- data.table::copy(variant_data)
    }

    # Filter for specified gene and add NA check
    gene_vars <- variant_data[gene == gene_name & !is.na(gene)]

    # Extract amino acid position from hgvs_p (primary) or hgvs_c (fallback
    # for splice variants where hgvs_p is NA).
    # For splice variants, the affected codon is derived from the coding position
    # in hgvs_c: c.NNN+d or c.NNN-d -> aa = ceiling(NNN / 3).
    gene_vars[, aa_pos := as.numeric(sub(".*?(\\d+).*", "\\1", hgvs_p))]
    if ("hgvs_c" %in% names(gene_vars)) {
        gene_vars[
            is.na(aa_pos) & !is.na(hgvs_c),
            aa_pos := {
                # Extract the coding nucleotide position (before +/- offset)
                cdna_pos <- suppressWarnings(
                    as.numeric(sub(".*?c\\.([0-9]+)[+\\-\\*]?.*", "\\1", hgvs_c))
                )
                ceiling(cdna_pos / 3L)
            }
        ]
    }

    # Initialize domain column
    gene_vars[,
        domain := fcase(
            is.na(aa_pos) , "Unknown" ,
            default = "Other"
        )
    ]

    # Assign domains based on position
    for (domain_name in names(domains)) {
        domain_range <- domains[[domain_name]]
        if (length(domain_range) != 2 || !is.numeric(domain_range)) {
            warning(sprintf(
                "Skipping domain '%s': invalid range specification",
                domain_name
            ))
            next
        }
        # Capitalize first letter of domain name for display
        display_name <- paste0(
            toupper(substring(domain_name, 1, 1)),
            substring(domain_name, 2)
        )
        # Replace hyphen/underscore with space and capitalize
        display_name <- gsub("_", "-", display_name)
        display_name <- gsub("-", " ", display_name)
        display_name <- tools::toTitleCase(display_name)

        # Use data.table syntax for efficient assignment
        gene_vars[
            !is.na(aa_pos) &
                aa_pos >= domain_range[1] &
                aa_pos <= domain_range[2],
            domain := display_name
        ]
    }

    # Add primary domain classification if specified
    if (!is.null(primary_domain)) {
        primary_display <- tools::toTitleCase(gsub("[_-]", " ", primary_domain))
        gene_vars[, `:=`(
            is_primary_domain = (domain == primary_display),
            is_non_primary = (domain != primary_display &
                domain != "Other" &
                domain != "Unknown")
        )]
    }

    # Add NA check - count variants with NA in critical fields
    n_na_gene <- sum(is.na(gene_vars$gene))
    n_na_hgvs <- sum(is.na(gene_vars$hgvs_p))
    n_na_aa_pos <- sum(is.na(gene_vars$aa_pos))

    if (n_na_gene > 0 || n_na_hgvs > 0 || n_na_aa_pos > 0) {
        message(sprintf(
            "NA values detected in %s variants: gene=%d, hgvs_p=%d, aa_pos=%d",
            gene_name,
            n_na_gene,
            n_na_hgvs,
            n_na_aa_pos
        ))
    }

    return(gene_vars)
}
