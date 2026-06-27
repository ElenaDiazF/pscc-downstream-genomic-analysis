# SAMPLE OUTLIER DETERMINATION
#
# Documentation
# - This script does not define custom functions.
# - It calculates tumor mutational burden (TMB), applies Tukey's upper
#   outlier cutoff, and lists/plots hypermutated samples.
#
# Run from 1 to 6 rmd files with filtered_Variants_outlier.
#Normalmente para detectar hypermutated outliers en TMB se usa el criterio de Tukey:
  # Upper cutoff = Q3 + 2.5 × IQR
#Q3 = percentil 75
#IQR = Q3 − Q1
#TMB > upper cutoff--> hypermutated

# Include premalignant samples in the TMB table.
#INCLUDING PREMALIGNANT SAMPLES
tmb_table_all <- tmb(maf = maf_tumor, captureSize = 43, logScale = FALSE)

# Attach clinical metadata to each TMB value.
tmb_with_clinical_all <- tmb_table_all %>%
  inner_join(clinical_for_maf_all, by = "Tumor_Sample_Barcode")

# Vector con TMB
tmb_values <- tmb_with_clinical_all$total_perMB

# Calcular quartiles
# Calculate the quartiles used by Tukey's outlier criterion.
Q1 <- quantile(tmb_values, 0.25, na.rm = TRUE)
Q3 <- quantile(tmb_values, 0.75, na.rm = TRUE)

# IQR
# Calculate the interquartile range.
IQR_value <- IQR(tmb_values, na.rm = TRUE)

# Tukey upper cutoff
# Define Tukey's upper cutoff for hypermutated outliers.
upper_cutoff <- Q3 + 2.5 * IQR_value

upper_cutoff


# Identify samples above the TMB outlier cutoff.
hypermutated_samples <- tmb_with_clinical_all %>%
  filter(total_perMB > upper_cutoff)

hypermutated_samples

# Print the sample IDs classified as hypermutated.
hypermutated_samples$Tumor_Sample_Barcode
# "AS9750"

# Plot TMB per sample with the hypermutation cutoff.
ggplot(tmb_with_clinical, aes(x = Tumor_Sample_Barcode, y = total_perMB)) +
  geom_point(size = 2) +
  geom_hline(
    yintercept = upper_cutoff,
    linetype = "dashed",
    color = "red"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )


#Hypermutated outlier samples were identified using the Tukey outlier criterion (Q3 + 2.5 × IQR).
