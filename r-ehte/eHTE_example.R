# #############################################################################
# eHTE Analysis Example Script
#
# This script demonstrates how to run the eHTE analysis in R.
# The eHTE package assesses heterogeneity of treatment effects between
# study arms in clinical trials.
#
# Written for eHTE R implementation
# Date: 2024
# #############################################################################

# Source the eHTE functions
source("eHTE.R")

# =============================================================================
# EXAMPLE 1: Using simulated data
# =============================================================================

# Create example data with:
# - TRT01P: Treatment name (must include "Placebo")
# - TRT01PN: Treatment code (1 = Placebo, 2+ = Active treatments)
# - CHG: Change from baseline (the outcome variable)

set.seed(42)

# Simulate a clinical trial with placebo and two active treatment arms
n_placebo <- 100
n_trt1 <- 80
n_trt2 <- 85

# Placebo arm: mean = 0, sd = 12
placebo_data <- data.frame(
  TRT01P = "Placebo",
  TRT01PN = 1L,
  CHG = rnorm(n_placebo, mean = 0, sd = 12)
)

# Treatment 1: mean = -5 (improvement), some heterogeneity
# Mix of responders and non-responders
trt1_responders <- rnorm(n_trt1 * 0.6, mean = -10, sd = 8)
trt1_nonresponders <- rnorm(n_trt1 * 0.4, mean = 2, sd = 8)
trt1_chg <- c(trt1_responders, trt1_nonresponders)

treatment1_data <- data.frame(
  TRT01P = "Treatment Low Dose",
  TRT01PN = 2L,
  CHG = trt1_chg
)

# Treatment 2: mean = -8 (more improvement), more heterogeneity
trt2_responders <- rnorm(n_trt2 * 0.5, mean = -15, sd = 10)
trt2_nonresponders <- rnorm(n_trt2 * 0.5, mean = 0, sd = 10)
trt2_chg <- c(trt2_responders, trt2_nonresponders)

treatment2_data <- data.frame(
  TRT01P = "Treatment High Dose",
  TRT01PN = 3L,
  CHG = trt2_chg
)

# Combine into single dataset
example_data <- rbind(placebo_data, treatment1_data, treatment2_data)

# View the data structure
cat("\n=== Example Data Structure ===\n")
str(example_data)
cat("\n=== Data Summary ===\n")
print(summary(example_data))

# =============================================================================
# Run eHTE Analysis
# =============================================================================

cat("\n\n=== Running eHTE Analysis ===\n")
cat("Note: Using 1000 permutations for faster demo. Use 10000 for production.\n\n")

# Run the analysis (use fewer permutations for demo)
results <- eHTE_p(example_data, n_perms = 1000)

# =============================================================================
# View Results
# =============================================================================

cat("\n\n=== Results Summary ===\n")
cat("\nAll Data Percentile Method:\n")
print(results$results_all)

cat("\nFixed 48 Percentile Method:\n")
print(results$results_48)

# =============================================================================
# Generate Plots
# =============================================================================

cat("\n\n=== Generating Plots ===\n")

# Individual plots
p_hist <- plot_hist(example_data, title = "Distribution of Change from Baseline")
p_cum <- plot_cumulative(example_data, title = "Cumulative Distribution by Treatment")
p_ite <- plot_ite(results$ite_df, title = "Individual Treatment Effects")

# Save plots
ggsave("histogram.png", p_hist, width = 8, height = 6)
ggsave("cumulative.png", p_cum, width = 8, height = 6)
ggsave("ite.png", p_ite, width = 6, height = 8)

cat("Plots saved: histogram.png, cumulative.png, ite.png\n")

# Combined 3-panel plot (similar to Python eHTEplot)
cat("\nGenerating combined 3-panel plot...\n")
combined <- eHTE_plot(example_data, results$ite_df, var_name = "CHG")

# Save combined plot
ggsave("eHTE_combined.png", combined, width = 12, height = 6)
cat("Combined plot saved: eHTE_combined.png\n")

# =============================================================================
# EXAMPLE 2: Reading data from file
# =============================================================================

cat("\n\n=== Example: Reading Data from File ===\n")
cat("
# If you have data in a CSV file:
# my_data <- read.csv('path/to/your/data.csv')

# If you have data in a SAS dataset (requires haven package):
# library(haven)
# my_data <- read_sas('path/to/your/data.sas7bdat')

# Ensure your data has the required columns:
# - TRT01P (character): Treatment name, must include 'Placebo'
# - TRT01PN (integer): Treatment code, 1 = Placebo, 2+ = Active
# - CHG (numeric): Change from baseline outcome

# Then run:
# results <- eHTE_p(my_data, n_perms = 10000)
")

# =============================================================================
# Interpretation Guide
# =============================================================================

cat("\n=== Interpretation Guide ===\n")
cat("
sigma: Standard deviation of Individual Treatment Effects (ITE)
       Higher values indicate more variability in treatment response

eHTE:  Standardized effect size = sigma / sigma_placebo
       Values > 1 suggest treatment response variability exceeds
       what would be expected from placebo variability alone

p-value: Permutation-based significance test
         H0: Treatment effect is homogeneous (no HTE)
         H1: Treatment effect is heterogeneous
         Small p-values (< 0.05) suggest significant heterogeneity

Plots:
- Histogram: Shows distribution of outcomes by treatment arm
- Cumulative: Shows cumulative distribution curves
- ITE scatter: Shows individual treatment effects by percentile rank
  Points far from 0 indicate responders (negative) or non-responders (positive)
")

cat("\n=== Analysis Complete ===\n")
