# R eHTE - Heterogeneous Treatment Effect Estimator

R implementation of the eHTE (estimated Heterogeneous Treatment Effect) package for assessing treatment response variability in clinical trials.

## Overview

The eHTE package assesses heterogeneity of treatment effects between study arms. It calculates the variance in the differences in cumulative responses between the placebo and active treatment arms across various percentiles, testing the null hypothesis (H₀) of homogeneous treatment effect against the alternative hypothesis of heterogeneous treatment effect.

## Files

1. **eHTE.R**: Core functions for eHTE analysis
2. **eHTE_example.R**: Example script demonstrating usage

## Requirements

### R Version
- R >= 4.0.0

### Required Packages
```r
install.packages(c("dplyr", "tidyr", "ggplot2", "purrr", "gridExtra"))
```

Optional (for reading SAS datasets):
```r
install.packages("haven")
```

## Usage

### Basic Usage

```r
# Source the eHTE functions
source("eHTE.R")

# Prepare your data with required columns:
# - TRT01P (character): Treatment name, must include "Placebo"
# - TRT01PN (integer): Treatment code (1 = Placebo, 2+ = Active treatments)
# - CHG (numeric): Change from baseline outcome

# Run the analysis
results <- eHTE_p(your_data, n_perms = 10000)

# Generate plots
eHTE_plot(your_data, results$ite_df, var_name = "CHG")
```

### Input Data Format

| Column | Type | Description |
|--------|------|-------------|
| TRT01P | character | Treatment name (must include "Placebo") |
| TRT01PN | integer | Treatment code: 1 = Placebo, 2+ = Active |
| CHG | numeric | Change from baseline (outcome variable) |

Example data structure:
```r
data <- data.frame(
  TRT01P = c("Placebo", "Placebo", "Treatment", "Treatment"),
  TRT01PN = c(1L, 1L, 2L, 2L),
  CHG = c(-2.5, 1.3, -8.4, -3.2)
)
```

### Reading Data from Files

```r
# From CSV
data <- read.csv("your_data.csv")

# From SAS dataset (requires haven package)
library(haven)
data <- read_sas("your_data.sas7bdat")

# Ensure proper data types
data$TRT01PN <- as.integer(data$TRT01PN)
data$TRT01P <- as.character(data$TRT01P)
```

## Functions

### Main Analysis

| Function | Description |
|----------|-------------|
| `eHTE_p(indf, n_perms)` | Main analysis function. Returns list with results tables and ITE data |
| `validate_input(df)` | Validates input data structure |

### Plotting

| Function | Description |
|----------|-------------|
| `eHTE_plot(indf, ite_df, var_name)` | Combined 3-panel visualization |
| `plot_hist(indf, title)` | Histogram with kernel density overlay |
| `plot_cumulative(indf, title)` | Cumulative distribution curves |
| `plot_ite(ite_df, title)` | Individual treatment effect scatter plot |

### Internal Functions

| Function | Description |
|----------|-------------|
| `sas_percentile(x, probs)` | SAS-compatible percentile calculation (PCTLDEF=5) |
| `calculate_percentile(x)` | Rank-based percentile assignment |
| `calculate_sigma(pbo_df, trt_df)` | Sigma calculation using actual data percentiles |
| `calculate_sigma_48(pbo_df, trt_df)` | Sigma calculation using fixed 48 percentiles |
| `gen_sim(...)` | Generate simulated data for permutation testing |
| `cal_pvalue(...)` | Calculate p-values from permutation distribution |

## Output

### Results Tables

The `eHTE_p()` function returns a list containing:

- **results_all**: Results using actual data percentiles
- **results_48**: Results using fixed 48 percentiles (3 to 97 by 2)
- **ite_df**: Individual Treatment Effect data for plotting
- **summary_stats**: Summary statistics by treatment arm

Each results table contains:

| Column | Description |
|--------|-------------|
| TRT01PN | Treatment code |
| TRT01P | Treatment name |
| sigma (σ) | Standard deviation of ITE |
| eHTE | Standardized effect size (σ / σ_placebo) |
| pvalue | Permutation-based p-value |

### Example Output

```
=== Fixed 48 percentiles ===
treatment: Treatment, sigma: 4.15006, eHTE: 0.335, p-value: 0.062
```

| TRT01PN | TRT01P | sigma | eHTE | pvalue |
|---------|--------|-------|------|--------|
| 2 | Treatment | 4.15006 | 0.33547 | 0.0615 |

### Plots

1. **Histogram**: Distribution of outcomes with kernel density curves
2. **Cumulative**: Cumulative distribution curves by treatment arm
3. **ITE Scatter**: Individual treatment effects plotted by percentile rank

## Interpretation

- **sigma (σ)**: Higher values indicate more variability in treatment response
- **eHTE**: Values > 1 suggest treatment response variability exceeds placebo variability
- **p-value**: Tests H₀ (homogeneous treatment effect). Small values (< 0.05) suggest significant heterogeneity

## Algorithm Details

1. **Data Preparation**: Separate placebo and active treatment arms
2. **Percentile Calculation**: Uses SAS-compatible PCTLDEF=5 method
3. **ITE Calculation**: Treatment response minus matched placebo response at each percentile
4. **Sigma Estimation**: Standard deviation of ITE across percentiles 3-97
5. **Permutation Testing**: 10,000 simulations using observed means but placebo SD
6. **P-value**: Proportion of simulated σ ≥ observed σ

## Compatibility

This R implementation produces results equivalent to:
- Python `ehte` package
- SAS eHTE macros

The percentile calculation uses the same method as SAS `PROC UNIVARIATE` with `PCTLDEF=5`.

## References

Siegel JS, Zhong J, Tomioka S, Ogirala A, Faraone SV, Szabo ST, Koblan KS, Hopkins SC. Estimating heterogeneity of treatment effect in psychiatric clinical trials. medRxiv [Preprint]. 2024 Apr 23:2024.04.23.24306211. doi: 10.1101/2024.04.23.24306211. PMID: 38712180; PMCID: PMC11071592.

## License

Apache License 2.0

## Contact

For inquiries, contact:
- Joshua.siegel@us.sumitomo-pharma.com
- sam.tomioka@us.sumitomo-pharma.com
