# #############################################################################
# Copyright 2023-2024 Sumitomo Pharma America
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# For inquiry about this package, contact Joshua.siegel@us.sumitomo-pharma.com,
# sam.tomioka@us.sumitomo-pharma.com
#
# R version of eHTE estimator
# Equivalent to Python ehte package and SAS eHTE macros
#
# #############################################################################

# Required packages
required_packages <- c("dplyr", "tidyr", "ggplot2", "purrr", "gridExtra")

# Function to check and install packages if needed
check_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Package", pkg, "is not installed. Please install with: install.packages('", pkg, "')"))
    }
  }
}

check_packages(required_packages)

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(gridExtra)

#' SAS-compatible percentile calculation (PCTLDEF=5)
#'
#' Implements the empirical distribution function with averaging method
#' which matches SAS PROC UNIVARIATE with PCTLDEF=5
#'
#' @param x Numeric vector of values
#' @param probs Numeric vector of probabilities (0-1)
#' @return Numeric vector of percentile values
sas_percentile <- function(x, probs) {
  n <- length(x)
  x <- sort(x)

  result <- sapply(probs, function(p) {
    np_value <- n * p
    j <- floor(np_value)
    g <- np_value - j

    if (j == 0) {
      # Average of first two elements
      return(mean(x[1:min(2, n)]))
    } else if (g == 0) {
      if (j < n) {
        return(mean(c(x[j], x[j + 1])))
      } else {
        return(x[n])
      }
    } else {
      if (j < n) {
        return(x[j + 1])
      } else {
        return(x[n])
      }
    }
  })

  return(result)
}

#' Calculate rank-based percentile for observations
#'
#' Assigns percentile ranks to observations (similar to PROC RANK group=101)
#'
#' @param x Numeric vector
#' @return Integer vector of percentile ranks (0-100)
calculate_percentile <- function(x) {
  ranked <- rank(x, ties.method = "average")
  percentiles <- as.integer((ranked / length(x)) * 100)
  return(percentiles)
}

#' Validate input data frame
#'
#' Checks that the input data frame has required columns and proper values
#'
#' @param df Data frame with TRT01P, TRT01PN, CHG columns
#' @return TRUE if valid, stops with error message otherwise
validate_input <- function(df) {
  required_columns <- c("TRT01P", "TRT01PN", "CHG")

  # Check required columns
  if (!all(required_columns %in% names(df))) {
    missing <- setdiff(required_columns, names(df))
    stop(paste("DataFrame is missing required columns:", paste(missing, collapse = ", ")))
  }

  # Check TRT01P is character/string
  if (!is.character(df$TRT01P)) {
    stop("TRT01P column should contain only strings/characters.")
  }

  # Check TRT01PN is numeric/integer
  if (!is.numeric(df$TRT01PN)) {
    stop("TRT01PN column should contain only integers.")
  }

  # Check for 'Placebo' in TRT01P
  if (!"Placebo" %in% df$TRT01P) {
    stop("TRT01P column should contain the value 'Placebo'.")
  }

  # Check TRT01PN is 1 when TRT01P is 'Placebo'
  placebo_rows <- df[df$TRT01P == "Placebo", ]
  if (!any(placebo_rows$TRT01PN == 1)) {
    stop("TRT01PN should be 1 when TRT01P is 'Placebo'.")
  }

  # Check CHG is numeric
  if (!is.numeric(df$CHG)) {
    stop("CHG column should contain only numeric values.")
  }

  return(TRUE)
}

#' Generate simulated data for permutation testing
#'
#' Simulates data using observed means but placebo standard deviation
#'
#' @param trts Vector of treatment codes
#' @param nobs_list Named list of observation counts per treatment
#' @param mean_list Named list of means per treatment
#' @param sd_placebo Placebo standard deviation
#' @param n_perms Number of permutations
#' @param seed Random seed
#' @return Data frame with simulated data
gen_sim <- function(trts, nobs_list, mean_list, sd_placebo, n_perms, seed = 123) {
  set.seed(seed)

  # Generate all combinations
  sim_data <- expand.grid(
    NPERMS = 1:n_perms,
    TRT01PN = trts
  )

  # For each row, generate the appropriate number of observations
  result <- do.call(rbind, lapply(1:nrow(sim_data), function(i) {
    perm <- sim_data$NPERMS[i]
    trt <- sim_data$TRT01PN[i]
    n <- nobs_list[[as.character(trt)]]
    m <- mean_list[[as.character(trt)]]

    data.frame(
      NPERMS = perm,
      TRT01PN = trt,
      PT = 1:n,
      M = m,
      CHG = rnorm(n, mean = m, sd = sd_placebo)
    )
  }))

  return(result)
}

#' Calculate sigma using actual data percentiles
#'
#' Uses all observations and their rank-based percentiles
#'
#' @param pbo_df Placebo data frame
#' @param trt_df Treatment data frame
#' @param interval95 Logical, whether to restrict to 3-97 percentile range
#' @return List with ehte_df (sigma results) and ite_df (individual ITE data)
calculate_sigma <- function(pbo_df, trt_df, interval95 = TRUE) {
  # Add NPERMS if not present
  if (!"NPERMS" %in% names(pbo_df)) {
    pbo_df$NPERMS <- 0
  }
  if (!"NPERMS" %in% names(trt_df)) {
    trt_df$NPERMS <- 0
  }

  # Calculate percentiles for placebo (0-100)
  pct_pcb <- pbo_df %>%
    group_by(NPERMS, TRT01PN) %>%
    summarise(
      percentile = 0:100,
      PCB_CHG = sas_percentile(CHG, seq(0, 1, length.out = 101)),
      .groups = "drop"
    )

  # Calculate percentile ranks for treatment observations
  pct_trt <- trt_df %>%
    group_by(NPERMS, TRT01PN) %>%
    mutate(percentile = calculate_percentile(CHG)) %>%
    ungroup()

  # Merge treatment with placebo on percentile
  ite_df <- pct_trt %>%
    left_join(pct_pcb %>% select(NPERMS, percentile, PCB_CHG),
              by = c("NPERMS", "percentile")) %>%
    mutate(ITE = CHG - PCB_CHG)

  # Filter to central 95% if requested
  if (interval95) {
    ite_df <- ite_df %>%
      filter(percentile >= 3 & percentile <= 97)
  }

  # Calculate sigma (SD of ITE) for each permutation and treatment
  ehte_df <- ite_df %>%
    group_by(NPERMS, TRT01PN) %>%
    summarise(sigma = sd(ITE, na.rm = TRUE), .groups = "drop")

  return(list(ehte_df = ehte_df, ite_df = ite_df))
}

#' Calculate sigma using fixed 48 percentiles
#'
#' Uses fixed percentiles from 3 to 97 by 2 (48 total)
#'
#' @param pbo_df Placebo data frame
#' @param trt_df Treatment data frame
#' @return List with ehte_df (sigma results) and ite_df (individual ITE data)
calculate_sigma_48 <- function(pbo_df, trt_df) {
  # Add NPERMS if not present
  if (!"NPERMS" %in% names(pbo_df)) {
    pbo_df$NPERMS <- 0
  }
  if (!"NPERMS" %in% names(trt_df)) {
    trt_df$NPERMS <- 0
  }

  # Fixed percentiles from 3 to 97 by 2
  pctiles <- seq(3, 97, by = 2)
  probs <- pctiles / 100

  # Calculate percentiles for placebo
  pct_pcb <- pbo_df %>%
    group_by(NPERMS, TRT01PN) %>%
    summarise(
      percentile = pctiles,
      PCB_CHG = sas_percentile(CHG, probs),
      .groups = "drop"
    )

  # Calculate percentiles for treatment
  pct_trt <- trt_df %>%
    group_by(NPERMS, TRT01PN) %>%
    summarise(
      percentile = pctiles,
      CHG = sas_percentile(CHG, probs),
      .groups = "drop"
    )

  # Merge and calculate ITE
  ite_df <- pct_trt %>%
    left_join(pct_pcb %>% select(NPERMS, percentile, PCB_CHG),
              by = c("NPERMS", "percentile")) %>%
    mutate(ITE = CHG - PCB_CHG)

  # Calculate sigma
  ehte_df <- ite_df %>%
    group_by(NPERMS, TRT01PN) %>%
    summarise(sigma = sd(ITE, na.rm = TRUE), .groups = "drop")

  return(list(ehte_df = ehte_df, ite_df = ite_df))
}

#' Calculate p-value from permutation distribution
#'
#' @param sim_ehte_df Data frame with simulated sigma values
#' @param obs_ehte_df Data frame with observed sigma values
#' @param sd_placebo Placebo standard deviation
#' @param trt_dict Named vector mapping TRT01PN to treatment names
#' @return Data frame with sigma, eHTE, and p-value for each treatment
cal_pvalue <- function(sim_ehte_df, obs_ehte_df, sd_placebo, trt_dict) {
  # Get observed sigmas
  obs_sigmas <- obs_ehte_df %>%
    select(TRT01PN, sigma) %>%
    rename(obs_sigma = sigma)

  # Merge simulated with observed
  count_df <- sim_ehte_df %>%
    rename(sim_sigma = sigma) %>%
    left_join(obs_sigmas, by = "TRT01PN") %>%
    mutate(ind = ifelse(sim_sigma >= obs_sigma, 1, 0))

  # Calculate p-values
  pvalues <- count_df %>%
    group_by(TRT01PN) %>%
    summarise(pvalue = mean(ind, na.rm = TRUE), .groups = "drop")

  # Build results table
  results <- obs_ehte_df %>%
    left_join(pvalues, by = "TRT01PN") %>%
    mutate(
      eHTE = sigma / sd_placebo,
      TRT01P = trt_dict[as.character(TRT01PN)]
    ) %>%
    select(TRT01PN, TRT01P, sigma, eHTE, pvalue)

  # Print results
  for (i in 1:nrow(results)) {
    cat(sprintf("treatment: %s, sigma: %.5f, eHTE: %.3f, p-value: %.3f\n",
                results$TRT01P[i], results$sigma[i], results$eHTE[i], results$pvalue[i]))
  }

  return(results)
}

#' Main eHTE analysis function
#'
#' Runs the complete eHTE analysis pipeline
#'
#' @param indf Input data frame with TRT01P, TRT01PN, CHG columns
#' @param n_perms Number of permutations for p-value calculation (default 10000)
#' @return List with results tables and individual ITE data
#' @export
eHTE_p <- function(indf, n_perms = 10000) {
  # Validate input
  validate_input(indf)

  # Create treatment dictionary
  trt_dict <- indf %>%
    select(TRT01PN, TRT01P) %>%
    distinct() %>%
    { setNames(.$TRT01P, as.character(.$TRT01PN)) }

  # Print summary statistics
  summary_stats <- indf %>%
    group_by(TRT01PN, TRT01P) %>%
    summarise(
      N = n(),
      Mean = mean(CHG, na.rm = TRUE),
      SD = sd(CHG, na.rm = TRUE),
      Min = min(CHG, na.rm = TRUE),
      Max = max(CHG, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_stats)

  # Separate placebo and treatment arms
  pbo_df <- indf %>% filter(tolower(TRT01P) == "placebo")
  trt_df <- indf %>% filter(tolower(TRT01P) != "placebo")

  # Get placebo statistics
  pbo_stats <- pbo_df %>%
    group_by(TRT01PN) %>%
    summarise(
      nobs = n(),
      m = mean(CHG, na.rm = TRUE),
      s = sd(CHG, na.rm = TRUE),
      .groups = "drop"
    )

  placebo_code <- pbo_stats$TRT01PN[1]
  sd_placebo <- pbo_stats$s[1]

  # Get treatment statistics
  trt_stats <- trt_df %>%
    group_by(TRT01PN) %>%
    summarise(
      nobs = n(),
      m = mean(CHG, na.rm = TRUE),
      .groups = "drop"
    )

  trt_codes <- sort(unique(trt_stats$TRT01PN))

  # Create named lists for simulation
  pbo_nobs_list <- setNames(list(pbo_stats$nobs[1]), as.character(placebo_code))
  pbo_mean_list <- setNames(list(pbo_stats$m[1]), as.character(placebo_code))

  trt_nobs_list <- setNames(as.list(trt_stats$nobs), as.character(trt_stats$TRT01PN))
  trt_mean_list <- setNames(as.list(trt_stats$m), as.character(trt_stats$TRT01PN))

  # Generate simulated data
  cat("\nGenerating simulated data for permutation testing...\n")
  simp_df <- gen_sim(placebo_code, pbo_nobs_list, pbo_mean_list, sd_placebo, n_perms, seed = 456)
  simt_df <- gen_sim(trt_codes, trt_nobs_list, trt_mean_list, sd_placebo, n_perms, seed = 456)

  # Calculate sigma using all records
  cat("\nCalculating sigma using actual data percentiles...\n")
  obs_result <- calculate_sigma(pbo_df, trt_df, interval95 = TRUE)
  sim_result <- calculate_sigma(simp_df, simt_df, interval95 = TRUE)

  # Calculate sigma using 48 fixed percentiles
  cat("\nCalculating sigma using fixed 48 percentiles...\n")
  obs_result_48 <- calculate_sigma_48(pbo_df, trt_df)
  sim_result_48 <- calculate_sigma_48(simp_df, simt_df)

  # Calculate p-values
  cat("\n=== 3-97% percentiles based on actual data ===\n")
  results_all <- cal_pvalue(sim_result$ehte_df, obs_result$ehte_df, sd_placebo, trt_dict)

  cat("\n=== Fixed 48 percentiles ===\n")
  results_48 <- cal_pvalue(sim_result_48$ehte_df, obs_result_48$ehte_df, sd_placebo, trt_dict)

  # Add TRT01P to ITE data for plotting
  ite_df <- obs_result$ite_df
  if ("TRT01P" %in% names(trt_df)) {
    trt_mapping <- trt_df %>% select(TRT01PN, TRT01P) %>% distinct()
    ite_df <- ite_df %>%
      left_join(trt_mapping, by = "TRT01PN")
  }

  return(list(
    ite_df = ite_df,
    results_all = results_all,
    results_48 = results_48,
    summary_stats = summary_stats
  ))
}

#' Generate eHTE plots
#'
#' Creates a 3-panel visualization: histogram, ECDF, and ITE scatter plot
#'
#' @param indf Original input data frame
#' @param ite_df ITE data frame from eHTE_p results
#' @param var_name Variable name for axis labels (default "CHG")
#' @return ggplot object with combined plots
#' @export
eHTE_plot <- function(indf, ite_df, var_name = "CHG") {
  # Color palette (matching Python/SAS)
  color_palette <- c("Placebo" = "#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F")

  # Get treatment names excluding placebo
  trt_names <- unique(indf$TRT01P[tolower(indf$TRT01P) != "placebo"])

  # Create color dictionary
  color_dict <- c("Placebo" = "#0072BD")
  for (i in seq_along(trt_names)) {
    color_dict[trt_names[i]] <- color_palette[i + 1]
  }

  # Plot 1: Histogram with KDE
  p1 <- ggplot(indf, aes(x = CHG, fill = TRT01P, color = TRT01P)) +
    geom_histogram(aes(y = after_stat(count)), bins = 40, alpha = 0.5, position = "dodge") +
    geom_density(aes(y = after_stat(count) * 2), linewidth = 1) +
    scale_fill_manual(values = color_dict) +
    scale_color_manual(values = color_dict) +
    labs(x = var_name, y = "Frequency") +
    theme_minimal() +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    )

  # Plot 2: ECDF
  p2 <- ggplot(indf, aes(x = CHG, color = TRT01P)) +
    stat_ecdf(linewidth = 1) +
    scale_color_manual(values = color_dict) +
    labs(x = var_name, y = "Cumulative Percentile") +
    ylim(0, 1) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank()
    )

  # Plot 3: ITE scatter plot
  p3 <- ggplot(ite_df, aes(x = ITE, y = percentile, color = TRT01P)) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    scale_color_manual(values = color_dict) +
    labs(x = paste("ITE", var_name), y = "Rank") +
    ylim(0, 100) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  # Combine plots
  combined_plot <- gridExtra::grid.arrange(p1, p2, p3, ncol = 3, widths = c(5, 3, 1))

  return(combined_plot)
}

#' Plot histogram with kernel density
#'
#' @param indf Data frame with TRT01PN and CHG columns
#' @param title Plot title (default "Overlapping Histograms of Score")
#' @return ggplot object
#' @export
plot_hist <- function(indf, title = "Overlapping Histograms of Score") {
  color_palette <- c("1" = "#0072BD", "2" = "#D95319", "3" = "#EDB120")

  ggplot(indf, aes(x = CHG, fill = factor(TRT01PN), color = factor(TRT01PN))) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, alpha = 0.5, position = "identity") +
    geom_density(linewidth = 1, alpha = 0) +
    scale_fill_manual(values = color_palette, name = "Treatment") +
    scale_color_manual(values = color_palette, name = "Treatment") +
    labs(title = title, x = "Score", y = "Density") +
    theme_minimal()
}

#' Plot cumulative distribution
#'
#' @param indf Data frame with TRT01PN and CHG columns
#' @param title Plot title (default "Cumulative Score Distribution")
#' @return ggplot object
#' @export
plot_cumulative <- function(indf, title = "Cumulative Score Distribution") {
  color_palette <- c("1" = "#0072BD", "2" = "#D95319", "3" = "#EDB120")

  ggplot(indf, aes(x = CHG, color = factor(TRT01PN))) +
    stat_ecdf(linewidth = 1) +
    scale_color_manual(values = color_palette, name = "Treatment") +
    labs(title = title, x = "Score", y = "Cumulative Percentile") +
    theme_minimal()
}

#' Plot ITE scatter plot
#'
#' @param ite_df ITE data frame from eHTE_p results
#' @param title Plot title (default "Individual Treatment Effects")
#' @return ggplot object
#' @export
plot_ite <- function(ite_df, title = "Individual Treatment Effects") {
  color_palette <- c("1" = "#0072BD", "2" = "#D95319", "3" = "#EDB120")

  # Add placebo reference line (ITE = 0 across all percentiles)
  pbo_ref <- data.frame(
    TRT01PN = 1,
    percentile = seq(3, 97, by = 2),
    ITE = 0
  )

  plot_data <- bind_rows(ite_df, pbo_ref)

  ggplot(plot_data, aes(x = ITE, y = percentile, color = factor(TRT01PN))) +
    geom_point(size = 2, alpha = 0.7) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    scale_color_manual(values = color_palette, name = "Treatment") +
    labs(title = title, x = "ITE", y = "Percentile") +
    ylim(0, 100) +
    theme_minimal()
}

# Print package info when sourced
cat("eHTE R Package loaded successfully.\n")
cat("Use eHTE_p(data) to run the analysis.\n")
cat("Use eHTE_plot(data, ite_df) to generate visualizations.\n")
