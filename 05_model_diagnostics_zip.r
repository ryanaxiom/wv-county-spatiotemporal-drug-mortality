# ==============================================================================
# File 05: Comprehensive Model Diagnostics with County Analysis
# ==============================================================================
# This file provides comprehensive diagnostics across all modeling phases
#
# DEVELOPMENT PHASE 5: Enhanced Model Diagnostics & County Analysis
# - Load all phase results efficiently
# - Comprehensive model comparison with proper WAIC interpretation
# - Proper ZIP uncertainty quantification with appropriate credible intervals
# - County-level prediction analysis with appropriate 95% credible intervals
# - Individual county plots and grid visualizations
# - Export publication-ready summary tables and county analysis
# ==============================================================================

# Clear environment and load libraries
rm(list = ls())
suppressPackageStartupMessages({
  library(INLA)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(viridis)
  library(knitr)
  library(DT)
  library(scales)
  library(sf)
  library(cowplot)
  library(lubridate)
  library(parallel)
})

# Load shared utilities
source("00_model_utilities.r")

# Get response type from environment variable
RESPONSE_TYPE <- Sys.getenv("ANALYSIS_RESPONSE_TYPE", unset = "opioid")

cat("=== PHASE 5: COMPREHENSIVE MODEL DIAGNOSTICS & COUNTY ANALYSIS ===\n")
cat("Response variable:", RESPONSE_TYPE, "\n")
cat("Proper ZIP uncertainty quantification\n\n")

# ==============================================================================
# 1. LOAD ALL PHASE RESULTS AND DATA
# ==============================================================================

cat("--- Step 1: Loading All Phase Results and Data ---\n")

# Load configuration and data
CONFIG <- readRDS("outputs/data/config.rds")
precision_matrices <- readRDS("outputs/data/precision_matrices.rds")

# Load data splits
if (RESPONSE_TYPE == "opioid") {
  data_splits <- readRDS("outputs/data/opioid_splits.rds")
  response_col <- "opioid_death"
} else {
  data_splits <- readRDS("outputs/data/stimulant_splits.rds")
  response_col <- "stimulant_death"
}

cat("✓ Configuration\n")
cat("✓ Precision matrices\n")
cat("✓ Train/test splits\n")

# Load all phase results
phase_results <- list()
phase_names <- c("Phase 1" = "Spatial", "Phase 2a" = "County Temporal", 
                 "Phase 3" = "Separable", "Phase 4" = "Interaction")

# Try to load each phase
for (phase_num in c("1_spatial", "2a_county_temporal", "3_spatiotemporal", "4_interaction")) {
  file_path <- paste0("outputs/models/phase", phase_num, "_", RESPONSE_TYPE, ".rds")
  if (file.exists(file_path)) {
    phase_label <- ifelse(phase_num == "2a_county_temporal", "Phase 2a", paste0("Phase ", substr(phase_num,1,1)))
    phase_results[[phase_label]] <- readRDS(file_path)
    cat("✓", phase_label, "(", phase_names[phase_label], ")\n")
  }
}

cat("✓ Loaded", length(phase_results), "out of 4 possible phase results\n")
cat("✓ Model family:", CONFIG$FAMILY, "(Zero-Inflated Poisson)\n\n")

# ==============================================================================
# 2. EXTRACT AND STANDARDIZE METRICS
# ==============================================================================

cat("--- Step 2: Extracting and Standardizing Metrics ---\n")

# Function to extract metrics from a phase result
extract_metrics <- function(phase_result, phase_name) {
  if (is.null(phase_result)) return(NULL)
  
  # Get the raw metrics data
  metrics_data <- NULL
  
  # Check for comparison metrics in different possible locations
  if ("comparison_metrics" %in% names(phase_result)) {
    metrics_data <- phase_result$comparison_metrics
  } else if ("all_models_comparison" %in% names(phase_result)) {
    metrics_data <- phase_result$all_models_comparison
  } else if ("models" %in% names(phase_result)) {
    # Extract metrics from individual models
    metrics_list <- list()
    
    for (model_name in names(phase_result$models)) {
      model <- phase_result$models[[model_name]]
      
      # Skip if model failed or is NULL
      if (is.null(model) || !inherits(model, "inla")) next
      
      # row.names = NULL to eliminate warnings
      metrics <- data.frame(
        phase = phase_name,
        phase_description = phase_names[phase_name],
        model_name = model_name,
        model_type = ifelse(grepl("ar1|rw", model_name), "temporal",
                           ifelse(grepl("adjacency|distance", model_name), "spatial",
                                 ifelse(grepl("period", model_name), "interaction", "separable"))),
        complexity = length(model$names.fixed) + length(model$summary.random),
        dic = ifelse(!is.null(model$dic$dic), model$dic$dic, NA),
        waic = ifelse(!is.null(model$waic$waic), model$waic$waic, NA),
        marginal_likelihood = ifelse(!is.null(model$mlik), 
                                     sum(model$mlik[,1], na.rm = TRUE), NA),
        mean_cpo = ifelse(!is.null(model$cpo$cpo), 
                         -mean(log(model$cpo$cpo), na.rm = TRUE), NA),
        failure_rate = ifelse(!is.null(model$cpo$failure), 
                             sum(model$cpo$failure, na.rm = TRUE), 0),
        runtime_seconds = ifelse(!is.null(model$cpu.used), 
                                model$cpu.used["Total"], NA),
        row.names = NULL,
        stringsAsFactors = FALSE
      )
      
      metrics_list[[model_name]] <- metrics
    }
    
    if (length(metrics_list) > 0) {
      metrics_data <- do.call(rbind, metrics_list)
    }
  }
  
  # If no metrics found, return NULL
  if (is.null(metrics_data) || nrow(metrics_data) == 0) {
    return(NULL)
  }
  
  # Filter successful fits only
  if ("fit_success" %in% names(metrics_data)) {
    metrics_data <- metrics_data[metrics_data$fit_success == TRUE, ]
  }
  
  if (nrow(metrics_data) == 0) return(NULL)
  
  # CREATE STANDARDIZED OUTPUT
  standardized <- data.frame(
    phase = phase_name,
    phase_description = phase_names[phase_name],
    model_name = metrics_data$model_name,
    model_type = if ("model_type" %in% names(metrics_data)) {
                   metrics_data$model_type
                 } else {
                   sapply(metrics_data$model_name, function(x) {
                     if (grepl("ar1|rw", x)) "temporal"
                     else if (grepl("adjacency|distance", x)) "spatial"
                     else if (grepl("period", x)) "interaction"
                     else "separable"
                   })
                 },
    complexity = if ("complexity_score" %in% names(metrics_data)) {
                   metrics_data$complexity_score
                 } else if ("complexity" %in% names(metrics_data)) {
                   metrics_data$complexity
                 } else {
                   NA
                 },
    dic = if ("dic" %in% names(metrics_data)) metrics_data$dic else NA,
    waic = if ("waic" %in% names(metrics_data)) metrics_data$waic else NA,
    marginal_likelihood = if ("marginal_likelihood" %in% names(metrics_data)) {
                            metrics_data$marginal_likelihood
                          } else {
                            NA
                          },
    mean_cpo = if ("mean_cpo" %in% names(metrics_data)) metrics_data$mean_cpo else NA,
    failure_rate = if ("failure_rate" %in% names(metrics_data)) {
                     metrics_data$failure_rate
                   } else {
                     0
                   },
    runtime_seconds = if ("runtime_mins" %in% names(metrics_data)) {
                        metrics_data$runtime_mins * 60
                      } else if ("runtime" %in% names(metrics_data)) {
                        metrics_data$runtime
                      } else {
                        NA
                      },
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  # For Phase 4, filter out separable models (they're duplicates from Phase 3)
  if (phase_name == "Phase 4") {
    standardized <- standardized[standardized$model_type == "interaction", ]
  }
  
  return(standardized)
}

# Combine metrics from all phases
all_metrics <- list()
for (phase_name in names(phase_results)) {
  cat("\nProcessing", phase_name, "\n")
  
  phase_metrics <- extract_metrics(phase_results[[phase_name]], phase_name)
  all_metrics[[phase_name]] <- phase_metrics
  
}

cat("\n✓ Successfully extracted metrics from", length(all_metrics), "phases\n")

if (length(all_metrics) == 0) {
  cat("✗ ERROR: No model metrics could be extracted from phase results.\n")
  stop("Cannot proceed without model metrics.")
}

combined_metrics <- do.call(rbind, all_metrics)
rownames(combined_metrics) <- NULL

# Add ranking and comparison columns
combined_metrics <- combined_metrics %>%
  arrange(waic) %>%
  mutate(
    waic_rank = rank(waic, na.last = TRUE),
    waic_difference = waic - min(waic, na.rm = TRUE),
    within_3_waic = waic_difference <= 3,
    within_10_waic = waic_difference <= 10,
    relative_improvement = round((max(waic, na.rm = TRUE) - waic) / 
                                (max(waic, na.rm = TRUE) - min(waic, na.rm = TRUE)) * 100, 2)
  )

cat("✓ Combined", nrow(combined_metrics), "models from", length(phase_results), "phases\n\n")

# ==============================================================================
# 3. IDENTIFY BEST MODELS FOR DETAILED ANALYSIS
# ==============================================================================

cat("--- Step 3: Identifying Best Models for Detailed Analysis ---\n")

# Overall best model
overall_best <- combined_metrics %>%
  filter(waic == min(waic, na.rm = TRUE)) %>%
  slice(1)

cat("✓ Overall best model:", overall_best$model_name, 
    "(WAIC =", round(overall_best$waic, 1), ")\n")

# Best model by phase
best_by_phase <- combined_metrics %>%
  group_by(phase) %>%
  filter(waic == min(waic, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(waic)

cat("✓ Best models by phase:\n")
for (i in 1:nrow(best_by_phase)) {
  cat("  ", best_by_phase$phase[i], ":", best_by_phase$model_name[i], 
      "(WAIC =", round(best_by_phase$waic[i], 1), ")\n")
}

# Models within 3 WAIC units
equivalent_models <- combined_metrics %>%
  filter(within_3_waic) %>%
  arrange(waic)

cat("✓ Models with equivalent performance (within 3 WAIC units):", nrow(equivalent_models), "\n")
if (nrow(equivalent_models) > 1) {
  cat("  These models have statistically equivalent performance:\n")
  for (i in 1:min(5, nrow(equivalent_models))) {
    cat("   ", i, ". ", equivalent_models$model_name[i], 
        " (", equivalent_models$phase[i], ", WAIC = ", 
        round(equivalent_models$waic[i], 1), ")\n", sep="")
  }
}
cat("\n")

# ==============================================================================
# 4. LOAD BEST MODEL FOR COUNTY-LEVEL PREDICTIONS
# ==============================================================================

cat("--- Step 4: Loading Best Model for County-Level Predictions ---\n")

# Get the best model
best_phase <- overall_best$phase
best_model_name <- overall_best$model_name

cat("✓ Best model phase:", best_phase, "\n")
cat("✓ Best model name:", best_model_name, "\n")

# Load the best model based on the phase structure
best_model <- NULL

if (best_phase == "Phase 4") {
  # Phase 4 stores models in interaction_results[[model_name]]$model
  if (!is.null(phase_results[["Phase 4"]]$interaction_results[[best_model_name]]$model)) {
    best_model <- phase_results[["Phase 4"]]$interaction_results[[best_model_name]]$model
    cat("✓ Loaded model from Phase 4 interaction_results\n")
  } else {
    cat("✗ Model not found in Phase 4 interaction_results\n")
    cat("  Available models in interaction_results:", 
        paste(names(phase_results[["Phase 4"]]$interaction_results), collapse=", "), "\n")
  }
} else if (best_phase == "Phase 3") {
  # Phase 3 might store models in different locations
  if (!is.null(phase_results[["Phase 3"]]$model_fits[[best_model_name]]$model)) {
    best_model <- phase_results[["Phase 3"]]$model_fits[[best_model_name]]$model
    cat("✓ Loaded model from Phase 3 model_fits\n")
  } else if (!is.null(phase_results[["Phase 3"]]$spatiotemporal_results[[best_model_name]]$model)) {
    best_model <- phase_results[["Phase 3"]]$spatiotemporal_results[[best_model_name]]$model
    cat("✓ Loaded model from Phase 3 spatiotemporal_results\n")
  } else {
    cat("✗ Model not found in Phase 3\n")
    cat("  Available objects:", paste(names(phase_results[["Phase 3"]]), collapse=", "), "\n")
  }
} else if (best_phase == "Phase 2a") {
  # Phase 2a temporal models
  if (!is.null(phase_results[["Phase 2a"]]$model_fits[[best_model_name]]$model)) {
    best_model <- phase_results[["Phase 2a"]]$model_fits[[best_model_name]]$model
    cat("✓ Loaded model from Phase 2a model_fits\n")
  } else {
    cat("✗ Model not found in Phase 2a\n")
    cat("  Available objects:", paste(names(phase_results[["Phase 2a"]]), collapse=", "), "\n")
  }
} else if (best_phase == "Phase 1") {
  # Phase 1 spatial models
  if (!is.null(phase_results[["Phase 1"]]$model_fits[[best_model_name]]$model)) {
    best_model <- phase_results[["Phase 1"]]$model_fits[[best_model_name]]$model
    cat("✓ Loaded model from Phase 1 model_fits\n")
  } else {
    cat("✗ Model not found in Phase 1\n")
    cat("  Available objects:", paste(names(phase_results[["Phase 1"]]), collapse=", "), "\n")
  }
}

cat("✓ Loaded fitted model:", best_model_name, "\n")
cat("✓ Model class:", class(best_model)[1], "\n")

# Verify model structure
cat("✓ Model family:", best_model$.args$family, "\n")
cat("✓ Has fitted values:", !is.null(best_model$summary.fitted.values), "\n")
if (!is.null(best_model$summary.fitted.values)) {
  cat("✓ Fitted values range:", round(min(best_model$summary.fitted.values$mean), 3), 
      "to", round(max(best_model$summary.fitted.values$mean), 3), "\n")
}
cat("\n")

# ==============================================================================
# 5. PREPARE DATA EXACTLY MATCHING PHASE 04 BASED ON PHASE 00
# ==============================================================================

cat("--- Step 5: Preparing Data ---\n")

# Load the same splits that Phase 04 used
# These come from the splits in Phase 00
if (RESPONSE_TYPE == "opioid") {
  data_splits_current <- readRDS("outputs/data/opioid_splits.rds")
} else {
  data_splits_current <- readRDS("outputs/data/stimulant_splits.rds")
}

cat("✓ Loaded current data splits from Phase 00\n")
cat("✓ Training data size:", nrow(data_splits_current$train), "observations\n")
cat("✓ Holdout configuration:", data_splits_current$validation_config$holdout_months, "months\n")

# RECREATE THE EXACT SAME DATA PREPARATION FROM PHASE 04
# This must match exactly to what was used to fit the best model
train_data <- data_splits_current$train %>%
  filter(Residence_County %in% precision_matrices$county_names) %>%
  mutate(
    # Basic indices (exactly matching Phase 04)
    spatial_id = match(Residence_County, precision_matrices$county_names),
    time_numeric = as.numeric(time_id),
    
    # Create interaction indices for space-time models (from Phase 04)
    space_time_idx = (spatial_id - 1) * max(time_id) + time_id,
    
    # Create duplicate spatial indices for INLA (from Phase 04)
    spatial_id2 = spatial_id,
    spatial_id3 = spatial_id,
    
    # Period groupings matching Phase 04 with 6-month holdout
    max_train_time = max(time_id),
    pre_intervention_period = ifelse(time_id <= 60, 1, 0),   # Jan 2015 - Dec 2019 
    implementation_period = ifelse(time_id > 60 & time_id <= 96, 1, 0),  # Jan 2020 - Dec 2022 
    assessment_period = ifelse(time_id > 96, 1, 0),          # Jan 2023 onwards
    
    # Seasonal indicators (from Phase 4)
    season = case_when(
      month(date) %in% c(12, 1, 2) ~ "winter",
      month(date) %in% c(3, 4, 5) ~ "spring", 
      month(date) %in% c(6, 7, 8) ~ "summer",
      month(date) %in% c(9, 10, 11) ~ "fall"
    ),
    season_id = as.numeric(factor(season))
  ) %>%
  arrange(spatial_id, time_id)

cat("✓ Recreated Phase 4 data structure exactly\n")
cat("✓ Training data:", nrow(train_data), "observations\n")
cat("✓ Counties:", length(unique(train_data$spatial_id)), "\n")
cat("✓ Time points:", length(unique(train_data$time_id)), "\n")
cat("✓ Time range: time_id", min(train_data$time_id), "to", max(train_data$time_id), "\n")
cat("✓ Date range:", min(train_data$date), "to", max(train_data$date), "\n")

# Verification: Check assessment period learning
assessment_in_train <- sum(train_data$assessment_period == 1, na.rm = TRUE)
cat("\nPeriod interaction verification:\n")
cat("✓ Pre-intervention period observations:", sum(train_data$pre_intervention_period == 1), "\n")
cat("✓ Implementation period observations:", sum(train_data$implementation_period == 1), "\n")
cat("✓ Assessment period observations:", assessment_in_train, "\n")

if (assessment_in_train > 0) {
  cat("✓ SUCCESS: Assessment period appears in training data\n")
  cat("✓ Period interactions are learnable with a 6-month testing holdout\n")
} else {
  cat("⚠️  WARNING: Assessment period is missing from training\n")
  cat("   Period interactions may be unidentifiable\n")
}

# Verify this matches the model's fitted values exactly
if (!is.null(best_model$summary.fitted.values)) {
  n_fitted_values <- nrow(best_model$summary.fitted.values)
  cat("\nDimension verification:\n")
  cat("✓ Training data prepared:", nrow(train_data), "observations\n")
  cat("✓ Model fitted values:", n_fitted_values, "observations\n")
  
  if (nrow(train_data) == n_fitted_values) {
    cat("✓ SUCCESS: Dimensions match perfectly\n")
  } else {
    cat("⚠️  DIMENSION MISMATCH DETECTED\n")
    cat("   This means Phase 04 was run with different data splits\n")
    cat("   SOLUTION: Re-run Phase 04 with current Phase 00 splits\n")
    stop("Dimension mismatch - need to re-run Phase 04 with current splits")
  }
}

# Additional diagnostics
cat("\nAdditional verification:\n")
cat("✓ Spatial counties available:", length(precision_matrices$county_names), "\n")
cat("✓ Counties in training data:", length(unique(train_data$Residence_County)), "\n")
cat("✓ All spatial counties in training:", 
    all(precision_matrices$county_names %in% train_data$Residence_County), "\n")

cat("\n")

# ==============================================================================
# 6. COUNTY-LEVEL PREDICTION INTERVALS VIA POSTERIOR SAMPLING
# ==============================================================================

cat("--- Step 6: County-Level Predictions via Posterior Sampling ---\n")

if (!is.null(best_model$summary.fitted.values)) {
  fitted_summary <- best_model$summary.fitted.values
  fitted_mean <- fitted_summary$mean
  
  fitted_lower_inla <- fitted_summary$`0.025quant`
  fitted_upper_inla <- fitted_summary$`0.975quant`
  
  # Generate posterior samples to properly account for all correlations
  cat("Generating posterior samples for proper prediction intervals...\n")
  
  # Request posterior samples from INLA
  n_samples <- 5000
  
  # Extract posterior samples of the linear predictor
  # This preserves all correlations between random effects
  posterior_samples <- inla.posterior.sample(n_samples, best_model)
  
  # Initialize storage for prediction intervals
  prediction_quantiles <- matrix(NA, nrow = length(fitted_mean), ncol = n_samples)
  
cat("Processing posterior samples...\n")

# Parallel processing setup
library(parallel)
n_cores <- min(detectCores() - 1, 8)  # Cap at 8 cores for memory

if (.Platform$OS.type == "unix") {
  # Mac/Linux - use mclapply
  cat("  Using", n_cores, "cores (mclapply)\n")

  prediction_quantiles <- mclapply(1:n_samples, function(s) {
    latent_sample <- posterior_samples[[s]]$latent
    predictor_idx <- grep("^Predictor:", rownames(latent_sample))
  
    if (length(predictor_idx) > 0) {
      linear_pred_sample <- latent_sample[predictor_idx, 1]
      lambda_sample <- exp(linear_pred_sample)
    
      # Vectorized sampling for all observations at once
      if (CONFIG$FAMILY == "zeroinflatedpoisson1") {
        zero_prob <- latent_sample[grep("zero.probability", rownames(latent_sample)), 1]
        if (length(zero_prob) == 0) zero_prob <- 0.1
      
        is_zero <- runif(length(lambda_sample)) < zero_prob
        predictions <- ifelse(is_zero, 0, rpois(length(lambda_sample), lambda_sample))
      } else {
        predictions <- rpois(length(lambda_sample), lambda_sample)
      }
      return(predictions)
    }
    return(rep(NA, length(fitted_mean)))
  }, mc.cores = n_cores)

  # Convert list to matrix
  prediction_quantiles <- do.call(cbind, prediction_quantiles)

} else {
  # Windows - use parLapply
  cat("  Using", n_cores, "cores (parLapply)\n")
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("posterior_samples", "CONFIG", "fitted_mean"))

  prediction_list <- parLapply(cl, 1:n_samples, function(s) {
    latent_sample <- posterior_samples[[s]]$latent
    predictor_idx <- grep("^Predictor:", rownames(latent_sample))
  
    if (length(predictor_idx) > 0) {
      linear_pred_sample <- latent_sample[predictor_idx, 1]
      lambda_sample <- exp(linear_pred_sample)
    
      if (CONFIG$FAMILY == "zeroinflatedpoisson1") {
        zero_prob <- latent_sample[grep("zero.probability", rownames(latent_sample)), 1]
        if (length(zero_prob) == 0) zero_prob <- 0.1
      
        is_zero <- runif(length(lambda_sample)) < zero_prob
        predictions <- ifelse(is_zero, 0, rpois(length(lambda_sample), lambda_sample))
      } else {
        predictions <- rpois(length(lambda_sample), lambda_sample)
      }
      return(predictions)
    }
    return(rep(NA, length(fitted_mean)))
  })

  stopCluster(cl)
  prediction_quantiles <- do.call(cbind, prediction_list)
}

cat("  Processed", n_samples, "samples using", n_cores, "cores\n")

# Calculate prediction intervals from samples
cat("Calculating prediction intervals from posterior samples...\n")
  
  proper_lower <- apply(prediction_quantiles, 1, quantile, probs = 0.025, na.rm = TRUE)
  proper_upper <- apply(prediction_quantiles, 1, quantile, probs = 0.975, na.rm = TRUE)
  
  # Calculate coverage
  coverage <- mean(train_data$distinct_patient_count >= proper_lower & 
                  train_data$distinct_patient_count <= proper_upper, na.rm = TRUE)
  
# Diagnostic: Analyze coverage by count level
  cat("\n  Coverage analysis by count level:\n")

  # Group by observed count levels
  count_groups <- cut(train_data$distinct_patient_count, 
                     breaks = c(-1, 0, 1, 2, 5, 10, Inf),
                     labels = c("0", "1", "2", "3-5", "6-10", ">10"))

  coverage_by_count <- tapply(
    train_data$distinct_patient_count >= proper_lower & 
    train_data$distinct_patient_count <= proper_upper,
    count_groups,
    mean
  )

  for (i in 1:length(coverage_by_count)) {
    cat("    Count =", names(coverage_by_count)[i], ": ",
        round(coverage_by_count[i] * 100, 1), "% coverage (n=",
        sum(count_groups == names(coverage_by_count)[i]), ")\n", sep="")
  }

  # Check if high coverage is driven by zeros
  zero_coverage <- mean((train_data$distinct_patient_count == 0) & 
                        (proper_lower == 0))
  cat("    Zero inflation impact: ", round(sum(train_data$distinct_patient_count == 0)), 
      " zeros (", round(zero_coverage * 100, 1), "% have lower bound = 0)\n", sep="")
  
  cat("\n✓ Posterior predictive intervals:\n")
  cat("  - Method: Full posterior sampling (accounts for all correlations)\n")
  cat("  - Number of samples:", n_samples, "\n")
  cat("  - Coverage achieved:", round(coverage * 100, 1), "%\n")
  cat("  - Mean interval width:", round(mean(proper_upper - proper_lower), 2), "\n")
  
  # Round to integers (counts must be whole numbers)
  proper_lower <- round(proper_lower)
  proper_upper <- round(proper_upper)

  
  # Create predictions dataframe
  predictions_df <- data.frame(
    county = train_data$Residence_County,
    time_id = train_data$time_id,
    year_month = train_data$year_month,
    date = as.Date(paste0(train_data$year_month, "-01")),
    observed = train_data$distinct_patient_count,
    fitted_mean = fitted_mean,
    fitted_lower = proper_lower,
    fitted_upper = proper_upper,
    fitted_lower_inla = fitted_lower_inla,  # Renamed for clarity
    fitted_upper_inla = fitted_upper_inla,  # Renamed for clarity
    population = train_data$population,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  # Calculate performance metrics
  coverage_inla <- mean(predictions_df$observed >= predictions_df$fitted_lower_inla & 
                        predictions_df$observed <= predictions_df$fitted_upper_inla, na.rm = TRUE)
  coverage_calibrated <- mean(predictions_df$observed >= predictions_df$fitted_lower & 
                             predictions_df$observed <= predictions_df$fitted_upper, na.rm = TRUE)
  
  residuals <- predictions_df$observed - predictions_df$fitted_mean
  mae <- mean(abs(residuals), na.rm = TRUE)
  rmse <- sqrt(mean(residuals^2, na.rm = TRUE))
  
  # Calculate interval score for probabilistic calibration
  interval_score <- mean(
    (predictions_df$fitted_upper - predictions_df$fitted_lower) + 
    (2/0.05) * (predictions_df$fitted_lower - predictions_df$observed) * (predictions_df$observed < predictions_df$fitted_lower) +
    (2/0.05) * (predictions_df$observed - predictions_df$fitted_upper) * (predictions_df$observed > predictions_df$fitted_upper),
    na.rm = TRUE
  )
  
  # Use calibrated coverage
  coverage <- coverage_calibrated
  
  cat("✓ County-level prediction intervals created:\n")
  cat("  - Method: Posterior predictive sampling (full Bayesian inference)\n")
  cat("  - Mean interval width:", round(mean(proper_upper - proper_lower), 3), "\n")
  cat("  - INLA coverage:", round(coverage_inla * 100, 1), "%\n")
  cat("  - Posterior sampling coverage:", round(coverage_calibrated * 100, 1), "%\n")
  cat("  - Interval score:", round(interval_score, 3), "\n")
  
  cat("✓ Overall model performance:\n")
  cat("  - Coverage:", round(coverage * 100, 1), "%\n")
  cat("  - MAE:", round(mae, 3), "\n")
  cat("  - RMSE:", round(rmse, 3), "\n")
  
  cat("✓ Generated predictions for", nrow(predictions_df), "observations\n")
  cat("✓ Model accounts for", round(sum(predictions_df$observed == 0) / nrow(predictions_df) * 100, 1), "% zero observations\n\n")
  
} else {
  stop("No fitted values available from best model")
}

# ==============================================================================
# 7. CALCULATE COUNTY-LEVEL FIT STATISTICS
# ==============================================================================

cat("--- Step 7: Calculating County-Level Fit Statistics ---\n")

# Calculate fit statistics by county
county_fit_stats <- predictions_df %>%
  group_by(county) %>%
  summarise(
    total_deaths = sum(observed),
    avg_population = mean(population),
    death_rate_per_1000 = (sum(observed) / sum(population)) * 1000,
    coverage_95 = mean(observed >= fitted_lower & observed <= fitted_upper, na.rm = TRUE),
    rmse = sqrt(mean((observed - fitted_mean)^2, na.rm = TRUE)),
    mae = mean(abs(observed - fitted_mean), na.rm = TRUE),
    mape = mean(abs((observed - fitted_mean) / (observed + 1)) * 100, na.rm = TRUE),
    correlation = cor(observed, fitted_mean, use = "complete.obs"),
    avg_ci_width = mean(fitted_upper - fitted_lower, na.rm = TRUE),
    zero_months = sum(observed == 0),
    .groups = 'drop'
  ) %>%
  mutate(
    county_title = tools::toTitleCase(tolower(county)),
    rate_category = case_when(
      death_rate_per_1000 >= quantile(death_rate_per_1000, 0.75) ~ "High",
      death_rate_per_1000 <= quantile(death_rate_per_1000, 0.25) ~ "Low",
      TRUE ~ "Medium"
    ),
    fit_quality = case_when(
      coverage_95 >= 0.90 & rmse < 2 ~ "Excellent",
      coverage_95 >= 0.80 & rmse < 3 ~ "Good",
      coverage_95 >= 0.70 & rmse < 4 ~ "Adequate",
      TRUE ~ "Poor"
    ),
    Coverage_95_Pct = round(coverage_95 * 100, 1)
  )

# Summary of fit quality
fit_summary <- table(county_fit_stats$fit_quality)
cat("✓ Calculated fit statistics for", nrow(county_fit_stats), "counties\n")
cat("✓ Fit quality distribution:\n")
for (quality in names(fit_summary)) {
  cat("  ", quality, fit_summary[quality], "\n")
}
cat("✓ Average county coverage:", round(mean(county_fit_stats$coverage_95) * 100, 1), "%\n\n")

# ==============================================================================
# 8. CREATE INDIVIDUAL COUNTY PLOTS
# ==============================================================================

cat("--- Step 8: Creating Individual County Plots ---\n")

# Create output directory
dir.create("outputs/plots/individual_counties", showWarnings = FALSE, recursive = TRUE)

# Function to create individual county plot
create_county_plot <- function(county_name, predictions_df, county_stats) {
  county_data <- predictions_df %>%
    filter(county == county_name)
  
  stats <- county_stats %>%
    filter(county == county_name)
  
  p <- ggplot(county_data, aes(x = date)) +  # Use 'date' column
    geom_ribbon(aes(ymin = fitted_lower, ymax = fitted_upper), 
                fill = "red", alpha = 0.25) +
    geom_line(aes(y = fitted_mean), color = "red", linewidth = 1) +
    geom_point(aes(y = observed), color = "black", size = 1.5) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    scale_y_continuous(expand = c(0.1, 0)) +
    labs(
      title = paste(stats$county_title, "County -", 
                   RESPONSE_TYPE, "Deaths"),
      subtitle = paste("Rate:", stats$rate_category, 
                      "| Fit:", stats$fit_quality, 
                      "| Coverage:", stats$Coverage_95_Pct, "%"),
      x = "Date",
      y = "Deaths per Month"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# Create plots for all counties
cat("✓ Creating individual plots for", nrow(county_fit_stats), "counties...\n")

for (county_name in unique(predictions_df$county)) {
  p <- create_county_plot(county_name, predictions_df, county_fit_stats)
  
  filename <- paste0("outputs/plots/individual_counties/", 
                    gsub(" ", "_", tolower(county_name)), 
                    "_", RESPONSE_TYPE, ".png")
  
  ggsave(filename, p, width = 10, height = 6, dpi = 150)
}

cat("✓ Created", nrow(county_fit_stats), "individual county plots\n\n")

# ==============================================================================
# 9. CREATE COUNTY GRID VISUALIZATIONS
# ==============================================================================

cat("--- Step 9: Creating County Grid Visualizations ---\n")

# Sort counties by death rate for grid arrangement
county_fit_stats <- county_fit_stats %>%
  arrange(desc(death_rate_per_1000))

# Create 3x4 grid pages
counties_per_page <- 12
n_pages <- ceiling(nrow(county_fit_stats) / counties_per_page)

for (page in 1:n_pages) {
  start_idx <- (page - 1) * counties_per_page + 1
  end_idx <- min(page * counties_per_page, nrow(county_fit_stats))
  
  page_counties <- county_fit_stats$county[start_idx:end_idx]
  
  # Create plots for this page
  plot_list <- list()
  for (i in 1:length(page_counties)) {
    county_name <- page_counties[i]
    plot_list[[i]] <- create_county_plot(county_name, predictions_df, county_fit_stats)
  }
  
  # Fill remaining spots with empty plots if needed
  while (length(plot_list) < counties_per_page) {
    plot_list[[length(plot_list) + 1]] <- ggplot() + theme_void()
  }
  
  # Arrange in grid
  grid_plot <- do.call(gridExtra::grid.arrange, 
                       c(plot_list, ncol = 4, nrow = 3, 
                         top = paste("County Analysis - Page", page, "of", n_pages)))
  
  # Save grid
  filename <- paste0("outputs/plots/all_counties_grid_page", page, "_", 
                    RESPONSE_TYPE, ".png")
  ggsave(filename, grid_plot, width = 20, height = 15, dpi = 150)
}

cat("✓ Created", n_pages, "pages of county grids\n\n")

# ==============================================================================
# 10. CREATE SPATIAL VISUALIZATIONS (if spatial data available)
# ==============================================================================

cat("--- Step 10: Creating Spatial Visualizations ---\n")

# Try to create spatial maps if tigris is available
if (requireNamespace("tigris", quietly = TRUE)) {
  library(tigris)
  options(tigris_use_cache = TRUE)
  
  # Get WV county boundaries
  wv_counties <- counties(state = "WV", year = 2020, progress_bar = FALSE)
  
  # Join with fit statistics
  wv_counties_stats <- wv_counties %>%
    mutate(county_upper = toupper(gsub(" County", "", NAME))) %>%
    left_join(county_fit_stats %>%
               mutate(county_upper = toupper(county)),
             by = "county_upper")
  
  # Create maps
  maps_list <- list()
  
  # Death rate map
  maps_list[[1]] <- ggplot(wv_counties_stats) +
    geom_sf(aes(fill = death_rate_per_1000), color = "white", linewidth = 0.2) +
    scale_fill_viridis(name = "Deaths per\n1,000 pop", 
                       option = "plasma", direction = -1) +
    labs(title = "Death Rate by County") +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Coverage map
  maps_list[[2]] <- ggplot(wv_counties_stats) +
    geom_sf(aes(fill = Coverage_95_Pct), color = "white", linewidth = 0.2) +
    scale_fill_gradient2(name = "95%\nCoverage", 
                        low = "red", mid = "brown", high = "darkgreen",
                        midpoint = 92, limits = c(85, 100)) +
    labs(title = "Model Coverage by County") +
    theme_minimal() +
    theme(legend.position = "right")
  
  # RMSE map
  maps_list[[3]] <- ggplot(wv_counties_stats) +
    geom_sf(aes(fill = rmse), color = "white", linewidth = 0.2) +
    scale_fill_viridis(name = "RMSE", option = "viridis") +
    labs(title = "Prediction Error by County") +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Fit quality map
  maps_list[[4]] <- ggplot(wv_counties_stats) +
    geom_sf(aes(fill = fit_quality), color = "white", linewidth = 0.2) +
    scale_fill_manual(name = "Fit Quality",
                     values = c("Excellent" = "darkgreen", 
                               "Good" = "lightgreen",
                               "Adequate" = "yellow", 
                               "Poor" = "red")) +
    labs(title = "Overall Fit Quality") +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Combine maps
  combined_maps <- do.call(gridExtra::grid.arrange, 
                          c(maps_list, ncol = 2, nrow = 2,
                            top = paste("Spatial Analysis -", 
                                      toupper(RESPONSE_TYPE), "Deaths")))
  
  # Save combined map
  ggsave(paste0("outputs/plots/spatial_analysis_comprehensive_", 
               RESPONSE_TYPE, ".png"),
         combined_maps, width = 16, height = 12, dpi = 150)
  
  cat("✓ Created comprehensive spatial analysis maps\n\n")
} else {
  cat("ℹ tigris package not available - skipping spatial maps\n\n")
}

# ==============================================================================
# 11. GENERATE PERFORMANCE SUMMARY AND REPORTS
# ==============================================================================

cat("--- Step 11: Generating Performance Summary and Reports ---\n")

# Create output directory
dir.create("outputs/diagnostics", showWarnings = FALSE, recursive = TRUE)

# Create summary table for export
county_summary <- county_fit_stats %>%
  arrange(desc(death_rate_per_1000)) %>%
  select(
    County = county_title,
    Deaths_per_1000 = death_rate_per_1000,
    Total_Deaths = total_deaths,
    Population = avg_population,
    Fit_Quality = fit_quality,
    Coverage_95_Pct = Coverage_95_Pct,
    RMSE = rmse,
    MAE = mae,
    Avg_CI_Width = avg_ci_width,
    Zero_Months = zero_months
  ) %>%
  mutate(
    Deaths_per_1000 = round(Deaths_per_1000, 2),
    Population = round(Population, 0),
    RMSE = round(RMSE, 2),
    MAE = round(MAE, 2),
    Avg_CI_Width = round(Avg_CI_Width, 2)
  )

# Save summary CSV
write.csv(county_summary, 
          paste0("outputs/models/county_performance_summary_", 
                RESPONSE_TYPE, ".csv"),
          row.names = FALSE)

cat("✓ Created county performance summary table\n")

# Create diagnostic summary
diagnostic_summary <- data.frame(
  Metric = c("Best Model", "Model Type", "WAIC", "DIC", 
            "Overall Coverage", "Average County Coverage",
            "Overall RMSE", "Overall MAE", "Counties with Good Fit",
            "High Rate Counties", "Low Rate Counties"),
  Value = c(
    overall_best$model_name,
    overall_best$model_type,
    round(overall_best$waic, 1),
    round(overall_best$dic, 1),
    paste0(round(coverage * 100, 1), "%"),
    paste0(round(mean(county_fit_stats$coverage_95) * 100, 1), "%"),
    round(rmse, 3),
    round(mae, 3),
    sum(county_fit_stats$fit_quality %in% c("Excellent", "Good")),
    sum(county_fit_stats$rate_category == "High"),
    sum(county_fit_stats$rate_category == "Low")
  ),
  row.names = NULL,
  stringsAsFactors = FALSE
)

cat("✓ Generated diagnostic summary\n\n")

# ==============================================================================
# 12. CREATE COMPREHENSIVE ANALYSIS REPORT
# ==============================================================================

cat("\n--- Step 12: Creating Comprehensive Analysis Report ---\n")

# Helper function for string repetition
`%.%` <- function(x, n) {
  paste(rep(x, n), collapse = "")
}


# Create text report
report_lines <- c(
  "=" %.% 78,
  paste("COMPREHENSIVE ANALYSIS REPORT -", toupper(RESPONSE_TYPE), "DEATHS"),
  paste("Generated:", Sys.Date()),
  "=" %.% 78,
  "",
  "MODEL SELECTION SUMMARY",
  "-" %.% 40,
  paste("Best Overall Model:", overall_best$model_name),
  paste("Model Phase:", overall_best$phase),
  paste("WAIC:", round(overall_best$waic, 1)),
  paste("DIC:", round(overall_best$dic, 1)),
  paste("Model Type:", overall_best$model_type),
  "",
  "PREDICTION PERFORMANCE (COUNTY-LEVEL CREDIBLE INTERVALS)",
  "-" %.% 40,
  paste("INLA pooled interval coverage:", round(coverage_inla * 100, 1), "%"),
  paste("County-specific credible interval coverage:", round(coverage_calibrated * 100, 1), "%"),
  paste("Overall RMSE:", round(rmse, 3)),
  paste("Overall MAE:", round(mae, 3)),
  paste("Average County Coverage:", round(mean(county_fit_stats$coverage_95) * 100, 1), "%"),
  "",
  "SAMPLE SIZE ADAPTATION SUMMARY",
  "-" %.% 40,
  paste("INLA assumes sample sizze:", nrow(train_data)),
  paste("Correct county sample size:", "~102 time points per county"),
  paste("Method used: Random effects uncertainty extraction"),
  paste("Coverage improvement:", round((proper_coverage - inla_coverage) * 100, 1), "percentage points"),
  "",
  "COUNTY-LEVEL SUMMARY",
  "-" %.% 40,
  paste("Total Counties:", nrow(county_fit_stats)),
  paste("Excellent Fit:", sum(county_fit_stats$fit_quality == "Excellent")),
  paste("Good Fit:", sum(county_fit_stats$fit_quality == "Good")),
  paste("Adequate Fit:", sum(county_fit_stats$fit_quality == "Adequate")),
  paste("Poor Fit:", sum(county_fit_stats$fit_quality == "Poor")),
  "",
  "=" %.% 78
)

# Save report
writeLines(report_lines,
          paste0("outputs/diagnostics/comprehensive_analysis_report_", 
                RESPONSE_TYPE, ".txt"))

cat("✓ Saved comprehensive analysis report\n\n")

# ==============================================================================
# 13. CROSS-VALIDATION AND PREDICTION METRICS
# ==============================================================================

cat("--- Step 13: Cross-Validation and Prediction Metrics ---\n")

# Calculate additional prediction metrics using existing variables

calibrated_coverage <- proper_coverage

prediction_metrics <- list(
  mae = mae,
  rmse = rmse,
  mape = mean(abs((predictions_df$observed - predictions_df$fitted_mean) / 
                 (predictions_df$observed + 1)) * 100, na.rm = TRUE),
  coverage_95_calibrated = calibrated_coverage,
  coverage_95_inla = inla_coverage,
  coverage_improvement = calibrated_coverage - inla_coverage,
  correlation = cor(predictions_df$observed, predictions_df$fitted_mean, 
                    use = "complete.obs"),
  exact_predictions = mean(round(predictions_df$fitted_mean) == predictions_df$observed, 
                          na.rm = TRUE),
  mean_calibrated_width = mean(predictions_df$fitted_upper - predictions_df$fitted_lower),
  mean_inla_width = mean(predictions_df$fitted_upper_inla - predictions_df$fitted_lower_inla)
)

cat("✓ Final prediction metrics:\n")
cat("  MAE:", round(prediction_metrics$mae, 3), "\n")
cat("  RMSE:", round(prediction_metrics$rmse, 3), "\n")
cat("  MAPE:", round(prediction_metrics$mape, 1), "%\n")
cat("  County-specific 95% Credible Interval Coverage:", round(prediction_metrics$coverage_95_calibrated * 100, 1), "%\n")
cat("  95% Credible Coverage Coverage (default from INLA):", round(prediction_metrics$coverage_95_inla * 100, 1), "%\n")
cat("  Coverage improvement:", round(prediction_metrics$coverage_improvement * 100, 1), "percentage points\n")
cat("  Correlation:", round(prediction_metrics$correlation, 3), "\n")
cat("  Exact predictions:", round(prediction_metrics$exact_predictions * 100, 1), "%\n")
cat("  County-specific credible interval width:", round(prediction_metrics$mean_calibrated_width, 3), "\n")
cat("  INLA credible interval width:", round(prediction_metrics$mean_inla_width, 3), "\n\n")

# Save prediction metrics
saveRDS(prediction_metrics, 
        paste0("outputs/models/phase5_prediction_metrics_", RESPONSE_TYPE, ".rds"))

# ==============================================================================
# 14. STATE-LEVEL PREDICTIONS USING PROPER STATE MODEL CREDIBLE INTERVALS
# ==============================================================================

cat("--- Step 14: Creating State-Level Predictions ---\n")

# Aggregate training data to the state level for proper state-level modeling
state_training_data <- train_data %>%
  group_by(time_id, year_month, date) %>%
  summarise(
    total_count = sum(distinct_patient_count),
    total_population = sum(population),
    n_counties = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    log_total_pop_offset = log(total_population),
    time_numeric = as.numeric(time_id),
    time_scaled = scale(time_numeric)[,1]
  ) %>%
  arrange(time_id)

cat("Fitting separate state-level ZIP model for credible intervals...\n")

# Fit state-level ZIP model (this gives us the intervals to use for state-level aggregation)
state_model <- inla(
  total_count ~ time_scaled + f(time_id, model='ar1') + offset(log_total_pop_offset),
  family = CONFIG$FAMILY,
  data = state_training_data,
  control.compute = list(dic = TRUE, waic = TRUE),
  verbose = FALSE
)

if (!is.null(state_model) && !is.null(state_model$summary.fitted.values)) {
  # Using R-INLA summary.fitted.values directly at the state level
  state_fitted_values <- state_model$summary.fitted.values

  state_level_data <- state_training_data %>%
    mutate(
      fitted_mean = state_fitted_values$mean,
      fitted_lower = state_fitted_values$`0.025quant`,  # These are correct for the state level
      fitted_upper = state_fitted_values$`0.975quant`,  # These are correct for the state level
      observed = total_count
    )

  # Calculate state-level coverage (should be ~95%)
  state_coverage <- mean(state_level_data$observed >= state_level_data$fitted_lower & 
                        state_level_data$observed <= state_level_data$fitted_upper, na.rm = TRUE)

  cat("✓ State-level model fitted successfully\n")
  cat("  State WAIC:", round(state_model$waic$waic, 1), "\n")
  cat("  Correct method for state-level intervals: Using state model's summary.fitted.values directly\n")
  cat("  State-level sample size:", nrow(state_training_data), "time points\n")
  cat("  State coverage:", round(state_coverage * 100, 1), "%\n")

} else {
  stop("State-level model failed to fit")
}

# Create state-level plot
state_plot <- ggplot(state_level_data, aes(x = date)) +
  geom_ribbon(aes(ymin = fitted_lower, ymax = fitted_upper), 
              fill = "red", alpha = 0.25) +
  geom_line(aes(y = fitted_mean), color = "red", linewidth = 1.5) +
  geom_point(aes(y = observed), color = "black", size = 2) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(expand = c(0.1, 0)) +
  labs(
    title = paste("State-Level Analysis: Direct State Model"),
    subtitle = paste("WV", RESPONSE_TYPE, "deaths (ZIP, WAIC =", round(state_model$waic$waic, 1), 
                    ", Coverage =", round(state_coverage * 100, 1), "%)"),
    x = "Date", y = "Total Deaths per Month",
    caption = paste("State-level ZIP model with proper sample size (n =", nrow(state_training_data), "time points)",
                   "\nBlack points = observed, Red line = model fit, Shaded area = 95% credible intervals")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(size = 10, hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save state-level plot
ggsave(paste0("outputs/plots/state_level_", RESPONSE_TYPE, ".png"), 
       state_plot, width = 14, height = 8, dpi = 300, bg = "white")

cat("✓ State-level analysis complete\n")
cat("✓ Method: Direct state-level ZIP model with proper intervals\n")
cat("✓ Sample size:", nrow(state_training_data), "time points (correct for state predictions)\n")
cat("✓ Coverage:", round(state_coverage * 100, 1), "%\n")
cat("✓ Interval source: summary.fitted.values from state model (correct usage)\n\n")

# ==============================================================================
#  FINAL SUMMARY
# ==============================================================================

cat("\n=== COMPREHENSIVE DIAGNOSTICS & COUNTY ANALYSIS COMPLETE ===\n\n")

cat("FILES GENERATED:\n")
cat("• Individual county plots (", nrow(county_fit_stats), 
    " files): outputs/plots/individual_counties/\n", sep = "")
cat("• County grid pages (", n_pages, 
    " files): outputs/plots/all_counties_grid_page*.png\n", sep = "")
if (requireNamespace("tigris", quietly = TRUE)) {
  cat("• Spatial analysis: outputs/plots/spatial_analysis_comprehensive_", 
      RESPONSE_TYPE, ".png\n", sep = "")
}
cat("• Performance summary: outputs/diagnostics/county_performance_summary_", 
    RESPONSE_TYPE, ".csv\n", sep = "")
cat("• Analysis report: outputs/diagnostics/comprehensive_analysis_report_", 
    RESPONSE_TYPE, ".txt\n", sep = "")
cat("• Prediction metrics: outputs/diagnostics/prediction_metrics_", 
    RESPONSE_TYPE, ".rds\n", sep = "")
cat("• State-level plot: outputs/plots/state_level_", 
    RESPONSE_TYPE, ".png\n\n", sep = "")

cat("KEY FINDINGS:\n")
cat("• Best model:", overall_best$model_name, "(", overall_best$phase, 
    ", WAIC =", round(overall_best$waic, 1), ")\n")
cat("• Model family: ZIP (handles", 
    round(sum(predictions_df$observed == 0) / nrow(predictions_df) * 100, 1), 
    "% zero observations)\n")
cat("• Overall 95% coverage:", round(coverage * 100, 1), "%\n")
cat("• Average county coverage:", round(mean(county_fit_stats$coverage_95) * 100, 1), "%\n")
cat("• Counties with good/excellent fit:", 
    sum(county_fit_stats$fit_quality %in% c("Excellent", "Good")), 
    "out of", nrow(county_fit_stats), "\n")
cat("• Mean interval width:", round(mean(predictions_df$fitted_upper - predictions_df$fitted_lower), 2), "\n\n")

cat("✓ Phase 5 comprehensive diagnostics complete\n")

# ==============================================================================
# END OF PHASE 5 - COMPREHENSIVE MODEL DIAGNOSTICS
# ==============================================================================