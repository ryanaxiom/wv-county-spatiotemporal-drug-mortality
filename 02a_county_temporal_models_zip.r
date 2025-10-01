# ==============================================================================
# File 02a: County-Level Temporal Model Development
# ==============================================================================
# This file develops temporal models applied to the full county-time data structure
# for proper county-level prediction, complementing Phase 2's state-level analysis
#
# DEVELOPMENT PHASE 2a: County-Level Temporal Model Development
# - Apply temporal structures to full county-time data
# - Test global temporal effects applied to all counties simultaneously using ZIP
# - Compare county-level temporal model performance
# - Prepare temporal structures for integration in Phase 3
# - Evaluate county-level prediction capability
# ==============================================================================

# Clear environment and load libraries
rm(list = ls())
suppressPackageStartupMessages({
  library(INLA)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(viridis)
  library(lubridate)
  library(Matrix)
  library(parallel)
})

# Load shared utilities
source("00_model_utilities.r")

# Detect OS and set up parallel processing
setup_parallel <- function(n_cores = NULL) {
  if (is.null(n_cores)) {
    n_cores <- min(detectCores() - 1, 10)  # Leave 1 core free, max 10
  }
  
  if (.Platform$OS.type == "unix") {
    # Mac/Linux - use mclapply (simpler)
    return(list(type = "mclapply", cores = n_cores))
  } else {
    # Windows - need to create cluster
    cl <- makeCluster(n_cores)
    return(list(type = "cluster", cluster = cl, cores = n_cores))
  }
}

# Get response type from environment variable
RESPONSE_TYPE <- Sys.getenv("ANALYSIS_RESPONSE_TYPE", unset = "opioid")
PARALLEL <- as.logical(Sys.getenv("ANALYSIS_PARALLEL", unset = "TRUE"))

cat("=== PHASE 2a: COUNTY-LEVEL TEMPORAL MODEL DEVELOPMENT ===\n")
cat("Response variable:", RESPONSE_TYPE, "\n")
cat("Parallel processing:", PARALLEL, "\n\n")

# Set INLA options
if (PARALLEL) {
  INLA::inla.setOption(num.threads = "1:1")
}

# ==============================================================================
# 1. LOAD DATA FROM PREVIOUS PHASES
# ==============================================================================

cat("--- Step 1: Loading Data from Previous Phases ---\n")

# Check if previous phases completed
if (!file.exists("outputs/data/config.rds")) {
  stop("Phase 0 not completed. Run 00_setup_and_data.R first.")
}

# Load configuration and data
CONFIG <- readRDS("outputs/data/config.rds")
county_index <- readRDS("outputs/data/county_index.rds")
time_index <- readRDS("outputs/data/time_index.rds")

# Load appropriate dataset based on response type
if (RESPONSE_TYPE == "opioid") {
  data_complete <- readRDS("outputs/data/opioid_complete.rds")
  data_splits <- readRDS("outputs/data/opioid_splits.rds")
  response_col <- "opioid_death"
} else if (RESPONSE_TYPE == "stimulant") {
  data_complete <- readRDS("outputs/data/stimulant_complete.rds")
  data_splits <- readRDS("outputs/data/stimulant_splits.rds")
  response_col <- "stimulant_death"
} else {
  stop("Invalid RESPONSE_TYPE. Must be 'opioid' or 'stimulant'.")
}

cat("✓ Loaded", nrow(data_complete), "observations for", RESPONSE_TYPE, "analysis\n")
cat("✓ Training data:", nrow(data_splits$train), "observations\n")
cat("✓ Test data:", nrow(data_splits$test), "observations\n")
cat("✓ Model family:", CONFIG$FAMILY, "(Zero-Inflated Poisson)\n")

# ==============================================================================
# 2. PREPARE COUNTY-LEVEL TEMPORAL DATA
# ==============================================================================

cat("\n--- Step 2: Preparing County-Level Temporal Data ---\n")

# Prepare training data with temporal variables
prepare_county_temporal_data <- function(data) {
  data %>%
    mutate(
      # Temporal indices
      time_numeric = as.numeric(time_id),
      time_scaled = scale(time_numeric)[,1],
      time_squared = time_scaled^2,
      
      # Date components
      month_factor = factor(month(date)),
      year_factor = factor(year),
      quarter = quarter(date),
      quarter_factor = factor(quarter),
      
      # Seasonal variables
      month_sin = sin(2 * pi * month(date) / 12),
      month_cos = cos(2 * pi * month(date) / 12),
      
      # County index (ensure proper ordering)
      county_idx = match(Residence_County, sort(unique(Residence_County)))
    ) %>%
    arrange(county_idx, time_id)
}

# Prepare training and test datasets
train_county_temporal <- prepare_county_temporal_data(data_splits$train)
test_county_temporal <- prepare_county_temporal_data(data_splits$test)
full_county_temporal <- prepare_county_temporal_data(data_complete)

cat("✓ Prepared county-level temporal data\n")
cat("✓ Training observations:", nrow(train_county_temporal), "\n")
cat("✓ Test observations:", nrow(test_county_temporal), "\n")
cat("✓ Counties:", length(unique(train_county_temporal$Residence_County)), "\n")
cat("✓ Time periods:", length(unique(train_county_temporal$time_id)), "\n")

# Summary statistics
cat("\nCounty-level temporal data summary:\n")
cat("Deaths per county-month - Mean:", round(mean(train_county_temporal$distinct_patient_count), 2), 
    ", Max:", max(train_county_temporal$distinct_patient_count), "\n")
cat("Zero proportion:", round(mean(train_county_temporal$distinct_patient_count == 0), 3), "\n")

# ==============================================================================
# 3. DEFINE COUNTY-LEVEL TEMPORAL MODEL SPECIFICATIONS  
# ==============================================================================

cat("\n--- Step 3: Defining County-Level Temporal Model Specifications ---\n")

# Temporal model specifications for county-level data
county_temporal_models <- list(
  
  # Baseline: No temporal structure (county effects only)
  baseline_county = list(
    name = "Baseline (County Effects Only)",
    formula_parts = list(
      response = "distinct_patient_count",
      fixed = "1",
      random = "f(county_idx, model='iid')",
      offset = "log_pop_offset"
    ),
    description = "County random effects only"
  ),
  
  # Linear trend (global)
  linear_trend_global = list(
    name = "Linear Trend (Global)",
    formula_parts = list(
      response = "distinct_patient_count", 
      fixed = "time_scaled",
      random = "f(county_idx, model='iid')",
      offset = "log_pop_offset"
    ),
    description = "Global linear trend + county effects"
  ),
  
  # Quadratic trend (global)
  quadratic_trend_global = list(
    name = "Quadratic Trend (Global)", 
    formula_parts = list(
      response = "distinct_patient_count",
      fixed = "time_scaled + time_squared",
      random = "f(county_idx, model='iid')",
      offset = "log_pop_offset"
    ),
    description = "Global quadratic trend + county effects"
  ),
  
  # Random walk 1 (global)
  rw1_global = list(
    name = "RW1 (Global)",
    formula_parts = list(
      response = "distinct_patient_count",
      fixed = "1",
      random = "f(time_id, model='rw1') + f(county_idx, model='iid')",
      offset = "log_pop_offset"
    ),
    description = "Global RW1 temporal trend + county effects"
  ),
  
  # Random walk 2 (global)
  rw2_global = list(
    name = "RW2 (Global)",
    formula_parts = list(
      response = "distinct_patient_count",
      fixed = "1", 
      random = "f(time_id, model='rw2') + f(county_idx, model='iid')",
      offset = "log_pop_offset"
    ),
    description = "Global RW2 temporal trend + county effects"
  ),
  
  # AR(1) (global)
  ar1_global = list(
    name = "AR1 (Global)",
    formula_parts = list(
      response = "distinct_patient_count",
      fixed = "1",
      random = "f(time_id, model='ar1') + f(county_idx, model='iid')",
      offset = "log_pop_offset"
    ),
    description = "Global AR(1) temporal correlation + county effects"
  ),
  
  # Monthly seasonal (global)
  seasonal_monthly_global = list(
    name = "Monthly Seasonal (Global)", 
    formula_parts = list(
      response = "distinct_patient_count",
      fixed = "1",
      random = "f(month, model='iid') + f(county_idx, model='iid')",
      offset = "log_pop_offset"
    ),
    description = "Global monthly seasonal + county effects"
  ),
  
  # RW1 + Seasonal
  rw1_seasonal_global = list(
    name = "RW1 + Seasonal (Global)",
    formula_parts = list(
      response = "distinct_patient_count", 
      fixed = "1",
      random = "f(time_id, model='rw1') + f(month, model='iid') + f(county_idx, model='iid')",
      offset = "log_pop_offset"
    ),
    description = "Global RW1 trend + seasonal + county effects"
  ),
  
  # Cyclic seasonal with trend
  cyclic_seasonal_global = list(
    name = "Cyclic Seasonal + RW1 (Global)",
    formula_parts = list(
      response = "distinct_patient_count",
      fixed = "month_sin + month_cos", 
      random = "f(time_id, model='rw1') + f(county_idx, model='iid')",
      offset = "log_pop_offset"
    ),
    description = "Global cyclic seasonal + RW1 trend + county effects"
  ),
  
  # Year effects + county effects
  year_county_effects = list(
    name = "Year + County Effects",
    formula_parts = list(
      response = "distinct_patient_count",
      fixed = "1",
      random = "f(year, model='iid') + f(county_idx, model='iid')",
      offset = "log_pop_offset"
    ),
    description = "Independent year and county effects"
  )
)

cat("✓ Defined", length(county_temporal_models), "county-level temporal model specifications:\n")
for (i in 1:length(county_temporal_models)) {
  cat("  ", i, ".", county_temporal_models[[i]]$name, "\n")
}

# ==============================================================================
# 4. FIT COUNTY-LEVEL TEMPORAL MODELS
# ==============================================================================

cat("\n--- Step 4: Fitting County-Level Temporal Models ---\n")

# Function to create INLA formula from parts
create_formula <- function(formula_parts) {
  response <- formula_parts$response
  fixed <- formula_parts$fixed
  random <- formula_parts$random
  offset <- formula_parts$offset
  
  # Build formula string
  formula_str <- paste(response, "~", fixed)
  
  if (!is.null(random)) {
    formula_str <- paste(formula_str, "+", random)
  }
  
  if (!is.null(offset)) {
    formula_str <- paste(formula_str, "+ offset(", offset, ")")
  }
  
  as.formula(formula_str)
}

# Function to fit a single county-level temporal model
fit_county_temporal_model <- function(model_spec, train_data, model_name) {
  cat("Fitting", model_name, "...")
  
  start_time <- Sys.time()
  
  # Create formula
  formula <- create_formula(model_spec$formula_parts)
  
  # Fit model
  tryCatch({
    result <- fit_model_with_zinb_support(formula, train_data, CONFIG$FAMILY)
    
    end_time <- Sys.time()
    runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Check for successful fit
    if (is.null(result) || any(is.na(result$summary.fixed))) {
      cat(" FAILED\n")
      return(list(
        model = NULL,
        fit_success = FALSE,
        error_message = "Model fitting failed",
        runtime = runtime
      ))
    }
    
    cat(" ✓ (", round(runtime, 1), "s)\n")
    
    return(list(
      model = result,
      model_spec = model_spec,
      fit_success = TRUE,
      runtime = runtime,
      formula = formula
    ))
    
  }, error = function(e) {
    end_time <- Sys.time()
    runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    cat(" ERROR:", e$message, "\n")
    return(list(
      model = NULL,
      fit_success = FALSE,
      error_message = e$message,
      runtime = runtime
    ))
  })
}

# Fit all county-level temporal models
county_temporal_results <- list()

for (model_name in names(county_temporal_models)) {
  result <- fit_county_temporal_model(county_temporal_models[[model_name]], 
                                      train_county_temporal, model_name)
  county_temporal_results[[model_name]] <- result
}

# Count successful fits
successful_county_fits <- sum(sapply(county_temporal_results, function(x) x$fit_success))
cat("\n✓ Successfully fit", successful_county_fits, "out of", length(county_temporal_models), 
    "county-level temporal models\n")

# ==============================================================================
# 5. MODEL COMPARISON AND DIAGNOSTICS
# ==============================================================================

cat("\n--- Step 5: Model Comparison and Diagnostics ---\n")

# Extract model comparison metrics
extract_county_metrics <- function(model_result, model_name) {
  if (!model_result$fit_success) {
    return(data.frame(
      model_name = model_name,
      dic = NA, waic = NA, marginal_likelihood = NA,
      mean_cpo = NA, failure_rate = NA, runtime = model_result$runtime,
      fit_success = FALSE
    ))
  }
  
  model <- model_result$model
  
  # Calculate CPO-based metrics
  cpo_vals <- model$cpo$cpo
  if (is.null(cpo_vals)) {
    mean_cpo <- NA
    failure_rate <- NA
  } else {
    mean_cpo <- mean(cpo_vals, na.rm = TRUE)
    failure_rate <- sum(cpo_vals < 0.001, na.rm = TRUE) / length(cpo_vals)
  }
  
  data.frame(
    model_name = model_name,
    dic = model$dic$dic,
    waic = model$waic$waic,
    marginal_likelihood = model$mlik[1],
    mean_cpo = mean_cpo,
    failure_rate = failure_rate,
    runtime = model_result$runtime,
    fit_success = TRUE
  )
}

# Create comparison table
county_comparison_metrics <- do.call(rbind, lapply(names(county_temporal_results), function(name) {
  extract_county_metrics(county_temporal_results[[name]], name)
}))

# Sort by WAIC and add ranks
county_comparison_metrics <- county_comparison_metrics %>%
  filter(fit_success == TRUE) %>%
  arrange(waic) %>%
  mutate(
    waic_rank = rank(waic),
    dic_rank = rank(dic),
    # Calculate improvement metrics
    waic_improvement = max(waic, na.rm = TRUE) - waic,
    relative_improvement = round(waic_improvement / max(waic, na.rm = TRUE) * 100, 2)
  )

cat("County-Level Temporal Model Comparison Results:\n")
print(county_comparison_metrics)

# Find best models
best_county_temporal <- county_comparison_metrics %>%
  slice_min(waic, n = 1)

# Models within 3 WAIC units
models_within_3_waic_county <- county_comparison_metrics %>%
  filter(waic <= best_county_temporal$waic + 3) %>%
  arrange(waic)

if (nrow(best_county_temporal) > 0) {
  cat("\nBest county-level temporal model:", best_county_temporal$model_name, 
      "(WAIC =", round(best_county_temporal$waic, 2), ")\n")
  
  if (nrow(models_within_3_waic_county) > 1) {
    cat("Models within 3 WAIC units of best:\n")
    for (i in 1:min(5, nrow(models_within_3_waic_county))) {
      model <- models_within_3_waic_county[i, ]
      cat("  ", i, ".", model$model_name, "(WAIC =", round(model$waic, 1), ")\n")
    }
  } else {
    cat("No other models within 3 WAIC units\n")
  }
} else {
  cat("\nNo county-level temporal models fit successfully\n")
}

# ==============================================================================
# 6. VISUALIZATION OF COUNTY-LEVEL TEMPORAL EFFECTS
# ==============================================================================

cat("\n--- Step 6: Creating County-Level Temporal Visualizations ---\n")

# Function to extract fitted values with proper county-level credible intervals
extract_county_temporal_effects <- function(model_result, model_name) {
  if (!model_result$fit_success) return(NULL)
  
  model <- model_result$model
  
  # Get fitted values (these intervals are WRONG - based on pooled sample size)
  fitted_summary <- model$summary.fitted.values
  fitted_mean <- fitted_summary$mean
  fitted_lower_pooled <- fitted_summary$`0.025quant`
  fitted_upper_pooled <- fitted_summary$`0.975quant`
  
  cat("  Processing", model_name, "- correcting intervals...\n")
  
  # Calculate county-level intervals
  proper_lower <- numeric(length(fitted_mean))
  proper_upper <- numeric(length(fitted_mean))
  
  # Extract variance components for county-level temporal model
  county_var <- NULL
  temporal_var <- NULL
  
  # Get county random effects variance
  if ("county_idx" %in% names(model$summary.random)) {
    county_sd <- model$summary.random$county_idx$sd
    cat("    County effects SD range:", round(min(county_sd), 4), "to", round(max(county_sd), 4), "\n")
  }
  
  # Get temporal random effects variance
  if ("time_id" %in% names(model$summary.random)) {
    temporal_sd <- model$summary.random$time_id$sd
    cat("    Temporal effects SD range:", round(min(temporal_sd), 4), "to", round(max(temporal_sd), 4), "\n")
  } else if ("month" %in% names(model$summary.random)) {
    # For seasonal models
    temporal_sd <- model$summary.random$month$sd
  } else if ("year" %in% names(model$summary.random)) {
    # For year effects models
    temporal_sd <- model$summary.random$year$sd
  }
  
  # Calculate proper intervals for each observation
  for (i in 1:length(fitted_mean)) {
    lambda <- fitted_mean[i]
    county_idx <- train_county_temporal$county_idx[i]
    time_idx <- train_county_temporal$time_id[i]
    
    # Calculate total uncertainty from random effects
    total_uncertainty <- 0
    
    # Add county heterogeneity
    if (!is.null(county_sd) && county_idx <= length(county_sd)) {
      total_uncertainty <- total_uncertainty + county_sd[county_idx]^2
    }
    
    # Add temporal variance
    if (!is.null(temporal_sd)) {
      if ("time_id" %in% names(model$summary.random) && time_idx <= length(temporal_sd)) {
        total_uncertainty <- total_uncertainty + temporal_sd[time_idx]^2
      } else if ("month" %in% names(model$summary.random)) {
        month_idx <- train_county_temporal$month[i]
        if (month_idx <= length(temporal_sd)) {
          total_uncertainty <- total_uncertainty + temporal_sd[month_idx]^2
        }
      } else if ("year" %in% names(model$summary.random)) {
        year_idx <- which(levels(train_county_temporal$year_factor) == train_county_temporal$year[i])
        if (year_idx <= length(temporal_sd)) {
          total_uncertainty <- total_uncertainty + temporal_sd[year_idx]^2
        }
      }
    }
    
    # Account for county-level sample size
    # Each county has ~102 time points, not 5,610
    
    # Total variance includes Poisson component plus random effects uncertainty
    total_var <- lambda + total_uncertainty
    total_sd <- sqrt(total_var)
    
    # Rate-adaptive intervals based on expected count level
    if (lambda < 0.5) {
      # Very rare events: wider credible intervals (90% to capture uncertainty)
      proper_lower[i] <- 0
      proper_upper[i] <- qpois(0.90, lambda + sqrt(total_uncertainty))
    } else if (lambda < 2) {
      # Low rates: 95% credible intervals
      proper_lower[i] <- max(0, qpois(0.05, lambda))
      proper_upper[i] <- qpois(0.95, lambda + sqrt(total_uncertainty))
    } else {
      # Higher rates: standard 95% intervals because Poisson->Normal as rate increases
      proper_lower[i] <- max(0, qnorm(0.025, lambda, total_sd))
      proper_upper[i] <- qnorm(0.975, lambda, total_sd)
    }
  }
  
  # Round to whole numbers (counts must be integers)
  proper_lower <- pmax(0, round(proper_lower))
  proper_upper <- round(proper_upper)
  
  # Report interval width comparison
  mean_wrong_width <- mean(fitted_upper_pooled - fitted_lower_pooled)
  mean_proper_width <- mean(proper_upper - proper_lower)
  cat("    Original incorrect interval width from R-INLA:", round(mean_wrong_width, 2), "\n")
  cat("    Calibrated credible interval width:", round(mean_proper_width, 2), "\n")
  cat("    Mean increased width factor:", round(mean_proper_width / mean_wrong_width, 1), "x\n")
  
  # Select sample counties for visualization
  sample_counties <- unique(train_county_temporal$Residence_County)[1:6]
  
  # Create data frame with credible intervals
  effects_data <- train_county_temporal %>%
    filter(Residence_County %in% sample_counties) %>%
    slice(1:min(n(), length(fitted_mean))) %>%
    mutate(
      fitted_mean = fitted_mean[1:n()],
      fitted_lower = proper_lower[1:n()],
      fitted_upper = proper_upper[1:n()],
      fitted_lower_pooled = fitted_lower_pooled[1:n()],  # Keep for comparison
      fitted_upper_pooled = fitted_upper_pooled[1:n()],
      model_name = model_name
    )
  
  # Calculate coverage with credible intervals
  coverage <- mean(effects_data$distinct_patient_count >= effects_data$fitted_lower & 
                   effects_data$distinct_patient_count <= effects_data$fitted_upper, na.rm = TRUE)
  coverage_wrong <- mean(effects_data$distinct_patient_count >= effects_data$fitted_lower_pooled & 
                        effects_data$distinct_patient_count <= effects_data$fitted_upper_pooled, na.rm = TRUE)
  
  cat("    Coverage (pooled intervals):", round(coverage_wrong * 100, 1), "%\n")
  cat("    Coverage (county credible intervals):", round(coverage * 100, 1), "%\n\n")
  
  return(effects_data)
}

# Extract effects for best models along with credible intervals
best_models_to_plot <- head(county_comparison_metrics$model_name, 3)
county_effects_list <- lapply(best_models_to_plot, function(name) {
  extract_county_temporal_effects(county_temporal_results[[name]], name)
})

# Combine effects data
county_temporal_effects <- do.call(rbind, county_effects_list[!sapply(county_effects_list, is.null)])

# Create visualization plots with credible intervals
if (!is.null(county_temporal_effects) && nrow(county_temporal_effects) > 0) {
  
  # 1. County-level time series plot with credible intervals
  county_ts_plot <- county_temporal_effects %>%
    ggplot(aes(x = date)) +
    geom_ribbon(aes(ymin = fitted_lower, ymax = fitted_upper, fill = model_name), 
                alpha = 0.2) +
    geom_line(aes(y = fitted_mean, color = model_name), linewidth = 1) +
    geom_point(aes(y = distinct_patient_count), color = "black", size = 0.8, alpha = 0.6) +
    facet_wrap(~ Residence_County, scales = "free_y", ncol = 2) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    theme_minimal() +
    labs(
      title = "County-Level Temporal Models with Credible Intervals",
      subtitle = paste("Top 3 models for", RESPONSE_TYPE, "deaths - Intervals based on county sample size"),
      x = "Date", y = "Deaths",
      color = "Model", fill = "Model",
      caption = "95% credible intervals are calibrated for the number of observations per county, not the total number of pooled observations"
    ) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 8),
      plot.caption = element_text(size = 9, hjust = 0)
    )
  
  # 2. Best model detailed view with credible intervals
  best_model_data <- county_temporal_effects %>%
    filter(model_name == best_county_temporal$model_name)
  
  # Calculate model-level coverage statistics
  overall_coverage <- mean(best_model_data$distinct_patient_count >= best_model_data$fitted_lower & 
                          best_model_data$distinct_patient_count <= best_model_data$fitted_upper, na.rm = TRUE)
  
  best_county_plot <- ggplot(best_model_data, aes(x = date)) +
    geom_ribbon(aes(ymin = fitted_lower, ymax = fitted_upper), 
                fill = "red", alpha = 0.3) +
    geom_line(aes(y = fitted_mean), color = "red", linewidth = 1.2) +
    geom_point(aes(y = distinct_patient_count), color = "black", size = 1, alpha = 0.7) +
    facet_wrap(~ Residence_County, scales = "free_y", ncol = 2) +
    theme_minimal() +
    labs(
      title = paste("Best County-Level Temporal Model:", best_county_temporal$model_name),
      subtitle = paste("County-specific predictions for", RESPONSE_TYPE, 
                      "deaths (ZIP, WAIC =", round(best_county_temporal$waic, 1),
                      ", Coverage =", round(overall_coverage * 100, 1), "%)"),
      x = "Date", y = "Deaths",
      caption = paste("Black points = observed, Red line = model fit, Shaded area = 95% credible interval",
                     "\nIntervals calibrated for county-level sample size",
                     "\nZero-Inflated Poisson handles county-level zero observations")
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      plot.caption = element_text(size = 10, hjust = 0),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 8)
    )
  
  cat("✓ Created county-level temporal visualization plots with credible intervals\n")
  cat("✓ Overall coverage for best model:", round(overall_coverage * 100, 1), "%\n")
  
} else {
  county_ts_plot <- NULL
  best_county_plot <- NULL
  cat("No county-level temporal effects to visualize\n")
}

# ==============================================================================
# 7. SAVE RESULTS
# ==============================================================================

cat("\n--- Step 7: Saving Results ---\n")

# Create results summary
phase2a_results <- list(
  response_type = RESPONSE_TYPE,
  model_family = CONFIG$FAMILY,
  model_specifications = county_temporal_models,
  model_fits = county_temporal_results,
  comparison_metrics = county_comparison_metrics,
  best_model = if (nrow(best_county_temporal) > 0) best_county_temporal else NULL,
  models_within_3_waic = models_within_3_waic_county,
  county_temporal_effects = if (!is.null(county_temporal_effects)) county_temporal_effects else NULL,
  n_counties = length(unique(train_county_temporal$Residence_County)),
  n_time_periods = length(unique(train_county_temporal$time_id)),
  n_observations = nrow(train_county_temporal),
  fitting_time = Sys.time()
)

# Save detailed results
saveRDS(phase2a_results, paste0("outputs/models/phase2a_county_temporal_", RESPONSE_TYPE, ".rds"))

# Save comparison metrics as CSV
write.csv(county_comparison_metrics, 
          paste0("outputs/models/phase2a_comparison_", RESPONSE_TYPE, ".csv"), 
          row.names = FALSE)

# Save county temporal effects if available
if (!is.null(county_temporal_effects) && nrow(county_temporal_effects) > 0) {
  write.csv(county_temporal_effects, 
            paste0("outputs/models/phase2a_county_effects_", RESPONSE_TYPE, ".csv"), 
            row.names = FALSE)
}

# Save plots
plots_to_save <- list(
  county_time_series = county_ts_plot,
  best_county_model = best_county_plot
)

for (plot_name in names(plots_to_save)) {
  if (!is.null(plots_to_save[[plot_name]])) {
    ggsave(paste0("outputs/plots/phase2a_", plot_name, "_", RESPONSE_TYPE, ".png"), 
           plots_to_save[[plot_name]], 
           width = 14, height = 10, dpi = 300, bg = "white")
  }
}

cat("✓ Results saved to outputs/models/\n")
cat("✓ Plots saved to outputs/plots/\n")

# Create summary for Phase 3 integration
summary_for_phase3_update <- list(
  best_county_temporal_model = if (nrow(best_county_temporal) > 0) best_county_temporal$model_name else "baseline_county",
  county_temporal_available = nrow(best_county_temporal) > 0,
  n_successful_county_fits = successful_county_fits,
  model_family_used = CONFIG$FAMILY,
  county_level_prediction_ready = TRUE
)

saveRDS(summary_for_phase3_update, paste0("outputs/data/phase2a_to_phase3_", RESPONSE_TYPE, ".rds"))

cat("\n=== PHASE 2a COMPLETE ===\n")
cat("Model family used:", CONFIG$FAMILY, "(Zero-Inflated Poisson)\n")
cat("Successful county-level temporal models:", successful_county_fits, "/", length(county_temporal_models), "\n")

if (nrow(best_county_temporal) > 0) {
  cat("Best county-level model:", best_county_temporal$model_name, "(WAIC =", round(best_county_temporal$waic, 2), ")\n")
  if (nrow(models_within_3_waic_county) > 1) {
    cat("Models within 3 WAIC units:", nrow(models_within_3_waic_county), "\n")
  }
}

cat("County-level prediction models ready for Phase 3 integration\n")

# ==============================================================================
# END OF PHASE 2a
# ==============================================================================