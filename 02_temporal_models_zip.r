# ==============================================================================
# File 02: Temporal Model Development for WV County Drug Death Analysis
# ==============================================================================
# This file develops and compares different temporal model structures using R-INLA
#
# DEVELOPMENT PHASE 2: Temporal Model Development
# - AR(1) temporal correlation models
# - Random walk smoothing (RW1, RW2) for temporal trends
# - Seasonal/cyclical components (yearly patterns)
# - Test temporal structures on spatially-aggregated data using Zero-Inflated Poisson
# - Evaluate temporal model fit and forecasting capability
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
  library(forecast)
  library(parallel)
})

# Load shared utilities
source("00_model_utilities.r")

# Get response type from environment variable (set by master script)
RESPONSE_TYPE <- Sys.getenv("ANALYSIS_RESPONSE_TYPE", unset = "opioid")
PARALLEL <- as.logical(Sys.getenv("ANALYSIS_PARALLEL", unset = "TRUE"))

cat("=== PHASE 2: TEMPORAL MODEL DEVELOPMENT ===\n")
cat("Response variable:", RESPONSE_TYPE, "\n")
cat("Parallel processing:", PARALLEL, "\n\n")

# Set INLA options
if (PARALLEL) {
  INLA::inla.setOption(num.threads = "1:1")
}

# ==============================================================================
# 1. LOAD DATA FROM PHASE 0
# ==============================================================================

cat("--- Step 1: Loading Data from Phase 0 ---\n")

# Check if Phase 0 completed
if (!file.exists("outputs/data/config.rds")) {
  stop("Phase 0 not completed. Run 00_setup_and_data.R first.")
}

# Load configuration and data
CONFIG <- readRDS("outputs/data/config.rds")
time_index <- readRDS("outputs/data/time_index.rds")
county_index <- readRDS("outputs/data/county_index.rds")

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
# 2. CREATE TEMPORALLY-AGGREGATED DATA
# ==============================================================================

cat("\n--- Step 2: Creating Temporally-Aggregated Data ---\n")

# Aggregate data across all counties by time to focus on temporal patterns
create_temporal_data <- function(data) {
  data %>%
    group_by(time_id, year_month, year, month, date) %>%
    summarise(
      total_count = sum(distinct_patient_count, na.rm = TRUE),
      total_population = sum(population, na.rm = TRUE),
      n_counties = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      log_pop_offset = log(total_population),
      time_numeric = as.numeric(time_id),
      year_factor = factor(year),
      month_factor = factor(month),
      quarter = quarter(date),
      quarter_factor = factor(quarter),
      # Create seasonal indices
      month_sin = sin(2 * pi * month / 12),
      month_cos = cos(2 * pi * month / 12),
      # Create trend variables
      time_scaled = scale(time_numeric)[,1],
      time_squared = time_scaled^2
    ) %>%
    arrange(time_id)
}

# Create aggregated training and test data
train_temporal <- create_temporal_data(data_splits$train)
test_temporal <- create_temporal_data(data_splits$test)
full_temporal <- create_temporal_data(data_complete)

cat("✓ Created temporal aggregated data\n")
cat("✓ Training periods:", nrow(train_temporal), "months\n")
cat("✓ Test periods:", nrow(test_temporal), "months\n")
cat("✓ Date range:", min(full_temporal$date), "to", max(full_temporal$date), "\n")

# Check for temporal patterns
cat("\nTemporal data summary:\n")
cat("Total deaths over time - Min:", min(full_temporal$total_count), 
    "Max:", max(full_temporal$total_count), 
    "Mean:", round(mean(full_temporal$total_count), 1), "\n")

# ==============================================================================
# 3. DEFINE TEMPORAL MODEL SPECIFICATIONS
# ==============================================================================

cat("\n--- Step 3: Defining Temporal Model Specifications ---\n")

# Model specifications to test
temporal_models <- list(
  
  # Baseline: No temporal structure
  baseline = list(
    name = "Baseline (No Temporal)",
    formula_parts = list(
      response = "total_count",
      fixed = "1",
      random = NULL,
      offset = "log_pop_offset"
    ),
    description = "Intercept only with population offset"
  ),
  
  # Linear trend
  linear_trend = list(
    name = "Linear Trend",
    formula_parts = list(
      response = "total_count",
      fixed = "time_scaled",
      random = NULL,
      offset = "log_pop_offset"
    ),
    description = "Linear time trend"
  ),
  
  # Quadratic trend
  quadratic_trend = list(
    name = "Quadratic Trend",
    formula_parts = list(
      response = "total_count",
      fixed = "time_scaled + time_squared",
      random = NULL,
      offset = "log_pop_offset"
    ),
    description = "Quadratic time trend"
  ),
  
  # Random walk order 1 (RW1)
  rw1 = list(
    name = "Random Walk 1 (RW1)",
    formula_parts = list(
      response = "total_count",
      fixed = "1",
      random = "f(time_id, model='rw1')",
      offset = "log_pop_offset"
    ),
    description = "Random walk order 1 - smooth trends"
  ),
  
  # Random walk order 2 (RW2)
  rw2 = list(
    name = "Random Walk 2 (RW2)",
    formula_parts = list(
      response = "total_count",
      fixed = "1",
      random = "f(time_id, model='rw2')",
      offset = "log_pop_offset"
    ),
    description = "Random walk order 2 - smoother trends"
  ),
  
  # AR(1) temporal correlation
  ar1 = list(
    name = "AR(1) Correlation",
    formula_parts = list(
      response = "total_count",
      fixed = "1",
      random = "f(time_id, model='ar1')",
      offset = "log_pop_offset"
    ),
    description = "AR(1) temporal correlation"
  ),
  
  # Seasonal (monthly) effects
  seasonal_monthly = list(
    name = "Monthly Seasonal",
    formula_parts = list(
      response = "total_count",
      fixed = "1",
      random = "f(month, model='iid')",
      offset = "log_pop_offset"
    ),
    description = "Monthly seasonal effects"
  ),
  
  # Seasonal with trend
  seasonal_rw1 = list(
    name = "RW1 + Monthly Seasonal",
    formula_parts = list(
      response = "total_count",
      fixed = "1",
      random = "f(time_id, model='rw1') + f(month, model='iid')",
      offset = "log_pop_offset"
    ),
    description = "Random walk with monthly seasonal effects"
  ),
  
  # Cyclic seasonal (smoother)
  cyclic_seasonal = list(
    name = "Cyclic Seasonal",
    formula_parts = list(
      response = "total_count",
      fixed = "month_sin + month_cos",
      random = "f(time_id, model='rw1')",
      offset = "log_pop_offset"
    ),
    description = "Smooth cyclic seasonal pattern with RW1 trend"
  ),
  
  # Year effects + month effects
  year_month_effects = list(
    name = "Year + Month Effects",
    formula_parts = list(
      response = "total_count",
      fixed = "1",
      random = "f(year, model='iid') + f(month, model='iid')",
      offset = "log_pop_offset"
    ),
    description = "Independent year and month effects"
  )
)

cat("✓ Defined", length(temporal_models), "temporal model specifications:\n")
for (i in 1:length(temporal_models)) {
  cat("  ", i, ".", temporal_models[[i]]$name, "\n")
}

# ==============================================================================
# 4. FIT TEMPORAL MODELS
# ==============================================================================

cat("\n--- Step 4: Fitting Temporal Models ---\n")

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

# Function to fit a single temporal model
fit_temporal_model <- function(model_spec, train_data, model_name) {
  cat("Fitting", model_name, "...\n")
  
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
      cat("  ❌ Failed to fit", model_name, "\n")
      return(list(
        model = NULL,
        fit_success = FALSE,
        error_message = "Model fitting failed",
        runtime = runtime
      ))
    }
    
    cat("  ✓ Completed", model_name, "in", round(runtime, 2), "seconds\n")
    
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
    
    cat("  ❌ Error fitting", model_name, ":", e$message, "\n")
    return(list(
      model = NULL,
      fit_success = FALSE,
      error_message = e$message,
      runtime = runtime
    ))
  })
}

# Fit all temporal models
temporal_results <- list()

for (model_name in names(temporal_models)) {
  result <- fit_temporal_model(temporal_models[[model_name]], train_temporal, model_name)
  temporal_results[[model_name]] <- result
}

# Count successful fits
successful_fits <- sum(sapply(temporal_results, function(x) x$fit_success))
cat("\n✓ Successfully fit", successful_fits, "out of", length(temporal_models), "temporal models\n")

# ==============================================================================
# 5. MODEL COMPARISON AND DIAGNOSTICS
# ==============================================================================

cat("\n--- Step 5: Model Comparison and Diagnostics ---\n")

# Extract model comparison metrics
extract_metrics <- function(model_result, model_name) {
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
comparison_metrics <- do.call(rbind, lapply(names(temporal_results), function(name) {
  extract_metrics(temporal_results[[name]], name)
}))

# Sort by WAIC (lower is better)
comparison_metrics <- comparison_metrics %>%
  arrange(waic) %>%
  mutate(
    waic_rank = ifelse(fit_success, rank(waic, na.last = TRUE), NA),
    dic_rank = ifelse(fit_success, rank(dic, na.last = TRUE), NA)
  )

cat("Temporal Model Comparison Results:\n")
print(comparison_metrics)

# Find best model
best_temporal <- comparison_metrics %>%
  filter(fit_success == TRUE) %>%
  slice_min(waic, n = 1)

if (nrow(best_temporal) > 0) {
  cat("\n✓ Best temporal model:", best_temporal$model_name, "(WAIC =", round(best_temporal$waic, 2), ")\n")
} else {
  cat("\n❌ No temporal models fit successfully\n")
}

# ==============================================================================
# 6. TEMPORAL EFFECTS VISUALIZATION
# ==============================================================================

cat("\n--- Step 6: Creating Temporal Effects Visualizations ---\n")

# Function to extract temporal effects
extract_temporal_effects <- function(model_result, model_name) {
  if (!model_result$fit_success) return(NULL)
  
  model <- model_result$model
  
  # Get fitted values
  fitted_values <- model$summary.fitted.values
  
  # Combine with temporal data
  effects_data <- train_temporal %>%
    mutate(
      fitted_mean = fitted_values$mean,
      fitted_lower = fitted_values$`0.025quant`,
      fitted_upper = fitted_values$`0.975quant`,
      model_name = model_name
    )
  
  # Extract random effects if they exist
  if ("time_id" %in% names(model$summary.random)) {
    time_effects <- model$summary.random$time_id %>%
      mutate(
        time_id = 1:nrow(.),
        temporal_effect = mean,
        temporal_lower = `0.025quant`,
        temporal_upper = `0.975quant`
      ) %>%
      select(time_id, temporal_effect, temporal_lower, temporal_upper)
    
    effects_data <- effects_data %>%
      left_join(time_effects, by = "time_id")
  }
  
  return(effects_data)
}

# Extract temporal effects for successful models
temporal_effects_list <- lapply(names(temporal_results), function(name) {
  if (temporal_results[[name]]$fit_success) {
    extract_temporal_effects(temporal_results[[name]], name)
  } else {
    NULL
  }
})

# Filter out NULL results and ensure consistent column structure
temporal_effects_list <- temporal_effects_list[!sapply(temporal_effects_list, is.null)]

# Check if we have any effects to combine
if (length(temporal_effects_list) > 0) {
  # Get all unique column names across all data frames
  all_cols <- unique(unlist(lapply(temporal_effects_list, names)))
  
  # Ensure each data frame has all columns (fill missing with NA)
  temporal_effects_list <- lapply(temporal_effects_list, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    if (length(missing_cols) > 0) {
      df[missing_cols] <- NA
    }
    return(df[all_cols])  # Reorder columns consistently
  })
  
  # Now combine with consistent column structure
  temporal_effects <- do.call(rbind, temporal_effects_list)
} else {
  temporal_effects <- NULL
}

# Create temporal visualization plots
if (!is.null(temporal_effects) && nrow(temporal_effects) > 0) {
  
  # 1. Time series plot of observed vs fitted
  ts_plot <- temporal_effects %>%
    filter(model_name %in% head(comparison_metrics$model_name[comparison_metrics$fit_success], 4)) %>%
    ggplot(aes(x = date)) +
    geom_point(aes(y = total_count), color = "black", size = 1, alpha = 0.6) +
    geom_line(aes(y = fitted_mean, color = model_name), linewidth = 1) +
    geom_ribbon(aes(ymin = fitted_lower, ymax = fitted_upper, fill = model_name), 
                alpha = 0.2) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    theme_minimal() +
    labs(
      title = paste("Temporal Models: Observed vs Fitted Values (ZIP)"),
      subtitle = paste("Response:", RESPONSE_TYPE, "deaths - Zero-Inflated Poisson"),
      x = "Date", y = "Total Deaths",
      color = "Model", fill = "Model"
    ) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # 2. Seasonal pattern plot (if seasonal effects exist)
  if (any(grepl("month", names(temporal_results)))) {
    seasonal_data <- full_temporal %>%
      group_by(month) %>%
      summarise(
        mean_count = mean(total_count),
        se_count = sd(total_count) / sqrt(n()),
        .groups = 'drop'
      )
    
    seasonal_plot <- ggplot(seasonal_data, aes(x = month, y = mean_count)) +
      geom_col(fill = "steelblue", alpha = 0.7) +
      geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count),
                    width = 0.2) +
      theme_minimal() +
      labs(
        title = "Seasonal Pattern in Drug Deaths",
        subtitle = paste("Average monthly", RESPONSE_TYPE, "deaths ± SE (ZIP model)"),
        x = "Month", y = "Average Deaths"
      ) +
      scale_x_continuous(breaks = 1:12, 
                        labels = month.abb)
  } else {
    seasonal_plot <- NULL
  }
  
  # 3. Best model detailed plot
  if (nrow(best_temporal) > 0) {
    best_model_data <- temporal_effects %>%
      filter(model_name == best_temporal$model_name)
    
    best_plot <- ggplot(best_model_data, aes(x = date)) +
      geom_point(aes(y = total_count), color = "black", size = 2, alpha = 0.7) +
      geom_line(aes(y = fitted_mean), color = "red", linewidth = 1.2) +
      geom_ribbon(aes(ymin = fitted_lower, ymax = fitted_upper), 
                  fill = "red", alpha = 0.3) +
      theme_minimal() +
      labs(
        title = paste("Best Temporal Model:", best_temporal$model_name),
        subtitle = paste("WV", RESPONSE_TYPE, "deaths over time (ZIP, WAIC =", round(best_temporal$waic, 1), ")"),
        x = "Date", y = "Total Deaths",
        caption = "Black points = observed, Red line = model fit, Shaded area = 95% credible interval\nZero-Inflated Poisson model handles temporal aggregation of county zeros"
      ) +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        plot.caption = element_text(size = 10, hjust = 0),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  } else {
    best_plot <- NULL
  }
  
  cat("✓ Created temporal visualization plots\n")
} else {
  ts_plot <- NULL
  seasonal_plot <- NULL  
  best_plot <- NULL
  cat("No temporal effects to visualize\n")
}

# ==============================================================================
# 7. SAVE RESULTS
# ==============================================================================

cat("\n--- Step 7: Saving Results ---\n")

# Create results summary
phase2_results <- list(
  response_type = RESPONSE_TYPE,
  model_family = CONFIG$FAMILY,
  model_specifications = temporal_models,
  model_fits = temporal_results,
  comparison_metrics = comparison_metrics,
  best_model = if (nrow(best_temporal) > 0) best_temporal else NULL,
  temporal_effects = if (!is.null(temporal_effects)) temporal_effects else NULL,
  temporal_data = list(
    train = train_temporal,
    test = test_temporal,
    full = full_temporal
  ),
  n_time_periods = nrow(train_temporal),
  date_range = c(min(full_temporal$date), max(full_temporal$date)),
  fitting_time = Sys.time()
)

# Save detailed results
saveRDS(phase2_results, paste0("outputs/models/phase2_temporal_", RESPONSE_TYPE, ".rds"))

# Save comparison metrics as CSV
write.csv(comparison_metrics, 
          paste0("outputs/models/phase2_comparison_", RESPONSE_TYPE, ".csv"), 
          row.names = FALSE)

# Save temporal effects if available
if (!is.null(temporal_effects) && nrow(temporal_effects) > 0) {
  write.csv(temporal_effects, 
            paste0("outputs/models/phase2_temporal_effects_", RESPONSE_TYPE, ".csv"), 
            row.names = FALSE)
}

# Save plots
plots_to_save <- list(
  time_series = ts_plot,
  seasonal = seasonal_plot,
  best_model = best_plot
)

for (plot_name in names(plots_to_save)) {
  if (!is.null(plots_to_save[[plot_name]])) {
    ggsave(paste0("outputs/plots/phase2_", plot_name, "_", RESPONSE_TYPE, ".png"), 
           plots_to_save[[plot_name]], 
           width = 12, height = 8, dpi = 300, bg = "white")
  }
}

cat("✓ Results saved to outputs/models/\n")
cat("✓ Plots saved to outputs/plots/\n")

# Create summary for next phase (maintains pipeline compatibility)
summary_for_phase3 <- list(
  best_temporal_model = if (nrow(best_temporal) > 0) best_temporal$model_name else "baseline",
  temporal_model_available = nrow(best_temporal) > 0,
  n_successful_fits = successful_fits,
  temporal_improvement = if (nrow(best_temporal) > 0) {
    baseline_waic <- comparison_metrics$waic[comparison_metrics$model_name == "baseline"]
    if (length(baseline_waic) > 0) {
      round((1 - best_temporal$waic/baseline_waic) * 100, 1)
    } else {
      NA
    }
  } else {
    NA
  },
  model_family_used = CONFIG$FAMILY,
  data_ready_for_phase3 = TRUE
)

saveRDS(summary_for_phase3, paste0("outputs/data/phase2_to_phase3_", RESPONSE_TYPE, ".rds"))

# Calculate response-scale improvement metrics for better interpretation
response_scale_improvement <- NULL
if (nrow(best_temporal) > 0) {
  # Extract fitted values for baseline and best model
  baseline_result <- temporal_results[["baseline"]]
  best_result <- temporal_results[[best_temporal$model_name]]
  
  if (!is.null(baseline_result$model) && !is.null(best_result$model)) {
    # Get fitted values
    baseline_fitted <- baseline_result$model$summary.fitted.values$mean
    best_fitted <- best_result$model$summary.fitted.values$mean
    observed <- train_temporal$total_count
    
    # Calculate prediction accuracy metrics
    baseline_mae <- mean(abs(observed - baseline_fitted))
    best_mae <- mean(abs(observed - best_fitted))
    
    baseline_rmse <- sqrt(mean((observed - baseline_fitted)^2))
    best_rmse <- sqrt(mean((observed - best_fitted)^2))
    
    # Percentage improvements in prediction accuracy
    mae_improvement <- round((1 - best_mae/baseline_mae) * 100, 1)
    rmse_improvement <- round((1 - best_rmse/baseline_rmse) * 100, 1)
    
    # WAIC difference for proper interpretation
    waic_diff <- comparison_metrics$waic[comparison_metrics$model_name == "baseline"] - best_temporal$waic
    
    response_scale_improvement <- list(
      mae_improvement = mae_improvement,
      rmse_improvement = rmse_improvement,
      waic_difference = waic_diff,
      baseline_mae = baseline_mae,
      best_mae = best_mae
    )
  }
}

cat("\n=== PHASE 2 COMPLETE ===\n")
cat("Model family used:", CONFIG$FAMILY, "(Zero-Inflated Poisson)\n")
cat("Successful temporal models:", successful_fits, "/", length(temporal_models), "\n")

if (nrow(best_temporal) > 0) {
  cat("Best model:", best_temporal$model_name, "\n")
  cat("Best WAIC:", round(best_temporal$waic, 2), "\n")
  
  # Report WAIC improvement properly
  baseline_waic <- comparison_metrics$waic[comparison_metrics$model_name == "baseline"]
  if (length(baseline_waic) > 0) {
    waic_improvement <- baseline_waic - best_temporal$waic
    cat("WAIC improvement:", round(waic_improvement, 1), "units (lower is better)\n")
    
    # Approximate interpretation: WAIC difference of 10 is strong evidence
    if (waic_improvement > 20) {
      cat("  Evidence strength: Very strong (>20 WAIC units)\n")
    } else if (waic_improvement > 10) {
      cat("  Evidence strength: Strong (>10 WAIC units)\n")
    } else if (waic_improvement > 6) {
      cat("  Evidence strength: Positive (>6 WAIC units)\n")
    } else if (waic_improvement > 2) {
      cat("  Evidence strength: Weak (>2 WAIC units)\n")
    }
  }
  
  # Report response-scale improvements
  if (!is.null(response_scale_improvement)) {
    cat("\nPrediction accuracy improvements (response scale):\n")
    cat("  Mean Absolute Error:", response_scale_improvement$mae_improvement, "% better\n")
    cat("  Root Mean Square Error:", response_scale_improvement$rmse_improvement, "% better\n")
    cat("  Baseline MAE:", round(response_scale_improvement$baseline_mae, 1), "deaths/month\n")
    cat("  Best model MAE:", round(response_scale_improvement$best_mae, 1), "deaths/month\n")
  }
  
  # Note about scale interpretation
  cat("\nNote: WAIC improvements are on log-likelihood scale, not directly interpretable %\n")
  cat("      Response-scale metrics (MAE, RMSE) show actual prediction improvements\n")
}

cat("\n✓ Ready for Phase 2a: County Temporal Models\n")

# ==============================================================================
# END OF PHASE 2
# ==============================================================================