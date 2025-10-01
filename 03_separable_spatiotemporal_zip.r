# ==============================================================================
# File 03: Separable Spatiotemporal Model Development
# ==============================================================================
# This file develops separable spatiotemporal models combining spatial and temporal effects
#
# DEVELOPMENT PHASE 3: Separable Spatiotemporal Models
# - Combine all spatial + temporal effect combinations (5 × 10 = 50+ models)
# - Test separable (additive) space × time interactions
# - Compare computational efficiency vs pure spatial/temporal models
# - Assess improvement from combining spatial and temporal structure using ZIP
# - Include model complexity considerations for practical decision-making
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
  library(knitr)
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

cat("=== PHASE 3: SEPARABLE SPATIOTEMPORAL MODEL DEVELOPMENT ===\n")
cat("Response variable:", RESPONSE_TYPE, "\n")
cat("Parallel processing:", PARALLEL, "\n\n")

# Set INLA options
if (PARALLEL) {
  INLA::inla.setOption(num.threads = "1:1")
}

# Enable caching for spatial data
options(tigris_use_cache = TRUE)

# ==============================================================================
# 1. LOAD DATA AND PREVIOUS RESULTS
# ==============================================================================

cat("--- Step 1: Loading Data and Previous Results ---\n")

# Check if previous phases completed
required_files <- c("outputs/data/config.rds", 
                   paste0("outputs/models/phase1_spatial_", RESPONSE_TYPE, ".rds"),
                   paste0("outputs/models/phase2a_county_temporal_", RESPONSE_TYPE, ".rds"))

missing_files <- !file.exists(required_files)
if (any(missing_files)) {
  stop("Missing required files from previous phases: ", 
       paste(required_files[missing_files], collapse = ", "))
}

# Load configuration and data
CONFIG <- readRDS("outputs/data/config.rds")
precision_matrices <- readRDS("outputs/data/precision_matrices.rds")
county_index <- readRDS("outputs/data/county_index.rds")
time_index <- readRDS("outputs/data/time_index.rds")

# Load previous phase results
phase1_results <- readRDS(paste0("outputs/models/phase1_spatial_", RESPONSE_TYPE, ".rds"))
phase2a_results <- readRDS(paste0("outputs/models/phase2a_county_temporal_", RESPONSE_TYPE, ".rds"))

# Load appropriate dataset
if (RESPONSE_TYPE == "opioid") {
  data_complete <- readRDS("outputs/data/opioid_complete.rds")
  data_splits <- readRDS("outputs/data/opioid_splits.rds")
} else {
  data_complete <- readRDS("outputs/data/stimulant_complete.rds")
  data_splits <- readRDS("outputs/data/stimulant_splits.rds")
}

cat("✓ Loaded data and previous phase results\n")
cat("✓ Phase 1 best spatial model:", phase1_results$best_model$model_name, "\n")
cat("✓ Phase 2a best temporal model:", phase2a_results$best_model$model_name, "\n")
cat("✓ Model family:", CONFIG$FAMILY, "(Zero-Inflated Poisson)\n")

# ==============================================================================
# 2. PREPARE SPATIOTEMPORAL DATA
# ==============================================================================

cat("\n--- Step 2: Preparing Spatiotemporal Data ---\n")

# Filter to counties with spatial information
train_data <- data_splits$train %>%
  filter(Residence_County %in% precision_matrices$county_names) %>%
  mutate(
    # Add spatial index
    spatial_id = match(Residence_County, precision_matrices$county_names),
    # Ensure proper indexing
    time_numeric = as.numeric(time_id),
    month_factor = factor(month(date)),
    year_factor = factor(year),
    # Create seasonal variables
    month_sin = sin(2 * pi * month(date) / 12),
    month_cos = cos(2 * pi * month(date) / 12)
  ) %>%
  arrange(spatial_id, time_id)

test_data <- data_splits$test %>%
  filter(Residence_County %in% precision_matrices$county_names) %>%
  mutate(
    spatial_id = match(Residence_County, precision_matrices$county_names),
    time_numeric = as.numeric(time_id),
    month_factor = factor(month(date)),
    year_factor = factor(year),
    month_sin = sin(2 * pi * month(date) / 12),
    month_cos = cos(2 * pi * month(date) / 12)
  ) %>%
  arrange(spatial_id, time_id)

cat("✓ Training data:", nrow(train_data), "observations\n")
cat("✓ Test data:", nrow(test_data), "observations\n")
cat("✓ Spatial units:", length(precision_matrices$county_names), "counties\n")
cat("✓ Temporal units:", max(train_data$time_id), "time periods\n")

# ==============================================================================
# 3. DEFINE SPATIOTEMPORAL MODEL COMBINATIONS
# ==============================================================================

cat("\n--- Step 3: Defining Spatiotemporal Model Combinations ---\n")

# Spatial model components (from Phase 1)
spatial_components <- list(
  none = list(formula = NULL, matrix = NULL, complexity = 0),
  adjacency = list(formula = "f(spatial_id, model='generic', Cmatrix=Q_adj, rankdef=1, constr=TRUE)", 
                   matrix = precision_matrices$adjacency, complexity = 1),
  exponential = list(formula = "f(spatial_id, model='generic', Cmatrix=Q_exp, rankdef=1, constr=TRUE)", 
                     matrix = precision_matrices$exponential, complexity = 1),
  gaussian = list(formula = "f(spatial_id, model='generic', Cmatrix=Q_gauss, rankdef=1, constr=TRUE)", 
                  matrix = precision_matrices$gaussian, complexity = 1),
  independent = list(formula = "f(spatial_id, model='iid')", 
                     matrix = NULL, complexity = 1)
)

# Temporal model components (from Phase 2a, adapted for full data)
temporal_components <- list(
  none = list(formula = NULL, complexity = 0),
  linear = list(formula = "time_numeric", complexity = 1),
  quadratic = list(formula = "time_numeric + I(time_numeric^2)", complexity = 2),
  rw1 = list(formula = "f(time_id, model='rw1')", complexity = 2),
  rw2 = list(formula = "f(time_id, model='rw2')", complexity = 3),
  ar1 = list(formula = "f(time_id, model='ar1')", complexity = 2),
  monthly = list(formula = "f(month_factor, model='iid')", complexity = 2),
  yearly = list(formula = "f(year_factor, model='iid')", complexity = 1),
  cyclic = list(formula = "month_sin + month_cos", complexity = 2)
)

# Create all model combinations
create_model_combinations <- function(spatial_components, temporal_components, max_models = 60) {
  combinations <- list()
  
  for (s_name in names(spatial_components)) {
    for (t_name in names(temporal_components)) {
      
      # Skip if both are "none" (that's just intercept only)
      if (s_name == "none" && t_name == "none") next
      
      model_name <- paste(s_name, t_name, sep = "_")
      spatial_comp <- spatial_components[[s_name]]
      temporal_comp <- temporal_components[[t_name]]
      
      # Build formula components
      fixed_parts <- c("1")
      random_parts <- c()
      
      # Add temporal fixed effects
      if (t_name %in% c("linear", "quadratic", "cyclic")) {
        fixed_parts <- c(fixed_parts, temporal_comp$formula)
      } else if (!is.null(temporal_comp$formula)) {
        random_parts <- c(random_parts, temporal_comp$formula)
      }
      
      # Add spatial random effects
      if (!is.null(spatial_comp$formula)) {
        random_parts <- c(random_parts, spatial_comp$formula)
      }
      
      # Calculate complexity score
      complexity <- spatial_comp$complexity + temporal_comp$complexity
      
      combinations[[model_name]] <- list(
        name = model_name,
        spatial_component = s_name,
        temporal_component = t_name,
        fixed_formula = paste(fixed_parts, collapse = " + "),
        random_formula = if(length(random_parts) > 0) paste(random_parts, collapse = " + ") else NULL,
        spatial_matrix = spatial_comp$matrix,
        complexity_score = complexity,
        description = paste("Spatial:", s_name, "| Temporal:", t_name)
      )
      
      # Limit total number of models for computational feasibility
      if (length(combinations) >= max_models) {
        cat("Limiting to", max_models, "model combinations for computational feasibility\n")
        break
      }
    }
    if (length(combinations) >= max_models) break
  }
  
  return(combinations)
}

# Generate model combinations
model_combinations <- create_model_combinations(spatial_components, temporal_components, max_models = 50)

cat("✓ Generated", length(model_combinations), "spatiotemporal model combinations\n")

# Show complexity distribution
complexity_table <- table(sapply(model_combinations, function(x) x$complexity_score))
cat("Model complexity distribution:\n")
print(complexity_table)

# ==============================================================================
# 4. FIT SPATIOTEMPORAL MODELS
# ==============================================================================

cat("\n--- Step 4: Fitting Spatiotemporal Models ---\n")

# Function to create formula from components
create_spatiotemporal_formula <- function(model_spec) {
  fixed_part <- model_spec$fixed_formula
  random_part <- model_spec$random_formula
  
  formula_str <- paste("distinct_patient_count ~", fixed_part)
  
  if (!is.null(random_part)) {
    formula_str <- paste(formula_str, "+", random_part)
  }
  
  formula_str <- paste(formula_str, "+ offset(log_pop_offset)")
  
  as.formula(formula_str)
}

# Function to fit a single spatiotemporal model
fit_spatiotemporal_model <- function(model_spec, train_data, model_name) {

  
  start_time <- Sys.time()
  
  # Set up precision matrices in global environment if needed
  if (!is.null(model_spec$spatial_matrix)) {
    if (model_spec$spatial_component == "adjacency") {
      assign("Q_adj", model_spec$spatial_matrix, envir = .GlobalEnv)
    } else if (model_spec$spatial_component == "exponential") {
      assign("Q_exp", model_spec$spatial_matrix, envir = .GlobalEnv)
    } else if (model_spec$spatial_component == "gaussian") {
      assign("Q_gauss", model_spec$spatial_matrix, envir = .GlobalEnv)
    }
  }
  
  # Create formula
  formula <- create_spatiotemporal_formula(model_spec)
  
  # Fit model
  tryCatch({
    result <- fit_model_with_zinb_support(formula, train_data, CONFIG$FAMILY)
  
    end_time <- Sys.time()
    runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
    # Check for successful fit
    if (is.null(result) || any(is.na(result$summary.fixed))) {
      cat("Fitting", model_name, "model (complexity:", model_spec$complexity_score, ") FAILED\n")
      return(list(
        model = NULL,
        fit_success = FALSE,
        error_message = "Model fitting failed",
        runtime = runtime
      ))
    }
  
    # Print complete success message
    cat("Fit", model_name, "model (complexity:", model_spec$complexity_score, ") ✓ (", round(runtime, 1), "s)\n")
  
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
  
    cat("Fitting", model_name, "model (complexity:", model_spec$complexity_score, ") ERROR:", e$message, "\n")
    return(list(
      model = NULL,
      fit_success = FALSE,
      error_message = e$message,
      runtime = runtime
    ))
  })
}
# Fit models in parallel
models_by_complexity <- model_combinations[order(sapply(model_combinations, function(x) x$complexity_score))]

cat("Fitting models in parallel...\n")
parallel_env <- setup_parallel() # Uses all cores minus 1
cat("Using", parallel_env$cores, "cores\n\n")

# Timing estimate
cat("Fitting", length(models_by_complexity), "models on", parallel_env$cores, "cores\n")

if (parallel_env$type == "mclapply") {
  # Mac/Linux
  spatiotemporal_results <- mclapply(names(models_by_complexity), function(model_name) {
    model_spec <- models_by_complexity[[model_name]]
    fit_spatiotemporal_model(model_spec, train_data, model_name)
  }, mc.cores = parallel_env$cores)
} else {
  # Windows
  clusterEvalQ(parallel_env$cluster, {
    library(INLA)
    library(dplyr)
    source("00_model_utilities.r")
  })
  clusterExport(parallel_env$cluster, c("train_data", "precision_matrices", "CONFIG", 
                                        "fit_spatiotemporal_model", "models_by_complexity"))
  
  spatiotemporal_results <- parLapply(parallel_env$cluster, names(models_by_complexity), function(model_name) {
    model_spec <- models_by_complexity[[model_name]]
    fit_spatiotemporal_model(model_spec, train_data, model_name)
  })
  
  stopCluster(parallel_env$cluster)
}

names(spatiotemporal_results) <- names(models_by_complexity)
successful_fits <- sum(sapply(spatiotemporal_results, function(x) x$fit_success))
total_models <- length(models_by_complexity)

cat("\n✓ Successfully fit", successful_fits, "out of", total_models, "spatiotemporal models\n")

# ==============================================================================
# 5. MODEL COMPARISON WITH COMPLEXITY CONSIDERATIONS
# ==============================================================================

cat("\n--- Step 5: Model Comparison with Complexity Analysis ---\n")

# Extract comprehensive model metrics
extract_detailed_metrics <- function(model_result, model_name) {
  if (!model_result$fit_success) {
    return(data.frame(
      model_name = model_name,
      spatial_component = model_result$model_spec$spatial_component %||% "none",
      temporal_component = model_result$model_spec$temporal_component %||% "none",
      complexity_score = model_result$model_spec$complexity_score %||% 0,
      dic = NA, waic = NA, marginal_likelihood = NA,
      mean_cpo = NA, failure_rate = NA, runtime = model_result$runtime,
      fit_success = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  
  model <- model_result$model
  spec <- model_result$model_spec
  
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
    spatial_component = spec$spatial_component %||% "none",
    temporal_component = spec$temporal_component %||% "none",
    complexity_score = spec$complexity_score,
    dic = model$dic$dic,
    waic = model$waic$waic,
    marginal_likelihood = model$mlik[1],
    mean_cpo = mean_cpo,
    failure_rate = failure_rate,
    runtime = model_result$runtime,
    fit_success = TRUE,
    stringsAsFactors = FALSE
  )
}

# Create comprehensive comparison table
comparison_metrics <- do.call(rbind, lapply(names(spatiotemporal_results), function(name) {
  extract_detailed_metrics(spatiotemporal_results[[name]], name)
}))

# Add derived metrics for model selection
comparison_metrics <- comparison_metrics %>%
  filter(fit_success == TRUE) %>%
  arrange(waic) %>%
  mutate(
	  model_type = case_when(
	    spatial_component == "none" & temporal_component != "none" ~ "single_component",
	    spatial_component != "none" & temporal_component == "none" ~ "single_component", 
	    spatial_component != "none" & temporal_component != "none" ~ "separable",
	    TRUE ~ "unknown"
	  ),
    waic_rank = rank(waic),
    dic_rank = rank(dic),
    # Complexity-adjusted metrics
    waic_per_complexity = waic / pmax(complexity_score, 1),
    # Improvement metrics
    waic_improvement = max(waic, na.rm = TRUE) - waic,
    relative_improvement = waic_improvement / max(waic, na.rm = TRUE) * 100,
    # Parsimony score (lower WAIC with lower complexity is better)
    parsimony_score = waic + 10 * complexity_score  # Penalty for complexity
  ) %>%
  arrange(waic)

cat("Top 10 models by WAIC:\n")
print(comparison_metrics %>% 
      select(model_name, spatial_component, temporal_component, complexity_score, 
             waic, dic, relative_improvement, runtime) %>%
      head(10))

# Find models representing different complexity-performance tradeoffs
best_overall <- comparison_metrics[1, ]
best_simple <- comparison_metrics %>% 
  filter(complexity_score <= 2) %>% 
  slice_min(waic, n = 1)
best_parsimony <- comparison_metrics %>% 
  slice_min(parsimony_score, n = 1)

# Find models within 3 WAIC units of the best model
models_within_3_waic <- comparison_metrics %>%
  filter(waic <= best_overall$waic + 3) %>%
  arrange(waic)

cat("\nModel selection summary:\n")
cat("Best overall WAIC:", best_overall$model_name, "(WAIC =", round(best_overall$waic, 1), 
    ", Complexity =", best_overall$complexity_score, ")\n")

if (nrow(models_within_3_waic) > 1) {
  cat("Models within 3 WAIC units (similar performance):\n")
  for (i in 1:min(5, nrow(models_within_3_waic))) {  # Show up to 5 similar models
    model <- models_within_3_waic[i, ]
    cat("  ", i, ".", model$model_name, "(WAIC =", round(model$waic, 1), 
        ", Complexity =", model$complexity_score, ")\n")
  }
} else {
  cat("No other models within 3 WAIC units of the best model\n")
}

cat("Best simple model (≤2 complexity):", best_simple$model_name, "(WAIC =", round(best_simple$waic, 1), 
    ", Complexity =", best_simple$complexity_score, ")\n")
cat("Best parsimony (complexity-adjusted):", best_parsimony$model_name, "(WAIC =", round(best_parsimony$waic, 1), 
    ", Complexity =", best_parsimony$complexity_score, ")\n")

# Compare with individual models on same data structure
# Find temporal-only and spatial-only models from our results
temporal_only_models <- comparison_metrics %>% filter(spatial_component == "none")
spatial_only_models <- comparison_metrics %>% filter(temporal_component == "none")

best_temporal_only <- temporal_only_models %>% slice_min(waic, n = 1)
best_spatial_only <- spatial_only_models %>% slice_min(waic, n = 1)

# Calculate meaningful improvements
if (nrow(best_temporal_only) > 0 && nrow(best_spatial_only) > 0) {
  temporal_improvement <- best_temporal_only$waic - best_overall$waic
  spatial_improvement <- best_spatial_only$waic - best_overall$waic
  
  cat("\nMeaningful improvement analysis (same dataset comparisons):\n")
  cat("Best temporal-only model:", best_temporal_only$model_name, "WAIC =", round(best_temporal_only$waic, 1), "\n")
  cat("Best spatial-only model:", best_spatial_only$model_name, "WAIC =", round(best_spatial_only$waic, 1), "\n")
  cat("Best combined model:", best_overall$model_name, "WAIC =", round(best_overall$waic, 1), "\n")
  cat("WAIC improvement over temporal-only:", round(temporal_improvement, 1), "units\n")
  cat("WAIC improvement over spatial-only:", round(spatial_improvement, 1), "units\n")
  
  # Calculate response-scale improvements
  response_scale_improvements <- NULL
  
  # Extract fitted values from the three models for response-scale comparison
  best_overall_result <- spatiotemporal_results[[best_overall$model_name]]
  temporal_only_result <- spatiotemporal_results[[best_temporal_only$model_name]]
  spatial_only_result <- spatiotemporal_results[[best_spatial_only$model_name]]
  
  if (!is.null(best_overall_result$model) && 
      !is.null(temporal_only_result$model) && 
      !is.null(spatial_only_result$model)) {
    
    # Get fitted values
    best_fitted <- best_overall_result$model$summary.fitted.values$mean
    temporal_fitted <- temporal_only_result$model$summary.fitted.values$mean
    spatial_fitted <- spatial_only_result$model$summary.fitted.values$mean
    observed <- train_data$distinct_patient_count
    
    # Calculate MAE for each model
    best_mae <- mean(abs(observed - best_fitted))
    temporal_mae <- mean(abs(observed - temporal_fitted))
    spatial_mae <- mean(abs(observed - spatial_fitted))
    
    # Calculate RMSE for each model
    best_rmse <- sqrt(mean((observed - best_fitted)^2))
    temporal_rmse <- sqrt(mean((observed - temporal_fitted)^2))
    spatial_rmse <- sqrt(mean((observed - spatial_fitted)^2))
    
    # Calculate percentage improvements
    mae_improvement_temporal <- round((1 - best_mae/temporal_mae) * 100, 1)
    mae_improvement_spatial <- round((1 - best_mae/spatial_mae) * 100, 1)
    rmse_improvement_temporal <- round((1 - best_rmse/temporal_rmse) * 100, 1)
    rmse_improvement_spatial <- round((1 - best_rmse/spatial_rmse) * 100, 1)
    
    response_scale_improvements <- list(
      mae_improvement_temporal = mae_improvement_temporal,
      mae_improvement_spatial = mae_improvement_spatial,
      rmse_improvement_temporal = rmse_improvement_temporal,
      rmse_improvement_spatial = rmse_improvement_spatial,
      best_mae = best_mae,
      temporal_mae = temporal_mae,
      spatial_mae = spatial_mae
    )
    
    cat("\nPrediction accuracy improvements (response scale):\n")
    cat("  Over temporal-only - MAE:", mae_improvement_temporal, "% better, RMSE:", 
        rmse_improvement_temporal, "% better\n")
    cat("  Over spatial-only - MAE:", mae_improvement_spatial, "% better, RMSE:", 
        rmse_improvement_spatial, "% better\n")
    cat("  Actual MAE values:\n")
    cat("    Best combined:", round(best_mae, 3), "deaths per county-month\n")
    cat("    Temporal-only:", round(temporal_mae, 3), "deaths per county-month\n")
    cat("    Spatial-only:", round(spatial_mae, 3), "deaths per county-month\n")
  }
  
  # Interpret WAIC differences
  cat("\nWAIC Evidence Strength:\n")
  if (temporal_improvement > 20) {
    cat("  vs. temporal-only: Very strong evidence (>20 WAIC units)\n")
  } else if (temporal_improvement > 10) {
    cat("  vs. temporal-only: Strong evidence (>10 WAIC units)\n")
  } else if (temporal_improvement > 6) {
    cat("  vs. temporal-only: Positive evidence (>6 WAIC units)\n")
  } else {
    cat("  vs. temporal-only: Weak evidence (<6 WAIC units)\n")
  }
  
  if (spatial_improvement > 20) {
    cat("  vs. spatial-only: Very strong evidence (>20 WAIC units)\n")
  } else if (spatial_improvement > 10) {
    cat("  vs. spatial-only: Strong evidence (>10 WAIC units)\n")
  } else if (spatial_improvement > 6) {
    cat("  vs. spatial-only: Positive evidence (>6 WAIC units)\n")
  } else {
    cat("  vs. spatial-only: Weak evidence (<6 WAIC units)\n")
  }
  
  cat("\nNote: WAIC improvements are on log-likelihood scale, not percentages\n")
  cat("      Response-scale metrics (MAE, RMSE) show actual prediction improvements\n")
  
} else {
  cat("\nSome comparison models not available in this run\n")
}

# ==============================================================================
# 6. VISUALIZATION OF RESULTS
# ==============================================================================

cat("\n--- Step 6: Creating Visualizations ---\n")

# DIAGNOSTIC: Verify data before plotting
cat("Data check before plotting:\n")
cat("  adjacency_ar1 WAIC in comparison_metrics:", 
    comparison_metrics$waic[comparison_metrics$model_name == "adjacency_ar1"], "\n")
cat("  Best overall WAIC:", min(comparison_metrics$waic), "\n")
cat("  Number of unique models:", length(unique(comparison_metrics$model_name)), "\n")

# Check for duplicates
if (any(duplicated(comparison_metrics$model_name))) {
  cat("  WARNING: Duplicate model names found!\n")
  dup_models <- comparison_metrics$model_name[duplicated(comparison_metrics$model_name)]
  cat("  Duplicates:", paste(dup_models, collapse=", "), "\n")
  
  # Remove duplicates, keeping best WAIC for each model
  comparison_metrics <- comparison_metrics %>%
    group_by(model_name) %>%
    slice_min(waic, n = 1) %>%
    ungroup()
}

# 1. Model component heatmap (green/brown/red style)
if (successful_fits > 10) {
# Create improved heatmap data with better WAIC scaling
heatmap_data <- comparison_metrics %>%
  group_by(spatial_component, temporal_component) %>%
  summarise(
    mean_waic = mean(waic, na.rm = TRUE),
    best_waic = min(waic, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    # Add WAIC labels for the heatmap
    waic_label = round(best_waic, 0)
  )

# Function to format axis labels
format_labels <- function(x) {
  formatted <- gsub("_", " ", x)
  tools::toTitleCase(formatted)
}

# Heatmap with continuous color scale
component_heatmap <- ggplot(heatmap_data, 
                           aes(x = temporal_component, y = spatial_component)) +
  # Use continuous WAIC values for fill
  geom_tile(aes(fill = best_waic), color = "white", linewidth = 0.5) +
  # Add WAIC values as text
  geom_text(aes(label = waic_label), color = "white", size = 3.5, fontface = "bold") +
  # Use scale_fill_viridis_c for continuous values
  scale_fill_gradientn(
          name = "WAIC\n(lower better)",
          colors = c("#1a7a1a",    # Dark green for best (~11000)
                     "#5aae61",    # Medium green (~11100)
					 "#8b7355",    # Medium brown (~11150)
                     "#cd9b45",    # Golden brown (~11200)
                     "#fdae61",    # Orange (~11300)
                     "#f46d43",    # Orange-red (~11500)
                     "#d73027"),   # Red for worst (12000+)
          values = c(0, 0.05, 0.1, 0.2, 0.4, 0.7, 1),  # Compress green end, expand red
          limits = c(min(heatmap_data$best_waic), max(heatmap_data$best_waic))
        ) +
  # Better axis labels
  scale_x_discrete(labels = format_labels) +
  scale_y_discrete(labels = format_labels) +
  labs(
    title = "Spatial × Temporal Component Performance",
    subtitle = paste("WAIC by component combination -", RESPONSE_TYPE, "deaths (ZIP)"),
    x = "Temporal Component", 
    y = "Spatial Component",
    caption = "Lower WAIC = better model • Numbers show exact WAIC values"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    panel.grid = element_blank()
  )
} else {
component_heatmap <- NULL
}

# 2. Top models comparison with compressed scale and WAIC thresholds
top_models <- comparison_metrics %>% head(15)  # Show more models to see differences

# Calculate WAIC differences from best model
best_waic <- min(top_models$waic)
top_models <- top_models %>%
  mutate(
    waic_diff = waic - best_waic,
    # Color categories based on WAIC differences
    waic_category = case_when(
      waic_diff == 0 ~ "Best Model",
      waic_diff <= 3 ~ "Equivalent (≤3 WAIC)",
      waic_diff <= 7 ~ "Moderate Diff (3-7 WAIC)", 
      waic_diff <= 20 ~ "Large Diff (7-20 WAIC)",
      TRUE ~ "Very Large Diff (>20 WAIC)"
    ),
    waic_category = factor(waic_category, levels = c("Best Model", "Equivalent (≤3 WAIC)", 
                                                    "Moderate Diff (3-7 WAIC)", 
                                                    "Large Diff (7-20 WAIC)", 
                                                    "Very Large Diff (>20 WAIC)"))
  )

# Dynamic scale compression - focus on the actual range of top models
waic_range <- max(top_models$waic) - min(top_models$waic)
y_min <- best_waic - (waic_range * 0.1)  # Add 10% buffer below
y_max <- max(top_models$waic) + (waic_range * 0.2)  # Add 20% buffer above for annotations

top_models_plot <- ggplot(top_models, aes(x = reorder(model_name, -waic), y = waic)) +
  geom_col(aes(fill = waic_category), alpha = 0.8, width = 0.7) +
  geom_text(aes(label = round(waic, 0)), vjust = -0.5, size = 3) +
  # Add horizontal reference lines for WAIC thresholds
  geom_hline(yintercept = best_waic + 3, linetype = "dashed", color = "darkgray", alpha = 0.7) +
  geom_hline(yintercept = best_waic + 7, linetype = "dashed", color = "darkgray", alpha = 0.7) +
  geom_hline(yintercept = best_waic + 20, linetype = "dashed", color = "darkgray", alpha = 0.7) +
  # Custom colors for different WAIC categories
  scale_fill_manual(
    name = "Model Performance",
    values = c("Best Model" = "#2166ac", 
               "Equivalent (≤3 WAIC)" = "#5aae61", 
               "Moderate Diff (3-7 WAIC)" = "#fee08b", 
               "Large Diff (7-20 WAIC)" = "#fdae61", 
               "Very Large Diff (>20 WAIC)" = "#d73027")
  ) +
  # Compressed y-axis scale
  coord_cartesian(ylim = c(y_min, y_max)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
  theme_minimal() +
  labs(
    title = "Top Spatiotemporal Models (ZIP) - Compressed Scale",
    subtitle = paste("WAIC differences emphasized for", RESPONSE_TYPE, "deaths"),
    x = "Model", y = "WAIC (lower is better)",
    caption = "Dashed lines: +3, +7, +20 WAIC thresholds • Colors indicate performance tiers"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "bottom",
    plot.caption = element_text(size = 9, hjust = 0),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  # Add threshold annotations
  annotate("text", x = 1, y = best_waic + 3, label = "+3 WAIC", 
           hjust = 0, vjust = -0.5, size = 3, color = "darkgray") +
  annotate("text", x = 1, y = best_waic + 7, label = "+7 WAIC", 
           hjust = 0, vjust = -0.5, size = 3, color = "darkgray") +
  annotate("text", x = 1, y = best_waic + 20, label = "+20 WAIC", 
           hjust = 0, vjust = -0.5, size = 3, color = "darkgray")

cat("✓ Created visualization plots\n")

# ==============================================================================
# 7. SAVE RESULTS
# ==============================================================================

cat("\n--- Step 7: Saving Results ---\n")

# Create comprehensive results summary
phase3_results <- list(
  response_type = RESPONSE_TYPE,
  model_family = CONFIG$FAMILY,
  model_combinations = model_combinations,
  model_fits = spatiotemporal_results,
  comparison_metrics = comparison_metrics,
  best_models = list(
    overall = best_overall,
    simple = best_simple,
    parsimony = best_parsimony,
    within_3_waic = models_within_3_waic
  ),
  improvement_analysis = list(
    temporal_only_best = if (nrow(best_temporal_only) > 0) best_temporal_only else NULL,
    spatial_only_best = if (nrow(best_spatial_only) > 0) best_spatial_only else NULL,
    combined_best = best_overall,
    temporal_improvement = if (exists("temporal_improvement")) temporal_improvement else NA,
    spatial_improvement = if (exists("spatial_improvement")) spatial_improvement else NA
  ),
  complexity_analysis = complexity_table,
  n_models_tested = total_models,
  n_successful_fits = successful_fits,
  fitting_time = Sys.time()
)

# Save detailed results
saveRDS(phase3_results, paste0("outputs/models/phase3_spatiotemporal_", RESPONSE_TYPE, ".rds"))

# Save comparison metrics as CSV
write.csv(comparison_metrics, 
          paste0("outputs/models/phase3_comparison_", RESPONSE_TYPE, ".csv"), 
          row.names = FALSE)

# Save plots
plots_to_save <- list(
  heatmap = component_heatmap,
  top_models = top_models_plot
)

for (plot_name in names(plots_to_save)) {
  if (!is.null(plots_to_save[[plot_name]])) {
    ggsave(paste0("outputs/plots/phase3_", plot_name, "_", RESPONSE_TYPE, ".png"), 
           plots_to_save[[plot_name]], 
           width = 12, height = 8, dpi = 300, bg = "white")
  }
}

# Create summary for next phase
summary_for_phase4 <- list(
  best_separable_model = best_overall$model_name,
  separable_model_available = TRUE,
  model_type = comparison_metrics$model_type,
  n_successful_separable_fits = successful_fits,
  temporal_improvement = if (exists("temporal_improvement")) round(temporal_improvement, 1) else NA,
  spatial_improvement = if (exists("spatial_improvement")) round(spatial_improvement, 1) else NA,
  complexity_vs_performance_available = TRUE,
  model_family_used = CONFIG$FAMILY,
  data_ready_for_phase4 = TRUE
)

saveRDS(summary_for_phase4, paste0("outputs/data/phase3_to_phase4_", RESPONSE_TYPE, ".rds"))

cat("✓ Results saved to outputs/models/\n")
cat("✓ Plots saved to outputs/plots/\n")

cat("\n=== PHASE 3 COMPLETE ===\n")
cat("Model family used:", CONFIG$FAMILY, "(Zero-Inflated Poisson)\n")
cat("Successful spatiotemporal models:", successful_fits, "/", total_models, "\n")
cat("Best overall model:", best_overall$model_name, "(WAIC =", round(best_overall$waic, 1), ")\n")

if (nrow(models_within_3_waic) > 1) {
  cat("Models within 3 WAIC units:", nrow(models_within_3_waic), "models\n")
} else {
  cat("No other models within 3 WAIC units\n")
}

cat("Best simple model:", best_simple$model_name, "(WAIC =", round(best_simple$waic, 1), ")\n")
if (exists("temporal_improvement") && exists("spatial_improvement")) {
  cat("Improvement from combination: temporal ", round(temporal_improvement, 1), 
      " units, spatial ", round(spatial_improvement, 1), " units\n")
} else {
  cat("Individual improvement components calculated above\n")
}
cat("Ready for Phase 4: Full Spatiotemporal Interaction Models\n")

# Clean up global environment
rm(list = ls(pattern = "^Q_"))

# ==============================================================================
# END OF PHASE 3
# ==============================================================================