# ==============================================================================
# File 01: Spatial Model Development for WV County Drug Death Analysis
# ==============================================================================
# This file develops and compares different spatial model structures using R-INLA
#
# DEVELOPMENT PHASE 1: Spatial Model Development
# - Binary adjacency spatial models (Queen/Rook contiguity)
# - Distance-weighted spatial models (exponential/Gaussian decay)
# - Multiple precision matrix formulations and spatial random effects
# - Test spatial structures independently using Zero-Inflated Poisson
# - Compare spatial model performance and create spatial effect visualizations
# ==============================================================================

# Clear environment and load libraries
rm(list = ls())
suppressPackageStartupMessages({
  library(INLA)
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(viridis)
  library(gridExtra)
  library(Matrix)
})

# Load shared utilities
source("00_model_utilities.r")

# Get response type from environment variable (set by master script)
RESPONSE_TYPE <- Sys.getenv("ANALYSIS_RESPONSE_TYPE", unset = "opioid")
PARALLEL <- as.logical(Sys.getenv("ANALYSIS_PARALLEL", unset = "TRUE"))

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
precision_matrices <- readRDS("outputs/data/precision_matrices.rds")
county_index <- readRDS("outputs/data/county_index.rds")
time_index <- readRDS("outputs/data/time_index.rds")
county_coords <- readRDS("outputs/data/county_coords.rds")

cat("=== PHASE 1: SPATIAL MODEL DEVELOPMENT ===\n")
cat("Response variable:", RESPONSE_TYPE, "\n")
cat("Model family:", CONFIG$FAMILY, "(Zero-Inflated Poisson for zero-inflation)\n")
cat("Parallel processing:", PARALLEL, "\n")
cat("Note: Spatial boundaries will be loaded from US Census via tigris package\n\n")

# Set INLA options
if (PARALLEL) {
  INLA::inla.setOption(num.threads = "1:1")  # Let INLA handle threading
}

# Enable caching for tigris to speed up boundary downloads
options(tigris_use_cache = TRUE)

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

# Check spatial data availability
if (length(precision_matrices$county_names) == 0) {
  stop("No spatial data available. Check county name matching in Phase 0.")
}
cat("✓ Spatial data available for", length(precision_matrices$county_names), "counties\n")

# ==============================================================================
# 2. PREPARE SPATIAL DATA FOR MODELING
# ==============================================================================

cat("\n--- Step 2: Preparing Spatial Data for Modeling ---\n")

# Filter data to counties with spatial information
train_data <- data_splits$train %>%
  filter(Residence_County %in% precision_matrices$county_names) %>%
  arrange(county_id, time_id)

test_data <- data_splits$test %>%
  filter(Residence_County %in% precision_matrices$county_names) %>%
  arrange(county_id, time_id)

cat("✓ Training data after spatial filtering:", nrow(train_data), "observations\n")
cat("✓ Test data after spatial filtering:", nrow(test_data), "observations\n")

# Create spatial index mapping for filtered data
spatial_counties <- precision_matrices$county_names
n_spatial <- length(spatial_counties)
spatial_index_map <- data.frame(
  county_name = spatial_counties,
  spatial_id = 1:n_spatial
)

# Add spatial indices to data
add_spatial_index <- function(data) {
  data %>%
    left_join(spatial_index_map, by = c("Residence_County" = "county_name")) %>%
    filter(!is.na(spatial_id))
}

train_data <- add_spatial_index(train_data)
test_data <- add_spatial_index(test_data)

cat("✓ Added spatial indices to data\n")

# ==============================================================================
# 3. DEFINE SPATIAL MODEL SPECIFICATIONS
# ==============================================================================

cat("\n--- Step 3: Defining Spatial Model Specifications ---\n")

# Model specifications to test
spatial_models <- list(
  
  # Baseline: No spatial structure
  baseline = list(
    name = "Baseline (No Spatial)",
    formula_parts = list(
      response = "distinct_patient_count",
      fixed = "1",
      random = NULL,
      offset = "log_pop_offset"
    ),
    precision_matrix = NULL,
    description = "Intercept only with population offset"
  ),
  
  # Binary adjacency (standard CAR)
  adjacency = list(
    name = "Binary Adjacency CAR",
    formula_parts = list(
      response = "distinct_patient_count",
      fixed = "1",
      random = "f(spatial_id, model='generic', Cmatrix=Q_adj, rankdef=1, constr=TRUE)",
      offset = "log_pop_offset"
    ),
    precision_matrix = precision_matrices$adjacency,
    description = "Conditional autoregressive model with binary adjacency"
  ),
  
  # Distance-weighted (exponential)
  exponential = list(
    name = "Exponential Distance Weights",
    formula_parts = list(
      response = "distinct_patient_count",
      fixed = "1", 
      random = "f(spatial_id, model='generic', Cmatrix=Q_exp, rankdef=1, constr=TRUE)",
      offset = "log_pop_offset"
    ),
    precision_matrix = precision_matrices$exponential,
    description = "CAR model with exponential distance decay"
  ),
  
  # Distance-weighted (Gaussian)
  gaussian = list(
    name = "Gaussian Distance Weights",
    formula_parts = list(
      response = "distinct_patient_count",
      fixed = "1",
      random = "f(spatial_id, model='generic', Cmatrix=Q_gauss, rankdef=1, constr=TRUE)",
      offset = "log_pop_offset"
    ),
    precision_matrix = precision_matrices$gaussian,
    description = "CAR model with Gaussian distance decay"
  ),
  
  # Independent spatial effects
  independent = list(
    name = "Independent Spatial Effects",
    formula_parts = list(
      response = "distinct_patient_count",
      fixed = "1",
      random = "f(spatial_id, model='iid')",
      offset = "log_pop_offset"
    ),
    precision_matrix = NULL,
    description = "Independent spatial random effects"
  )
)

cat("✓ Defined", length(spatial_models), "spatial model specifications:\n")
for (i in 1:length(spatial_models)) {
  cat("  ", i, ".", spatial_models[[i]]$name, "\n")
}

# ==============================================================================
# 4. FIT SPATIAL MODELS
# ==============================================================================

cat("\n--- Step 4: Fitting Spatial Models ---\n")

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

# Function to fit a single spatial model
fit_spatial_model <- function(model_spec, train_data, model_name) {
  cat("Fitting", model_name, "...\n")
  
  start_time <- Sys.time()
  
  # Set up precision matrix in global environment if needed
  if (!is.null(model_spec$precision_matrix)) {
    matrix_name <- paste0("Q_", gsub("[^a-zA-Z]", "", tolower(model_name)))
    assign(matrix_name, model_spec$precision_matrix, envir = .GlobalEnv)
    
    # Update formula to use the correct matrix name
    model_spec$formula_parts$random <- gsub("Q_\\w+", matrix_name, model_spec$formula_parts$random)
  }
  
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

# Fit all spatial models
spatial_results <- list()

for (model_name in names(spatial_models)) {
  result <- fit_spatial_model(spatial_models[[model_name]], train_data, model_name)
  spatial_results[[model_name]] <- result
}

# Count successful fits
successful_fits <- sum(sapply(spatial_results, function(x) x$fit_success))
cat("\n✓ Successfully fit", successful_fits, "out of", length(spatial_models), "spatial models\n")

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
comparison_metrics <- do.call(rbind, lapply(names(spatial_results), function(name) {
  extract_metrics(spatial_results[[name]], name)
}))

# Sort by WAIC (lower is better)
comparison_metrics <- comparison_metrics %>%
  arrange(waic) %>%
  mutate(
    waic_rank = ifelse(fit_success, rank(waic, na.last = TRUE), NA),
    dic_rank = ifelse(fit_success, rank(dic, na.last = TRUE), NA)
  )

cat("Model Comparison Results:\n")
print(comparison_metrics)

# Find best model
best_spatial <- comparison_metrics %>%
  filter(fit_success == TRUE) %>%
  slice_min(waic, n = 1)

if (nrow(best_spatial) > 0) {
  cat("\n✓ Best spatial model:", best_spatial$model_name, "(WAIC =", round(best_spatial$waic, 2), ")\n")
} else {
  cat("\n❌ No spatial models fit successfully\n")
}

# ==============================================================================
# 6. SPATIAL EFFECT VISUALIZATION
# ==============================================================================

cat("\n--- Step 6: Creating Spatial Effect Visualizations ---\n")

# Function to extract spatial effects
extract_spatial_effects <- function(model_result, model_name) {
  if (!model_result$fit_success) return(NULL)
  
  model <- model_result$model
  
  # Check if model has spatial effects
  if (!"spatial_id" %in% names(model$summary.random)) {
    return(data.frame(
      model_name = model_name,
      spatial_id = 1:n_spatial,
      county_name = spatial_counties,
      spatial_effect = 0,
      spatial_lower = 0,
      spatial_upper = 0
    ))
  }
  
  spatial_summary <- model$summary.random$spatial_id
  
  data.frame(
    model_name = model_name,
    spatial_id = 1:nrow(spatial_summary),
    county_name = spatial_counties,
    spatial_effect = spatial_summary$mean,
    spatial_lower = spatial_summary$`0.025quant`,
    spatial_upper = spatial_summary$`0.975quant`
  )
}

# Extract spatial effects for all models
spatial_effects_list <- lapply(names(spatial_results), function(name) {
  extract_spatial_effects(spatial_results[[name]], name)
})

# Combine and filter out NULL results
spatial_effects <- do.call(rbind, spatial_effects_list[!sapply(spatial_effects_list, is.null)])

# Create spatial effects plot
if (nrow(spatial_effects) > 0) {
  # Add coordinates for mapping
  spatial_effects <- spatial_effects %>%
    left_join(county_coords, by = "county_name")
  
  # Load WV boundaries for plotting context (do this once, outside the loop)
  wv_state <- NULL
  wv_counties_sf <- NULL
  
  cat("Loading WV boundaries for enhanced plotting...\n")
  tryCatch({
    if (require(tigris, quietly = TRUE)) {
      wv_state <- tigris::states(year = 2020) %>%
        filter(NAME == "West Virginia")
      wv_counties_sf <- tigris::counties(state = "WV", year = 2020) %>%
        mutate(county_name = tolower(NAME))
      cat("✓ WV boundaries loaded successfully\n")
    }
  }, error = function(e) {
    cat("Note: Could not load state boundaries\n")
  })
  
  # Create plot for each successful model
  plot_list <- list()
  successful_models <- unique(spatial_effects$model_name)
  
  for (model_name in successful_models) {
    model_effects <- spatial_effects %>% filter(model_name == !!model_name)
    
    p <- ggplot()
    
    # Add state and county boundaries if available
    if (!is.null(wv_state) && !is.null(wv_counties_sf)) {
      p <- p + 
        geom_sf(data = wv_state, fill = "transparent", color = "darkgray", linewidth = 0.8) +
        geom_sf(data = wv_counties_sf, fill = "transparent", color = "lightgray", linewidth = 0.3) +
        coord_sf(expand = FALSE)
    }
    
    # Add spatial effects points
    p <- p +
      geom_point(data = model_effects, aes(x = lon, y = lat, color = spatial_effect), 
                 size = 3, alpha = 0.8) +
      scale_color_viridis_c(name = "Spatial\nEffect") +
      theme_minimal() +
      labs(
        title = paste("Spatial Effects:", model_name),
        subtitle = paste("Response:", RESPONSE_TYPE, "using ZIP"),
        x = "Longitude", y = "Latitude"
      ) +
      theme(
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 10),
        plot.subtitle = element_text(size = 9),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)
      )
    
    plot_list[[model_name]] <- p
  }
  
  # Combine plots
  if (length(plot_list) > 1) {
    combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 2))
  } else {
    combined_plot <- plot_list[[1]]
  }
  
  cat("✓ Created spatial effects visualization\n")
  
  # Create and save separate plot for best spatial model
  if (nrow(best_spatial) > 0) {
    best_model_name <- best_spatial$model_name
    best_model_effects <- spatial_effects %>% filter(model_name == best_model_name)
    
    cat("Creating separate plot for best model:", best_model_name, "\n")
    
    # Create enhanced plot for best model
    best_plot <- ggplot()
    
    # Add state and county boundaries if available
    if (!is.null(wv_state) && !is.null(wv_counties_sf)) {
      best_plot <- best_plot + 
        geom_sf(data = wv_state, fill = "transparent", color = "darkgray", linewidth = 1.0) +
        geom_sf(data = wv_counties_sf, fill = "transparent", color = "lightgray", linewidth = 0.4) +
        coord_sf(expand = FALSE)
    }
    
    # Add spatial effects points with enhanced styling
    best_plot <- best_plot +
      geom_point(data = best_model_effects, aes(x = lon, y = lat, color = spatial_effect), 
                 size = 4, alpha = 0.9, stroke = 0.2) +
      scale_color_viridis_c(name = "Spatial\nEffect", option = "plasma") +
      theme_minimal() +
      labs(
        title = paste("Best Spatial Model:", best_model_name),
        subtitle = paste("WV", RESPONSE_TYPE, "deaths - Spatial random effects (Zero-Inflated Poisson, WAIC =", round(best_spatial$waic, 1), ")"),
        x = "Longitude", y = "Latitude",
        caption = paste("• Yellow/Pink: Higher than expected rates  • Blue/Purple: Lower than expected rates",
                       "\n• ZIP model handles", round(sum(train_data$distinct_patient_count == 0)/nrow(train_data)*100, 1), "% zero observations")
      ) +
      theme(
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        plot.caption = element_text(size = 10, hjust = 0, margin = margin(t = 10)),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        panel.grid = element_line(linewidth = 0.2, color = "lightgray"),
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    # Save best model plot
    ggsave(paste0("outputs/plots/phase1_best_spatial_model_", RESPONSE_TYPE, ".png"), 
           best_plot, width = 12, height = 9, dpi = 300, bg = "white")
    
    cat("✓ Best spatial model plot saved as: phase1_best_spatial_model_", RESPONSE_TYPE, ".png\n")
    
    # === ADD CHOROPLETH MAP ===
    cat("Creating choropleth map for best spatial model...\n")
    
    tryCatch({
      if (!is.null(wv_counties_sf)) {
        # Join spatial effects with county boundaries
        wv_counties_enhanced <- wv_counties_sf %>%
          left_join(best_model_effects, by = "county_name")
        
        # Create choropleth map
        choropleth_plot <- ggplot(wv_counties_enhanced) +
          # Fill counties with spatial effects
          geom_sf(aes(fill = spatial_effect), color = "white", size = 0.3) +
          # Add state outline
          geom_sf(data = wv_state, fill = "transparent", color = "black", size = 1.0) +
          # Use plasma color scheme (better for choropleth)
          scale_fill_viridis_c(
            name = "Spatial\nEffect",
            option = "plasma",
            na.value = "lightgray"
          ) +
          # Remove axes and grid for clean map
          theme_void() +
          # Add titles and styling
          labs(
            title = paste("Spatial Effects by County:", best_model_name),
            subtitle = paste("WV", RESPONSE_TYPE, "deaths - ZIP model spatial effects"),
            caption = paste("WAIC =", round(best_spatial$waic, 1),
                           "• Purple: Lower than expected • Yellow/Pink: Higher than expected")
          ) +
          theme(
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 14, hjust = 0.5),
            plot.caption = element_text(size = 11, hjust = 0.5, margin = margin(t = 10)),
            legend.position = "right",
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            plot.background = element_rect(fill = "white", color = NA)
          )
        
        # Save the choropleth map
        ggsave(paste0("outputs/plots/phase1_choropleth_spatial_", RESPONSE_TYPE, ".png"), 
               choropleth_plot, 
               width = 12, height = 9, dpi = 300, bg = "white")
        
        cat("✓ Choropleth map saved as: phase1_choropleth_spatial_", RESPONSE_TYPE, ".png\n")
        
        # Print spatial effects summary
        cat("\nSpatial Effects Summary:\n")
        cat("Highest spatial effect:", round(max(best_model_effects$spatial_effect, na.rm = TRUE), 3), 
            "in", best_model_effects$county_name[which.max(best_model_effects$spatial_effect)], "County\n")
        cat("Lowest spatial effect:", round(min(best_model_effects$spatial_effect, na.rm = TRUE), 3), 
            "in", best_model_effects$county_name[which.min(best_model_effects$spatial_effect)], "County\n")
        
      } else {
        cat("County boundaries not available for choropleth map\n")
      }
    }, error = function(e) {
      cat("Error creating choropleth map:", e$message, "\n")
    })
  }
  
} else {
  cat("No spatial effects to visualize\n")
  combined_plot <- NULL
}

# ==============================================================================
# 7. SAVE RESULTS
# ==============================================================================

cat("\n--- Step 7: Saving Results ---\n")

# Create results summary
phase1_results <- list(
  response_type = RESPONSE_TYPE,
  model_family = CONFIG$FAMILY,
  model_specifications = spatial_models,
  model_fits = spatial_results,
  comparison_metrics = comparison_metrics,
  best_model = if (nrow(best_spatial) > 0) best_spatial else NULL,
  spatial_effects = if (!is.null(spatial_effects) && nrow(spatial_effects) > 0) spatial_effects else NULL,
  n_counties = n_spatial,
  n_train_obs = nrow(train_data),
  n_test_obs = nrow(test_data),
  fitting_time = Sys.time()
)

# Save detailed results
saveRDS(phase1_results, paste0("outputs/models/phase1_spatial_", RESPONSE_TYPE, ".rds"))
cat("✓ Phase 1 results (.rds) saved\n")

# Save comparison metrics as CSV for easy viewing
write.csv(comparison_metrics, 
          paste0("outputs/models/phase1_comparison_", RESPONSE_TYPE, ".csv"), 
          row.names = FALSE)
cat("✓ Comparison metrics CSV saved\n")

# Save spatial effects if available
if (!is.null(spatial_effects) && nrow(spatial_effects) > 0) {
  write.csv(spatial_effects, 
            paste0("outputs/models/phase1_spatial_effects_", RESPONSE_TYPE, ".csv"), 
            row.names = FALSE)
}

# Save plots
if (!is.null(combined_plot)) {
  ggsave(paste0("outputs/plots/phase1_spatial_effects_", RESPONSE_TYPE, ".png"), 
         combined_plot, width = 15, height = 10, dpi = 300)
}

cat("✓ Results saved to outputs/models/\n")
cat("✓ Plots saved to outputs/plots/\n")

# Create summary for next phase
summary_for_phase2 <- list(
  best_spatial_model = if (nrow(best_spatial) > 0) best_spatial$model_name else "baseline",
  spatial_model_available = nrow(best_spatial) > 0,
  n_successful_fits = successful_fits,
  data_ready_for_phase2 = TRUE,
  model_family_used = CONFIG$FAMILY
)

saveRDS(summary_for_phase2, paste0("outputs/data/phase1_to_phase2_", RESPONSE_TYPE, ".rds"))

cat("\n=== PHASE 1 COMPLETE ===\n")
cat("Model family used:", CONFIG$FAMILY, "(Zero-Inflated Poisson)\n")
cat("Successful spatial models:", successful_fits, "/", length(spatial_models), "\n")
if (nrow(best_spatial) > 0) {
  cat("Best model:", best_spatial$model_name, "\n")
  cat("Best WAIC:", round(best_spatial$waic, 2), "\n")
}
cat("✓ Ready for Phase 2: Temporal Model Development\n")

# Clean up global environment
rm(list = ls(pattern = "^Q_"))

# ==============================================================================
# END OF PHASE 1
# ==============================================================================