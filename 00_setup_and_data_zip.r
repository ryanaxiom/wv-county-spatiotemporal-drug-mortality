# ==============================================================================
# File 00: Setup with K-Months Ahead Validation
# ==============================================================================
# This file creates both temporal holdout and k-months ahead validation splits
# for proper model selection and final evaluation
#
# VALIDATION FRAMEWORK:
# 1. Temporal holdout: Last 12 months for final evaluation  
# 2. K-months ahead: Rolling windows within training period for model selection
# ==============================================================================

# Clear environment and set options
rm(list = ls())
options(stringsAsFactors = FALSE)

# Load required libraries
suppressPackageStartupMessages({
  library(INLA)
  library(sf)
  library(spdep)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(maps)
  library(Matrix)
  library(readr)
  library(zoo)
})

# Enable caching for tigris to speed up boundary downloads
options(tigris_use_cache = TRUE)

# Create directory structure for outputs
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/data", showWarnings = FALSE)
dir.create("outputs/models", showWarnings = FALSE)
dir.create("outputs/diagnostics", showWarnings = FALSE)
dir.create("outputs/plots", showWarnings = FALSE)

cat("=== PHASE 0: DATA INFRASTRUCTURE WITH K-MONTHS VALIDATION ===\n")

# ==============================================================================
# FIXED CONFIGURATION PARAMETERS
# ==============================================================================

CONFIG <- list(
  # Data parameters
  START_YEAR = 2015,
  END_YEAR = 2023,
  START_MONTH = 1,
  END_MONTH = 12,
  
  # Response variable settings
  PRIMARY_RESPONSE = "opioid_death",
  COUNT_VAR = "distinct_patient_count",
  POPULATION_VAR = "population",
  
  # Model settings
  FAMILY = "zeroinflatedpoisson1",
  USE_OFFSET = TRUE,
  
  # Validation parameters (6-month holdout to allow final assessment period learning)
  HOLDOUT_MONTHS = 6,  # 6 months
  K_AHEAD = c(1, 3, 6, 12),  # K-months ahead validation
  MIN_TRAIN_MONTHS = 36,  # Minimum training months before validation starts
  VALIDATION_STEP = 6,  # Step size for rolling windows (months)
  
  # Cross-validation
  CV_FOLDS = 5,
  
  # Computational parameters
  PARALLEL = TRUE,
  N_CORES = min(parallel::detectCores() - 1, 8),
  
  # Random seed for reproducibility
  SEED = 12345
)

# Set seed
set.seed(CONFIG$SEED)

# Save configuration
saveRDS(CONFIG, "outputs/data/config.rds")
cat("✓ Configuration saved with k-months validation framework\n")

# ==============================================================================
# 1. GET WV COUNTY SPATIAL DATA
# ==============================================================================

cat("\n--- Step 1: Obtaining WV County Spatial Data ---\n")

# Initialize boundary tracking variable
USING_REAL_BOUNDARIES <- FALSE

# Try to get real WV county boundaries
tryCatch({
  if (require(tigris, quietly = TRUE)) {
    cat("Downloading WV county boundaries from US Census...\n")
    wv_counties_sf <- tigris::counties(state = "WV", year = 2020)
    
    county_coords <- wv_counties_sf %>%
      mutate(
        county_name = tolower(NAME),
        centroid = st_centroid(geometry)
      ) %>%
      mutate(
        lon = st_coordinates(centroid)[,1],
        lat = st_coordinates(centroid)[,2]
      ) %>%
      st_drop_geometry() %>%
      select(county_name, lon, lat) %>%
      arrange(county_name)
    
    wv_nb <- spdep::poly2nb(wv_counties_sf, queen = TRUE)
    adj_matrix_contiguity <- spdep::nb2mat(wv_nb, style = "B")
    rownames(adj_matrix_contiguity) <- colnames(adj_matrix_contiguity) <- 
      tolower(wv_counties_sf$NAME)
    
    USING_REAL_BOUNDARIES <- TRUE
    cat("✓ Downloaded", nrow(county_coords), "WV counties from US Census\n")
    
  } else {
    stop("tigris package not available")
  }
  
}, error = function(e) {
  cat("Could not download real boundaries, using simplified test data...\n")
  
  # Use simplified county coordinates in case of failure to download the real WV county boundaries
  county_coords <- data.frame(
    county_name = c("barbour", "berkeley", "boone", "braxton", "brooke", "cabell", 
                    "calhoun", "clay", "doddridge", "fayette", "gilmer", "grant", 
                    "greenbrier", "hampshire", "hancock", "hardy", "harrison", "jackson",
                    "jefferson", "kanawha", "lewis", "lincoln", "logan", "marion",
                    "marshall", "mason", "mcdowell", "mercer", "mineral", "mingo",
                    "monongalia", "monroe", "morgan", "nicholas", "ohio", "pendleton",
                    "pleasants", "pocahontas", "preston", "putnam", "raleigh", "randolph",
                    "ritchie", "roane", "summers", "taylor", "tucker", "tyler", 
                    "upshur", "wayne", "webster", "wetzel", "wirt", "wood", "wyoming"),
    lon = c(-80.1, -77.9, -81.7, -80.6, -80.5, -82.3, -81.1, -81.1, -80.4, -81.1,
            -80.8, -79.1, -80.4, -78.8, -80.5, -78.9, -80.3, -81.8, -77.8, -81.6,
            -80.4, -82.1, -81.9, -80.1, -80.7, -82.0, -81.5, -81.1, -78.9, -82.3,
            -79.9, -80.5, -78.2, -80.9, -80.7, -79.3, -81.2, -80.1, -79.7, -81.9,
            -81.2, -79.9, -80.9, -81.4, -80.8, -80.2, -79.6, -80.9, -80.2, -82.4,
            -80.4, -80.7, -81.4, -81.4, -81.6),
    lat = c(39.1, 39.4, 37.9, 38.7, 40.3, 38.4, 38.8, 38.5, 39.4, 38.0, 38.9, 39.3,
            37.8, 39.3, 40.5, 39.0, 39.3, 38.9, 39.3, 38.4, 39.1, 38.2, 37.8, 39.5,
            39.8, 38.7, 37.4, 37.4, 39.4, 37.7, 39.6, 37.5, 39.6, 38.2, 40.1, 38.7,
            39.4, 38.1, 39.3, 38.4, 37.8, 38.9, 39.0, 38.5, 37.7, 39.3, 39.1, 39.2,
            38.9, 38.2, 38.4, 39.6, 39.1, 39.3, 37.6)
  )
  
  # Create spatial relationships using distance threshold
  dist_matrix <- as.matrix(dist(county_coords[, c("lon", "lat")]))
  rownames(dist_matrix) <- colnames(dist_matrix) <- county_coords$county_name
  
  distance_threshold <- 1.2
  nb_list <- spdep::dnearneigh(as.matrix(county_coords[, c("lon", "lat")]), 0, distance_threshold)
  adj_matrix_contiguity <- spdep::nb2mat(nb_list, style = "B")
  rownames(adj_matrix_contiguity) <- colnames(adj_matrix_contiguity) <- county_coords$county_name
  
  cat("Could not download real boundaries, using simplified test data...\n ✓ Using simplified spatial data for", nrow(county_coords), "counties\n")
})

# Distance matrix for all spatial calculations
dist_matrix <- as.matrix(dist(county_coords[, c("lon", "lat")]))
rownames(dist_matrix) <- colnames(dist_matrix) <- county_coords$county_name

# Save spatial data info
spatial_data_info <- list(
  using_real_boundaries = USING_REAL_BOUNDARIES,
  n_counties = nrow(county_coords),
  spatial_source = if (USING_REAL_BOUNDARIES) "US Census TIGER" else "Simplified coordinates"
)
saveRDS(spatial_data_info, "outputs/data/spatial_data_info.rds")

# ==============================================================================
# 2-4. LOAD, PROCESS DATA AND CREATE COMPLETE GRID
# ==============================================================================

cat("\n--- Steps 2-4: Loading and Processing Drug Death Data ---\n")

# Load datasets
opioid_data <- read_csv("monthly_opioid_county_death.csv", show_col_types = FALSE)
stimulant_data <- read_csv("monthly_stimulant_county_death.csv", show_col_types = FALSE)

cat("Loaded", nrow(opioid_data), "opioid observations\n")
cat("Loaded", nrow(stimulant_data), "stimulant observations\n")

# Create complete time sequence
all_months <- seq(from = as.Date(paste(CONFIG$START_YEAR, CONFIG$START_MONTH, "01", sep = "-")),
                  to = as.Date(paste(CONFIG$END_YEAR, CONFIG$END_MONTH, "01", sep = "-")),
                  by = "month")

year_month_seq <- format(all_months, "%Y-%m")
counties_in_data <- unique(c(opioid_data$Residence_County, stimulant_data$Residence_County))

# Create complete grid
complete_grid <- expand.grid(
  Residence_County = counties_in_data,
  year_month = year_month_seq,
  stringsAsFactors = FALSE
) %>%
  mutate(
    year = as.integer(substr(year_month, 1, 4)),
    month = as.integer(substr(year_month, 6, 7)),
    date = as.Date(paste(year_month, "01", sep = "-"))
  ) %>%
  arrange(Residence_County, date)

# Merge and fill zeros
merge_and_fill <- function(data, response_col) {
  complete_grid %>%
    left_join(data, by = c("Residence_County", "year", "year_month")) %>%
    mutate(
      distinct_patient_count = ifelse(is.na(distinct_patient_count), 0, distinct_patient_count),
      !!response_col := ifelse(is.na(.data[[response_col]]), 0, .data[[response_col]]),
      population = ifelse(is.na(population), NA, population)
    ) %>%
    group_by(Residence_County) %>%
    arrange(date) %>%
    mutate(population = zoo::na.fill(population, "extend")) %>%
    ungroup()
}

opioid_complete <- merge_and_fill(opioid_data, "opioid_death")
stimulant_complete <- merge_and_fill(stimulant_data, "stimulant_death")

cat("✓ Opioid complete data:", nrow(opioid_complete), "observations\n")
cat("✓ Stimulant complete data:", nrow(stimulant_complete), "observations\n")

# ==============================================================================
# 5. ADD SPATIAL AND TEMPORAL INDICES 
# ==============================================================================

county_index <- data.frame(
  county_name = sort(unique(complete_grid$Residence_County)),
  county_id = 1:length(unique(complete_grid$Residence_County))
)

time_index <- data.frame(
  year_month = sort(unique(complete_grid$year_month)),
  time_id = 1:length(unique(complete_grid$year_month))
)

add_indices <- function(data) {
  data %>%
    left_join(county_index, by = c("Residence_County" = "county_name")) %>%
    left_join(time_index, by = "year_month") %>%
    mutate(
      space_time_id = (county_id - 1) * nrow(time_index) + time_id,
      log_pop_offset = ifelse(population > 0, log(population), 0)
    ) %>%
    arrange(county_id, time_id)
}

opioid_complete <- add_indices(opioid_complete)
stimulant_complete <- add_indices(stimulant_complete)

# ==============================================================================
# 6. CREATE VALIDATION FRAMEWORK
# ==============================================================================

cat("\n--- Step 6: Creating Proper Validation Framework ---\n")

create_proper_splits <- function(data, CONFIG) {
  max_time_id <- max(time_index$time_id)
  
  # Use 6-month holdout to keep some 2023 data for learning assessment period
  # This allows the model to learn what assessment_period=1 means
  train_time_threshold <- max_time_id - 6
  
  final_train <- data %>% filter(time_id <= train_time_threshold)
  final_test <- data %>% filter(time_id > train_time_threshold)
  
  # Verify that assessment_period appears in training data
  if ("assessment_period" %in% names(final_train)) {
    assessment_in_train <- sum(final_train$assessment_period == 1, na.rm = TRUE)
    cat("  Assessment period observations in training:", assessment_in_train, "\n")
  }
  
  # K-MONTHS AHEAD VALIDATION: Rolling windows within training period  
  min_train_time <- CONFIG$MIN_TRAIN_MONTHS
  validation_splits <- list()
  
  # Create rolling windows for model selection
  window_starts <- seq(min_train_time, train_time_threshold - max(CONFIG$K_AHEAD), 
                      by = CONFIG$VALIDATION_STEP)
  
  cat("  Creating", length(window_starts), "validation windows for k-months ahead validation\n")
  
  for (i in seq_along(window_starts)) {
    train_end <- window_starts[i]
    
    # For each k-ahead value, create validation set
    for (k in CONFIG$K_AHEAD) {
      val_start <- train_end + 1
      val_end <- min(train_end + k, train_time_threshold)
      
      if (val_end > val_start) {  # Ensure we have validation data
        validation_splits[[paste0("window_", i, "_k", k)]] <- list(
          train = data %>% filter(time_id <= train_end),
          validation = data %>% filter(time_id >= val_start & time_id <= val_end),
          train_end_time = train_end,
          validation_start_time = val_start,
          validation_end_time = val_end,
          k_ahead = k,
          window_id = i
        )
      }
    }
  }
  
  # Return comprehensive validation structure
  list(
    # Final evaluation splits (6-month holdout)
    train = final_train,
    test = final_test,
    holdout_start_time = train_time_threshold + 1,
    n_train = nrow(final_train),
    n_test = nrow(final_test),
    
    # K-months ahead validation splits for model selection
    k_months_validation = validation_splits,
    n_validation_windows = length(validation_splits),
    
    # Metadata
    validation_config = list(
      k_ahead_values = CONFIG$K_AHEAD,
      min_train_months = CONFIG$MIN_TRAIN_MONTHS,
      validation_step = CONFIG$VALIDATION_STEP,
      n_windows = length(window_starts),
      holdout_months = 6  # 6 months instead of 12, see above
    )
  )
}

# Create proper splits for both datasets
opioid_splits <- create_proper_splits(opioid_complete, CONFIG)
stimulant_splits <- create_proper_splits(stimulant_complete, CONFIG)

cat("✓ Validation framework created:\n")
cat("  Final evaluation - Opioid train/test:", opioid_splits$n_train, "/", opioid_splits$n_test, "\n")
cat("  Final evaluation - Stimulant train/test:", stimulant_splits$n_train, "/", stimulant_splits$n_test, "\n")
cat("  K-months validation windows:", opioid_splits$n_validation_windows, "\n")
cat("  K-ahead values:", paste(CONFIG$K_AHEAD, collapse = ", "), "months\n")
cat("  Holdout period: 6 months (allows first 6 months of assessment period learning)\n")

# ==============================================================================
# 7. CREATE PRECISION MATRICES
# ==============================================================================

cat("\n--- Step 7: Creating Precision Matrices ---\n")

data_counties <- sort(unique(opioid_complete$Residence_County))
spatial_counties <- rownames(adj_matrix_contiguity)
common_counties <- intersect(data_counties, spatial_counties)

county_match_idx <- match(common_counties, spatial_counties)
adj_matrix_matched <- adj_matrix_contiguity[county_match_idx, county_match_idx]
dist_matrix_matched <- dist_matrix[county_match_idx, county_match_idx]

create_precision_matrices <- function(adj_mat, dist_mat) {
  n <- nrow(adj_mat)
  
  # Binary adjacency
  Q_adj <- diag(rowSums(adj_mat)) - adj_mat
  
  # Distance-weighted (exponential)
  W_exp <- exp(-dist_mat / median(dist_mat[dist_mat > 0]))
  diag(W_exp) <- 0
  W_exp[adj_mat == 0] <- 0
  Q_exp <- diag(rowSums(W_exp)) - W_exp
  
  # Distance-weighted (Gaussian)
  W_gauss <- exp(-(dist_mat^2) / (2 * (median(dist_mat[dist_mat > 0]))^2))
  diag(W_gauss) <- 0
  W_gauss[adj_mat == 0] <- 0
  Q_gauss <- diag(rowSums(W_gauss)) - W_gauss
  
  # Ensure positive definiteness
  epsilon <- 1e-6
  Q_adj <- Q_adj + epsilon * diag(n)
  Q_exp <- Q_exp + epsilon * diag(n)
  Q_gauss <- Q_gauss + epsilon * diag(n)
  
  list(
    adjacency = as(Q_adj, "Matrix"),
    exponential = as(Q_exp, "Matrix"),
    gaussian = as(Q_gauss, "Matrix"),
    county_names = common_counties
  )
}

precision_matrices <- create_precision_matrices(adj_matrix_matched, dist_matrix_matched)

# ==============================================================================
# 8. SAVE ALL OUTPUTS WITH VALIDATION
# ==============================================================================

cat("\n--- Step 8: Saving Enhanced Outputs ---\n")

# Save datasets
saveRDS(opioid_complete, "outputs/data/opioid_complete.rds")
saveRDS(stimulant_complete, "outputs/data/stimulant_complete.rds")
saveRDS(opioid_splits, "outputs/data/opioid_splits.rds")
saveRDS(stimulant_splits, "outputs/data/stimulant_splits.rds")

# Save spatial data
saveRDS(precision_matrices, "outputs/data/precision_matrices.rds")
saveRDS(county_index, "outputs/data/county_index.rds")
saveRDS(time_index, "outputs/data/time_index.rds")
saveRDS(county_coords, "outputs/data/county_coords.rds")

# Enhanced summary statistics
summary_stats <- list(
  opioid = list(
    n_counties = length(unique(opioid_complete$Residence_County)),
    n_months = length(unique(opioid_complete$year_month)),
    total_observations = nrow(opioid_complete),
    non_zero_observations = sum(opioid_complete$distinct_patient_count > 0),
    zero_proportion = sum(opioid_complete$distinct_patient_count == 0) / nrow(opioid_complete),
    mean_count = mean(opioid_complete$distinct_patient_count),
    max_count = max(opioid_complete$distinct_patient_count)
  ),
  stimulant = list(
    n_counties = length(unique(stimulant_complete$Residence_County)),
    n_months = length(unique(stimulant_complete$year_month)),
    total_observations = nrow(stimulant_complete),
    non_zero_observations = sum(stimulant_complete$distinct_patient_count > 0),
    zero_proportion = sum(stimulant_complete$distinct_patient_count == 0) / nrow(stimulant_complete),
    mean_count = mean(stimulant_complete$distinct_patient_count),
    max_count = max(stimulant_complete$distinct_patient_count)
  ),
  spatial = list(
    n_counties_with_spatial = length(precision_matrices$county_names),
    spatial_structure_types = names(precision_matrices)[1:3]
  ),
  # Validation framework summary
  validation = list(
    k_ahead_values = CONFIG$K_AHEAD,
    n_validation_windows = opioid_splits$n_validation_windows,
    validation_framework = "Enhanced k-months ahead + temporal holdout"
  )
)

saveRDS(summary_stats, "outputs/data/summary_stats.rds")

cat("\n=== ENHANCED DATA SETUP COMPLETE ===\n")
cat("VALIDATION FRAMEWORK:\n")
cat("✓ Temporal holdout: Last", CONFIG$HOLDOUT_MONTHS, "months for final evaluation\n")
cat("✓ K-months ahead:", length(opioid_splits$k_months_validation), "validation windows\n")
cat("✓ K-ahead values:", paste(CONFIG$K_AHEAD, collapse = ", "), "months\n")
cat("✓ Model family: Zero-Inflated Poisson\n")
cat("✓ Zero proportions - Opioid:", round(summary_stats$opioid$zero_proportion, 3), 
    "| Stimulant:", round(summary_stats$stimulant$zero_proportion, 3), "\n")

# Create enhanced status report
data_status <- data.frame(
  Component = c("Enhanced Configuration", "Opioid Data", "Stimulant Data", 
                "Final Train/Test Splits", "K-Months Validation", "Spatial Matrices", 
                "County Indices", "Time Indices", "Summary Stats"),
  Status = rep("✓ Complete", 9),
  File = c("config.rds", "opioid_complete.rds", "stimulant_complete.rds", 
           "splits contain final train/test", "splits contain k_months_validation", 
           "precision_matrices.rds", "county_index.rds", "time_index.rds", "summary_stats.rds"),
  Enhancement = c("K-months validation config", "", "", "", 
                 "Rolling windows", "", "", "", "Added validation summary"),
  stringsAsFactors = FALSE
)

write.csv(data_status, "outputs/data/enhanced_setup_status.csv", row.names = FALSE)

cat("\n✓ Enhanced validation framework ready for phases 1-4\n")
cat("✓ Phases use k-months validation for model selection\n")
cat("✓ Final evaluation uses 6-month temporal holdout\n")

# ==============================================================================
# END OF PHASE 0
# ==============================================================================