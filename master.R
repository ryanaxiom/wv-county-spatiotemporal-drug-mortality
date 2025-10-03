#Sys.setenv(ANALYSIS_RESPONSE_TYPE = "stimulant") # Default: opioid, this option changes the analysis to stimulant
#Sys.setenv(ANALYSIS_PARALLEL = "FALSE") # Default: True, this option changes the parallelization
source("00_setup_and_data_zip.R")
source("01_spatial_models_zip.R")
source("02_temporal_models_zip.R")
source("02a_county_temporal_models_zip.R")
source("03_separable_spatiotemporal_zip.R")
source("04_full_interaction_models_zip.R")
source("05_model_diagnostics_zip.R")
source("06_visualization_zip.R")