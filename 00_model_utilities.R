# Shared utilities for model fitting
fit_model_with_zinb_support <- function(formula, data, family, verbose = FALSE) {
  if (family == "zeroinflatednbinomial1") {
    cat("  [ZINB] Fitting NB first for initialization...\n")
    
    # Try NB first for starting values
    nb_model <- inla(
      formula = formula,
      family = "nbinomial",
      data = data,
      control.compute = list(dic = TRUE),
      verbose = FALSE
    )
    
    overdispersion_init <- nb_model$summary.hyperpar["size for nbinomial", "mean"]
    cat("  [ZINB] Overdispersion from NB:", round(overdispersion_init, 3), "\n")
    cat("  [ZINB] Now fitting Zero-Inflated Negative Binomial model...\n")
    
    hyper_zinb <- list(
      size = list(
        initial = log(max(overdispersion_init, 0.1)),
        prior = "loggamma",
        param = c(1, 0.01)
      )
    )
    
    result <- inla(
      formula = formula,
      family = family,
      data = data,
      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
      control.family = list(hyper = hyper_zinb),
      control.inla = list(
        strategy = 'adaptive',
        int.strategy = 'eb',
        tolerance = 1e-5,
        stupid.search = TRUE
      ),
      verbose = verbose
    )
    cat("  [ZINB] Model fitted successfully\n")
    
  } else {
    
    # ZIP or regular Poisson
    result <- inla(
      formula = formula,
      family = family,
      data = data,
      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
      control.inla = list(strategy = 'adaptive', int.strategy = 'auto'),
      verbose = verbose
    )

  }
  
  return(result)
}

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