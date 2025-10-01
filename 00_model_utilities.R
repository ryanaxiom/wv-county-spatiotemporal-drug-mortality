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