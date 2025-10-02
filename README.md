# wv-county-spatiotemporal-drug-mortality
County-level and state-level spatiotemporal analysis of opioid and stimulant mortality in West Virginia (2015-2023) using Zero-Inflated Poisson models in R-INLA
## Overview
The goal of this project was to develop a single model using the R programming language to predict both county-level and state-level death counts due to opioid and stimulant overdose in WV.  This required careful consideration of model type, development, selection, and error propagation at various levels. Overall, the project was exceptionally effective in both modeling the training set and predicting the testing set with excellent Bayesian credible interval coverage for both the individual counties as well as at the state level.

INLA provides a fast and efficient alternative to Markov Chain Monte Carlo (MCMC) methods for Bayesian inference in latent Gaussian models. INLA utilizes Integrated Nested Laplace Approximation to deterministically approximate the marginal posterior distributions.  INLA is widely applicable to latent Gaussian models and only requires that the latent field (random effects, spatial effects, etc) must have a joint multivariate Gaussian distribution, though the observations themselves can follow any distribution in the exponential family (normal, Poisson, binomial, negative binomial, gamma, beta, etc).  In that case, INLA is much faster (100-1000x+) and more computationally efficient than MCMC as well as being easier to implement (no need for convergence diagnostics or tuning) and generally very accurate.  However, it is not a "silver bullet" solution as it cannot model non-Gaussian latent fields (eg models with t-distributed random effects), complex non-linear relationships that can't be approximated within LGM, discrete latent variables, random forests, nerual networks.  In contrast, MCMC works with non-Gaussian latent fields as well and on very small datasets (N<100) as approximate methods will be less accurate.  It can also model complex non-linear relationships and will give the full joint posterior distributions (unlike ILNA which approximates the posterior marginal distributions). This comes at the expense of speed (hours or days for complex spatial-temporal models), convergence issues, diagnostics burdens, high memory usage, and tuning random sampling methods. Choosing INLA over MCMC likely saved months of computational time in developing this predictive model while giving results that are essentially identical.

The data for this study is monthly death counts for each of the 55 counties from January 2015 through December 2023.

This exploratory data analysis (EDA) involved fitting 70+ zero-inflated Poisson (ZIP) models using the R-INLA procedure through a progression of increasing complexity where the direction of development was decided by the data. This would be impossible using MCMC since it would take weeks instead of about 30 minutes or less. The progression of models led to being able to determine the added benefits of sequential fits and a selection was made both by model fit statistics as well as favoring parsimony.

**Manuscript in preparation**

## Key Features

This code performs spatio-temporal modeling of county-level drug mortality by following an EDA approach where increasingly complex models were fitted depending on the previous output of earlier phases.  This provides a way to perform separate analyses for opioid or stimulant mortality.  However, though automatic selection is built in to this code, it is not a replacement for experience as the operator must pay close attention to when a more parsimonious model would be a better choice and modify the code to go in that direction.

The model complexities start with single component spatial and temporal models, then proceeds to separable spatio-temporal models, then onto several potentially promising interaction models.  The separable spatio-temporal models assume that the spatial and temporal effects are independent of each other while the interaction models allow for the temporal effects to change based on the spatial levels.

This methodology handles zero-inflated count data very effectively.

## Key Findings

### Model Comparison: Interaction vs Separable
The analysis provides **strong statistical evidence** for space-time interactions in opioid mortality patterns:

- **Best interaction model:** adjacency_ar1_period (WAIC: 10935.3)  
- **Best separable model:** adjacency_ar1 (WAIC: 11014.4)
- **WAIC improvement:** 79.1 units (decisive evidence for interactions)
- **MAE improvement:** 0.63 → 0.58 deaths/county-month (8% reduction)
- **RMSE improvement:** 1.01 → 0.94 deaths/county-month (7% reduction)

### Epidemiological Interpretation
The strong interaction effects indicate that:
- Spatial patterns of opioid deaths evolved differently across time periods
- Policy interventions (2020-2022) had county-specific impacts
- The epidemic's spatial structure changed, not just its intensity
- County-targeted interventions may be more effective than uniform state-level approaches

### Computational Trade-off
- Interaction models: 13 minutes runtime, complexity score 8
- Separable models: 2 minutes runtime, complexity score 3
- Overall recommendation: The 8% accuracy improvement justifies the additional computational cost

### Software

This code utilizes R (>= 4.0.0) as well as the INLA package to implement the R-INLA fitting procedures. The R-INLA package can be downloaded from www.r-inla.org . Following this, the package must be installed manually by adding the R-INLA repository to R's options and then using the install.packages() command.

### R Packages

## Installation

# Install R-INLA

# Install other dependencies
