# wv-county-spatiotemporal-drug-mortality
County-level and state-level spatiotemporal analysis of opioid and stimulant mortality in West Virginia (2015-2023) using Zero-Inflated Poisson models in R-INLA
## Overview
The goal of this project was to develop a single model using the R programming language to predict both county-level and state-level death counts due to opioid and stimulant overdose in WV.  This required careful consideration of model type, development, selection, and error propagation at various levels. Overall, the project was exceptionally effective in both modeling the training set and predicting the testing set with excellent Bayesian credible interval coverage for both the individual counties as well as at the state level.

INLA provides a fast and efficient alternative to Markov Chain Monte Carlo (MCMC) methods for Bayesian inference in latent Gaussian models. INLA utilizes Integrated Nested Laplace Approximation to deterministically approximate the posterior distributions.  INLA is widely applicable to a variety of models (GLM, GAM, spatio-temporal, survival models, Poisson, negative binomial, ZIP, ZINB, mixed effects, spline smoothing, heirarchical Bayesian, ...). INLA is much faster and more computationally efficient than MCMC and is easier to implement as well as being generally very accurate  In contrast, MCMC 

The data for this study is monthly death counts for each of the 55 counties from January 2015 through December 2023.

This exploratory data analysis (EDA) involved fitting a variety of zero-inflated Poisson (ZIP) models using the R-INLA procedure through a progression of increasing complexity where the direction of development was decided by the data. The progression of models led to being able to determine the added benefits of sequential fits and a selection was made both by model fit statistics as well as favoring parsimony.

**Manuscript in preparation**

## Key Features

### Software

This utilizes R (>= 4.0.0) as well as the INLA package to implement the R-INLA fitting procedures. The R-INLA package can be downloaded from www.r-inla.org . Following this, the package must be installed manually by adding the R-INLA repository to R's options and then using the install.packages() command.
