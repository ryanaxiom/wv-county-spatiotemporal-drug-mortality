# wv-county-spatiotemporal-drug-mortality
County-level and state-level spatiotemporal analysis of opioid and stimulant mortality in West Virginia (2015-2023) using Zero-Inflated Poisson models in R-INLA

## Overview
The goal of this project was to develop a single model using the R programming language to predict both county-level and state-level death counts due to opioid and stimulant overdose in WV. This required careful consideration of model type, development, selection, and error propagation at various levels. Overall, the project was exceptionally effective in both modeling the training set and predicting the testing set with excellent Bayesian credible interval coverage (99%) for both the individual counties as well as at the state level, while properly handling the 58.7% zero-inflation in the mortality data.

INLA provides a fast and efficient alternative to Markov Chain Monte Carlo (MCMC) methods for Bayesian inference in latent Gaussian models. INLA utilizes Integrated Nested Laplace Approximation to deterministically approximate the marginal posterior distributions. INLA is widely applicable to latent Gaussian models and only requires that the latent field (random effects, spatial effects, etc) must have a joint multivariate Gaussian distribution, though the observations themselves can follow any distribution in the exponential family (normal, Poisson, binomial, negative binomial, gamma, beta, etc). In that case, INLA is much faster (100-1000x+) and more computationally efficient than MCMC as well as being easier to implement (no need for convergence diagnostics or tuning) and generally very accurate. However, it is not a "silver bullet" solution as it cannot model non-Gaussian latent fields (eg models with t-distributed random effects), complex non-linear relationships that can't be approximated within LGM, discrete latent variables, random forests, neural networks. In contrast, MCMC works with non-Gaussian latent fields as well and on very small datasets (N<100) as approximate methods will be less accurate. It can also model complex non-linear relationships and will give the full joint posterior distributions (unlike INLA which approximates the posterior marginal distributions). This comes at the expense of speed (hours or days for complex spatial-temporal models), convergence issues, diagnostics burdens, high memory usage, and tuning random sampling methods. Choosing INLA over MCMC likely saved months of computational time in developing this predictive model while giving results that are essentially identical.

The data for this study is monthly death counts for each of the 55 counties from January 2015 through December 2023.

This exploratory data analysis (EDA) involved fitting 70+ zero-inflated Poisson (ZIP) models using the R-INLA procedure through a progression of increasing complexity where the direction of development was decided by the data. This would be impossible using MCMC since it would take weeks instead of about 30 minutes or less. The progression of models led to being able to determine the added benefits of sequential fits and a selection was made both by model fit statistics as well as favoring parsimony.

**Manuscript in preparation**

## Key Features

This code performs spatio-temporal modeling of county-level drug mortality by following an EDA approach where increasingly complex models were fitted depending on the previous output of earlier phases. This provides a way to perform separate analyses for opioid or stimulant mortality. However, though automatic selection is built in to this code, it is not a replacement for experience as the operator must pay close attention to when a more parsimonious model would be a better choice and modify the code to go in that direction.

The model complexities start with single component spatial and temporal models, then proceeds to separable spatio-temporal models, then onto several potentially promising interaction models. The separable spatio-temporal models assume that the spatial and temporal effects are independent of each other while the interaction models allow for the temporal effects to change based on the spatial levels.

This methodology handles zero-inflated count data very effectively. ZINB could possibly be better but there were convergence issues that could not be resolved despite many attempts, both simple and complex, to improve convergence.

### Methodological Innovations
- **Corrected credible intervals for county-level predictions**: R-INLA defaults to using the full pooled sample size for uncertainty quantification, which underestimates uncertainty at the county level. This analysis implements proper county-specific sample size adjustments for accurate 95% credible intervals.
- **Comprehensive model selection framework**: Beyond automated WAIC selection, the pipeline generates diagnostic visualizations across all phases to support informed decisions about model parsimony and practical deployment.
- **Parsimony-aware selection**: The final model (adjacency-based spatial effects) was chosen not solely on WAIC but also considering that it performed equivalently to more complex Gaussian and exponential spatial weightings while being simpler to implement and interpret.

### Diagnostic Outputs by Phase
- **Phase 1**: Spatial effect comparisons, maps of county-specific random effects
- **Phase 2/2a**: Temporal trend visualizations, county-specific vs state-level patterns  
- **Phase 3**: Separable model performance heatmaps, component contribution analysis
- **Phase 4**: Interaction strength assessment, period-specific spatial patterns
- **Phase 5**: County-level and state-level predictions with properly calibrated intervals, model coverage validation
- **Phase 6**: Executive dashboard synthesizing all model comparisons

## Key Findings

### Model Performance Progression
The analysis demonstrates clear improvements at each modeling stage:
- **Baseline → Single Component**: 28% MAE reduction
- **Single Component → Separable**: 20.6% MAE reduction  
- **Separable → Interaction**: 1.7% MAE reduction (52.8% of counties improved)
- **Overall (Baseline → Best)**: 44% MAE reduction

**Final model comparison:**
- **Best interaction model:** adjacency_ar1_period (WAIC: 10935.3)  
- **Best separable model:** adjacency_ar1 (WAIC: 11014.4)
- **Statistical evidence:** 79.1 WAIC units (decisive for interactions)

**Technical achievements:**
- ZIP models effectively handle 58.7% zero observations in the data
- Best model achieves 99% coverage with properly calibrated credible intervals

While the final interaction improvement is modest (1.7%), the 79.1 WAIC unit difference provides decisive statistical evidence for space-time interactions.

### Epidemiological Interpretation
The interaction effects indicate that:
- Spatial patterns of opioid deaths evolved somewhat differently across time periods
- Policy interventions (2020-2022) had county-specific impacts
- The epidemic's spatial and time structure changed, not just its intensity
- County-targeted interventions may be more effective than uniform state-level approaches

### Computational Trade-off

When comparing separable versus interaction models:
- **Runtime increase:** 1.5× (2.61s → 3.79s)
- **Complexity increase:** 167% (score 3 → 8)
- **Accuracy gain:** 1.7% over separable, 22.3% over single component
- **Recommendation:** The modest accuracy gain and strong statistical evidence justify the computational cost for research/planning purposes as well as real-time applications since the penalty is small (only 50% longer time to fit).

## Requirements

### Software
- R (>= 4.0.0)
- R-INLA package from www.r-inla.org

### R Packages
- INLA
- dplyr (>= 1.0.0)
- ggplot2 (>= 3.3.0)
- ggrepel
- cowplot
- sf
- sn
- spdep
- tidyr
- maps
- readr
- zoo
- DT
- tigris
- viridis
- gridExtra
- forecast
- scales
- knitr
- lubridate
- Matrix
- parallel

## Installation
```r
# Install R-INLA
install.packages("INLA", repos = c(getOption("repos"), 
                INLA = "https://inla.r-inla-download.org/R/stable"))

# Install other dependencies
install.packages(c("dplyr", "ggplot2", "ggrepel", "cowplot", "sf", "sn", "spdep", "tidyr", "maps", "readr", "zoo", "DT", "tigris", "viridis", "gridExtra", "forecast", "scales", "knitr", "lubridate", "Matrix", "parallel"))
```
## Data Requirements
Note: Original data cannot be shared due to confidentiality agreements

### Input Data Structure

Required csv files:

1. `monthly_opioid_county_death.csv` - Primary analysis data for one substance

2. `monthly_stimulant_county_death.csv` - Required for pipeline execution

**If you only have one substance data**, create a dummy stimulant file:
```r
   # Create dummy stimulant CSV that will survive train/test splits
months <- seq(as.Date("2015-01-01"), as.Date("2023-12-01"), by = "month")
dummy_stimulant <- data.frame(
  year = year(months),
  year_month = format(months, "%Y-%m"),
  Residence_County = rep(c("BERKELEY", "JEFFERSON", "MONONGALIA"), length.out = length(months)),
  distinct_patient_count = 0,
  population = 10000,
  stimulant_death = 0
)
write.csv(dummy_stimulant, "monthly_stimulant_county_death.csv", row.names = FALSE)
```

### Column Requirements

Both csv files must contain these exact columns:
   - `year` (integer): 2015-2023
   - `year_month` (string): "YYYY-MM" format
   - `Residence_County` (string): County name
   - `distinct_patient_count` (integer): Monthly deaths
   - `population` (integer): County population
   - `[drug]_death` (float): Rate per population


### Data Access

Researchers interested in similar data should contact the West Virginia Department of Health and Human Resources.

## Usage

### Expected Runtime
  Total pipeline runtime varies by hardware:
   - Apple M1 Max MacBook Pro: ~40-45 minutes
   - Modern desktop (Intel i9-14900KF): ~17 minutes
   
  Individual phase times on M1 Max:
   - Phase 0-1: ~2 minutes
   - Phase 2-2a: ~5 minutes  
   - Phase 3: ~8 minutes
   - Phase 4: ~25 minutes
   - Phase 5-6: ~5 minutes

  Individual phase times on M1 Max:
   - Phase 0-1: ~22 seconds
   - Phase 2-2a: ~40 seconds
   - Phase 3: ~40 seconds
   - Phase 4: ~14 minutes
   - Phase 5-6: ~1.5 minutes

### Running the Full Pipeline

By default, this analysis runs for opioid mortality. To analyze other outcomes:
```bash
# Default (opioid analysis)
Rscript 00_setup_and_data_zip.R

# See the csv requirements in the Input Data Structure section
export ANALYSIS_RESPONSE_TYPE="stimulant"
Rscript 00_setup_and_data_zip.R
```

```bash
# Set response type (opioid or stimulant)
export ANALYSIS_RESPONSE_TYPE="opioid"

# Run phases sequentially
Rscript 00_setup_and_data_zip.R
Rscript 01_spatial_models_zip.R
Rscript 02_temporal_models_zip.R
Rscript 02a_county_temporal_models_zip.R
Rscript 03_separable_spatiotemporal_zip.R
Rscript 04_full_interaction_models_zip.R
Rscript 05_model_diagnostics_zip.R
Rscript 06_visualization_zip.R
```

### Project Structure
```
├── 00_model_utilities.R              # Shared functions and utilities
├── 00_setup_and_data_zip.R           # Data preparation and splits
├── 01_spatial_models_zip.R           # Spatial-only models
├── 02_temporal_models_zip.R          # State-level temporal exploration
├── 02a_county_temporal_models_zip.R  # County-specific temporal
├── 03_separable_spatiotemporal_zip.R # Separable space-time models
├── 04_full_interaction_models_zip.R  # Space-time interaction models
├── 05_model_diagnostics_zip.R        # Diagnostics and county predictions
└── 06_visualization_zip.R            # Executive dashboard generation
```

## Important Note on Model Selection

While this pipeline automatically selects the model with the lowest WAIC, **this may not always be the optimal choice for your application**. In this analysis, the best WAIC model (adjacency_ar1_period) happened to also be reasonably parsimonious, but this alignment is not guaranteed.

### Considerations Beyond WAIC:
- **Computational resources**: A model with marginally better WAIC but 73× runtime may not be practical for operational use
- **Interpretability**: Simpler models may be preferred for policy communication even with slightly higher WAIC
- **Stability**: Complex models may overfit despite good WAIC scores
- **Deployment needs**: Real-time prediction systems may require faster, simpler models

### Recommended Workflow:
1. Run the full pipeline to generate all model comparisons
2. Review models within 10 WAIC units of the best
3. Examine the diagnostic plots from each phase
4. Consider your specific context (computational budget, interpretability needs, deployment constraints)
5. Make an informed selection that may differ from the automated choice

**Example**: If the exponential distance-weighted spatial model had WAIC of 10934 vs adjacency's 10935, the adjacency model would likely be preferred for its simplicity and interpretability despite marginally worse fit.

## Expected Outputs

### Critical Technical Outputs
- **County-level prediction plots**: 55 individual plots with properly calibrated 95% credible intervals (corrected for county-specific sample sizes)
- **Model selection diagnostics**: Comprehensive visualizations supporting parsimony decisions when multiple models have similar WAIC
- **Spatial effect comparisons**: Visual evidence that adjacency, Gaussian, and exponential weightings perform equivalently, justifying the simpler adjacency approach

### Directories Created

`outputs/data/` - Processed data objects
`outputs/diagnostics` - comprehensive analysis report, executive summary, and 
`outputs/models/` - Fitted model objects (.rds files)
`outputs/plots/` - Visualizations (.png files)
`outputs/models/` - Performance and phase summaries (.csv and .txt files)

### Key Output Files

`phase6_executive_dashboard.png` - Main results summary
`comprehensive_analysis_report_[opioid/stimulant].txt` - Detailed report
`phase6_all_models_comparison_[opioid/stimulant].csv` - Model metrics
State-level and county-level prediction plots with respective credible intervals

## Citation
Hansen, Ryan. (2025). West Virginia County-Level Spatiotemporal Drug 
Mortality Analysis. GitHub repository. 
https://github.com/ryanaxiom/wv-county-spatiotemporal-drug-mortality
DOI: [pending]

## License
MIT License - see LICENSE file for details

## Contact
Ryan Hansen - ryan.hansen@mail.wvu.edu
ORCID: https://orcid.org/0000-0003-1350-4921
West Virginia University
Project Link: https://github.com/ryanaxiom/wv-county-spatiotemporal-drug-mortality

## Acknowledgments

West Virginia Department of Health and Human Resources for data access
Mohammad Al-Mamun (mohammad.almamun@hsc.wvu.edu) for data cleaning and provision
Thanks to West Virginia Clinical and Translational Science Institute for funding this project
R-INLA development team

## Disclaimer
This analysis is for research purposes only. Results should not be used for clinical decision-making without consultation with appropriate health professionals. Individual-level predictions or county-specific results require appropriate data use agreements.
