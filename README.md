MoCo: **Mo**tion-**Co**ntrolled Brain Phenotype Differences Between Groups <img src="fig/MoCo.png" width="120" align="right"/>
===================================================

MoCo is an R package designed to remove motion artifacts in brain phenotype analysis. Please note that MoCo is still under development.

## Installation

Ensure you have the following R packages installed:

-   SuperLearner
-   haldensify
-   devtools

You can install them by running the following code:

``` r
if(!require(c("SuperLearner","haldensify", "devtools"))){
    install.packages(c("SuperLearner","haldensify", "devtools"))
}
```

Then, you can install MoCo from GitHub using the following code:

```{r}
library(devtools)
install_github("thebrisklab/MoCo")

library(MoCo)
```

## Tutorial

The `moco()` function serves as the main function of MoCo. The input and output of the function are illustrated in the figure below.

<img src="fig/input_output.png" width="600" align="center"/>

### Explanation of Arguments  

```
moco(
  X, Z, A, M, Y, 
  Delta_M, 
  thresh = NULL,
  Delta_Y,
  SL_library = c("SL.earth","SL.glmnet","SL.gam","SL.glm", "SL.glm.interaction", "SL.step","SL.step.interaction","SL.xgboost","SL.ranger","SL.mean"),
  SL_library_customize = list(
    gA = NULL, 
    gDM = NULL,
    gDY_AX = NULL,
    gDY_AXZ = NULL,
    mu_AMXZ = NULL,
    eta_AXZ = NULL,
    eta_AXM = NULL,
    xi_AX = NULL
  ), 
  glm_formula = list(gA = NULL, 
                     gDM = NULL,
                     gDY_AX = NULL,
                     gDY_AXZ = NULL,
                     mu_AMXZ = NULL,
                     eta_AXZ = NULL,
                     eta_AXM = NULL,
                     xi_AX = NULL,
                     pMX = NULL,
                     pMXZ = NULL),
  HAL_pMX = TRUE,
  HAL_pMXZ = TRUE,
  HAL_options = list(max_degree = 3, lambda_seq = exp(seq(-1, -10, length = 100)), num_knots = c(1000, 500, 250)),
  cross_fit = TRUE,
  cv_folds = 5,
  test = TRUE,
  fwer = 0.05, 
  seed_rgn = 1, 
  ...
)
```

- `A`: A binary vector of length n (number of participants), serving as a group indicator, such as diagnosis group or control group.
- `X`: A dataframe or matrix containing demographic confounders that would ideally be balanced in a randomized controlled trial.
- `Z`: A dataframe or matrix of covariates representing brain phenotypes.
- `M`: A numeric vector of length n representing continuous motion values for each participant.
- `Y`: A matrix of dimension n $\times$ p, where n is the number of participants, and p is the number of regions of interest.
          If it represents seed-based association measures:
            Each (i, j) element denotes participant i's association measure between the seed region and region j.
            The column representing the association measure of the seed region with itself should be filled with NA values to indicate its position.
          If it represents other types of association measures:
            Each (i, j) element denotes participant i's association measure between two brain regions of interest, such as the upper diagonal part of the functional connectivity matrix.
            No NA values are allowed in \code{Y} in this case.
- Delta_M Binary vector of length n indicating data usability, based on the motion value.
                It corresponds to a binary variable indicating whether motion is available and meets inclusion criteria for conventional analyses.
                If motion meets inclusion criteria for analysis, set \code{Delta_M = 1}; otherwise, set \code{Delta_M = 0}.                
- thresh Value used to threshold M to produce Delta_M. One can specify either Delta_M or thresh.
- Delta_Y Binary vector indicating the non-missingness of brain image data \code{Y}. 
                Set \code{Delta_Y = 1} if \code{Y} is available; otherwise, set \code{Delta_Y = 0}.
 
- SL_library SuperLearner library for estimating nuisance regressions. 
                   Defaults to c("SL.earth","SL.glmnet","SL.gam","SL.glm", "SL.glm.interaction", "SL.step","SL.step.interaction","SL.xgboost","SL.ranger","SL.mean") if not specified.
- SL_library_customize Customize SuperLearner library for estimating each nuisance regression. 
                    - \code{gA}: SuperLearner library for estimating the propensity score.
                    - \code{gDM}: SuperLearner library for estimating the probability P(Delta_M = 1 | A, X).
                    - \code{gDY_AX}: SuperLearner library for estimating the probability P(Delta_Y = 1 | A, X).
                    - \code{gDY_AXZ}: SuperLearner library for estimating the probability P(Delta_Y = 1 | A, X, Z).
                    - \code{mu_AMXZ}: SuperLearner library for estimating the outcome regression E(Y | Delta_Y = 1, A, M, X, Z).
                    - \code{eta_AXZ}: SuperLearner library for estimating E(mu_AMXZ pMXD / pMXZD | A, X, Z, Delta_M = 1).                
                    - \code{eta_AXM}: SuperLearner library for estimating E(mu_AMXZ pMX/pMXZ gDY_AX/gDY_AXZ | A, M, X, Delta_Y = 1).
                    - \code{xi_AX}: SuperLearner library for estimating E(eta_AXZ | A, X).
                    
- glm_formula All glm formulas default to NULL, indicating SuperLearner will be used for nuisance regressions.
                    - \code{gA}: GLM formula for estimating the propensity score.
                    - \code{gDM}: GLM formula for estimating the probability P(Delta_M = 1 | A, X).
                    - \code{gDY_AX}: GLM formula for estimating the probability P(Delta_Y = 1 | A, X).
                    - \code{gDY_AXZ}: GLM formula for estimating the probability P(Delta_Y = 1 | A, X, Z).
                    - \code{mu_AMXZ}: GLM formula for estimating the outcome regression E(Y | Delta_Y = 1, A, M, X, Z).
                    - \code{eta_AXZ}: GLM formula for estimating E(mu_AMXZ pMXD / pMXZD | A, X, Z, Delta_M = 1).                
                    - \code{eta_AXM}: GLM formula for estimating E(mu_AMXZ pMX/pMXZ gDY_AX/gDY_AXZ | A, M, X, Delta_Y = 1).
                    - \code{xi_AX}: GLM formula for estimating E(eta_AXZ | A, X).
                    - \code{pMX}: GLM formula for estimating p(m | a, x, Delta_Y = 1) and p(m | a, x, Delta_M = 1), assuming M follows a log normal distribution.
                    - \code{pMXZ}: GLM formula for estimating p(m | a, x, z, Delta_Y = 1) and p(m | a, x, z, Delta_M = 1), assuming M follows a log normal distribution.
                    
- HAL_pMX Specifies whether to estimate p(m | a, x, Delta_Y = 1) and p(m | a, x, Delta_M=1) using the highly adaptive lasso conditional density estimation method. 
 Defaults to \code{TRUE}. If set to \code{FALSE}, please specify the \code{pMX} option in \code{glm_formula}, such as \code{pMX = "."}.
- HAL_pMXZ Specifies whether to estimate p(m | a, x, z, Delta_Y = 1) and p(m | a, x, z, Delta_M=1) using the highly adaptive lasso conditional density estimation method. 
 Defaults to \code{TRUE}. If set to \code{FALSE}, please specify the \code{pMXZ} option in \code{glm_formula}, such as \code{pMXZ = "."}.
  
 
- HAL_options Additional options for highly adaptive lasso (HAL) method.
                   - \code{max_degree}: The highest order of interaction terms for generating basis functions (passed to \code{haldensify}).
                   - \code{lambda_seq}: A numeric sequence of values for the regularization parameter of Lasso regression (passed to \code{haldensify}).
                   - \code{num_knots}: The maximum number of knot points (i.e., bins) for any covariate for generating basis functions (passed to \code{haldensify}).
 
- cross_fit Logical indicating whether to develop the estimator based on cross-fitting. Defaults to TRUE.
- test Logical indicating whether to conduct hypothesis testing based on simultaneous confidence band. Defaults to TRUE.
- fwer A vector of family-wise error rates (FWER) to control for multiple hypothesis testing. Defaults to 0.05. Set to NULL if hypo_test is FALSE.
- seed_rgn Specifies the value of seed(s) for nuisance regression calculation using super learner. Can be a vector. Defaults to value 1. 


