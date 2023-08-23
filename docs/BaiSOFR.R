## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
# devtools::install_github("YunanGao/BaiSOFR")
library(BaiSOFR)

## ----cache=TRUE---------------------------------------------------------------
set.seed(1234) # set random seed to ensure duplication of the results

# simulate the SOFR data
data = simulate_SOFR(n = 1000, # number of observations
                     p = 3, # number of functional covariates
                     SNR = 5, # Signal-to-noise ratio
                     num_z_continuous = 5, # number of continuous covariates in z
                     num_z_dummy = 2, # number of dummy covariates in z
                     beta_types = c('spiky','smooth','stepwise'), # the ground truth shape of beta
                     domain_X = c(0, 1), # maximal domain over which X are being observed
                     max_obs = 200, # maximal number of discrete points to observe X
                     subject_specific_domain = TRUE # whether to generate subject-specific domain
                     )

# View the simulated dataset
names(data)

## ---- fig.width=6, fig.height=4-----------------------------------------------

# Visualize the ground truth beta functions
plot_beta(beta_list = data$groundTruth$beta, x.grids = data$t.grid)

## ----cache=TRUE---------------------------------------------------------------
# Create the z-variables matrix
z = cbind(scale(data$Z_continuous), data$Z_dummy)

# Fit the BASOFR model
BASOFR_fit <- fitBASOFR(y = data$y, X = data$X, 
                        t.grid = data$t.grid,  # seq(0, 1, length.out=200)
                        z = z,
                        burnin = 1000,
                        S = 1000)


## -----------------------------------------------------------------------------
names(BASOFR_fit)

## ----fig.width=6, fig.height=4------------------------------------------------

# Obtain the ground truth beta
groundTruth.beta = matrix(unlist(data$groundTruth$beta), byrow=TRUE, nrow=3)

# Visualize the fit of Bayesian SOFR
plot_BASOFR_fit(BASOFR_fit = BASOFR_fit,
                groundTruth = groundTruth.beta # optional
                )

## ----fig.width=5, fig.height=3.75, cache=TRUE---------------------------------
BayesSOFR <- BASOFR_fit # It can be any Bayesian SOFR model, does not limit to our BASOFR model.

SOFR_Interpret <- Interpret_SOFR(# the posterior samples
                                 beta.post = BayesSOFR$beta.post, 
                                 Alpha.post = BayesSOFR$Alpha.post,
                                 intercept.post = BayesSOFR$intercept.post,
                                 sigma2.post = BayesSOFR$sigma.2.post,
                                 # the data observations
                                 y = data$y, 
                                 X = data$X, 
                                 z = z,
                                 t.grid = data$t.grid)

## -----------------------------------------------------------------------------
names(SOFR_Interpret)

