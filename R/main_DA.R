#' Interpretation of Scalar-on-Function Regression Models
#'
#' This function extracts locally constant point estimates from a given Bayesian scalar-on-function regression (SOFR) model.
#' It does not have to be specific to the BASOFR approach. The output includes plots and a list of results.
#'
#' @param beta.post Posterior draws of the beta parameter from the Bayesian SOFR model.
#' @param Alpha.post Posterior draws of the Alpha parameter from the Bayesian SOFR model. Can be NULL if there are no scalar covariates.
#' @param intercept.post Posterior draws of the intercept from the Bayesian SOFR model.
#' @param sigma2.post Posterior draws of the sigma-squared parameter from the Bayesian SOFR model.
#' @param y Vector of observations for the dependent variable.
#' @param X Matrix of functional covariate observations.
#' @param z Matrix of scalar covariate observations. Can be NULL.
#' @param t.grid A numeric vector of the discrete points on which the functions X are observed.
#' @param plot_all_locally_constant Logical. If TRUE, all locally constant estimates will be plotted.
#' @param K Number of partitions for the aggregated trajectories. Defaults to 20.
#' @param eps_level A numeric value for selecting an acceptable family. It determines the threshold for predictive performance of models under the Bayesian SOFR model. Defaults to 0.1.
#'
#' @return A list containing:
#' \itemize{
#'   \item {simplest_acc_family: A list of vectors representing the simplest member in the acceptable family for each functional regression coefficient.}
#'   \item {locally_constant: A list of matrices representing all of the locally constant estimates in the solution path of the fused lasso problem.}
#'   \item {empirical_MSE: A list of vectors containing the empirical MSE of all of the locally constant models in the fused lasso solution path.}
#'   \item {predictive_MSE: A list of matrices representing the predictive MSE of all of the locally constant models in the fused lasso solution path.}
#' }
#'
#' @examples
#' \dontrun{
#' # simulate the data for BASOFR
#' data = simulate_SOFR(n = 1500, p = 3, SNR = 5, num_z_continuous = 5, num_z_dummy = 2,
#' beta_types = c('spiky','smooth','stepwise'), domain_X = c(0, 1),
#'               max_obs = 105,
#'               subject_specific_domain = TRUE)
#'
#' # Scale the continuous z-variables to have std 1 and mean 0, and create the z-variables matrix
#' z = cbind(scale(data$Z_continuous), data$Z_dummy)
#'
#' # Fit the BASOFR model
#' BASOFR_fit <- fitBASOFR(y=data$y, X=data$X, t.grid = data$t.grid,
#'                         z_continuous = data$Z_continuous, z_dummy = data$Z_dummy,
#'                         S = 1000, burnin = 100)
#'
#' plot_BASOFR_fit(BASOFR_fit)
#'
#' # Decision Analysis
#' BayesSOFR <- BASOFR_fit # it does not have to be a BASOFR model
#'
#' # Extract more interpretable locally constant estimates from the Bayesian SOFR model:
#' SOFR_Interpret <- Interpret_SOFR(beta.post = beta.post, Alpha.post = Alpha.post, intercept.post = intercept.post, sigma2.post = sigma.2.post,
#'                y = data$y, X = data$X, z=z, t.grid=data$t.grid)
#'
#' names(SOFR_Interpret)
#' }
#' @export
Interpret_SOFR <- function(beta.post, Alpha.post=NULL, intercept.post, sigma2.post,
                           y, X, z=NULL, t.grid, plot_all_locally_constant = FALSE, K=20, eps_level = 0.1){
  n = length(y)

  # stored iterations in the MCMC
  S = dim(beta.post)[2]

  p = dim(beta.post)[1]

  # check if S is greater than 1000
  if(S > 1000) {
    cat('The posterior distribution has more than 1000 samples. We will reduce this to 1000 by evenly sampling across the distribution.\n')
    # create an evenly spaced sequence of 1000 indices spanning the range of S
    indices <- round(seq(from = 1, to = 1001, length.out = 1000))

    # use these indices to select the corresponding slices from each array
    beta.post <- beta.post[, indices, ]
    if(!is.null(Alpha.post)){
      Alpha.post <- Alpha.post[indices,]
    }
    intercept.post <- intercept.post[indices]
    sigma2.post <- as.matrix(sigma2.post[indices,])
  }

  cat("Start calculating the aggregated trajectories of the functional covariates...")
  # 1. obtain the aggregated trajectories
  x.subdomain.aver.K <- get_aggregated_trajectory(X = X, t.grid = t.grid, K=K)
  gc()


  cat("Start to solve the 'fit-to-the-fit' to obtain the locally constant estimates...")
  # 2. obtain the "fit-to-the-fit"

  # 2.1 obtain the target variable of the "fit-to-the-fit" using pseudo-data
  h_bar <- get_h_bar(X = X,
                     beta.post = beta.post,
                     t.grid = t.grid)


  # 2.2 obtain the locally constant estimates
  locallyConstant_fitToAFit <- get_locallyconstant(h_bar=h_bar, x.subdomain.aver.K=x.subdomain.aver.K,
                                                   beta.post=beta.post, Alpha.post=Alpha.post, intercept.post=intercept.post,sigma2.post = sigma2.post,
                                                   t.grid=t.grid, y=y, X=X, z=z,
                                                   plot_all_locally_constant=plot_all_locally_constant)

  cat("Start to extract the acceptable family... the simplest member in the acceptable family will be plotted automatically.")
  # 3. Extract the acceptable family
  acc_family <- model_selection(locallyConstant_fitToAFit = locallyConstant_fitToAFit,
                                eps_level=eps_level, t.grid = t.grid, beta.post = beta.post)

  DA_result = suppressWarnings(list())
  # a length p list, each component is the simplest member in the acceptable family.
  # (The function is estimated on the same discrete points that the functional covariates are observed on)
  DA_result$simplest_acc_family = lapply(acc_family, function(x) x$simplest_accFam)

  # a length p list of matrices, each matrix represents all of the locally constant estimates in the solution path of the fusedlasso problem
  # (The locally constant model is estimated on the same discrete points that the functional covariates are observed on)
  DA_result$locally_constant = locallyConstant_fitToAFit$LocallyConstant.est

  # a length p list of vectors, each vector contains the empirical MSE of all of the locally constant model in the fused lasso solution path
  DA_result$empirical_MSE = locallyConstant_fitToAFit$point.MSE

  # a length p list of matrices, each matrix (ncol=S, nrow = # of locally constant model in the fused lasso solution path) represents the predictive MSE of all of the locally constant model in the fused lasso solution
  DA_result$predictive_MSE = locallyConstant_fitToAFit$predictive.MSE

  return(DA_result)
}
