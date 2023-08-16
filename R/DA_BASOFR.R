#' Interpretation of Bayesian Additive Scalar-on-Function Regression Models
#'
#' This function is designed to simplify the user's interaction when dealing with the BASOFR approach.
#' It extracts locally constant point estimates from a given Bayesian scalar-on-function regression (SOFR) model.
#' The output includes plots and a list of results.
#'
#' @param BASOFR_fit A list obtained from the `fitBASOFR` function. It should contain the following components:
#'                   - `beta.post`: Posterior draws of the beta parameter from the Bayesian SOFR model.
#'                   - `Alpha_continuous_post`: Posterior draws of the Alpha parameter for continuous covariates.
#'                                                Can be NULL if there are no continuous scalar covariates.
#'                   - `Alpha_dummy_post`: Posterior draws of the Alpha parameter for dummy covariates.
#'                                         Can be NULL if there are no dummy scalar covariates.
#'                   - `intercept_post`: Posterior draws of the intercept from the Bayesian SOFR model.
#'                   - `sigma.2.post`: Posterior draws of the sigma-squared parameter from the Bayesian SOFR model.
#' @param y Vector of observations for the dependent variable.
#' @param X Matrix of functional covariate observations.
#' @param z Matrix of scalar covariate observations. Can be NULL.
#' @param t.grid A numeric vector of the discrete points on which the functions X are observed.
#' @param plot_all_locally_constant Logical. If TRUE, all locally constant estimates will be plotted. Defaults to TRUE.
#' @param K Number of partitions for the aggregated trajectories. Defaults to 20.
#' @param eps_level A numeric value for selecting an acceptable family. It determines the threshold for predictive performance
#'                  of models under the Bayesian SOFR model. Defaults to 0.1.
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
#' data = simulate_BASOFR(n = 1000, p = 3, SNR = 5, num_z_continuous = 0, num_z_dummy = 1,
#'                        x.grids = seq(0, 1, length.out = 100), beta_types = c('spiky','smooth', 'stepwise'))
#'
#' # Fit the BASOFR model
#' BASOFR_fit <- fitBASOFR(y=data$y, X=data$X, t.grid = data$t.grid,
#'                         z_continuous = data$Z_continuous, z_dummy = data$Z_dummy,
#'                         S = 1000, burnin = 100)
#'
#' # Plot the BASOFR fit
#' plot_BASOFR_fit(BASOFR_fit)
#'
#' # Check for presence of scalar covariates and create the appropriate `z` matrix
#' if (!is.null(data$Z_continuous) && !is.null(data$Z_dummy)) {
#'   z <- cbind(data$Z_continuous, data$Z_dummy)
#' } else if (is.null(data$Z_continuous) && !is.null(data$Z_dummy)) {
#'   z <- data$Z_dummy
#' } else if (!is.null(data$Z_continuous) && is.null(data$Z_dummy)) {
#'   z <- data$Z_continuous
#' } else if (is.null(data$Z_continuous) && is.null(data$Z_dummy)) {
#'   z <- NULL
#' }
#'
#' # Use Interpret_BASOFR function
#' Interpret_BASOFR_result <- Interpret_BASOFR(BASOFR_fit = BASOFR_fit,
#'                                             y=data$y, X=data$X, z=z,
#'                                             t.grid=data$t.grid)
#' }
#' @export
Interpret_BASOFR <- function(BASOFR_fit, y, X, z=NULL, t.grid, plot_all_locally_constant=TRUE, K=20, eps_level=0.1){
  cat("Start calculating the aggregated trajectories of the functional covariates...")
  # 1. obtain the aggregated trajectories
  x.subdomain.aver.K <- get_aggregated_trajectory(X = X, t.grid = t.grid, K=K)

  cat("Start to solve the 'fit-to-the-fit' to obtain the locally constant estimates...")
  # 2. obtain the "fit-to-the-fit"

  # 2.1 obtain the target variable of the "fit-to-the-fit" using pseudo-data
  h_bar <- get_h_bar(X = X, beta.post = BASOFR_fit$beta.post, t.grid = t.grid)

  # 2.2 obtain the locally constant estimates
  if (!is.null(BASOFR_fit$Alpha_continuous_post) && !is.null(BASOFR_fit$Alpha_dummy_post)) {
    Alpha.post <- cbind(BASOFR_fit$Alpha_continuous_post, BASOFR_fit$Alpha_dummy_post)
  } else if (is.null(BASOFR_fit$Alpha_continuous_post) && !is.null(BASOFR_fit$Alpha_dummy_post)) {
    Alpha.post <- BASOFR_fit$Alpha_dummy_post
  } else if (!is.null(BASOFR_fit$Alpha_continuous_post) && is.null(BASOFR_fit$Alpha_dummy_post)) {
    Alpha.post <- BASOFR_fit$Alpha_continuous_post
  } else if (is.null(BASOFR_fit$Alpha_continuous_post) && is.null(BASOFR_fit$Alpha_dummy_post)) {
    Alpha.post <- NULL
  }
  locallyConstant_fitToAFit <- get_locallyconstant(h_bar=h_bar, x.subdomain.aver.K=x.subdomain.aver.K,
                                                               beta.post=BASOFR_fit$beta.post, Alpha.post=Alpha.post,
                                                               intercept.post=BASOFR_fit$intercept_post, sigma2.post = BASOFR_fit$sigma.2.post,
                                                               t.grid=t.grid, y=y, X=X, z=z,
                                                               plot_all_locally_constant=plot_all_locally_constant)

  cat("Start to extract the acceptable family... \nThe model selection process, and the simplest member in the acceptable family will be plotted automatically.")
  # 3. Extract the acceptable family
  acc_family <- model_selection(locallyConstant_fitToAFit = locallyConstant_fitToAFit,
                                            eps_level=eps_level, t.grid = t.grid, beta.post = BASOFR_fit$beta.post)

  DA_result = suppressWarnings(list())
  # a length p list, each component is the simplest member in the acceptable family.
  # (The function is estimated on the same discrete points as the that the functional covariates are observed on)
  DA_result$simplest_acc_family = lapply(acc_family, function(x) x$simplest_accFam)

  # a length p list of matrices, each matrix represents all of the locally constant estimates in the solution path of the fusedlasso problem
  # (The locally constant model is estimated on the same discrete points as the that the functional covariates are observed on)
  DA_result$locally_constant = locallyConstant_fitToAFit$LocallyConstant.est

  # a length p list of vectors, each vector contains the empirical MSE of all of the locally constant model in the fused lasso solution path
  DA_result$empirical_MSE = locallyConstant_fitToAFit$point.MSE

  # a length p list of matrices, each matrix (ncol=S, nrow = # of locally constant model in the fused lasso solution path) represents the predictive MSE of all of the locally constant model in the fused lasso solution
  DA_result$predictive_MSE = locallyConstant_fitToAFit$predictive.MSE

  return(DA_result)



}
