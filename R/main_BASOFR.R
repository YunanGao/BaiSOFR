#' fit a BASOFR Model
#'
#' This function fits a Bayesian Adaptive Scalar-on-Function Regression (BASOFR) model
#' using a Gibbs sampler.
#'
#' @param y A nx1 matrix or an array of the response variable.
#' @param X A list of matrices, each of dimension n x length(t.grid), representing observations of functional covariates. Missing values in the rows of X are allowed.
#' @param t.grid A numeric vector of the discrete points at which the functions X are observed.
#' @param z (Optional) A matrix of scalar covariates. Default is NULL. Note that we do not do any scaling or centering to preprocess this matrix. Therefore, it is advised that any necessary preprocessing, such as normalization or standardization, is performed on this matrix before passing it as an argument.
#' @param S An integer representing the number of iterations in the MCMC process. Default is 10000.
#' @param burnin An integer representing the number of initial iterations to discard (the "burn-in" period). Default is 10000.
#' @param save (Optional) Boolean to decide whether to save the MCMC samples. Default is FALSE.
#' @param save_file (Optional) File name to save the MCMC samples. Default is NULL.
#' @param num_knots (Optional) Number of knots for the B-spline basis. Default depends on the number of observations.
#' @param boundary_reg A logical flag indicating whether to implement boundary regularization for the regression function betas. Default is TRUE.
#' @param boundary_prior (Optional) A vector of length two representing the Inverse-Gamma prior parameters for the boundary basis coefficients variance. The first one is the shape parameter, and the second one is the rate parameter. Default is c(0.01, 0.01).
#' @param z_prior (Optional) A vector of length two representing the Inverse-Gamma prior parameters for the variance of the z-variables. The first one is the shape parameter, and the second one is the rate parameter. Default is c(0.01, 0.01).
#' @param sigma2_prior (Optional) A vector of length two representing the Inverse-Gamma prior parameters for the variance of the error term. The first element is the shape parameter, and the second is the rate parameter. Default is c(0.01, 0.01).
#' @param plot_basisExpansionX (Optional) Boolean to decide whether to plot the basis expansion of X. Default is FALSE. The B-spline basis selected by default typically captures the trend of the functional covariates while smoothing out irregular features. However, when the BASOFR model does not yield satisfactory results, it may be beneficial to examine the basis expansion of X by setting this parameter to TRUE. For each unique domain in the X covariates, a random observation is selected and the basis expansion for all functional covariates (where the coefficients are derived from ordinary least squares estimates) is displayed. If the B-spline basis expansion fails to adequately represent the original functional observations, consider increasing the number of knots of the B-spline using the `num_knots` argument. Conversely, if the B-spline basis expansion appears to overfit the original functional observations, consider decreasing the number of knots of the B-spline using the `num_knots` argument.
#'
#' @return Returns a list containing the following components:
#' \itemize{
#' \item{beta.post:}{An array of dimension p x S x length(t.grid), containing the posterior draws of the functional regression coefficients, where p is the number of functional covariates.}
#' \item{intercept_post:}{A length S vector of posterior draws of the intercept term.}
#' \item{Alpha.post:}{If z is not NULL, a S x p_scalar matrix of posterior draws of the scalar regression coefficients, where p_scalar is the number of scalar covariates.}
#' \item{sigma.2.post:}{A S x 1 array of posterior draws of the error term variance.}
#' \item{beta.post.mean:}{A p x length(t.grid) matrix, containing the posterior mean of the functional regression coefficients.}
#' \item{beta.95CI:}{An array of dimension p x length(t.grid) x 2, containing the 95% posterior credible interval of the functional regression coefficients.}
#' \item{beta.50CI:}{An array of dimension p x length(t.grid) x 2, containing the 50% posterior credible interval of the functional regression coefficients.}
#' \item{t.grid:}{A vector of discrete points where the functional regression coefficients are estimated on. Same as the argument t.grid.}
#' \item{other_parameters:}{A list containing the following components:
#' \itemize{
#' \item{DHS.post:}{A list containing the posterior draws of the DHS mu, phi, lambda terms}
#' \item{B.post:}{An array of dimension p x S x K containing the posterior draws of the B-spline basis coefficients of the functional regression coefficients, where K is the dimension of the B-spline basis we used to expand the regression functions and the functional covariates.}
#' \item{z.sigma.2.post:}{If z is not NULL, a S x p_scalar matrix containing the posterior draws of the variance of the scalar regression coefficients.}}}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate data for BASOFR
#' data = simulate_SOFR(n = 1500, p = 3, SNR = 5, num_z_continuous = 5, num_z_dummy = 2,
#'                      beta_types = c('spiky','smooth','stepwise'), domain_X = c(0, 1),
#'                      max_obs = 105,
#'                      subject_specific_domain = TRUE)
#'
#' # Scale the continuous z-variables to have std 1 and mean 0, and create the z-variables matrix
#' z = cbind(scale(data$Z_continuous), data$Z_dummy)
#'
#' # Store original means and standard deviations of z-variables for later use
#' z_continuous_means = colMeans(data$Z_continuous, na.rm = TRUE)
#' z_continuous_sds = apply(data$Z_continuous, 2, sd, na.rm = TRUE)
#'
#' # Fit the BASOFR model
#' BASOFR_fit <- fitBASOFR(y = data$y, X = data$X, t.grid = data$t.grid, z = z,
#'                         S = 1000, burnin = 1000)
#'
#' # Examine the items in the BASOFR_fit
#' str(BASOFR_fit)
#'
#' # Obtain the ground truth beta
#' groundTruth = matrix(unlist(data$groundTruth$beta), byrow=TRUE, nrow=3)
#'
#' # Plot the fit of the BASOFR model
#' plot_BASOFR_fit(BASOFR_fit, groundTruth)
#' }
fitBASOFR <- function(y, X, t.grid,
                      z = NULL,
                      S = 10000, burnin=10000,
                      num_knots=NULL,
                      boundary_reg = TRUE,
                      save = FALSE, save_file=NULL,
                      boundary_prior=c(0.01,0.01),
                      z_prior=c(0.01,0.01),
                      sigma2_prior = c(0.01,0.01),
                      plot_basisExpansionX = FALSE){

  # Extract n, p,
  n = length(y) # number of observation
  p = length(X) # number of functional covarites

  if(is.null(num_knots)){
    if(n > 100){
      num_knots = 51
    }else{
      num_knots = 30
    }
  }

  # 1. create B-spline
  knots <- seq(min(t.grid), max(t.grid), (max(t.grid)-min(t.grid))/(num_knots-1))
  Bspline_basis <- get_Bspline_basis(knots)

  # 2. apply basis expansion, scale the Z_variables
  basisExpansion_results <- basisExpansion_X(y = y, X = X,
                                             t.grid = t.grid,
                                             z = z,
                                             Bspline_basis = Bspline_basis,
                                             plot = plot_basisExpansionX)

  # 3. run the MCMC
  BASOFR_fit <- MCMC_BASOFR(y = basisExpansion_results$y,
                            X_feature = basisExpansion_results$X_feature,
                            Z = basisExpansion_results$z,
                            Bspline_basis = Bspline_basis,
                            t.grid = t.grid,
                            S = S, burnin = burnin,
                            save = save,
                            save_file = save_file,
                            boundary_reg = boundary_reg,
                            boundary_prior = boundary_prior,
                            z_prior = z_prior,
                            sigma2_prior = sigma2_prior)

  # Save the MCMC
  if(save == TRUE){
    save(BASOFR_fit, file = save_file)
  }

  return(BASOFR_fit)


}
