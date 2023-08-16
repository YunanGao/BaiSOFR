#' Simulate functional predictors from Gaussian processes
#'
#' This function simulates functional predictors with or without seasonality
#' from Gaussian processes.
#'
#' @param n Number of observations
#' @param l Characteristic length-scale of the process, default is 0.01
#' @param x A vector of values where the functional predictors are observed,
#' same as the x.grids parameter in the simulate_BASOFR function
#' @param sigma_x Standard deviation for the Gaussian process kernel, default is 0.7
#' @param seasonality Logical indicating if seasonality should be included in the predictors, default is FALSE
#' @param P A parameter controlling the period of the sinusoidal component. Only effective when seasonality == TRUE, default is 365/(295-1)
#' @param A Amplitude of the sinusoidal component. Only effective when seasonality == TRUE, default is 1
#'
#' @return A matrix of simulated functional predictors
#' @import MASS
FunctionalPredictors <- function(n, l=0.01, x=seq(0,1,length.out = 201), sigma_x = 0.7, seasonality=FALSE,  P=365/(295-1), A=1 ) {

  ## phase parameters
  if(seasonality == TRUE){
    phi.phase = as.matrix(runif(n, min=0, max=2*pi))

    ## mean functions for all n subjects
    mu.x.matrix = apply(phi.phase, 1, function(phi) A*sin((2*pi*x/P) + phi))
  }

  ## Kernels
  obs.x = length(x)
  x_diff= abs(matrix(rep(x, obs.x), nrow=obs.x) - t(matrix(rep(x, obs.x), nrow=obs.x)))
  Kernel = sigma_x^2 * exp(-x_diff^2/(2*l^2))

  ## functional predictors
  if(seasonality == TRUE){
    x.sim = t(apply(mu.x.matrix, 2, function(mu.x) mvrnorm(n=1, mu=mu.x, Sigma=Kernel)))
    return(x.sim)
  }else{
    mu = rep(0, obs.x)
    x.sim = mvrnorm(n=n, mu=mu, Sigma=Kernel)
    return(x.sim)}

}

#' Plot the true beta functions
#'
#' @param beta_list A list of vectors representing the true beta functions.
#' @param x.grids A vector of values representing the points at which the
#'                beta functions are evaluated.
#'
#' @return A plot showing all beta functions, each in a different color,
#'         against the corresponding values of x.grids.
#' @importFrom graphics plot lines legend
#' @export
#' @examples
#' \dontrun{
#' sim_data = simulate_SOFR(n=1000, p=3, SNR = 7,
#'                            num_z_continuous = 0, num_z_dummy = 0,
#'                            x.grids =seq(0,1,length.out = 201),
#'                            beta_types = c('spiky', 'smooth', 'spiky'))
#' plot_beta(beta_list = sim_data$groundTruth$beta, x.grids = sim_data$groundTruth$x.grid)
#' }
plot_beta <- function(beta_list, x.grids) {
  p <- length(beta_list)

  # Define colors for each beta function
  colors <- rainbow(p)

  # Initialize an empty plot
  plot(1, 1, xlim = range(x.grids), ylim = range(sapply(beta_list, range)),
       type = "n", xlab = "t", ylab = bquote(beta(t)), main = expression(beta ~ 'function(s)'))

  # Add each beta function to the plot
  for (j in 1:p) {
    lines(x.grids, beta_list[[j]], col = colors[j], lwd = 2)
  }

  # Add a legend
  legend("topright", legend = sapply(1:p, function(i) bquote(beta[.(i)](t))),
                                     col = colors, lty = 1, lwd = 2)
}
