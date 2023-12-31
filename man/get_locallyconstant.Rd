% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getLocallyConstant.R
\name{get_locallyconstant}
\alias{get_locallyconstant}
\title{Fit-of-the-fit: Locally constant estimates and mean square errors}
\usage{
get_locallyconstant(
  h_bar,
  x.subdomain.aver.K,
  beta.post,
  Alpha.post,
  intercept.post,
  sigma2.post,
  t.grid,
  y,
  X,
  z,
  plot_all_locally_constant
)
}
\arguments{
\item{h_bar}{A list of target variables of the "Fit-to-the-Fit" using Pseudo-data. This should have the same length as X and is a returned object from get_h_bar() function.}

\item{x.subdomain.aver.K}{A list of matrices with the same length as X. Each matrix is of size n x K, corresponding to the aggregated trajectories of each of the functional covariates. This parameter is a returned object from the get_aggregated_trajectory() function.}

\item{beta.post}{A 3D array of the posterior samples of the functional regression coefficients.}

\item{Alpha.post}{A matrix of the posterior samples of the coefficients for z.}

\item{intercept.post}{A vector of the posterior samples of the intercepts.}

\item{sigma2.post}{A vector of the posterior samples of the error term variance.}

\item{t.grid}{A numeric vector of the discrete points on which the functions X are observed.}

\item{y}{A numeric vector of the response variable.}

\item{X}{A list of matrices, each of which represents a functional covariate.}

\item{z}{A matrix or NULL. If it is a matrix, it represents the control variable(s) z.}

\item{plot_all_locally_constant}{A logical value indicating whether to plot the locally constant estimates.}
}
\value{
A list containing four elements:
- predictive.MSE: the predictive mean square errors,
- point.MSE: the empirical mean square errors,
- LocallyConstant.est: the locally constant estimates,
- model.sum: the unique locally constant models for each functional covariate.
}
\description{
This function performs the 'fit-of-the-fit' process to obtain locally constant estimates
of different complexity. It also calculates the empirical and predictive mean square errors.
}
\note{
This function is a part of Interpret_SOFR and not typically called by the user directly.
}
\examples{
# The function is usually not called directly, so examples are not provided.

}
