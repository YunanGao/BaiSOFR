% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotBASOFR.R
\name{plot_BASOFR_fit}
\alias{plot_BASOFR_fit}
\title{Plot the Posterior Inference of a BASOFR Fit}
\usage{
plot_BASOFR_fit(BASOFR_fit, groundTruth = NULL)
}
\arguments{
\item{BASOFR_fit}{A list object containing the result of a BASOFR model fit using the \code{fitBASOFR()} function. It is expected to contain the following components: \code{t.grid}, \code{beta.post.mean}, \code{beta.95CI}, \code{beta.50CI}, and \code{B.post}.}

\item{groundTruth}{An optional p x length(BASOFR_fit$t.grid) matrix, with each row representing a ground truth regression function. If provided, the ground truth will be plotted along with the posterior inference.}
}
\value{
A plot for each functional covariate's posterior inference including 95\% and 50\% credible intervals. The posterior mean is represented by a black line, while the 95\% and 50\% CI are represented by gray and light gray areas, respectively. If the ground truth is provided, it will be added to the plot as a red dashed line.
}
\description{
This function generates plots for the posterior inference of a Bayesian Adaptive Scalar-on-Function Regression (BASOFR) model fit.
}
\details{
This function provides a way to visualize the results of a BASOFR model fit. For each functional covariate in the model, it produces a plot of the posterior mean and credible intervals over the range of the function's input (t.grid). If the ground truth functions are provided, they are also included in the plot for comparison. This can be useful for assessing the quality of the model fit and for understanding the relationship between the functional covariates and the response variable.
}
\seealso{
\code{fitBASOFR()}
}
