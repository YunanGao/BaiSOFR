#' @title Plot the Posterior Inference of a BASOFR Fit
#'
#' @description This function generates plots for the posterior inference of a Bayesian Adaptive Scalar-on-Function Regression (BASOFR) model fit.
#'
#' @param BASOFR_fit A list object containing the result of a BASOFR model fit using the `fitBASOFR()` function. It is expected to contain the following components: `t.grid`, `beta.post.mean`, `beta.95CI`, `beta.50CI`, and `B.post`.
#' @param groundTruth An optional p x length(BASOFR_fit$t.grid) matrix, with each row representing a ground truth regression function. If provided, the ground truth will be plotted along with the posterior inference.
#'
#' @return A plot for each functional covariate's posterior inference including 95% and 50% credible intervals. The posterior mean is represented by a black line, while the 95% and 50% CI are represented by gray and light gray areas, respectively. If the ground truth is provided, it will be added to the plot as a red dashed line.
#'
#' @details This function provides a way to visualize the results of a BASOFR model fit. For each functional covariate in the model, it produces a plot of the posterior mean and credible intervals over the range of the function's input (t.grid). If the ground truth functions are provided, they are also included in the plot for comparison. This can be useful for assessing the quality of the model fit and for understanding the relationship between the functional covariates and the response variable.
#'
#' @seealso `fitBASOFR()`
#' @export
plot_BASOFR_fit <- function(BASOFR_fit, groundTruth = NULL){
  t.grid = BASOFR_fit$t.grid
  beta.post.mean = BASOFR_fit$beta.post.mean
  beta.post.95CI = BASOFR_fit$beta.95CI
  beta.post.50CI = BASOFR_fit$beta.50CI

  p = dim(BASOFR_fit$beta.post)[1]

  par(mfrow=c(1,1))
  for (j in 1:p){
    if(!is.null(groundTruth)){
      ylim = range(c(beta.post.95CI[j,,], beta.post.50CI[j,,], groundTruth[j,]))
    }else{
      ylim = range(c(beta.post.95CI[j,,], beta.post.50CI[j,,]))
    }

    plot(t.grid, beta.post.mean[j, ], type = "l",
         main = bquote(beta[.(j)]~"Posterior Inference"),
         xlab = "t", ylab = bquote(beta[.(j)](t)),
         ylim = ylim)
    abline(h = 0, lwd = 2, lty = 2, col = "lightgray")

    # Add 95% and 50% CI
    polygon(c(t.grid, rev(t.grid)), c(beta.post.95CI[j,,2], rev(beta.post.95CI[j,,1])), col = "gray", border = NA)
    polygon(c(t.grid, rev(t.grid)), c(beta.post.50CI[j,,2], rev(beta.post.50CI[j,,1])), col = "gray80", border = NA)

    lines(t.grid, beta.post.mean[j,], lwd = 2)

    if(!is.null(groundTruth)){
      # add red line as ground truth
      lines(t.grid, groundTruth[j,], col = "red", lty = 3, lwd=3)
      # Add a legend
      legend("topright", legend = c("95% CI", "50% CI", "Posterior Mean", "Ground Truth"),
             fill = c("gray", "gray80", "black", "red"), bty = "n", cex = 0.8)
    }

    # Add a legend
    legend("topright", legend = c("95% CI", "50% CI", "Posterior Mean"),
           fill = c("gray", "gray80", "black"), bty = "n", cex = 0.8)


  }
}


