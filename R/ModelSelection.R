#' Model Selection and Visualization for Locally Constant Estimates
#'
#' This function performs the following:
#' 1. Selects the acceptable family of locally constant estimates for each functional regression coefficient.
#' 2. Plots the simplest member of the acceptable family for each functional regression coefficient.
#'
#' It is a part of Interpret_SOFR function, so typically it's not called directly by the user.
#'
#' @param locallyConstant_fitToAFit An object returned by the `get_locallyconstant()` function.
#' @param eps_level A numeric value for selecting an acceptable family. It determines the threshold for predictive performance of models under the Bayesian SOFR model.
#' @param t.grid A numeric vector of the discrete points on which the functions X are observed.
#' @param beta.post A 3D array of the posterior samples of the functional regression coefficients.
#'
#' @return A list of model summaries for each functional regression coefficient. Each summary includes:
#'         - ell_min_in: Best model by in-sample MSE,
#'         - post.dmse: Percent difference in predictive MSE relative to the "best" model,
#'         - dmse: Percent difference in MSE relative to the "best" model,
#'         - ci_dmse: 100(1 - 2*eps_level)% confidence interval of dmse,
#'         - simplest_accFam_index: Index of the simplest models in the acceptable family,
#'         - simplest_accFam: Matrix of simplest models in the acceptable family.
#'
#' @note This function is a part of Interpret_SOFR and not typically called by the user directly.
#'
#' @examples
#' # This function is usually not called directly, so examples are not provided.
model_selection <- function(locallyConstant_fitToAFit, eps_level = 0.1,
                            t.grid, beta.post){
  p = length(locallyConstant_fitToAFit$predictive.MSE)
  eta_level = 0.0 # % percent of best model
  model.summary = list()

  beta.post.mean = array(NA, dim=c(p, dim(beta.post)[3]))
  beta.post.95CI = array(NA,  dim=c(p, dim(beta.post)[3], 2))
  for (j in 1:p) {
    beta.post.mean[j, ] = apply(beta.post[j,,], 2, mean)
    beta.post.95CI[j,,] = t(apply(beta.post[j,,], 2, function(x) quantile(x, c(0.025,0.975))))
  }

  for (j in 1:p) {
    model.summary.j = list()
    MSE.point = locallyConstant_fitToAFit$point.MSE[[j]]
    MSE.predictive= locallyConstant_fitToAFit$predictive.MSE[[j]]
    model.sum = locallyConstant_fitToAFit$model.sum[[j]]
    locallyConstant_est = locallyConstant_fitToAFit$LocallyConstant.est[[j]]
    beta.post.mean.j = beta.post.mean[j,]
    beta.post.95CI.j = beta.post.95CI[j,,]

    # Best model by in-sample MSE:
    ell_min_in = which.min(MSE.point)
    model.summary.j$ell_min_in = ell_min_in

    # Percent difference in predictive MSE relative to the "best" model:
    post.dmse = apply(t(MSE.predictive), 2 , function(pmse)
      100*(pmse- t(MSE.predictive)[,ell_min_in])/t(MSE.predictive)[,ell_min_in])
    model.summary.j$post.dmse = post.dmse

    # Percent difference in MSE relative to the "best" model:
    dmse = 100*(MSE.point - MSE.point[ell_min_in])/MSE.point[ell_min_in]
    model.summary.j$dmse = dmse

    # 100(1 - 2*eps_level)%
    ci_dmse = t(apply(post.dmse, 2, quantile, c(eps_level, 1 - eps_level)))
    model.summary.j$ci_dmse = ci_dmse

    # one of the Smallest model(s) for which the lower CI includes zero:
    ell_eps_in = min(which(colMeans(post.dmse <= eta_level) >= eps_level))

    # num of Cps in the ell_eps_in-th model
    ncp.ell_eps_in = model.sum$num.cp[ell_eps_in]

    # all the models' index with the same num of Cps
    same.size.index.in = which(model.sum$num.cp == ncp.ell_eps_in)

    # the set of models we are intrested in
    models.int.in = same.size.index.in[same.size.index.in >= ncp.ell_eps_in]

    # the smallest model(s) for which the lower CI includes zero:
    ell_eps_in_set = intersect(which(colMeans(post.dmse <= eta_level) >= eps_level),
                               models.int.in)
    model.summary.j$simplest_accFam_index = ell_eps_in_set

    model.summary.j$simplest_accFam = as.matrix(locallyConstant_est[ell_eps_in_set,])

    model.index = 1:dim(model.sum)[1]

    main = bquote(paste('Acceptable Family Selection for ',beta[.(j)]))
    par(mfrow=c(1,1))
    plot(model.index,
         colMeans(post.dmse),type='p', ylim = range(ci_dmse, 0), lwd=3,
         xlab = 'Model Index', ylab = 'Difference in predictive loss (%)',main=main,
         cex.lab=1, cex.axis= 1)
    abline(h = eta_level, lwd=3, lty=6)
    abline(v = ell_eps_in_set[1], lwd = 8, col="grey50")
    abline(v = ell_min_in, lwd=6, col="grey80", lty=2)
    suppressWarnings(
      arrows(model.index[-ell_min_in], ci_dmse[-ell_min_in,1],
             model.index[-ell_min_in], ci_dmse[-ell_min_in,2],
             length=0.05, angle=90, code=3, lwd=5)
    )

    # Add legend
    legend("topright",
           legend = c("95% CI of % increase in\nPredictive MSE compared to\nempirical MSE minimizer",
                      "Empirical MSE minimizer",
                      "Simplest member in acceptable family"),
           fill = c("black", "grey80", "grey50"),
           bty = "n",
           cex = 0.6, y.intersp = 0.5)


    par(mfrow=c(1,length(ell_eps_in_set)))
    for (i in (1:length(ell_eps_in_set))){
      main = bquote(paste(beta[.(j)], ': Simplest Member in Acceptable Family '))
      plot(t.grid, beta.post.mean.j, col="blue",lty=6,lwd=3, type="l",xlab="t",
           ylab = bquote(beta[.(j)](t)), main=main, ylim = c(min(beta.post.95CI.j[,1]),max(beta.post.95CI.j[,2])))

      polygon(c(t.grid, rev(t.grid)), c(beta.post.95CI[j,,2], rev(beta.post.95CI[j,,1])), col = "gray", border = NA)


      lines(t.grid, beta.post.mean.j, lwd = 2)

      # Add locally constant estimate
      lines(t.grid, model.summary.j$simplest_accFam[,i], lwd=3, col='blue', lty=1)

      legend("topright", legend = c("95% CI", "Posterior Mean", "Locally Constant Estimate"),
             fill = c("gray", "black", "blue"), bty = "n", cex = 0.8)
      abline(h=0, lty=6, lwd=1)

    }
    model.summary[[j]] = model.summary.j
  }
  return(model.summary)
}
