#' Fit-of-the-fit: Locally constant estimates and mean square errors
#'
#' This function performs the 'fit-of-the-fit' process to obtain locally constant estimates
#' of different complexity. It also calculates the empirical and predictive mean square errors.
#'
#' @param h_bar A list of target variables of the "Fit-to-the-Fit" using Pseudo-data. This should have the same length as X and is a returned object from get_h_bar() function.
#' @param x.subdomain.aver.K A list of matrices with the same length as X. Each matrix is of size n x K, corresponding to the aggregated trajectories of each of the functional covariates. This parameter is a returned object from the get_aggregated_trajectory() function.
#' @param beta.post A 3D array of the posterior samples of the functional regression coefficients.
#' @param Alpha.post A matrix of the posterior samples of the coefficients for z.
#' @param intercept.post A vector of the posterior samples of the intercepts.
#' @param sigma2.post A vector of the posterior samples of the error term variance.
#' @param t.grid A numeric vector of the discrete points on which the functions X are observed.
#' @param y A numeric vector of the response variable.
#' @param X A list of matrices, each of which represents a functional covariate.
#' @param z A matrix or NULL. If it is a matrix, it represents the control variable(s) z.
#' @param plot_all_locally_constant A logical value indicating whether to plot the locally constant estimates.
#'
#' @return A list containing four elements:
#'         - predictive.MSE: the predictive mean square errors,
#'         - point.MSE: the empirical mean square errors,
#'         - LocallyConstant.est: the locally constant estimates,
#'         - model.sum: the unique locally constant models for each functional covariate.
#'
#' @note This function is a part of Interpret_SOFR and not typically called by the user directly.
#'
#' @importFrom genlasso fusedlasso
#'
#' @examples
#' # The function is usually not called directly, so examples are not provided.
#'
get_locallyconstant <- function(h_bar, x.subdomain.aver.K,
                                beta.post, Alpha.post, intercept.post, sigma2.post,
                                t.grid, y, X, z, plot_all_locally_constant ){

  n =  dim(x.subdomain.aver.K[[1]])[1]
  p = length(x.subdomain.aver.K)
  S = dim(beta.post)[2]
  K = dim(x.subdomain.aver.K[[1]])[2]

  grid.length = (max(t.grid) - min(t.grid))/ (length(t.grid)-1)
  partition.K = seq(min(t.grid), max(t.grid), length.out=(K + 1))
  X = lapply(X, function(Xp) swapNAwith0(Xp))

  # Check if both Alpha.post and z are non-null
  if (!is.null(Alpha.post) && !is.null(z)) {
    # calculate Alpha.post.mean
    Alpha.post.mean = apply(Alpha.post, 2, mean)
  } else if ((as.numeric(is.null(Alpha.post)) + as.numeric(is.null(z))) == 1) {
    # if one of Alpha.post, z is null and the other is non-null, then there is something wrong going on
    stop("Alpha.post and z should be either both non-null or both null.")
  }

  intercept.post.mean = mean(intercept.post)
  beta.post.mean = array(NA, dim=c(p, dim(beta.post)[3]))
  beta.95CI = array(NA,  dim=c(p, dim(beta.post)[3], 2))
  for (j in 1:p) {
    beta.post.mean[j, ] = apply(beta.post[j,,], 2, mean)
    beta.95CI[j,,] = t(apply(beta.post[j,,], 2, function(x) quantile(x, c(0.025,0.975))))
  }


  domain_K = matrix(NA, nrow = K, ncol = 2)
  for (i in 1:K){
    left_boundary = partition.K[i]
    right_boundary= partition.K[i+1]
    if(i != K){
      domain.k = which(t.grid >= left_boundary & t.grid < right_boundary)
    }else{domain.k = which(t.grid >= left_boundary)}
    domain_K[i,] = c(min(domain.k), max(domain.k))

  }

  # penalization matrix
  D1 = diff(diag(K), diff=1)

  fusedlasso_output = list()

  # first fit fusedlasso
  for (j in 1:p){fusedlasso_output[[j]] = genlasso::fusedlasso(y = h_bar[[j]], X = x.subdomain.aver.K[[j]], D=D1); gc()}

  # store the unique models for each functional covariates
  model.sum = lapply(fusedlasso_output, function(fusedlasso_output_j) get_unique_locallyConstant(fusedlasso_output_j, K))

  # then refit and store the point and predictive MSE
  posterior_predictive_dist <- get_posterior_predictive(X, t.grid, beta.post, sigma2.post, add_noise = TRUE)

  predictive.MSE = list()
  point.MSE = list()
  stepwise.est = list()


  for (j in 1:p){
    cat(paste("Start of 'fit-of-the-fit' process for beta", j, " to obtain locally constant estimates.\n"))
    model.sum.j = model.sum[[j]]
    x.subdomain.aver.K.j = x.subdomain.aver.K[[j]]
    h_bar.j = h_bar[[j]]

    posterior_predictive_dist.j = posterior_predictive_dist[[j]]

    if(p > 1){
      X_j_else = X[-j]
      beta.post.mean.j_else = beta.post.mean[-j,]
      func.pred.j_else = rep(0,n)
      for (k in 1:(p-1)){func.pred.j_else = func.pred.j_else + X_j_else[[k]] %*% beta.post.mean.j_else[k,] * grid.length}
      if(!is.null(z) && !is.null(Alpha.post.mean)){
        point.MSE.target.j = y - z %*% Alpha.post.mean - func.pred.j_else - intercept.post.mean
      }else if (is.null(z)){
        point.MSE.target.j = y - func.pred.j_else - intercept.post.mean
      }
    }else if(p == 1){
      if(!is.null(z) && !is.null(Alpha.post.mean)){
        point.MSE.target.j = y - z %*% Alpha.post.mean - intercept.post.mean
      }else if (is.null(z)){
        point.MSE.target.j = y - intercept.post.mean
      }
    }


    predictive.MSE.j <- matrix(NA, nrow=nrow(model.sum.j), ncol=S)
    point.MSE.j <- rep(NA,nrow(model.sum.j))
    stepwise.j <- matrix(NA, nrow=nrow(model.sum.j), ncol=length(t.grid))
    par(mfrow=c(2,2))
    for(i in 1:nrow(model.sum.j)){

      model.i = model.sum.j[i,]
      n.cp.i = model.i$num.cp
      domain.model.i =matrix(sort(c(1, which(model.i[-1]!=0), which(model.i[-1]!=0)+1, K)), nrow = n.cp.i+1, byrow=TRUE)
      x.i = apply(domain.model.i, 1, function(x) if(x[1] < x[2]){rowSums(x.subdomain.aver.K.j[, (x[1]:x[2])])}else if(x[1]==x[2]){x.subdomain.aver.K.j[, x[1]]}   )
      lm.i = lm(h_bar.j~x.i-1)
      coeff.i = lm.i$coefficients
      times = apply(t(apply(domain.model.i,1, function(x) c(min(domain_K[(x[1]:x[2]),]), max(domain_K[(x[1]:x[2]),])))),1, function(region) region[2] -region[1] + 1)
      if(sum(times) != length(t.grid)){
        print('warning! approx beta length!= obs x length')
      }
      stepwise.j[i, ] = rep(coeff.i, times)
      main = paste('beta', j, ': ncp', n.cp.i, ', model', i)

      # please fill in:
      if(plot_all_locally_constant == TRUE){
        # plot the locally constant estimate

        main = bquote(paste(beta[.(j)], " Locally Constant Model ", .(i), ", ", .(n.cp.i), " changing pts"))

        plot(t.grid, beta.post.mean[j, ], type = "l", lty=3,
             main = main,
             xlab = "t", ylab = bquote(beta[.(j)](t))
             , ylim = range(beta.95CI[j,,]))
        abline(h = 0, lwd = 2, lty = 2, col = "lightgray")

        # Add 95% CI
        polygon(c(t.grid, rev(t.grid)), c(beta.95CI[j,,2], rev(beta.95CI[j,,1])), col = "gray", border = NA)


        lines(t.grid, beta.post.mean[j,], lwd = 2)

        # Add locally constant estimate
        lines(t.grid, stepwise.j[i,], lwd=3, col='blue', lty=1)

        legend("topright", legend = c("95% CI", "Posterior Mean", "Locally Constant Estimate"),
               fill = c("gray", "black", "blue"), bty = "n", cex = 0.8)

      }


      predictive.MSE.j[i,] = apply( posterior_predictive_dist.j, 2, function(x)  mean((x-lm.i$fitted.values)^2))
      point.MSE.j[i] = mean((    point.MSE.target.j - x.i %*% lm.i$coefficients)^2)
    }


    predictive.MSE[[j]] = predictive.MSE.j
    point.MSE[[j]] = point.MSE.j
    stepwise.est[[j]] = stepwise.j
    cat(paste("End of 'fit-of-the-fit' process for beta", j, " to obtain locally constant estimates.\n"))
  }
  return(list(predictive.MSE = predictive.MSE,
              point.MSE = point.MSE,
              LocallyConstant.est = stepwise.est,
              model.sum = model.sum))

}
