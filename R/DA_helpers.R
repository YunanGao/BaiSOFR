#' @title Substitute NA with 0 in a Matrix
#' @description This function is a helper function for get_h_bar(). It is used to replace all NA values in a given matrix with 0.
#' @param Xp A matrix with potential NA values.
#' @return The input matrix with all NA values replaced by 0.
#' @noRd
swapNAwith0 <- function(Xp){
  Xp.0 = Xp
  Xp.0[is.na(Xp)] <- 0
  return(Xp.0)
}

#' @title Posterior Predictive Distribution Calculation for Functional Regression Terms
#' @description This function calculates the posterior distribution of key functional regression terms. It is part of the Interpret_SOFR() function, and therefore not intended for direct use by users. The output of this function is used to calculate the predictive mean square error of the locally constant estimates.
#' @param X A list of functional covariate matrices, each matrix is n x length(t.grid), corresponding to a functional covariate.
#' @param t.grid A numeric vector of the discrete points on which the functions X are observed.
#' @param beta.post An array of dimension p x S x length(t.grid), containing the posterior draws of the functional regression coefficients.
#' @param sigma2.post A vector representing the posterior distribution of the error variance.
#' @param add_noise (Optional) A Boolean flag indicating whether to add a noise term to the calculation. Default is TRUE.
#' @return Returns a list of the same length as X, with each component being a n x S matrix. This list represents the posterior distribution of the key functional regression term.
#' @noRd
get_posterior_predictive <- function(X, t.grid, beta.post, sigma2.post, add_noise = TRUE){
  p = length(X)
  n = dim(X[[1]])[1]
  grid.length = (max(t.grid) - min(t.grid))/ (length(t.grid)-1)
  X.func.with0 <- lapply(X, function(Xp) swapNAwith0(Xp))

  posterior_predictive = list()
  for (j in 1:p){
    X.func.p <- X.func.with0[[j]]
    beta.post.p <- beta.post[j,,]
    posterior_predictive[[j]] = X.func.p %*% t(beta.post.p) * grid.length
    gc()
  }

  if(add_noise == TRUE){
    noise_term = apply(sigma2.post, 1, function(x) rnorm(n, mean=0, sd=sqrt(x)))
    posterior_predictive.1 <- lapply(posterior_predictive, function(dist) dist + noise_term)
    gc()
    return(posterior_predictive.1)
  }else{
    return(posterior_predictive)
  }
}

#' Get Unique Locally Constant Functions
#'
#' This function keeps the models with different numbers of change points (cps),
#' or cps at different locations. It is a helper function of the `get_locallyconstant()`
#' function, which is part of the `Interpret_SOFR()` function, and is not expected to
#' be used directly by the users.
#'
#' @param fusedlasso_output_j An object returned by `genlasso::fusedlasso()`.
#' @param K The maximum number of locally constant pieces.
#'
#' @return A matrix indicating the unique locally constant functions. The output matrix
#' contains the number of change points and their locations.
#'
#' @noRd
get_unique_locallyConstant <- function(fusedlasso_output_j, K){
  num.cp <- apply(diff(diag(K), diff=1) %*% fusedlasso_output_j$beta, 2, function(x) sum(x!=0))
  cp.location <- (diff(diag(K), diff=1) %*% fusedlasso_output_j$beta != 0) * 1
  model.sum.each <- unique(as.data.frame(cbind(num.cp, t(cp.location))))
  model.sum.each <- model.sum.each[order(model.sum.each$num.cp),]
  return(model.sum.each)
}
