#' @title Calculate the Target Variable of the "Fit-to-the-Fit" using Pseudo-data
#' @description This function obtains the target variable of the "fit-to-the-fit" using pseudo-data.
#'   It is part of the Interpret_SOFR() function. This function is not expected to be used directly by users.
#'   The function returns the target variable for each "fit-of-a-fit". It returns a list with length p,
#'   each component is a nx1 matrix/array.
#' @param X A list of functional covariate matrices, each matrix is n x length(t.grid), corresponding to a functional covariate.
#' @param beta.post An array of dimension p x S x length(t.grid), containing the posterior draws of the functional regression coefficients.
#' @param t.grid A numeric vector of the discrete points on which the functions X are observed.
#' @return A list with length p, each component is a nx1 matrix/array which is the target variable of the "fit-to-the-fit".
#' @noRd
get_h_bar <- function(X, beta.post, t.grid){
  grid.length = (max(t.grid) - min(t.grid))/ (length(t.grid)-1)
  p = length(X)
  h_bar = list()
  for (j in 1:p){
    X.j <- swapNAwith0(X[[j]])
    beta.j.pm <- apply(beta.post[j,,], 2, mean)
    h_bar[[j]]=X.j %*% as.matrix(beta.j.pm) * grid.length
  }
  return(h_bar)
}
