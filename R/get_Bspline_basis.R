#' Create a B-spline basis
#'
#' This function creates a B-spline basis for a given set of knots.
#'
#' @param knots A numeric vector of knots for the B-spline basis.
#'        Defaults to a sequence of 100 equally spaced points spanning
#'        the range of `t.grid` in the global `data` object.
#'
#' @return A basis object for B-splines. The output is an object of class `basisfd`.
#'
#' @examples
#' knots <- seq(0, 1, length.out = 100)
#' basis <- get_Bspline_basis(knots)
#'
#' @importFrom fda create.bspline.basis
#' @noRd
get_Bspline_basis <- function(knots = seq(min(data$t.grid), max(data$t.grid), (max(data$t.grid)-min(data$t.grid))/100)){
  norder = 4
  nknots = length(knots)
  nbasis = length(knots) + norder - 2

  basis = fda::create.bspline.basis(c(min(knots),max(knots)), nbasis=nbasis,
                                    norder=norder, breaks=knots)
  return(basis)
}
