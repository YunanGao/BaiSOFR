#' Basis Expansion for Functional Covariates
#'
#' @description
#' This function B-spline performs basis expansion on functional covariates X using the given B-spline basis, with the
#' option of scaling the Z covariates. The result is a reordered
#' response and covariates according to the observed domain of each subject.
#'
#' @param y A matrix/array of the response variable.
#' @param X A list of functional covariates matrices.
#' @param t.grid A numeric vector of the discrete points on which the functions are observed.
#' @param z A matrix/array of scalar covariates, default to NULL.
#' @param Bspline_basis A basis object created by get_Bspline_basis() function
#' @param plot A logical value indicating whether to plot the basis expansion.
#'
#' @return
#' A list containing:
#' - X_feature: a list of basis expanded functional covariates
#' - y: reordered y
#' - z: reordered z
#' - ID: the reordered indices
#'
#' @importFrom fda eval.penalty eval.basis
#' @importFrom Matrix rankMatrix
#'
#'
#' @details
#'
#' This function is mainly used for data preprocessing in the fitBASOFR() function, and it is not
#' intended to be used directly by most users. However, in the case where fitBASOFR() does not yield satisfactory
#' results, users are advised to utilize basisExpansion_X() to inspect the basis expansion on X.
#' If the basis expansion fails to capture the trend, consider increasing the number of knots in the num_knots argument in the fitBASOFR() function.
#' If the basis expansion overfits the X trend, consider decreasing the number of knots in the num_knots argument in the fitBASOFR() function.
#'
#' @examples
#' \dontrun{
#' # simulate data
#' data = simulate_SOFR(n = 100, p = 2, SNR = 5, num_z_continuous = 1, num_z_dummy = 1,
#'                        x.grids = seq(0, 1, length.out = 100), beta_types = c('spiky', 'smooth'))
#' # create B-spline
#' num_knots=50
#' knots <- seq(min(data$t.grid), max(data$t.grid),
#'              (max(data$t.grid)-min(data$t.grid))/(num_knots-1))
#' Bspline_basis <- get_Bspline_basis(knots)
#'
#' # use basisExpansion_X to inspect the basis expansion
#' basisExpansion_X(y=data$y, X=data$X, t.grid=data$t.grid, z_continuous=data$z_continuous,
#'                  z_dummy= data$z_dummy, Bspline_basis=Bspline_basis, plot= TRUE, scale=TRUE)
#' }
#'
#' @export
basisExpansion_X <- function(y = y,
                             X = X, t.grid = t.grid,
                             z = z,
                             #z_continuous = z_continuous, z_dummy = z_dummy,
                             Bspline_basis = Bspline_basis,
                             plot = visualize_basisExpansion_X
                             #,scale = scale_z
){
  n  = dim(X[[1]])[1]
  p = length(X)

  # Obtain the subject specific domain
  domain_counts = sapply(X, function(x) apply(x, 1, function(row) sum(!is.na(row))))
  subject_specific_domain = domain_counts[,1]
  subject_domain_dist <- table(subject_specific_domain)

  # Check if the functional covariates have the same missingness
  if(!all(domain_counts == domain_counts[,1])){
    stop('Functional covariates do not have the same missingness pattern at certain rows.')
  }

  basismat = fda::eval.basis(t.grid, Bspline_basis)

  # define id_reorder and X_func_feature as empty lists
  X_func_feature = list()
  id_reorder = list()

  for(domain in unique(subject_specific_domain)){

    # the basismat
    basismat.domain = basismat[1:domain,]
    if(domain < length(t.grid)){
      basismat.domain.fullrank = basismat.domain[,1:Matrix::rankMatrix(basismat.domain)]
    }else{
      basismat.domain.fullrank = basismat
    }


    # the ID of X with the same domain
    domain_ID <- which(subject_specific_domain == domain)


    # Subset of X for a specific domain
    X_domain = lapply(X, function(x) x[domain_ID, 1:domain])

    # Perform basis expansion
    if(length(domain_ID) > 1){
      basis.coeff.domain = lapply(X_domain, function (xp) t(tcrossprod(tcrossprod(solve(crossprod(basismat.domain.fullrank)),
                                                                                  basismat.domain.fullrank), xp)))
    }else{basis.coeff.domain = lapply(X_domain, function (xp) t(tcrossprod(tcrossprod(solve(crossprod(basismat.domain.fullrank)),
                                                                                      basismat.domain.fullrank), t(as.matrix(xp)))))}


    if(plot == TRUE){
      # Visualize the basis expansion
      reproduce.domain = lapply(basis.coeff.domain, function(x_star_p) as.matrix(tcrossprod(x_star_p, basismat.domain.fullrank)))

      if(length(domain_ID) > 1){
        original.i = lapply(X, function(xp) t(as.matrix(t(xp[domain_ID, 1:domain] ))))

      }else{
        original.i = lapply(X, function(xp) (as.matrix(t(xp[domain_ID, 1:domain] ))))
      }
      par(mfrow=c(2,2))
      for (j in 1:p){
        index = sample(1:length(domain_ID), 1)
        main = paste0('Basis expansion reproduce curve, \n ', 'domain: [', t.grid[1], ',' ,t.grid[domain], ']', ' -- X',j )
        plot(t.grid[1:domain], original.i[[j]][index,], col='red', type='l', main=main, xlab='t', ylab = 'X and reproduced X' )
        lines(t.grid[1:domain], reproduce.domain[[j]][index,])
      }
    }

    # Impute 0 for parts of basis coefficients that do not exist
    basis.coeff.domain.imput0 = lapply(basis.coeff.domain, function(xp) {
      if(ncol(xp) < ncol(basismat)) {
        xp = cbind(xp, matrix(0, nrow=nrow(xp), ncol=(ncol(basismat)-ncol(xp))))
      }
      return(xp)
    })

    # obtain domain-specific J matrix
    J.i = fda::eval.penalty(basisobj = Bspline_basis, rng=c(t.grid[1], t.grid[domain]))

    # obtain basis coefficient  %*% J.i
    X.func.feature.domain = lapply(basis.coeff.domain.imput0, function(xp_star) tcrossprod(xp_star, J.i))


    # Store the domain_ID and X.func.feature.domain in lists
    id_reorder[[as.character(domain)]] = domain_ID
    X_func_feature[[as.character(domain)]] = X.func.feature.domain


  }

  # Combine each component in the id_reorder list head to toe
  id_reorder_all = unlist(id_reorder)

  # Rbind all of the X_func_feature components
  X_func_feature_all = lapply(1:p, function(i) do.call(rbind, lapply(X_func_feature, `[[`, i)))

  # Reorder y
  if (is.null(y)) {
    stop('Response variable y is NULL. Cannot proceed with BASOFR.')
  } else {
    y <- y[id_reorder_all, , drop = FALSE]
  }

  # Variables to store the original mean and sd
  mean_continuous = NULL
  sd_continuous = NULL
  mean_dummy = NULL
  sd_dummy = NULL

  # # Reorder z_continuous if not NULL
  # if (!is.null(z_continuous)) {
  #   z_continuous <- z_continuous[id_reorder_all, , drop = FALSE]
  #   if(scale) {
  #     mean_continuous <- colMeans(z_continuous, na.rm = TRUE)
  #     sd_continuous <- apply(z_continuous, 2, sd, na.rm = TRUE)
  #     z_continuous <- scale(z_continuous, center = mean_continuous, scale = sd_continuous)
  #   }
  # }
  #
  # # Reorder z_dummy if not NULL
  # if (!is.null(z_dummy)) {
  #   z_dummy <- z_dummy[id_reorder_all, , drop = FALSE]
  #   if(scale) {
  #     mean_dummy <- colMeans(z_dummy, na.rm = TRUE)
  #     sd_dummy <- apply(z_dummy, 2, sd, na.rm = TRUE)
  #     z_dummy <- scale(z_dummy, center = mean_dummy, scale = sd_dummy)
  #   }
  # }

  # Reorder z if not NULL
  if(!is.null(z)){
    z <- z[id_reorder_all, , drop = FALSE]
  }


  # Return a list of reordered elements
  return(list(X_feature = X_func_feature_all,
              y = y,
              z = z,
              # z_continuous = z_continuous,
              # z_dummy = z_dummy,
              #mean_continuous = mean_continuous, sd_continuous = sd_continuous,
              #mean_dummy = mean_dummy, sd_dummy = sd_dummy,
              ID = id_reorder_all))
}
