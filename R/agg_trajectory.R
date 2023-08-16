#' Get Aggregated Trajectories of Functional Covariates
#'
#' This function obtains the aggregated trajectories of the functional covariates X
#' over a given partition. The partition is defined by K and t.grid: it is going to be
#' evenly divided into K regions from min(t.grid) and max(t.grid).
#' Note: This is a helper function for the Interpret_SOFR() function and may not be
#' intended for direct use by package users.
#'
#' @param X A list of functional covariate matrices, each matrix is n x length(t.grid),
#' corresponding to a functional covariate.
#' @param t.grid A numeric vector of the discrete points on which the functions X are observed.
#' @param K The number of regions in the partition, indicating a partition with K even regions.
#'
#' @return A list of matrices with the same length as X. Each matrix is n x K, corresponding
#' to the aggregated trajectories of one of the functional covariates.
#'
#' @noRd
get_aggregated_trajectory <- function(X, t.grid,  K){
  p = length(X)
  n = dim(X[[1]])[1]
  partition.K <- seq(min(t.grid), max(t.grid), length.out=(K + 1))
  x.subdomain.aver.K <- lapply(X, function(xp) aggregate_for_one_func(xp, partition.K, t.grid) )
  return(x.subdomain.aver.K)
}

#' @title Aggregate For One Function
#' @description This is a helper function for `get_aggregated_trajectory()`.
#' It obtains the aggregated trajectories of a given functional covariate (Xp.func) over a given partition (partition.K).
#' This function is not expected to be directly used by users, it's mainly for internal usage within the package.
#'
#' @param Xp.func A numeric matrix of size n x length(t.grid). Represents one of the functional covariate observations.
#' @param partition.K A numeric vector indicating the partition of [min(t.grid), max(t.grid)].
#' @param t.grid A numeric vector representing the grid on which the functions are observed.
#'
#' @return A numeric matrix of size n x (length(partition.K) - 1), representing the aggregated trajectories for the given functional covariate over the given partition.
#'
#' @noRd
aggregate_for_one_func <- function(Xp.func, partition.K, t.grid){
  K = length(partition.K) - 1
  xp.subdomain.aver.K = matrix(NA, nrow=nrow(Xp.func), ncol=K)
  domain_K = matrix(NA, nrow = K, ncol = 2)
  grid.length = (max(t.grid) - min(t.grid))/ (length(t.grid)-1)
  for (i in 1:K){
    left_boundary = partition.K[i]
    right_boundary= partition.K[i+1]
    if(i != K){
      domain.k = which(t.grid >= left_boundary & t.grid < right_boundary)
    }else{domain.k = which(t.grid >= left_boundary)}

    domain_K[i,] = c(min(domain.k), max(domain.k))
    xp.subdomain.aver.K[,i] = apply(Xp.func[,domain.k],1, function(x) sum(x, na.rm=TRUE) * grid.length)
  }
  return(xp.subdomain.aver.K)

}




