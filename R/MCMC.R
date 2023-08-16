#' @title Run the MCMC for the BASOFR Model
#'
#' @description This function runs the MCMC for the Bayesian Adaptive Scalar-on-function regression (BASOFR) model.
#' This function is designed to be used internally by fitBASOFR() function and it is not expected to be called by users directly.
#'
#' @param y The response variable, already processed by the basisExpansion_X() function.
#' @param X_feature The functional covariates after basis expansion and multiplication by the J matrix. This list is obtained by the basisExpansion_X() function
#' @param Z Scalar covariates, already processed by the basisExpansion_X() function.
#' @param Bspline_basis B-spline basis matrix.
#' @param t.grid A sequence of discrete points at which the functional covariates are observed.
#' @param boundary_reg A boolean flag to indicate whether boundary regularization should be implemented.
#' @param S The number of iterations in the MCMC chain.
#' @param burnin The number of initial MCMC iterations to discard as "burn-in".
#' @param save A boolean flag to indicate whether the MCMC output should be saved.
#' @param save_file The name of the file where the MCMC output should be saved.
#' @param boundary_prior The prior for the boundary parameter.
#' @param z_prior The prior for the z variables.
#' @param sigma2_prior The prior for the sigma2 parameter.
#'
#' @return A list containing the MCMC outputs.
#'
#' @importFrom fda eval.penalty eval.basis
#'
#' @details This function fits the Bayesian Adaptive Scalar-on-Function Regression (BASOFR) model via MCMC. The function receives inputs that have already been processed by basisExpansion_X(). The output is a list of MCMC chains for the parameters of interest.
#' @seealso \code{\link{fitBASOFR}} for an external function that uses this function internally.
#' @examples
MCMC_BASOFR <- function(y,
                        X_feature,
                        Z, #z_continuous, #z_dummy,
                        Bspline_basis,
                        t.grid,
                        S = 10000,
                        burnin = 10000,
                        save = FALSE,
                        save_file,
                        boundary_reg = TRUE,
                        boundary_prior = c(0.01,0.01),
                        z_prior = c(0.01,0.01),
                        sigma2_prior = c(0.01, 0.01)
){

  n  = dim(X_feature[[1]])[1]
  p = length(X_feature)

  # Inverse gamma prior on the boundary basis coefficients
  if(boundary_reg == TRUE){
    a.lambda = boundary_prior[1]  # shape parameter
    b.lambda = boundary_prior[2]  # rate parameter
  }

  # Inverse gamma prior on the alpha's variance
  aj0 = z_prior[1] # shape parameter
  bj0 = z_prior[2] # rate parameter

  # Inverse gamma prior on sigma^2
  a0 = sigma2_prior[1] # shape parameter
  b0 = sigma2_prior[2] # rate parameter

  # B-spline basis estimated on the discrete points we collect X
  basismat <- fda::eval.basis(t.grid, Bspline_basis)
  K = ncol(basismat)

  # 2nd-order difference matrix
  D_d = diff(diag(K), diff=2)

  # # create the Z variables, which should be a column bind of both of the z_continous and z_dummy
  # Z = NULL; p_scalar = 0; intercept_Z = matrix(rep(1, n), ncol = 1)
  # if (!is.null(z_continuous) && !is.null(z_dummy)) {
  #   Z = cbind(z_continuous, z_dummy); p_scalar = ncol(Z); intercept_Z = cbind(intercept_Z, Z)
  # } else if (!is.null(z_continuous)) {
  #   Z = z_continuous; p_scalar = ncol(Z); intercept_Z = cbind(intercept_Z, Z)
  # } else if (!is.null(z_dummy)) {
  #   Z = z_dummy; p_scalar = ncol(Z); intercept_Z = cbind(intercept_Z, Z)
  # }


  intercept_Z = matrix(rep(1, n), ncol = 1)
  if(!is.null(Z)){
    p_scalar = ncol(Z); intercept_Z = cbind(intercept_Z, Z)
  }else{
    p_scalar = 0;
  }

  # Model: y = mu (intercept) + Z %*% alpha
  # + X_feature[[1]] %*% B[[1]] + X_feature[[2]] %*% B[[p]] + ... + X_feature[[p]] %*% B[[p]] + epsilon,
  #  each element in epsilon ~ iid N(0,sigma^2)

  # Create the "design matrix", which should be a column bind of (1, Z (likely to be NULL), X_feature[[1]], ..., X_feature[[p]] )
  # Start with intercept column (mu)
  design_matrix = matrix(rep(1, n), ncol = 1)

  # Adding Z variables to design_matrix if it exists
  if(!is.null(Z)){
    design_matrix = cbind(design_matrix, Z)
  }

  # Adding X_features to design_matrix
  for(i in 1:p){
    design_matrix = cbind(design_matrix, X_feature[[i]])
  }

  # Obtain the OLS estimator of the regression coefficients
  diag_val = 1e-7
  success = FALSE
  while(!success) {
    tryCatch({
      coeff.OLS = crossprod(t(solve(crossprod(design_matrix) + diag(diag_val, dim(crossprod(design_matrix))[1]))), crossprod(design_matrix, y))
      success = TRUE
    },
    error = function(e) {
      message("Error in OLS computation, increasing diagonal value and trying again")
      diag_val <<- diag_val * 10; print(diag_val)
    })
  }

  # Initialize the MCMC (Gibbs sampler)
  cat("Initializing the MCMC...\n")

  # Store the MCMC
  if(boundary_reg == TRUE){
    posterior.draws = list(
      B.post = array(NA, dim = c(p, S + burnin, K))[ , , , drop = FALSE],
      Alpha.post = array(NA, dim = c(S + burnin, p_scalar + 1)),
      sigma.2.post = array(NA, dim = c(S + burnin, 1)),
      lambda.sq.post = array(NA, dim = c(p, S + burnin, K)),
      mu.post = array(NA, dim = c(S + burnin, p)),
      phi.post = array(NA, dim = c(S + burnin, p))
    )
    if (!is.null(Z)) {
      posterior.draws$sigma.2.j.post = array(NA, dim = c(S + burnin, p_scalar))
    }
    D.matrix = rbind(c(1, rep(0, ncol(D_d) - 1)), D_d, c(rep(0, ncol(D_d) - 1), 1))
  } else {
    posterior.draws = list(
      B.post = array(NA, dim = c(p, S + burnin, K))[ , , , drop = FALSE],
      Alpha.post = array(NA, dim = c(S + burnin, p_scalar + 1)),
      sigma.2.post = array(NA, dim = c(S + burnin, 1)),
      lambda.sq.post = array(NA, dim = c(p, S + burnin, K - 2)),
      mu.post = array(NA, dim = c(S + burnin, p)),
      phi.post = array(NA, dim = c(S + burnin, p))
    )
    if (!is.null(Z)) {
      posterior.draws$sigma.2.j.post = array(NA, dim = c(S + burnin, p_scalar))
    }
  }


  # Include t.grid in the posterior.draws
  posterior.draws$t.grid = t.grid

  # initialize the regression coefficients as OLS estimator
  coeff = coeff.OLS

  B.all = matrix(tail(coeff.OLS, p * K), nrow = p, byrow = TRUE)
  rownames(B.all) = paste0(rep('beta', p), 1:p)
  Alpha = head(coeff.OLS, p_scalar+1) # include the intercept term in the scalar regression coeff


  # initialize the log-variance/ log-vol of the second-differenced B-spline coeffs
  # evolParams: a list of length p
  evolParams = apply(B.all, 1, function(B)  initDHS(omega = as.matrix(tcrossprod(D_d, t(as.matrix(B))))))

  # h_k: a matrix of p_func col
  h_k = sapply(evolParams, function(evol_p) as.matrix(evol_p$ht))
  lambda_k.sq = exp(h_k)

  if(boundary_reg == TRUE){
    lambda0.sq =  1/rgamma(1,a.lambda, b.lambda)
    # lambda_k.sq.full: n.basis x p
    lambda_k.sq.full = apply(lambda_k.sq, 2, function(lambda_k.sq) c(lambda0.sq, lambda_k.sq, lambda0.sq))
  }

  mu.DHS = phi.DHS = rep(NA,p)

  # initialize the Alpha's variance parameter
  if (!is.null(Z)) {
    Sigma2_Alpha = rep(1, p_scalar)
  }

  # inform users that we start running the Gibbs sampler
  cat("Starting the Gibbs sampler...\n")

  # obtain the 10%,..., 90% percentile of S + burnin
  percentiles = round(quantile(1:(S + burnin), probs = seq(0.1, 0.9, by = 0.1)))



  # Gibbs sampler
  for (s in (1: (burnin + S))){

    # update sigma.2
    a.new = a0 + n/2
    b.new = b0 + sum((y-design_matrix%*%coeff)^2/2)
    sigma.2 = 1/rgamma(1, shape=a.new, rate = b.new)

    # update B-spline coeff for beta
    for (j in (1:p)) {
      if(boundary_reg == TRUE){
        Q.B.j = crossprod(X_feature[[j]])/sigma.2  + crossprod(D.matrix, 1/lambda_k.sq.full[,j] * D.matrix)
      }else{
        Q.B.j = crossprod(X_feature[[j]])/sigma.2 + crossprod(D_d, 1/lambda_k.sq[,j] * D_d)}

      # make a design matrix
      if(p == 1){
        l.B.j = crossprod(X_feature[[j]], (y - intercept_Z %*% Alpha))/sigma.2
      }else{
        design.X.func.j = cbind(intercept_Z,
                                get_func_pred(X_feature[-j]))
        coeff.func.j = c(Alpha, c(t(B.all[-j,])))
        l.B.j = crossprod(X_feature[[j]], (y-design.X.func.j %*% as.matrix(coeff.func.j)))/sigma.2}

      ch.Q.Bj = chol(Q.B.j)
      B.all[j, ] =  backsolve(ch.Q.Bj, forwardsolve(t(ch.Q.Bj), l.B.j) + rnorm(length(l.B.j)))
    }

    # Update mu and Alpha
    if(!is.null(Z)){
      Q.Alpha =  crossprod(intercept_Z)/sigma.2 + rbind(0,cbind(0,diag(1/Sigma2_Alpha, nrow=length(Sigma2_Alpha))))
    }else{
      Q.Alpha = crossprod(intercept_Z)/sigma.2
    }
    design.X.scalar = cbind(get_func_pred(X_feature))
    coeff.Alpha = c(t(B.all))
    l.Alpha = crossprod(intercept_Z, (y- design.X.scalar %*% coeff.Alpha))/sigma.2
    ch.Q.Alpha =chol(Q.Alpha)
    Alpha = backsolve(ch.Q.Alpha, forwardsolve(t(ch.Q.Alpha), l.Alpha) + rnorm(length(l.Alpha)))


    # update all coefficients
    coeff = as.matrix(c(Alpha, c(t(B.all))))

    # update the dynamic shrinkage process  {lambda_k.sq}
    Omega = apply(B.all, 1, function(B) tcrossprod(D_d, t(B))) # Matrix: (n.basis-2) x n_func
    for (j in 1:p) {
      sample.j = sampleDSP(omega = Omega[,j], evolParams = evolParams[[j]])
      h_k.j = sample.j$ht
      lambda_k.sq[,j] = as.numeric(exp(h_k.j))
      mu.DHS[j] = sample.j$dhs_mean; phi.DHS[j] = sample.j$dhs_phi
      evolParams[[j]] = list(
        ht = h_k.j,
        dhs_mean = sample.j$dhs_mean,
        dhs_phi = sample.j$dhs_phi,
        sigma_eta_t  = sample.j$sigma_eta_t,
        sigma_eta_0  = sample.j$sigma_eta_0,
        dhs_mean0 = sample.j$dhs_mean0
      )
    }

    # update Sigma2_Alpha
    if(!is.null(Z)){
      a.j = aj0 + 1/2
      b.j = bj0 + (Alpha[-1]^2)/2
      Sigma2_Alpha = 1/rgamma(p_scalar, shape=a.j, rate=b.j)
    }

    # update the boundary regularization parameters
    if(boundary_reg == TRUE){
      lambda0.sq = apply(B.all,1, function(B) 1/rgamma(1, a.lambda+1, b.lambda + (B[1]^2  + B[length(B)]^2)/2)) # vector of length p_func
      lambda_k.sq.full = rbind(lambda0.sq,
                               lambda_k.sq,
                               lambda0.sq)
    }


    # storage
    posterior.draws$B.post[,s,] = B.all
    posterior.draws$Alpha.post[s,] = Alpha
    posterior.draws$sigma.2.post[s] = sigma.2
    if(boundary_reg == TRUE){
      posterior.draws$lambda.sq.post[,s,] = t(lambda_k.sq.full)
    }else{
      posterior.draws$lambda.sq.post[,s,] = t(lambda_k.sq)
    }
    if(!is.null(Z)){
      posterior.draws$sigma.2.j.post[s,] = Sigma2_Alpha
    }
    posterior.draws$mu.post[s,] = mu.DHS
    posterior.draws$phi.post[s,] = phi.DHS


    # Inform the user about the progress
    if(s %in% percentiles){
      cat(paste("Gibbs sampler ", round(which(percentiles == s)*10, digits = 1), "% complete\n"))
      print(Sys.time())
      beta.post = list()
      for (j in 1:p) {
        beta.post[[j]] = tcrossprod(posterior.draws$B.post[j,(1:s),], basismat)
      }
      posterior.draws$beta.post = beta.post
      rm(beta.post)
      gc()
      if(save == TRUE){
        save(posterior.draws, file = save_file)
      }
    }
  }
  # inform the user that Gibbs sampler is done
  cat("Gibbs sampler is done!\n")

  # discard the burn-in period
  posterior.draws$B.post = posterior.draws$B.post[, ( (burnin+1): (S+burnin)), , drop = FALSE]
  posterior.draws$Alpha.post = tail(posterior.draws$Alpha.post, S)
  posterior.draws$sigma.2.post = tail(posterior.draws$sigma.2.post,S)
  posterior.draws$lambda.sq.post = posterior.draws$lambda.sq.post[,  ( (burnin+1): (S+burnin)),]
  if(!is.null(Z)){
    posterior.draws$sigma.2.j.post = tail(posterior.draws$sigma.2.j.post, S)
  }
  posterior.draws$mu.post = tail(posterior.draws$mu.post,S)
  posterior.draws$phi.post = tail(posterior.draws$phi.post,S)
  gc()

  # obtain posterior draws for intercept and Alpha:
  posterior.draws$intercept.post = posterior.draws$Alpha.post[,1]
  if(!is.null(Z)){
    posterior.draws$Alpha.post = posterior.draws$Alpha.post[,-1]
  }else{
    posterior.draws$Alpha.post = NULL
  }

  # obtain posterior draws for beta(t)
  posterior.draws$beta.post = array(NA, dim=c(p, S, dim(basismat)[1]))
  posterior.draws$beta.post.mean = array(NA, dim=c(p, dim(basismat)[1]))
  posterior.draws$beta.95CI = array(NA,  dim=c(p, dim(basismat)[1], 2))
  posterior.draws$beta.50CI = array(NA,  dim=c(p, dim(basismat)[1], 2))
  for (j in 1:p){
    posterior.draws$beta.post[j,,] = tcrossprod(posterior.draws$B.post[j,,] , basismat)
    posterior.draws$beta.post.mean[j, ] = apply(posterior.draws$beta.post[j,,], 2, mean)
    posterior.draws$beta.95CI[j,,] = t(apply(posterior.draws$beta.post[j,,], 2, function(x) quantile(x, c(0.025,0.975))))
    posterior.draws$beta.50CI[j,,] =  t(apply(posterior.draws$beta.post[j,,], 2, function(x) quantile(x, c(0.25,0.75))))
  }

  BASOFR.fit <-list()
  BASOFR.fit$beta.post <- posterior.draws$beta.post; posterior.draws$beta.post = NULL; gc()
  BASOFR.fit$intercept.post <- posterior.draws$intercept.post; posterior.draws$intercetp.post = NULL; gc()
  BASOFR.fit$Alpha.post <- posterior.draws$Alpha.post; posterior.draws$Alpha.post = NULL; gc()
  BASOFR.fit$sigma.2.post <- posterior.draws$sigma.2.post;  posterior.draws$sigma.2.post = NULL; gc()
  BASOFR.fit$beta.post.mean <- posterior.draws$beta.post.mean
  BASOFR.fit$beta.95CI <- posterior.draws$beta.95CI
  BASOFR.fit$beta.50CI <- posterior.draws$beta.50CI
  BASOFR.fit$t.grid <- posterior.draws$t.grid

  BASOFR.fit$other_parameters <- list(DHS.post = list(mu.post = posterior.draws$mu.post,
                                                      phi.post = posterior.draws$phi.post,
                                                      lambda.sq.post = posterior.draws$lambda.sq.post),
                                      B.post = posterior.draws$B.post)

  if(!is.null(posterior.draws$sigma.2.j.post)){
    BASOFR.fit$other_parameters$z.sigma.2.post = posterior.draws$sigma.2.j.post
  }

  rm(posterior.draws); gc()

  # if(save == TRUE){
  #   save(posterior.draws, file = save_file)
  # }

  if(save == TRUE){
    save(BASOFR.fit, file = save_file)
  }

  return(BASOFR.fit)

}

