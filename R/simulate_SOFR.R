#' Simulate data for Scalar-on-Function Regression (SOFR) model
#'
#' This function simulates data for a SOFR model, which includes one or more functional covariates,
#' one or more control variables (both continuous and dummy), and a scalar response. The function
#' also calculates the signal-to-noise ratio for each functional covariate.
#'
#' @param n The number of observations (default: 1000)
#' @param p The number of functional covariates (default: 3)
#' @param SNR The signal-to-noise ratio (default: 7)
#' @param num_z_continuous The number of continuous control variables (default: 0)
#' @param num_z_dummy The number of dummy control variables (default: 0)
#' @param domain_X The maximal domain over which the functional covariates X may be observed (default: c(0,1)). The discrete points we get to observe functional X is given by t.grid = seq(domain_X[1], domain_X[2], length.out = max_obs)
#' @param max_obs The maximal number of discrete points to observe the functional covariates on domain_X (default: 201)
#' @param subject_specific_domain Whether to generate subject-specific domain (default: FALSE)
#' @param beta_types The shape of each beta (default: c('spiky', 'smooth', 'stepwise'))
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item t.grid: The time grid
#'   \item X: The functional covariates
#'   \item y: The response variable
#'   \item Z_continuous: The continuous control variables
#'   \item Z_dummy: The dummy control variables
#'   \item groundTruth: A list containing the true coefficients, signal-to-noise ratio and SNR for each functional covariate
#'   \item subject_specific_domain: The subject-specific domain vector if `subject_specific_domain` is TRUE
#' }
#'
#' @examples
#' # Simulate the data for a SOFR model with specific parameters
#' data = simulate_SOFR(n = 500, p = 2, SNR = 5, num_z_continuous = 1, num_z_dummy = 1,
#'                      domain_X = c(0,2), max_obs = 150, subject_specific_domain = TRUE,
#'                      beta_types = c('spiky', 'smooth'))
#'
#' # Check the output
#' str(data)
#'
#' @importFrom MASS mvrnorm
#' @export
simulate_SOFR <- function(n=1000, # the number of observations (xi, yi, zi(optional))
                          p=3, #the number of functional covariates,
                          SNR = 7, # the signal to noise ratio
                          num_z_continuous =  0, # the number of continuous control variables z
                          num_z_dummy = 0, # the number of dummy control variables z
                          domain_X = c(0,1), # the maximal domain over which the functional covariates X may be observed
                          max_obs = 201, # the maximal number of discrete points we get to observe the functional covariates on domain_X
                          subject_specific_domain = FALSE,
                          beta_types = c('spiky', 'smooth', 'stepwise') # a vector same length as p, indicating the shape of each beta.
){
  # the discrete points we get to observe functional X
  t.grid = seq(domain_X[1], domain_X[2], length.out = max_obs)
  # setup beta
  beta_list = list()
  for (j in 1:p) {
    if (beta_types[j] == 'spiky') {
      beta_list[[j]] = 20*(0.5 * (2+exp(20-60*t.grid) + exp(60*t.grid-20))^(-1) - 0.5*(2 + exp(40-60*t.grid) + exp(60*t.grid-40))^(-1))
    } else if (beta_types[j] == 'smooth') {
      beta_list[[j]] = 20*(0.15* exp(-20*(t.grid-0.5)^2) - 0.15*exp(-20*(t.grid-0.25)^2))
    } else if (beta_types[j] == 'stepwise') {
      beta_list[[j]] = 20*(-0.05 * (t.grid < 0.3) + 0.03 *(t.grid>=0.3 & t.grid< 0.7) - 0.04 * (t.grid>=0.7))
    } else {
      stop("Invalid beta type. Available options are 'spiky', 'smooth', 'stepwise'")
    }
  }

  # simulate X
  # use the FunctionalPredictors() function to simulate the X
  # it is a list of length p
  X = list()
  for(j in 1:p){
    X[[j]] = FunctionalPredictors(n=n, x=t.grid)
  }

  # if subject_specific_domain = TRUE, we generate the subject-specific domain
  if(subject_specific_domain == TRUE){
    if(max_obs > 200){
      # there are 10 unique subject-specific domains
      unique_subject_specific_domain = round(seq(round(max_obs * 0.7),   max_obs, length.out=10))
    }else{
      # there are 5 unique subject-specific domains
      unique_subject_specific_domain = round(seq(round(max_obs * 0.7),   max_obs, length.out=5))
    }

    # Additionally, ensure the max_obs is part of the unique_subject_specific_domain vector
    unique_subject_specific_domain[length(unique_subject_specific_domain)] <- max_obs

    # generate a length n vector, with each number in unique_subject_specific_domain appears about the same time:
    # Determine how many times each number should appear for approximate equal representation
    times_each_number = n %/% length(unique_subject_specific_domain)

    # Create a new vector by repeating each number 'times_each_number' times
    subject_specific_vector = rep(unique_subject_specific_domain, each = times_each_number)

    # If n is not a multiple of the length of unique_subject_specific_domain, there will be a remainder
    # Randomly sample the remaining elements from unique_subject_specific_domain
    remainder = n %% length(unique_subject_specific_domain)
    if (remainder > 0) {
      subject_specific_vector = c(subject_specific_vector, sample(unique_subject_specific_domain, size = remainder))
    }

    # Shuffle the resulting vector to ensure randomness
    subject_specific_vector = sample(subject_specific_vector)

    # Loop through each matrix in the list
    for(j in 1:p){
      # Get the current matrix
      current_matrix <- X[[j]]

      # Loop through each row in the current matrix
      for(i in 1:nrow(current_matrix)){
        # Replace values with NA after subject_specific_vector[i] in the i-th row
        # Check if subject_specific_vector[i] is less than ncol(current_matrix)
        if(subject_specific_vector[i] < ncol(current_matrix)) {
          current_matrix[i, (subject_specific_vector[i]+1):ncol(current_matrix)] <- NA
        }
      }

      # Replace the matrix in the list with the updated matrix
      X[[j]] <- current_matrix
    }
  }

  # simulate Z
  z_continuous = NULL
  z_dummy = NULL
  if(num_z_continuous != 0){
    A <- matrix(runif(num_z_continuous^2, min=0, max=0.3), ncol=num_z_continuous)
    Sigma <- t(A) %*% A + diag(1, nrow=num_z_continuous, ncol=num_z_continuous)
    z_continuous = MASS::mvrnorm(n =n, mu = rep(0,num_z_continuous), Sigma = Sigma)
  }

  if(num_z_dummy != 0){
    prob_z = runif(num_z_dummy) # generates random probabilities
    z_dummy = matrix(nrow = n, ncol = num_z_dummy)
    for(j in 1:num_z_dummy){
      z_dummy[,j] = rbinom(n, 1, prob = prob_z[j])
    }
  }

  # simulate alpha
  alpha = NULL
  if (num_z_continuous != 0 || num_z_dummy != 0) {
    p_Z = num_z_continuous + num_z_dummy
    alpha <- runif(p_Z, min=-1, max=1)
    # set ~20% of alpha to be zero
    alpha[sample(1:p_Z, size = floor(0.2 * p_Z))] <- 0
  }


  # simulate y
  Z = NULL
  if (num_z_continuous != 0 && num_z_dummy != 0) {
    Z = cbind(z_continuous, z_dummy)
  } else if (num_z_continuous != 0) {
    Z = z_continuous
  } else if (num_z_dummy != 0) {
    Z = z_dummy
  }

  if (is.null(alpha) || is.null(Z)) {
    y.hat <- rep(0, n)  # Create a vector of zeroes with length n
  } else {
    y.hat <- Z %*% alpha
  }

  grid_length = 1 / length(t.grid)
  if(subject_specific_domain != TRUE){
    for (j in 1:p) {
      y.hat = y.hat + X[[j]] %*% beta_list[[j]] * grid_length
    }
  }else{
    # Initialize an empty list to store the modified matrices
    X_0 = vector("list", length = p)

    # Loop through each matrix in X and replace NAs with zeros
    for (j in 1:p) {
      X_0[[j]] = replace(X[[j]], is.na(X[[j]]), 0)
    }

    # Now use the filled matrices for computation
    for (j in 1:p) {
      y.hat = y.hat + X_0[[j]] %*% beta_list[[j]] * grid_length
    }
  }

  noise.y = rnorm(n, 0, sd(y.hat) / SNR)
  y = y.hat + noise.y

  # calculate the signal-to-noise ratio for each X[[j]]
  SNR_X = rep(NA, p)
  if(subject_specific_domain != TRUE){
    for (j in 1:p) {
      SNR_X[j] = sd(X[[j]] %*% beta_list[[j]] * grid_length) / sd(y)
    }
  }else{
    for (j in 1:p) {
      SNR_X[j] = sd(X_0[[j]] %*% beta_list[[j]] * grid_length) / sd(y)
    }
  }


  # return
  data = list()
  data$t.grid = t.grid
  data$X = X
  data$y = y
  data$Z_continuous = z_continuous
  data$Z_dummy = z_dummy
  data$groundTruth = list(beta = beta_list,
                          alpha = alpha,
                          SNR = SNR,
                          SNR_X = SNR_X)
  if(subject_specific_domain == TRUE){
    data$subject_specific_domain = subject_specific_vector
  }

  return(data)
}
