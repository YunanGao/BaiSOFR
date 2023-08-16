#' Internal function to process functional features
#'
#' This function processes the functional features by combining them into one matrix.
#' It's intended for internal use within the package only.
#'
#' @param X.func.feature A list of functional features. Each element of the list is expected
#' to be a matrix or data frame representing one functional feature.
#'
#' @return A matrix combining all the functional features.
#'
#' @keywords internal
get_func_pred <- function(X.func.feature){
  p_func = length(X.func.feature)
  for (i in 1:p_func){
    if(i  == 1){
      X.func = X.func.feature[[i]]
    }else{
      X.func = cbind(X.func, X.func.feature[[i]])
    }
  }
  return(X.func)
}
