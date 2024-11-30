#' Extract sample weights
#'
#' A generic function to extract sample weights from objects of class 'balancing_weights'
#'
#' @param object An objects of class 'balancing_weights'.
#' @return Empirically normalized sample weights.
#' @examples
#' data(demo)
#' output1 <- balancing.weights(S, Z, X, method = "IC", naturalGroupProp)
#' get_weights(output1)
#' @export
#'
get_weights <- function(object) {
  if (!inherits(object, "balancing_weights")) {
    stop("Object is not of class 'balancing_weights'")
  }
  return(object$wt.v)
}


#' Extract percentage sample ESS
#'
#' A generic function to extract percentage sample ESS for different object classes
#'
#' @param object An objects of class 'balancing_weights' or 'causal_estimates'.
#' @return Percentage sample effective sample size (ESS) for the pseudo-population.
#'
#' @examples
#' data(demo)
#' output1 <- balancing.weights(S, Z, X, method = "IC", naturalGroupProp)
#' percentESS(output1)
#'
#' output2 <- causal.estimate(S, Z, X, Y, B = 5, method = "IC", naturalGroupProp)
#' percentESS(output2)
#'
#' @export
#'
# Generic function
percentESS <- function(object) {
  if (!inherits(object, "balancing_weights")) {
    if (!inherits(object, "causal_estimates")) {
      stop("Object is not of class 'causal_estimates'")
    }
  }

  return(object$percentESS)

}


#' Extract causal estimates (mean differences)
#'
#'
#' A generic function to extract mean differences for objects of class 'causal_estimates'.
#'
#' @param object An objects of class 'causal_estimates'.
#' @return The mean differences between two groups
#'
#' @examples
#' data(demo)
#' output2 <- causal.estimate(S, Z, X, Y, B = 5, method = "IC", naturalGroupProp)
#' mean_diff(output2)
#'
#' @export
mean_diff <- function(object) {
  if (!inherits(object, "causal_estimates")) {
    stop("Object is not of class 'causal_estimates'")
  }

  return(round(object$otherFeatures.v,2))
}

#' Extract sigma ratios
#'
#' A generic function to extract the ratios of standard deviations for objects of class 'causal_estimates'.
#'
#' @param object An objects of class 'causal_estimates'.
#' @return The ratios of standard deviations between two groups
#'
#' @examples
#' data(demo)
#' output2 <- causal.estimate(S, Z, X, Y, B = 5, method = "IC", naturalGroupProp)
#' mean_diff(output2)
#'
#' @export
#'
sigma_ratio <- function(object) {
  if (!inherits(object, "causal_estimates")) {
    stop("Object is not of class 'causal_estimates")
  }

  num_outcomes = length(object$otherFeatures.v)
  sigma_ratio_est = rep(NA, num_outcomes)
  for (i in 1:num_outcomes) {
    sigma_ratio_est[i] = object$moments.ar[,,i][2,1]/object$moments.ar[,,i][2,2]
  }

  return(round(sigma_ratio_est,2))
}
