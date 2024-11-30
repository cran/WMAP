#' Summary method for objects of class 'balancing_weights'
#'
#' @param object An object of class 'balancing_weights'
#' @param ... Additional arguments affecting the summary produced (so far no additional arguments are needed, so leave blank).
#' @return Printed summary of the 'balancing_weights' object, including:
#'
#' \describe{
#'   \item{Weight length}{The total number of weights.}
#'   \item{Weight distribution}{Statistical summary of weight values.}
#'   \item{Percentage sample ESS}{Percentage sample effective sample size (ESS) for the pseudo-population.}
#' }
#'
#' @examples
#' data(demo)
#' output1 <- balancing.weights(S, Z, X, method = "IC", naturalGroupProp)
#' summary(output1)
#'
#' @export
summary.balancing_weights <- function(object,...) {
  cat("Summary of object of class 'balancing_weights':\n")

  cat("-------- Weight length --------\n")
  print(length(object$wt.v))

  cat("\n-------- Weight distribution --------\n")
  print(summary(object$wt.v))

  cat("\n-------- Percentage sample ESS --------\n")
  print(object$percentESS)
}

#' Summary method for objects of class 'causal_estimates'
#'
#' @param object An object of class 'causal_estimates'
#' @param ... Additional arguments affecting the summary produced (so far no additional arguments are needed, so leave blank).
#'
#' @return Printed summary of the 'causal_estimates' object, including:
#'
#' \describe{
#'   \item{Percentage sample ESS}{Percentage sample effective sample size (ESS) for the pseudo-population.}
#'   \item{Mean differences with 95% CI (if \code{B > 0})}{The mean differences between two groups with their corresponding 95% confidence intervals.}
#'   \item{Sigma ratios with 95% CI (if \code{B > 0})}{The ratios of standard deviations between two groups with their corresponding 95% confidence intervals.}
#'   \item{Mean differences (if \code{B = 0})}{The mean differences between two groups.}
#'   \item{Sigma ratios (if \code{B = 0})}{The ratios of standard deviations between two groups.}
#' }
#'
#' @examples
#' data(demo)
#' set.seed(1)
#' output2 <- causal.estimate(S, Z, X, Y, B = 5, method = "IC", naturalGroupProp)
#' summary(output2)
#'
#'
#' @export
summary.causal_estimates <- function(object,...) {
  cat("Summary of object of class 'causal_estimates':\n")

  cat("-------- Percentage sample ESS --------\n")
  print(object$percentESS)

  if(!is.null(object$collatedESS)){
    cat("\n-------- Mean differences with 95% CI --------\n")
    CI.f = round(t(apply(object$collatedOtherFeatures.mt, 1, function(x) quantile(x,probs = c(0.025,0.975)))),2)
    print(write_res(round(object$otherFeatures.v,2), CI.f))

    cat("\n-------- Sigma ratios with 95% CI --------\n")
    print(write_sigma_ratio(object))
  } else if(is.null(object$collatedESS)){
    cat("\n-------- Mean differences --------\n")
    print(round(object$otherFeatures.v,2))

    cat("\n-------- Sigma ratios --------\n")
      num_outcomes = length(object$otherFeatures.v)
      sigma_ratio_est = rep(NA, num_outcomes)
      for (i in 1:num_outcomes) {
        sigma_ratio_est[i] = object$moments.ar[,,i][2,1]/object$moments.ar[,,i][2,2]
      }
    print(round(sigma_ratio_est,2))
  }


}
