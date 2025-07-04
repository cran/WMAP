#' Summary method for objects of class 'balancing_weights'
#'
#' @param object An object of class 'balancing_weights'
#' @param ... Additional arguments affecting the summary produced (so far no additional arguments are needed, so leave blank).
#' @return Printed summary of the 'balancing_weights' object, including:
#'
#' \describe{
#'   \item{Weight length}{The total number of weights.}
#'   \item{Weight distribution}{Statistical summary of weight values.}
#'   \item{percentESS}{Percentage sample effective sample size (ESS) for the pseudo-population.}
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
  print(summary(output1$wt.v))

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
#'   \item{Mean differences with 95% CI}{The mean differences between two groups with their corresponding 95% confidence intervals.}
#'   \item{Sigma ratios with 95% CI}{The ratios of standard deviations between two groups with their corresponding 95% confidence intervals.}
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

  cat("\n-------- Mean differences with 95% CI --------\n")
  CI.f = round(t(apply(object$collatedOtherFeatures.mt, 1, function(x) quantile(x,probs = c(0.025,0.975)))),2)
  print(write_res(round(object$otherFeatures.v,2), CI.f))

  cat("\n-------- Sigma ratios with 95% CI --------\n")
  print(write_sigma_ratio(object))
}
