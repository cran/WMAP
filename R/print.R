#' Print method for objects of class 'balancing_weights'
#'
#' @param x An object of class 'balancing_weights'
#' @param ... Additional arguments affecting the printed results (so far no additional arguments are needed, so leave blank).
#' @return Print values of the 'balancing_weights' object, including:
#'
#' \describe{
#'   \item{Weight length}{The total number of weights.}
#'   \item{percentESS}{Percentage sample effective sample size (ESS) for the pseudo-population.}
#' }
#'
#' @examples
#' data(demo)
#' output1 <- balancing.weights(S, Z, X, method = "IC", naturalGroupProp)
#' print(output1)
#'
#' @export
#'
#'
print.balancing_weights <- function(x, ...) {
  cat("An object of class 'balancing_weights'\n")
  cat("-------- Weight length --------\n")
  print(length(x$wt.v))

  cat("\n-------- Percentage sample ESS --------\n")
  print(x$percentESS)

  invisible(x)
}


#' Print method for objects of class 'causal_estimates'
#'
#' @param x An object of class 'causal_estimates'
#' @param ... Additional arguments affecting the printed results (so far no additional arguments are needed, so leave blank).
#' @return Print values of the 'causal_estimates' object, including:
#'
#' \describe{
#'   \item{Percentage sample ESS}{Percentage sample effective sample size (ESS) for the pseudo-population.}
#'   \item{Mean differences}{The mean differences between two groups}
#'   \item{Sigma ratios}{The ratios of standard deviations between two groups}
#' }
#'
#' @examples
#' data(demo)
#' output2 <- causal.estimate(S, Z, X, Y, B = 5, method = "IC", naturalGroupProp)
#' print(output2)
#'
#' @export
#'
#'
print.causal_estimates <- function(x, ...) {
  cat("An object of class 'causal_estimates'\n")

  cat("-------- Percentage sample ESS --------\n")
  print(x$percentESS)

  cat("\n-------- Mean differences --------\n")
  print(round(x$otherFeatures.v,2))

  cat("\n-------- Sigma ratios --------\n")

  num_outcomes = length(x$otherFeatures.v)
  sigma_ratio_est = rep(NA, num_outcomes)
  for (i in 1:num_outcomes) {
    sigma_ratio_est[i] = x$moments.ar[,,i][2,1]/x$moments.ar[,,i][2,2]
  }

  print(round(sigma_ratio_est,2))


  invisible(x)

}
