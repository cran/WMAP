#' Compute balancing weights using FLEXOR or other methods
#'
#' This function calculates balancing weights based on the specified pseudo-population method.
#'
#' @param S Vector of factor levels representing the study memberships. Takes values in \{1, ..., J\}.
#' @param Z Vector of factor levels representing the group memberships. Takes values in \{1, ..., K\}.
#' @param X Covariate matrix of \eqn{N} rows and \eqn{p} columns.
#' @param method Pseudo-population method, i.e., weighting method. Take values in \code{FLEXOR}, \code{IC}, or \code{IGO}.
#' @param seed Seed for random number generation. Default is \code{NULL}.
#' @param naturalGroupProp Relevant only for FLEXOR method: a fixed user-specified probability vector \eqn{\theta}.
#' @param num.random Relevant only for FLEXOR method: number of random starting points of \eqn{\gamma} in the two-step iterative procedure. Default is 40.
#' @param gammaMin Relevant only for FLEXOR method: Lower bound for each \eqn{\gamma_s} in the two-step iterative procedure. Default is 0.001.
#' @param gammaMax Relevant only for FLEXOR method: Upper bound for each \eqn{\gamma_s} in the two-step iterative procedure. Default is 0.999.
#' @param verbose Logical; Relevant only for FLEXOR method: if \code{TRUE} (default), displays progress messages during computation to the console. Set to \code{FALSE} to suppress these messages.
#'
#' @return An S3 list object with the following components:
#'
#' \describe{
#'   \item{wt.v}{\eqn{N} empirically normalized sample weights.}
#'   \item{percentESS}{Percentage sample effective sample size (ESS) for the pseudo-population.}
#' }
#'
#' @importFrom zeallot %<-%
#' @examples
#' data(demo)
#' balancing.weights(S, Z, X, method = "IC", naturalGroupProp)
#'
#' @export
balancing.weights <- function(S, Z, X, method, naturalGroupProp, num.random=40, gammaMin=1e-3, gammaMax = (1-1e-3), seed=NULL, verbose = TRUE)
{

  flexorFlag = method %in% c("FLEXOR", "flexor")

  fn.userInputChecks_1(method, flexorFlag, gammaMin, gammaMax, num.random, naturalGroupProp)

  if (!is.null(seed)){
    if (!is.natural(seed))
      stop("seed must be a natural number.")

    set.seed(seed)
  }

  data = fn.createListData(S, Z, X)

  text = fn.userInputChecks_2(data)

  c(return.parm, percentESS) %<-% pkgcond::suppress_conditions(fn.Weights(data, method, flexorFlag, num.random, naturalGroupProp, gammaMin, gammaMax, verbose))

  output1 = NULL
  output1$wt.v = return.parm$wt.v
  output1$percentESS =  percentESS

  class(output1) <- "balancing_weights"
  return(output1)
}

