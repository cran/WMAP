#' Estimate causal effects using FLEXOR or other methods
#'
#' This function estimates causal effects based on the specified pseudo-population method.
#' The FLEXOR method involves an iterative two-step procedure.
#'
#' @param S Vector of factor levels representing the study memberships. Takes values in \{1, ..., J\}.
#' @param Z Vector of factor levels representing the group memberships. Takes values in \{1, ..., K\}.
#' @param X Covariate matrix of \eqn{N} rows and \eqn{p} columns.
#' @param Y Matrix of \eqn{L} outcomes, with dimensions \eqn{N \times L}.
#' @param B Number of bootstrap samples for variance estimation. Default is 100.
#' @param method Pseudo-population method, i.e., weighting method. Take values in \code{FLEXOR}, \code{IC}, or \code{IGO}.
#' @param seed Seed for random number generation. Default is \code{NULL}.
#' @param naturalGroupProp Relevant only for FLEXOR method: a fixed user-specified probability vector \eqn{\theta}.
#' @param num.random Relevant only for FLEXOR method: number of random starting points of \eqn{\gamma} in the two-step iterative procedure. Default is 40.
#' @param gammaMin Relevant only for FLEXOR method: Lower bound for each \eqn{\gamma_s} in the two-step iterative procedure. Default is 0.001.
#' @param gammaMax Relevant only for FLEXOR method: Upper bound for each \eqn{\gamma_s} in the two-step iterative procedure. Default is 0.999.
#' @param verbose Logical; if \code{TRUE} (default), displays progress messages during computation to the console. Set to \code{FALSE} to suppress these messages.
#'
#' @return An S3 list object with the following components:
#'
#' \describe{
#'   \item{percentESS}{Percentage sample effective sample size (ESS) of the pseudo-population.}
#'   \item{moments.ar}{An array of dimension \eqn{3 \times K \times L}, containing:
#'   \itemize{
#'     \item Estimated means, standard deviations (SDs), and medians (dimension 1),
#'     \item For \eqn{K} groups (dimension 2),
#'     \item And \eqn{L} counterfactual outcomes (dimension 3).
#'   }}
#'   \item{otherFeatures.v}{Estimated mean group differences for \eqn{L} outcomes.}
#'   \item{collatedMoments.ar}{An array of dimension \eqn{3 \times K \times L \times B}, containing:
#'   \itemize{
#'     \item \code{moments.ar} of the \eqn{b}th bootstrap sample (dimensions 1–3),
#'     \item For \eqn{B} bootstrap samples (dimension 4).
#'   }}
#'   \item{collatedOtherFeatures.mt}{A matrix of dimension \eqn{L \times B} containing:
#'   \itemize{
#'     \item \code{otherFeatures.v} of the \eqn{b}th bootstrap sample (dimension 1),
#'     \item For \eqn{B} bootstrap samples (dimension 2).
#'   }}
#'   \item{collatedESS}{A vector of length \eqn{B}} containing percentage sample ESS for \eqn{B} bootstrap samples.
#'   \item{method}{Pseudo-population method, i.e., weighting method.}
#' }
#'
#' @examples
#' data(demo)
#' set.seed(1)
#' causal.estimate(S, Z, X, Y, B = 5, method = "IC", naturalGroupProp)
#'
#' @export
causal.estimate <- function(S, Z, X, Y, B=100, method, naturalGroupProp=NULL, num.random=40, gammaMin=1e-3, gammaMax = (1-1e-3), seed=NULL, verbose = TRUE)
{

  #preamble_path <- function(){
  #  system.file("R", "preamble.R", package = "FLEXOR")
  #}
#
  #if (!("fn.FLEXOR" %in% objects())){
  #  #source("R/preamble.R", echo=FALSE)
  #  #source("preamble.R", echo=FALSE)
  #  source(preamble_path())
  #}

  Y = matrix(Y, nrow=nrow(Y))

  output1 = balancing.weights(S, Z, X, method, naturalGroupProp, num.random, gammaMin, gammaMax, seed, verbose)

  data = fn.createListData(S, Z, X, Y)

  fn.userInputChecks_3(data)

  flexorFlag = method %in% c("FLEXOR", "flexor")

  c(moments.ar, otherFeatures.v) %<-% fn.Stage2(local_data=data, local_wt.v=output1$wt.v)

  output2 = NULL

  if (B>0){

    c(collatedMoments.ar, collatedOtherFeatures.mt,collatedESS) %<-% bootStrap(data, moments.ar, otherFeatures.v, B, method=method, method_wt.v=output1$wt.v, naturalGroupProp, num.random.b = num.random, gammaMin, gammaMax, seed, verbose)

    output2$collatedMoments.ar = collatedMoments.ar
    output2$collatedOtherFeatures.mt = collatedOtherFeatures.mt

    output2$collatedESS = collatedESS
  }


  output2$percentESS = output1$percentESS
  output2$moments.ar = moments.ar
  output2$otherFeatures.v = otherFeatures.v
  output2$method = method

  class(output2) <- "causal_estimates"
  output2
}

