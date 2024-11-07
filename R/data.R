#' Demo Dataset
#' @title Demo Dataset
#' @description A dataset containing example data for demonstration purposes.
#'
#' @docType data
#' @name demo
#' @aliases S X Y Z groupNames naturalGroupProp
#' @usage data(demo)
#' @format An rda object, with 450 observations and the following variables:
#' \describe{
#'   \item{S}{A vector of factor levels, representing the study memberships.}
#'   \item{Z}{A vector of factor levels, representing the group memberships.}
#'   \item{X}{A covariate matrix.}
#'   \item{Y}{An outcome matrix.}
#'   \item{naturalGroupProp}{The relative group prevalences of the larger natural population. Necessary only for FLEXOR weights; it should be skipped for IC and IGO weights.}
#'   \item{groupNames}{Disease subtype names "IDC" or "ILC"}
#' }
#' @examples
#' data(demo)
data(demo, envir = environment())





