% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{demo}
\alias{demo}
\alias{S}
\alias{X}
\alias{Y}
\alias{Z}
\alias{groupNames}
\alias{naturalGroupProp}
\title{Demo Dataset}
\format{
An rda object, with 450 observations and the following variables:
\describe{
\item{S}{A vector of factor levels, representing the study memberships.}
\item{Z}{A vector of factor levels, representing the group memberships.}
\item{X}{A covariate matrix.}
\item{Y}{An outcome matrix.}
\item{naturalGroupProp}{The relative group prevalences of the larger natural population. Necessary only for FLEXOR weights; it should be skipped for IC and IGO weights.}
\item{groupNames}{Disease subtype names "IDC" or "ILC"}
}
}
\usage{
data(demo)
}
\description{
A dataset containing example data for demonstration purposes.
}
\details{
Demo Dataset
}
\examples{
data(demo)
}
