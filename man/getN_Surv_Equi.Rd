% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/39 getN_Surv_Equi.R
\name{getN_Surv_Equi}
\alias{getN_Surv_Equi}
\title{Title}
\usage{
getN_Surv_Equi(delta, cut, alpha, beta, N, r, maxN = 1e+06)
}
\arguments{
\item{maxN}{}
}
\description{
Title
}
\examples{
(v <- getN_Surv_Equi(
  delta = log(1.1), cut = log(1.2), alpha = 0.025,
  beta = 0.2, N = NA, r = 1
))
getN_Surv_Equi(
  delta = log(1.1), cut = log(1.2), alpha = 0.025,
  beta = NA, N = v$N, r = 1
)
}
