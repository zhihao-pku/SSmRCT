% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/58 getN_Surv_Noninf_JM1.R
\name{getN_Surv_Noninf_JM1}
\alias{getN_Surv_Noninf_JM1}
\title{Title}
\usage{
getN_Surv_Noninf_JM1(delta_j, delta_nj, pi, cut, beta1, N, r, direct = -1)
}
\arguments{
\item{direct}{}
}
\description{
Title
}
\examples{
getN_Surv_Noninf_JM1(
  delta_j = log(1.1), delta_nj = log(1.0),
  pi = 0.5, cut = log(1.3), beta1 = 0.2, N = seq(400, 800, 200),
  r = 1, direct = -1
)
}
