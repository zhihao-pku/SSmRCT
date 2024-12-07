% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/b3_getNe_Surv_Super_JM1.R
\name{getNe_Surv_JM1}
\alias{getNe_Surv_JM1}
\alias{getNe_Surv_Super_JM1}
\alias{getNe_Surv_Noninf_JM1}
\alias{getNe_Surv_Equi_JM1}
\title{Regional number of events allocation using Japan's Method 1 for survival endpoints}
\usage{
getNe_Surv_Super_JM1(
  delta_a,
  delta_j,
  pi = 0.5,
  alpha = NA,
  beta = NA,
  beta1 = 0.2,
  Ne = NA,
  r = 1,
  criterion = 1
)

getNe_Surv_Noninf_JM1(
  delta_a,
  delta_j,
  pi = 0.5,
  cut,
  alpha = NA,
  beta = NA,
  beta1 = 0.2,
  Ne = NA,
  r = 1,
  criterion = 1,
  direct = 1
)

getNe_Surv_Equi_JM1(
  delta_a,
  delta_j,
  pi = 0.5,
  cut,
  alpha = NA,
  beta = NA,
  beta1 = 0.2,
  Ne = NA,
  r = 1,
  criterion = 1,
  maxNe = 1e+06
)
}
\arguments{
\item{delta_a}{A vector. log(HR) between treatment and control groups globally.}

\item{delta_j}{A vector. log(HR) between treatment and control groups in target region.}

\item{pi}{A vector. Proportion of global efficacy to retain. Default value is 0.5, which means retaining half of the efficacy.}

\item{alpha}{A vector. One-sided type I error rate for global success, which is used to calculate global number of events only when \code{Ne} is \code{NA}. Default value is 0.025.}

\item{beta}{A vector. Type II error rate for global success, which is used to calculate global number of events only when \code{Ne} is \code{NA}.}

\item{beta1}{A vector. Type II error rate for efficacy consistency between target region and globally. Default value is 0.2.}

\item{Ne}{A vector. Global number of events. When \code{Ne} is \code{NA} and \code{alpha} and \code{beta} are not \code{NA}, \code{Ne} will be calculated automatically.}

\item{r}{A vector. Ratio of the number of events of the treatment group to the control group. Default value is 1.}

\item{criterion}{A vector. If \code{criterion = 1}, the consistency criterion defined on the log(HR) scale will be used. If \code{criterion = 2}, the consistency criterion defined on the HR scale will be used. See \code{details} for more information.}

\item{cut}{A vector. Positive value for non-inferiority or equivalence margin. For example, if the non-inferiority margin for HR is 0.6, then \code{cut = -log(0.6)}. If the non-inferiority margin for HR is 1.3, then \code{cut = log(1.3)}.}

\item{direct}{\code{direct = 1} indicates that a larger HR is preferable, while \code{direct = -1} indicates that a smaller HR is preferable.}

\item{maxNe}{Maximum possible number of events (\code{Ne}) in equivalence design. Default value is 1e+06.}
}
\value{
A data frame where \code{f} is required proportion of number of events allocated to the target region, and \code{Nej} is required number of events for the target region, calculated as \code{Nej = Ne * f}.
}
\description{
Based on Japan's Method 1, given the global number of events and marginal probability (power) of efficacy consistency between target region and globally, calculate the required number of events allocated to the target region, in clinical trials using superiority, non-inferiority, and equivalence designs with survival endpoints.
}
\details{
The global success criterion and the efficacy consistency criterion between target region and globally could be found in \code{\link{getPwr_Surv_Super_JM1}}.
}
\examples{
getNe_Surv_Super_JM1(
  delta_a = log(1.4),
  delta_j = log(1.3),
  pi = 0.5, beta1 = 0.2, Ne = 200, r = 1, criterion = 1
)

# Global number of events will be calculated based on alpha and beta.
getNe_Surv_Noninf_JM1(
  delta_a = log(1.0),
  delta_j = log(1.1),
  pi = 0.5, cut = log(1.4),
  alpha = 0.025, beta = 0.2, beta1 = 0.2, Ne = NA, r = 1, criterion = 2,
  direct = -1
)
}
\references{
1. Quan H, Li M, Chen J, et al. Assessment of Consistency of Treatment Effects in Multiregional Clinical Trials. Drug Information J. 2010;44(5):617-632. doi:10.1177/009286151004400509

2. Liao JJZ, Yu Z, Li Y. Sample size allocation in multiregional equivalence studies. Pharm Stat. 2018;17(5):570-577. doi:10.1002/pst.1871
}
\seealso{
\code{\link{getPwr_Surv_Super_JM1}}, \code{\link{getNe_Surv_Super}}.
}
