#' Combine vectors or matrices
#'
#' Combine two vectors to one vector or two matrices to one matrix
#'
#' @name combine
#' @param A A vector or matrix
#' @param B A vector or matrix
#'
#' @return A vector or matrix
#'
#' @export
#'
#' @examples
#' combine(c(1, 2, 3), c(4, 5, 6))
#'
#' combine(
#'   matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE),
#'   matrix(c(5, 6, 7, 8), nrow = 2, byrow = TRUE)
#' )
combine <- function(A, B) {
  if (class(A)[1] != class(B)[1] | (!is.vector(A) & !is.matrix(A)) | (!is.vector(B) & !is.matrix(B))) {
    stop("The data types of A and B should be consistent, either as vectors or data frames.")
  }
  if (is.matrix(A)) {
    dim <- nrow(A) * 2
    combined_matrix <- matrix(0, nrow = dim, ncol = dim)
    for (o in 1:nrow(A)) {
      combined_matrix[(2 * o - 1), seq(1, dim, by = 2)] <- A[o, ]
      combined_matrix[(2 * o), seq(1, dim, by = 2)] <- B[o, ]
      combined_matrix[(2 * o - 1), seq(2, dim, by = 2)] <- A[o, ]
      combined_matrix[(2 * o), seq(2, dim, by = 2)] <- B[o, ]
    }
    return(combined_matrix)
  }
  if (is.vector(A)) {
    combined_vector <- numeric(length(A) + length(B))
    combined_vector[seq(1, length(combined_vector), by = 2)] <- A
    combined_vector[seq(2, length(combined_vector), by = 2)] <- B
    return(combined_vector)
  }
}

#' Transform number of events to sample size
#'
#' Calculate required sample size according to required number of events, given hazard of event and follow-up parameters in survival analysis. Assuming uniform enrollment of subject and the event time and dropout time follow an exponential distribution.
#'
#' @name Ne_to_N
#' @param Ne A vector. Number of events
#' @param r A vector. Ratio of sample sizes of treatment group to control group. Default value is 1.
#' @param lambda0 A vector. Hazard of control group.
#' @param lambda1 A vector. Hazard of treatment group.
#' @param dropoutRate A vector. Dropout rate within the time interval specified by \code{dropoutTime} parameter.
#' @param dropoutTime A vector. Time interval for dropout rate.
#' @param a A vector. Accrual time for trial, which is only used when \code{follow_up = 'until_end'}.
#' @param f A vector. Follow up time for trial, which is only used when \code{follow_up = 'until_end'}.
#' @param l A vector. Fixed follow up period for each subject, which is only used when \code{follow_up = 'fixed_period'}.
#' @param follow_up A vector. If \code{follow_up = 'until_end'}, subjects will be followed up until the end of trial. If \code{follow_up = 'fixed_period'}, each subject will be followed up a fixed period.
#'
#' @return A data frame containing input parameters and returned event rate and required sample size.
#' \itemize{
#'   \item{\code{eventRate0 }}{Event rate for control group.}
#'   \item{\code{eventRate1 }}{Event rate for treatment group.}
#'   \item{\code{eventRate }}{Event rate for trial.}
#'   \item{\code{N0 }}{Required sample size for control group.}
#'   \item{\code{N1 }}{Required sample size for treatment group.}
#'   \item{\code{N }}{Required sample size for trial.}
#'   }
#'
#' @references
#' 1. Quan H, Zhao PL, Zhang J, Roessner M, Aizawa K. Sample size considerations for Japanese patients in a multi-regional trial based on MHLW guidance. Pharm Stat. 2010;9(2):100-112. doi:10.1002/pst.380
#'
#' @export
#'
#' @examples
#' # Median survival time in control group is 20 months, HR = 0.75, and annual dropout rate is 5%.
#' # Accrual time is 18 months, and follow-up time is 18 months.
#' # Each subject is followed up until the end of trial.
#' Ne_to_N(
#'   Ne = 100, r = 1, lambda0 = log(2) / 20, lambda1 = log(2) / 20 * 0.75,
#'   dropoutRate = 0.05, dropoutTime = 12,
#'   a = 18, f = 18, follow_up = "until_end"
#' )
#'
#' # Each subject is followed for 18 months after enrollment.
#' Ne_to_N(
#'   Ne = 100, r = 1, lambda0 = log(2) / 20, lambda1 = log(2) / 20 * 0.75,
#'   dropoutRate = 0.05, dropoutTime = 12,
#'   l = 18, follow_up = "fixed_period"
#' )
Ne_to_N <- function(Ne = NA, r = 1, lambda0, lambda1, dropoutRate, dropoutTime = 1, a = NA, f = NA, l = NA, follow_up = "until_end") {
  eg <- as.data.frame(expand.grid(Ne = Ne, r = r, lambda0 = lambda0, lambda1 = lambda1, dropoutRate = dropoutRate, dropoutTime = dropoutTime, a = a, f = f, l = l, follow_up = follow_up))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    Ne <- R$Ne
    r <- R$r
    lambda0 <- R$lambda0
    lambda1 <- R$lambda1
    dropoutRate <- R$dropoutRate
    dropoutTime <- R$dropoutTime
    a <- R$a
    f <- R$f
    l <- R$l
    follow_up <- R$follow_up
    if (!follow_up %in% c("until_end", "fixed_period")) {
      stop("Parameter follow_up should be one of `until_end` and `fixed_period`.")
    }
    if (follow_up == "until_end" & (is.na(a) | is.na(f))) {
      stop("When patients are followed up until the end of study, Paramters a and f can not be NA.")
    }
    if (follow_up == "until_end" & (!is.na(l))) {
      warning("When patients are followed up until the end of study, Paramter l will be ignored.")
    }
    if (follow_up == "fixed_period" & (is.na(l))) {
      stop("When patients are followed up a fixed period, Paramter l can not be NA.")
    }
    if (follow_up == "fixed_period" & (!is.na(a) | !is.na(f))) {
      warning("When patients are followed up a fixed period, Paramter a and b will be ignored.")
    }
    dR <- 1 - exp(log(1 - dropoutRate) / dropoutTime)
    if (follow_up == "until_end") {
      eR <- c(lambda0, lambda1) / (c(lambda0, lambda1) + dR) * (a - exp(-(c(lambda0, lambda1) + dR) * (a + f)) / (c(lambda0, lambda1) + dR) * (exp((c(lambda0, lambda1) + dR) * a) - 1)) / a
      eventRate <- sum(eR * c(1 / (r + 1), r / (r + 1)))
    }
    if (follow_up == "fixed_period") {
      eR <- c(lambda0, lambda1) / (c(lambda0, lambda1) + dR) * (1 - exp(-(c(lambda0, lambda1) + dR) * l))
      eventRate <- sum(eR * c(1 / (r + 1), r / (r + 1)))
    }
    data.frame(Ne = Ne, r = r, lambda0 = lambda0, lambda1 = lambda1, dropoutRate = dropoutRate, dropoutTime = dropoutTime, a = a, f = f, l = l, follow_up = follow_up, eventRate0 = eR[1], eventRate1 = eR[2], eventRate = eventRate, N0 = Ne / eventRate / (1 + r), N1 = Ne / eventRate * r / (1 + r), N = Ne / eventRate)
  }, .options = furrr::furrr_options(seed = TRUE))
  res
}
