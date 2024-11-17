#' Number of events and power for survival endpoints
#'
#' Calculating required number of events when given power or power when given number of events for survival endpoints.
#'
#' @rdname getNe_Surv
#'
#' @name getNe_Surv
#'
#' @param delta A vector. log(HR) between treatment and control groups.
#' @param alpha A vector. One-sided type I error rate. Default value is 0.025.
#' @param cut A vector. Positive value for non-inferiority or equivalence margin. For example, if the non-inferiority margin for HR is 0.6, then \code{cut = -log(0.6)}. If the non-inferiority margin for HR is 1.3, then \code{cut = log(1.3)}.
#' @param beta A vector. Type II error rate. When \code{beta} is \code{NA} and \code{Ne} is not \code{NA}, the power will be returned.
#' @param Ne A vector. Number of events. When \code{Ne} is \code{NA} and \code{beta} is not \code{NA}, the number of events will be returned.
#' @param r A vector. Ratio of number of events of the treatment group to the control group. Default value is 1 which under H0 assumption.
#' @param direct \code{direct = 1} indicates that a larger HR is preferable, while \code{direct = -1} indicates that a smaller HR is preferable.
#' @param maxNe Maximum possible number of events (\code{Ne}) in equivalence design. Default value is 1e+06.

#' @return A data frame containing input parameters and returned number of events or power.
#'
#' @details
#' Taking the larger HR is preferable as an example. Number of events calculation is based on the following Z test, with the success criterion:
#'
#' in superiority design:
#' \deqn{Z = \frac{\hat \delta}{\sqrt{\frac{1}{N_e / (r + 1)} + \frac{1}{N_e r / (r + 1)}}} > \Phi^{-1}(1 - \alpha)}
#'
#' in non-inferiority design:
#' \deqn{Z = \frac{\hat \delta + \Delta}{\sqrt{\frac{1}{N_e / (r + 1)} + \frac{1}{N_e r / (r + 1)}}} > \Phi^{-1}(1 - \alpha)}
#'
#' in equivalence design:
#' \deqn{Z_u = \frac{\hat \delta + \Delta}{\sqrt{\frac{1}{N_e / (r + 1)} + \frac{1}{N_e r / (r + 1)}}} > \Phi^{-1}(1 - \alpha) \text{ and } Z_l = \frac{\hat \delta - \Delta}{\sqrt{\frac{1}{N_e / (r + 1)} + \frac{1}{N_e r / (r + 1)}}} < \Phi^{-1}(\alpha)}
#'
#' Where \eqn{\hat \delta = log({\hat{HR}})} between treatment and control groups, and \eqn{\Delta} is the non-inferiority or equivalence margin (\code{cut}).
#'
#' @export
#'
#' @examples
#' (v <- getNe_Surv_Super(
#'   delta = log(1.2),
#'   alpha = 0.025, beta = 0.2, Ne = NA, r = 1
#' ))
#' getNe_Surv_Super(
#'   delta = log(1.2),
#'   alpha = 0.025, beta = NA, Ne = v$Ne, r = 1
#' )
#'
#' (v <- getNe_Surv_Noninf(
#'   delta = log(1.1),
#'   cut = log(1.3),
#'   alpha = 0.025, beta = 0.2, Ne = NA, r = 1, direct = -1
#' ))
#' getNe_Surv_Noninf(
#'   delta = log(1.1),
#'   cut = log(1.3),
#'   alpha = 0.025, beta = NA, Ne = v$Ne, r = 1, direct = -1
#' )
getNe_Surv_Super <- function(delta, alpha = 0.025, beta = NA, Ne = NA, r = 1) {
  eg <- as.data.frame(expand.grid(delta = delta, alpha = alpha, beta = beta, Ne = Ne, r = r, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta <- R$delta
    alpha <- R$alpha
    beta <- R$beta
    Ne <- R$Ne
    r <- R$r
    if (is.na(beta) & is.na(Ne)) {
      stop("Beta and Ne cannot be NA simultaneously.")
    }
    if (!is.na(beta) & (!is.na(Ne))) {
      stop("Set one of beta and Ne to NA.")
    }
    getPwr <- function(ne0) {
      ne1 <- r * ne0
      z <- delta / sqrt(1 / ne1 + 1 / ne0)
      dplyr::if_else(delta > 0, 1 - pnorm(q = qnorm(1 - alpha), mean = z, sd = 1), pnorm(q = -qnorm(1 - alpha), mean = z, sd = 1))
    }
    if (is.na(Ne) & (!is.na(beta))) {
      ne0 <- (qnorm(1 - alpha) + qnorm(1 - beta))^2 * (1 + 1 / r) / delta^2
      ne0 <- ceiling(ne0)
      ne1 <- r * ne0
      Ne <- ne1 + ne0
      pwr <- getPwr(ne0)
    }
    if (is.na(beta) & (!is.na(Ne))) {
      ne0 <- Ne / (r + 1)
      ne1 <- r * ne0
      pwr <- getPwr(ne0)
    }
    data.frame(delta, alpha, beta, Ne, ne1, ne0, r, pwr)
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getNe_Surv
#' @export
getNe_Surv_Noninf <- function(delta, cut, alpha = 0.025, beta = NA, Ne = NA, r = 1, direct = 1) {
  eg <- as.data.frame(expand.grid(delta = delta, cut = cut, alpha = alpha, beta = beta, Ne = Ne, r = r, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta <- R$delta
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    Ne <- R$Ne
    r <- R$r
    if (cut < 0) {
      warning("Parameter cut should be a positive value.")
    }
    if (is.na(beta) & is.na(Ne)) {
      stop("Beta and Ne cannot be NA simultaneously.")
    }
    if (!is.na(beta) & (!is.na(Ne))) {
      stop("Set one of beta and Ne to NA.")
    }
    if (!direct %in% c(-1, 1)) {
      stop("Parameter direct should be one of `1` or `-1`.")
    }
    getPwr <- function(ne0) {
      ne1 <- r * ne0
      z <- dplyr::if_else(direct == 1, (delta + cut) / sqrt(1 / ne1 + 1 / ne0), (delta - cut) / sqrt(1 / ne1 + 1 / ne0))
      dplyr::if_else(direct == 1, 1 - pnorm(q = qnorm(1 - alpha), mean = z, sd = 1), pnorm(q = -qnorm(1 - alpha), mean = z, sd = 1))
    }
    if (is.na(Ne) & (!is.na(beta))) {
      ne0 <- dplyr::if_else(direct == 1, (qnorm(1 - alpha) + qnorm(1 - beta))^2 * (1 + 1 / r) / (delta + cut)^2, (qnorm(1 - alpha) + qnorm(1 - beta))^2 * (1 + 1 / r) / (delta - cut)^2)
      ne0 <- ceiling(ne0)
      ne1 <- r * ne0
      Ne <- ne1 + ne0
      pwr <- getPwr(ne0)
    }
    if (is.na(beta)) {
      ne0 <- Ne / (r + 1)
      ne1 <- r * ne0
      pwr <- getPwr(ne0)
    }
    data.frame(delta, cut, alpha, beta, direct, Ne, ne1, ne0, r, pwr)
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getNe_Surv
#' @export
getNe_Surv_Equi <- function(delta, cut, alpha = 0.025, beta = NA, Ne = NA, r = 1, maxNe = 1e+06) {
  eg <- as.data.frame(expand.grid(delta = delta, cut = cut, alpha = alpha, beta = beta, Ne = Ne, r = r, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta <- R$delta
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    Ne <- R$Ne
    r <- R$r
    if (cut < 0) {
      warning("Parameter cut should be a positive value.")
    }
    if (is.na(beta) & is.na(Ne)) {
      stop("Beta and Ne cannot be NA simultaneously.")
    }
    if (!is.na(beta) & (!is.na(Ne))) {
      stop("Set one of beta and Ne to NA.")
    }
    getPwr <- function(ne0) {
      ne1 <- r * ne0
      z1 <- (delta + cut) / sqrt(1 / ne1 + 1 / ne0)
      z2 <- (delta - cut) / sqrt(1 / ne1 + 1 / ne0)
      mvtnorm::pmvnorm(lower = c(qnorm(1 - alpha), -Inf), upper = c(Inf, -qnorm(1 - alpha)), mean = c(z1, z2), sigma = matrix(1, nrow = 2, ncol = 2))
    }
    if (is.na(Ne) & (!is.na(beta))) {
      ne0 <- tryCatch(
        {
          uniroot(f = function(ne0) {
            getPwr(ne0) - (1 - beta)
          }, interval = c(1e-06, maxNe / (r + 1)))$root
        },
        error = function(e) {
          warning(sprintf("The calculated `Ne` exceeds %s", maxNe))
          NA
        }
      )
      ne0 <- ceiling(ne0)
      ne1 <- r * ne0
      Ne <- ne1 + ne0
      pwr <- getPwr(ne0)
    }
    if (is.na(beta) & (!is.na(Ne))) {
      ne0 <- Ne / (r + 1)
      ne1 <- r * ne0
      pwr <- getPwr(ne0)
    }
    data.frame(delta, cut, alpha, beta, Ne, ne1, ne0, r, pwr)
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}
