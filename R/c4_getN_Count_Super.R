#' Sample size and power for count endpoints
#'
#' Calculating required sample size when given power or power when given sample size for count endpoints.
#'
#' @rdname getN_Count
#'
#' @name getN_Count
#'
#' @param delta A vector. log(RR) between treatment and control groups.
#' @param lambda0 A vector. Baseline hazard of control group.
#' @param t A vector. Average exposure time.
#' @param k A vector. The over-dispersion parameter (k > 0) for negative binomial distribution, which is 0 for poisson distribution. Default value is 0.
#' @param cut A vector. Positive value for non-inferiority or equivalence margin. For example, if the non-inferiority margin for RR is 0.6, then \code{cut = -log(0.6)}. If the non-inferiority margin for RR is 1.3, then \code{cut = log(1.3)}.
#' @param alpha A vector. One-sided type I error rate. Default value is 0.025.
#' @param beta A vector. Type II error rate. When \code{beta} is \code{NA} and \code{N} is not \code{NA}, the power will be returned.
#' @param N A vector. Sample size. When \code{N} is \code{NA} and \code{beta} is not \code{NA}, the sample size will be returned.
#' @param r A vector. Ratio of sample sizes of treatment group to control group. Default value is 1.
#' @param direct \code{direct = 1} indicates that a larger RR is preferable, while \code{direct = -1} indicates that a smaller RR is preferable.
#' @param maxN Maximum possible sample size (\code{N}) in equivalence design. Default value is 1e+06.
#'
#' @return A data frame containing input parameters and returned sample size or power.
#'
#' @details
#' Taking the larger RR is preferable as an example. Sample size calculation is based on the following Z test, with the success criterion:
#'
#' in superiority design:
#' \deqn{Z = \frac{\hat \delta}{\sqrt{\frac{\hat \sigma_0^2}{N / (r + 1)} + \frac{\hat \sigma_1^2}{Nr / (r + 1)}}} > \Phi^{-1}(1 - \alpha)}
#'
#' in non-inferiority design:
#' \deqn{Z = \frac{\hat \delta + \Delta}{\sqrt{\frac{\hat \sigma_0^2}{N / (r + 1)} + \frac{\hat \sigma_1^2}{Nr / (r + 1)}}} > \Phi^{-1}(1 - \alpha)}
#'
#' in equivalence design:
#' \deqn{Z_u = \frac{\hat \delta + \Delta}{\sqrt{\frac{\hat \sigma_0^2}{N / (r + 1)} + \frac{\hat \sigma_1^2}{Nr / (r + 1)}}} > \Phi^{-1}(1 - \alpha)\text{ and } Z_l = \frac{\hat \delta - \Delta}{\sqrt{\frac{\hat \sigma_0^2}{N / (r + 1)} + \frac{\hat \sigma_1^2}{Nr / (r + 1)}}} < \Phi^{-1}(\alpha)}
#'
#' Where \eqn{\hat \delta = log(\hat{RR})} between treatment and control groups, \eqn{\hat \sigma_0^2 = \frac{1}{\hat \lambda_0 t} + \hat k, \sigma_1^2 = \frac{1}{e^{\hat \delta}\hat \lambda_0 t} + \hat k}, and \eqn{\Delta} is the non-inferiority or equivalence margin (\code{cut}).
#'
#' @export
#'
#' @examples
#' (v <- getN_Count_Super(
#'   delta = log(1.2),
#'   lambda0 = 0.5, t = 5, k = 0, alpha = 0.025, beta = 0.2, N = NA, r = 1
#' ))
#' getN_Count_Super(
#'   delta = log(1.2),
#'   lambda0 = 0.5, t = 5, k = 0, alpha = 0.025, beta = NA, N = v$N, r = 1
#' )
#'
#' (v <- getN_Count_Noninf(
#'   delta = log(1.1),
#'   lambda0 = 0.1, t = 5, k = 1, cut = log(1.4),
#'   alpha = 0.025, beta = 0.2, N = NA, r = 1, direct = -1
#' ))
#' getN_Count_Noninf(
#'   delta = log(1.1),
#'   lambda0 = 0.1, t = 5, k = 1, cut = log(1.4),
#'   alpha = 0.025, beta = NA, N = v$N, r = 1, direct = -1
#' )
getN_Count_Super <- function(delta, lambda0, t, k = 0, alpha = 0.025, beta = NA, N = NA, r = 1) {
  eg <- as.data.frame(expand.grid(delta = delta, lambda0 = lambda0, t = t, k = k, alpha = alpha, beta = beta, N = N, r = r, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta <- R$delta
    lambda0 <- R$lambda0
    t <- R$t
    k <- R$k
    alpha <- R$alpha
    beta <- R$beta
    N <- R$N
    r <- R$r
    if (is.na(beta) & is.na(N)) {
      stop("Beta and N cannot be NA simultaneously.")
    }
    if (!is.na(beta) & (!is.na(N))) {
      stop("Set one of beta and N to NA.")
    }
    lambda1 <- exp(delta) * lambda0
    sigma1 <- sqrt(1 / lambda1 / t + k)
    sigma0 <- sqrt(1 / lambda0 / t + k)
    getPwr <- function(n0) {
      n1 <- r * n0
      z <- delta / sqrt(sigma1^2 / n1 + sigma0^2 / n0)
      dplyr::if_else(delta > 0, 1 - pnorm(q = qnorm(1 - alpha), mean = z, sd = 1), pnorm(q = -qnorm(1 - alpha), mean = z, sd = 1))
    }
    if (is.na(N) & (!is.na(beta))) {
      n0 <- (qnorm(1 - alpha) + qnorm(1 - beta))^2 * (sigma1^2 + sigma0^2 / r) / delta^2
      n0 <- ceiling(n0)
      n1 <- r * n0
      N <- n1 + n0
      pwr <- getPwr(n0)
    }
    if (is.na(beta) & (!is.na(N))) {
      n0 <- N / (r + 1)
      n1 <- r * n0
      pwr <- getPwr(n0)
    }
    data.frame(delta, lambda0, t, k, alpha, beta, N, n1, n0, r, pwr)
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getN_Count
#' @export
getN_Count_Noninf <- function(delta, lambda0, t, k = 0, cut, alpha = 0.025, beta = NA, N = NA, r = 1, direct = 1) {
  eg <- as.data.frame(expand.grid(delta = delta, lambda0 = lambda0, t = t, k = k, cut = cut, alpha = alpha, beta = beta, N = N, r = r, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta <- R$delta
    lambda0 <- R$lambda0
    t <- R$t
    k <- R$k
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    N <- R$N
    r <- R$r
    if (cut < 0) {
      warning("Parameter cut should be a positive value.")
    }
    if (is.na(beta) & is.na(N)) {
      stop("Beta and N cannot be NA simultaneously.")
    }
    if (!is.na(beta) & (!is.na(N))) {
      stop("Set one of beta and N to NA.")
    }
    if (!direct %in% c(-1, 1)) {
      stop("Parameter direct should be one of `1` or `-1`.")
    }
    lambda1 <- exp(delta) * lambda0
    sigma1 <- sqrt(1 / lambda1 / t + k)
    sigma0 <- sqrt(1 / lambda0 / t + k)
    getPwr <- function(n0) {
      n1 <- r * n0
      z <- dplyr::if_else(direct == 1, (delta + cut) / sqrt(sigma1^2 / n1 + sigma0^2 / n0), (delta - cut) / sqrt(sigma1^2 / n1 + sigma0^2 / n0))
      dplyr::if_else(direct == 1, 1 - pnorm(q = qnorm(1 - alpha), mean = z, sd = 1), pnorm(q = -qnorm(1 - alpha), mean = z, sd = 1))
    }
    if (is.na(N) & (!is.na(beta))) {
      n0 <- dplyr::if_else(direct == 1, (qnorm(1 - alpha) + qnorm(1 - beta))^2 * (sigma1^2 + sigma0^2 / r) / (delta + cut)^2, (qnorm(1 - alpha) + qnorm(1 - beta))^2 * (sigma1^2 + sigma0^2 / r) / (delta - cut)^2)
      n0 <- ceiling(n0)
      n1 <- r * n0
      N <- n1 + n0
      pwr <- getPwr(n0)
    }
    if (is.na(beta) & (!is.na(N))) {
      n0 <- N / (r + 1)
      n1 <- r * n0
      pwr <- getPwr(n0)
    }
    data.frame(delta, lambda0, t, k, cut, alpha, beta, direct, N, n1, n0, r, pwr)
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getN_Count
#' @export
getN_Count_Equi <- function(delta, lambda0, t, k = 0, cut, alpha = 0.025, beta = NA, N = NA, r = 1, maxN = 1e+06) {
  eg <- as.data.frame(expand.grid(delta = delta, lambda0 = lambda0, t = t, k = k, cut = cut, alpha = alpha, beta = beta, N = N, r = r, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta <- R$delta
    lambda0 <- R$lambda0
    t <- R$t
    k <- R$k
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    N <- R$N
    r <- R$r
    if (cut < 0) {
      warning("Parameter cut should be a positive value.")
    }
    if (is.na(beta) & is.na(N)) {
      stop("Beta and N cannot be NA simultaneously.")
    }
    if (!is.na(beta) & (!is.na(N))) {
      stop("Set one of beta and N to NA.")
    }
    lambda1 <- exp(delta) * lambda0
    sigma1 <- sqrt(1 / lambda1 / t + k)
    sigma0 <- sqrt(1 / lambda0 / t + k)
    getPwr <- function(n0) {
      n1 <- r * n0
      z1 <- (delta + cut) / sqrt(sigma1^2 / n1 + sigma0^2 / n0)
      z2 <- (delta - cut) / sqrt(sigma1^2 / n1 + sigma0^2 / n0)
      mvtnorm::pmvnorm(lower = c(qnorm(1 - alpha), -Inf), upper = c(Inf, -qnorm(1 - alpha)), mean = c(z1, z2), sigma = matrix(1, nrow = 2, ncol = 2))
    }
    if (is.na(N) & (!is.na(beta))) {
      n0 <- tryCatch(
        {
          uniroot(f = function(n0) {
            getPwr(n0) - (1 - beta)
          }, interval = c(1e-06, maxN / (r + 1)))$root
        },
        error = function(e) {
          warning(sprintf("The calculated `N` exceeds %s", maxN))
          NA
        }
      )
      n0 <- ceiling(n0)
      n1 <- r * n0
      N <- n1 + n0
      pwr <- getPwr(n0)
    }
    if (is.na(beta) & (!is.na(N))) {
      n0 <- N / (r + 1)
      n1 <- r * n0
      pwr <- getPwr(n0)
    }
    data.frame(delta, lambda0, t, k, cut, alpha, beta, N, n1, n0, r, pwr)
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}
