#' Sample size and power for continuous endpoints
#'
#' Calculating required sample size when given power or power when given sample size for continuous endpoints.
#'
#' @rdname getN_Con
#'
#' @name getN_Con
#'
#' @param delta A vector. Mean difference between treatment and control groups.
#' @param sigma A vector. Common standard deviation.
#' @param cut A vector. Positive value for non-inferiority or equivalence margin.
#' @param alpha A vector. One-sided type I error rate. Default value is 0.025.
#' @param beta A vector. Type II error rate. When \code{beta} is \code{NA} and \code{N} is not \code{NA}, the power will be returned.
#' @param N A vector. Sample size. When \code{N} is \code{NA} and \code{beta} is not \code{NA}, the sample size will be returned.
#' @param r A vector. Ratio of sample sizes of treatment group to control group. Default value is 1.
#' @param direct \code{direct = 1} indicates that a larger mean difference is preferable, while \code{direct = -1} indicates that a smaller mean difference is preferable.
#' @param maxN Maximum possible sample size (\code{N}) in equivalence design. Default value is 1e+06.
#'
#' @return A data frame containing input parameters and returned sample size or power.
#'
#' @details
#' Taking the larger mean difference is preferable as an example. Sample size calculation is based on the following Z test, with the success criterion:
#'
#' in superiority design:
#' \deqn{Z = \frac{\hat \delta}{\hat \sigma\sqrt{\frac{1}{N / (r + 1)} + \frac{1}{Nr / (r + 1)}}} > \Phi^{-1}(1 - \alpha)}
#'
#' in non-inferiority design:
#' \deqn{Z = \frac{\hat \delta + \Delta}{\hat \sigma\sqrt{\frac{1}{N / (r + 1)} + \frac{1}{Nr / (r + 1)}}} > \Phi^{-1}(1 - \alpha)}
#'
#' in equivalence design:
#' \deqn{Z_u = \frac{\hat \delta + \Delta}{\hat \sigma\sqrt{\frac{1}{N / (r + 1)} + \frac{1}{Nr / (r + 1)}}} > \Phi^{-1}(1 - \alpha)\text{ and } Z_l = \frac{\hat \delta - \Delta}{\hat \sigma\sqrt{\frac{1}{N / (r + 1)} + \frac{1}{Nr / (r + 1)}}} < \Phi^{-1}(\alpha)}
#'
#' Where \eqn{\hat \delta} is the mean difference between treatment and control groups, and \eqn{\Delta} is the non-inferiority or equivalence margin (\code{cut}).
#'
#' @export
#'
#' @examples
#' (v <- getN_Con_Super(delta = 1.5, sigma = 4, alpha = 0.025, beta = 0.2, N = NA, r = 1))
#' getN_Con_Super(delta = 1.5, sigma = 4, alpha = 0.025, beta = NA, N = v$N, r = 1)
#'
#' (v <- getN_Con_Noninf(
#'   delta = 1, sigma = 4, cut = 2, alpha = 0.025, beta = 0.2, N = NA, r = 1, direct = -1
#' ))
#' getN_Con_Noninf(
#'   delta = 1, sigma = 4, cut = 2, alpha = 0.025, beta = NA, N = v$N, r = 1, direct = -1
#' )
getN_Con_Super <- function(delta, sigma, alpha = 0.025, beta = NA, N = NA, r = 1) {
  eg <- as.data.frame(expand.grid(delta = delta, sigma = sigma, alpha = alpha, beta = beta, N = N, r = r, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta <- R$delta
    sigma <- R$sigma
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
    getPwr <- function(n0) {
      n1 <- r * n0
      z <- delta / (sigma * sqrt(1 / n1 + 1 / n0))
      dplyr::if_else(delta > 0, 1 - pnorm(q = qnorm(1 - alpha), mean = z, sd = 1), pnorm(q = -qnorm(1 - alpha), mean = z, sd = 1))
    }
    if (is.na(N) & (!is.na(beta))) {
      n0 <- (qnorm(1 - alpha) + qnorm(1 - beta))^2 * sigma^2 * (1 + 1 / r) / delta^2
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
    data.frame(delta, sigma, alpha, beta, N, n1, n0, r, pwr)
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getN_Con
#' @export
getN_Con_Noninf <- function(delta, sigma, cut, alpha = 0.025, beta = NA, N = NA, r = 1, direct = 1) {
  eg <- as.data.frame(expand.grid(delta = delta, sigma = sigma, cut = cut, alpha = alpha, beta = beta, N = N, r = r, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta <- R$delta
    sigma <- R$sigma
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
    getPwr <- function(n0) {
      n1 <- r * n0
      z <- dplyr::if_else(direct == 1, (delta + cut) / (sigma * sqrt(1 / n1 + 1 / n0)), (delta - cut) / (sigma * sqrt(1 / n1 + 1 / n0)))
      dplyr::if_else(direct == 1, 1 - pnorm(q = qnorm(1 - alpha), mean = z, sd = 1), pnorm(q = -qnorm(1 - alpha), mean = z, sd = 1))
    }
    if (is.na(N) & (!is.na(beta))) {
      n0 <- dplyr::if_else(direct == 1, (qnorm(1 - alpha) + qnorm(1 - beta))^2 * sigma^2 * (1 + 1 / r) / (delta + cut)^2, (qnorm(1 - alpha) + qnorm(1 - beta))^2 * sigma^2 * (1 + 1 / r) / (delta - cut)^2)
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
    data.frame(delta, sigma, cut, alpha, beta, direct, N, n1, n0, r, pwr)
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getN_Con
#' @export
getN_Con_Equi <- function(delta, sigma, cut, alpha = 0.025, beta = NA, N = NA, r = 1, maxN = 1e+06) {
  eg <- as.data.frame(expand.grid(delta = delta, sigma = sigma, cut = cut, alpha = alpha, beta = beta, N = N, r = r, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta <- R$delta
    sigma <- R$sigma
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
    getPwr <- function(n0) {
      n1 <- r * n0
      z1 <- (delta + cut) / (sigma * sqrt(1 / n1 + 1 / n0))
      z2 <- (delta - cut) / (sigma * sqrt(1 / n1 + 1 / n0))
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
    data.frame(delta, sigma, cut, alpha, beta, N, n1, n0, r, pwr)
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}
