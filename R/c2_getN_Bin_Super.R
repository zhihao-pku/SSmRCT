#' Sample size and power for binary endpoints
#'
#' Calculating required sample size when given power or power when given sample size for binary endpoints.
#'
#' @rdname getN_Bin
#'
#' @name getN_Bin
#'
#' @param p1 A vector. Rate of treatment group.
#' @param p0 A vector. Rate of control group.
#' @param cut A vector. Positive value for non-inferiority or equivalence margin. For RD, the margin is on the original scale. For RR and OR, the margin is on the log scale. For example, if the non-inferiority margins for RD, RR, OR are 0.2, 0.6, and 1.3, then the \code{cut = 0.2}, \code{cut = -log(0.6)}, and \code{cut = log(1.3)}, respectively.
#' @param alpha A vector. One-sided type I error rate. Default value is 0.025.
#' @param beta A vector. Type II error rate. When \code{beta} is \code{NA} and \code{N} is not \code{NA}, the power will be returned.
#' @param N A vector. Sample size. When \code{N} is \code{NA} and \code{beta} is not \code{NA}, the sample size will be returned.
#' @param r A vector. Ratio of sample sizes of treatment group to control group. Default value is 1.
#' @param scale A vector. Optional values are "RD" for rate difference, "RR" for relative risk, and "OR" for odds ratio. Default value is "RD".
#' @param direct If \code{direct = 1}, larger values of RD, RR, and OR are preferable. If \code{direct = -1}, smaller values of RD, RR, and OR are preferable.
#' @param maxN Maximum possible sample size (\code{N}) in equivalence design. Default value is 1e+06.
#'
#' @return A data frame containing input parameters and returned sample size or power.
#'
#' @details
#' Taking the larger RD, RR, and OR are preferable as an example. Sample size calculation is based on the following Z test, with the success criterion:
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
#' Where \eqn{\hat \delta = \hat p_1 - \hat p_0, \hat \sigma_0^2 = \hat p_0(1 - \hat p_0), \hat \sigma_1^2 = \hat p_1(1 - \hat p_1)} for RD; \eqn{\hat \delta = log(\frac{\hat p_1}{\hat p_0}), \hat \sigma_0^2 = \frac{(1 - \hat p_0)}{\hat p_0}, \hat \sigma_1^2 = \frac{(1 - \hat p_1)}{\hat p_1}} for RR; \eqn{\hat \delta = log(\frac{\hat p_1 / (1 - \hat p_1)}{\hat p_0 / (1 - \hat p_0)}), \hat \sigma_0^2 = \frac{1}{\hat p_0(1 - \hat p_0)}, \hat \sigma_1^2 = \frac{1}{\hat p_1(1 - \hat p_1)}} for OR; and \eqn{\Delta} is the non-inferiority or equivalence margin (\code{cut}).
#'
#' @export
#'
#' @examples
#' (v <- getN_Bin_Super(
#'   p1 = 0.6, p0 = 0.4, alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RD"
#' ))
#' getN_Bin_Super(
#'   p1 = 0.6, p0 = 0.4, alpha = 0.025, beta = NA, N = v$N, r = 1, scale = "RD"
#' )
#'
#' (v <- getN_Bin_Noninf(
#'   p1 = 0.6, p0 = 0.5, cut = log(1.4),
#'   alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RR", direct = -1
#' ))
#' getN_Bin_Noninf(
#'   p1 = 0.6, p0 = 0.5, cut = log(1.4),
#'   alpha = 0.025, beta = NA, N = v$N, r = 1, scale = "RR", direct = -1
#' )
getN_Bin_Super <- function(p1, p0, alpha = 0.025, beta = NA, N = NA, r = 1, scale = "RD") {
  eg <- as.data.frame(expand.grid(p1 = p1, p0 = p0, alpha = alpha, beta = beta, N = N, r = r, scale = scale, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    p1 <- R$p1
    p0 <- R$p0
    alpha <- R$alpha
    beta <- R$beta
    N <- R$N
    r <- R$r
    scale <- R$scale
    if (is.na(beta) & is.na(N)) {
      stop("Beta and N cannot be NA simultaneously.")
    }
    if (!is.na(beta) & (!is.na(N))) {
      stop("Set one of beta and N to NA.")
    }
    if (!scale %in% c("RD", "RR", "OR")) {
      stop("Parameter scale should be one of RD, RR, and OR")
    }
    getPwr <- function(n0, delta, sigma1, sigma0) {
      n1 <- r * n0
      z <- delta / sqrt(sigma1^2 / n1 + sigma0^2 / n0)
      dplyr::if_else(delta > 0, 1 - pnorm(q = qnorm(1 - alpha), mean = z, sd = 1), pnorm(q = -qnorm(1 - alpha), mean = z, sd = 1))
    }
    if (is.na(N) & (!is.na(beta))) {
      if (scale == "RD") {
        delta <- p1 - p0
        sigma1 <- sqrt(p1 * (1 - p1))
        sigma0 <- sqrt(p0 * (1 - p0))
      }
      if (scale == "RR") {
        delta <- log(p1 / p0)
        sigma1 <- sqrt((1 - p1) / p1)
        sigma0 <- sqrt((1 - p0) / p0)
      }
      if (scale == "OR") {
        delta <- log(p1 / (1 - p1) / p0 * (1 - p0))
        sigma1 <- sqrt(1 / (1 - p1) / p1)
        sigma0 <- sqrt(1 / (1 - p0) / p0)
      }
      n0 <- (qnorm(1 - alpha) + qnorm(1 - beta))^2 * (sigma1^2 / r + sigma0^2) / delta^2
      n0 <- ceiling(n0)
      n1 <- r * n0
      N <- n1 + n0
      pwr <- getPwr(n0, delta, sigma1, sigma0)
    }
    if (is.na(beta) & (!is.na(N))) {
      if (scale == "RD") {
        delta <- p1 - p0
        sigma1 <- sqrt(p1 * (1 - p1))
        sigma0 <- sqrt(p0 * (1 - p0))
      }
      if (scale == "RR") {
        delta <- log(p1 / p0)
        sigma1 <- sqrt((1 - p1) / p1)
        sigma0 <- sqrt((1 - p0) / p0)
      }
      if (scale == "OR") {
        delta <- log(p1 / (1 - p1) / p0 * (1 - p0))
        sigma1 <- sqrt(1 / (1 - p1) / p1)
        sigma0 <- sqrt(1 / (1 - p0) / p0)
      }
      n0 <- N / (r + 1)
      n1 <- r * n0
      pwr <- getPwr(n0, delta, sigma1, sigma0)
    }
    data.frame(p1, p0, scale, delta, alpha, beta, N, n1, n0, r, pwr)
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getN_Bin
#' @export
getN_Bin_Noninf <- function(p1, p0, cut, alpha = 0.025, beta = NA, N = NA, r = 1, scale = "RD", direct = 1) {
  eg <- as.data.frame(expand.grid(p1 = p1, p0 = p0, cut = cut, alpha = alpha, beta = beta, N = N, r = r, scale = scale, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    p1 <- R$p1
    p0 <- R$p0
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    N <- R$N
    r <- R$r
    scale <- R$scale
    if (cut < 0) {
      warning("Parameter cut should be a positive value.")
    }
    if (is.na(beta) & is.na(N)) {
      stop("Beta and N cannot be NA simultaneously.")
    }
    if (!is.na(beta) & (!is.na(N))) {
      stop("Set one of beta and N to NA.")
    }
    if (!scale %in% c("RD", "RR", "OR")) {
      stop("Parameter scale should be one of RD, RR, and OR")
    }
    if (!direct %in% c(-1, 1)) {
      stop("Parameter direct should be one of `1` or `-1`.")
    }
    getPwr <- function(n0, delta, sigma1, sigma0) {
      n1 <- r * n0
      z <- dplyr::if_else(direct == 1, (delta + cut) / sqrt(sigma1^2 / n1 + sigma0^2 / n0), (delta - cut) / sqrt(sigma1^2 / n1 + sigma0^2 / n0))
      dplyr::if_else(direct == 1, 1 - pnorm(q = qnorm(1 - alpha), mean = z, sd = 1), pnorm(q = -qnorm(1 - alpha), mean = z, sd = 1))
    }
    if (is.na(N) & (!is.na(beta))) {
      if (scale == "RD") {
        delta <- p1 - p0
        sigma1 <- sqrt(p1 * (1 - p1))
        sigma0 <- sqrt(p0 * (1 - p0))
      }
      if (scale == "RR") {
        delta <- log(p1 / p0)
        sigma1 <- sqrt((1 - p1) / p1)
        sigma0 <- sqrt((1 - p0) / p0)
      }
      if (scale == "OR") {
        delta <- log(p1 / (1 - p1) / p0 * (1 - p0))
        sigma1 <- sqrt(1 / (1 - p1) / p1)
        sigma0 <- sqrt(1 / (1 - p0) / p0)
      }
      n0 <- dplyr::if_else(direct == 1, (qnorm(1 - alpha) + qnorm(1 - beta))^2 * (sigma1^2 / r + sigma0^2) / (delta + cut)^2, (qnorm(1 - alpha) + qnorm(1 - beta))^2 * (sigma1^2 / r + sigma0^2) / (delta - cut)^2)
      n0 <- ceiling(n0)
      n1 <- r * n0
      N <- n1 + n0
      pwr <- getPwr(n0, delta, sigma1, sigma0)
    }
    if (is.na(beta) & (!is.na(N))) {
      if (scale == "RD") {
        delta <- p1 - p0
        sigma1 <- sqrt(p1 * (1 - p1))
        sigma0 <- sqrt(p0 * (1 - p0))
      }
      if (scale == "RR") {
        delta <- log(p1 / p0)
        sigma1 <- sqrt((1 - p1) / p1)
        sigma0 <- sqrt((1 - p0) / p0)
      }
      if (scale == "OR") {
        delta <- log(p1 / (1 - p1) / p0 * (1 - p0))
        sigma1 <- sqrt(1 / (1 - p1) / p1)
        sigma0 <- sqrt(1 / (1 - p0) / p0)
      }
      n0 <- N / (r + 1)
      n1 <- r * n0
      pwr <- getPwr(n0, delta, sigma1, sigma0)
    }
    data.frame(p1, p0, scale, delta, cut, alpha, beta, direct, N, n1, n0, r, pwr)
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}


#' @rdname getN_Bin
#' @export
getN_Bin_Equi <- function(p1, p0, cut, alpha = 0.025, beta = NA, N = NA, r = 1, scale = "RD", maxN = 1e+06) {
  eg <- as.data.frame(expand.grid(p1 = p1, p0 = p0, cut = cut, alpha = alpha, beta = beta, N = N, r = r, scale = scale, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    p1 <- R$p1
    p0 <- R$p0
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    N <- R$N
    r <- R$r
    scale <- R$scale
    if (cut < 0) {
      warning("Parameter cut should be a positive value.")
    }
    if (is.na(beta) & is.na(N)) {
      stop("Beta and N cannot be NA simultaneously.")
    }
    if (!is.na(beta) & (!is.na(N))) {
      stop("Set one of beta and N to NA.")
    }
    if (!scale %in% c("RD", "RR", "OR")) {
      stop("Parameter scale should be one of RD, RR, and OR")
    }
    getPwr <- function(n0, delta, sigma1, sigma0) {
      n1 <- r * n0
      z1 <- (delta + cut) / sqrt(sigma1^2 / n1 + sigma0^2 / n0)
      z2 <- (delta - cut) / sqrt(sigma1^2 / n1 + sigma0^2 / n0)
      mvtnorm::pmvnorm(lower = c(qnorm(1 - alpha), -Inf), upper = c(Inf, -qnorm(1 - alpha)), mean = c(z1, z2), sigma = matrix(1, nrow = 2, ncol = 2))
    }
    if (is.na(N) & (!is.na(beta))) {
      if (scale == "RD") {
        delta <- p1 - p0
        sigma1 <- sqrt(p1 * (1 - p1))
        sigma0 <- sqrt(p0 * (1 - p0))
      }
      if (scale == "RR") {
        delta <- log(p1 / p0)
        sigma1 <- sqrt((1 - p1) / p1)
        sigma0 <- sqrt((1 - p0) / p0)
      }
      if (scale == "OR") {
        delta <- log(p1 / (1 - p1) / p0 * (1 - p0))
        sigma1 <- sqrt(1 / (1 - p1) / p1)
        sigma0 <- sqrt(1 / (1 - p0) / p0)
      }
      n0 <- tryCatch(
        {
          uniroot(f = function(n0) {
            getPwr(n0, delta, sigma1, sigma0) - (1 - beta)
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
      pwr <- getPwr(n0, delta, sigma1, sigma0)
    }
    if (is.na(beta) & (!is.na(N))) {
      if (scale == "RD") {
        delta <- p1 - p0
        sigma1 <- sqrt(p1 * (1 - p1))
        sigma0 <- sqrt(p0 * (1 - p0))
      }
      if (scale == "RR") {
        delta <- log(p1 / p0)
        sigma1 <- sqrt((1 - p1) / p1)
        sigma0 <- sqrt((1 - p0) / p0)
      }
      if (scale == "OR") {
        delta <- log(p1 / (1 - p1) / p0 * (1 - p0))
        sigma1 <- sqrt(1 / (1 - p1) / p1)
        sigma0 <- sqrt(1 / (1 - p0) / p0)
      }
      n0 <- N / (r + 1)
      n1 <- r * n0
      pwr <- getPwr(n0, delta, sigma1, sigma0)
    }
    data.frame(p1, p0, scale, delta, cut, alpha, beta, N, n1, n0, r, pwr)
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}
