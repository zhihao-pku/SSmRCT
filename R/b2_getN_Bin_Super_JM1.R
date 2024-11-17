#' Regional sample size allocation using Japan's Method 1 for binary endpoints
#'
#' Based on Japan's Method 1, given the global sample size and marginal probability (power) of efficacy consistency between target region and globally, calculate the required sample size allocated to the target region, in clinical trials using superiority, non-inferiority, and equivalence designs with binary endpoints.
#'
#' @rdname getN_Bin_JM1
#'
#' @name getN_Bin_JM1
#'
#' @param p1_a A vector. Rate of treatment group globally.
#' @param p0_a A vector. Rate of the control group globally.
#' @param p1_j A vector. Rate of treatment group in target region.
#' @param p0_j A vector. Rate of control group in target region.
#' @param pi A vector. Proportion of global efficacy to retain. Default value is 0.5, which means retaining half of the efficacy.
#' @param cut A vector. Positive value for non-inferiority or equivalence margin. For RD, the margin is on the original scale. For RR and OR, the margin is on the log scale. For example, if the non-inferiority margins for RD, RR, OR are 0.2, 0.6, and 1.3, then the \code{cut = 0.2}, \code{cut = -log(0.6)}, and \code{cut = log(1.3)}, respectively.
#' @param alpha A vector. One-sided type I error rate for global success, which is used to calculate global sample size only when \code{N} is \code{NA}. Default value is 0.025.
#' @param beta A vector. Type II error rate for global success, which is used to calculate global sample size only when \code{N} is missing.
#' @param beta1 A vector. Type II error rate for efficacy consistency between target region and globally. Default value is 0.2.
#' @param N A vector. Global sample size. When \code{N} is \code{NA} and \code{alpha} and \code{beta} are not \code{NA}, \code{N} will be calculated automatically.
#' @param r A vector. Ratio of sample sizes of treatment group to control group. Default value is 1.
#' @param scale A vector. Optional values are "RD" for rate difference, "RR" for relative risk, and "OR" for odds ratio. Default value is "RD".
#' @param direct If \code{direct = 1}, larger values of RD, RR, and OR are preferable. If \code{direct = -1}, smaller values of RD, RR, and OR are preferable.
#' @param maxN Maximum possible sample size (\code{N}) in equivalence design. Default value is 1e+06.
#'
#' @return A data frame where \code{f} is required proportion of sample size allocated to the target region, and \code{Nj} is required sample size for the target region, calculated as \code{Nj = N * f}.
#'
#' @details
#' The global success criterion and the efficacy consistency criterion between target region and globally could be found in \code{\link{getPwr_Bin_Super_JM1}}.
#'
#' @references
#' 1. Quan H, Li M, Chen J, et al. Assessment of Consistency of Treatment Effects in Multiregional Clinical Trials. Drug Information J. 2010;44(5):617-632. doi:10.1177/009286151004400509
#'
#' 2. Liao JJZ, Yu Z, Li Y. Sample size allocation in multiregional equivalence studies. Pharm Stat. 2018;17(5):570-577. doi:10.1002/pst.1871
#'
#' @seealso
#' \code{\link{getPwr_Bin_Super_JM1}}, \code{\link{getN_Bin_Super}}.
#'
#' @export
#'
#' @examples
#' getN_Bin_Super_JM1(
#'   p1_a = 0.75, p0_a = 0.5, p1_j = 0.7, p0_j = 0.5, pi = 0.5, beta1 = 0.2,
#'   N = 200, r = 1, scale = "RD"
#' )
#'
#' # Global sample size will be calculated based on alpha and beta.
#' getN_Bin_Noninf_JM1(
#'   p1_a = 0.5, p0_a = 0.5, p1_j = 0.6, p0_j = 0.5, pi = 0.5, cut = log(1.6),
#'   alpha = 0.025, beta = 0.2, beta1 = 0.2, N = NA, r = 1, scale = "RR",
#'   direct = -1
#' )
getN_Bin_Super_JM1 <- function(p1_a, p0_a, p1_j, p0_j, pi = 0.5, alpha = NA, beta = NA, beta1 = 0.2, N = NA, r = 1, scale = "RD") {
  eg <- as.data.frame(expand.grid(p1_a = p1_a, p0_a = p0_a, p1_j = p1_j, p0_j = p0_j, pi = pi, alpha = alpha, beta = beta, beta1 = beta1, N = N, r = r, scale = scale, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    p1_a <- R$p1_a
    p0_a <- R$p0_a
    p1_j <- R$p1_j
    p0_j <- R$p0_j
    pi <- R$pi
    alpha <- R$alpha
    beta <- R$beta
    beta1 <- R$beta1
    N <- R$N
    r <- R$r
    scale <- R$scale
    if (pi < 0 | pi > 1) {
      warning("Parameter pi generally is between 0 and 1.")
    }
    if ((is.na(alpha) | is.na(beta)) & is.na(N)) {
      stop("The combination of alpha and beta, and N, cannot both be NA.")
    }
    if ((!is.na(alpha) | (!is.na(beta))) & (!is.na(N))) {
      warning("When both alpha and beta are not NA, N will be calculated automatically.")
    }
    if (!scale %in% c("RD", "RR", "OR")) {
      stop("Parameter scale should be one of `RD`, `RR`, and `OR`")
    }
    if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getN_Bin_Super(p1 = p1_a, p0 = p0_a, alpha = alpha, beta = beta, N = NA, r = r, scale = scale)$N
    }
    getPwr <- function(f) {
      Nj <- N * f
      if (scale == "RD") {
        delta_a <- p1_a - p0_a
        delta_j <- p1_j - p0_j
        var_j <- p1_j * (1 - p1_j) / (r * Nj / (1 + r)) + p0_j * (1 - p0_j) / (Nj / (1 + r))
        var_a <- p1_a * (1 - p1_a) / (r * N / (1 + r)) + p0_a * (1 - p0_a) / (N / (1 + r))
      }
      if (scale == "RR") {
        delta_a <- log(p1_a / p0_a)
        delta_j <- log(p1_j / p0_j)
        var_j <- (1 - p1_j) / p1_j / (r * Nj / (1 + r)) + (1 - p0_j) / p0_j / (Nj / (1 + r))
        var_a <- (1 - p1_a) / p1_a / (r * N / (1 + r)) + (1 - p0_a) / p0_a / (N / (1 + r))
      }
      if (scale == "OR") {
        delta_a <- log(p1_a / (1 - p1_a) / p0_a * (1 - p0_a))
        delta_j <- log(p1_j / (1 - p1_j) / p0_j * (1 - p0_j))
        var_j <- 1 / (1 - p1_j) / p1_j / (r * Nj / (1 + r)) + 1 / (1 - p0_j) / p0_j / (Nj / (1 + r))
        var_a <- 1 / (1 - p1_a) / p1_a / (r * N / (1 + r)) + 1 / (1 - p0_a) / p0_a / (N / (1 + r))
      }
      sej <- sqrt(var_j + pi^2 * var_a - 2 * pi * sqrt(f) * sqrt(var_j * var_a))
      uj <- (delta_j - pi * delta_a) / sej
      uj <- dplyr::if_else(delta_a < 0, (-1) * uj, uj)
      mvtnorm::pmvnorm(lower = c(0), upper = c(Inf), mean = c(uj), sigma = 1)
    }
    f <- tryCatch(
      {
        uniroot(f = function(f) {
          getPwr(f) - (1 - beta1)
        }, interval = c(1e-06, 1 - 1e-06))$root
      },
      error = function(e) {
        warning("Cannot find an f value that meets the conditions.")
        NA
      }
    )
    p1_nj <- p0_nj <- pwr <- NA
    if (!is.na(f)) {
      p1_nj <- (p1_a - f * p1_j) / (1 - f)
      p0_nj <- (p0_a - f * p0_j) / (1 - f)
      pwr <- getPwr(f)
    }
    data.frame(p1_a, p0_a, p1_j, p0_j, p1_nj, p0_nj, pi, alpha, beta, N, r, scale, pwr, beta1, f, Nj = N * f) %>%
      dplyr::do(magrittr::set_rownames(., 1:nrow(.)))
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getN_Bin_JM1
#' @export
getN_Bin_Noninf_JM1 <- function(p1_a, p0_a, p1_j, p0_j, pi = 0.5, cut, alpha = NA, beta = NA, beta1 = 0.2, N = NA, r = 1, scale = "RD", direct = 1) {
  eg <- as.data.frame(expand.grid(p1_a = p1_a, p0_a = p0_a, p1_j = p1_j, p0_j = p0_j, pi = pi, cut = cut, alpha = alpha, beta = beta, beta1 = beta1, N = N, r = r, scale = scale, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    p1_a <- R$p1_a
    p0_a <- R$p0_a
    p1_j <- R$p1_j
    p0_j <- R$p0_j
    pi <- R$pi
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    beta1 <- R$beta1
    N <- R$N
    r <- R$r
    scale <- R$scale
    if (pi < 0 | pi > 1) {
      warning("Parameter pi generally is between 0 and 1.")
    }
    if (cut < 0) {
      warning("Parameter cut should be a positive value.")
    }
    if ((is.na(alpha) | is.na(beta)) & is.na(N)) {
      stop("The combination of alpha and beta, and N, cannot both be NA.")
    }
    if ((!is.na(alpha) | (!is.na(beta))) & (!is.na(N))) {
      warning("When both alpha and beta are not NA, N will be calculated automatically.")
    }
    if (!scale %in% c("RD", "RR", "OR")) {
      stop("Parameter scale should be one of `RD`, `RR`, and `OR`")
    }
    if (!direct %in% c(-1, 1)) {
      stop("Parameter direct should be one of `1` or `-1`.")
    }
    if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getN_Bin_Noninf(p1 = p1_a, p0 = p0_a, cut = cut, alpha = alpha, beta = beta, N = NA, r = r, scale = scale, direct = direct)$N
    }
    getPwr <- function(f) {
      Nj <- N * f
      if (scale == "RD") {
        delta_a <- p1_a - p0_a
        delta_j <- p1_j - p0_j
        var_j <- p1_j * (1 - p1_j) / (r * Nj / (1 + r)) + p0_j * (1 - p0_j) / (Nj / (1 + r))
        var_a <- p1_a * (1 - p1_a) / (r * N / (1 + r)) + p0_a * (1 - p0_a) / (N / (1 + r))
      }
      if (scale == "RR") {
        delta_a <- log(p1_a / p0_a)
        delta_j <- log(p1_j / p0_j)
        var_j <- (1 - p1_j) / p1_j / (r * Nj / (1 + r)) + (1 - p0_j) / p0_j / (Nj / (1 + r))
        var_a <- (1 - p1_a) / p1_a / (r * N / (1 + r)) + (1 - p0_a) / p0_a / (N / (1 + r))
      }
      if (scale == "OR") {
        delta_a <- log(p1_a / (1 - p1_a) / p0_a * (1 - p0_a))
        delta_j <- log(p1_j / (1 - p1_j) / p0_j * (1 - p0_j))
        var_j <- 1 / (1 - p1_j) / p1_j / (r * Nj / (1 + r)) + 1 / (1 - p0_j) / p0_j / (Nj / (1 + r))
        var_a <- 1 / (1 - p1_a) / p1_a / (r * N / (1 + r)) + 1 / (1 - p0_a) / p0_a / (N / (1 + r))
      }
      sej <- sqrt(var_j + var_a - 2 * sqrt(f) * sqrt(var_j * var_a))
      uj <- dplyr::if_else(direct == 1, (delta_j - delta_a + pi * cut) / sej, (delta_j - delta_a - pi * cut) / sej)
      uj <- dplyr::if_else(direct == -1, (-1) * uj, uj)
      mvtnorm::pmvnorm(lower = c(0), upper = c(Inf), mean = c(uj), sigma = 1)
    }
    f <- tryCatch(
      {
        uniroot(f = function(f) {
          getPwr(f) - (1 - beta1)
        }, interval = c(1e-06, 1 - 1e-06))$root
      },
      error = function(e) {
        warning("Cannot find an f value that meets the conditions.")
        NA
      }
    )
    p1_nj <- p0_nj <- pwr <- NA
    if (!is.na(f)) {
      p1_nj <- (p1_a - f * p1_j) / (1 - f)
      p0_nj <- (p0_a - f * p0_j) / (1 - f)
      pwr <- getPwr(f)
    }
    data.frame(p1_a, p0_a, p1_j, p0_j, p1_nj, p0_nj, pi, cut, alpha, beta, N, r, scale, direct, pwr, beta1, f, Nj = N * f) %>%
      dplyr::do(magrittr::set_rownames(., 1:nrow(.)))
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getN_Bin_JM1
#' @export
getN_Bin_Equi_JM1 <- function(p1_a, p0_a, p1_j, p0_j, pi = 0.5, cut, alpha = NA, beta = NA, beta1 = 0.2, N = NA, r = 1, scale = "RD", maxN = 1e+06) {
  eg <- as.data.frame(expand.grid(p1_a = p1_a, p0_a = p0_a, p1_j = p1_j, p0_j = p0_j, pi = pi, cut = cut, alpha = alpha, beta = beta, beta1 = beta1, N = N, r = r, scale = scale, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    p1_a <- R$p1_a
    p0_a <- R$p0_a
    p1_j <- R$p1_j
    p0_j <- R$p0_j
    pi <- R$pi
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    beta1 <- R$beta1
    N <- R$N
    r <- R$r
    scale <- R$scale
    if (pi < 0 | pi > 1) {
      warning("Parameter pi generally is between 0 and 1.")
    }
    if (cut < 0) {
      warning("Parameter cut should be a positive value.")
    }
    if ((is.na(alpha) | is.na(beta)) & is.na(N)) {
      stop("The combination of alpha and beta, and N, cannot both be NA.")
    }
    if ((!is.na(alpha) | (!is.na(beta))) & (!is.na(N))) {
      warning("When both alpha and beta are not NA, N will be calculated automatically.")
    }
    if (!scale %in% c("RD", "RR", "OR")) {
      stop("Parameter scale should be one of `RD`, `RR`, and `OR`")
    }
    if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getN_Bin_Equi(p1 = p1_a, p0 = p0_a, cut = cut, alpha = alpha, beta = beta, N = NA, r = r, scale = scale, maxN = maxN)$N
    }
    getPwr <- function(f) {
      Nj <- N * f
      if (scale == "RD") {
        delta_a <- p1_a - p0_a
        delta_j <- p1_j - p0_j
        var_j <- p1_j * (1 - p1_j) / (r * Nj / (1 + r)) + p0_j * (1 - p0_j) / (Nj / (1 + r))
        var_a <- p1_a * (1 - p1_a) / (r * N / (1 + r)) + p0_a * (1 - p0_a) / (N / (1 + r))
      }
      if (scale == "RR") {
        delta_a <- log(p1_a / p0_a)
        delta_j <- log(p1_j / p0_j)
        var_j <- (1 - p1_j) / p1_j / (r * Nj / (1 + r)) + (1 - p0_j) / p0_j / (Nj / (1 + r))
        var_a <- (1 - p1_a) / p1_a / (r * N / (1 + r)) + (1 - p0_a) / p0_a / (N / (1 + r))
      }
      if (scale == "OR") {
        delta_a <- log(p1_a / (1 - p1_a) / p0_a * (1 - p0_a))
        delta_j <- log(p1_j / (1 - p1_j) / p0_j * (1 - p0_j))
        var_j <- 1 / (1 - p1_j) / p1_j / (r * Nj / (1 + r)) + 1 / (1 - p0_j) / p0_j / (Nj / (1 + r))
        var_a <- 1 / (1 - p1_a) / p1_a / (r * N / (1 + r)) + 1 / (1 - p0_a) / p0_a / (N / (1 + r))
      }
      sej <- sqrt(var_j + var_a - 2 * sqrt(f) * sqrt(var_j * var_a))
      uj1 <- (delta_j - delta_a + pi * cut) / sej
      uj2 <- (delta_j - delta_a - pi * cut) / sej
      mvtnorm::pmvnorm(lower = c(0, -Inf), upper = c(Inf, 0), mean = c(uj1, uj2), sigma = matrix(1, nrow = 2, ncol = 2))
    }
    getPwr(0.1)
    f <- tryCatch(
      {
        uniroot(f = function(f) getPwr(f) - (1 - beta1), interval = c(1e-06, 1 - 1e-06))$root
      },
      error = function(e) {
        warning("Cannot find an f value that meets the conditions.")
        NA
      }
    )
    p1_nj <- p0_nj <- pwr <- NA
    if (!is.na(f)) {
      p1_nj <- (p1_a - f * p1_j) / (1 - f)
      p0_nj <- (p0_a - f * p0_j) / (1 - f)
      pwr <- getPwr(f)
    }
    data.frame(p1_a, p0_a, p1_j, p0_j, p1_nj, p0_nj, pi, cut, alpha, beta, N, r, scale, pwr, beta1, f, Nj = N * f) %>%
      dplyr::do(magrittr::set_rownames(., 1:nrow(.)))
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}
