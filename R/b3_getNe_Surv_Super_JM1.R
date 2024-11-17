#' Regional number of events allocation using Japan's Method 1 for survival endpoints
#'
#' Based on Japan's Method 1, given the global number of events and marginal probability (power) of efficacy consistency between target region and globally, calculate the required number of events allocated to the target region, in clinical trials using superiority, non-inferiority, and equivalence designs with survival endpoints.
#'
#' @rdname getNe_Surv_JM1
#'
#' @name getNe_Surv_JM1
#'
#' @param delta_a A vector. log(HR) between treatment and control groups globally.
#' @param delta_j A vector. log(HR) between treatment and control groups in target region.
#' @param pi A vector. Proportion of global efficacy to retain. Default value is 0.5, which means retaining half of the efficacy.
#' @param cut A vector. Positive value for non-inferiority or equivalence margin. For example, if the non-inferiority margin for HR is 0.6, then \code{cut = -log(0.6)}. If the non-inferiority margin for HR is 1.3, then \code{cut = log(1.3)}.
#' @param alpha A vector. One-sided type I error rate for global success, which is used to calculate global number of events only when \code{Ne} is \code{NA}. Default value is 0.025.
#' @param beta A vector. Type II error rate for global success, which is used to calculate global number of events only when \code{Ne} is \code{NA}.
#' @param beta1 A vector. Type II error rate for efficacy consistency between target region and globally. Default value is 0.2.
#' @param Ne A vector. Global number of events. When \code{Ne} is \code{NA} and \code{alpha} and \code{beta} are not \code{NA}, \code{Ne} will be calculated automatically.
#' @param r A vector. Ratio of the number of events of the treatment group to the control group. Default value is 1.
#' @param criterion A vector. If \code{criterion = 1}, the consistency criterion defined on the log(HR) scale will be used. If \code{criterion = 2}, the consistency criterion defined on the HR scale will be used. See \code{details} for more information.
#' @param direct \code{direct = 1} indicates that a larger HR is preferable, while \code{direct = -1} indicates that a smaller HR is preferable.
#' @param maxNe Maximum possible number of events (\code{Ne}) in equivalence design. Default value is 1e+06.
#'
#' @return A data frame where \code{f} is required proportion of number of events allocated to the target region, and \code{Nej} is required number of events for the target region, calculated as \code{Nej = Ne * f}.
#'
#' @details
#' The global success criterion and the efficacy consistency criterion between target region and globally could be found in \code{\link{getPwr_Surv_Super_JM1}}.
#'
#' @references
#' 1. Quan H, Li M, Chen J, et al. Assessment of Consistency of Treatment Effects in Multiregional Clinical Trials. Drug Information J. 2010;44(5):617-632. doi:10.1177/009286151004400509
#'
#' 2. Liao JJZ, Yu Z, Li Y. Sample size allocation in multiregional equivalence studies. Pharm Stat. 2018;17(5):570-577. doi:10.1002/pst.1871
#'
#' @seealso
#' \code{\link{getPwr_Surv_Super_JM1}}, \code{\link{getNe_Surv_Super}}.
#'
#' @export
#'
#' @examples
#' getNe_Surv_Super_JM1(
#'   delta_a = log(1.4),
#'   delta_j = log(1.3),
#'   pi = 0.5, beta1 = 0.2, Ne = 200, r = 1, criterion = 1
#' )
#'
#' # Global number of events will be calculated based on alpha and beta.
#' getNe_Surv_Noninf_JM1(
#'   delta_a = log(1.0),
#'   delta_j = log(1.1),
#'   pi = 0.5, cut = log(1.4),
#'   alpha = 0.025, beta = 0.2, beta1 = 0.2, Ne = NA, r = 1, criterion = 2,
#'   direct = -1
#' )
getNe_Surv_Super_JM1 <- function(delta_a, delta_j, pi = 0.5, alpha = NA, beta = NA, beta1 = 0.2, Ne = NA, r = 1, criterion = 1) {
  eg <- as.data.frame(expand.grid(delta_a = delta_a, delta_j = delta_j, pi = pi, alpha = alpha, beta = beta, beta1 = beta1, Ne = Ne, r = r, criterion = criterion, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta_a <- R$delta_a
    delta_j <- R$delta_j
    pi <- R$pi
    alpha <- R$alpha
    beta <- R$beta
    beta1 <- R$beta1
    Ne <- R$Ne
    r <- R$r
    criterion <- R$criterion
    if (pi < 0 | pi > 1) {
      warning("Parameter pi generally is between 0 and 1.")
    }
    if ((is.na(alpha) | is.na(beta)) & is.na(Ne)) {
      stop("The combination of alpha and beta, and Ne, cannot both be NA.")
    }
    if ((!is.na(alpha) | (!is.na(beta))) & (!is.na(Ne))) {
      warning("When both alpha and beta are not NA, Ne will be calculated automatically.")
    }
    if (!criterion %in% c(1, 2)) {
      stop("Parameter criterion should be one of `1` or `2`.")
    }
    gr <- 2 + r + 1 / r
    if (is.na(Ne) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getNe_Surv_Super(delta = delta_a, alpha = alpha, beta = beta, Ne = NA, r = r)$Ne
    }
    getPwr <- function(f) {
      Nej <- Ne * f
      if (criterion == 1) {
        sej <- sqrt(gr / Nej + pi^2 * gr / Ne - 2 * pi * sqrt(f) * sqrt(gr / Nej * gr / Ne))
        uj <- (delta_j - pi * delta_a) / sej
      }
      if (criterion == 2) {
        sej <- sqrt(gr / Nej * exp(delta_j)^2 + pi^2 * gr / Ne * exp(delta_a)^2 - 2 * pi * gr / Ne * exp(delta_j) * exp(delta_a))
        uj <- -(1 - exp(delta_j) - pi * (1 - exp(delta_a))) / sej
      }
      uj <- dplyr::if_else(delta_a < 0, (-1) * uj, uj)
      mvtnorm::pmvnorm(lower = c(0), upper = c(Inf), mean = c(uj), sigma = 1)
    }
    f <- tryCatch(
      {
        uniroot(f = function(f) getPwr(f) - (1 - beta1), interval = c(1e-06, 1 - 1e-06))$root
      },
      error = function(e) {
        warning("Cannot find an f value that meets the conditions.")
        NA
      }
    )
    delta_nj <- pwr <- NA
    if (!is.na(f)) {
      delta_nj <- (delta_a - f * delta_j) / (1 - f)
      pwr <- getPwr(f)
    }
    data.frame(delta_a, delta_j, delta_nj, pi, alpha, beta, Ne, r, criterion, pwr, beta1, f, Nej = Ne * f) %>%
      dplyr::do(magrittr::set_rownames(., 1:nrow(.)))
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getNe_Surv_JM1
#' @export
getNe_Surv_Noninf_JM1 <- function(delta_a, delta_j, pi = 0.5, cut, alpha = NA, beta = NA, beta1 = 0.2, Ne = NA, r = 1, criterion = 1, direct = 1) {
  eg <- as.data.frame(expand.grid(delta_a = delta_a, delta_j = delta_j, pi = pi, cut = cut, alpha = alpha, beta = beta, beta1 = beta1, Ne = Ne, r = r, criterion = criterion, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta_a <- R$delta_a
    delta_j <- R$delta_j
    pi <- R$pi
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    beta1 <- R$beta1
    Ne <- R$Ne
    r <- R$r
    criterion <- R$criterion
    if (pi < 0 | pi > 1) {
      warning("Parameter pi generally is between 0 and 1.")
    }
    if (cut < 0) {
      warning("Parameter cut should be a positive value.")
    }
    if ((is.na(alpha) | is.na(beta)) & is.na(Ne)) {
      stop("The combination of alpha and beta, and Ne, cannot both be NA.")
    }
    if ((!is.na(alpha) | (!is.na(beta))) & (!is.na(Ne))) {
      warning("When both alpha and beta are not NA, Ne will be calculated automatically.")
    }
    if (!criterion %in% c(1, 2)) {
      stop("Parameter criterion should be one of `1` or `2`.")
    }
    if (!direct %in% c(-1, 1)) {
      stop("Parameter direct should be one of `1` or `-1`.")
    }
    gr <- 2 + r + 1 / r
    if (is.na(Ne) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getNe_Surv_Noninf(delta = delta_a, cut = cut, alpha = alpha, beta = beta, Ne = NA, r = r, direct = direct)$Ne
    }
    getPwr <- function(f) {
      Nej <- Ne * f
      if (criterion == 1) {
        sej <- sqrt(gr / Nej + gr / Ne - 2 * sqrt(f) * sqrt(gr / Nej * gr / Ne))
        uj <- dplyr::if_else(direct == 1, (delta_j - delta_a + pi * cut) / sej, (delta_j - delta_a - pi * cut) / sej)
      }
      if (criterion == 2) {
        sej <- sqrt(gr / Nej * exp(delta_j)^2 + pi^2 * gr / Ne * exp(delta_a)^2 - 2 * pi * gr / Ne * exp(delta_j) * exp(delta_a))
        uj <- dplyr::if_else(direct == 1, -(1 / exp(cut) - exp(delta_j) - pi * (1 / exp(cut) - exp(delta_a))) / sej, -(exp(cut) - exp(delta_j) - pi * (exp(cut) - exp(delta_a))) / sej)
      }
      uj <- dplyr::if_else(direct == -1, (-1) * uj, uj)
      mvtnorm::pmvnorm(lower = c(0), upper = c(Inf), mean = c(uj), sigma = 1)
    }
    f <- tryCatch(
      {
        uniroot(f = function(f) getPwr(f) - (1 - beta1), interval = c(1e-06, 1 - 1e-06))$root
      },
      error = function(e) {
        warning("Cannot find an f value that meets the conditions.")
        NA
      }
    )
    delta_nj <- pwr <- NA
    if (!is.na(f)) {
      delta_nj <- (delta_a - f * delta_j) / (1 - f)
      pwr <- getPwr(f)
    }
    data.frame(delta_a, delta_j, delta_nj, pi, cut, alpha, beta, Ne, r, criterion, direct, pwr, beta1, f, Nej = Ne * f) %>%
      dplyr::do(magrittr::set_rownames(., 1:nrow(.)))
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getNe_Surv_JM1
#' @export
getNe_Surv_Equi_JM1 <- function(delta_a, delta_j, pi = 0.5, cut, alpha = NA, beta = NA, beta1 = 0.2, Ne = NA, r = 1, criterion = 1, maxNe = 1e+06) {
  eg <- as.data.frame(expand.grid(delta_a = delta_a, delta_j = delta_j, pi = pi, cut = cut, alpha = alpha, beta = beta, beta1 = beta1, Ne = Ne, r = r, criterion = criterion, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta_a <- R$delta_a
    delta_j <- R$delta_j
    pi <- R$pi
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    beta1 <- R$beta1
    Ne <- R$Ne
    r <- R$r
    criterion <- criterion
    if (pi < 0 | pi > 1) {
      warning("Parameter pi generally is between 0 and 1.")
    }
    if (cut < 0) {
      warning("Parameter cut should be a positive value.")
    }
    if ((is.na(alpha) | is.na(beta)) & is.na(Ne)) {
      stop("The combination of alpha and beta, and Ne, cannot both be NA.")
    }
    if ((!is.na(alpha) | (!is.na(beta))) & (!is.na(Ne))) {
      warning("When both alpha and beta are not NA, Ne will be calculated automatically.")
    }
    if (!criterion %in% c(1, 2)) {
      stop("Parameter criterion should be one of `1` or `2`.")
    }
    gr <- 2 + r + 1 / r
    if (is.na(Ne) & (!is.na(alpha)) & (!is.na(beta))) {
      Ne <- getNe_Surv_Equi(delta = delta_a, cut = cut, alpha = alpha, beta = beta, Ne = NA, r = r, maxNe = maxNe)$Ne
    }
    getPwr <- function(f) {
      Nej <- Ne * f
      if (criterion == 1) {
        sej <- sqrt(gr / Nej + gr / Ne - 2 * sqrt(f) * sqrt(gr / Nej * gr / Ne))
        uj1 <- (delta_j - delta_a + pi * cut) / sej
        uj2 <- (delta_j - delta_a - pi * cut) / sej
      }
      if (criterion == 2) {
        sej <- sqrt(gr / Nej * exp(delta_j)^2 + pi^2 * gr / Ne * exp(delta_a)^2 - 2 * pi * gr / Ne * exp(delta_j) * exp(delta_a))
        uj1 <- -(1 / exp(cut) - exp(delta_j) - pi * (1 / exp(cut) - exp(delta_a))) / sej
        uj2 <- -(exp(cut) - exp(delta_j) - pi * (exp(cut) - exp(delta_a))) / sej
      }
      mvtnorm::pmvnorm(lower = c(0, -Inf), upper = c(Inf, 0), mean = c(uj1, uj2), sigma = matrix(1, nrow = 2, ncol = 2))
    }
    f <- tryCatch(
      {
        uniroot(f = function(f) getPwr(f) - (1 - beta1), interval = c(1e-06, 1 - 1e-06))$root
      },
      error = function(e) {
        warning("Cannot find an f value that meets the conditions.")
        NA
      }
    )
    delta_nj <- pwr <- NA
    if (!is.na(f)) {
      delta_nj <- (delta_a - f * delta_j) / (1 - f)
      pwr <- getPwr(f)
    }
    data.frame(delta_a, delta_j, delta_nj, pi, cut, alpha, beta, Ne, r, criterion, pwr, beta1, f, Nej = Ne * f) %>%
      dplyr::do(magrittr::set_rownames(., 1:nrow(.)))
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}
