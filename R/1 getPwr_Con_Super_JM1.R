#' Title
#'
#' @param delta_j
#' @param delta_nj
#' @param sigma
#' @param f
#' @param pi
#' @param alpha
#' @param beta
#' @param N
#' @param r
#' @param sim
#' @param nsim
#' @param seed
#' @param numcore
#'
#' @return
#' @export
#'
#' @examples
#' getPwr_Con_Super_JM1(
#'   delta_j = 0.5, delta_nj = 0.7, sigma = 1,
#'   f = seq(0.1, 0.9, 0.1),
#'   pi = 0.5, alpha = 0.025, beta = NA, N = 100, r = 1, sim = FALSE
#' )
#' getPwr_Con_Super_JM1(
#'   delta_j = 0.5, delta_nj = 0.7, sigma = 1,
#'   f = seq(0.1, 0.9, 0.1),
#'   pi = 0.5, alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = FALSE
#' )
getPwr_Con_Super_JM1 <- function(delta_j, delta_nj, sigma, f, pi = 0.5, alpha = 0.025, beta = NA, N, r = 1, sim = FALSE, nsim = 1000, seed = 0, numcore = 2) {
  if (!sim) {
    eg <- as.data.frame(expand.grid(
      delta_j = delta_j,
      delta_nj = delta_nj,
      sigma = sigma,
      f = f,
      pi = pi,
      alpha = alpha,
      beta = beta,
      N = N,
      r = r
    ))
    res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      sigma <- R$sigma
      f <- R$f
      pi <- R$pi
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      gr <- 2 + r + 1 / r
      delta <- delta_j * f + delta_nj * (1 - f)
      if (is.na(N)) {
        N <- getN_Con_Super(
          delta, sigma, alpha,
          beta, NA, r
        )$N
      }
      Nj <- N * f
      sej <- sqrt(gr * sigma^2 / Nj +
        pi^2 * gr * sigma^2 / N -
        2 * pi * sqrt(f) *
          sqrt(gr * sigma^2 / Nj * gr * sigma^2 / N))
      se <- sqrt(gr * sigma^2 / N)
      uj <- (delta_j - pi * delta) / sej
      u <- delta / se
      cov <- sqrt(f) * sqrt(gr * sigma^2 / N * gr * sigma^2 / Nj) -
        pi * gr * sigma^2 / N
      corr <- cov / sqrt(sej * se)
      M <- matrix(c(1, corr, corr, 1), nrow = 2, byrow = T)
      if (delta < 0) {
        uj <- (-1) * uj
        u <- (-1) * u
      }
      p1 <- pmvnorm(
        lower = c(qnorm(1 - alpha), -Inf),
        upper = c(Inf, Inf),
        mean = c(u, uj),
        corr = M
      )
      p2 <- pmvnorm(
        lower = c(-Inf, 0),
        upper = c(Inf, Inf),
        mean = c(u, uj),
        corr = M
      )
      p3 <- pmvnorm(
        lower = c(qnorm(1 - alpha), 0),
        upper = c(Inf, Inf),
        mean = c(u, uj),
        corr = M
      )
      p4 <- p3 / p1
      data.frame(
        delta_j, delta_nj, delta,
        sigma, f, pi,
        alpha, beta, N, r,
        p1, p2, p3, p4
      )
    }, .options = furrr_options(seed = TRUE))
  }
  if (sim) {
    eg <- as.data.frame(expand.grid(
      delta_j = delta_j,
      delta_nj = delta_nj,
      sigma = sigma,
      f = f,
      pi = pi,
      alpha = alpha,
      beta = beta,
      N = N,
      r = r,
      nsim = 1:nsim
    ))
    if (numcore >= 2) {
      plan(multisession, workers = numcore)
    }
    ss <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      sigma <- R$sigma
      f <- R$f
      pi <- R$pi
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      gr <- 2 + r + 1 / r
      delta <- delta_j * f + delta_nj * (1 - f)
      if (is.na(N)) {
        N <- getN_Con_Super(
          delta, sigma, alpha,
          beta, NA, r
        )$N
      }
      Nj <- N * f
      set.seed(i + seed)
      xt_j <- rnorm(n = Nj * r / (1 + r), mean = delta_j, sd = sigma)
      xc_j <- rnorm(n = Nj / (1 + r), mean = 0, sd = sigma)
      xt_nj <- rnorm(n = (N - Nj) * r / (1 + r), mean = delta_nj, sd = sigma)
      xc_nj <- rnorm(n = (N - Nj) / (1 + r), mean = 0, sd = sigma)
      xt <- c(xt_j, xt_nj)
      xc <- c(xc_j, xc_nj)
      za <- (mean(xt) - mean(xc)) /
        (sigma * sqrt(1 / length(xt) + 1 / length(xc)))
      zj <- mean(xt_j) - mean(xc_j) - pi * (mean(xt) - mean(xc))
      if (delta < 0) {
        za <- (-1) * za
        zj <- (-1) * zj
      }
      succ_a <- if_else(za > qnorm(1 - alpha), 1, 0)
      succ_j <- if_else(zj > 0, 1, 0)
      data.frame(
        delta_j, delta_nj, delta,
        sigma, f, pi, alpha, beta,
        N, r,
        succ_a, succ_j
      )
    }, .progress = TRUE, .options = furrr_options(seed = TRUE))
    if (numcore >= 2) {
      plan(sequential)
    }
    res <- ss %>%
      group_by(
        delta_j, delta_nj, delta, sigma, f, pi, alpha,
        beta, N, r
      ) %>%
      summarise(
        p1 = mean(succ_a),
        p2 = mean(succ_j),
        p3 = mean(succ_a & succ_j),
        p4 = mean(succ_j[succ_a == 1])
      ) %>%
      arrange(f) %>%
      as.data.frame()
  }
  return(res)
}
