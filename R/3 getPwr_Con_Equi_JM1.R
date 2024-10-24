#' Title
#'
#' @param delta_j
#' @param delta_nj
#' @param sigma
#' @param f
#' @param pi
#' @param cut
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
#' getPwr_Con_Equi_JM1(
#'   delta_j = -0.2, delta_nj = -0.1, sigma = 1,
#'   f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.4, alpha = 0.025, beta = NA,
#'   N = 400, r = 1, sim = FALSE
#' )
#' getPwr_Con_Equi_JM1(
#'   delta_j = -0.2, delta_nj = -0.1, sigma = 1,
#'   f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.4, alpha = 0.025, beta = 0.2,
#'   N = NA, r = 1, sim = FALSE
#' )
getPwr_Con_Equi_JM1 <- function(delta_j, delta_nj, sigma, f, pi = 0.5, cut, alpha = 0.025, beta = NA, N, r = 1, sim = FALSE, nsim = 1000, seed = 0, numcore = 2) {
  if (!sim) {
    eg <- as.data.frame(expand.grid(
      delta_j = delta_j,
      delta_nj = delta_nj,
      sigma = sigma,
      f = f,
      pi = pi,
      cut = cut,
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
      cut <- R$cut
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      gr <- 2 + r + 1 / r
      delta <- delta_j * f + delta_nj * (1 - f)
      if (is.na(N)) {
        N <- getN_Con_Equi(
          delta, sigma, cut, alpha,
          beta, NA, r
        )$N
      }
      Nj <- N * f
      se <- sqrt(gr * sigma^2 / N)
      u1 <- (delta + cut) / se
      u2 <- (delta - cut) / se
      sej <- sqrt(gr * sigma^2 / Nj +
        gr * sigma^2 / N -
        2 * sqrt(f) *
          sqrt(gr * sigma^2 / Nj * gr * sigma^2 / N))
      uj1 <- (delta_j - delta + pi * cut) / sej
      uj2 <- (delta_j - delta - pi * cut) / sej
      cov <- sqrt(f) * sqrt(gr * sigma^2 / N * gr * sigma^2 / Nj) -
        gr * sigma^2 / N
      corr <- cov / sqrt(sej * se)
      M <- matrix(c(
        1, 1, corr, corr,
        1, 1, corr, corr,
        corr, corr, 1, 1,
        corr, corr, 1, 1
      ), nrow = 4, byrow = T)
      p1 <- pmvnorm(
        lower = c(qnorm(1 - alpha), -Inf, -Inf, -Inf),
        upper = c(Inf, -qnorm(1 - alpha), Inf, Inf),
        mean = c(u1, u2, uj1, uj2),
        corr = M
      )
      p2 <- pmvnorm(
        lower = c(-Inf, -Inf, 0, -Inf),
        upper = c(Inf, Inf, Inf, 0),
        mean = c(u1, u2, uj1, uj2),
        corr = M
      )
      p3 <- pmvnorm(
        lower = c(qnorm(1 - alpha), -Inf, 0, -Inf),
        upper = c(Inf, -qnorm(1 - alpha), Inf, 0),
        mean = c(u1, u2, uj1, uj2),
        corr = M
      )
      p4 <- p3 / p1
      data.frame(
        delta_j, delta_nj, delta,
        sigma, f, pi, cut,
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
      cut = cut,
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
      cut <- R$cut
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      delta <- delta_j * f + delta_nj * (1 - f)
      if (is.na(N)) {
        N <- getN_Con_Equi(
          delta, sigma, cut, alpha,
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
      za1 <- (mean(xt) - mean(xc) + cut) /
        (sigma * sqrt(1 / length(xt) + 1 / length(xc)))
      za2 <- (mean(xt) - mean(xc) - cut) /
        (sigma * sqrt(1 / length(xt) + 1 / length(xc)))
      zj1 <- mean(xt_j) - mean(xc_j) - (mean(xt) - mean(xc)) + pi * cut
      zj2 <- mean(xt_j) - mean(xc_j) - (mean(xt) - mean(xc)) - pi * cut
      succ_a <- if_else(za1 > qnorm(1 - alpha) & za2 < -qnorm(1 - alpha), 1, 0)
      succ_j <- if_else(zj1 > 0 & zj2 < 0, 1, 0)
      data.frame(
        delta_j, delta_nj, delta,
        sigma, f, pi, cut,
        alpha, beta, N, r,
        succ_a, succ_j
      )
    }, .progress = TRUE, .options = furrr_options(seed = TRUE))
    if (numcore >= 2) {
      plan(sequential)
    }
    res <- ss %>%
      group_by(
        delta_j, delta_nj, delta, sigma,
        f, pi, cut, alpha, beta, N, r,
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
