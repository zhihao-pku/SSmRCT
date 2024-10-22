#' Title
#'
#' @param delta_j
#' @param delta_nj
#' @param sigma
#' @param f
#' @param pi
#' @param cut
#' @param alpha
#' @param N
#' @param r
#' @param direct
#' @param sim
#' @param nsim
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
#' getPwr_Con_Noninf_JM1(
#'   delta_j = -0.2, delta_nj = -0.1, sigma = 1,
#'   f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.4, alpha = 0.025,
#'   N = 400, r = 1, direct = 1, sim = FALSE
#' )
getPwr_Con_Noninf_JM1 <- function(delta_j, delta_nj, sigma, f, pi, cut, alpha, N, r, direct = 1, sim = FALSE, nsim = 1000, seed = 0) {
  if (!sim) {
    eg <- as.data.frame(expand.grid(
      delta_j = delta_j,
      delta_nj = delta_nj,
      sigma = sigma,
      f = f,
      pi = pi,
      cut = cut,
      alpha = alpha,
      N = N,
      r = r,
      direct = direct
    ))
    res <- map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      sigma <- R$sigma
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      N <- R$N
      r <- R$r
      direct <- R$direct
      Nj <- N * f
      gr <- 2 + r + 1 / r
      delta <- delta_j * f + delta_nj * (1 - f)
      se <- sqrt(gr * sigma^2 / N)
      u <- if_else(direct == 1,
        (delta + cut) / se,
        (delta - cut) / se
      )
      sej <- sqrt(gr * sigma^2 / Nj +
        gr * sigma^2 / N -
        2 * sqrt(f) *
          sqrt(gr * sigma^2 / Nj * gr * sigma^2 / N))
      uj <- if_else(direct == 1,
        (delta_j - delta + pi * cut) / sej,
        (delta_j - delta - pi * cut) / sej
      )
      cov <- sqrt(f) * sqrt(gr * sigma^2 / N * gr * sigma^2 / Nj) -
        gr * sigma^2 / N
      corr <- cov / sqrt(sej * se)
      M <- matrix(c(1, corr, corr, 1), nrow = 2, byrow = T)
      if (direct == -1) {
        u <- (-1) * u
        uj <- (-1) * uj
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
        sigma, f, pi, cut,
        alpha, N, r, direct,
        p1, p2, p3, p4
      )
    })
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
      N = N,
      r = r,
      direct = direct,
      nsim = 1:nsim
    ))
    ss <- map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      sigma <- R$sigma
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      N <- R$N
      r <- R$r
      direct <- R$direct
      Nj <- N * f
      delta <- delta_j * f + delta_nj * (1 - f)
      set.seed(i + seed)
      xt_j <- rnorm(n = Nj * r / (1 + r), mean = delta_j, sd = sigma)
      xc_j <- rnorm(n = Nj / (1 + r), mean = 0, sd = sigma)
      xt_nj <- rnorm(n = (N - Nj) * r / (1 + r), mean = delta_nj, sd = sigma)
      xc_nj <- rnorm(n = (N - Nj) / (1 + r), mean = 0, sd = sigma)
      xt <- c(xt_j, xt_nj)
      xc <- c(xc_j, xc_nj)
      za <- if_else(direct == 1,
        (mean(xt) - mean(xc) + cut) /
          (sigma * sqrt(1 / length(xt) + 1 / length(xc))),
        (mean(xt) - mean(xc) - cut) /
          (sigma * sqrt(1 / length(xt) + 1 / length(xc)))
      )
      zj <- if_else(direct == 1,
        mean(xt_j) - mean(xc_j) - (mean(xt) - mean(xc)) + pi * cut,
        mean(xt_j) - mean(xc_j) - (mean(xt) - mean(xc)) - pi * cut
      )
      if (direct == -1) {
        za <- (-1) * za
        zj <- (-1) * zj
      }
      succ_a <- if_else(za > qnorm(1 - alpha), 1, 0)
      succ_j <- if_else(zj > 0, 1, 0)
      data.frame(
        delta_j, delta_nj, delta,
        sigma, f, pi, cut,
        alpha,
        N = N, r = r, direct,
        succ_a, succ_j
      )
    })
    res <- ss %>%
      group_by(
        delta_j, delta_nj, delta, sigma,
        f, pi, cut, alpha, N, r, direct,
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
