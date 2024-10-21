#' Title
#'
#' @param p1_j
#' @param p0_j
#' @param p1_nj
#' @param p0_nj
#' @param f
#' @param pi
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
#' library(ggplot2)
#' dat1 <- getPwr_Bin_Super_JM1(
#'   p1_j = 0.35, p0_j = 0.5, p1_nj = 0.25, p0_nj = 0.5,
#'   f = seq(0.1, 0.9, 0.1),
#'   pi = 0.5, alpha = 0.025, N = 200, r = 1, direct = -1, sim = F
#' )
#' dat1$M <- "calc"
#' dat2 <- getPwr_Bin_Super_JM1(
#'   p1_j = 0.35, p0_j = 0.5, p1_nj = 0.25, p0_nj = 0.5,
#'   f = seq(0.1, 0.9, 0.1),
#'   pi = 0.5, alpha = 0.025, N = 200, r = 1, direct = -1, sim = T
#' )
#' dat2$M <- "sim"
#'
#' dat <- bind_rows(dat1, dat2)
#' dat <- gather(data = dat, key = "pwr class", value = "pwr", p1, p2, p3, p4)
#'
#' ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
#'   geom_point() +
#'   geom_line() +
#'   facet_wrap(vars(`pwr class`), nrow = 2) +
#'   scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
#'   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
#'   labs(x = "allocation ratio", linetype = "")
getPwr_Bin_Super_JM1 <- function(p1_j, p0_j, p1_nj, p0_nj, f, pi, alpha, N, r, direct, sim = F, nsim = 1000, seed = 0) {
  if (!sim) {
    eg <- as.data.frame(expand.grid(
      p1_j = p1_j,
      p0_j = p0_j,
      p1_nj = p1_nj,
      p0_nj = p0_nj,
      f = f,
      pi = pi,
      alpha = alpha,
      N = N,
      r = r,
      direct = direct
    ))
    res <- map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      p1_j <- R$p1_j
      p0_j <- R$p0_j
      p1_nj <- R$p1_nj
      p0_nj <- R$p0_nj
      f <- R$f
      pi <- R$pi
      alpha <- R$alpha
      N <- R$N
      r <- R$r
      direct <- R$direct
      Nj <- N * f
      gr <- 2 + r + 1 / r
      p0_a <- p0_j * f + p0_nj * (1 - f)
      p1_a <- p1_j * f + p1_nj * (1 - f)
      p_j <- p0_j / (1 + r) + p1_j * r / (1 + r)
      p_a <- p0_a / (1 + r) + p1_a * r / (1 + r)
      sigma_j <- sqrt(p_j * (1 - p_j))
      sigma_a <- sqrt(p_a * (1 - p_a))
      delta_j <- p1_j - p0_j
      delta <- p1_a - p0_a
      sej <- sqrt(gr * sigma_j^2 / Nj +
        pi^2 * gr * sigma_a^2 / N -
        2 * pi * sqrt(f) *
          sqrt(gr * sigma_j^2 / Nj * gr * sigma_a^2 / N))
      se <- sqrt(gr * sigma_a^2 / N)
      uj <- (delta_j - pi * delta) / sej
      u <- delta / se
      cov <- sqrt(f) * sqrt(gr * sigma_a^2 / N * gr * sigma_j^2 / Nj) -
        pi * gr * sigma_a^2 / N
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
        p1_j, p0_j,
        p1_nj, p0_nj,
        p1_a, p0_a,
        f, pi,
        alpha, N, r, direct,
        p1, p2, p3, p4
      )
    })
  }
  if (sim) {
    eg <- as.data.frame(expand.grid(
      p1_j = p1_j,
      p0_j = p0_j,
      p1_nj = p1_nj,
      p0_nj = p0_nj,
      f = f,
      pi = pi,
      alpha = alpha,
      N = N,
      r = r,
      direct = direct,
      nsim = 1:nsim
    ))
    ss <- map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      p1_j <- R$p1_j
      p0_j <- R$p0_j
      p1_nj <- R$p1_nj
      p0_nj <- R$p0_nj
      f <- R$f
      pi <- R$pi
      alpha <- R$alpha
      N <- R$N
      r <- R$r
      direct <- R$direct
      Nj <- N * f
      p1_a <- p1_j * f + p1_nj * (1 - f)
      p0_a <- p0_j * f + p0_nj * (1 - f)
      set.seed(i + seed)
      xt_j <- rbinom(n = Nj * r / (1 + r), size = 1, prob = p1_j)
      xc_j <- rbinom(n = Nj / (1 + r), size = 1, prob = p0_j)
      xt_nj <- rbinom(n = (N - Nj) * r / (1 + r), size = 1, prob = p1_nj)
      xc_nj <- rbinom(n = (N - Nj) / (1 + r), size = 1, prob = p0_nj)
      xt <- c(xt_j, xt_nj)
      xc <- c(xc_j, xc_nj)
      p_a <- mean(c(xt, xc))
      p_j <- mean(c(xt_j, xc_j))
      sigma_j <- sqrt(p_j * (1 - p_j))
      sigma_a <- sqrt(p_a * (1 - p_a))
      za <- (mean(xt) - mean(xc)) /
        (sigma_a * sqrt(1 / length(xt) + 1 / length(xc)))
      zj <- mean(xt_j) - mean(xc_j) - pi * (mean(xt) - mean(xc))
      if (direct == -1) {
        za <- (-1) * za
        zj <- (-1) * zj
      }
      succ_a <- if_else(za > qnorm(1 - alpha), 1, 0)
      succ_j <- if_else(zj > 0, 1, 0)
      data.frame(
        p1_j, p0_j,
        p1_nj, p0_nj,
        p1_a, p0_a,
        f, pi, alpha,
        N, r, direct,
        succ_a, succ_j
      )
    })
    res <- ss %>%
      group_by(
        p1_j, p0_j, p1_nj, p0_nj, p1_a, p0_a, f, pi, alpha,
        N, r, direct
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
