#' Title
#'
#' @param p1_j
#' @param p0_j
#' @param p1_nj
#' @param p0_nj
#' @param f
#' @param pi
#' @param cut
#' @param alpha
#' @param N
#' @param r
#' @param sim
#' @param nsim
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
#' library(ggplot2)
#' dat1 <- getPwr_Bin_Equi_JM1(
#'   p1_j = 0.55, p0_j = 0.65, p1_nj = 0.65, p0_nj = 0.65,
#'   f = seq(0.1, 0.9, 0.1), pi = 0.5,
#'   cut = 0.2, alpha = 0.025, N = 400, r = 1, sim = F
#' )
#' dat1$M <- "calc"
#' dat2 <- getPwr_Bin_Equi_JM1(
#'   p1_j = 0.55, p0_j = 0.65, p1_nj = 0.65, p0_nj = 0.65,
#'   f = seq(0.1, 0.9, 0.1), pi = 0.5,
#'   cut = 0.2, alpha = 0.025, N = 400, r = 1, sim = T
#' )
#' dat2$M <- "sim"
#'
#' dat <- bind_rows(dat1, dat2)
#' dat <- gather(
#'   data = dat, key = "pwr class",
#'   value = "pwr", p1, p2, p3, p4
#' )
#'
#' ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
#'   geom_point() +
#'   geom_line() +
#'   facet_wrap(vars(`pwr class`), nrow = 2) +
#'   scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
#'   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
#'   labs(x = "allocation ratio", color = "", linetype = "")
getPwr_Bin_Equi_JM1 <- function(p1_j, p0_j, p1_nj, p0_nj, f, pi, cut, alpha, N, r, sim = F, nsim = 1000, seed = 0) {
  if (!sim) {
    eg <- as.data.frame(expand.grid(
      p1_j = p1_j,
      p0_j = p0_j,
      p1_nj = p1_nj,
      p0_nj = p0_nj,
      f = f,
      pi = pi,
      cut = cut,
      alpha = alpha,
      N = N,
      r = r
    ))
    res <- map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      p1_j <- R$p1_j
      p0_j <- R$p0_j
      p1_nj <- R$p1_nj
      p0_nj <- R$p0_nj
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      N <- R$N
      r <- R$r
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
      se <- sqrt(gr * sigma_a^2 / N)
      u1 <- (delta + cut) / se
      u2 <- (delta - cut) / se
      sej <- sqrt(gr * sigma_j^2 / Nj +
        gr * sigma_a^2 / N -
        2 * sqrt(f) *
          sqrt(gr * sigma_j^2 / Nj * gr * sigma_a^2 / N))
      uj1 <- (delta_j - delta + pi * cut) / sej
      uj2 <- (delta_j - delta - pi * cut) / sej
      cov <- sqrt(f) * sqrt(gr * sigma_a^2 / N * gr * sigma_j^2 / Nj) -
        gr * sigma_a^2 / N
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
        p1_j, p0_j,
        p1_nj, p0_nj,
        p1_a, p0_a,
        f, pi, cut,
        alpha, N, r,
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
      cut = cut,
      alpha = alpha,
      N = N,
      r = r,
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
      cut <- R$cut
      alpha <- R$alpha
      N <- R$N
      r <- R$r
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
      p_a <- mean(c(xc, xt))
      p_j <- mean(c(xt_j, xc_j))
      sigma_j <- sqrt(p_j * (1 - p_j))
      sigma_a <- sqrt(p_a * (1 - p_a))
      za1 <- (mean(xt) - mean(xc) + cut) /
        (sigma_a * sqrt(1 / length(xt) + 1 / length(xc)))
      za2 <- (mean(xt) - mean(xc) - cut) /
        (sigma_a * sqrt(1 / length(xt) + 1 / length(xc)))
      zj1 <- mean(xt_j) - mean(xc_j) - (mean(xt) - mean(xc)) + pi * cut
      zj2 <- mean(xt_j) - mean(xc_j) - (mean(xt) - mean(xc)) - pi * cut
      succ_a <- if_else(za1 > qnorm(1 - alpha) & za2 < -qnorm(1 - alpha), 1, 0)
      succ_j <- if_else(zj1 > 0 & zj2 < 0, 1, 0)
      data.frame(
        p1_j, p0_j,
        p1_nj, p0_nj,
        p1_a, p0_a,
        f, pi, cut,
        alpha, N, r,
        succ_a, succ_j
      )
    })
    res <- ss %>%
      group_by(
        p1_j, p0_j, p1_nj, p0_nj, p1_a, p0_a, f,
        pi, cut, alpha, N, r
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
