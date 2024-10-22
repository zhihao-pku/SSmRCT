#' Title
#'
#' @param delta_j
#' @param delta_nj
#' @param f
#' @param pi
#' @param alpha
#' @param N
#' @param r
#' @param criterion
#' @param direct
#' @param lambda0_j
#' @param lambda0_nj
#' @param sim
#' @param nsim
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
#' getPwr_Surv_Super_JM1(
#'   delta_j = log(0.8), delta_nj = log(0.6),
#'   f = seq(0.1, 0.9, 0.1),
#'   pi = 0.5, alpha = 0.025, N = 200, r = 1,
#'   criterion = 1, direct = -1, sim = FALSE
#' )
getPwr_Surv_Super_JM1 <- function(delta_j, delta_nj, f, pi, alpha, N, r, criterion, direct = -1, lambda0_j = 1, lambda0_nj = 1, sim = FALSE, nsim = 1000, seed = 0) {
  if (!sim) {
    eg <- as.data.frame(expand.grid(
      delta_j = delta_j,
      delta_nj = delta_nj,
      f = f,
      pi = pi,
      alpha = alpha,
      N = N,
      r = r,
      criterion = criterion,
      direct = direct
    ))
    res <- map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      f <- R$f
      pi <- R$pi
      alpha <- R$alpha
      N <- R$N
      r <- R$r
      criterion <- R$criterion
      direct <- R$direct
      Nj <- N * f
      gr <- 2 + r + 1 / r
      delta <- delta_j * f + delta_nj * (1 - f)
      if (criterion == 1) {
        se <- sqrt(gr / N)
        u <- delta / se
        sej <- sqrt(gr / Nj +
          pi^2 * gr / N -
          2 * pi * sqrt(f) *
            sqrt(gr / Nj * gr / N))
        uj <- (delta_j - pi * delta) / sej
        cov <- sqrt(f) * sqrt(gr / N * gr / Nj) -
          pi * gr / N
        corr <- cov / sqrt(sej * se)
      }
      if (criterion == 2) {
        se <- sqrt(gr / N * exp(delta)^2)
        u <- -(1 - exp(delta)) / se
        sej <- sqrt(gr / Nj * exp(delta_j)^2 +
          pi^2 * gr / N * exp(delta)^2 -
          2 * pi * gr / N * exp(delta_j) * exp(delta))
        uj <- -(1 - exp(delta_j) - pi * (1 - exp(delta))) / sej
        cov <- gr / N * exp(delta_j) * exp(delta) - pi * gr / N * exp(delta)^2
        corr <- cov / sqrt(sej * se)
      }
      if (direct == -1) {
        u <- (-1) * u
        uj <- (-1) * uj
      }
      M <- matrix(c(1, corr, corr, 1), nrow = 2, byrow = T)
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
        f, pi,
        alpha, N, r, criterion, direct,
        p1, p2, p3, p4
      )
    })
  }
  if (sim) {
    eg <- as.data.frame(expand.grid(
      delta_j = delta_j,
      delta_nj = delta_nj,
      lambda0_j = lambda0_j,
      lambda0_nj = lambda0_nj,
      f = f,
      pi = pi,
      alpha = alpha,
      N = N,
      r = r,
      criterion = criterion,
      direct = direct,
      nsim = 1:nsim
    ))
    ss <- map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      lambda0_j <- R$lambda0_j
      lambda0_nj <- R$lambda0_nj
      f <- R$f
      pi <- R$pi
      alpha <- R$alpha
      N <- R$N
      r <- R$r
      criterion <- R$criterion
      direct <- R$direct
      Nj <- N * f
      delta <- delta_j * f + delta_nj * (1 - f)
      set.seed(i + seed)
      xt_j <- rexp(n = Nj * r / (1 + r), rate = lambda0_j * exp(delta_j))
      xc_j <- rexp(n = Nj / (1 + r), rate = lambda0_j)
      xt_nj <- rexp(n = (N - Nj) * r / (1 + r), rate = lambda0_nj * exp(delta_nj))
      xc_nj <- rexp(n = (N - Nj) / (1 + r), rate = lambda0_nj)
      xt <- c(xt_j, xt_nj)
      xc <- c(xc_j, xc_nj)
      dat_a <- data.frame(
        time = c(xt, xc), status = 1,
        trt = c(rep(1, length(xt)), rep(0, length(xc)))
      )
      dat_j <- data.frame(
        time = c(xt_j, xc_j), status = 1,
        trt = c(rep(1, length(xt_j)), rep(0, length(xc_j)))
      )
      fit_j <- coxph(Surv(time, status) ~ trt, dat = dat_j)
      coef_j <- coef(fit_j)
      fit_a <- coxph(Surv(time, status) ~ trt, dat = dat_a)
      coef_a <- coef(fit_a)
      za <- summary(fit_a)$coefficients[4]
      if (criterion == 1) {
        zj <- coef_j - pi * coef_a
      }
      if (criterion == 2) {
        zj <- -(1 - exp(coef_j) - pi * (1 - exp(coef_a)))
      }
      if (direct == -1) {
        za <- (-1) * za
        zj <- (-1) * zj
      }
      succ_a <- if_else(za > qnorm(1 - alpha), 1, 0)
      succ_j <- if_else(zj > 0, 1, 0)
      data.frame(
        delta_j, delta_nj, delta,
        lambda0_j, lambda0_nj,
        f, pi, alpha,
        N, r, criterion, direct,
        succ_a, succ_j
      )
    })
    res <- ss %>%
      group_by(
        delta_j, delta_nj, delta,
        lambda0_j, lambda0_nj,
        f, pi, alpha, N, r, criterion, direct
      ) %>%
      summarise(
        p1 = mean(succ_a, na.rm = T),
        p2 = mean(succ_j, na.rm = T),
        p3 = mean(succ_a & succ_j, na.rm = T),
        p4 = mean(succ_j[succ_a == 1], na.rm = T)
      ) %>%
      arrange(f) %>%
      as.data.frame()
  }
  return(res)
}
