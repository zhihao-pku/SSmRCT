#' Title
#'
#' @param delta_j
#' @param delta_nj
#' @param f
#' @param pi
#' @param cut
#' @param alpha
#' @param N
#' @param r
#' @param direct
#'
#' @return
#' @export
#'
#' @examples
#' library(ggplot2)
#' dat1 <- getPwr_Surv_Noninf_JM1(
#'   delta_j = log(1.1), delta_nj = log(1.0),
#'   f = seq(0.1, 0.9, 0.1), cut = log(1.3),
#'   pi = 0.5, alpha = 0.025, N = 400, r = 1,
#'   direct = -1, sim = F
#' )
#' dat1$M <- "calc"
#' dat2 <- getPwr_Surv_Noninf_JM1(
#'   delta_j = log(1.1), delta_nj = log(1.0),
#'   f = seq(0.1, 0.9, 0.1), cut = log(1.3),
#'   pi = 0.5, alpha = 0.025, N = 400, r = 1,
#'   direct = -1, sim = T
#' )
#' dat2$M <- "sim"
#'
#' dat <- bind_rows(dat1, dat2)
#' dat <- gather(
#'   data = dat, key = "pwr class", value = "pwr",
#'   p1, p2, p3, p4
#' )
#'
#' ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
#'   geom_point() +
#'   geom_line() +
#'   facet_wrap(vars(`pwr class`), nrow = 2) +
#'   scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
#'   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
#'   labs(x = "allocation ratio", color = "", linetype = "")
getPwr_Surv_Noninf_JM1 <- function(delta_j, delta_nj, f, pi, cut, alpha, N, r, direct = -1, lambda0_j = 1, lambda0_nj = 1, sim = F, nsim = 1000, seed = 0) {
  if (!sim) {
    eg <- as.data.frame(expand.grid(
      delta_j = delta_j,
      delta_nj = delta_nj,
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
      se <- sqrt(gr / N)
      u <- if_else(direct == 1,
        (delta + cut) / se,
        (delta - cut) / se
      )
      sej <- sqrt(gr / Nj +
        gr / N -
        2 * sqrt(f) *
          sqrt(gr / Nj * gr / N))
      uj <- if_else(direct == 1,
        (delta_j - delta + pi * cut) / sej,
        (delta_j - delta - pi * cut) / sej
      )
      cov <- sqrt(f) * sqrt(gr / N * gr / Nj) - gr / N
      corr <- cov / sqrt(sej * se)
      M <- matrix(c(1, corr, corr, 1), nrow = 2, byrow = T)
      if (direct == -1) {
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
        f, pi, cut,
        alpha, N, r, direct,
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
      cut = cut,
      alpha = alpha,
      N = N,
      r = r,
      direct = direct,
      nsim = 1:nsim
    ))
    ss <- map_dfr(.x = 1:nsim, .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      lambda0_j <- R$lambda0_j
      lambda0_nj <- R$lambda0_nj
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
      ss <- summary(fit_a)$coefficients
      za_u <- ss[, 1] + qnorm(1 - alpha) * ss[, 3]
      za_l <- ss[, 1] - qnorm(1 - alpha) * ss[, 3]
      zj <- if_else(direct == 1,
        coef_j - coef_a + pi * cut,
        coef_j - coef_a - pi * cut
      )
      if (direct == 1) {
        succ_a <- if_else(za_l > -cut, 1, 0)
        succ_j <- if_else(zj > 0, 1, 0)
      }
      if (direct == -1) {
        succ_a <- if_else(za_u < cut, 1, 0)
        succ_j <- if_else(zj < 0, 1, 0)
      }
      data.frame(
        delta_j, delta_nj, delta,
        lambda0_j, lambda0_nj,
        f, pi, cut, alpha,
        N, r, direct,
        succ_a, succ_j
      )
    })
    res <- ss %>%
      group_by(
        delta_j, delta_nj, delta,
        lambda0_j, lambda0_nj,
        f, pi, cut, alpha, N, r, direct
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
