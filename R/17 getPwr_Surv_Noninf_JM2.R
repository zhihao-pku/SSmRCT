#' Title
#'
#' @param delta_i
#' @param fi
#' @param alpha
#' @param N
#' @param r
#' @param direct
#' @param lambda0_i
#' @param sim
#' @param nsim
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
#' f_set <- seq(0.1, 0.9, 0.1)
#' map_dfr(.x = 1:length(f_set), .f = function(i) {
#'   f <- f_set[i]
#'   res <- getPwr_Surv_Noninf_JM2(
#'     delta_i = c(log(1.1), log(1.0)),
#'     fi = c(f, 1 - f), cut = log(1.3),
#'     alpha = 0.025, N = 300, r = 1, direct = -1, sim = FALSE
#'   )$overall
#'   res$M <- "calc"
#'   res$f <- f
#'   res
#' })
getPwr_Surv_Noninf_JM2 <- function(delta_i, fi, cut, alpha, N, r, direct = -1, lambda0_i = 1, sim = FALSE, nsim = 1000, seed = 0) {
  if (!sim) {
    Ni <- N * fi
    gr <- 2 + r + 1 / r
    num <- length(delta_i)
    delta <- sum(delta_i * fi)
    sei <- sqrt(gr / Ni)
    ui <- if_else(direct == rep(1, num),
      (delta_i + cut) / sei,
      (delta_i - cut) / sei
    )
    se <- sqrt(gr / N)
    u <- if_else(direct == 1,
      (delta + cut) / se,
      (delta - cut) / se
    )
    M <- diag(num + 1)
    M[1, 2:(num + 1)] <- sqrt(fi)
    M[2:(num + 1), 1] <- sqrt(fi)
    if (direct == -1) {
      ui <- (-1) * ui
      u <- (-1) * u
    }
    p1 <- 1 - pnorm(q = qnorm(1 - alpha), mean = u, sd = 1)
    p2 <- prod(1 - pnorm(q = 0, mean = ui, sd = 1))
    p3 <- pmvnorm(
      lower = c(qnorm(1 - alpha), rep(0, num)),
      upper = rep(Inf, num + 1),
      mean = c(u, ui),
      sigma = M
    )
    p4 <- p3 / p1
    res <- data.frame(
      delta = delta,
      p1 = p1, p2 = p2, p3 = p3, p4 = p4
    )
    p_margin <- 1 - pnorm(q = 0, mean = ui, sd = 1)
    p_joint <- c()
    for (k in 1:num) {
      p_joint_ <- pmvnorm(
        lower = c(qnorm(1 - alpha), 0),
        upper = c(Inf, Inf),
        mean = c(u, ui[k]),
        sigma = M[c(1, k + 1), c(1, k + 1)]
      )
      p_joint <- c(p_joint, p_joint_)
    }
    p_conditional <- p_joint / p1
    L <- list(
      overall = res,
      p_margin = p_margin,
      p_joint = p_joint,
      p_conditional = p_conditional
    )
  }
  if (sim) {
    Ni <- N * fi
    delta <- sum(delta_i * fi)
    num <- length(delta_i)
    da <- data.frame()
    di <- NULL
    for (j in 1:nsim) {
      set.seed(j + seed)
      Xt <- matrix(NA, nrow = N * r / (1 + r), ncol = num)
      Xc <- matrix(NA, nrow = N / (1 + r), ncol = num)
      zi <- c()
      for (k in 1:num) {
        nt_k <- round(Ni[k] * r / (1 + r))
        nc_k <- round(Ni[k] / (1 + r))
        if (length(lambda0_i) == 1) {
          lambda0_i <- rep(lambda0_i, num)
        }
        Xt[1:nt_k, k] <- rexp(n = nt_k, rate = lambda0_i[k] * exp(delta_i[k]))
        Xc[1:nc_k, k] <- rexp(n = nc_k, rate = lambda0_i[k])
        dat_i <- data.frame(
          time = c(Xt[1:nt_k, k], Xc[1:nc_k, k]),
          status = 1,
          trt = c(rep(1, nt_k), rep(0, nc_k))
        )
        fit_i <- coxph(Surv(time, status) ~ trt, dat = dat_i)
        coef_i <- coef(fit_i)
        zi <- c(zi, coef_i)
      }
      Xt <- as.numeric(Xt)
      Xt <- Xt[!is.na(Xt)]
      Xc <- as.numeric(Xc)
      Xc <- Xc[!is.na(Xc)]
      dat_a <- data.frame(
        time = c(Xt, Xc),
        status = 1,
        trt = c(rep(1, length(Xt)), rep(0, length(Xc)))
      )
      fit_a <- coxph(Surv(time, status) ~ trt, data = dat_a)
      ss <- summary(fit_a)$coefficients
      za_u <- ss[, 1] + qnorm(1 - alpha) * ss[, 3]
      za_l <- ss[, 1] - qnorm(1 - alpha) * ss[, 3]
      if (direct == 1) {
        succ_a <- if_else(za_l > -cut, 1, 0)
        succ_i <- if_else(all(zi > -cut), 1, 0)
        succ_i_ <- if_else(zi > -cut, 1, 0)
      }
      if (direct == -1) {
        succ_a <- if_else(za_u < cut, 1, 0)
        succ_i <- if_else(all(zi < cut), 1, 0)
        succ_i_ <- if_else(zi < cut, 1, 0)
      }
      da <- bind_rows(da, c(succ_a = succ_a, succ_i = succ_i))
      di <- rbind(di, succ_i_)
    }
    res <- da %>%
      summarise(
        p1 = mean(succ_a, na.rm = T),
        p2 = mean(succ_i, na.rm = T),
        p3 = mean(succ_a & succ_i, na.rm = T),
        p4 = mean(succ_i[succ_a == 1], na.rm = T)
      ) %>%
      mutate(
        delta = delta
      ) %>%
      select(delta, p1, p2, p3, p4)
    p_margin <- colMeans(di)
    p_joint <- c()
    for (k in 1:num) {
      p_joint_ <- mean(da$succ_a & di[, k])
      p_joint <- c(p_joint, p_joint_)
    }
    p_conditional <- colMeans(di[da$succ_a == 1, ])
    L <- list(
      overall = res,
      p_margin = p_margin,
      p_joint = p_joint,
      p_conditional = p_conditional
    )
  }
  return(L)
}
