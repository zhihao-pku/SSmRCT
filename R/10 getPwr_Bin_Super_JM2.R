#' Title
#'
#' @param pt_i
#' @param pc_i
#' @param fi
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
#' f_set <- seq(0.1, 0.9, 0.1)
#' map_dfr(.x = 1:length(f_set), .f = function(i) {
#'   f <- f_set[i]
#'   res <- getPwr_Bin_Super_JM2(
#'     pt_i = c(0.3, 0.4),
#'     pc_i = c(0.6, 0.6),
#'     fi = c(f, 1 - f),
#'     alpha = 0.025, N = 100, r = 1, direct = -1, sim = FALSE
#'   )$overall
#'   res$M <- "calc"
#'   res$f <- f
#'   res
#' })
getPwr_Bin_Super_JM2 <- function(pt_i, pc_i, fi, alpha, N, r, direct = 1, sim = FALSE, nsim = 1000, seed = 0) {
  if (!sim) {
    Ni <- N * fi
    gr <- 2 + r + 1 / r
    delta_i <- pt_i - pc_i
    pt <- sum(pt_i * fi)
    pc <- sum(pc_i * fi)
    delta <- pt - pc
    p_i <- pt_i * r / (1 + r) + pc_i / (1 + r)
    p_a <- pt * r / (1 + r) + pc / (1 + r)
    sigma_i <- sqrt(p_i * (1 - p_i))
    sigma_a <- sqrt(p_a * (1 - p_a))
    sei <- sqrt(gr * sigma_i^2 / Ni)
    ui <- delta_i / sei
    se <- sqrt(gr * sigma_a^2 / N)
    u <- delta / se
    num <- length(delta_i)
    M <- diag(num + 1)
    M[1, ] <- c(1, sqrt(fi))
    M[, 1] <- c(1, sqrt(fi))
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
    gr <- 2 + r + 1 / r
    delta_i <- pt_i - pc_i
    pt <- sum(pt_i * fi)
    pc <- sum(pc_i * fi)
    delta <- pt - pc
    num <- length(delta_i)
    da <- data.frame()
    di <- NULL
    for (j in 1:nsim) {
      set.seed(j + seed)
      Xt <- matrix(NA, nrow = N * r / (1 + r), ncol = num)
      Xc <- matrix(NA, nrow = N / (1 + r), ncol = num)
      for (k in 1:num) {
        nt_k <- round(Ni[k] * r / (1 + r))
        nc_k <- round(Ni[k] / (1 + r))
        Xt[1:nt_k, k] <- rbinom(n = nt_k, size = 1, prob = pt_i[k])
        Xc[1:nc_k, k] <- rbinom(n = nc_k, size = 1, prob = pc_i[k])
      }
      ut <- colMeans(Xt, na.rm = T)
      uc <- colMeans(Xc, na.rm = T)
      sigma_a <- sqrt(mean(c(as.numeric(Xt), as.numeric(Xc)), na.rm = T) *
        (1 - mean(c(as.numeric(Xt), as.numeric(Xc)), na.rm = T)))
      zi <- ut - uc
      za <- (mean(as.numeric(Xt), na.rm = T) -
        mean(as.numeric(Xc), na.rm = T)) / sqrt(gr * sigma_a^2 / N)
      if (direct == -1) {
        zi <- (-1) * zi
        za <- (-1) * za
      }
      succ_a <- if_else(za > qnorm(1 - alpha), 1, 0)
      succ_i <- if_else(all(zi > 0), 1, 0)
      succ_i_ <- if_else(zi > 0, 1, 0)
      da <- bind_rows(da, c(succ_a = succ_a, succ_i = succ_i))
      di <- rbind(di, succ_i_)
    }
    res <- da %>%
      summarise(
        p1 = mean(succ_a),
        p2 = mean(succ_i),
        p3 = mean(succ_a & succ_i),
        p4 = mean(succ_i[succ_a == 1])
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
