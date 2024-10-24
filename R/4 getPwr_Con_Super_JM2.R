#' Title
#'
#' @param delta_i
#' @param sigma
#' @param fi
#' @param alpha
#' @param beta
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
# f_set <- seq(0.1, 0.9, 0.1)
# map_dfr(.x = 1:length(f_set), .f = function(i) {
#   f <- f_set[i]
#   res <- getPwr_Con_Super_JM2(
#     delta_i = c(1, 0.8), sigma = 4,
#     fi = c(f, 1 - f),
#     alpha = 0.025, beta = NA, N = 200, r = 1, sim = FALSE
#   )$overall
#   res$M <- "calc"
#   res$f <- f
#   res
# })
# f_set <- seq(0.1, 0.9, 0.1)
# map_dfr(.x = 1:length(f_set), .f = function(i) {
#   f <- f_set[i]
#   res <- getPwr_Con_Super_JM2(
#     delta_i = c(1, 0.8), sigma = 4,
#     fi = c(f, 1 - f),
#     alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = FALSE
#   )$overall
#   res$M <- "calc"
#   res$f <- f
#   res
# })
getPwr_Con_Super_JM2 <- function(delta_i, sigma, fi, alpha = 0.025, beta = NA, N, r = 1, sim = FALSE, nsim = 1000, seed = 0) {
  if (!sim) {
    gr <- 2 + r + 1 / r
    delta <- sum(delta_i * fi)
    if (is.na(N)) {
      N <- getN_Con_Super(
        delta, sigma, alpha,
        beta, NA, r
      )$N
    }
    Ni <- N * fi
    sei <- sqrt(gr * sigma^2 / Ni)
    ui <- delta_i / sei
    se <- sqrt(gr * sigma^2 / N)
    u <- delta / se
    num <- length(delta_i)
    M <- diag(num + 1)
    M[1, ] <- c(1, sqrt(fi))
    M[, 1] <- c(1, sqrt(fi))
    if (delta < 0) {
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
      delta = delta, N = N,
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
    gr <- 2 + r + 1 / r
    delta <- sum(delta_i * fi)
    if (is.na(N)) {
      N <- getN_Con_Super(
        delta, sigma, alpha,
        beta, NA, r
      )$N
    }
    Ni <- N * fi
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
        Xt[1:nt_k, k] <- rnorm(n = nt_k, mean = delta_i[k], sd = sigma)
        Xc[1:nc_k, k] <- rnorm(n = nc_k, mean = 0, sd = sigma)
      }
      ut <- colMeans(Xt, na.rm = T)
      uc <- colMeans(Xc, na.rm = T)
      zi <- ut - uc
      za <- (mean(as.numeric(Xt), na.rm = T) -
        mean(as.numeric(Xc), na.rm = T)) / sqrt(gr * sigma^2 / N)
      if (delta < 0) {
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
        delta = delta,
        N = N,
      ) %>%
      select(delta, N, p1, p2, p3, p4)
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
