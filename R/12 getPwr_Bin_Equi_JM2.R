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
#'   res <- getPwr_Bin_Equi_JM2(
#'     pt_i = c(0.5, 0.6),
#'     pc_i = c(0.6, 0.6),
#'     fi = c(f, 1 - f), cut = 0.3,
#'     alpha = 0.025, N = 100, r = 1, sim = FALSE
#'   )$overall
#'   res$M <- "calc"
#'   res$f <- f
#'   res
#' })
getPwr_Bin_Equi_JM2 <- function(pt_i, pc_i, fi, cut, alpha, N, r, sim = FALSE, nsim = 1000, seed = 0) {
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
    num <- length(delta_i)
    sei <- sqrt(gr * sigma_i^2 / Ni)
    ui1 <- (delta_i + cut) / sei
    ui2 <- (delta_i - cut) / sei
    se <- sqrt(gr * sigma_a^2 / N)
    u1 <- (delta + cut) / se
    u2 <- (delta - cut) / se
    M <- diag(num + 1)
    M[1, ] <- c(1, sqrt(fi))
    M[, 1] <- c(1, sqrt(fi))
    M1 <- combine_matrix(M, M)
    p1 <- pmvnorm(
      lower = c(qnorm(1 - alpha), -Inf),
      upper = c(Inf, -qnorm(1 - alpha)),
      mean = c(u1, u2),
      corr = matrix(1, nrow = 2, ncol = 2)
    )
    p_margin <- c()
    for (k in 1:num) {
      p_margin_ <- pmvnorm(
        lower = c(0, -Inf),
        upper = c(Inf, 0),
        mean = c(ui1[k], ui2[k]),
        corr = matrix(1, nrow = 2, ncol = 2)
      )
      p_margin <- c(p_margin, p_margin_)
    }
    p2 <- prod(p_margin)
    p3 <- pmvnorm(
      lower = c(qnorm(1 - alpha), -Inf, rep(c(0, -Inf), num)),
      upper = c(Inf, -qnorm(1 - alpha), rep(c(Inf, 0), num)),
      mean = c(u1, u2, combine_vector(ui1, ui2)),
      sigma = M1
    )
    p4 <- p3 / p1
    res <- data.frame(
      delta = delta,
      p1 = p1, p2 = p2, p3 = p3, p4 = p4
    )
    p_joint <- c()
    for (k in 1:num) {
      M <- diag(2)
      M[1, ] <- c(1, sqrt(fi[k]))
      M[, 1] <- c(1, sqrt(fi[k]))
      M2 <- combine_matrix(M, M)
      p_joint_ <- pmvnorm(
        lower = c(qnorm(1 - alpha), -Inf, 0, -Inf),
        upper = c(Inf, -qnorm(1 - alpha), Inf, 0),
        mean = c(u1, u2, ui1[k], ui2[k]),
        sigma = M2
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
      zi1 <- ut - uc + cut
      zi2 <- ut - uc - cut
      za1 <- (mean(as.numeric(Xt), na.rm = T) -
        mean(as.numeric(Xc), na.rm = T) + cut) / sqrt(gr * sigma_a^2 / N)
      za2 <- (mean(as.numeric(Xt), na.rm = T) -
        mean(as.numeric(Xc), na.rm = T) - cut) / sqrt(gr * sigma_a^2 / N)
      succ_a <- if_else(za1 > qnorm(1 - alpha) & za2 < -qnorm(1 - alpha), 1, 0)
      succ_i <- if_else(all(zi1 > 0) & all(zi2 < 0), 1, 0)
      succ_i_ <- if_else(zi1 > 0 & zi2 < 0, 1, 0)
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
