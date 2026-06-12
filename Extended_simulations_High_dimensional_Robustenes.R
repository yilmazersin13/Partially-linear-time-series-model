
options(stringsAsFactors = FALSE)
set.seed(20260611)



OUT_DIR <- "outputs_reviewer_extension"
TABLE_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

USE_PARALLEL <- TRUE
N_CORES <- 3
B <- 20

BW_SEQ <- c(seq(0.005, 0.05, by = 0.005), seq(0.06, 0.25, by = 0.02))
K_SEQ  <- c(0, 0.001, 0.005, 0.01, 0.025, 0.05, 0.10, 0.25,
            0.50, 1, 2, 5, 10, 20, 50, 100)
CRITERIA <- c("BIC")

HIGH_DIM_GRID <- expand.grid(
  n = c(150, 300),
  p = c(50, 100),
  delta = c(0.90, 0.99),
  rho = 0.60,
  f_id = 1,
  err_case = "AR1_Gauss",
  stringsAsFactors = FALSE
)

ROBUST_GRID <- expand.grid(
  n = c(150, 300),
  p = 6,
  delta = c(0.90, 0.99),
  rho = 0.60,
  f_id = 1,
  err_case = c("AR1_Gauss", "ARMA11_Gauss", "AR1_Hetero"),
  stringsAsFactors = FALSE
)

# Put GRTK first intentionally to highlight it in the figures.
ESTIMATORS <- c("GRTK", "RTK", "KS")
EST_COLS <- c(GRTK = "#1B9E77", RTK = "#D95F02", KS = "#7570B3")

FIG_DEVICE <- "png"   # use "png" or "jpeg"; no PDF output
SAVE_FIGURES <- TRUE



gaussian_kernel <- function(u) exp(-0.5 * u^2)

NW_smoother_matrix <- function(z, bw) {
  n <- length(z)
  D <- outer(z, z, "-") / bw
  K <- gaussian_kernel(D)
  rs <- rowSums(K) + 1e-12
  K / rs
}

AR1_corr_matrix <- function(n, rho) {
  rho <- max(min(as.numeric(rho), 0.999), -0.999)
  toeplitz(rho^(0:(n - 1)))
}

safe_chol_inv <- function(A) {
  R <- try(chol(A), silent = TRUE)
  if (inherits(R, "try-error")) return(solve(A))
  chol2inv(R)
}

safe_solve <- function(A, b = NULL, ridge = 1e-8) {
  out <- try(if (is.null(b)) solve(A) else solve(A, b), silent = TRUE)
  if (!inherits(out, "try-error") && all(is.finite(out))) return(out)
  
  p <- nrow(A)
  out <- try(if (is.null(b)) solve(A + ridge * diag(p)) else solve(A + ridge * diag(p), b), silent = TRUE)
  if (!inherits(out, "try-error") && all(is.finite(out))) return(out)
  
  if (is.null(b)) return(qr.solve(A + ridge * diag(p), diag(p)))
  qr.solve(A + ridge * diag(p), b)
}

estimate_rho_ar1 <- function(res, lag = 1, eta = 1e-3) {
  ac <- stats::acf(as.numeric(res), plot = FALSE, lag.max = lag)$acf
  rho_hat <- as.numeric(ac[2])
  min(1 - eta, max(-1 + eta, rho_hat))
}

toeplitz_sigma <- function(p, delta) {
  outer(1:p, 1:p, function(i, j) delta^abs(i - j))
}

gen_X_multicollinear <- function(n, p, delta) {
  Sigma <- toeplitz_sigma(p, delta)
  U <- chol(Sigma)
  Z <- matrix(rnorm(n * p), n, p)
  Z %*% U
}

make_beta <- function(p) {
  beta <- rep(0, p)
  beta[1:min(3, p)] <- c(1.0, 0.8, 0.5)[1:min(3, p)]
  beta
}

f_true <- function(z, which = 1) {
  if (which == 1) return(2 * sin(2 * pi * z))
  if (which == 2) return(4 * (z - 0.5)^2 - 1)
  stop("f(which) must be 1 or 2")
}

open_device <- function(file, width = 8, height = 5) {
  if (FIG_DEVICE == "png") {
    grDevices::png(filename = file, width = width, height = height, units = "in", res = 300)
  } else if (FIG_DEVICE == "jpeg") {
    grDevices::jpeg(filename = file, width = width, height = height, units = "in", res = 300, quality = 95)
  } else {
    stop("Use only PNG or JPEG output in this revision.")
  }
}

fig_file <- function(stem) {
  file.path(FIG_DIR, paste0(stem, ".", FIG_DEVICE))
}

pretty_err_label <- function(x) {
  out <- x
  out[out == "AR1_Gauss"] <- "AR(1), Gaussian"
  out[out == "ARMA11_Gauss"] <- "ARMA(1,1), Gaussian"
  out[out == "AR1_Hetero"] <- "AR(1), heteroscedastic"
  out
}

scenario_label_hd <- function(n, p, delta) paste0("n=", n, ", p=", p, ", delta=", delta)
scenario_label_rb <- function(n, delta, err_case) paste0("n=", n, ", delta=", delta, "\n", pretty_err_label(err_case))



gen_errors_AR1 <- function(n, rho, innov = c("gaussian", "t5"), burn = 300, sigma = 1) {
  innov <- match.arg(innov)
  N <- n + burn
  u <- switch(
    innov,
    gaussian = rnorm(N, 0, sigma),
    t5 = {
      df <- 5
      scale <- sqrt((df - 2) / df) * sigma
      stats::rt(N, df = df) * scale
    }
  )
  e <- as.numeric(stats::filter(u, filter = rho, method = "recursive", init = 0))
  e[(burn + 1):N]
}

gen_errors_ARMA11 <- function(n, ar = 0.60, ma = 0.40, burn = 500, sigma = 1) {
  N <- n + burn
  u <- rnorm(N, 0, sigma)
  e <- numeric(N)
  for (t in 2:N) {
    e[t] <- ar * e[t - 1] + u[t] + ma * u[t - 1]
  }
  e[(burn + 1):N]
}

gen_errors_AR1_hetero <- function(n, rho = 0.60, z = NULL, burn = 300) {
  if (is.null(z)) z <- (1:n) / n
  N <- n + burn
  z_full <- seq(0, 1, length.out = N)
  sigma_t <- 0.6 + 0.8 * z_full
  u <- rnorm(N, 0, sigma_t)
  e <- as.numeric(stats::filter(u, filter = rho, method = "recursive", init = 0))
  e[(burn + 1):N]
}

generate_dataset <- function(n, p, delta, rho, f_id = 1,
                             err_case = c("AR1_Gauss", "AR1_t5", "ARMA11_Gauss", "AR1_Hetero")) {
  err_case <- match.arg(err_case)
  
  z <- (1:n) / n
  x <- gen_X_multicollinear(n, p, delta)
  f <- f_true(z, which = f_id)
  beta <- make_beta(p)
  
  eps <- switch(
    err_case,
    AR1_Gauss    = gen_errors_AR1(n, rho, innov = "gaussian"),
    AR1_t5       = gen_errors_AR1(n, rho, innov = "t5"),
    ARMA11_Gauss = gen_errors_ARMA11(n, ar = rho, ma = 0.40),
    AR1_Hetero   = gen_errors_AR1_hetero(n, rho = rho, z = z)
  )
  
  y <- as.vector(x %*% beta + f + eps)
  
  list(y = y, x = x, z = z, f = f, beta = beta, eps = eps,
       n = n, p = p, delta = delta, rho = rho, f_id = f_id, err_case = err_case)
}



fit_plm_once <- function(x, z, y, bw, k, Rinv, W = NULL) {
  n <- length(y)
  p <- ncol(x)
  if (is.null(W)) W <- NW_smoother_matrix(z, bw)
  Iw <- diag(n) - W
  
  xtil <- Iw %*% x
  ytil <- Iw %*% y
  
  A <- crossprod(xtil, Rinv %*% xtil) + k * diag(p)
  b <- crossprod(xtil, Rinv %*% ytil)
  beta_hat <- as.vector(safe_solve(A, b))
  
  yhat <- as.vector(W %*% y + (Iw %*% x) %*% beta_hat)
  fhat <- as.vector(W %*% (y - x %*% beta_hat))
  
  list(beta_hat = beta_hat, fhat = fhat, yhat = yhat,
       W = W, Iw = Iw, xtil = xtil, A = A)
}

df_fast <- function(W, Iw, xtil, A, Rinv) {
  trW <- sum(diag(W))
  C <- crossprod(xtil, Rinv %*% (Iw %*% xtil))
  tr2 <- sum(diag(safe_solve(A, C)))
  as.numeric(trW + tr2)
}

H_matrix <- function(W, Iw, x, xtil, A, Rinv) {
  B <- safe_solve(A, crossprod(xtil, Rinv %*% Iw))
  H2 <- (Iw %*% x) %*% B
  W + H2
}

Selection_bw_k_safe <- function(x, z, y,
                                bw_seq, k_seq,
                                lag = 1,
                                criteria = c("BIC"),
                                R_mode = c("AR1", "I"),
                                k_pilot = NULL,
                                eta_rho = 1e-3) {
  R_mode <- match.arg(R_mode)
  n <- length(y)
  p <- ncol(x)
  
  if (is.null(k_pilot)) k_pilot <- median(k_seq[k_seq > 0])
  if (length(k_pilot) == 0L || is.na(k_pilot) || !is.finite(k_pilot)) k_pilot <- 0.1
  
  rho_hat <- 0
  Rinv <- diag(n)
  
  if (R_mode == "AR1") {
    bw_pilot <- median(bw_seq)
    pilot <- fit_plm_once(x, z, y, bw = bw_pilot, k = k_pilot, Rinv = diag(n))
    res_pilot <- y - pilot$yhat
    rho_hat <- estimate_rho_ar1(res_pilot, lag = lag, eta = eta_rho)
    Rinv <- safe_chol_inv(AR1_corr_matrix(n, rho_hat))
  }
  
  L <- length(bw_seq)
  M <- length(k_seq)
  W_cache <- vector("list", L)
  Iw_cache <- vector("list", L)
  xtil_cache <- vector("list", L)
  ytil_cache <- vector("list", L)
  Xt_Rinv_cache <- vector("list", L)
  
  for (i in seq_len(L)) {
    W_cache[[i]] <- NW_smoother_matrix(z, bw_seq[i])
    Iw_cache[[i]] <- diag(n) - W_cache[[i]]
    xtil_cache[[i]] <- Iw_cache[[i]] %*% x
    ytil_cache[[i]] <- Iw_cache[[i]] %*% y
    Xt_Rinv_cache[[i]] <- crossprod(xtil_cache[[i]], Rinv)
  }
  
  make_score_matrix <- function(criterion) {
    S <- matrix(NA_real_, L, M)
    
    for (i in seq_len(L)) {
      W <- W_cache[[i]]
      Iw <- Iw_cache[[i]]
      xtil <- xtil_cache[[i]]
      ytil <- ytil_cache[[i]]
      Xt_Rinv <- Xt_Rinv_cache[[i]]
      
      for (j in seq_len(M)) {
        k <- k_seq[j]
        A <- Xt_Rinv %*% xtil + k * diag(p)
        b <- Xt_Rinv %*% ytil
        beta_hat <- as.vector(safe_solve(A, b))
        yhat <- as.vector(W %*% y + (Iw %*% x) %*% beta_hat)
        RSS <- sum((y - yhat)^2)
        df <- df_fast(W, Iw, xtil, A, Rinv)
        df_res <- n - df
        if (!is.finite(df_res) || df_res <= 5) df_res <- 5
        
        if (criterion == "GCV") {
          S[i, j] <- (RSS / n) / (1 - df / n)^2
        }
        if (criterion == "AICc") {
          denom <- 1 - (df + 2) / n
          if (denom <= 0.01) denom <- 0.01
          S[i, j] <- log(RSS / n) + (1 + df / n) / denom
        }
        if (criterion == "BIC") {
          S[i, j] <- (RSS / n) * exp(df * log(n) / n)
        }
        if (criterion == "RECP") {
          sigma2_hat <- RSS / df_res
          H <- H_matrix(W, Iw, x, xtil, A, Rinv)
          trHHt <- sum(H * H)
          S[i, j] <- (RSS + sigma2_hat * trHHt) / df_res
        }
      }
    }
    S
  }
  
  pairs <- data.frame(Criterion = criteria, bw_hat = NA_real_, k_hat = NA_real_)
  score_list <- list()
  
  for (m in seq_along(criteria)) {
    S <- make_score_matrix(criteria[m])
    score_list[[criteria[m]]] <- S
    idx <- which(S == min(S, na.rm = TRUE), arr.ind = TRUE)[1, ]
    pairs$bw_hat[m] <- bw_seq[idx[1]]
    pairs$k_hat[m] <- k_seq[idx[2]]
  }
  
  list(pairs = pairs, rho_hat = rho_hat, scores = score_list, R_mode = R_mode)
}

fit_KS <- function(x, z, y, bw) {
  fit_plm_once(x, z, y, bw = bw, k = 0, Rinv = diag(length(y)))
}

fit_RTK <- function(x, z, y, bw, k) {
  fit_plm_once(x, z, y, bw = bw, k = k, Rinv = diag(length(y)))
}

fit_GRTK <- function(x, z, y, bw, k, lag = 1, bw_pilot = 0.15, k_pilot = 0.10) {
  n <- length(y)
  pilot <- fit_plm_once(x, z, y, bw = bw_pilot, k = k_pilot, Rinv = diag(n))
  rho_hat <- estimate_rho_ar1(y - pilot$yhat, lag = lag)
  Rinv <- safe_chol_inv(AR1_corr_matrix(n, rho_hat))
  fit <- fit_plm_once(x, z, y, bw = bw, k = k, Rinv = Rinv)
  fit$rho_hat <- rho_hat
  fit
}



compute_metrics <- function(fit, dat) {
  mu_true <- as.vector(dat$x %*% dat$beta + dat$f)
  beta_error <- as.vector(fit$beta_hat - dat$beta)
  active <- which(dat$beta != 0)
  inactive <- which(dat$beta == 0)
  
  beta_smsd_total <- sum(beta_error^2)
  beta_smsd_per_coef <- beta_smsd_total / length(dat$beta)
  active_mse <- mean(beta_error[active]^2)
  inactive_mse <- if (length(inactive) > 0) mean(beta_error[inactive]^2) else NA_real_
  
  list(
    beta_smsd_total = beta_smsd_total,
    beta_smsd_per_coef = beta_smsd_per_coef,
    active_beta_mse = active_mse,
    inactive_beta_mse = inactive_mse,
    f_mse = mean((fit$fhat - dat$f)^2),
    pred_mse = mean((fit$yhat - mu_true)^2)
  )
}

run_one_replication <- function(sc, rep_id, criteria = CRITERIA) {
  dat <- generate_dataset(
    n = sc$n, p = sc$p, delta = sc$delta, rho = sc$rho,
    f_id = sc$f_id, err_case = sc$err_case
  )
  
  out <- list()
  idx <- 1L
  bw_pilot <- median(BW_SEQ)
  k_pilot <- median(K_SEQ[K_SEQ > 0])
  
  for (crit in criteria) {
    for (est in ESTIMATORS) {
      fit <- NULL
      bw_hat <- NA_real_
      k_hat <- NA_real_
      rho_hat <- NA_real_
      ok <- TRUE
      msg <- NA_character_
      
      tmp <- try({
        if (est == "KS") {
          sel <- Selection_bw_k_safe(dat$x, dat$z, dat$y, BW_SEQ, k_seq = 0,
                                     criteria = crit, R_mode = "I")
          bw_hat <- sel$pairs$bw_hat[1]
          k_hat <- 0
          fit <- fit_KS(dat$x, dat$z, dat$y, bw = bw_hat)
        }
        if (est == "RTK") {
          sel <- Selection_bw_k_safe(dat$x, dat$z, dat$y, BW_SEQ, K_SEQ,
                                     criteria = crit, R_mode = "I")
          bw_hat <- sel$pairs$bw_hat[1]
          k_hat <- sel$pairs$k_hat[1]
          fit <- fit_RTK(dat$x, dat$z, dat$y, bw = bw_hat, k = k_hat)
        }
        if (est == "GRTK") {
          sel <- Selection_bw_k_safe(dat$x, dat$z, dat$y, BW_SEQ, K_SEQ,
                                     criteria = crit, R_mode = "AR1")
          bw_hat <- sel$pairs$bw_hat[1]
          k_hat <- sel$pairs$k_hat[1]
          fit <- fit_GRTK(dat$x, dat$z, dat$y, bw = bw_hat, k = k_hat,
                          bw_pilot = bw_pilot, k_pilot = k_pilot)
          rho_hat <- fit$rho_hat
        }
      }, silent = TRUE)
      
      if (inherits(tmp, "try-error")) {
        ok <- FALSE
        msg <- as.character(tmp)
        metrics <- list(beta_smsd_total = NA_real_, beta_smsd_per_coef = NA_real_,
                        active_beta_mse = NA_real_, inactive_beta_mse = NA_real_,
                        f_mse = NA_real_, pred_mse = NA_real_)
      } else {
        metrics <- compute_metrics(fit, dat)
      }
      
      out[[idx]] <- data.frame(
        rep_id = rep_id,
        n = sc$n, p = sc$p, delta = sc$delta, rho = sc$rho,
        f_id = sc$f_id, err_case = sc$err_case,
        criterion = crit, estimator = est,
        bw_hat = bw_hat, k_hat = k_hat, rho_hat = rho_hat,
        beta_smsd_total = metrics$beta_smsd_total,
        beta_smsd_per_coef = metrics$beta_smsd_per_coef,
        active_beta_mse = metrics$active_beta_mse,
        inactive_beta_mse = metrics$inactive_beta_mse,
        f_mse = metrics$f_mse,
        pred_mse = metrics$pred_mse,
        ok = ok,
        message = msg,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  
  do.call(rbind, out)
}

run_grid <- function(grid, label) {
  tasks <- vector("list", nrow(grid) * B)
  z <- 1L
  for (i in seq_len(nrow(grid))) {
    for (b in seq_len(B)) {
      tasks[[z]] <- list(sc = grid[i, ], rep_id = b)
      z <- z + 1L
    }
  }
  
  cat("\nRunning", label, "simulation with", length(tasks), "tasks.\n")
  start_time <- Sys.time()
  
  if (USE_PARALLEL && N_CORES > 1L) {
    cl <- parallel::makeCluster(N_CORES)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterSetRNGStream(cl, 20260611)
    
    export_names <- c(
      "B", "BW_SEQ", "K_SEQ", "CRITERIA", "ESTIMATORS",
      "gaussian_kernel", "NW_smoother_matrix", "AR1_corr_matrix",
      "safe_chol_inv", "safe_solve", "estimate_rho_ar1",
      "toeplitz_sigma", "gen_X_multicollinear", "make_beta", "f_true",
      "gen_errors_AR1", "gen_errors_ARMA11", "gen_errors_AR1_hetero",
      "generate_dataset", "fit_plm_once", "df_fast", "H_matrix",
      "Selection_bw_k_safe", "fit_KS", "fit_RTK", "fit_GRTK",
      "compute_metrics", "run_one_replication"
    )
    parallel::clusterExport(cl, varlist = export_names, envir = .GlobalEnv)
    
    results <- vector("list", length(tasks))
    batch_size <- max(1L, N_CORES)
    n_batches <- ceiling(length(tasks) / batch_size)
    
    for (bb in seq_len(n_batches)) {
      i1 <- (bb - 1L) * batch_size + 1L
      i2 <- min(bb * batch_size, length(tasks))
      results[i1:i2] <- parallel::parLapply(cl, tasks[i1:i2], function(task) {
        run_one_replication(task$sc, task$rep_id)
      })
      if (bb %% 5 == 0 || bb == n_batches) {
        cat(sprintf("  %s progress: %d/%d batches\n", label, bb, n_batches))
      }
    }
  } else {
    results <- lapply(tasks, function(task) run_one_replication(task$sc, task$rep_id))
  }
  
  res <- do.call(rbind, results)
  rownames(res) <- NULL
  cat(label, "simulation completed in", format(Sys.time() - start_time), "\n")
  res
}


summarise_results <- function(df) {
  group_vars <- c("n", "p", "delta", "rho", "f_id", "err_case", "criterion", "estimator")
  split_key <- interaction(df[group_vars], drop = TRUE)
  blocks <- split(df, split_key)
  
  res <- lapply(blocks, function(d) {
    base <- d[1, group_vars, drop = FALSE]
    data.frame(
      base,
      reps = sum(d$ok, na.rm = TRUE),
      fail = sum(!d$ok, na.rm = TRUE),
      bw_mean = mean(d$bw_hat, na.rm = TRUE),
      k_mean = mean(d$k_hat, na.rm = TRUE),
      rho_hat_mean = mean(d$rho_hat, na.rm = TRUE),
      beta_smsd_total_mean = mean(d$beta_smsd_total, na.rm = TRUE),
      beta_smsd_total_sd = stats::sd(d$beta_smsd_total, na.rm = TRUE),
      beta_smsd_per_coef_mean = mean(d$beta_smsd_per_coef, na.rm = TRUE),
      beta_smsd_per_coef_sd = stats::sd(d$beta_smsd_per_coef, na.rm = TRUE),
      active_beta_mse_mean = mean(d$active_beta_mse, na.rm = TRUE),
      inactive_beta_mse_mean = mean(d$inactive_beta_mse, na.rm = TRUE),
      f_mse_mean = mean(d$f_mse, na.rm = TRUE),
      f_mse_sd = stats::sd(d$f_mse, na.rm = TRUE),
      pred_mse_mean = mean(d$pred_mse, na.rm = TRUE),
      pred_mse_sd = stats::sd(d$pred_mse, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out[order(out$n, out$p, out$delta, out$err_case, match(out$estimator, ESTIMATORS)), ]
}

grtk_gain_table <- function(summary_df) {
  split_key <- interaction(summary_df$n, summary_df$p, summary_df$delta, summary_df$rho,
                           summary_df$err_case, summary_df$criterion, drop = TRUE)
  blocks <- split(summary_df, split_key)
  out <- lapply(blocks, function(d) {
    if (!all(c("GRTK", "RTK", "KS") %in% d$estimator)) return(NULL)
    dg <- d[d$estimator == "GRTK", ][1, ]
    dr <- d[d$estimator == "RTK", ][1, ]
    dk <- d[d$estimator == "KS", ][1, ]
    data.frame(
      n = dg$n,
      p = dg$p,
      delta = dg$delta,
      rho = dg$rho,
      err_case = dg$err_case,
      criterion = dg$criterion,
      gain_pred_vs_RTK_pct = 100 * (dr$pred_mse_mean - dg$pred_mse_mean) / dr$pred_mse_mean,
      gain_pred_vs_KS_pct  = 100 * (dk$pred_mse_mean - dg$pred_mse_mean) / dk$pred_mse_mean,
      gain_fmse_vs_RTK_pct = 100 * (dr$f_mse_mean - dg$f_mse_mean) / dr$f_mse_mean,
      gain_fmse_vs_KS_pct  = 100 * (dk$f_mse_mean - dg$f_mse_mean) / dk$f_mse_mean,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  if (!is.null(out)) rownames(out) <- NULL
  out
}

fmt_num <- function(x, digits = 4) ifelse(is.na(x), "--", formatC(x, format = "f", digits = digits))

write_highdim_latex_table <- function(tab, file) {
  lines <- c(
    "\\begin{table}[h!]",
    "\\centering",
    "\\caption{High-dimensional simulation results under severe multicollinearity using BIC-selected tuning parameters.}",
    "\\label{tab:highdim-extension}",
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{rrrrlrrrr}",
    "\\toprule",
    "$n$ & $p$ & $\\delta$ & $\\rho$ & Estimator & $\\bar{k}$ & SMSD$/p$ & MSE($\\widehat{f}$) & Pred. MSE \\",
    "\\midrule"
  )
  for (i in seq_len(nrow(tab))) {
    lines <- c(lines, paste(
      tab$n[i], "&", tab$p[i], "&", fmt_num(tab$delta[i], 2), "&", fmt_num(tab$rho[i], 2), "&",
      tab$estimator[i], "&", fmt_num(tab$k_mean[i], 3), "&",
      fmt_num(tab$beta_smsd_per_coef_mean[i], 4), "&",
      fmt_num(tab$f_mse_mean[i], 4), "&",
      fmt_num(tab$pred_mse_mean[i], 4), "\\\\"
    ))
  }
  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\begin{flushleft}",
    "\\footnotesize Note: SMSD$/p$ denotes the coefficient-wise scaled squared deviation of $\\widehat{\\boldsymbol{\\beta}}$ from the true coefficient vector. Pred. MSE is computed relative to the true conditional mean.",
    "\\end{flushleft}",
    "\\end{table}"
  )
  writeLines(lines, file)
}

write_robust_latex_table <- function(tab, file) {
  lines <- c(
    "\\begin{table}[h!]",
    "\\centering",
    "\\caption{Robustness simulation results for alternative error structures using BIC-selected tuning parameters.}",
    "\\label{tab:robustness-extension}",
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{rrrlrrrr}",
    "\\toprule",
    "$n$ & $\\delta$ & $\\rho$ & Error & Estimator & SMSD$/p$ & MSE($\\widehat{f}$) & Pred. MSE \\",
    "\\midrule"
  )
  for (i in seq_len(nrow(tab))) {
    lines <- c(lines, paste(
      tab$n[i], "&", fmt_num(tab$delta[i], 2), "&", fmt_num(tab$rho[i], 2), "&",
      tab$err_case[i], "&", tab$estimator[i], "&",
      fmt_num(tab$beta_smsd_per_coef_mean[i], 4), "&",
      fmt_num(tab$f_mse_mean[i], 4), "&",
      fmt_num(tab$pred_mse_mean[i], 4), "\\\\"
    ))
  }
  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\begin{flushleft}",
    "\\footnotesize Note: ARMA11\\_Gauss denotes an ARMA(1,1) error with moving-average coefficient 0.40; AR1\\_Hetero denotes an AR(1) error with time-varying innovation variance.",
    "\\end{flushleft}",
    "\\end{table}"
  )
  writeLines(lines, file)
}


plot_hd_boxpanels <- function(df, metric, file, ylab) {
  if (!SAVE_FIGURES) return(invisible(NULL))
  open_device(file, width = 11, height = 7)
  oldpar <- par(no.readonly = TRUE)
  on.exit({par(oldpar); dev.off()}, add = TRUE)
  
  scenarios <- unique(df[c("n", "p", "delta")])
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  for (i in seq_len(nrow(scenarios))) {
    ss <- scenarios[i, ]
    dd <- subset(df, n == ss$n & p == ss$p & delta == ss$delta)
    vals <- split(dd[[metric]], factor(dd$estimator, levels = ESTIMATORS))
    boxplot(vals, col = EST_COLS[ESTIMATORS], border = EST_COLS[ESTIMATORS],
            las = 1, outline = FALSE, main = scenario_label_hd(ss$n, ss$p, ss$delta),
            ylab = ylab, cex.axis = 0.9, cex.lab = 0.95)
    med <- sapply(vals, median, na.rm = TRUE)
    points(seq_along(med), med, pch = 19, cex = 0.9)
  }
  mtext("High-dimensional comparison across estimators", outer = TRUE, cex = 1.1, font = 2)
}

plot_hd_strip_panels <- function(df, metric, file, ylab) {
  if (!SAVE_FIGURES) return(invisible(NULL))
  open_device(file, width = 11, height = 7)
  oldpar <- par(no.readonly = TRUE)
  on.exit({par(oldpar); dev.off()}, add = TRUE)
  
  scenarios <- unique(df[c("n", "p", "delta")])
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  for (i in seq_len(nrow(scenarios))) {
    ss <- scenarios[i, ]
    dd <- subset(df, n == ss$n & p == ss$p & delta == ss$delta)
    grp <- factor(dd$estimator, levels = ESTIMATORS)
    stripchart(dd[[metric]] ~ grp, vertical = TRUE, method = "jitter", jitter = 0.15,
               pch = 16, cex = 0.65, col = grDevices::adjustcolor(EST_COLS[grp], alpha.f = 0.40),
               las = 1, main = scenario_label_hd(ss$n, ss$p, ss$delta), ylab = ylab)
    boxplot(dd[[metric]] ~ grp, add = TRUE, outline = FALSE, axes = FALSE,
            col = grDevices::adjustcolor(EST_COLS[ESTIMATORS], alpha.f = 0.20), border = EST_COLS[ESTIMATORS])
  }
  mtext("High-dimensional jittered distributions", outer = TRUE, cex = 1.1, font = 2)
}

plot_summary_profile_hd <- function(summary_df, metric, file, ylab) {
  if (!SAVE_FIGURES) return(invisible(NULL))
  
  metric_col <- metric
  if (!metric_col %in% names(summary_df)) {
    candidate <- paste0(metric, "_mean")
    if (candidate %in% names(summary_df)) metric_col <- candidate
  }
  if (!metric_col %in% names(summary_df)) {
    warning("Skipping profile figure: metric column not found: ", metric)
    return(invisible(NULL))
  }
  
  summary_df <- summary_df[is.finite(summary_df[[metric_col]]), , drop = FALSE]
  if (nrow(summary_df) == 0L) {
    warning("Skipping profile figure: no finite values for ", metric_col)
    return(invisible(NULL))
  }
  
  open_device(file, width = 10, height = 6)
  oldpar <- par(no.readonly = TRUE)
  on.exit({par(oldpar); dev.off()}, add = TRUE)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  
  for (pp in sort(unique(summary_df$p))) {
    dd <- subset(summary_df, p == pp)
    yr <- range(dd[[metric_col]], na.rm = TRUE)
    if (!all(is.finite(yr))) yr <- c(0, 1)
    if (diff(yr) == 0) yr <- yr + c(-0.05, 0.05) * max(1, abs(yr[1]))
    
    plot(NA, xlim = range(dd$n), ylim = yr,
         xlab = "Sample size n", ylab = ylab, main = paste0("p = ", pp), xaxt = "n")
    axis(1, at = sort(unique(dd$n)))
    
    for (est in ESTIMATORS) {
      de <- subset(dd, estimator == est & delta == 0.90)
      if (nrow(de) > 0L) lines(de$n, de[[metric_col]], type = "b", lwd = 2, pch = 16, col = EST_COLS[est])
      de2 <- subset(dd, estimator == est & delta == 0.99)
      if (nrow(de2) > 0L) lines(de2$n, de2[[metric_col]], type = "b", lwd = 2, lty = 2, pch = 1, col = EST_COLS[est])
    }
    
    legend("topright", legend = c("GRTK, delta=0.90", "GRTK, delta=0.99",
                                  "RTK, delta=0.90", "RTK, delta=0.99",
                                  "KS, delta=0.90", "KS, delta=0.99"),
           col = c(EST_COLS["GRTK"], EST_COLS["GRTK"], EST_COLS["RTK"], EST_COLS["RTK"], EST_COLS["KS"], EST_COLS["KS"]),
           lty = c(1, 2, 1, 2, 1, 2), pch = c(16, 1, 16, 1, 16, 1), cex = 0.8, bty = "n")
  }
  mtext("Mean performance profiles across sample sizes", outer = TRUE, cex = 1.1, font = 2)
}

plot_robust_boxpanels <- function(df, metric, file, ylab) {
  if (!SAVE_FIGURES) return(invisible(NULL))
  open_device(file, width = 12, height = 7)
  oldpar <- par(no.readonly = TRUE)
  on.exit({par(oldpar); dev.off()}, add = TRUE)
  
  scenarios <- unique(df[c("n", "delta", "err_case")])
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  for (i in seq_len(nrow(scenarios))) {
    ss <- scenarios[i, ]
    dd <- subset(df, n == ss$n & delta == ss$delta & err_case == ss$err_case)
    vals <- split(dd[[metric]], factor(dd$estimator, levels = ESTIMATORS))
    boxplot(vals, col = EST_COLS[ESTIMATORS], border = EST_COLS[ESTIMATORS],
            las = 1, outline = FALSE, main = scenario_label_rb(ss$n, ss$delta, ss$err_case),
            ylab = ylab, cex.axis = 0.85, cex.main = 0.85)
    med <- sapply(vals, median, na.rm = TRUE)
    points(seq_along(med), med, pch = 19, cex = 0.9)
  }
  mtext("Robustness comparison across error structures", outer = TRUE, cex = 1.1, font = 2)
}

plot_robust_strip_panels <- function(df, metric, file, ylab) {
  if (!SAVE_FIGURES) return(invisible(NULL))
  open_device(file, width = 12, height = 7)
  oldpar <- par(no.readonly = TRUE)
  on.exit({par(oldpar); dev.off()}, add = TRUE)
  
  scenarios <- unique(df[c("n", "delta", "err_case")])
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  for (i in seq_len(nrow(scenarios))) {
    ss <- scenarios[i, ]
    dd <- subset(df, n == ss$n & delta == ss$delta & err_case == ss$err_case)
    grp <- factor(dd$estimator, levels = ESTIMATORS)
    stripchart(dd[[metric]] ~ grp, vertical = TRUE, method = "jitter", jitter = 0.15,
               pch = 16, cex = 0.65, col = grDevices::adjustcolor(EST_COLS[grp], alpha.f = 0.40),
               las = 1, main = scenario_label_rb(ss$n, ss$delta, ss$err_case), ylab = ylab,
               cex.main = 0.85)
    boxplot(dd[[metric]] ~ grp, add = TRUE, outline = FALSE, axes = FALSE,
            col = grDevices::adjustcolor(EST_COLS[ESTIMATORS], alpha.f = 0.20), border = EST_COLS[ESTIMATORS])
  }
  mtext("Robustness jittered distributions", outer = TRUE, cex = 1.1, font = 2)
}

make_example_fits <- function(sc, criteria = "BIC") {
  dat <- generate_dataset(n = sc$n, p = sc$p, delta = sc$delta, rho = sc$rho, f_id = sc$f_id, err_case = sc$err_case)
  out <- list(dat = dat)
  
  sel_ks <- Selection_bw_k_safe(dat$x, dat$z, dat$y, BW_SEQ, k_seq = 0, criteria = criteria, R_mode = "I")
  fit_ks <- fit_KS(dat$x, dat$z, dat$y, bw = sel_ks$pairs$bw_hat[1])
  out$KS <- list(fit = fit_ks, bw = sel_ks$pairs$bw_hat[1], k = 0)
  
  sel_rtk <- Selection_bw_k_safe(dat$x, dat$z, dat$y, BW_SEQ, K_SEQ, criteria = criteria, R_mode = "I")
  fit_rtk <- fit_RTK(dat$x, dat$z, dat$y, bw = sel_rtk$pairs$bw_hat[1], k = sel_rtk$pairs$k_hat[1])
  out$RTK <- list(fit = fit_rtk, bw = sel_rtk$pairs$bw_hat[1], k = sel_rtk$pairs$k_hat[1])
  
  sel_grtk <- Selection_bw_k_safe(dat$x, dat$z, dat$y, BW_SEQ, K_SEQ, criteria = criteria, R_mode = "AR1")
  fit_grtk <- fit_GRTK(dat$x, dat$z, dat$y, bw = sel_grtk$pairs$bw_hat[1], k = sel_grtk$pairs$k_hat[1])
  out$GRTK <- list(fit = fit_grtk, bw = sel_grtk$pairs$bw_hat[1], k = sel_grtk$pairs$k_hat[1], rho_hat = fit_grtk$rho_hat)
  
  out
}

plot_example_f_curves <- function(example_obj, file, main_title) {
  if (!SAVE_FIGURES) return(invisible(NULL))
  dat <- example_obj$dat
  z <- dat$z
  true_f <- dat$f
  ord <- order(z)
  z <- z[ord]
  true_f <- true_f[ord]
  
  open_device(file, width = 11, height = 4)
  oldpar <- par(no.readonly = TRUE)
  on.exit({par(oldpar); dev.off()}, add = TRUE)
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  
  for (est in ESTIMATORS) {
    fit <- example_obj[[est]]$fit
    partial_res <- dat$y - as.vector(dat$x %*% fit$beta_hat)
    fhat <- fit$fhat[ord]
    plot(z, partial_res[ord], pch = 16, cex = 0.65,
         col = grDevices::adjustcolor("grey40", alpha.f = 0.35),
         xlab = "z", ylab = "Partial residual / f(z)", main = est)
    lines(z, true_f, lwd = 2, lty = 2, col = "black")
    lines(z, fhat, lwd = 2.5, col = EST_COLS[est])
    legend("topright", legend = c("Partial residuals", "True f", "Estimated f"),
           pch = c(16, NA, NA), lty = c(NA, 2, 1), lwd = c(NA, 2, 2.5),
           col = c(grDevices::adjustcolor("grey40", alpha.f = 0.35), "black", EST_COLS[est]),
           cex = 0.75, bty = "n")
  }
  mtext(main_title, outer = TRUE, cex = 1.1, font = 2)
}

#run simulations--------------

cat("Reviewer extension simulation script (v4)\n")
cat("B =", B, "| criteria =", paste(CRITERIA, collapse = ", "), "| cores =", ifelse(USE_PARALLEL, N_CORES, 1), "| figure device =", FIG_DEVICE, "\n")

highdim_raw <- run_grid(HIGH_DIM_GRID, label = "High-dimensional")
robust_raw <- run_grid(ROBUST_GRID, label = "Robustness")

write.csv(highdim_raw, file.path(TABLE_DIR, "raw_highdim_results.csv"), row.names = FALSE)
write.csv(robust_raw, file.path(TABLE_DIR, "raw_robustness_results.csv"), row.names = FALSE)

highdim_sum <- summarise_results(highdim_raw)
robust_sum <- summarise_results(robust_raw)

write.csv(highdim_sum, file.path(TABLE_DIR, "summary_highdim_all.csv"), row.names = FALSE)
write.csv(robust_sum, file.path(TABLE_DIR, "summary_robustness_all.csv"), row.names = FALSE)

hd_table <- subset(highdim_sum, criterion == "BIC",
                   select = c(n, p, delta, rho, estimator, reps, fail,
                              bw_mean, k_mean, rho_hat_mean,
                              beta_smsd_per_coef_mean, beta_smsd_total_mean,
                              active_beta_mse_mean, inactive_beta_mse_mean,
                              f_mse_mean, pred_mse_mean))
hd_table <- hd_table[order(hd_table$n, hd_table$p, hd_table$delta, match(hd_table$estimator, ESTIMATORS)), ]

rob_table <- subset(robust_sum, criterion == "BIC",
                    select = c(n, p, delta, rho, err_case, estimator, reps, fail,
                               bw_mean, k_mean, rho_hat_mean,
                               beta_smsd_per_coef_mean, beta_smsd_total_mean,
                               active_beta_mse_mean, inactive_beta_mse_mean,
                               f_mse_mean, pred_mse_mean))
rob_table <- rob_table[order(rob_table$n, rob_table$delta, rob_table$err_case, match(rob_table$estimator, ESTIMATORS)), ]

write.csv(hd_table, file.path(TABLE_DIR, "Table_R3_2_highdim_BIC.csv"), row.names = FALSE)
write.csv(rob_table, file.path(TABLE_DIR, "Table_R3_3_robustness_BIC.csv"), row.names = FALSE)
write_highdim_latex_table(hd_table, file.path(TABLE_DIR, "Table_R3_2_highdim_BIC.tex"))
write_robust_latex_table(rob_table, file.path(TABLE_DIR, "Table_R3_3_robustness_BIC.tex"))

hd_gain <- grtk_gain_table(subset(highdim_sum, criterion == "BIC"))
rb_gain <- grtk_gain_table(subset(robust_sum, criterion == "BIC"))
write.csv(hd_gain, file.path(TABLE_DIR, "Table_R3_2_highdim_GRTK_gain.csv"), row.names = FALSE)
write.csv(rb_gain, file.path(TABLE_DIR, "Table_R3_3_robustness_GRTK_gain.csv"), row.names = FALSE)

# Improved figures, PNG or JPEG only.
plot_hd_boxpanels(highdim_raw, "pred_mse",
                  fig_file("Figure_R3_2_highdim_pred_mse_boxpanels"),
                  ylab = "Prediction MSE")
plot_hd_strip_panels(highdim_raw, "pred_mse",
                     fig_file("Figure_R3_2_highdim_pred_mse_strip"),
                     ylab = "Prediction MSE")
plot_hd_boxpanels(highdim_raw, "beta_smsd_per_coef",
                  fig_file("Figure_R3_2_highdim_beta_smsd_boxpanels"),
                  ylab = "SMSD / p")
plot_summary_profile_hd(subset(highdim_sum, criterion == "BIC"), "pred_mse",
                        fig_file("Figure_R3_2_highdim_pred_profile"),
                        ylab = "Mean prediction MSE")

plot_robust_boxpanels(robust_raw, "pred_mse",
                      fig_file("Figure_R3_3_robustness_pred_mse_boxpanels"),
                      ylab = "Prediction MSE")
plot_robust_strip_panels(robust_raw, "f_mse",
                         fig_file("Figure_R3_3_robustness_f_mse_strip"),
                         ylab = expression(MSE(hat(f))))
plot_robust_boxpanels(robust_raw, "beta_smsd_per_coef",
                      fig_file("Figure_R3_3_robustness_beta_smsd_boxpanels"),
                      ylab = "SMSD / p")

# Illustrative smooth-function plots.
example_hd <- make_example_fits(list(n = 300, p = 100, delta = 0.99, rho = 0.60, f_id = 1, err_case = "AR1_Gauss"))
plot_example_f_curves(example_hd,
                      fig_file("Figure_R3_2_highdim_example_fits"),
                      main_title = "Illustrative fitted nonparametric components in a high-dimensional setting")

example_rb <- make_example_fits(list(n = 300, p = 6, delta = 0.99, rho = 0.60, f_id = 1, err_case = "ARMA11_Gauss"))
plot_example_f_curves(example_rb,
                      fig_file("Figure_R3_3_robustness_example_fits"),
                      main_title = "Illustrative fitted nonparametric components under ARMA(1,1) errors")

