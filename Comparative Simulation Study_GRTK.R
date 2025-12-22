###########---------------------
# simulation_study.R (by E. Yilmaz)
# Purpose:
#   - Full simulation study for GRTK paper with comprehensive outputs
#   - Metrics: Bias, Variance, SMSD, MSE(f), Prediction MSE
#   - Figures: Bias grids, fitted curves, 3D surfaces, boxplots, comparison plots
#   - Tables: Tuning parameters, performance metrics, robustness results
##------------------------------

options(stringsAsFactors = FALSE)
set.seed(20251218)

source("tuning_selection.R")

#-Parallel setup -------------------------------------------

USE_PARALLEL <- TRUE
N_CORES <- parallel::detectCores() - 1
if (N_CORES < 1) N_CORES <- 1

if (USE_PARALLEL) {
  library(parallel)
  cat("Using parallel processing with", N_CORES, "cores\n")
}

# - Output directory -----------------------------------------

OUT_DIR <- "outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------ Simulation grid ------------------------------------------

SIM_N     <- c(40, 150, 300)
SIM_DELTA <- c(0.80, 0.90, 0.99)
SIM_RHO   <- c(0.00, 0.60, 0.90)

P <- 6

# - replications
B <- 50  # 1000

# Tuning grids (Be careful when deciding them :)--------

BW_SEQ <- c(seq(0.001, 0.01, by = 0.005),
            seq(0.015, 0.1, by = 0.01),
            seq(0.11, 0.25, by = 0.03))

K_SEQ  <- c(0, 0.001, 0.005, 0.01, 0.025, 0.05, 0.10, 0.25, 
            0.50, 1, 2, 5, 10, 20)

CRITERIA <- c("GCV", "AICc", "BIC", "RECP")

# Truth: f(z) and beta -------------------------------------

f_true <- function(z, which = 1) {
  if (which == 1) return(2 * sin(2 * pi * z))
  if (which == 2) return(4 * (z - 0.5)^2 - 1)
  stop("f(which) must be 1 or 2")
}

beta_true <- c(1.0, 0.8, 0.5, 0, 0, 0)

#  Multicollinear X generation ------------------------------

toeplitz_sigma <- function(p, delta) {
  outer(1:p, 1:p, function(i, j) delta^abs(i - j))
}

gen_X_multicollinear <- function(n, p, delta) {
  Sigma <- toeplitz_sigma(p, delta)
  U <- chol(Sigma)
  Z <- matrix(rnorm(n * p), n, p)
  X <- Z %*% U
  X
}

#  Error generators -----------------------------------------

gen_errors_AR1 <- function(n, rho, innov = c("gaussian", "t5"), burn = 300, sigma = 1) {
  innov <- match.arg(innov)
  N <- n + burn
  
  u <- switch(innov,
              gaussian = rnorm(N, 0, sigma),
              t5 = {
                df <- 5
                scale <- sqrt((df - 2) / df) * sigma
                rt(N, df = df) * scale
              })
  
  e <- as.numeric(stats::filter(u, filter = rho, method = "recursive", init = 0))
  e[(burn + 1):N]
}

gen_errors_AR2 <- function(n, rho1 = 0.50, rho2 = 0.25,
                           innov = c("gaussian", "t5"), burn = 500, sigma = 1) {
  innov <- match.arg(innov)
  N <- n + burn
  
  u <- switch(innov,
              gaussian = rnorm(N, 0, sigma),
              t5 = {
                df <- 5
                scale <- sqrt((df - 2) / df) * sigma
                rt(N, df = df) * scale
              })
  
  e <- as.numeric(stats::filter(u, filter = c(rho1, rho2), method = "recursive", init = c(0, 0)))
  e[(burn + 1):N]
}

#  Dataset generator ----------------------------------------

generate_dataset <- function(n, p, delta, rho, f_id = 1,
                             err_case = c("AR1_Gauss", "AR1_t5", "AR2_Gauss")) {
  
  err_case <- match.arg(err_case)
  
  z <- (1:n) / n
  x <- gen_X_multicollinear(n, p, delta)
  f <- f_true(z, which = f_id)
  
  eps <- switch(err_case,
                AR1_Gauss = gen_errors_AR1(n, rho, innov = "gaussian"),
                AR1_t5    = gen_errors_AR1(n, rho, innov = "t5"),
                AR2_Gauss = gen_errors_AR2(n, rho1 = 0.50, rho2 = 0.25, innov = "gaussian"))
  
  y <- as.vector(x %*% beta_true + f + eps)
  
  list(y = y, x = x, z = z, f = f, beta = beta_true, eps = eps, err_case = err_case)
}

# ------------------ Fit wrappers ---------------------------------------------

fit_KS <- function(x, z, y, bw) {
  n <- length(y)
  fit_plm_once(x, z, y, bw = bw, k = 0, Rinv = diag(n))
}

fit_RTK <- function(x, z, y, bw, k) {
  n <- length(y)
  fit_plm_once(x, z, y, bw = bw, k = k, Rinv = diag(n))
}

fit_GRTK <- function(x, z, y, bw, k, lag = 1, bw_pilot = 0.15, k_pilot = 0.10) {
  n <- length(y)
  
  pilot <- fit_plm_once(x, z, y, bw = bw_pilot, k = k_pilot, Rinv = diag(n))
  rho_hat <- estimate_rho_ar1(y - pilot$yhat, lag = lag)
  
  R <- AR1_corr_matrix(n, rho_hat)
  Rinv <- safe_chol_inv(R)
  
  fit <- fit_plm_once(x, z, y, bw = bw, k = k, Rinv = Rinv)
  fit$rho_hat <- rho_hat
  fit
}

# ------------------ Tuning selection wrapper ---------------------------------

select_pair <- function(x, z, y, bw_seq, k_seq, criterion, R_mode = c("AR1", "I"), lag = 1, out_dir = "outputs") {
  R_mode <- match.arg(R_mode)
  sel <- Selection_bw_k(
    x = x, z = z, y = y,
    bw_seq = bw_seq, k_seq = k_seq,
    lag = lag,
    criteria = criterion,
    R_mode = R_mode,
    make_plots = FALSE,
    make_contours = FALSE,
    save_outputs = FALSE,
    out_dir = out_dir,
    prefix = "tmp"
  )
  sel$pairs[1, c("bw_hat", "k_hat")]
}

# metrics ------------------------------------

compute_metrics <- function(beta_hat, f_hat, yhat, beta_true, f_true_vec, x) {
  mu_true <- as.vector(x %*% beta_true + f_true_vec)
  n <- length(yhat)
  p <- length(beta_true)
  
  # Individual beta metrics
  beta_bias <- beta_hat - beta_true
  beta_bias_sq <- beta_bias^2
  
  # Aggregated metrics
  total_bias_sq <- sum(beta_bias_sq)         # Sum of squared biases
  total_bias <- sqrt(total_bias_sq)           # L2 norm of bias
  
  # f metrics
  f_ise <- mean((f_hat - f_true_vec)^2)       # Integrated squared error
  f_bias <- mean(f_hat - f_true_vec)          # Mean bias of f
  
  # Prediction metrics
  pred_mse <- mean((yhat - mu_true)^2)
  
  list(
    beta_hat = beta_hat,
    beta_bias = beta_bias,
    beta_bias_sq = beta_bias_sq,
    total_bias_sq = total_bias_sq,
    total_bias = total_bias,
    f_ise = f_ise,
    f_bias = f_bias,
    pred_mse = pred_mse,
    f_hat = f_hat
  )
}

# (1) TUNING DEMO
cat("\n=== Running tuning demonstration ===\n")

demo_data <- generate_dataset(n = 300, p = P, delta = 0.90, rho = 0.60, f_id = 1, err_case = "AR1_Gauss")

demo_prefix <- "tuning_demo_n300_delta090_rho060_f1_AR1Gauss"

demo_sel <- Selection_bw_k(
  x = demo_data$x, z = demo_data$z, y = demo_data$y,
  bw_seq = BW_SEQ, k_seq = K_SEQ,
  lag = 1,
  criteria = CRITERIA,
  R_mode = "AR1",
  make_plots = TRUE,
  make_contours = TRUE,
  save_outputs = TRUE,
  out_dir = OUT_DIR,
  prefix = demo_prefix
)

cat("Tuning demo complete. Selected pairs:\n")
print(demo_sel$pairs)

# (2) FULL SIMULATION 

scenarios <- expand.grid(
  n = SIM_N,
  delta = SIM_DELTA,
  rho = SIM_RHO,
  f_id = c(1, 2),
  err_case = c("AR1_Gauss", "AR1_t5"),
  stringsAsFactors = FALSE
)

cat("\n=== Running full simulation study ===\n")
cat("Scenarios:", nrow(scenarios), "| Replications:", B, "\n")

# Function to run one replication
run_one_replication <- function(sc, bw_seq, k_seq, criteria, p, beta_true) {
  
  dat <- generate_dataset(
    n = sc$n, p = p, delta = sc$delta, rho = sc$rho,
    f_id = sc$f_id, err_case = sc$err_case
  )
  
  results_list <- vector("list", length(criteria) * 3)
  idx <- 1L
  
  bw_pilot <- median(bw_seq)
  k_pilot <- median(k_seq[k_seq > 0])
  
  for (crit in criteria) {
    
    # ---- KS ----
    bw_ks <- select_pair(dat$x, dat$z, dat$y, bw_seq, k_seq = 0, criterion = crit, R_mode = "I")
    fitKS <- fit_KS(dat$x, dat$z, dat$y, bw = bw_ks$bw_hat)
    mKS <- compute_metrics(fitKS$beta_hat, fitKS$fhat, fitKS$yhat, dat$beta, dat$f, dat$x)
    
    results_list[[idx]] <- data.frame(
      n = sc$n, delta = sc$delta, rho = sc$rho, f_id = sc$f_id, err_case = sc$err_case,
      criterion = crit, estimator = "KS",
      bw_hat = bw_ks$bw_hat, k_hat = 0, rho_hat = NA_real_,
      # Individual beta biases
      bias_beta1 = mKS$beta_bias[1], bias_beta2 = mKS$beta_bias[2], bias_beta3 = mKS$beta_bias[3],
      bias_beta4 = mKS$beta_bias[4], bias_beta5 = mKS$beta_bias[5], bias_beta6 = mKS$beta_bias[6],
      # Aggregated metrics
      total_bias_sq = mKS$total_bias_sq,
      f_ise = mKS$f_ise,
      pred_mse = mKS$pred_mse,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
    
    # ---- RTK ----
    bk_rtk <- select_pair(dat$x, dat$z, dat$y, bw_seq, k_seq, criterion = crit, R_mode = "I")
    fitRTK <- fit_RTK(dat$x, dat$z, dat$y, bw = bk_rtk$bw_hat, k = bk_rtk$k_hat)
    mRTK <- compute_metrics(fitRTK$beta_hat, fitRTK$fhat, fitRTK$yhat, dat$beta, dat$f, dat$x)
    
    results_list[[idx]] <- data.frame(
      n = sc$n, delta = sc$delta, rho = sc$rho, f_id = sc$f_id, err_case = sc$err_case,
      criterion = crit, estimator = "RTK",
      bw_hat = bk_rtk$bw_hat, k_hat = bk_rtk$k_hat, rho_hat = NA_real_,
      bias_beta1 = mRTK$beta_bias[1], bias_beta2 = mRTK$beta_bias[2], bias_beta3 = mRTK$beta_bias[3],
      bias_beta4 = mRTK$beta_bias[4], bias_beta5 = mRTK$beta_bias[5], bias_beta6 = mRTK$beta_bias[6],
      total_bias_sq = mRTK$total_bias_sq,
      f_ise = mRTK$f_ise,
      pred_mse = mRTK$pred_mse,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
    
    # ---- GRTK ----
    bk_grtk <- select_pair(dat$x, dat$z, dat$y, bw_seq, k_seq, criterion = crit, R_mode = "AR1")
    fitGRTK <- fit_GRTK(dat$x, dat$z, dat$y, bw = bk_grtk$bw_hat, k = bk_grtk$k_hat,
                        lag = 1, bw_pilot = bw_pilot, k_pilot = k_pilot)
    mGRTK <- compute_metrics(fitGRTK$beta_hat, fitGRTK$fhat, fitGRTK$yhat, dat$beta, dat$f, dat$x)
    
    results_list[[idx]] <- data.frame(
      n = sc$n, delta = sc$delta, rho = sc$rho, f_id = sc$f_id, err_case = sc$err_case,
      criterion = crit, estimator = "GRTK",
      bw_hat = bk_grtk$bw_hat, k_hat = bk_grtk$k_hat, rho_hat = fitGRTK$rho_hat,
      bias_beta1 = mGRTK$beta_bias[1], bias_beta2 = mGRTK$beta_bias[2], bias_beta3 = mGRTK$beta_bias[3],
      bias_beta4 = mGRTK$beta_bias[4], bias_beta5 = mGRTK$beta_bias[5], bias_beta6 = mGRTK$beta_bias[6],
      total_bias_sq = mGRTK$total_bias_sq,
      f_ise = mGRTK$f_ise,
      pred_mse = mGRTK$pred_mse,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
  
  do.call(rbind, results_list)
}

# Build task list
tasks <- vector("list", nrow(scenarios) * B)
task_idx <- 1L
for (s in seq_len(nrow(scenarios))) {
  for (b in seq_len(B)) {
    tasks[[task_idx]] <- list(scenario = scenarios[s, ], rep_id = b, scenario_id = s)
    task_idx <- task_idx + 1L
  }
}

total_tasks <- length(tasks)
cat("Total tasks:", total_tasks, "\n")

start_time <- Sys.time()

if (USE_PARALLEL && N_CORES > 1) {
  
  cl <- makeCluster(N_CORES)
  
  clusterExport(cl, c("BW_SEQ", "K_SEQ", "CRITERIA", "P", "beta_true",
                      "generate_dataset", "f_true", "toeplitz_sigma", "gen_X_multicollinear",
                      "gen_errors_AR1", "gen_errors_AR2",
                      "fit_KS", "fit_RTK", "fit_GRTK", "fit_plm_once",
                      "select_pair", "Selection_bw_k", "compute_metrics",
                      "NW_smoother_matrix", "gaussian_kernel", "AR1_corr_matrix",
                      "safe_chol_inv", "estimate_rho_ar1", "df_fast", "H_matrix",
                      "run_one_replication", "OUT_DIR"))
  
  clusterSetRNGStream(cl, 20251218)
  
  cat("Running parallel simulation...\n")
  
  BATCH_SIZE <- max(10, N_CORES * 2)
  n_batches <- ceiling(total_tasks / BATCH_SIZE)
  results_list <- vector("list", total_tasks)
  
  for (batch in seq_len(n_batches)) {
    start_idx <- (batch - 1) * BATCH_SIZE + 1
    end_idx <- min(batch * BATCH_SIZE, total_tasks)
    batch_tasks <- tasks[start_idx:end_idx]
    
    batch_results <- parLapply(cl, batch_tasks, function(task) {
      run_one_replication(task$scenario, BW_SEQ, K_SEQ, CRITERIA, P, beta_true)
    })
    
    results_list[start_idx:end_idx] <- batch_results
    
    completed <- end_idx
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    rate <- completed / elapsed
    eta <- (total_tasks - completed) / rate
    
    cat(sprintf("\r  Progress: %d/%d tasks (%.1f%%) | Elapsed: %.0fs | ETA: %.0fs     ",
                completed, total_tasks, 100 * completed / total_tasks, elapsed, eta))
    flush.console()
  }
  cat("\n")
  
  stopCluster(cl)
  
} else {
  
  cat("Running sequential simulation...\n")
  results_list <- vector("list", total_tasks)
  
  for (i in seq_along(tasks)) {
    if (i %% 10 == 0 || i == 1) {
      cat(sprintf("\rTask %d/%d (%.1f%%)", i, total_tasks, 100 * i / total_tasks))
      flush.console()
    }
    results_list[[i]] <- run_one_replication(tasks[[i]]$scenario, BW_SEQ, K_SEQ, CRITERIA, P, beta_true)
  }
  cat("\n")
}

end_time <- Sys.time()
cat("\nSimulation completed in:", format(end_time - start_time), "\n")

res_df <- do.call(rbind, results_list)
rownames(res_df) <- NULL

write.csv(res_df, file.path(OUT_DIR, "simulation_raw_results.csv"), row.names = FALSE)

# (3) COMPUTE  STATISTICS - SMSD CALCULATION

cat("\n=== Computing summary statistics (including SMSD) ===\n")

# For SMSD: SMSD = Bias^2 + Variance (in case of)

compute_smsd_table <- function(df, groupvars) {
  agg <- aggregate(
    cbind(bias_beta1, bias_beta2, bias_beta3, total_bias_sq, f_ise, pred_mse, bw_hat, k_hat) ~ 
      n + delta + rho + f_id + err_case + criterion + estimator,
    data = df,
    FUN = function(x) {
      c(mean = mean(x, na.rm = TRUE), 
        var = var(x, na.rm = TRUE),
        sd = sd(x, na.rm = TRUE))
    }
  )
  
  # Flatten and compute SMSD
  result <- agg[, 1:7]
  

  for (b in 1:3) {
    col <- paste0("bias_beta", b)
    result[[paste0("bias_beta", b, "_mean")]] <- agg[[col]][, "mean"]
    result[[paste0("bias_beta", b, "_var")]] <- agg[[col]][, "var"]
    result[[paste0("smsd_beta", b)]] <- agg[[col]][, "mean"]^2 + agg[[col]][, "var"]
  }
  
  # Total SMSD (sum of individual SMSDs)
  result$smsd_total <- result$smsd_beta1 + result$smsd_beta2 + result$smsd_beta3
  
  # Alternative: use total_bias_sq directly
  result$total_bias_sq_mean <- agg$total_bias_sq[, "mean"]
  result$total_bias_sq_var <- agg$total_bias_sq[, "var"]
  
  # f and prediction metrics
  result$f_ise_mean <- agg$f_ise[, "mean"]
  result$f_ise_sd <- agg$f_ise[, "sd"]
  result$pred_mse_mean <- agg$pred_mse[, "mean"]
  result$pred_mse_sd <- agg$pred_mse[, "sd"]
  
  # Tuning parameters
  result$bw_mean <- agg$bw_hat[, "mean"]
  result$k_mean <- agg$k_hat[, "mean"]
  
  result
}

summary_smsd <- compute_smsd_table(res_df)
write.csv(summary_smsd, file.path(OUT_DIR, "simulation_summary_with_SMSD.csv"), row.names = FALSE)

# (4) CREATE TABLES

cat("\n=== Creating paper-style tables ===\n")

# TABLE 1: Optimal (lambda, k) pairs by criteria (like paper Table 1)
table1 <- subset(summary_smsd, 
                 f_id == 1 & err_case == "AR1_Gauss" & estimator == "GRTK",
                 select = c(n, delta, rho, criterion, bw_mean, k_mean))
table1 <- table1[order(table1$n, table1$delta, table1$rho, table1$criterion), ]
write.csv(table1, file.path(OUT_DIR, "Table1_tuning_pairs.csv"), row.names = FALSE)

# TABLE 2: Bias, Variance, SMSD for beta (like paper Table 2)
table2 <- subset(summary_smsd, 
                 f_id == 1 & err_case == "AR1_Gauss" & estimator == "GRTK" & rho == 0.60,
                 select = c(n, delta, criterion, 
                            bias_beta1_mean, bias_beta1_var, smsd_beta1,
                            bias_beta2_mean, bias_beta2_var, smsd_beta2,
                            smsd_total))
write.csv(table2, file.path(OUT_DIR, "Table2_beta_bias_var_SMSD.csv"), row.names = FALSE)

# TABLE 3: MSE for nonparametric component (like paper Table 3)
table3 <- subset(summary_smsd,
                 f_id == 1 & err_case == "AR1_Gauss" & rho == 0.60,
                 select = c(n, delta, criterion, estimator, f_ise_mean, f_ise_sd))
table3_wide <- reshape(table3, 
                       idvar = c("n", "delta", "criterion"),
                       timevar = "estimator",
                       direction = "wide")
write.csv(table3_wide, file.path(OUT_DIR, "Table3_f_MSE.csv"), row.names = FALSE)

# TABLE 4: Comprehensive comparison across estimators
table4 <- subset(summary_smsd,
                 f_id == 1 & err_case == "AR1_Gauss" & criterion == "BIC",
                 select = c(n, delta, rho, estimator, smsd_total, f_ise_mean, pred_mse_mean))
write.csv(table4, file.path(OUT_DIR, "Table4_estimator_comparison.csv"), row.names = FALSE)

cat("Tables 1-4 saved.\n")

# TABLE 5: Comprehensive Criteria Comparison (GRTK estimator)
# Shows performance of each criterion across all metrics
table5_data <- subset(summary_smsd,
                      f_id == 1 & err_case == "AR1_Gauss" & estimator == "GRTK")
table5 <- aggregate(
  cbind(smsd_total, f_ise_mean, pred_mse_mean, bw_mean, k_mean) ~ 
    n + delta + rho + criterion,
  data = table5_data,
  FUN = mean
)
table5 <- table5[order(table5$n, table5$delta, table5$rho, table5$criterion), ]
write.csv(table5, file.path(OUT_DIR, "Table5_criteria_comparison.csv"), row.names = FALSE)

# TABLE 6: Best criterion by configuration (which criterion wins?)
# For each (n, delta, rho), find which criterion gives lowest SMSD
find_best_criterion <- function(df) {
  configs <- unique(df[, c("n", "delta", "rho")])
  results <- data.frame()
  
  for (i in 1:nrow(configs)) {
    cfg_data <- subset(df, n == configs$n[i] & 
                         abs(delta - configs$delta[i]) < 1e-10 & 
                         abs(rho - configs$rho[i]) < 1e-10)
    
    if (nrow(cfg_data) > 0) {
      best_smsd <- cfg_data[which.min(cfg_data$smsd_total), ]
      best_fise <- cfg_data[which.min(cfg_data$f_ise_mean), ]
      best_pred <- cfg_data[which.min(cfg_data$pred_mse_mean), ]
      
      results <- rbind(results, data.frame(
        n = configs$n[i],
        delta = configs$delta[i],
        rho = configs$rho[i],
        best_SMSD_crit = best_smsd$criterion,
        best_SMSD_val = best_smsd$smsd_total,
        best_fMSE_crit = best_fise$criterion,
        best_fMSE_val = best_fise$f_ise_mean,
        best_Pred_crit = best_pred$criterion,
        best_Pred_val = best_pred$pred_mse_mean
      ))
    }
  }
  results
}

table6 <- find_best_criterion(subset(summary_smsd, 
                                     f_id == 1 & err_case == "AR1_Gauss" & estimator == "GRTK"))
write.csv(table6, file.path(OUT_DIR, "Table6_best_criterion_by_config.csv"), row.names = FALSE)

cat("Tables 5-6 saved.\n")

# (5)  FIGURES

cat("\n=== Generating figures ===\n")

est_levels <- c("KS", "RTK", "GRTK")
est_cols   <- c("black", "red", "blue")
est_pchs   <- c(16, 17, 15)
est_ltys   <- c(1, 2, 4)

crit_cols  <- c("black", "red", "blue", "darkgreen")
crit_names <- c("GCV", "AICc", "BIC", "RECP")
crit_pchs  <- c(16, 17, 15, 18)
crit_ltys  <- c(1, 2, 3, 4)

# -----------------------------------------------------------------------------
# FIGURE 1: Grid of Bias plots for individual beta coefficients
# -----------------------------------------------------------------------------

png(file.path(OUT_DIR, "fig1_beta_bias_grid.png"), width = 1400, height = 1000, res = 120)
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

for (delta_val in SIM_DELTA) {
  for (beta_idx in 1:3) {
    bias_col <- paste0("bias_beta", beta_idx, "_mean")
    
    plot_data <- subset(summary_smsd, 
                        abs(delta - delta_val) < 1e-10 & 
                          f_id == 1 & err_case == "AR1_Gauss" & 
                          criterion == "BIC" & rho == 0.60)
    
    if (nrow(plot_data) > 0) {
      ylim_range <- range(plot_data[[bias_col]], na.rm = TRUE)
      ylim_range <- ylim_range + c(-0.1, 0.1) * max(abs(ylim_range), 0.1)
      
      plot(NA, xlim = range(SIM_N), ylim = ylim_range,
           xlab = "n", ylab = bquote("Bias(" * beta[.(beta_idx)] * ")"),
           main = bquote(delta == .(delta_val)))
      abline(h = 0, lty = 3, col = "gray50")
      
      for (i in seq_along(est_levels)) {
        dd <- plot_data[plot_data$estimator == est_levels[i], ]
        dd <- dd[order(dd$n), ]
        if (nrow(dd) > 0) {
          lines(dd$n, dd[[bias_col]], type = "b", 
                pch = est_pchs[i], col = est_cols[i], lty = est_ltys[i], lwd = 1.5)
        }
      }
      if (delta_val == SIM_DELTA[1] && beta_idx == 1) {
        legend("topright", legend = est_levels, col = est_cols, pch = est_pchs, 
               lty = est_ltys, lwd = 1.5, bty = "n", cex = 0.9)
      }
    }
  }
}
mtext("Bias of Beta Coefficients by Sample Size (BIC, rho=0.60, f1)", outer = TRUE, line = 0.5, cex = 1.2)
dev.off()
cat("Figure 1 saved.\n")

# -----------------------------------------------------------------------------
# FIGURE 2: Fitted curves comparison (KS, RTK, GRTK) - 3x3 grid
# -----------------------------------------------------------------------------

png(file.path(OUT_DIR, "fig2_fitted_curves_grid.png"), width = 1400, height = 1000, res = 120)
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

config_list <- list(
  list(n = 40, delta = 0.80, rho = 0.00),
  list(n = 40, delta = 0.90, rho = 0.60),
  list(n = 40, delta = 0.99, rho = 0.90),
  list(n = 150, delta = 0.80, rho = 0.00),
  list(n = 150, delta = 0.90, rho = 0.60),
  list(n = 150, delta = 0.99, rho = 0.90),
  list(n = 300, delta = 0.80, rho = 0.00),
  list(n = 300, delta = 0.90, rho = 0.60),
  list(n = 300, delta = 0.99, rho = 0.90)
)

set.seed(12345)
for (cfg in config_list) {
  dat <- generate_dataset(n = cfg$n, p = P, delta = cfg$delta, rho = cfg$rho, 
                          f_id = 1, err_case = "AR1_Gauss")
  
  bw_ks <- select_pair(dat$x, dat$z, dat$y, BW_SEQ, k_seq = 0, criterion = "BIC", R_mode = "I")
  bk_rtk <- select_pair(dat$x, dat$z, dat$y, BW_SEQ, K_SEQ, criterion = "BIC", R_mode = "I")
  bk_grtk <- select_pair(dat$x, dat$z, dat$y, BW_SEQ, K_SEQ, criterion = "BIC", R_mode = "AR1")
  
  fitKS <- fit_KS(dat$x, dat$z, dat$y, bw = bw_ks$bw_hat)
  fitRTK <- fit_RTK(dat$x, dat$z, dat$y, bw = bk_rtk$bw_hat, k = bk_rtk$k_hat)
  fitGRTK <- fit_GRTK(dat$x, dat$z, dat$y, bw = bk_grtk$bw_hat, k = bk_grtk$k_hat,
                      bw_pilot = median(BW_SEQ), k_pilot = median(K_SEQ[K_SEQ > 0]))
  
  ylim_range <- range(c(dat$f, fitKS$fhat, fitRTK$fhat, fitGRTK$fhat), na.rm = TRUE)
  
  plot(dat$z, dat$f, type = "l", lwd = 2.5, col = "gray40",
       xlab = "z", ylab = "f(z)", 
       main = bquote(n == .(cfg$n) ~ ", " ~ delta == .(cfg$delta) ~ ", " ~ rho == .(cfg$rho)),
       ylim = ylim_range)
  lines(dat$z, fitKS$fhat, col = est_cols[1], lty = est_ltys[1], lwd = 1.5)
  lines(dat$z, fitRTK$fhat, col = est_cols[2], lty = est_ltys[2], lwd = 1.5)
  lines(dat$z, fitGRTK$fhat, col = est_cols[3], lty = est_ltys[3], lwd = 1.5)
}
mtext("Fitted Curves: True f(z) vs Estimates (BIC)", outer = TRUE, line = 0.5, cex = 1.2)

# Add legend to last panel
legend("bottomright", legend = c("True", est_levels), 
       col = c("gray40", est_cols), lty = c(1, est_ltys), lwd = c(2.5, rep(1.5, 3)), bty = "n")
dev.off()
cat("Figure 2 saved.\n")

# -----------------------------------------------------------------------------
# FIGURE 3: SMSD vs rho (effect of autocorrelation)
# -----------------------------------------------------------------------------

png(file.path(OUT_DIR, "fig3_SMSD_vs_rho.png"), width = 1200, height = 800, res = 120)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

for (n_val in SIM_N) {
  for (delta_val in c(0.90, 0.99)) {
    plot_data <- subset(summary_smsd, 
                        n == n_val & abs(delta - delta_val) < 1e-10 & 
                          f_id == 1 & err_case == "AR1_Gauss" & criterion == "BIC")
    
    if (nrow(plot_data) > 0) {
      ylim_range <- range(plot_data$smsd_total, na.rm = TRUE)
      ylim_range <- ylim_range + c(-0.05, 0.1) * diff(ylim_range)
      
      plot(NA, xlim = range(SIM_RHO), ylim = ylim_range,
           xlab = expression(rho), ylab = "SMSD",
           main = bquote(n == .(n_val) ~ ", " ~ delta == .(delta_val)))
      
      for (i in seq_along(est_levels)) {
        dd <- plot_data[plot_data$estimator == est_levels[i], ]
        dd <- dd[order(dd$rho), ]
        if (nrow(dd) > 0) {
          lines(dd$rho, dd$smsd_total, type = "b",
                pch = est_pchs[i], col = est_cols[i], lty = est_ltys[i], lwd = 1.5)
        }
      }
      legend("topleft", legend = est_levels, col = est_cols, pch = est_pchs, 
             lty = est_ltys, lwd = 1.5, bty = "n", cex = 0.8)
    }
  }
}
mtext("SMSD vs Autocorrelation Level (BIC, f1, AR1-Gaussian)", outer = TRUE, line = 0.5, cex = 1.2)
dev.off()
cat("Figure 3 saved.\n")

# -----------------------------------------------------------------------------
# FIGURE 4: f ISE vs delta (effect of multicollinearity) 
# -----------------------------------------------------------------------------

png(file.path(OUT_DIR, "fig4_fISE_vs_delta.png"), width = 1200, height = 800, res = 120)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

for (n_val in SIM_N) {
  for (rho_val in c(0.00, 0.90)) {
    plot_data <- subset(summary_smsd,
                        n == n_val & abs(rho - rho_val) < 1e-10 &
                          f_id == 1 & err_case == "AR1_Gauss" & criterion == "BIC")
    
    if (nrow(plot_data) > 0) {
      ylim_range <- range(plot_data$f_ise_mean, na.rm = TRUE)
      ylim_range <- ylim_range + c(-0.05, 0.1) * diff(ylim_range)
      
      plot(NA, xlim = range(SIM_DELTA), ylim = ylim_range,
           xlab = expression(delta), ylab = "MSE(f)",
           main = bquote(n == .(n_val) ~ ", " ~ rho == .(rho_val)))
      
      for (i in seq_along(est_levels)) {
        dd <- plot_data[plot_data$estimator == est_levels[i], ]
        dd <- dd[order(dd$delta), ]
        if (nrow(dd) > 0) {
          lines(dd$delta, dd$f_ise_mean, type = "b",
                pch = est_pchs[i], col = est_cols[i], lty = est_ltys[i], lwd = 1.5)
        }
      }
      legend("topleft", legend = est_levels, col = est_cols, pch = est_pchs,
             lty = est_ltys, lwd = 1.5, bty = "n", cex = 0.8)
    }
  }
}
mtext("MSE of Nonparametric Component vs Multicollinearity (BIC, f1)", outer = TRUE, line = 0.5, cex = 1.2)
dev.off()
cat("Figure 4 saved.\n")

# -----------------------------------------------------------------------------
# FIGURE 5: Criteria comparison bar plots
# -----------------------------------------------------------------------------

png(file.path(OUT_DIR, "fig5_criteria_comparison.png"), width = 1400, height = 800, res = 120)
par(mfrow = c(2, 3), mar = c(5, 4, 3, 1), oma = c(0, 0, 2, 0))

for (n_val in SIM_N) {
  for (delta_val in c(0.90, 0.99)) {
    plot_data <- subset(summary_smsd,
                        n == n_val & abs(delta - delta_val) < 1e-10 &
                          abs(rho - 0.60) < 1e-10 & f_id == 1 & 
                          err_case == "AR1_Gauss" & estimator == "GRTK")
    
    if (nrow(plot_data) >= 4) {
      plot_data <- plot_data[order(match(plot_data$criterion, crit_names)), ]
      
      barplot(plot_data$smsd_total, names.arg = plot_data$criterion,
              col = crit_cols, ylim = c(0, max(plot_data$smsd_total) * 1.2),
              ylab = "SMSD",
              main = bquote(n == .(n_val) ~ ", " ~ delta == .(delta_val)))
    }
  }
}
mtext("SMSD by Selection Criterion (GRTK, rho=0.60, f1)", outer = TRUE, line = 0.5, cex = 1.2)
dev.off()
cat("Figure 5 saved.\n")

# -----------------------------------------------------------------------------
# FIGURE 6: Boxplots of beta estimates (like paper Figure 3)
# -----------------------------------------------------------------------------

png(file.path(OUT_DIR, "fig6_beta_boxplots.png"), width = 1400, height = 1000, res = 120)
par(mfrow = c(2, 3), mar = c(5, 4, 3, 1), oma = c(0, 0, 2, 0))

configs_box <- list(
  list(n = 40, delta = 0.90), list(n = 150, delta = 0.90), list(n = 300, delta = 0.90),
  list(n = 40, delta = 0.99), list(n = 150, delta = 0.99), list(n = 300, delta = 0.99)
)

for (cfg in configs_box) {
  box_data <- subset(res_df, 
                     n == cfg$n & abs(delta - cfg$delta) < 1e-10 &
                       abs(rho - 0.60) < 1e-10 & f_id == 1 & 
                       err_case == "AR1_Gauss" & criterion == "BIC")
  
  if (nrow(box_data) > 0) {
    # Create boxplot data for beta1 and beta2 for each estimator
    bp_list <- list()
    bp_cols <- c()
    bp_names <- c()
    
    for (est in est_levels) {
      est_data <- box_data[box_data$estimator == est, ]
      if (nrow(est_data) > 0) {
        bp_list[[paste0(est, "_b1")]] <- est_data$bias_beta1
        bp_list[[paste0(est, "_b2")]] <- est_data$bias_beta2
        bp_names <- c(bp_names, paste0(substr(est, 1, 1), "1"), paste0(substr(est, 1, 1), "2"))
        bp_cols <- c(bp_cols, est_cols[which(est_levels == est)], est_cols[which(est_levels == est)])
      }
    }
    
    if (length(bp_list) > 0) {
      boxplot(bp_list, names = bp_names,
              col = bp_cols,
              ylab = "Bias", 
              main = bquote(n == .(cfg$n) ~ ", " ~ delta == .(cfg$delta)),
              las = 2, cex.axis = 0.9)
      abline(h = 0, lty = 2, col = "gray50", lwd = 1.5)
    }
  } else {
    # Empty plot with message if no data
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    text(0.5, 0.5, "No data", cex = 1.2)
    title(main = bquote(n == .(cfg$n) ~ ", " ~ delta == .(cfg$delta)))
  }
}
mtext("Boxplots of Beta Bias: K=KS, R=RTK, G=GRTK (BIC, rho=0.60, f1)", outer = TRUE, line = 0.5, cex = 1.2)
dev.off()
cat("Figure 6 saved.\n")

# -----------------------------------------------------------------------------
# FIGURE 7: Robustness - Gaussian vs t5 errors
# -----------------------------------------------------------------------------

png(file.path(OUT_DIR, "fig7_robustness_comparison.png"), width = 1200, height = 500, res = 120)
par(mfrow = c(1, 3), mar = c(5, 4, 3, 1), oma = c(0, 0, 2, 0))

robust_data <- subset(summary_smsd, 
                      n == 300 & abs(delta - 0.90) < 1e-10 & abs(rho - 0.60) < 1e-10 &
                        f_id == 1 & criterion == "BIC" & estimator == "GRTK")

if (nrow(robust_data) >= 2) {
  robust_data <- robust_data[order(robust_data$err_case), ]
  
  barplot(robust_data$smsd_total, names.arg = c("Gaussian", "t5"),
          col = c("steelblue", "salmon"),
          ylab = "SMSD", main = "Beta SMSD")
  
  barplot(robust_data$f_ise_mean, names.arg = c("Gaussian", "t5"),
          col = c("steelblue", "salmon"),
          ylab = "MSE(f)", main = "Nonparametric MSE")
  
  barplot(robust_data$pred_mse_mean, names.arg = c("Gaussian", "t5"),
          col = c("steelblue", "salmon"),
          ylab = "Pred MSE", main = "Prediction MSE")
  
  mtext("Robustness: Gaussian vs Heavy-Tailed Errors (GRTK, BIC, n=300)", outer = TRUE, line = 0.5, cex = 1.2)
} else {
  # Create placeholder plots if data is insufficient
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.5, 0.5, "Insufficient data\nfor robustness comparison", cex = 1.2)
  
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.5, 0.5, "Run with both\nAR1_Gauss and AR1_t5", cex = 1.2)
  
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.5, 0.5, "error cases", cex = 1.2)
  
  mtext("Robustness Comparison (Data Pending)", outer = TRUE, line = 0.5, cex = 1.2)
}
dev.off()
cat("Figure 7 saved.\n")

# -----------------------------------------------------------------------------
# FIGURE 8: SMSD Heatmap (n vs delta)
# -----------------------------------------------------------------------------

png(file.path(OUT_DIR, "fig8_SMSD_heatmap.png"), width = 1000, height = 800, res = 120)

heatmap_data <- subset(summary_smsd,
                       abs(rho - 0.60) < 1e-10 & f_id == 1 & 
                         err_case == "AR1_Gauss" & criterion == "BIC" & estimator == "GRTK")

if (nrow(heatmap_data) >= length(SIM_N) * length(SIM_DELTA)) {
  mat <- matrix(NA, nrow = length(SIM_N), ncol = length(SIM_DELTA))
  rownames(mat) <- as.character(SIM_N)
  colnames(mat) <- as.character(SIM_DELTA)
  
  for (i in seq_along(SIM_N)) {
    for (j in seq_along(SIM_DELTA)) {
      val <- heatmap_data$smsd_total[
        heatmap_data$n == SIM_N[i] & abs(heatmap_data$delta - SIM_DELTA[j]) < 1e-10
      ]
      if (length(val) > 0) mat[i, j] <- val[1]
    }
  }
  
  col_palette <- colorRampPalette(c("white", "lightyellow", "orange", "red"))(100)
  
  par(mar = c(5, 5, 4, 2))
  image(1:length(SIM_DELTA), 1:length(SIM_N), t(mat),
        col = col_palette, axes = FALSE,
        xlab = expression(delta ~ "(Multicollinearity)"), 
        ylab = "n (Sample Size)",
        main = "SMSD Heatmap (GRTK, BIC, rho=0.60)")
  axis(1, at = 1:length(SIM_DELTA), labels = SIM_DELTA)
  axis(2, at = 1:length(SIM_N), labels = SIM_N)
  box()
  
  for (i in seq_along(SIM_DELTA)) {
    for (j in seq_along(SIM_N)) {
      if (!is.na(mat[j, i])) {
        text(i, j, sprintf("%.3f", mat[j, i]), cex = 1.2, font = 2)
      }
    }
  }
}
dev.off()
cat("Figure 8 saved.\n")

# -----------------------------------------------------------------------------
# FIGURE 9: Comprehensive Criteria Comparison (GRTK only)
# Shows all 4 criteria across different configurations
# -----------------------------------------------------------------------------

png(file.path(OUT_DIR, "fig9_criteria_performance_grid.png"), width = 1600, height = 1200, res = 120)
par(mfrow = c(3, 4), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))

# For each sample size, show SMSD, fMSE, PredMSE by criteria
for (n_val in SIM_N) {
  # SMSD by criteria
  crit_data <- subset(summary_smsd,
                      n == n_val & abs(delta - 0.90) < 1e-10 & abs(rho - 0.60) < 1e-10 &
                        f_id == 1 & err_case == "AR1_Gauss" & estimator == "GRTK")
  
  if (nrow(crit_data) >= 4) {
    crit_data <- crit_data[order(match(crit_data$criterion, crit_names)), ]
    
    barplot(crit_data$smsd_total, names.arg = crit_data$criterion,
            col = crit_cols, ylab = "SMSD",
            main = bquote("SMSD: n=" * .(n_val)))
    
    barplot(crit_data$f_ise_mean, names.arg = crit_data$criterion,
            col = crit_cols, ylab = "MSE(f)",
            main = bquote("MSE(f): n=" * .(n_val)))
    
    barplot(crit_data$pred_mse_mean, names.arg = crit_data$criterion,
            col = crit_cols, ylab = "Pred MSE",
            main = bquote("Pred MSE: n=" * .(n_val)))
    
    # Tuning parameter k selected
    barplot(crit_data$k_mean, names.arg = crit_data$criterion,
            col = crit_cols, ylab = expression(hat(k)),
            main = bquote("Ridge k: n=" * .(n_val)))
  } else {
    for (j in 1:4) {
      plot.new()
      plot.window(xlim = c(0, 1), ylim = c(0, 1))
      text(0.5, 0.5, "No data", cex = 1.2)
    }
  }
}
mtext("Selection Criteria Performance Comparison (GRTK, delta=0.90, rho=0.60)", 
      outer = TRUE, line = 1, cex = 1.3)
dev.off()
cat("Figure 9 saved.\n")

# -----------------------------------------------------------------------------
# FIGURE 10: Criteria Performance Lines - SMSD across all conditions
# -----------------------------------------------------------------------------

png(file.path(OUT_DIR, "fig10_criteria_SMSD_lines.png"), width = 1400, height = 1000, res = 120)
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

# For each (delta, rho), show SMSD vs n for all 4 criteria
for (delta_val in SIM_DELTA) {
  for (rho_val in SIM_RHO) {
    crit_data <- subset(summary_smsd,
                        abs(delta - delta_val) < 1e-10 & abs(rho - rho_val) < 1e-10 &
                          f_id == 1 & err_case == "AR1_Gauss" & estimator == "GRTK")
    
    if (nrow(crit_data) > 0) {
      ylim_range <- range(crit_data$smsd_total, na.rm = TRUE)
      ylim_range <- ylim_range + c(-0.05, 0.1) * diff(ylim_range)
      
      plot(NA, xlim = range(SIM_N), ylim = ylim_range,
           xlab = "n", ylab = "SMSD",
           main = bquote(delta == .(delta_val) ~ ", " ~ rho == .(rho_val)))
      
      for (i in seq_along(crit_names)) {
        dd <- crit_data[crit_data$criterion == crit_names[i], ]
        dd <- dd[order(dd$n), ]
        if (nrow(dd) > 0) {
          lines(dd$n, dd$smsd_total, type = "b",
                pch = crit_pchs[i], col = crit_cols[i], lty = i, lwd = 1.5)
        }
      }
      
      if (delta_val == SIM_DELTA[1] && rho_val == SIM_RHO[1]) {
        legend("topright", legend = crit_names, col = crit_cols, pch = crit_pchs,
               lty = 1:4, lwd = 1.5, bty = "n", cex = 0.8)
      }
    } else {
      plot.new()
      plot.window(xlim = c(0, 1), ylim = c(0, 1))
      text(0.5, 0.5, "No data", cex = 1.2)
      title(main = bquote(delta == .(delta_val) ~ ", " ~ rho == .(rho_val)))
    }
  }
}
mtext("SMSD by Selection Criterion vs Sample Size (GRTK)", outer = TRUE, line = 0.5, cex = 1.2)
dev.off()
cat("Figure 10 saved.\n")

# -----------------------------------------------------------------------------
# FIGURE 11: Estimator x Criterion Interaction Heatmap
# -----------------------------------------------------------------------------

png(file.path(OUT_DIR, "fig11_estimator_criterion_heatmap.png"), width = 1200, height = 800, res = 120)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2), oma = c(0, 0, 2, 0))

# Create heatmaps for each metric: SMSD, fMSE, PredMSE
# Rows = Estimators (KS, RTK, GRTK), Cols = Criteria (GCV, AICc, BIC, RECP)

for (metric in c("smsd_total", "f_ise_mean", "pred_mse_mean")) {
  heat_data <- subset(summary_smsd,
                      n == 300 & abs(delta - 0.90) < 1e-10 & abs(rho - 0.60) < 1e-10 &
                        f_id == 1 & err_case == "AR1_Gauss")
  
  if (nrow(heat_data) >= 12) {
    mat <- matrix(NA, nrow = 3, ncol = 4)
    rownames(mat) <- est_levels
    colnames(mat) <- crit_names
    
    for (i in seq_along(est_levels)) {
      for (j in seq_along(crit_names)) {
        val <- heat_data[[metric]][heat_data$estimator == est_levels[i] & 
                                     heat_data$criterion == crit_names[j]]
        if (length(val) > 0) mat[i, j] <- val[1]
      }
    }
    
    col_palette <- colorRampPalette(c("darkgreen", "yellow", "red"))(100)
    
    metric_label <- switch(metric,
                           "smsd_total" = "SMSD",
                           "f_ise_mean" = "MSE(f)",
                           "pred_mse_mean" = "Pred MSE")
    
    image(1:4, 1:3, t(mat), col = col_palette, axes = FALSE,
          xlab = "Criterion", ylab = "Estimator", main = metric_label)
    axis(1, at = 1:4, labels = crit_names)
    axis(2, at = 1:3, labels = est_levels, las = 1)
    box()
    
    # Add values
    for (i in 1:4) {
      for (j in 1:3) {
        if (!is.na(mat[j, i])) {
          text(i, j, sprintf("%.3f", mat[j, i]), cex = 1.0, font = 2)
        }
      }
    }
  } else {
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    text(0.5, 0.5, "Insufficient data", cex = 1.2)
  }
}
mtext("Estimator x Criterion Performance (n=300, delta=0.90, rho=0.60)", outer = TRUE, line = 0.5, cex = 1.2)
dev.off()
cat("Figure 11 saved.\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n" , paste(rep("=", 60), collapse = ""), "\n")
cat("SIMULATION COMPLETE\n")

