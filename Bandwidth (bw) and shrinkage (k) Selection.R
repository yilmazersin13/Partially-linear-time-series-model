################################################
# tuning_selection.R (by Ersin Yılmaz)
# 
#   - Select (bw, k) by GCV / AICc / BIC / RECP for semiparametric PLM
#   - Produce: table of selected (bw,k), 3D tuning surfaces, optional contours
#
# Notes:
#   - Pure base-R implementation (no external packages).
#   - Gaussian kernel Nadaraya–Watson smoothing matrix W(bw).
#   - AR(1) working correlation estimated from pilot residuals.
##-----------------------------------------------------------

options(stringsAsFactors = FALSE)

# utilities --------------------------------------

gaussian_kernel <- function(u) exp(-0.5 * u^2)

NW_smoother_matrix <- function(z, bw) {
  
  # Nadaraya–Watson smoothing matrix at observed points z
  n <- length(z)
  D <- outer(z, z, "-") / bw
  K <- gaussian_kernel(D)
  rs <- rowSums(K) + 1e-12
  W <- K / rs
  W
}

AR1_corr_matrix <- function(n, rho) {
  rho <- max(min(as.numeric(rho), 0.999), -0.999)
  toeplitz(rho^(0:(n - 1)))
}

safe_chol_inv <- function(A) {
  # stable inverse for SPD matrices
  R <- try(chol(A), silent = TRUE)
  if (inherits(R, "try-error")) {
    return(solve(A))
  }
  chol2inv(R)
}

estimate_rho_ar1 <- function(res, lag = 1, eta = 1e-3) {
  ac <- acf(res, plot = FALSE, lag.max = lag)$acf
  rho_hat <- as.numeric(ac[2])  # lag-1
  # truncation for numerical stability
  rho_hat <- min(1 - eta, max(-1 + eta, rho_hat))
  rho_hat
}

#  Core fit: PLM with kernel + ridge + GLS --------------

fit_plm_once <- function(x, z, y, bw, k, Rinv, W = NULL) {
  # If W is provided, use it; otherwise compute
  n <- length(y)
  p <- ncol(x)
  
  if (is.null(W)) {
    W <- NW_smoother_matrix(z, bw)
  }
  Iw <- diag(n) - W
  
  xtil <- Iw %*% x
  ytil <- Iw %*% y
  
  A <- crossprod(xtil, Rinv %*% xtil) + k * diag(p)
  b <- crossprod(xtil, Rinv %*% ytil)
  
  beta_hat <- solve(A, b)
  
  # linear smoother form:
  # yhat = W y + (I - W) X beta_hat
  yhat <- as.vector(W %*% y + (Iw %*% x) %*% beta_hat)
  fhat <- as.vector(W %*% (y - x %*% beta_hat))
  
  list(beta_hat = beta_hat, fhat = fhat, yhat = yhat,
       W = W, Iw = Iw, xtil = xtil, A = A)
}

# df(H) without building full H:
# H = W + (I-W) X A^{-1} xtil' Rinv (I-W)
# df = tr(H) = tr(W) + tr( A^{-1} * (xtil' Rinv (I-W)^2 X) )
# We compute C = xtil' Rinv (I-W) xtil = xtil' Rinv (I-W)^2 X, since xtil=(I-W)X.
df_fast <- function(W, Iw, xtil, A, Rinv) {
  trW <- sum(diag(W))
  C   <- crossprod(xtil, Rinv %*% (Iw %*% xtil))   # p x p
  tr2 <- sum(diag(solve(A, C)))
  trW + tr2
}

# Build full H only when needed (RECP / tr(HH')):
H_matrix <- function(W, Iw, x, xtil, A, Rinv) {
  # H = W + (I-W) X A^{-1} xtil' Rinv (I-W)
  B  <- solve(A, crossprod(xtil, Rinv %*% Iw))  # p x n
  H2 <- (Iw %*% x) %*% B                        # n x n
  W + H2
}

# ---------------------- Main tuning selection --------------------------------

Selection_bw_k <- function(x, z, y,
                           bw_seq, k_seq,
                           lag = 1,
                           criteria = c("GCV", "AICc", "BIC", "RECP"),
                           R_mode = c("AR1", "I"),
                           k_pilot = NULL,
                           eta_rho = 1e-3,
                           make_plots = TRUE,
                           make_contours = TRUE,
                           save_outputs = TRUE,
                           out_dir = "outputs",
                           prefix = "tuning_demo") {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  R_mode <- match.arg(R_mode)
  n <- length(y)
  p <- ncol(x)
  
  if (is.null(k_pilot)) k_pilot <- median(k_seq)
  
  # ---- Determine working inverse correlation
  rho_hat <- 0
  Rinv <- diag(n)
  
  if (R_mode == "AR1") {
    bw_pilot <- median(bw_seq)
    pilot <- fit_plm_once(x, z, y, bw = bw_pilot, k = k_pilot, Rinv = diag(n))
    res_pilot <- y - pilot$yhat
    rho_hat <- estimate_rho_ar1(res_pilot, lag = lag, eta = eta_rho)
    
    R <- AR1_corr_matrix(n, rho_hat)
    Rinv <- safe_chol_inv(R)
  }
  
  L <- length(bw_seq)
  M <- length(k_seq)
  
  GCV  <- matrix(NA_real_, L, M)
  AICc <- matrix(NA_real_, L, M)
  BIC  <- matrix(NA_real_, L, M)
  RECP <- matrix(NA_real_, L, M)
  
  need_RECP <- ("RECP" %in% criteria)
  
  # Pre-compute and cache all W matrices for speed
  W_cache <- vector("list", L)
  Iw_cache <- vector("list", L)
  xtil_cache <- vector("list", L)
  ytil_cache <- vector("list", L)
  Xt_Rinv_cache <- vector("list", L)
  
  for (i in seq_len(L)) {
    bw <- bw_seq[i]
    W_cache[[i]]  <- NW_smoother_matrix(z, bw)
    Iw_cache[[i]] <- diag(n) - W_cache[[i]]
    xtil_cache[[i]] <- Iw_cache[[i]] %*% x
    ytil_cache[[i]] <- Iw_cache[[i]] %*% y
    Xt_Rinv_cache[[i]] <- crossprod(xtil_cache[[i]], Rinv)  # p x n
  }
  
  for (i in seq_len(L)) {
    W  <- W_cache[[i]]
    Iw <- Iw_cache[[i]]
    xtil <- xtil_cache[[i]]
    ytil <- ytil_cache[[i]]
    Xt_Rinv <- Xt_Rinv_cache[[i]]
    
    for (j in seq_len(M)) {
      k <- k_seq[j]
      
      A <- Xt_Rinv %*% xtil + k * diag(p)
      b <- Xt_Rinv %*% ytil
      
      beta_hat <- solve(A, b)
      yhat <- as.vector(W %*% y + (Iw %*% x) %*% beta_hat)
      
      RSS <- sum((y - yhat)^2)
      
      # df and residual df
      df <- df_fast(W, Iw, xtil, A, Rinv)
      df_res <- n - df
      if (!is.finite(df_res) || df_res <= 5) df_res <- 5
      

      
      if ("GCV" %in% criteria) {
        # Standard GCV: (RSS/n) / (1 - df/n)^2
        GCV[i, j] <- (RSS / n) / (1 - df / n)^2
      }
      
      if ("AICc" %in% criteria) {
        # Hurvich-Simonoff-Tsai AICc for smoothing:
        # AICc = log(RSS/n) + (1 + df/n) / (1 - (df+2)/n)
        denom <- 1 - (df + 2) / n
        if (denom <= 0.01) denom <- 0.01  # numerical safeguard
        AICc[i, j] <- log(RSS / n) + (1 + df / n) / denom
      }
      
      if ("BIC" %in% criteria) {
        # Standard BIC: n*log(RSS/n) + df*log(n)
        # For minimization on comparable scale to GCV, use exponential form:
        # BIC_scale = (RSS/n) * exp(df * log(n) / n)
        BIC[i, j] <- (RSS / n) * exp(df * log(n) / n)
      }
      
      if (need_RECP) {
        sigma2_hat <- RSS / df_res
        H <- H_matrix(W, Iw, x, xtil, A, Rinv)
        trHHt <- sum(H * H)  # tr(H H')
        RECP[i, j] <- (RSS + sigma2_hat * trHHt) / df_res
      }
    }
  }
  
  get_min_pair <- function(S) {
    idx <- which(S == min(S, na.rm = TRUE), arr.ind = TRUE)[1, ]
    c(bw = bw_seq[idx[1]], k = k_seq[idx[2]])
  }
  
  pairs <- data.frame(
    Criterion = criteria,
    bw_hat = NA_real_,
    k_hat  = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (m in seq_along(criteria)) {
    cc <- criteria[m]
    if (cc == "GCV")  v <- get_min_pair(GCV)
    if (cc == "AICc") v <- get_min_pair(AICc)
    if (cc == "BIC")  v <- get_min_pair(BIC)
    if (cc == "RECP") v <- get_min_pair(RECP)
    pairs$bw_hat[m] <- v["bw"]
    pairs$k_hat[m]  <- v["k"]
  }
  
  if (save_outputs) {
    write.csv(pairs, file.path(out_dir, paste0(prefix, "_selected_pairs.csv")), row.names = FALSE)
    
    if (make_plots) {
      png(file.path(out_dir, paste0(prefix, "_3D_surfaces.png")), width = 1000, height = 800, res = 120)
      par(mfrow = c(2, 2), mar = c(3, 3, 2, 1))
      
      if ("GCV" %in% criteria) {
        persp(bw_seq, k_seq, GCV, main = "GCV", xlab = "bw", ylab = "k",
              zlab = "score", theta = 25, phi = 20, shade = 0.2, expand = 0.6)
      } else plot.new()
      
      if ("AICc" %in% criteria) {
        persp(bw_seq, k_seq, AICc, main = "AICc", xlab = "bw", ylab = "k",
              zlab = "score", theta = 25, phi = 20, shade = 0.2, expand = 0.6)
      } else plot.new()
      
      if ("BIC" %in% criteria) {
        persp(bw_seq, k_seq, BIC, main = "BIC", xlab = "bw", ylab = "k",
              zlab = "score", theta = 25, phi = 20, shade = 0.2, expand = 0.6)
      } else plot.new()
      
      if ("RECP" %in% criteria) {
        persp(bw_seq, k_seq, RECP, main = "RECP", xlab = "bw", ylab = "k",
              zlab = "score", theta = 25, phi = 20, shade = 0.2, expand = 0.6)
      } else plot.new()
      
      dev.off()
    }
    
    if (make_contours) {
      png(file.path(out_dir, paste0(prefix, "_contours.png")), width = 1000, height = 800, res = 120)
      par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
      
      if ("GCV" %in% criteria)  contour(bw_seq, k_seq, GCV,  main = "GCV (contour)",  xlab = "bw", ylab = "k")
      else plot.new()
      if ("AICc" %in% criteria) contour(bw_seq, k_seq, AICc, main = "AICc (contour)", xlab = "bw", ylab = "k")
      else plot.new()
      if ("BIC" %in% criteria)  contour(bw_seq, k_seq, BIC,  main = "BIC (contour)",  xlab = "bw", ylab = "k")
      else plot.new()
      if ("RECP" %in% criteria) contour(bw_seq, k_seq, RECP, main = "RECP (contour)", xlab = "bw", ylab = "k")
      else plot.new()
      
      dev.off()
    }
  }
  
  sel <- new.env()
  sel$pairs   <- pairs
  sel$GCV     <- GCV
  sel$AICc    <- AICc
  sel$BIC     <- BIC
  sel$RECP    <- RECP
  sel$rho_hat <- rho_hat
  sel$R_mode  <- R_mode
  
  return(sel)
}
