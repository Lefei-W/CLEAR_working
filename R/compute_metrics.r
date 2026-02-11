#### 1. CLEAR core specificity ####

safe_l2 <- function(x) { # on the rows (peaks)
  s <- sum(x^2)
  if (s == 0) x else x / sqrt(s)  # l2 euclidean normalization avoiding row sum of 0 (happen with the Human-scATAC-Corpus)
}

add_rank <- function(M, rn, cn) { # ranking assign to cell types (columns) which row names and column names
  R <- apply(M, 2, function(x) rank(x, ties.method = "min")) # rank 1 is the smallest, and tied values as the minimal rank
  dimnames(R) <- list(rn, cn)
  R # ranking matrix 
}

compute_metrics <- function(mat, rn, cn) { # generally use the PeakMatrix from ArchR rn are the peaks and cn the cell types
  mat_raw <- mat # keep a copy of the raw matrix for composite and logged metrics)
  dimnames(mat_raw) <- list(rn, cn)
  mat_qn <- preprocessCore::normalize.quantiles(mat)  # quantile normalisation based on CHEERS (keep empirical distri across cell types)
  dimnames(mat_qn) <- list(rn, cn)
  mat <- mat_qn # use this normalised matrix for all other specificity calculations 

  # composite matrix calculations
  log_raw <- log2(mat_raw + 1)
  dimnames(log_raw) <- list(rn, cn)
  
  scaled_mat  <- t(apply(mat, 1, scale));       dimnames(scaled_mat)  <- list(rn, cn) # Z scaling
  l2norm_mat  <- t(apply(mat, 1, safe_l2));     dimnames(l2norm_mat)  <- list(rn, cn) # l2 CHEERS norm
  compsc1_mat <- l2norm_mat + log2(mat_raw + 1);    dimnames(compsc1_mat) <- list(rn, cn) # Composite calculations 
  compsc2_mat <- l2norm_mat * log2(mat_raw + 1);    dimnames(compsc2_mat) <- list(rn, cn)
  
  # Tau calclulations 
  calc_tau <- function(x) if (max(x) == 0) 0 else sum(1 - (x / max(x))) / (length(x) - 1)
  tau_vec  <- apply(mat, 1, calc_tau)
  tau_mat  <- matrix(
    rep(tau_vec, ncol(mat)),
    nrow = nrow(mat), ncol = ncol(mat),
    dimnames = list(rn, cn)
  )
  
  compsc1_tau <- tau_mat + log2(mat_raw + 1); dimnames(compsc1_tau) <- list(rn, cn)
  compsc2_tau <- tau_mat * log2(mat_raw + 1); dimnames(compsc2_tau) <- list(rn, cn)
  
  list(
    zscaled             = scaled_mat,
    logged              = log_raw,
    quantile_normalized = mat_qn,
    l2_norm             = l2norm_mat,
    comp1               = compsc1_mat,
    comp2               = compsc2_mat,
    spec_rank           = add_rank(l2norm_mat,  rn, cn),
    comp1_rank          = add_rank(compsc1_mat, rn, cn),
    comp2_rank          = add_rank(compsc2_mat, rn, cn),
    tau                 = tau_mat,
    tau_rank            = add_rank(tau_mat, rn, cn),
    compsc1_tau         = compsc1_tau,
    compsc1_tau_rank    = add_rank(compsc1_tau, rn, cn),
    compsc2_tau         = compsc2_tau,
    compsc2_tau_rank    = add_rank(compsc2_tau, rn, cn)
  )
}

DERIVED_ASSAYS <- c( # names of assays generated and to be teseted for enrichment 
  "zscaled", "logged", "quantile_normalized", "l2_norm","comp1","comp2","spec_rank","comp1_rank","comp2_rank",
  "tau","tau_rank","compsc1_tau","compsc1_tau_rank","compsc2_tau","compsc2_tau_rank"
)

