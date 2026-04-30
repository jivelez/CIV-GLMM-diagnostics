# =============================================================
# CIV vs White IME test for Poisson GLMMs (random intercept)
# Scenarios
#   H0: Correct specification
#   M0: Fixed-effects misspecification (quadratic term omitted)
#   M1: Random-effects distribution misspecification (t random intercept)
#   M2: Random-effects variance heterogeneity across clusters
#
# Maintained fitted model in all cases:
#   y_ij | b_i ~ Poisson(exp(beta0 + beta1*x_ij + b_i))
#   b_i ~ Normal(0, sigma^2)
# =============================================================
suppressPackageStartupMessages({
  required <- c("lme4", "numDeriv", "mvtnorm", "MASS", "ggplot2", "dplyr", "tidyr", "purrr", "readr")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
  library(lme4)
  library(numDeriv)
  library(mvtnorm)
  library(MASS)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
})

# -----------------------------
# User settings
# -----------------------------
set.seed(123)
sample_sizes <- c(50,100,200,400,600,800,1000,2000,3000)      # increase as needed
cluster_size <- 4                    # repeated measures per cluster
R_reps <- 1000                        # Monte Carlo repetitions per design cell
alpha <- 0.05
quadrature_points <- 15              # Gauss-Hermite points for marginal LL
ridge <- 1e-8
output_dir <- "civ_white_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# True / maintained parameters
beta_true <- c(0.5, 0.7)             # beta0, beta1 in maintained model
sigma_true <- 0.6                    # random-intercept SD
beta_quad <- 0.5                     # used in M0
nu_t <- 5                            # df for t random effects in M1
hetero_gamma <- 0.6                  # variance heterogeneity in M2

# -----------------------------
# Utilities
# -----------------------------
vech_mat <- function(A) A[lower.tri(A, diag = TRUE)]
expit <- function(x) 1 / (1 + exp(-x))

# Gauss-Hermite nodes/weights for Normal(0, sigma^2)
get_gh <- function(m = 15) {
  gh <- statmod::gauss.quad(m, kind = "hermite")
  list(nodes = gh$nodes, weights = gh$weights)
}
if (!requireNamespace("statmod", quietly = TRUE)) {
  install.packages("statmod", repos = "https://cloud.r-project.org")
}
library(statmod)
GH <- get_gh(quadrature_points)

# -----------------------------
# Data generation
# -----------------------------
simulate_dataset <- function(n_clusters, cluster_size, scenario = c("H0", "M0", "M1", "M2")) {
  scenario <- match.arg(scenario)
  id <- rep(seq_len(n_clusters), each = cluster_size)
  t_idx <- rep(seq_len(cluster_size), times = n_clusters)
  x <- rnorm(n_clusters * cluster_size)

  z_cluster <- rnorm(n_clusters)  # used in M2 for variance heterogeneity

  if (scenario == "H0") {
    b <- rnorm(n_clusters, mean = 0, sd = sigma_true)
    eta <- beta_true[1] + beta_true[2] * x + b[id]
  }

  if (scenario == "M0") {
    b <- rnorm(n_clusters, mean = 0, sd = sigma_true)
    eta <- beta_true[1] + beta_true[2] * x + beta_quad * x^2 + b[id]
  }

  if (scenario == "M1") {
    b <- sqrt((nu_t - 2) / nu_t) * rt(n_clusters, df = nu_t) * sigma_true
    eta <- beta_true[1] + beta_true[2] * x + b[id]
  }

  if (scenario == "M2") {
    sigma_i <- sigma_true * exp(hetero_gamma * scale(z_cluster)[, 1])
    b <- rnorm(n_clusters, mean = 0, sd = sigma_i)
    eta <- beta_true[1] + beta_true[2] * x + b[id]
  }

  mu <- exp(eta)
  y <- rpois(length(mu), lambda = mu)

  data.frame(
    id = factor(id),
    t = t_idx,
    x = x,
    y = y,
    z_cluster = z_cluster[id]
  )
}

# -----------------------------
# Fit maintained model
# -----------------------------
fit_maintained_glmm <- function(dat) {
  fit <- suppressWarnings(
    glmer(
      y ~ x + (1 | id),
      data = dat,
      family = poisson(link = "log"),
      nAGQ = 10,
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
  )
  fit
}

extract_theta <- function(fit) {
  beta_hat <- fixef(fit)
  sigma_hat <- as.numeric(attr(VarCorr(fit)$id, "stddev"))
  c(beta_hat[[1]], beta_hat[[2]], log(sigma_hat))
}

# -----------------------------
# Cluster marginal log-likelihood for the maintained model
# params = (beta0, beta1, log_sigma)
# -----------------------------
cluster_loglik_poisson_ri <- function(params, y_i, x_i, gh = GH) {
  beta0 <- params[1]
  beta1 <- params[2]
  sigma <- exp(params[3])

  nodes <- gh$nodes
  weights <- gh$weights

  # For N(0, sigma^2): b = sqrt(2) * sigma * node
  b_vals <- sqrt(2) * sigma * nodes

  # log integrand excluding normal density because GH integrates exp(-x^2)
  log_terms <- vapply(seq_along(b_vals), function(k) {
    eta <- beta0 + beta1 * x_i + b_vals[k]
    sum(dpois(y_i, lambda = exp(eta), log = TRUE))
  }, numeric(1))

  # stable weighted log-sum-exp
  a <- max(log_terms)
  val <- (1 / sqrt(pi)) * sum(weights * exp(log_terms - a))
  ll <- log(val) + a
  as.numeric(ll)
}

# -----------------------------
# Clusterwise scores and Hessians
# -----------------------------
cluster_derivatives <- function(dat, theta_hat) {
  split_dat <- split(dat, dat$id)

  out <- lapply(split_dat, function(df_i) {
    y_i <- df_i$y
    x_i <- df_i$x

    f <- function(par) cluster_loglik_poisson_ri(par, y_i = y_i, x_i = x_i)

    s_i <- tryCatch(numDeriv::grad(f, x = theta_hat), error = function(e) rep(NA_real_, length(theta_hat)))
    Hobs_i <- tryCatch(numDeriv::hessian(f, x = theta_hat), error = function(e) matrix(NA_real_, length(theta_hat), length(theta_hat)))

    list(score = s_i, hessian = Hobs_i)
  })

  scores <- lapply(out, `[[`, "score")
  hessians <- lapply(out, `[[`, "hessian")
  list(scores = scores, hessians = hessians)
}

# -----------------------------
# White's IME test
# Using Delta_i = vech(s_i s_i' + d2 l_i)
# -----------------------------
white_ime_test <- function(scores, hessians, ridge = 1e-8) {
  n <- length(scores)
  p <- length(scores[[1]])

  bad <- any(vapply(scores, function(s) any(!is.finite(s)), logical(1))) ||
    any(vapply(hessians, function(H) any(!is.finite(H)), logical(1)))
  if (bad) {
    return(list(statistic = NA_real_, df = NA_integer_, p_value = NA_real_, failed = TRUE, ginv_used = NA))
  }

  Delta <- do.call(rbind, lapply(seq_len(n), function(i) {
    vech_mat(tcrossprod(scores[[i]]) + hessians[[i]])
  }))

  Delta_bar <- colMeans(Delta)
  Omega_hat <- cov(Delta)
  r <- length(Delta_bar)

  if (any(!is.finite(Omega_hat))) {
    return(list(statistic = NA_real_, df = r, p_value = NA_real_, failed = TRUE, ginv_used = NA))
  }

  eig_min <- tryCatch(min(eigen(Omega_hat, symmetric = TRUE, only.values = TRUE)$values), error = function(e) NA_real_)
  if (!is.finite(eig_min) || eig_min < ridge) {
    Omega_hat <- Omega_hat + diag(ridge, r)
  }

  ginv_used <- FALSE
  Omega_inv <- tryCatch(solve(Omega_hat), error = function(e) {
    ginv_used <<- TRUE
    MASS::ginv(Omega_hat)
  })

  T_IME <- as.numeric(n * t(Delta_bar) %*% Omega_inv %*% Delta_bar)
  p_val <- pchisq(T_IME, df = r, lower.tail = FALSE)

  list(statistic = T_IME, df = r, p_value = p_val, failed = FALSE, ginv_used = ginv_used)
}

# -----------------------------
# CIV statistics
# -----------------------------
civ_statistics <- function(scores, hessians, ridge = 1e-8) {
  n <- length(scores)
  p <- length(scores[[1]])

  bad <- any(vapply(scores, function(s) any(!is.finite(s)), logical(1))) ||
    any(vapply(hessians, function(H) any(!is.finite(H)), logical(1)))
  if (bad) {
    return(list(Ttr = NA_real_, Tld = NA_real_, failed = TRUE))
  }

  Jn <- Reduce(`+`, lapply(scores, function(s) tcrossprod(s))) / n
  Hn <- Reduce(`+`, lapply(hessians, function(H) -H)) / n

  if (any(!is.finite(Jn)) || any(!is.finite(Hn))) {
    return(list(Ttr = NA_real_, Tld = NA_real_, failed = TRUE))
  }

  eig_min <- tryCatch(min(eigen(Hn, symmetric = TRUE, only.values = TRUE)$values), error = function(e) NA_real_)
  if (!is.finite(eig_min) || eig_min < ridge) {
    Hn <- Hn + diag(ridge, nrow(Hn))
  }

  Hn_inv <- tryCatch(solve(Hn), error = function(e) MASS::ginv(Hn))
  Rn <- Hn_inv %*% Jn

  Ttr <- n * (sum(diag(Rn)) / p - 1)^2

  det_obj <- tryCatch(determinant(Rn, logarithm = TRUE), error = function(e) NULL)
  if (is.null(det_obj) || det_obj$sign <= 0) {
    Tld <- NA_real_
  } else {
    Tld <- n * (as.numeric(det_obj$modulus))^2
  }

  list(Ttr = Ttr, Tld = Tld, failed = FALSE)
}

# -----------------------------
# One Monte Carlo replication
# -----------------------------
run_one_rep <- function(n_clusters, scenario) {
  t0 <- proc.time()[3]

  dat <- simulate_dataset(n_clusters = n_clusters, cluster_size = cluster_size, scenario = scenario)

  fit <- tryCatch(fit_maintained_glmm(dat), error = function(e) NULL)
  if (is.null(fit)) {
    return(tibble(
      n = n_clusters,
      scenario = scenario,
      fit_failed = TRUE,
      white_stat = NA_real_, white_p = NA_real_, white_ginv = NA,
      civ_tr = NA_real_, civ_ld = NA_real_,
      elapsed = proc.time()[3] - t0
    ))
  }

  theta_hat <- tryCatch(extract_theta(fit), error = function(e) rep(NA_real_, 3))
  if (any(!is.finite(theta_hat))) {
    return(tibble(
      n = n_clusters,
      scenario = scenario,
      fit_failed = TRUE,
      white_stat = NA_real_, white_p = NA_real_, white_ginv = NA,
      civ_tr = NA_real_, civ_ld = NA_real_,
      elapsed = proc.time()[3] - t0
    ))
  }

  derivs <- cluster_derivatives(dat, theta_hat)
  white <- white_ime_test(derivs$scores, derivs$hessians, ridge = ridge)
  civ <- civ_statistics(derivs$scores, derivs$hessians, ridge = ridge)

  tibble(
    n = n_clusters,
    scenario = scenario,
    fit_failed = FALSE,
    white_stat = white$statistic,
    white_p = white$p_value,
    white_ginv = white$ginv_used,
    civ_tr = civ$Ttr,
    civ_ld = civ$Tld,
    elapsed = proc.time()[3] - t0
  )
}

# -----------------------------
# Simulation runner
# -----------------------------
run_simulations <- function(sample_sizes, scenarios, R_reps) {
  res <- vector("list", length(sample_sizes) * length(scenarios))
  idx <- 1

  for (n in sample_sizes) {
    for (sc in scenarios) {
      message(sprintf("Running scenario=%s, n=%s", sc, n))
      reps <- vector("list", R_reps)
      for (r in seq_len(R_reps)) {
        reps[[r]] <- run_one_rep(n_clusters = n, scenario = sc)
        if (r %% 20 == 0) message(sprintf("  rep %d/%d", r, R_reps))
      }
      res[[idx]] <- bind_rows(reps)
      idx <- idx + 1
    }
  }
  bind_rows(res)
}

# -----------------------------
# Calibrate CIV critical values under H0
# -----------------------------
calibrate_civ_critical_values <- function(raw_results, alpha = 0.05) {
  raw_results %>%
    filter(scenario == "H0") %>%
    group_by(n) %>%
    summarise(
      crit_tr = quantile(civ_tr, probs = 1 - alpha, na.rm = TRUE),
      crit_ld = quantile(civ_ld, probs = 1 - alpha, na.rm = TRUE),
      .groups = "drop"
    )
}

# -----------------------------
# Summaries
# -----------------------------
make_summary_tables <- function(raw_results, civ_crit, alpha = 0.05) {
  out <- raw_results %>%
    left_join(civ_crit, by = "n") %>%
    mutate(
      reject_white = white_p < alpha,
      reject_civ_tr = civ_tr > crit_tr,
      reject_civ_ld = civ_ld > crit_ld,
      numerical_failure = fit_failed | is.na(white_stat) | is.na(civ_tr)
    ) %>%
    group_by(scenario, n) %>%
    summarise(
      reps = n(),
      white_reject = mean(reject_white, na.rm = TRUE),
      civ_tr_reject = mean(reject_civ_tr, na.rm = TRUE),
      civ_ld_reject = mean(reject_civ_ld, na.rm = TRUE),
      white_ginv_rate = mean(white_ginv, na.rm = TRUE),
      fit_fail_rate = mean(fit_failed, na.rm = TRUE),
      numerical_failure_rate = mean(numerical_failure, na.rm = TRUE),
      mean_runtime_sec = mean(elapsed, na.rm = TRUE),
      .groups = "drop"
    )

  out
}

# -----------------------------
# Figures
# -----------------------------
make_plots <- function(summary_tbl, output_dir) {
    
    x_breaks_fun <- scales::pretty_breaks(n = 6)
    
    common_x_theme <- theme(
        axis.text.x = element_text(angle = 90, hjust = 1)
    )
    
    size_tbl <- summary_tbl %>%
        filter(scenario == "H0") %>%
        select(n, white_reject, civ_tr_reject, civ_ld_reject) %>%
        pivot_longer(-n, names_to = "test", values_to = "rejection")
    
    # Empirical Type I error under correct specification
    p_size <- ggplot(size_tbl, aes(x = n, y = rejection, color = test)) +
        geom_line(linewidth = 0.8) +
        geom_point(size = 2) +
        geom_hline(yintercept = alpha, linetype = 2) +
        scale_x_continuous(
            breaks = x_breaks_fun,
            labels = scales::label_number(big.mark = ",")
        ) +
        labs(
            x = "Number of clusters (n)",
            y = "Empirical Type I error",
            color = "Test",
            title = ""
        ) +
        theme_bw(base_size = 12) +
        common_x_theme
    
    ggsave(
        file.path(output_dir, "figure_size.png"),
        p_size, width = 7, height = 4.2, dpi = 300
    )
    
    power_tbl <- summary_tbl %>%
        filter(scenario != "H0") %>%
        select(scenario, n, white_reject, civ_tr_reject, civ_ld_reject) %>%
        pivot_longer(
            cols = c(white_reject, civ_tr_reject, civ_ld_reject),
            names_to = "test",
            values_to = "rejection"
        )
    
    # Empirical power under misspecification
    p_power <- ggplot(power_tbl, aes(x = n, y = rejection, color = test)) +
        geom_line(linewidth = 0.8) +
        geom_point(size = 2) +
        facet_wrap(~ scenario, nrow = 1) +
        scale_x_continuous(
            breaks = x_breaks_fun,
            labels = scales::label_number(big.mark = ",")
        ) +
        labs(
            x = "Number of clusters (n)",
            y = "Empirical power",
            color = "Test",
            title = ""
        ) +
        theme_bw(base_size = 12) +
        common_x_theme
    
    ggsave(
        file.path(output_dir, "figure_power.png"),
        p_power, width = 10, height = 4.2, dpi = 300
    )
    
    fail_tbl <- summary_tbl %>%
        select(scenario, n, white_ginv_rate, numerical_failure_rate) %>%
        pivot_longer(
            cols = c(white_ginv_rate, numerical_failure_rate),
            names_to = "metric",
            values_to = "value"
        )
    
    # Numerical stability diagnostics
    p_fail <- ggplot(fail_tbl, aes(x = n, y = value, color = metric)) +
        geom_line(linewidth = 0.8) +
        geom_point(size = 2) +
        facet_wrap(~ scenario, nrow = 1) +
        scale_x_continuous(
            breaks = x_breaks_fun,
            labels = scales::label_number(big.mark = ",")
        ) +
        labs(
            x = "Number of clusters (n)",
            y = "Rate",
            color = "Metric",
            title = ""
        ) +
        theme_bw(base_size = 12) +
        common_x_theme
    
    ggsave(
        file.path(output_dir, "figure_numerical_stability.png"),
        p_fail, width = 10, height = 4.2, dpi = 300
    )
}

# -----------------------------
# Main execution
# -----------------------------
if(!require(doParallel)) install.packages('doParallel')
require(doParallel)

## simulations
scenarios <- c("H0", "M0", "M1", "M2")
raw_results <- run_simulations(sample_sizes = sample_sizes, scenarios = scenarios, R_reps = R_reps)
write_csv(raw_results, file.path(output_dir, "raw_results.csv"))

## calibracion
civ_crit <- calibrate_civ_critical_values(raw_results, alpha = alpha)
write_csv(civ_crit, file.path(output_dir, "civ_critical_values.csv"))

## tables
summary_tbl <- make_summary_tables(raw_results, civ_crit, alpha = alpha)
write_csv(summary_tbl, file.path(output_dir, "summary_table.csv"))

## plots
make_plots(summary_tbl, output_dir)

message("\nDone. Files saved in: ", normalizePath(output_dir))
message("Generated files:")
print(list.files(output_dir))
