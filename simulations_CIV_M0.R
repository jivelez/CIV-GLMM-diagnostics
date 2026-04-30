## ============================================================
## CIV Monte Carlo COMPLETO + PARALELIZADO + SALIDAS POR CASO
## n = 20,100,200,400,600,1000,2000,3000
##
## Para cada caso M0/M1/M2 produce archivos independientes:
##  - results_M0_summary.csv
##  - results_M0_replications.csv
##  - results_M0_log.txt
##  - results_M0_TypeI_Power_only.csv
## (análogos para M1 y M2)
##
## ============================================================
## REQUERIMIENTOS:
## - lme4, statmod, numDeriv, MASS, mvtnorm
## - future.apply, data.table
## ============================================================


## ------------------ OPTIONS ------------------
options(warn = -1)

## ------------------ PACKAGES ------------------
pkgs <- c("lme4", "statmod", "numDeriv", "MASS", "mvtnorm",
          "future.apply", "data.table")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(lme4)
library(statmod)
library(numDeriv)
library(MASS)
library(mvtnorm)
library(future.apply)
library(data.table)

## ------------------ UTILITIES ------------------
logsumexp <- function(a) { m <- max(a); m + log(sum(exp(a - m))) }

safe_hessian <- function(f, x0) {
    H <- hessian(f, x0)
    0.5 * (H + t(H))
}

ridge_invert <- function(H, J, ridge = 1e-6) {
    Hs <- 0.5 * (H + t(H))
    lam <- ridge * mean(diag(Hs))
    Hreg <- Hs + diag(lam, nrow(Hs))
    solve(Hreg, J)
}

logdet_spd <- function(A, eps = 1e-12) {
    As <- 0.5 * (A + t(A))
    ev <- eigen(As, symmetric = TRUE, only.values = TRUE)$values
    ev <- pmax(ev, eps)
    sum(log(ev))
}

## ============================================================
## RANDOM EFFECTS
## ============================================================
r_t_re <- function(n, df = 3) {
    stopifnot(df > 2)
    s <- sqrt((df - 2) / df)  # scale so Var=1
    s * rt(n, df = df)
}

G_from_params <- function(sd1, sd2, rho) {
    matrix(c(sd1^2, rho*sd1*sd2,
             rho*sd1*sd2, sd2^2), 2, 2, byrow = TRUE)
}

## ============================================================
## SIMULATORS: H0/H1 for each case M0/M1/M2
## ============================================================
simulate_norm_intercept <- function(n, m, beta0, beta1, sigma) {
    id <- rep(1:n, each = m)
    x  <- runif(n*m, -1, 1)
    b  <- rnorm(n, 0, sigma)
    eta <- beta0 + beta1 * x + b[id]
    y <- rpois(n*m, exp(eta))
    data.frame(y=y, x=x, id=factor(id))
}

## M0: omitted quadratic
simulate_M0_H0 <- function(n, m, beta0, beta1, sigma) {
    simulate_norm_intercept(n, m, beta0, beta1, sigma)
}
simulate_M0_H1 <- function(n, m, beta0, beta1, beta2, sigma) {
    id <- rep(1:n, each=m)
    x  <- runif(n*m, -1, 1)
    b  <- rnorm(n, 0, sigma)
    eta <- beta0 + beta1*x + beta2*(x^2) + b[id]
    y <- rpois(n*m, exp(eta))
    data.frame(y=y, x=x, id=factor(id))
}

## M1: t random intercept vs normal fit
simulate_M1_H0 <- function(n, m, beta0, beta1, sigma) {
    simulate_norm_intercept(n, m, beta0, beta1, sigma)
}
simulate_M1_H1 <- function(n, m, beta0, beta1, df) {
    id <- rep(1:n, each=m)
    x  <- runif(n*m, -1, 1)
    b  <- r_t_re(n, df=df)
    eta <- beta0 + beta1*x + b[id]
    y <- rpois(n*m, exp(eta))
    data.frame(y=y, x=x, id=factor(id))
}

## M2: extra random effect in FIT
simulate_M2_H0 <- function(n, m, beta0, beta1, sd1, sd2, rho) {
    id <- rep(1:n, each=m)
    x  <- runif(n*m, -1, 1)
    z2 <- rnorm(n*m)
    
    G <- G_from_params(sd1, sd2, rho)
    bi <- MASS::mvrnorm(n, mu=c(0,0), Sigma=G)
    b1 <- bi[id, 1]
    b2 <- bi[id, 2]
    
    eta <- beta0 + beta1*x + b1 + b2*z2
    y <- rpois(n*m, exp(eta))
    data.frame(y=y, x=x, z2=z2, id=factor(id))
}
simulate_M2_H1 <- function(n, m, beta0, beta1, sd1) {
    id <- rep(1:n, each=m)
    x  <- runif(n*m, -1, 1)
    z2 <- rnorm(n*m)
    b1 <- rnorm(n, 0, sd1)
    eta <- beta0 + beta1*x + b1[id]
    y <- rpois(n*m, exp(eta))
    data.frame(y=y, x=x, z2=z2, id=factor(id))
}

## ============================================================
## CIV (1D) - Adaptive GH cluster loglik for random intercept
## ============================================================
post_mode_scale_1d <- function(beta0, beta1, sigma, y, x, maxit = 40) {
    b <- 0
    for (it in 1:maxit) {
        eta <- beta0 + beta1*x + b
        mu  <- exp(eta)
        g   <- sum(y - mu) - b/(sigma^2)
        H   <- -sum(mu) - 1/(sigma^2)
        step <- g / H
        b_new <- b - step
        if (abs(b_new - b) < 1e-10) { b <- b_new; break }
        b <- b_new
    }
    s <- sqrt(1 / (-H))
    list(bhat=b, scale=s)
}

cluster_loglik_GH1_adapt <- function(phi, data_i, gh_nodes, gh_weights) {
    beta0 <- phi[1]
    beta1 <- phi[2]
    sigma <- exp(phi[3])
    
    y <- data_i$y
    x <- data_i$x
    
    ms <- post_mode_scale_1d(beta0, beta1, sigma, y, x)
    bh <- ms$bhat
    sc <- ms$scale
    
    b_k <- bh + sc * sqrt(2) * gh_nodes
    
    logLik_k <- vapply(b_k, function(b) {
        eta <- beta0 + beta1*x + b
        mu  <- exp(eta)
        sum(dpois(y, mu, log=TRUE)) + dnorm(b, 0, sigma, log=TRUE)
    }, numeric(1))
    
    log_terms <- log(gh_weights) + logLik_k + (gh_nodes^2)
    log(sc * sqrt(2)) + logsumexp(log_terms)
}

compute_CIV_1d <- function(dat, phi_hat, K=25, ridge=1e-6) {
    dat$id <- as.factor(dat$id)
    ids <- levels(dat$id)
    ncl <- length(ids)
    
    gh <- statmod::gauss.quad(K, kind="hermite")
    
    scores <- vector("list", ncl)
    hessians <- vector("list", ncl)
    
    for (ii in seq_along(ids)) {
        di <- dat[dat$id == ids[ii], c("y","x")]
        li_fun <- function(phi) cluster_loglik_GH1_adapt(phi, di, gh$nodes, gh$weights)
        si <- grad(li_fun, phi_hat)
        Hi <- -safe_hessian(li_fun, phi_hat)
        scores[[ii]] <- si
        hessians[[ii]] <- Hi
    }
    
    p <- length(phi_hat)
    Jn <- Reduce(`+`, lapply(scores, function(s) tcrossprod(s))) / ncl
    Hn <- Reduce(`+`, hessians) / ncl
    Rn <- ridge_invert(Hn, Jn, ridge=ridge)
    
    Ttr <- ncl * ((sum(diag(Rn))/p) - 1)^2
    Tld <- ncl * (logdet_spd(Rn))^2
    list(Ttr=Ttr, Tld=Tld)
}

bootstrap_CIV_1d <- function(dat, B=199, Kgh=25, stat=c("Ttr","Tld"), ridge=1e-6) {
    stat <- match.arg(stat)
    ctrl <- glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
    
    fit <- glmer(y ~ x + (1|id), data=dat, family=poisson(link="log"), control=ctrl)
    beta_hat <- fixef(fit)
    sigma_hat <- as.numeric(sqrt(VarCorr(fit)$id[1,1])); sigma_hat <- max(sigma_hat, 1e-6)
    phi_hat <- c(beta_hat[1], beta_hat[2], log(sigma_hat))
    
    T0 <- compute_CIV_1d(dat, phi_hat, K=Kgh, ridge=ridge)[[stat]]
    
    ids <- levels(as.factor(dat$id))
    ncl <- length(ids)
    Tstar <- rep(NA_real_, B)
    
    for (b in 1:B) {
        bi <- rnorm(ncl, 0, sigma_hat); names(bi) <- ids
        eta <- phi_hat[1] + phi_hat[2]*dat$x + bi[as.character(dat$id)]
        dat_star <- dat
        dat_star$y <- rpois(nrow(dat), exp(eta))
        
        fit_star <- suppressWarnings(try(
            glmer(y ~ x + (1|id), data=dat_star, family=poisson(link="log"), control=ctrl),
            silent=TRUE
        ))
        if (inherits(fit_star, "try-error")) next
        
        beta_s <- fixef(fit_star)
        sigma_s <- as.numeric(sqrt(VarCorr(fit_star)$id[1,1])); sigma_s <- max(sigma_s, 1e-6)
        phi_s <- c(beta_s[1], beta_s[2], log(sigma_s))
        Tstar[b] <- compute_CIV_1d(dat_star, phi_s, K=Kgh, ridge=ridge)[[stat]]
    }
    
    Tstar <- Tstar[is.finite(Tstar)]
    pval <- (1 + sum(Tstar >= T0)) / (length(Tstar) + 1)
    
    list(pval=pval, beta_hat=beta_hat, sigma_hat=sigma_hat, B_eff=length(Tstar))
}

## ============================================================
## CIV (2D) - Adaptive GH cluster loglik for M2 fitted model (1 + z2 | id)
## ============================================================
post_mode_2d <- function(beta0, beta1, sd1, sd2, rho, y, x, z2, maxit=60) {
    G <- G_from_params(sd1, sd2, rho)
    Gi <- solve(G)
    b <- c(0,0)
    
    for (it in 1:maxit) {
        eta <- beta0 + beta1*x + b[1] + b[2]*z2
        mu  <- exp(eta)
        r   <- (y - mu)
        
        g1 <- sum(r * 1)  - as.numeric(Gi[1,] %*% b)
        g2 <- sum(r * z2) - as.numeric(Gi[2,] %*% b)
        g  <- c(g1, g2)
        
        S11 <- sum(mu)
        S12 <- sum(mu * z2)
        S22 <- sum(mu * z2 * z2)
        H <- -matrix(c(S11, S12, S12, S22), 2, 2, byrow=TRUE) - Gi
        
        step <- solve(H, g)
        b_new <- b - step
        if (max(abs(b_new - b)) < 1e-10) { b <- b_new; break }
        b <- b_new
    }
    
    cov <- solve(-H)
    list(bhat=b, cov=cov)
}

cluster_loglik_GH2_adapt <- function(phi, data_i, gh_nodes, gh_weights) {
    beta0 <- phi[1]
    beta1 <- phi[2]
    sd1   <- exp(phi[3])
    sd2   <- exp(phi[4])
    rho   <- tanh(phi[5])
    
    y  <- data_i$y
    x  <- data_i$x
    z2 <- data_i$z2
    
    ms <- post_mode_2d(beta0, beta1, sd1, sd2, rho, y, x, z2)
    bh <- ms$bhat
    C  <- ms$cov
    L  <- chol(C)  # upper
    
    grid <- expand.grid(t1=gh_nodes, t2=gh_nodes)
    w2 <- as.vector(outer(gh_weights, gh_weights))
    tmat <- cbind(grid$t1, grid$t2)
    
    bmat <- matrix(rep(bh, each=nrow(tmat)), ncol=2) + sqrt(2) * tmat %*% t(L)
    
    G <- G_from_params(sd1, sd2, rho)
    
    logLik_k <- apply(bmat, 1, function(bb) {
        b1 <- bb[1]; b2 <- bb[2]
        eta <- beta0 + beta1*x + b1 + b2*z2
        mu  <- exp(eta)
        sum(dpois(y, mu, log=TRUE)) +
            mvtnorm::dmvnorm(c(b1,b2), mean=c(0,0), sigma=G, log=TRUE)
    })
    
    quad_adj <- rowSums(tmat^2)
    log_jac <- log(2) + log(det(L))  # det(sqrt(2) L') = 2 det(L)
    
    log_terms <- log(w2) + logLik_k + quad_adj
    log_jac + logsumexp(log_terms)
}

compute_CIV_2d <- function(dat, phi_hat, K=11, ridge=1e-6) {
    dat$id <- as.factor(dat$id)
    ids <- levels(dat$id)
    ncl <- length(ids)
    
    gh <- statmod::gauss.quad(K, kind="hermite")
    
    scores <- vector("list", ncl)
    hessians <- vector("list", ncl)
    
    for (ii in seq_along(ids)) {
        di <- dat[dat$id == ids[ii], c("y","x","z2")]
        li_fun <- function(phi) cluster_loglik_GH2_adapt(phi, di, gh$nodes, gh$weights)
        si <- grad(li_fun, phi_hat)
        Hi <- -safe_hessian(li_fun, phi_hat)
        scores[[ii]] <- si
        hessians[[ii]] <- Hi
    }
    
    p <- length(phi_hat)
    Jn <- Reduce(`+`, lapply(scores, function(s) tcrossprod(s))) / ncl
    Hn <- Reduce(`+`, hessians) / ncl
    Rn <- ridge_invert(Hn, Jn, ridge=ridge)
    
    Ttr <- ncl * ((sum(diag(Rn))/p) - 1)^2
    Tld <- ncl * (logdet_spd(Rn))^2
    list(Ttr=Ttr, Tld=Tld)
}

bootstrap_CIV_2d <- function(dat, B=99, Kgh=11, stat=c("Ttr","Tld"), ridge=1e-6) {
    stat <- match.arg(stat)
    ctrl <- glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
    
    fit <- glmer(y ~ x + (1 + z2 | id), data=dat, family=poisson(link="log"), control=ctrl)
    beta_hat <- fixef(fit)
    
    vc <- VarCorr(fit)$id
    sd1 <- attr(vc, "stddev")[1]; sd1 <- max(sd1, 1e-6)
    sd2 <- attr(vc, "stddev")[2]; sd2 <- max(sd2, 1e-6)
    rho <- attr(vc, "correlation")[1,2]
    rho <- max(min(rho, 0.999), -0.999)
    
    phi_hat <- c(beta_hat[1], beta_hat[2], log(sd1), log(sd2), atanh(rho))
    T0 <- compute_CIV_2d(dat, phi_hat, K=Kgh, ridge=ridge)[[stat]]
    
    ids <- levels(as.factor(dat$id))
    ncl <- length(ids)
    Ghat <- G_from_params(sd1, sd2, rho)
    
    Tstar <- rep(NA_real_, B)
    
    for (b in 1:B) {
        bi <- MASS::mvrnorm(ncl, mu=c(0,0), Sigma=Ghat)
        rownames(bi) <- ids
        b1 <- bi[as.character(dat$id), 1]
        b2 <- bi[as.character(dat$id), 2]
        
        eta <- phi_hat[1] + phi_hat[2]*dat$x + b1 + b2*dat$z2
        dat_star <- dat
        dat_star$y <- rpois(nrow(dat), exp(eta))
        
        fit_star <- suppressWarnings(try(
            glmer(y ~ x + (1 + z2 | id), data=dat_star, family=poisson(link="log"), control=ctrl),
            silent=TRUE
        ))
        if (inherits(fit_star, "try-error")) next
        
        beta_s <- fixef(fit_star)
        vc_s <- VarCorr(fit_star)$id
        sd1_s <- attr(vc_s, "stddev")[1]; sd1_s <- max(sd1_s, 1e-6)
        sd2_s <- attr(vc_s, "stddev")[2]; sd2_s <- max(sd2_s, 1e-6)
        rho_s <- attr(vc_s, "correlation")[1,2]
        rho_s <- max(min(rho_s, 0.999), -0.999)
        
        phi_s <- c(beta_s[1], beta_s[2], log(sd1_s), log(sd2_s), atanh(rho_s))
        Tstar[b] <- compute_CIV_2d(dat_star, phi_s, K=Kgh, ridge=ridge)[[stat]]
    }
    
    Tstar <- Tstar[is.finite(Tstar)]
    pval <- (1 + sum(Tstar >= T0)) / (length(Tstar) + 1)
    
    list(pval=pval, beta_hat=beta_hat, sd1_hat=sd1, sd2_hat=sd2, rho_hat=rho, B_eff=length(Tstar))
}

## ============================================================
## PARALELIZACIÓN: plan
## ============================================================
n_workers <- 6
future::plan(future::multisession, workers = n_workers)

## ============================================================
## ONE REPLICATION (parallel unit)
## ============================================================
one_replication <- function(case, side, n, m,
                            beta0, beta1, sigma,
                            beta2, df,
                            sd1, sd2, rho,
                            Bboot1, Kgh1,
                            Bboot2, Kgh2,
                            stat, ridge,
                            seed_rep) {
    
    set.seed(seed_rep)
    
    dat <- switch(paste0(case, "_", side),
                  "M0_H0" = simulate_M0_H0(n, m, beta0, beta1, sigma),
                  "M0_H1" = simulate_M0_H1(n, m, beta0, beta1, beta2, sigma),
                  "M1_H0" = simulate_M1_H0(n, m, beta0, beta1, sigma),
                  "M1_H1" = simulate_M1_H1(n, m, beta0, beta1, df),
                  "M2_H0" = simulate_M2_H0(n, m, beta0, beta1, sd1, sd2, rho),
                  "M2_H1" = simulate_M2_H1(n, m, beta0, beta1, sd1),
                  stop("Combinación case/side no válida.")
    )
    
    if (case %in% c("M0", "M1")) {
        out <- try(bootstrap_CIV_1d(dat, B = Bboot1, Kgh = Kgh1, stat = stat, ridge = ridge), silent = TRUE)
        if (inherits(out, "try-error")) {
            return(list(ok = FALSE, pval = NA_real_,
                        beta0_hat=NA_real_, beta1_hat=NA_real_, sigma_hat=NA_real_))
        }
        return(list(ok = TRUE, pval = out$pval,
                    beta0_hat = out$beta_hat[1], beta1_hat = out$beta_hat[2], sigma_hat = out$sigma_hat))
    } else {
        out <- try(bootstrap_CIV_2d(dat, B = Bboot2, Kgh = Kgh2, stat = stat, ridge = ridge), silent = TRUE)
        if (inherits(out, "try-error")) {
            return(list(ok = FALSE, pval = NA_real_,
                        beta0_hat=NA_real_, beta1_hat=NA_real_,
                        sd1_hat=NA_real_, sd2_hat=NA_real_, rho_hat=NA_real_))
        }
        return(list(ok = TRUE, pval = out$pval,
                    beta0_hat=out$beta_hat[1], beta1_hat=out$beta_hat[2],
                    sd1_hat=out$sd1_hat, sd2_hat=out$sd2_hat, rho_hat=out$rho_hat))
    }
}

## ============================================================
## MONTE CARLO for one (case, side, n) in parallel
## ============================================================
mc_run_parallel <- function(case, side, n, m, R, alpha,
                            beta0, beta1, sigma,
                            beta2, df,
                            sd1, sd2, rho,
                            Bboot1, Kgh1,
                            Bboot2, Kgh2,
                            stat, ridge,
                            seed_base = 12345) {
    
    seeds <- seed_base + seq_len(R) +
        100000 * match(case, c("M0","M1","M2")) +
        1000   * (side=="H1") + n
    
    reps <- future_lapply(
        X = seq_len(R),
        FUN = function(r) {
            one_replication(case, side, n, m,
                            beta0, beta1, sigma,
                            beta2, df,
                            sd1, sd2, rho,
                            Bboot1, Kgh1, Bboot2, Kgh2,
                            stat, ridge,
                            seed_rep = seeds[r])
        },
        future.seed = TRUE
    )
    
    dt <- rbindlist(lapply(seq_along(reps), function(r) {
        as.list(c(list(case=case, side=side, n=n, m=m, r=r), reps[[r]]))
    }), fill = TRUE)
    
    dt[, reject := as.integer(pval <= alpha)]
    rej_rate <- dt[ok == TRUE, mean(reject, na.rm = TRUE)]
    n_ok <- dt[, sum(ok)]
    
    ## estimation metrics
    if (case %in% c("M0","M1")) {
        true <- c(beta0 = beta0, beta1 = beta1, sigma = sigma)
        est_ok <- dt[ok == TRUE, .(beta0_hat, beta1_hat, sigma_hat)]
        mean_est <- colMeans(est_ok, na.rm = TRUE)
        bias <- mean_est - true
        var_emp <- apply(est_ok, 2, var, na.rm = TRUE)
        rmse <- sqrt(colMeans((t(t(as.matrix(est_ok)) - true))^2, na.rm = TRUE))
        
        summ <- data.table(
            case=case, side=side, n=n, m=m, R=R, n_ok=n_ok, alpha=alpha, stat=stat,
            reject_rate=rej_rate,
            bias_beta0=bias["beta0_hat"], bias_beta1=bias["beta1_hat"], bias_sigma=bias["sigma_hat"],
            var_beta0=var_emp["beta0_hat"], var_beta1=var_emp["beta1_hat"], var_sigma=var_emp["sigma_hat"],
            rmse_beta0=rmse[1], rmse_beta1=rmse[2], rmse_sigma=rmse[3]
        )
    } else {
        true <- if (side=="H0") c(beta0=beta0, beta1=beta1, sd1=sd1, sd2=sd2, rho=rho) else c(beta0=beta0, beta1=beta1, sd1=sd1, sd2=0, rho=0)
        est_ok <- dt[ok == TRUE, .(beta0_hat, beta1_hat, sd1_hat, sd2_hat, rho_hat)]
        mean_est <- colMeans(est_ok, na.rm = TRUE)
        bias <- mean_est - true
        var_emp <- apply(est_ok, 2, var, na.rm = TRUE)
        rmse <- sqrt(colMeans((t(t(as.matrix(est_ok)) - true))^2, na.rm = TRUE))
        
        summ <- data.table(
            case=case, side=side, n=n, m=m, R=R, n_ok=n_ok, alpha=alpha, stat=stat,
            reject_rate=rej_rate,
            bias_beta0=bias["beta0_hat"], bias_beta1=bias["beta1_hat"],
            bias_sd1=bias["sd1_hat"], bias_sd2=bias["sd2_hat"], bias_rho=bias["rho_hat"],
            var_beta0=var_emp["beta0_hat"], var_beta1=var_emp["beta1_hat"],
            var_sd1=var_emp["sd1_hat"], var_sd2=var_emp["sd2_hat"], var_rho=var_emp["rho_hat"],
            rmse_beta0=rmse[1], rmse_beta1=rmse[2], rmse_sd1=rmse[3], rmse_sd2=rmse[4], rmse_rho=rmse[5]
        )
    }
    
    list(replications=dt, summary=summ)
}

## ============================================================
## RUN ONE CASE and WRITE independent files
## ============================================================
run_case_and_write <- function(case, sizes,
                               out_dir="Resultados_CIV",
                               m=10, R=200, alpha=0.05,
                               beta0=0.5, beta1=1.0, sigma=1.0,
                               beta2=0.3, df=3,
                               sd1=1.0, sd2=0.5, rho=0.2,
                               Bboot1=199, Kgh1=25,
                               Bboot2=99,  Kgh2=11,
                               stat="Ttr",
                               ridge=1e-6,
                               seed_base=12345) {
    
    dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
    
    f_sum   <- file.path(out_dir, paste0("results_", case, "_summary.csv"))
    f_rep   <- file.path(out_dir, paste0("results_", case, "_replications.csv"))
    f_log   <- file.path(out_dir, paste0("results_", case, "_log.txt"))
    f_final <- file.path(out_dir, paste0("results_", case, "_TypeI_Power_only.csv"))
    
    if (file.exists(f_sum)) file.remove(f_sum)
    if (file.exists(f_rep)) file.remove(f_rep)
    if (file.exists(f_log)) file.remove(f_log)
    if (file.exists(f_final)) file.remove(f_final)
    
    cat("LOG -", case, "\n", file=f_log, append=TRUE)
    cat("Workers:", n_workers, "\n", file=f_log, append=TRUE)
    
    all_summ <- list()
    
    for (n in sizes) {
        for (side in c("H0","H1")) {
            
            cat(sprintf("[%s] START case=%s side=%s n=%d\n", Sys.time(), case, side, n),
                file=f_log, append=TRUE)
            
            res <- mc_run_parallel(
                case=case, side=side, n=n, m=m, R=R, alpha=alpha,
                beta0=beta0, beta1=beta1, sigma=sigma,
                beta2=beta2, df=df,
                sd1=sd1, sd2=sd2, rho=rho,
                Bboot1=Bboot1, Kgh1=Kgh1,
                Bboot2=Bboot2, Kgh2=Kgh2,
                stat=stat, ridge=ridge,
                seed_base=seed_base
            )
            
            fwrite(res$replications, f_rep, append=file.exists(f_rep))
            fwrite(res$summary, f_sum, append=file.exists(f_sum))
            all_summ[[length(all_summ)+1]] <- res$summary
            
            cat(sprintf("[%s] END case=%s side=%s n=%d | reject=%.4f | n_ok=%d\n",
                        Sys.time(), case, side, n, res$summary$reject_rate, res$summary$n_ok),
                file=f_log, append=TRUE)
        }
    }
    
    sum_dt <- rbindlist(all_summ, fill=TRUE)
    typeI  <- sum_dt[side=="H0", .(TypeI=reject_rate), by=.(case,n,m,R,alpha,stat)]
    power  <- sum_dt[side=="H1", .(Power=reject_rate), by=.(case,n,m,R,alpha,stat)]
    final  <- merge(typeI, power, by=c("case","n","m","R","alpha","stat"), all=TRUE)
    
    fwrite(final, f_final)
    invisible(list(summary_by_side=sum_dt, typeI_power=final))
}

## ============================================================
## USER SETTINGS (edit here)
## ============================================================
## sampel size
sizes <- c(20, 100, 200, 400, 600, 1000, 2000, 3000)

## WARNING:
## If you truly need m=10000, M2 (2D) will be extremely heavy.
## Start with m=10 to validate; then scale carefully.
m <- 10

R <- 200
alpha <- 0.05

beta0 <- 0.5
beta1 <- 1.0
sigma <- 1.0

beta2 <- 0.3
df    <- 3

## M2-H0 parameters (when FIT is correct)
sd1 <- 1.0
sd2 <- 0.5
rho <- 0.2

## Bootstrap & GH
Bboot1 <- 199
Kgh1   <- 25

Bboot2 <- 99
Kgh2   <- 11

stat  <- "Ttr"  # "Tld" if you prefer
ridge <- 1e-6

out_dir <- "Resultados_CIV"

## ============================================================
## RUN ALL CASES (each produces independent files)
## ============================================================
cat('running M0 models...', '\n')
run_case_and_write("M0", sizes, out_dir,
                   m, R, alpha,
                   beta0, beta1, sigma,
                   beta2, df,
                   sd1, sd2, rho,
                   Bboot1, Kgh1,
                   Bboot2, Kgh2,
                   stat, ridge,
                   seed_base=50000)
cat('done!', '\n')


## ============================================================
## FIN
## ============================================================