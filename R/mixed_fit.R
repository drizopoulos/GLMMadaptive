mixed_fit <- function (y, X, Z, id, offset, family, initial_values, Funs, control,
                       penalized) {
    # Create lists of y, X, and Z per id
    y <- unattr(y); X <- unattr(X); Z <- unattr(Z); offset <- unattr(offset)
    id_unq <- unique(id)
    y_lis <- if (NCOL(y) == 2) lapply(id_unq, function (i) y[id == i, , drop = FALSE]) else split(y, id)
    N <- if (NCOL(y) == 2) y[, 1] + y[, 2]
    N_lis <- if (NCOL(y) == 2) split(N, id)
    X_lis <- lapply(id_unq, function (i) X[id == i, , drop = FALSE])
    Z_lis <- lapply(id_unq, function (i) Z[id == i, , drop = FALSE])
    offset_lis <- if (!is.null(offset)) split(offset, id)
    Zty_fun <- function (z, y) {
        if (NCOL(y) == 2) crossprod(z, y[, 1]) else crossprod(z, y)
    }
    Zty_lis <- lapply(mapply(Zty_fun, Z_lis, y_lis, SIMPLIFY = FALSE), drop)
    Xty <- drop(if (NCOL(y) == 2) crossprod(X, y[, 1]) else crossprod(X, y))
    # Functions
    log_dens <- Funs$log_dens
    mu_fun <- Funs$mu_fun
    var_fun <- Funs$var_fun
    mu.eta_fun <- Funs$mu.eta_fun
    score_eta_fun <- Funs$score_eta_fun
    score_phis_fun <- Funs$score_phis_fun
    canonical <- !is.null(family$family) &&
        ((family$family == "binomial" && family$link == "logit") ||
             (family$family == "poisson" && family$link == "log"))
    known_families <- c("binomial", "poisson", "negative binomial")
    user_defined <- !family$family %in% known_families
    numer_deriv <- if (control$numeric_deriv == "fd") fd else cd
    numer_deriv_vec <- if (control$numeric_deriv == "fd") fd_vec else cd_vec
    # dimensions
    n <- length(id_unq)
    ncx <- ncol(X)
    ncz <- ncol(Z)
    nAGQ <- control$nAGQ
    nAGQ_cartesian <- nAGQ^ncz
    # initial values
    betas <- unname(initial_values[["betas"]])
    D <- unname(initial_values[["D"]])
    diag_D <- !is.matrix(D)
    if (diag_D) {
        D <- diag(D, ncz)
    }
    phis <- unname(initial_values[["phis"]])
    has_phis <- !is.null(phis)
    nparams <- length(betas) + length(if (diag_D) diag(D) else D[lower.tri(D, TRUE)]) + length(phis)
    post_modes <- matrix(0.0, n, ncz)
    # penalized components
    pen_mu <- if (penalized$penalized) rep(penalized$pen_mu, ncx - 1)
    pen_invSigma <- if (penalized$penalized) diag(rep(1 / penalized$pen_sigma^2, ncx - 1), 
                                                  ncx - 1)
    pen_df <- if (penalized$penalized) penalized$pen_df
    penalized <- penalized$penalized
    # set up EM algorithm
    iter_EM <- control$iter_EM
    update_GH <- seq(0, iter_EM, control$update_GH_every)
    tol1 <- control$tol1; tol2 <- control$tol2; tol3 <- control$tol3
    converged <- FALSE
    if (iter_EM > 0) {
        Params <- matrix(0.0, iter_EM, nparams)
        GH <- GHfun(post_modes, y_lis, N_lis, X_lis, Z_lis, offset_lis, betas, solve(D), phis,
                    nAGQ, ncz, canonical, user_defined, Zty_lis, log_dens, mu_fun, var_fun,
                    mu.eta_fun, score_eta_fun, score_phis_fun)
        b <- GH$b
        b2 <- GH$b2
        Ztb <- GH$Ztb
        wGH <- GH$wGH
        dets <- GH$dets
        post_modes <- GH$post_modes
        lgLik <- numeric(iter_EM)
        for (it in seq_len(iter_EM)) {
            if (it %in% update_GH) {
                # calculate adaptive GH points and weights
                GH <- GHfun(post_modes, y_lis, N_lis, X_lis, Z_lis, offset_lis, betas, solve(D), phis,
                            nAGQ, ncz, canonical, user_defined, Zty_lis, log_dens, mu_fun,
                            var_fun, mu.eta_fun, score_eta_fun, score_phis_fun)
                b <- GH$b
                b2 <- GH$b2
                Ztb <- GH$Ztb
                wGH <- GH$wGH
                dets <- GH$dets
                post_modes <- GH$post_modes
            }
            # save parameters
            Params[it, ] <- c(betas, if (diag_D) diag(D) else D[lower.tri(D, TRUE)], phis)
            ##
            # calculate posterior distribution of the random effects
            log_p_yb <- rowsum(log_dens(y, as.vector(X %*% betas) + Ztb, mu_fun, phis), 
                               id, reorder = FALSE)
            log_p_b <- matrix(dmvnorm(b, rep(0, ncz), D, TRUE), n, nAGQ^ncz, byrow = TRUE)
            p_yb <- exp(log_p_yb + log_p_b)
            if (any(zero_ind <- p_yb == 0.0)) {
                p_yb[zero_ind] <- 1e-30
            }
            p_y <- c(p_yb %*% wGH)
            p_by <- p_yb / p_y
            t_p_by <- t(p_by)
            post_b <- apply(b, 2, function (b_k)
                colSums(t_p_by * matrix(b_k, nAGQ_cartesian, n) * wGH))
            post_b2 <- apply(b2, 2, function (b_k)
                colSums(t_p_by * matrix(b_k, nAGQ_cartesian, n) * wGH))
            post_vb <- post_b2 - if (ncz > 1) t(apply(post_b, 1, function (x) x %o% x)) else
                as.matrix(apply(post_b, 1, function (x) x %o% x))
            # calculate log-likelihood
            log_p_y <- log(p_y * dets)
            lgLik[it] <- sum(log_p_y[is.finite(log_p_y)], na.rm = TRUE)
            if (penalized) {
                lgLik[it] <- lgLik[it] + dmvt(betas[-1L], mu = pen_mu, invSigma = pen_invSigma,
                                              df = pen_df)
            }
            # check convergence
            if (it > 4 && lgLik[it] > lgLik[it - 1]) {
                thets1 <- Params[it - 1, ]
                thets2 <- Params[it, ]
                check1 <- max(abs(thets2 - thets1) / (abs(thets1) + tol1)) < tol2
                check2 <- (lgLik[it] - lgLik[it - 1]) < tol3 * (abs(lgLik[it - 1]) + tol3)
                if (check1 || check2) {
                    converged <- TRUE
                    if (control$verbose)
                        cat("\n\nconverged!\ncalculating Hessian...\n")
                    break
                }
            }
            # print results on screen
            if (control$verbose) {
                cat("\n\niter:", it, "\n")
                cat("log-likelihood:", lgLik[it], "\n")
                cat("betas:", round(betas, 4), "\n")
                if (has_phis)
                    cat("phis:", round(phis, 4), "\n")
                cat("D:", round(if (diag_D) diag(D) else D[lower.tri(D, TRUE)], 4), "\n")
            }
            ############################
            # update parameters
            Dn <- matrix(colMeans(post_b2, na.rm = TRUE), ncol(Z), ncol(Z))
            D <- 0.5 * (Dn + t(Dn))
            if (diag_D) {
                D <- diag(diag(D), ncz)
            }
            if (has_phis) {
                Hphis <- numer_deriv_vec(phis, score_phis, y = y, X = X, betas = betas,
                                         Ztb = Ztb, offset = offset, id = id, p_by = p_by,
                                         log_dens = log_dens, mu_fun = mu_fun, wGH = wGH,
                                         score_phis_fun = score_phis_fun)
                Hphis <- nearPD(Hphis)
                scphis <- score_phis(phis, y, X, betas, Ztb, offset, id, p_by, log_dens,
                                     mu_fun, wGH, score_phis_fun)
                phis <- phis - drop(solve(Hphis, scphis))
            }
            Hbetas <- numer_deriv_vec(betas, score_betas, y = y, N = N, X = X, id = id,
                                      offset = offset, phis = phis, Ztb = Ztb,
                                      p_by = p_by, wGH = wGH, canonical = canonical,
                                      user_defined = user_defined, Xty = Xty,
                                      log_dens = log_dens, mu_fun = mu_fun, var_fun = var_fun,
                                      mu.eta_fun = mu.eta_fun,
                                      score_eta_fun = score_eta_fun,
                                      score_phis_fun = score_phis_fun, 
                                      penalized = penalized, pen_mu = pen_mu, 
                                      pen_invSigma = pen_invSigma, pen_df = pen_df)
            Hbetas <- nearPD(Hbetas)
            scbetas <- score_betas(betas, y, N, X, id, offset, phis, Ztb, p_by, wGH,
                                   canonical, user_defined, Xty, log_dens, mu_fun, var_fun,
                                   mu.eta_fun, score_eta_fun, score_phis_fun,
                                   penalized, pen_mu, pen_invSigma, pen_df)
            betas <- betas - drop(solve(Hbetas, scbetas))
        }
    }
    list_thetas <- list(betas = betas, D = if (diag_D) log(diag(D)) else chol_transf(D),
                        phis = if (is.null(phis)) NA else phis)
    tht <- unlist(as.relistable(list_thetas))
    tht <- tht[!is.na(tht)]
    if (!converged) {
        # start quasi-Newton iterations
        if (control$verbose) {
            cat("\nStart quasi-Newton iterations...\n\n")
        }
        length_notNA <- function (x) length(x[!is.na(x)])
        ctrl <- list(maxit = control$iter_qN, trace = 10 * control$verbose,
                     reltol = control$tol3,
                     parscale = rep(c(control$parscale_betas, control$parscale_D,
                                      control$parscale_phis),
                                    sapply(list_thetas, length_notNA)))
        for (it in seq_len(control$iter_qN_outer)) {
            GH <- GHfun(post_modes, y_lis, N_lis, X_lis, Z_lis, offset_lis, betas, solve(D), phis,
                        nAGQ, ncz, canonical, user_defined, Zty_lis, log_dens, mu_fun,
                        var_fun, mu.eta_fun, score_eta_fun, score_phis_fun)
            # penalty include here
            opt <- optim(tht, logLik_mixed, score_mixed, method = control$optim_method,
                         control = ctrl, id = id, y = y, N = N, X = X, Z = Z, offset = offset,
                         phis = phis, Ztb = Ztb, GH = GH, canonical = canonical,
                         user_defined = user_defined, Xty = Xty, log_dens = log_dens,
                         mu_fun = mu_fun, var_fun = var_fun, mu.eta_fun = mu.eta_fun,
                         score_eta_fun = score_eta_fun, score_phis_fun = score_phis_fun,
                         list_thetas = list_thetas, diag_D = diag_D, 
                         penalized = penalized, pen_mu = pen_mu, 
                         pen_invSigma = pen_invSigma, pen_df = pen_df)
            new_pars <- relist(opt$par, skeleton = list_thetas)
            betas <- new_pars$betas
            phis <- new_pars$phis
            D <- if (diag_D) diag(exp(new_pars$D), length(new_pars$D)) else chol_transf(new_pars$D)
            post_modes <- GH$post_modes
            if (opt$convergence == 0) {
                converged <- TRUE
                break
            }
            ctrl$maxit <- ctrl$maxit + control$iter_qN_incr
            if (control$verbose) cat("\n")
        }
    }
    list_thetas <- list(betas = betas, D = if (diag_D) log(diag(D)) else chol_transf(D),
                        phis = if (is.null(phis)) NA else phis)
    tht <- unlist(as.relistable(list_thetas))
    tht <- tht[!is.na(tht)]
    GH <- GHfun(post_modes, y_lis, N_lis, X_lis, Z_lis, offset_lis, betas, solve(D), phis,
                nAGQ, ncz, canonical, user_defined, Zty_lis, log_dens, mu_fun, var_fun,
                mu.eta_fun, score_eta_fun, score_phis_fun)
    logLik <- - logLik_mixed(tht, id, y, N, X, Z, offset, phis, Ztb, GH, canonical,
                             user_defined, Xty, log_dens, mu_fun, var_fun, mu.eta_fun,
                             score_eta_fun, score_phis_fun, list_thetas, diag_D,
                             penalized, pen_mu, pen_invSigma, pen_df)
    Hessian <- cd_vec(tht, score_mixed, id = id, y = y, N = N, X = X, Z = Z, offset = offset,
                      phis = phis, Ztb = Ztb, GH = GH, list_thetas = list_thetas,
                      canonical = canonical, user_defined = user_defined, Xty = Xty,
                      log_dens = log_dens, mu_fun = mu_fun, var_fun = var_fun,
                      mu.eta_fun = mu.eta_fun, score_eta_fun = score_eta_fun,
                      score_phis_fun = score_phis_fun, diag_D = diag_D,
                      penalized = penalized, pen_mu = pen_mu, 
                      pen_invSigma = pen_invSigma, pen_df = pen_df)
    list(coefficients = betas, phis = if (has_phis) phis, D = D,
         post_modes = GH$post_modes, post_vars = GH$post_vars,
         logLik = logLik, Hessian = Hessian,
         converged = converged)
}

