mixed_fit <- function (y, X, Z, X_zi, Z_zi, id, offset, offset_zi, family, 
                       initial_values, Funs, control, penalized) {
    # Create lists of y, X, and Z per id
    y <- unattr(y); X <- unattr(X); Z <- unattr(Z); offset <- unattr(offset)
    X_zi <- unattr(X_zi); Z_zi <- unattr(Z_zi); offset_zi <- unattr(offset_zi)
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
    X_zi_lis <- lapply(id_unq, function (i) X_zi[id == i, , drop = FALSE])
    Z_zi_lis <- lapply(id_unq, function (i) Z_zi[id == i, , drop = FALSE])
    offset_zi_lis <- if (!is.null(offset_zi)) split(offset_zi, id)
    # Functions
    log_dens <- Funs$log_dens
    mu_fun <- Funs$mu_fun
    var_fun <- Funs$var_fun
    mu.eta_fun <- Funs$mu.eta_fun
    score_eta_fun <- Funs$score_eta_fun
    score_eta_zi_fun <- Funs$score_eta_zi_fun
    score_phis_fun <- Funs$score_phis_fun
    canonical <- !is.null(family$family) &&
        ((family$family == "binomial" && family$link == "logit") ||
             (family$family == "poisson" && family$link == "log"))
    known_families <- c("binomial", "poisson")
    user_defined <- !family$family %in% known_families
    numer_deriv <- if (control$numeric_deriv == "fd") fd else cd
    numer_deriv_vec <- if (control$numeric_deriv == "fd") fd_vec else cd_vec
    # dimensions
    n <- length(id_unq)
    ncx <- ncol(X)
    ncz <- ncol(Z)
    ncx_zi <- ncol(X_zi)
    ncz_zi <- ncol(Z_zi)
    nRE <- if (!is.null(Z_zi)) ncz + ncz_zi else ncz
    nAGQ <- control$nAGQ
    nAGQ_cartesian <- nAGQ^nRE
    ind_Z <- seq_len(ncol(Z))
    # initial values
    betas <- unname(initial_values[["betas"]])
    D <- unname(initial_values[["D"]])
    diag_D <- !is.matrix(D)
    if (diag_D) {
        D <- diag(D, nRE)
    }
    phis <- unname(initial_values[["phis"]])
    gammas <- unname(initial_values[["gammas"]])
    has_phis <- !is.null(phis)
    nparams <- length(betas) + length(if (diag_D) diag(D) else D[lower.tri(D, TRUE)]) + 
        length(phis) + length(gammas)
    post_modes <- matrix(0.0, n, nRE)
    # penalized components
    pen_mu <- if (penalized$penalized) rep(penalized$pen_mu, length.out = ncx)
    pen_invSigma <- if (penalized$penalized) 
        diag(rep(1 / penalized$pen_sigma^2, length.out = ncx), ncx)
    pen_df <- if (penalized$penalized) penalized$pen_df
    penalized <- penalized$penalized
    # set up EM algorithm
    iter_EM <- control$iter_EM
    update_GH <- seq(0, iter_EM, control$update_GH_every)
    tol1 <- control$tol1; tol2 <- control$tol2; tol3 <- control$tol3
    converged <- FALSE
    err_mgs <- paste("A large coefficient value has been detected during the optimization.\n",
                     "Please re-scale you covariates and/or try setting the control argument\n", 
                     "'iter_EM = 0'. Alternatively, this may due to a\n",
                     "divergence of the optimization algorithm, indicating that an overly\n",
                     "complex model is fitted to the data. For example, this could be\n",
                     "caused when including random-effects terms (e.g., in the\n", 
                     "zero-inflated part) that you do not need. Otherwise, adjust the\n",
                     "'max_coef_value' control argument.\n")
    large_shape_mgs <- paste("A value greater than 22000 has been detected for the shape/size\n",
                             "parameter of the negative binomial distribution. This typically\n",
                             "indicates that the Poisson model would be better. Otherwise,\n",
                             "adjust the 'max_phis_value' control argument.")
    if (iter_EM > 0) {
        Params <- matrix(0.0, iter_EM, nparams)
        GH <- GHfun(post_modes, y_lis, N_lis, X_lis, Z_lis, offset_lis, X_zi_lis, Z_zi_lis, 
                    offset_zi_lis, betas, solve(D), phis, gammas,
                    nAGQ, nRE, canonical, user_defined, Zty_lis, log_dens, mu_fun, var_fun,
                    mu.eta_fun, score_eta_fun, score_phis_fun, score_eta_zi_fun)
        b <- GH$b
        b2 <- GH$b2
        Ztb <- GH$Ztb
        Z_zitb <- GH$Z_zitb
        wGH <- GH$wGH
        dets <- GH$dets
        post_modes <- GH$post_modes
        lgLik <- numeric(iter_EM)
        for (it in seq_len(iter_EM)) {
            if (it %in% update_GH) {
                # calculate adaptive GH points and weights
                GH <- GHfun(post_modes, y_lis, N_lis, X_lis, Z_lis, offset_lis, X_zi_lis, 
                            Z_zi_lis, offset_zi_lis, betas, solve(D), phis, gammas,
                            nAGQ, nRE, canonical, user_defined, Zty_lis, log_dens, mu_fun, 
                            var_fun, mu.eta_fun, score_eta_fun, score_phis_fun, 
                            score_eta_zi_fun)
                b <- GH$b
                b2 <- GH$b2
                Ztb <- GH$Ztb
                Z_zitb <- GH$Z_zitb
                wGH <- GH$wGH
                dets <- GH$dets
                post_modes <- GH$post_modes
            }
            # save parameters
            Params[it, ] <- c(betas, if (diag_D) diag(D) else D[lower.tri(D, TRUE)], phis,
                              gammas)
            ##
            # calculate posterior distribution of the random effects
            eta_y <- as.vector(X %*% betas) + Ztb
            if (!is.null(offset))
                eta_y <- eta_y + offset
            eta_zi <- if (!is.null(X_zi)) as.vector(X_zi %*% gammas)
            if (!is.null(Z_zi))
                eta_zi <- eta_zi + Z_zitb
            if (!is.null(offset_zi))
                eta_zi <- eta_zi + offset_zi
            log_p_yb <- rowsum(log_dens(y, eta_y, mu_fun, phis, eta_zi), 
                               id, reorder = FALSE)
            log_p_b <- matrix(dmvnorm(b, rep(0, nRE), D, TRUE), n, nAGQ^nRE, byrow = TRUE)
            p_yb <- exp(log_p_yb + log_p_b)
            if (any(zero_ind <- p_yb == 0.0)) {
                p_yb[zero_ind] <- 1e-300
            }
            p_y <- c(p_yb %*% wGH)
            p_by <- p_yb / p_y
            t_p_by <- t(p_by)
            post_b <- apply(b, 2, function (b_k)
                colSums(t_p_by * matrix(b_k, nAGQ_cartesian, n) * wGH))
            post_b2 <- apply(b2, 2, function (b_k)
                colSums(t_p_by * matrix(b_k, nAGQ_cartesian, n) * wGH))
            # calculate log-likelihood
            log_p_y <- log(p_y * dets)
            lgLik[it] <- sum(log_p_y[is.finite(log_p_y)], na.rm = TRUE)
            if (penalized) {
                lgLik[it] <- lgLik[it] + dmvt(betas, mu = pen_mu, invSigma = pen_invSigma,
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
                    attr(converged, "during_EM") <- TRUE
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
                if (!is.null(gammas))
                    cat("gammas:", round(gammas, 4), "\n")
                cat("D:", round(if (diag_D) diag(D) else D[lower.tri(D, TRUE)], 4), "\n")
            }
            ############################
            # update parameters
            Dn <- matrix(colMeans(post_b2, na.rm = TRUE), nRE, nRE)
            D <- 0.5 * (Dn + t(Dn))
            if (diag_D) {
                D <- diag(diag(D), nRE)
            }
            if (has_phis) {
                Hphis <- numer_deriv_vec(phis, score_phis, y = y, X = X, betas = betas,
                                         Ztb = Ztb, offset = offset, eta_zi = eta_zi, 
                                         id = id, p_by = p_by, log_dens = log_dens, 
                                         mu_fun = mu_fun, wGH = wGH,
                                         score_phis_fun = score_phis_fun)
                Hphis <- nearPD(Hphis)
                scphis <- score_phis(phis, y, X, betas, Ztb, offset, eta_zi, id, p_by, 
                                     log_dens, mu_fun, wGH, score_phis_fun)
                phis <- phis - drop(solve(Hphis, scphis))
            }
            Hbetas <- numer_deriv_vec(betas, score_betas, y = y, N = N, X = X, id = id,
                                      offset = offset, phis = phis, Ztb = Ztb, eta_zi = eta_zi,
                                      p_by = p_by, wGH = wGH, canonical = canonical,
                                      user_defined = user_defined, Xty = Xty,
                                      log_dens = log_dens, mu_fun = mu_fun, var_fun = var_fun,
                                      mu.eta_fun = mu.eta_fun,
                                      score_eta_fun = score_eta_fun,
                                      score_phis_fun = score_phis_fun, 
                                      penalized = penalized, pen_mu = pen_mu, 
                                      pen_invSigma = pen_invSigma, pen_df = pen_df)
            Hbetas <- nearPD(Hbetas)
            scbetas <- score_betas(betas, y, N, X, id, offset, phis, Ztb, eta_zi, p_by, wGH,
                                   canonical, user_defined, Xty, log_dens, mu_fun, var_fun,
                                   mu.eta_fun, score_eta_fun, score_phis_fun,
                                   penalized, pen_mu, pen_invSigma, pen_df)
            betas <- betas - drop(solve(Hbetas, scbetas))
            if (!is.null(gammas)) {
                Hgammas <- numer_deriv_vec(gammas, score_gammas, y, X, betas, Ztb, offset, 
                                           X_zi, Z_zi, Z_zitb, offset_zi, log_dens,
                                           score_eta_zi_fun, phis, mu_fun, p_by, wGH, id)
                Hgammas <- nearPD(Hgammas)
                scgammas <- score_gammas(gammas, y, X, betas, Ztb, offset, X_zi, Z_zi, 
                                         Z_zitb, offset_zi, log_dens, score_eta_zi_fun, 
                                         phis, mu_fun, p_by, wGH, id)
                gammas <- gammas - drop(solve(Hgammas, scgammas))
            }
            if (any(abs(betas[-1L]) > control$max_coef_value) || 
                (!is.null(gammas) && any(abs(gammas) > control$max_coef_value))) {
                stop(err_mgs)
            }
            if (family$family %in% c("zero-inflated negative binomial", "negative binomial") &&
                exp(phis) > control$max_phis_value) {
                stop(large_shape_mgs)
            }
        }
    }
    list_thetas <- list(betas = betas, D = if (diag_D) log(diag(D)) else chol_transf(D))
    if (!is.null(phis)) {
        list_thetas <- c(list_thetas, list(phis = phis))
    }
    if (!is.null(gammas)) {
        list_thetas <- c(list_thetas, list(gammas = gammas))
    }
    tht <- unlist(as.relistable(list_thetas))
    if (!converged && control$iter_qN_outer > 0) {
        # start quasi-Newton iterations
        if (control$verbose) {
            cat("\nStart quasi-Newton iterations...\n\n")
        }
        length_notNA <- function (x) length(x[!is.na(x)])
        ns <- c("betas" = 0, "D" = 0, "phis" = 0, "gammas" = 0)
        lng <- sapply(list_thetas, length)
        ns[names(lng)] <- lng
        parscale <- rep(c(control$parscale_betas, control$parscale_D,
                          control$parscale_phis, control$parscale_gammas), ns)
       if (control$optimizer == "optimParallel") {
            cl <- parallel::makeCluster(2)
            parallel::setDefaultCluster(cl = cl)
            parallel::clusterExport(cl = cl, envir = environment(), 
                                    varlist = list("chol_transf", "deriv_D", "jacobian2",
                                                   "dmvt", "dmvnorm"))
        }
        for (it in seq_len(control$iter_qN_outer)) {
            GH <- GHfun(post_modes, y_lis, N_lis, X_lis, Z_lis, offset_lis, X_zi_lis, 
                        Z_zi_lis, offset_zi_lis, betas, solve(D), phis, gammas, nAGQ, nRE, 
                        canonical, user_defined, Zty_lis, log_dens, mu_fun, var_fun, 
                        mu.eta_fun, score_eta_fun, score_phis_fun, score_eta_zi_fun)
            opt <- optFun(tht, logLik_mixed, score_mixed, parscale = parscale,
                          control = control, id = id, y = y, N = N, X = X, Z = Z, 
                          offset = offset, X_zi = X_zi, Z_zi = Z_zi, offset_zi = offset_zi, 
                          GH = GH, canonical = canonical, user_defined = user_defined, 
                          Xty = Xty, log_dens = log_dens, mu_fun = mu_fun, 
                          var_fun = var_fun, mu.eta_fun = mu.eta_fun, 
                          score_eta_fun = score_eta_fun, score_eta_zi_fun = score_eta_zi_fun, 
                          score_phis_fun = score_phis_fun, list_thetas = list_thetas, 
                          diag_D = diag_D, penalized = penalized, pen_mu = pen_mu, 
                          pen_invSigma = pen_invSigma, pen_df = pen_df)
            tht <- opt$par
            new_pars <- relist(tht, skeleton = list_thetas)
            betas <- new_pars$betas
            phis <- new_pars$phis
            gammas <- new_pars$gammas
            D <- if (diag_D) diag(exp(new_pars$D), length(new_pars$D)) else chol_transf(new_pars$D)
            post_modes <- GH$post_modes
            if (any(abs(betas[-1L]) > control$max_coef_value) || 
                (!is.null(gammas) && any(abs(gammas) > control$max_coef_value))) {
                stop(err_mgs)
            }
            if (family$family %in% c("zero-inflated negative binomial", "negative binomial") &&
                exp(phis) > control$max_phis_value) {
                stop(large_shape_mgs)
            }
            if (opt$convergence == 0) {
                converged <- TRUE
                break
            }
            control$iter_qN <- control$iter_qN + control$iter_qN_incr
            if (control$verbose) cat("\n")
        }
        if (control$optimizer == "optimParallel") {
            parallel::stopCluster(cl)
        }
    }
    list_thetas <- list(betas = betas, D = if (diag_D) log(diag(D)) else chol_transf(D))
    if (!is.null(phis)) {
        list_thetas <- c(list_thetas, list(phis = phis))
    }
    if (!is.null(gammas)) {
        list_thetas <- c(list_thetas, list(gammas = gammas))
    }
    tht <- unlist(as.relistable(list_thetas))
    GH <- GHfun(post_modes, y_lis, N_lis, X_lis, Z_lis, offset_lis, X_zi_lis, Z_zi_lis, 
                offset_zi_lis, betas, solve(D), phis, gammas,
                nAGQ, nRE, canonical, user_defined, Zty_lis, log_dens, mu_fun, var_fun,
                mu.eta_fun, score_eta_fun, score_phis_fun, score_eta_zi_fun)
    logLik <- - logLik_mixed(tht, id, y, N, X, Z, offset, X_zi, Z_zi, offset_zi, GH, 
                             canonical, user_defined, Xty, log_dens, mu_fun, var_fun, 
                             mu.eta_fun, score_eta_fun, score_eta_zi_fun, score_phis_fun, 
                             list_thetas, diag_D, penalized, pen_mu, pen_invSigma, pen_df)
    Hessian <- cd_vec(tht, score_mixed, id = id, y = y, N = N, X = X, Z = Z, 
                      offset = offset, X_zi = X_zi, Z_zi = Z_zi, offset_zi = offset_zi, 
                      GH = GH, canonical = canonical, user_defined = user_defined, 
                      Xty = Xty, log_dens = log_dens, mu_fun = mu_fun, 
                      var_fun = var_fun, mu.eta_fun = mu.eta_fun, 
                      score_eta_fun = score_eta_fun, score_eta_zi_fun = score_eta_zi_fun, 
                      score_phis_fun = score_phis_fun, list_thetas = list_thetas, 
                      diag_D = diag_D, penalized = penalized, pen_mu = pen_mu, 
                      pen_invSigma = pen_invSigma, pen_df = pen_df)
    score_vect_contributions <- score_mixed(tht, id, y, N, X, Z, offset, X_zi, Z_zi, offset_zi, GH, 
                                            canonical, user_defined, Xty, log_dens, mu_fun, var_fun, 
                                            mu.eta_fun, score_eta_fun, score_eta_zi_fun, score_phis_fun, 
                                            list_thetas, diag_D, penalized, pen_mu, pen_invSigma, pen_df,
                                            i_contributions = TRUE)
    list(coefficients = betas, phis = if (has_phis) phis, D = D, gammas = gammas,
         post_modes = GH$post_modes, post_vars = GH$post_vars,
         logLik = logLik, Hessian = Hessian, 
         score_vect_contributions = score_vect_contributions,
         converged = converged)
}

optFun <- function (start, objective, gradient, parscale, control, ...) {
    if (control$optimizer == "optim") {
        optim(start, objective, gradient, method = control$optim_method,
              control = list(maxit = control$iter_qN, trace = 10 * control$verbose,
                             reltol = control$tol3, parscale = parscale), ...)
    } else if (control$optimizer == "optimParallel") {
        optimParallel::optimParallel(start, objective, gradient, method = control$optim_method,
                                     control = list(maxit = control$iter_qN, trace = 10 * control$verbose,
                                                    reltol = control$tol3, parscale = parscale), ...)
        
    } else {
        nlminb(start, objective, gradient, scale = 1 / parscale, 
               control = list(iter.max = control$iter_qN, trace = 10 * control$verbose,
                              rel.tol = control$tol3), ...)
    }
}

