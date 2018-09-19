logLik_mixed <- function (thetas, id, y, N, X, Z, offset, X_zi, Z_zi, offset_zi, GH, 
                          canonical, user_defined, Xty, log_dens, mu_fun, var_fun, 
                          mu.eta_fun, score_eta_fun, score_eta_zi_fun, score_phis_fun, 
                          list_thetas, diag_D, penalized, pen_mu, pen_invSigma, pen_df) {
    thetas <- relist(thetas, skeleton = list_thetas)
    betas <- thetas$betas
    phis <- thetas$phis
    gammas <- thetas$gammas
    D <- if (diag_D) diag(exp(thetas$D), length(thetas$D)) else chol_transf(thetas$D)
    nRE <- ncol(D)
    ##
    b <- GH$b
    Ztb <- GH$Ztb
    Z_zitb <- GH$Z_zitb
    wGH <- GH$wGH
    dets <- GH$dets
    ##
    eta_y <- as.vector(X %*% betas) + Ztb
    if (!is.null(offset))
        eta_y <- eta_y + offset
    eta_zi <- if (!is.null(X_zi)) as.vector(X_zi %*% gammas)
    if (!is.null(Z_zi))
        eta_zi <- eta_zi + Z_zitb
    if (!is.null(offset_zi))
        eta_zi <- eta_zi + offset_zi
    log_p_yb <- rowsum(log_dens(y, eta_y, mu_fun, phis, eta_zi), id, reorder = FALSE)
    log_p_b <- matrix(dmvnorm(b, rep(0, nRE), D, TRUE),
                      nrow(log_p_yb), ncol(log_p_yb), byrow = TRUE)
    # log penalty include here dmvt(betas, pen_mean, invSigma = pen_invsds, df = pen_df)
    p_yb <- exp(log_p_yb + log_p_b)
    if (any(zero_ind <- p_yb == 0.0)) {
        p_yb[zero_ind] <- 1e-300
    }
    p_y <- c(p_yb %*% wGH) * dets
    out <- - sum(log(p_y), na.rm = TRUE)
    if (penalized)
        out <- out - dmvt(betas[-1L], mu = pen_mu, invSigma = pen_invSigma, df = pen_df)
    out
}

score_mixed <- function (thetas, id, y, N, X, Z, offset, X_zi, Z_zi, offset_zi, GH, 
                         canonical, user_defined, Xty, log_dens, mu_fun, var_fun, 
                         mu.eta_fun, score_eta_fun, score_eta_zi_fun, score_phis_fun, 
                         list_thetas, diag_D, penalized, pen_mu, pen_invSigma, pen_df) {
    thetas <- relist(thetas, skeleton = list_thetas)
    betas <- thetas$betas
    phis <- thetas$phis
    gammas <- thetas$gammas
    D <- if (diag_D) diag(exp(thetas$D), length(thetas$D)) else chol_transf(thetas$D)
    nRE <- ncol(D)
    ##
    b <- GH$b
    b2 <- GH$b2
    Ztb <- GH$Ztb
    Z_zitb <- GH$Z_zitb
    wGH <- GH$wGH
    ##
    eta_y <- as.vector(X %*% betas) + Ztb
    if (!is.null(offset))
        eta_y <- eta_y + offset
    eta_zi <- if (!is.null(X_zi)) as.vector(X_zi %*% gammas)
    if (!is.null(Z_zi))
        eta_zi <- eta_zi + Z_zitb
    if (!is.null(offset_zi))
        eta_zi <- eta_zi + offset_zi
    log_Lik <- log_dens(y, eta_y, mu_fun, phis, eta_zi)
    log_p_yb <- rowsum(log_Lik, id, reorder = FALSE)
    log_p_b <- matrix(dmvnorm(b, rep(0, nRE), D, TRUE),
                      nrow(log_p_yb), ncol(log_p_yb), byrow = TRUE)
    # log penalty include here dmvt(betas, pen_mean, invSigma = pen_invsds, df = pen_df)
    p_yb <- exp(log_p_yb + log_p_b)
    if (any(zero_ind <- p_yb == 0.0)) {
        p_yb[zero_ind] <- 1e-300
    }
    p_y <- c(p_yb %*% wGH)
    p_by <- p_yb / p_y
    t_p_by <- t(p_by)
    n <- length(p_y)
    post_b <- apply(b, 2, function (b_k) colSums(t_p_by * matrix(b_k, ncol(Ztb), n) * wGH))
    post_b2 <- apply(b2, 2, function (b_k) colSums(t_p_by * matrix(b_k, ncol(Ztb), n) * wGH))
    post_vb <- post_b2 - if (nRE > 1) t(apply(post_b, 1, function (x) x %o% x)) else
        as.matrix(apply(post_b, 1, function (x) x %o% x))
    ###
    mu_y <- if (!is.null(attr(log_Lik, "mu_y"))) attr(log_Lik, "mu_y") else mu_fun(eta_y)
    score.betas <- if (user_defined) {
        if (!is.null(score_eta_fun)) {
            z <- score_eta_fun(y, mu_y, phis, eta_zi)
            ncx <- ncol(X)
            sc <- numeric(ncx)
            for (l in seq_len(ncx)) {
                cc <- drop(rowsum(X[, l] * z, id, reorder = FALSE))
                sc[l] <- sum(c((cc * p_by) %*% wGH), na.rm = TRUE)
            }
            - sc
        } else {
            l1 <- log_dens(y, eta_y + 1e-04, mu_fun, phis, eta_zi)
            l2 <- log_dens(y, eta_y - 1e-04, mu_fun, phis, eta_zi)
            z <- (l1 - l2) / (2 * 1e-04)
            ncx <- ncol(X)
            sc <- numeric(ncx)
            for (l in seq_len(ncx)) {
                cc <- drop(rowsum(X[, l] * z, id, reorder = FALSE))
                sc[l] <- sum(c((cc * p_by) %*% wGH), na.rm = TRUE)
            }
            - sc
        }
    } else {
        if (canonical) {
            if (!is.null(N))
                mu_y <- mu_y * N
            ncx <- ncol(X)
            sc <- numeric(ncx)
            for (l in seq_len(ncx)) {
                cc <- drop(rowsum(X[, l] * mu_y, id, reorder = FALSE))
                sc[l] <- sum(c((cc * p_by) %*% wGH), na.rm = TRUE)
            }
            - Xty + sc
        } else {
            var <- var_fun(mu_y)
            deriv <- mu.eta_fun(eta_y)
            z <- if (!is.null(N)) (y[, 1] - N * mu_y) * deriv / var else (y - mu_y) * deriv / var
            ncx <- ncol(X)
            sc <- numeric(ncx)
            for (l in seq_len(ncx)) {
                cc <- drop(rowsum(X[, l] * z, id, reorder = FALSE))
                sc[l] <- sum(c((cc * p_by) %*% wGH))
            }
            - sc
        }
    }
    if (penalized) {
        pen_invSigma_betas <- betas[-1L] * diag(pen_invSigma) / pen_df
        fact <- (pen_df + ncx) / c(1 + crossprod(betas[-1L], pen_invSigma_betas))
        score.betas <- score.betas + c(0, pen_invSigma_betas * fact)
    }
    ###
    score.phis <- if (!is.null(phis)) {
        if (is.null(score_phis_fun)) {
            n_phis <- length(phis)
            sc <- numeric(n_phis)
            for (i in seq_len(n_phis)) {
                phis1 <- phis2 <- phis
                phis1[i] <- phis[i] + 1e-03
                phis2[i] <- phis[i] - 1e-03
                l1 <- log_dens(y, eta_y, mu_fun, phis1, eta_zi)
                l2 <- log_dens(y, eta_y, mu_fun, phis2, eta_zi)
                z <- (l1 - l2) / (phis1[i] - phis2[i])
                sc[i] <- sum(c((rowsum(z, id, reorder = FALSE) * p_by) %*% wGH), na.rm = TRUE)
            }
            - sc
        } else {
            z <- score_phis_fun(y, mu_y, phis, eta_zi)
            -sum(c((rowsum(z, id, reorder = FALSE) * p_by) %*% wGH), na.rm = TRUE)
        }
    }
    score.gammas <- if (!is.null(X_zi)) {
        z <- if (!is.null(score_eta_zi_fun)) {
            score_eta_zi_fun(y, mu_y, phis, eta_zi)
        } else {
            l1 <- log_dens(y, eta_y, mu_fun, phis, eta_zi + 1e-03)
            l2 <- log_dens(y, eta_y, mu_fun, phis, eta_zi - 1e-03)
            (l1 - l2) / (2 * 1e-03)
        }
        ncx_zi <- ncol(X_zi)
        sc <- numeric(ncx_zi)
        for (l in seq_len(ncx_zi)) {
            cc <- drop(rowsum(X_zi[, l] * z, id, reorder = FALSE))
            sc[l] <- sum(c((cc * p_by) %*% wGH), na.rm = TRUE)
        }
        - sc
    }
    ###
    score.D <- if (diag_D) {
        D <- diag(D)
        svD <- 1/D
        svD2 <- svD^2
        cS.postVB <- colSums(as.matrix(post_vb), na.rm = TRUE)
        dim(cS.postVB) <- c(nRE, nRE)
        D * 0.5 * (n * svD - diag(cS.postVB) * svD2 - 
                       colSums(as.matrix(post_b^2), na.rm = TRUE) * svD2)
    } else {
        svD <- solve(D)
        dD <- deriv_D(D)
        ndD <- length(dD)
        D1 <- sapply(dD, function (x) sum(svD * x))
        D2 <- t(sapply(dD, function (x) c(svD %*% x %*% svD)))
        cS.postVB <- colSums(as.matrix(post_vb), na.rm = TRUE)
        out <- numeric(ndD)
        for (i in seq_along(dD)) {
            D.mat <- D2[i, ]
            dim(D.mat) <- c(nRE, nRE)
            out[i] <- sum(D2[i, ] * cS.postVB, na.rm = TRUE) +
                sum((post_b %*% D.mat) * post_b, na.rm = TRUE)
        }
        J <- jacobian2(attr(D, "L"), nRE)
        drop(0.5 * (n * D1 - out) %*% J)
    }
    ###
    c(score.betas, score.D, score.phis, score.gammas)
}

score_betas <- function (betas, y, N, X, id, offset, phis, Ztb, eta_zi, p_by, wGH, canonical,
                         user_defined, Xty, log_dens, mu_fun, var_fun, mu.eta_fun,
                         score_eta_fun, score_phis_fun, penalized, pen_mu, pen_invSigma, 
                         pen_df) {
    eta_y <- as.vector(X %*% betas) + Ztb
    if (!is.null(offset))
        eta_y <- eta_y + offset
    mu_y <- mu_fun(eta_y)
    ncx <- ncol(X)
    sc <- numeric(ncx)
    out <- if (user_defined) {
        if (!is.null(score_eta_fun)) {
            z <- score_eta_fun(y, mu_y, phis, eta_zi)
            for (l in seq_len(ncx)) {
                cc <- drop(rowsum(X[, l] * z, id, reorder = FALSE))
                sc[l] <- sum(c((cc * p_by) %*% wGH))
            }
            - sc
        } else {
            l1 <- log_dens(y, eta_y + 1e-05, mu_fun, phis, eta_zi)
            l2 <- log_dens(y, eta_y - 1e-05, mu_fun, phis, eta_zi)
            z <- (l1 - l2) / (2 * 1e-05)
            for (l in seq_len(ncx)) {
                cc <- drop(rowsum(X[, l] * z, id, reorder = FALSE))
                sc[l] <- sum(c((cc * p_by) %*% wGH))
            }
            - sc
        }
    } else {
        if (canonical) {
            if (!is.null(N))
                mu_y <- N * mu_y
            sc <- numeric(ncx)
            for (l in seq_len(ncx)) {
                cc <- rowsum(X[, l] * mu_y, id, reorder = FALSE)
                sc[l] <- sum(c((cc * p_by) %*% wGH))
            }
            - Xty + sc
        } else {
            var <- var_fun(mu_y)
            deriv <- mu.eta_fun(eta_y)
            z <- if (!is.null(N)) (y[, 1] - N * mu_y) * deriv / var else (y - mu_y) * deriv / var
            for (l in seq_len(ncx)) {
                cc <- rowsum(X[, l] * z, id, reorder = FALSE)
                sc[l] <- sum(c((cc * p_by) %*% wGH), na.rm = TRUE)
            }
            - sc
        }
    }
    if (penalized) {
        pen_invSigma_betas <- betas[-1L] * diag(pen_invSigma) / pen_df
        fact <- (pen_df + ncx) / c(1 + crossprod(betas[-1L], pen_invSigma_betas))
        out <- out + c(0, pen_invSigma_betas * fact)
    }
    out
}

score_phis <- function (phis, y, X, betas, Ztb, offset, eta_zi, id, p_by,
                        log_dens, mu_fun, wGH, score_phis_fun) {
    eta_y <- as.vector(X %*% betas) + Ztb
    if (!is.null(offset))
        eta_y <- eta_y + offset
    if (is.null(score_phis_fun)) {
        n_phis <- length(phis)
        sc <- numeric(n_phis)
        for (i in seq_len(n_phis)) {
            phis1 <- phis2 <- phis
            phis1[i] <- phis1[i] + 1e-03
            phis2[i] <- phis2[i] - 1e-03
            l1 <- log_dens(y, eta_y, mu_fun, phis1, eta_zi)
            l2 <- log_dens(y, eta_y, mu_fun, phis2, eta_zi)
            z <- (l1 - l2) / (phis1[i] - phis2[i])
            sc[i] <- sum(c((rowsum(z, id, reorder = FALSE) * p_by) %*% wGH), na.rm = TRUE)
        }
        - sc
    } else {
        mu_y <- mu_fun(eta_y)
        z <- score_phis_fun(y, mu_y, phis, eta_zi)
        -sum(c((rowsum(z, id, reorder = FALSE) * p_by) %*% wGH), na.rm = TRUE)
    }
}

score_gammas <- function (gammas, y, X, betas, Ztb, offset, X_zi, Z_zi, Z_zitb, offset_zi,
                          log_dens, score_eta_zi_fun, phis, mu_fun, p_by, wGH, id) {
    eta_y <- as.vector(X %*% betas) + Ztb
    if (!is.null(offset))
        eta_y <- eta_y + offset
    eta_zi <- as.vector(X_zi %*% gammas)
    if (!is.null(Z_zi))
        eta_zi <- eta_zi + Z_zitb
    if (!is.null(offset_zi))
        eta_zi <- eta_zi + offset_zi
    mu_y <- mu_fun(eta_y)
    z <- if (!is.null(score_eta_zi_fun)) {
        score_eta_zi_fun(y, mu_y, phis, eta_zi)
    } else {
        l1 <- log_dens(y, eta_y, mu_fun, phis, eta_zi + 1e-03)
        l2 <- log_dens(y, eta_y, mu_fun, phis, eta_zi - 1e-03)
        (l1 - l2) / (2 * 1e-03)
    }
    ncx_zi <- ncol(X_zi)
    sc <- numeric(ncx_zi)
    for (l in seq_len(ncx_zi)) {
        cc <- drop(rowsum(X_zi[, l] * z, id, reorder = FALSE))
        sc[l] <- sum(c((cc * p_by) %*% wGH), na.rm = TRUE)
    }
    - sc
}

binomial_log_dens = function (y, eta, mu_fun, phis, eta_zi) {
    mu_y <- mu_fun(eta)
    out <- if (NCOL(y) == 2) {
        dbinom(y[, 1], y[, 1] + y[, 2], mu_y, TRUE)
    } else {
        dbinom(y, 1, mu_y, TRUE)
    }
    attr(out, "mu_y") <- mu_y
    out
}

poisson_log_dens = function (y, eta, mu_fun, phis, eta_zi) {
    mu_y <- mu_fun(eta)
    out <- y * log(mu_y) - mu_y - lgamma(y + 1)
    attr(out, "mu_y") <- mu_y
    out
}

negative.binomial_log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
    phis <- exp(phis)
    mu <- mu_fun(eta)
    log_mu_phis <- log(mu + phis)
    comp1 <- lgamma(y + phis) - lgamma(phis) - lgamma(y + 1)
    comp2 <- phis * log(phis) - phis * log_mu_phis
    comp3 <- y * log(mu) - y * log_mu_phis
    out <- comp1 + comp2 + comp3
    attr(out, "mu_y") <- mu
    out
}

negative.binomial <- function () {
    stats <- make.link("log")
    stats <- make.link(link = "log")
    log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
        # the log density function
        phis <- exp(phis)
        mu <- mu_fun(eta)
        log_mu_phis <- log(mu + phis)
        comp1 <- lgamma(y + phis) - lgamma(phis) - lgamma(y + 1)
        comp2 <- phis * log(phis) - phis * log_mu_phis
        comp3 <- y * log(mu) - y * log_mu_phis
        out <- comp1 + comp2 + comp3
        attr(out, "mu_y") <- mu
        out
    }
    score_eta_fun <- function (y, mu, phis, eta_zi) {
        # the derivative of the log density w.r.t. mu
        phis <- exp(phis)
        mu_phis <- mu + phis
        comp2 <- - phis / mu_phis
        comp3 <- y / mu - y / mu_phis
        # the derivative of mu w.r.t. eta (this depends on the chosen link function)
        mu.eta <- mu
        (comp2 + comp3) * mu.eta
    }
    score_phis_fun <- function (y, mu, phis, eta_zi) {
        # the derivative of the log density w.r.t. phis
        phis <- exp(phis)
        mu_phis <- mu + phis
        comp1 <- digamma(y + phis) - digamma(phis)
        comp2 <- log(phis) + 1 - log(mu_phis) - phis / mu_phis
        comp3 <- - y / mu_phis
        (comp1 + comp2 + comp3) * phis
    }
    structure(list(family = "negative binomial", link = stats$name, 
                   linkfun = stats$linkfun, linkinv = stats$linkinv, log_dens = log_dens,
                   score_eta_fun = score_eta_fun, score_phis_fun = score_phis_fun),
              class = "family")
}

zi.poisson <- function () {
    stats <- make.link(link = "log")
    log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
        # the log density function
        ind_y0 <- y == 0
        ind_y1 <- y > 0
        mu <- as.matrix(mu_fun(eta))
        lambda <- as.matrix(exp(eta_zi))
        mu0 <- mu[ind_y0, ]
        lambda0 <- lambda[ind_y0, ]
        mu1 <- mu[ind_y1, ]
        lambda1 <- lambda[ind_y1, ]
        out <- as.matrix(eta)
        out[ind_y0, ] <- log(lambda0 + exp(-mu0))
        out[ind_y1, ] <- y[ind_y1] * log(mu1) - mu1 - lgamma(y[ind_y1] + 1)
        out <- out - log(1 + drop(lambda))
        attr(out, "mu_y") <- mu
        out
    }
    score_eta_fun <- function (y, mu, phis, eta_zi) {
        ind_y0 <- y == 0
        ind_y1 <- y > 0
        mu <- as.matrix(mu)
        lambda <- as.matrix(exp(eta_zi))
        mu0 <- mu[ind_y0, ]
        lambda0 <- lambda[ind_y0, ]
        mu1 <- mu[ind_y1, ]
        out <- mu
        out[ind_y0, ] <- - mu0 / (lambda0 * exp(mu0) + 1)
        out[ind_y1, ] <- y[ind_y1] - mu1
        out
    }
    score_eta_zi_fun <- function (y, mu, phis, eta_zi) {
        ind_y0 <- y == 0
        mu <- as.matrix(mu)
        eta_zi <- as.matrix(eta_zi)
        lambda <- exp(eta_zi)
        lambda0_exp_mu0 <- exp(eta_zi[ind_y0, ] + mu[ind_y0, ])
        lambda0_exp_mu0[lambda0_exp_mu0 == Inf] <- 1e200
        out <- matrix(- lambda / (1 + lambda), nrow = nrow(mu), ncol = ncol(mu))
        out[ind_y0, ] <- out[ind_y0, ] + 
            drop(lambda0_exp_mu0 / (lambda0_exp_mu0 + 1))
        out
    }
    structure(list(family = "zero-inflated poisson", link = stats$name, 
                   linkfun = stats$linkfun, linkinv = stats$linkinv, log_dens = log_dens,
                   score_eta_fun = score_eta_fun,
                   score_eta_zi_fun = score_eta_zi_fun),
              class = "family")
}

zi.negative.binomial <- function () {
    stats <- make.link(link = "log")
    log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
        # the log density function
        # NB part
        phis <- exp(phis)
        mu <- mu_fun(eta)
        log_mu_phis <- log(mu + phis)
        comp1 <- lgamma(y + phis) - lgamma(phis) - lgamma(y + 1)
        comp2 <- phis * log(phis) - phis * log_mu_phis
        comp3 <- y * log(mu) - y * log_mu_phis
        out <- as.matrix(comp1 + comp2 + comp3)
        # ZI part
        ind_y0 <- y == 0
        ind_y1 <- y > 0
        pis <- as.matrix(plogis(eta_zi))
        # combined
        out[ind_y0, ] <- log(pis[ind_y0, ] + (1 - pis[ind_y0, ]) * exp(out[ind_y0, ]))
        out[ind_y1, ] <- log(1 - pis[ind_y1, ]) + out[ind_y1, ]
        attr(out, "mu_y") <- mu
        out
    }
    score_eta_fun <- function (y, mu, phis, eta_zi) {
        # NB part
        phis <- exp(phis)
        mu <- as.matrix(mu)
        mu_phis <- mu + phis
        comp2 <- - phis / mu_phis
        comp3 <- y / mu - y / mu_phis
        mu.eta <- mu
        out <- (comp2 + comp3) * mu.eta
        # ZI part
        ind_y0 <- y == 0
        lambda <- as.matrix(exp(eta_zi))
        t <- phis / (phis + mu[ind_y0, ])
        den <- (lambda[ind_y0, ] + t^phis) * (phis + mu[ind_y0, ])^2
        out[ind_y0, ] <- - phis^2 * t^(phis - 1) * mu[ind_y0, ] / den
        out
    }
    score_eta_zi_fun <- function (y, mu, phis, eta_zi) {
        phis <- exp(phis)
        ind_y0 <- y == 0
        ind_y1 <- y > 0
        # NB part
        mu <- as.matrix(mu)
        lambda <- as.matrix(exp(eta_zi))
        out <- mu
        out[ind_y1, ] <- - lambda[ind_y1, ] / (1 + lambda[ind_y1, ])
        # ZI part
        t <- phis / (phis + mu[ind_y0, ])
        out[ind_y0, ] <- lambda[ind_y0, ] / (lambda[ind_y0, ] + t^phis) - 
            lambda[ind_y0, ] / (1 + lambda[ind_y0, ])
        out
    }
    score_phis_fun <- function (y, mu, phis, eta_zi) {
        # NB part
        phis <- exp(phis)
        mu <- as.matrix(mu)
        mu_phis <- mu + phis
        comp1 <- digamma(y + phis) - digamma(phis)
        comp2 <- log(phis) + 1 - log(mu_phis) - phis / mu_phis
        comp3 <- - y / mu_phis
        out <- (comp1 + comp2 + comp3) * phis
        # ZI part
        ind_y0 <- y == 0
        lambda <- as.matrix(exp(eta_zi))
        t <- phis / (phis + mu[ind_y0, ])
        t_phis <- t^phis
        out[ind_y0, ] <- t_phis * (log(t) + 1 - t) * phis / (lambda[ind_y0, ] + t_phis)
        out
    }
    structure(list(family = "zero-inflated negative binomial", link = stats$name, 
                   linkfun = stats$linkfun, linkinv = stats$linkinv, log_dens = log_dens,
                   score_eta_fun = score_eta_fun,
                   score_eta_zi_fun = score_eta_zi_fun,
                   score_phis_fun = score_phis_fun),
              class = "family")
}

hurdle.poisson <- function () {
    stats <- make.link("log")
    log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
        # binary indicator for y > 0
        ind <- y > 0
        # non-zero part
        eta <- as.matrix(eta)
        mu <- mu_fun(eta)
        eta_zi <- as.matrix(eta_zi)
        out <- eta
        out[ind, ] <- plogis(eta_zi[ind, ], lower.tail = FALSE, log.p = TRUE) - 
            mu[ind, ] + y[ind] * eta[ind, ] - log(- expm1(-mu[ind, ])) - lgamma(y[ind] + 1)
        # zero part
        out[!ind, ] <- plogis(eta_zi[!ind, ], log.p = TRUE)
        attr(out, "mu_y") <- mu
        out
    }
    score_eta_fun <- function (y, mu, phis, eta_zi) {
        # binary indicator for y > 0
        ind <- y > 0
        # non-zero part
        mu <- as.matrix(mu)
        mu_ind <- mu[ind, ]
        out <- mu
        out[!ind, ] <- 0
        out[ind, ] <- - mu_ind + y[ind] + (exp(-mu_ind) * mu_ind) / expm1(-mu_ind) 
        out
    }
    score_eta_zi_fun <- function (y, mu, phis, eta_zi) {
        ind <- y > 0
        probs <- plogis(as.matrix(eta_zi))
        out <- 1 - probs
        out[ind, ] <- - probs[ind, ]
        out
    }
    simulate <- function (n, mu, phis, eta_zi) {
        y <- qpois(runif(n, ppois(0, mu), 1), mu)
        y[as.logical(rbinom(n, 1, plogis(eta_zi)))] <- 0
        y
    }
    structure(list(family = "hurdle poisson", link = stats$name, 
                   linkfun = stats$linkfun, linkinv = stats$linkinv, log_dens = log_dens,
                   score_eta_fun = score_eta_fun, score_eta_zi_fun = score_eta_zi_fun,
                   simulate = simulate),
              class = "family")
}

hurdle.negative.binomial <- function () {
    stats <- make.link("log")
    log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
        phis <- exp(phis)
        # binary indicator for y > 0
        ind <- y > 0
        # non-zero part
        eta <- as.matrix(eta)
        mu <- mu_fun(eta)
        log_mu_phis <- log(mu + phis)
        eta_zi <- as.matrix(eta_zi)
        out <- eta
        comp1 <- lgamma(y + phis) - lgamma(phis) - lgamma(y + 1)
        comp2 <- phis * log(phis) - phis * log_mu_phis
        comp3 <- y * log(mu) - y * log_mu_phis
        log_g <- comp1 + comp2 + comp3
        comp4 <- log(1 - (1 + mu / phis)^(-phis))
        out[ind, ] <- plogis(eta_zi[ind, ], lower.tail = FALSE, log.p = TRUE) + 
            log_g[ind, ] - comp4[ind, ]
        # zero part
        out[!ind, ] <- plogis(eta_zi[!ind, ], log.p = TRUE)
        attr(out, "mu_y") <- mu
        out
    }
    score_eta_fun <- function (y, mu, phis, eta_zi) {
        phis <- exp(phis)
        # binary indicator for y > 0
        ind <- y > 0
        # non-zero part
        mu <- as.matrix(mu)
        mu_phis <- mu + phis
        comp2 <- - phis / mu_phis
        comp3 <- y / mu - y / mu_phis
        k <- (1 + mu / phis)
        comp4 <- k^(- phis - 1) / (1 - k^(-phis))
        mu.eta <- mu
        out <- (comp2 + comp3 - comp4) * mu.eta
        out[!ind, ] <- 0
        out
    }
    score_eta_zi_fun <- function (y, mu, phis, eta_zi) {
        ind <- y > 0
        probs <- plogis(as.matrix(eta_zi))
        out <- 1 - probs
        out[ind, ] <- - probs[ind, ]
        out
    }
    score_phis_fun <- function (y, mu, phis, eta_zi) {
        ind_y0 <- y == 0
        phis <- exp(phis)
        mu <- as.matrix(mu)
        mu_phis <- mu + phis
        comp1 <- digamma(y + phis) - digamma(phis)
        comp2 <- log(phis) + 1 - log(mu_phis) - phis / mu_phis
        comp3 <- - y / mu_phis
        k <- mu / phis
        k1 <- 1 + k
        comp4 <- k1^(-phis) * (k / k1 - log(k1)) / (1 - k1^(-phis))
        out <- (comp1 + comp2 + comp3 + comp4) * phis
        out[ind_y0, ] <- 0
        out
    }
    simulate <- function (n, mu, phis, eta_zi) {
        y <- qnbinom(runif(n, pnbinom(0, mu = mu, size = exp(phis)), 1), 
                     mu = mu, size = exp(phis))
        y[as.logical(rbinom(n, 1, plogis(eta_zi)))] <- 0
        y
    }
    structure(list(family = "hurdle negative binomial", link = stats$name, 
                   linkfun = stats$linkfun, linkinv = stats$linkinv, log_dens = log_dens,
                   score_eta_fun = score_eta_fun, score_eta_zi_fun = score_eta_zi_fun,
                   score_phis_fun = score_phis_fun,
                   simulate = simulate),
              class = "family")
}

hurdle.lognormal <- function () {
    stats <- make.link("identity")
    log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
        sigma <- exp(phis)
        # binary indicator for y > 0
        ind <- y > 0
        # non-zero part
        eta <- as.matrix(eta)
        eta_zi <- as.matrix(eta_zi)
        out <- eta
        out[ind, ] <- plogis(eta_zi[ind, ], lower.tail = FALSE, log.p = TRUE) + 
            dnorm(x = log(y[ind]), mean = eta[ind, ], sd = sigma, log = TRUE)
        # zero part
        out[!ind, ] <- plogis(eta_zi[!ind, ], log.p = TRUE)
        attr(out, "mu_y") <- eta
        out
    }
    score_eta_fun <- function (y, mu, phis, eta_zi) {
        sigma <- exp(phis)
        # binary indicator for y > 0
        ind <- y > 0
        # non-zero part
        eta <- as.matrix(mu)
        out <- eta
        out[!ind, ] <- 0
        out[ind, ] <- (log(y[ind]) - eta[ind, ]) / sigma^2
        out
    }
    score_eta_zi_fun <- function (y, mu, phis, eta_zi) {
        ind <- y > 0
        probs <- plogis(as.matrix(eta_zi))
        out <- 1 - probs
        out[ind, ] <- - probs[ind, ]
        out
    }
    score_phis_fun <- function (y, mu, phis, eta_zi) {
        sigma <- exp(phis)
        # binary indicator for y > 0
        ind <- y > 0
        # non-zero part
        eta <- as.matrix(mu)
        out <- eta
        out[!ind, ] <- 0
        out[ind, ] <- - 1 + (log(y[ind]) - eta[ind, ])^2 / sigma^2
        out
    }
    simulate <- function (n, mu, phis, eta_zi) {
        y <- rnorm(n = n, mean = mu, sd = exp(phis))
        y[as.logical(rbinom(n, 1, plogis(eta_zi)))] <- 0
        y
    }
    structure(list(family = "hurdle log-normal", link = stats$name, 
                   linkfun = stats$linkfun, linkinv = stats$linkinv, log_dens = log_dens,
                   score_eta_fun = score_eta_fun, score_eta_zi_fun = score_eta_zi_fun,
                   score_phis_fun = score_phis_fun, simulate = simulate),
              class = "family")
}

beta.fam <- function () {
    stats <- make.link("logit")
    log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
        # the log density function
        phi <- exp(phis)
        mu <- mu_fun(eta)
        mu_phi <- mu * phi
        comp1 <- lgamma(phi) - lgamma(mu_phi)
        comp2 <- (mu_phi - 1) * log(y) - lgamma(phi - mu_phi)
        comp3 <- (phi - mu_phi - 1) * log(1 - y)
        out <- comp1 + comp2 + comp3
        attr(out, "mu_y") <- mu
        out
    }
    score_eta_fun <- function (y, mu, phis, eta_zi) {
        # the derivative of the log density w.r.t. mu
        phi <- exp(phis)
        mu_phi <- mu * phi
        comp1 <- - digamma(mu_phi) * phi
        comp2 <- phi * (log(y) + digamma(phi - mu_phi))
        comp3 <- - phi * log(1 - y)
        # the derivative of mu w.r.t. eta (this depends on the chosen link function)
        mu.eta <- mu - mu * mu
        (comp1 + comp2 + comp3) * mu.eta
    }
    score_phis_fun <- function (y, mu, phis, eta_zi) {
        phi <- exp(phis)
        mu_phi <- mu * phi
        mu1 <- 1 - mu
        comp1 <- digamma(phi) - digamma(mu_phi) * mu
        comp2 <- mu * log(y) - digamma(phi - mu_phi) * mu1
        comp3 <- log(1 - y) * mu1
        (comp1 + comp2 + comp3) * phi
    }
    structure(list(family = "beta", link = stats$name, linkfun = stats$linkfun,
                   linkinv = stats$linkinv, log_dens = log_dens, 
                   score_eta_fun = score_eta_fun, score_phis_fun = score_phis_fun),
              class = "family")
}
