logLik_mixed <- function (thetas, id, y, N, X, Z, offset, X_zi, Z_zi, offset_zi, GH, 
                          canonical, user_defined, Xty, Xty_weights, log_dens, mu_fun, var_fun, 
                          mu.eta_fun, score_eta_fun, score_eta_zi_fun, score_phis_fun, 
                          list_thetas, diag_D, penalized, pen_mu, pen_invSigma, pen_df,
                          weights) {
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
    log_wGH <- rep(log(wGH), each = length(unique(id)))
    #dets <- GH$dets
    log_dets <- GH$log_dets
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
    log_p_yb <- unname(rowsum(log_Lik, id, reorder = FALSE))
    log_p_b <- matrix(dmvnorm(b, rep(0, nRE), D, TRUE),
                      nrow(log_p_yb), ncol(log_p_yb), byrow = TRUE)
    #p_yb <- exp(log_p_yb + log_p_b)
    #if (any(zero_ind <- p_yb == 0.0)) {
    #    p_yb[zero_ind] <- 1e-300
    #}
    #p_y <- c(p_yb %*% wGH) * dets
    log_p_y <- rowLogSumExps(log_p_yb + log_p_b + log_wGH) + log_dets
    out <- - sum(if (is.null(weights)) log_p_y else weights * log_p_y, na.rm = TRUE)
    if (penalized)
        out <- out - dmvt(betas, mu = pen_mu, invSigma = pen_invSigma, df = pen_df)
    out
}

score_mixed <- function (thetas, id, y, N, X, Z, offset, X_zi, Z_zi, offset_zi, GH, 
                         canonical, user_defined, Xty, Xty_weights, log_dens, mu_fun, var_fun, 
                         mu.eta_fun, score_eta_fun, score_eta_zi_fun, score_phis_fun, 
                         list_thetas, diag_D, penalized, pen_mu, pen_invSigma, pen_df,
                         i_contributions = FALSE, weights) {
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
    log_wGH <- rep(log(wGH), each = length(unique(id)))
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
    log_p_yb <- unname(rowsum(log_Lik, id, reorder = FALSE))
    log_p_b <- matrix(dmvnorm(b, rep(0, nRE), D, TRUE),
                      nrow(log_p_yb), ncol(log_p_yb), byrow = TRUE)
    # log penalty include here dmvt(betas, pen_mean, invSigma = pen_invsds, df = pen_df)
    #p_yb <- exp(log_p_yb + log_p_b)
    #if (any(zero_ind <- p_yb == 0.0)) {
    #    p_yb[zero_ind] <- 1e-300
    #}
    #p_y <- c(p_yb %*% wGH)
    #p_by <- p_yb / p_y
    log_p_yb_b <- log_p_yb + log_p_b
    log_p_y <- rowLogSumExps(log_p_yb_b + log_wGH)
    p_by <- exp(log_p_yb_b - log_p_y)
    t_p_by <- t(p_by)
    n <- length(log_p_y)
    NN <- if (NCOL(y) == 2) nrow(y) else length(y)
    post_b <- apply(b, 2, function (b_k) colSums(t_p_by * matrix(b_k, ncol(Ztb), n) * wGH))
    post_b2 <- apply(b2, 2, function (b_k) colSums(t_p_by * matrix(b_k, ncol(Ztb), n) * wGH))
    if (!is.null(weights)) {
        post_b <- weights * post_b
        post_b2 <- weights * post_b2
    }
    post_vb <- post_b2 - if (nRE > 1) t(apply(post_b, 1, function (x) x %o% x)) else
        as.matrix(apply(post_b, 1, function (x) x %o% x))
    ###
    mu_y <- if (!is.null(attr(log_Lik, "mu_y"))) attr(log_Lik, "mu_y") else mu_fun(eta_y)
    score.betas <- if (user_defined) {
        ncx <- ncol(X)
        sc <- if (i_contributions) matrix(0.0, NN, ncx) else numeric(ncx)
        if (!is.null(score_eta_fun)) {
            z <- score_eta_fun(y, mu_y, phis, eta_zi)
            for (l in seq_len(ncx)) {
                cc <- if (is.null(weights)) {
                    drop(rowsum(X[, l] * z, id, reorder = FALSE))
                } else {
                    weights * drop(rowsum(X[, l] * z, id, reorder = FALSE))
                } 
                if (i_contributions) {
                    sc[, l] <- c((X[, l] * z * p_by[id, , drop = FALSE]) %*% wGH)
                } else {
                    sc[l] <- sum(c((cc * p_by) %*% wGH), na.rm = TRUE)
                }
            }
            - sc
        } else {
            l1 <- log_dens(y, eta_y + 1e-04, mu_fun, phis, eta_zi)
            l2 <- log_dens(y, eta_y - 1e-04, mu_fun, phis, eta_zi)
            z <- (l1 - l2) / (2 * 1e-04)
            for (l in seq_len(ncx)) {
                cc <- if (is.null(weights)) {
                    drop(rowsum(X[, l] * z, id, reorder = FALSE))
                } else {
                    weights * drop(rowsum(X[, l] * z, id, reorder = FALSE))
                }
                if (i_contributions) {
                    sc[, l] <- c((X[, l] * z * p_by[id, , drop = FALSE]) %*% wGH)
                } else {
                    sc[l] <- sum(c((cc * p_by) %*% wGH), na.rm = TRUE)
                }
            }
            - sc
        }
    } else {
        ncx <- ncol(X)
        sc <- if (i_contributions) matrix(0.0, NN, ncx) else numeric(ncx)
        if (canonical) {
            if (!is.null(N))
                mu_y <- mu_y * N
            for (l in seq_len(ncx)) {
                cc <- if (is.null(weights)) {
                    drop(rowsum(X[, l] * mu_y, id, reorder = FALSE))
                } else {
                    weights * drop(rowsum(X[, l] * mu_y, id, reorder = FALSE))
                }
                if (i_contributions) {
                    sc[, l] <- c((X[, l] * mu_y * p_by[id, , drop = FALSE]) %*% wGH)
                } else {
                    sc[l] <- sum(c((cc * p_by) %*% wGH), na.rm = TRUE)
                }
           }
            if (i_contributions) {
                - (X * if (NCOL(y) == 2) y[, 1] else y) + sc
            } else {
                if (is.null(weights)) - Xty + sc else - Xty_weights + sc
            }
        } else {
            var <- var_fun(mu_y)
            deriv <- mu.eta_fun(eta_y)
            z <- if (!is.null(N)) (y[, 1] - N * mu_y) * deriv / var else (y - mu_y) * deriv / var
            for (l in seq_len(ncx)) {
                cc <- if (is.null(weights)) {
                    drop(rowsum(X[, l] * z, id, reorder = FALSE))
                } else {
                    weights * drop(rowsum(X[, l] * z, id, reorder = FALSE))
                }
                if (i_contributions) {
                    sc[, l] <- c((X[, l] * z * p_by[id, , drop = FALSE]) %*% wGH)
                } else {
                    sc[l] <- sum(c((cc * p_by) %*% wGH), na.rm = TRUE)
                }
            }
            - sc
        }
    }
    if (penalized) {
        pen_invSigma_betas <- betas * diag(pen_invSigma) / pen_df
        fact <- (pen_df + ncx) / c(1 + crossprod(betas, pen_invSigma_betas))
        score.betas <- score.betas + pen_invSigma_betas * fact
    }
    ###
    score.phis <- if (!is.null(phis)) {
        if (is.null(score_phis_fun)) {
            n_phis <- length(phis)
            sc <- if (i_contributions) matrix(0.0, NN, n_phis) else numeric(n_phis)
            for (i in seq_len(n_phis)) {
                phis1 <- phis2 <- phis
                phis1[i] <- phis[i] + 1e-03
                phis2[i] <- phis[i] - 1e-03
                l1 <- log_dens(y, eta_y, mu_fun, phis1, eta_zi)
                l2 <- log_dens(y, eta_y, mu_fun, phis2, eta_zi)
                z <- (l1 - l2) / (phis1[i] - phis2[i])
                if (i_contributions) {
                    sc[, i] <- c((z * p_by[id, , drop = FALSE]) %*% wGH)
                } else {
                    cc <- if (is.null(weights)) {
                        c((rowsum(z, id, reorder = FALSE) * p_by) %*% wGH)
                    } else {
                        weights * c((rowsum(z, id, reorder = FALSE) * p_by) %*% wGH)
                    }
                    sc[i] <- sum(cc, na.rm = TRUE)
                }
            }
            - sc
        } else {
            z <- score_phis_fun(y, mu_y, phis, eta_zi)
            if (i_contributions) {
                -c((z * p_by[id, , drop = FALSE]) %*% wGH)
            } else {
                cc <- if (is.null(weights)) {
                    c((rowsum(z, id, reorder = FALSE) * p_by) %*% wGH)
                } else {
                    weights * c((rowsum(z, id, reorder = FALSE) * p_by) %*% wGH)
                }
                -sum(cc, na.rm = TRUE)
            }
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
        sc <- if (i_contributions) matrix(0.0, NN, ncx_zi) else numeric(ncx_zi)
        for (l in seq_len(ncx_zi)) {
            cc <- if (is.null(weights)) {
                drop(rowsum(X_zi[, l] * drop(z), id, reorder = FALSE))
            } else {
                weights * drop(rowsum(X_zi[, l] * drop(z), id, reorder = FALSE))
            }
            if (i_contributions) {
                sc[, l] <- c((X_zi[, l] * drop(z) * p_by[id, , drop = FALSE]) %*% wGH)
            } else {
                sc[l] <- sum(c((cc * p_by) %*% wGH), na.rm = TRUE)
            }
        }
        - sc
    }
    ###
    score.D <- if (diag_D) {
        D <- diag(D)
        svD <- 1/D
        svD2 <- svD^2
        if (i_contributions) {
            NA
        } else {
            cS.postVB <- colSums(as.matrix(post_vb), na.rm = TRUE)
            dim(cS.postVB) <- c(nRE, nRE)
            D * 0.5 * (n * svD - diag(cS.postVB) * svD2 - 
                           colSums(as.matrix(post_b^2), na.rm = TRUE) * svD2)
        }
    } else {
        svD <- solve(D)
        dD <- deriv_D(D)
        ndD <- length(dD)
        D1 <- sapply(dD, function (x) sum(svD * x))
        D2 <- t(sapply(dD, function (x) c(svD %*% x %*% svD)))
        if (i_contributions) {
            rr <- matrix(0.0, n, ndD)
            for (j in seq_len(n)) {
                cS.postVB <- colSums(as.matrix(post_vb)[j, , drop = FALSE], na.rm = TRUE)
                out <- numeric(ndD)
                for (i in seq_along(dD)) {
                    D.mat <- D2[i, ]
                    dim(D.mat) <- c(nRE, nRE)
                    out[i] <- sum(D2[i, ] * cS.postVB, na.rm = TRUE) +
                        sum((post_b[j, , drop = FALSE] %*% D.mat) * post_b[j, , drop = FALSE], na.rm = TRUE)
                }
                J <- jacobian2(attr(D, "L"), nRE)
                rr[j, ] <- drop(0.5 * (D1 - out) %*% J)
            }
            rr
        } else {
            cS.postVB <- colSums(as.matrix(post_vb), na.rm = TRUE)
            out <- numeric(ndD)
            for (i in seq_along(dD)) {
                D.mat <- D2[i, ]
                dim(D.mat) <- c(nRE, nRE)
                out[i] <- sum(D2[i, ] * cS.postVB, na.rm = TRUE) +
                    sum((post_b %*% D.mat) * post_b, na.rm = TRUE)
            }
            J <- jacobian2(attr(D, "L"), nRE)
            if (is.null(weights)) {
                drop(0.5 * (n * D1 - out) %*% J)
            } else {
                drop(0.5 * (sum(weights) * D1 - out) %*% J)
            }
        }
    }
    ###
    if (i_contributions)
        list(score.betas = score.betas, score.D = score.D, score.phis = score.phis, 
             score.gammas = score.gammas)
    else
        c(score.betas, score.D, score.phis, score.gammas)
}

score_betas <- function (betas, y, N, X, id, offset, weights, phis, Ztb, eta_zi, p_by, wGH, canonical,
                         user_defined, Xty, Xty_weights, log_dens, mu_fun, var_fun, mu.eta_fun,
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
                cc <- if (is.null(weights)) {
                    drop(rowsum(X[, l] * z, id, reorder = FALSE))
                } else {
                    weights * drop(rowsum(X[, l] * z, id, reorder = FALSE))
                }
                sc[l] <- sum(c((cc * p_by) %*% wGH), na.rm = TRUE)
            }
            - sc
        } else {
            l1 <- log_dens(y, eta_y + 1e-05, mu_fun, phis, eta_zi)
            l2 <- log_dens(y, eta_y - 1e-05, mu_fun, phis, eta_zi)
            z <- (l1 - l2) / (2 * 1e-05)
            for (l in seq_len(ncx)) {
                cc <- if (is.null(weights)) {
                    drop(rowsum(X[, l] * z, id, reorder = FALSE))
                } else {
                    weights * drop(rowsum(X[, l] * z, id, reorder = FALSE))
                }
                sc[l] <- sum(c((cc * p_by) %*% wGH), na.rm = TRUE)
            }
            - sc
        }
    } else {
        if (canonical) {
            if (!is.null(N))
                mu_y <- N * mu_y
            sc <- numeric(ncx)
            for (l in seq_len(ncx)) {
                cc <- if (is.null(weights)) {
                    rowsum(X[, l] * mu_y, id, reorder = FALSE)
                } else {
                    weights * rowsum(X[, l] * mu_y, id, reorder = FALSE)
                }
                sc[l] <- sum(c((cc * p_by) %*% wGH), na.rm = TRUE)
            }
            if (is.null(weights)) - Xty + sc else - Xty_weights + sc
        } else {
            var <- var_fun(mu_y)
            deriv <- mu.eta_fun(eta_y)
            z <- if (!is.null(N)) (y[, 1] - N * mu_y) * deriv / var else (y - mu_y) * deriv / var
            for (l in seq_len(ncx)) {
                cc <- if (is.null(weights)) {
                    rowsum(X[, l] * z, id, reorder = FALSE)
                } else {
                    weights * rowsum(X[, l] * z, id, reorder = FALSE)
                }
                sc[l] <- sum(c((cc * p_by) %*% wGH), na.rm = TRUE)
            }
            - sc
        }
    }
    if (penalized) {
        pen_invSigma_betas <- betas * diag(pen_invSigma) / pen_df
        fact <- (pen_df + ncx) / c(1 + crossprod(betas, pen_invSigma_betas))
        out <- out + pen_invSigma_betas * fact
    }
    out
}

score_phis <- function (phis, y, X, betas, Ztb, offset, weights, eta_zi, id, p_by,
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
            cc <- if (is.null(weights)) {
                c((rowsum(z, id, reorder = FALSE) * p_by) %*% wGH)
            } else {
                weights * c((rowsum(z, id, reorder = FALSE) * p_by) %*% wGH)
            }
            sc[i] <- sum(cc, na.rm = TRUE)
        }
        - sc
    } else {
        mu_y <- mu_fun(eta_y)
        z <- score_phis_fun(y, mu_y, phis, eta_zi)
        -sum(c((rowsum(z, id, reorder = FALSE) * p_by) %*% wGH), na.rm = TRUE)
    }
}

score_gammas <- function (gammas, y, X, betas, Ztb, offset, weights, X_zi, Z_zi, Z_zitb, offset_zi,
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
        cc <- if (is.null(weights)) {
            drop(rowsum(X_zi[, l] * z, id, reorder = FALSE))
        } else {
            weights * drop(rowsum(X_zi[, l] * z, id, reorder = FALSE))
        }
        sc[l] <- sum(c((cc * p_by) %*% wGH), na.rm = TRUE)
    }
    - sc
}

binomial_log_dens = function (y, eta, mu_fun, phis, eta_zi) {
    mu_y <- mu_fun(eta)
    out <- if (NCOL(y) == 2L) {
        dbinom(y[, 1L], y[, 1L] + y[, 2L], mu_y, TRUE)
    } else {
        dbinom(y, 1L, mu_y, TRUE)
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
    comp3 <- y * (log(mu) - log_mu_phis)
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
        comp3 <- y * (log(mu) - log_mu_phis)
        out <- comp1 + comp2 + comp3
        attr(out, "mu_y") <- mu
        out
    }
    score_eta_fun <- function (y, mu, phis, eta_zi) {
        # the derivative of the log density w.r.t. mu
        phis <- exp(phis)
        #mu_phis <- mu + phis
        #comp2 <- - phis / mu_phis
        #comp3 <- y / mu - y / mu_phis
        ## the derivative of mu w.r.t. eta (this depends on the chosen link function)
        #mu.eta <- mu
        #(comp2 + comp3) * mu.eta
        mu.mu_phis <- mu / (mu + phis)
        - phis * mu.mu_phis + y * (1 - mu.mu_phis) 
    }
    score_phis_fun <- function (y, mu, phis, eta_zi) {
        # the derivative of the log density w.r.t. phis
        phis <- exp(phis)
        mu_phis <- mu + phis
        #comp1 <- digamma(y + phis) - digamma(phis)
        #comp2 <- log(phis) + 1 - log(mu_phis) - phis / mu_phis
        #comp3 <- - y / mu_phis
        #(comp1 + comp2 + comp3) * phis
        y_phis <- y + phis
        comp1 <- log(phis) + 1 - digamma(phis)
        comp2 <- digamma(y_phis)
        comp3 <- - log(mu_phis) - y_phis / mu_phis
        (comp1 + comp2 + comp3) * phis
    }
    structure(list(family = "negative binomial", link = stats$name, 
                   linkfun = stats$linkfun, linkinv = stats$linkinv, log_dens = log_dens,
                   variance = function (mu, theta) mu + mu^2 / theta,
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
        comp3 <- y * (log(mu) - log_mu_phis)
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
        mu.mu_phis <- mu / (mu + phis)
        out <- - phis * mu.mu_phis + y * (1 - mu.mu_phis) 
        # ZI part
        ind_y0 <- y == 0
        lambda <- exp(as.matrix(eta_zi)[ind_y0, ])
        mu0 <- mu[ind_y0, ]
        t <- phis / (phis + mu0)
        den <- (lambda + t^phis) * (mu0 + phis)^2
        out[ind_y0, ] <- - phis^2 * t^(phis - 1) * mu0 / den
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
        lambda0 <- lambda[ind_y0, ]
        out[ind_y0, ] <- lambda0 / (lambda0 + t^phis) - lambda0 / (1 + lambda0)
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
                   variance = function (mu) (mu + mu^2)/(1 - exp(-mu)) - mu^2/((1 - exp(-mu))^2),
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
        y <- rlnorm(n = n, meanlog = mu, sdlog = exp(phis))
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
    simulate <- function (n, mu, phis, eta_zi) {
        phi <- exp(phis)
        rbeta(n, shape1 = mu * phi, shape2 = phi * (1 - mu))
    }
    structure(list(family = "beta", link = stats$name, linkfun = stats$linkfun,
                   linkinv = stats$linkinv, variance = function (mu) mu * (1 - mu), 
                   log_dens = log_dens, simulate = simulate,
                   score_eta_fun = score_eta_fun, score_phis_fun = score_phis_fun),
              class = "family")
}

hurdle.beta.fam <- function () {
    stats <- make.link("logit")
    log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
        phi <- exp(phis)
        # binary indicator for y > 0
        ind <- y > 0
        # non-zero part
        eta <- as.matrix(eta)
        eta_zi <- as.matrix(eta_zi)
        out <- eta
        mu <- mu_fun(eta)
        mu_phi <- mu * phi
        comp1 <- lgamma(phi) - lgamma(mu_phi)
        comp2 <- (mu_phi - 1) * log(y) - lgamma(phi - mu_phi)
        comp3 <- (phi - mu_phi - 1) * log(1 - y)
        out[ind, ] <- plogis(eta_zi[ind, ], lower.tail = FALSE, log.p = TRUE) + 
            comp1[ind, ] + comp2[ind, ] + comp3[ind, ]
        # zero part
        out[!ind, ] <- plogis(eta_zi[!ind, ], log.p = TRUE)
        attr(out, "mu_y") <- eta
        out
    }
    score_eta_fun <- function (y, mu, phis, eta_zi) {
        phi <- exp(phis)
        # binary indicator for y > 0
        ind <- y > 0
        # non-zero part
        mu <- as.matrix(mu)
        mu_phi <- mu * phi
        out <- mu
        comp1 <- - digamma(mu_phi) * phi
        comp2 <- phi * (log(y) + digamma(phi - mu_phi))
        comp3 <- - phi * log(1 - y)
        mu.eta <- mu - mu * mu
        out[!ind, ] <- 0
        out[ind, ] <- (comp1[ind, ] + comp2[ind, ] + comp3[ind]) * mu.eta[ind, ]
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
        phi <- exp(phis)
        # binary indicator for y > 0
        ind <- y > 0
        # non-zero part
        mu <- as.matrix(mu)
        mu_phi <- mu * phi
        mu1 <- 1 - mu
        out <- mu
        comp1 <- digamma(phi) - digamma(mu_phi) * mu
        comp2 <- mu * log(y) - digamma(phi - mu_phi) * mu1
        comp3 <- log(1 - y) * mu1
        out[ind, ] <- (comp1[ind, ] + comp2[ind, ] + comp3[ind, ]) * phi
        out[!ind, ] <- 0
        out
    }
    simulate <- function (n, mu, phis, eta_zi) {
        phi <- exp(phis)
        y <- rbeta(n, shape1 = mu * phi, shape2 = phi * (1 - mu))
        y[as.logical(rbinom(n, 1, plogis(eta_zi)))] <- 0
        y
    }
    structure(list(family = "hurdle beta", link = stats$name, 
                   linkfun = stats$linkfun, linkinv = stats$linkinv, log_dens = log_dens,
                   score_eta_fun = score_eta_fun, score_eta_zi_fun = score_eta_zi_fun,
                   score_phis_fun = score_phis_fun, simulate = simulate),
              class = "family")
}

students.t <- function (df = stop("'df' must be specified"), link = "identity") {
    .df <- df
    env <- new.env(parent = .GlobalEnv)
    assign(".df", df, envir = env)
    stats <- make.link(link)
    log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
        # the log density function
        sigma <- exp(phis)
        out <- dt(x = (y - eta) / sigma, df = .df, log = TRUE) - log(sigma)
        attr(out, "mu_y") <- eta
        out
    }
    score_eta_fun <- function (y, mu, phis, eta_zi) {
        # the derivative of the log density w.r.t. mu
        sigma2 <- exp(phis)^2
        y_mu <- y - mu
        (y_mu * (.df + 1) / (.df * sigma2)) / (1 + y_mu^2 / (.df * sigma2))
    }
    score_phis_fun <- function (y, mu, phis, eta_zi) {
        sigma <- exp(phis)
        y_mu2_df <- (y - mu)^2 / .df
        (.df + 1) * y_mu2_df * sigma^{-2} / (1 + y_mu2_df / sigma^2) - 1
    }
    simulate <- function (n, mu, phis, eta_zi) {
        phi <- exp(phis)
        mu + phi * rt(n, df = .df)
    }
    environment(log_dens) <- environment(score_eta_fun) <- env
    environment(score_phis_fun) <- environment(simulate) <- env
    structure(list(family = "Student's-t", link = stats$name, linkfun = stats$linkfun,
                   linkinv = stats$linkinv, log_dens = log_dens, 
                   variance = function (mu) rep.int(1, length(mu)),
                   score_eta_fun = score_eta_fun, score_phis_fun = score_phis_fun,
                   simulate = simulate),
              class = "family")
}

compoisson <- function (max = 100) {
    stats <- make.link("log")
    .max <- max
    env <- new.env(parent = .GlobalEnv)
    assign(".max", max, envir = env)
    log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
        # the log density function
        phis <- exp(phis)
        mu <- mu_fun(eta)
        Z <- function (lambda, nu, sumTo) {
            out <- lambda
            j <- seq(1, sumTo)
            log_lambda <- log(lambda)
            nu_log_factorial <- nu * cumsum(log(j))
            for (i in seq_along(out)) {
                out[i] <- 1 + sum(exp(j * log_lambda[i] - nu_log_factorial))
            }
            out
        }
        out <- y * log(mu) - phis * lgamma(y + 1) - log(Z(mu, phis, .max))
        attr(out, "mu_y") <- mu
        out
    }
    score_eta_fun <- function (y, mu, phis, eta_zi) {
        phi <- exp(phis)
        Y <- function (lambda, nu, sumTo) {
            out <- lambda
            j <- seq(1, sumTo)
            log_lambda <- log(lambda)
            log_j <- log(j)
            nu_log_factorial <- nu * cumsum(log_j)
            for (i in seq_along(out)) {
                F1 <- j * log_lambda[i] - nu_log_factorial
                num <- sum(exp(log_j + F1))
                den <- 1 + sum(exp(F1))
                out[i] <- num / den
            }
            out
        }
        y - Y(mu, phi, .max)
    }
    score_phis_fun <- function (y, mu, phis, eta_zi) {
        phi <- exp(phis)
        W <- function (lambda, nu, sumTo) {
            out <- lambda
            j <- seq(1, sumTo)
            log_lambda <- log(lambda)
            log_factorial <- cumsum(log(j))
            log_log_factorial <- log(log_factorial)
            nu_log_factorial <- nu * log_factorial
            for (i in seq_along(out)) {
                num <- sum(exp(j * log_lambda[i] + log_log_factorial - nu_log_factorial))
                den <- 1 + sum(exp(j * log_lambda[i] - nu_log_factorial))
                out[i] <- num / den
            }
            out
        }
        (- lgamma(y + 1) + W(mu, phi, .max)) * phi
    }
    simulate <- function (n, mu, phis, eta_zi) {
        phi <- exp(phis)
    }
    structure(list(family = "Conway Maxwell Poisson", link = stats$name, 
                   linkfun = stats$linkfun, linkinv = stats$linkinv, log_dens = log_dens,
                   score_eta_fun = score_eta_fun, score_phis_fun = score_phis_fun),
              class = "family")
}

unit.lindley <- function () {
    stop("currently the 'unit.lindley()' family is unavailable.")
    stats <- make.link("logit")
    log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
        # the log density function
        # you link logit(mu) to covariates
        # where mu = (1 / (1 + theta))
        mu <- as.matrix(mu_fun(eta))
        theta <- 1 / mu - 1
        comp1 <- 2 * log(theta) - log(1 + theta)
        comp2 <- - 3 * log(1 - y) 
        comp3 <- - (theta * y) / (1 - y)
        out <- comp1 + comp2 + comp3
        attr(out, "mu_y") <- mu
        out
    }
    score_eta_fun <- function (y, mu, phis, eta_zi) {
        mu <- as.matrix(mu)
        theta <- 1 / mu - 1
        # the derivative of the log density w.r.t. theta
        comp1 <- 2 / theta - 1 / (1 + theta)
        comp3 <- - y / (1 - y)
        # the derivative of theta w.r.t mu
        tht_mu <- - 1 / mu^2
        # the derivative of mu w.r.t. eta
        mu_eta <- mu - mu * mu
        (comp1 + comp3) * tht_mu * mu_eta
    }
    simulate <- function (n, mu, phis, eta_zi) {
        NA
    }
    structure(list(family = "unit Lindley", link = stats$name, 
                   linkfun = stats$linkfun, linkinv = stats$linkinv, 
                   log_dens = log_dens, score_eta_fun = score_eta_fun,
                   simulate = simulate,
                   variance = function (mu) mu * (1 - mu)),
              class = "family")
}

compoisson2 <- function (max = 100) {
    stats <- make.link("log")
    .max <- max
    env <- new.env(parent = .GlobalEnv)
    assign(".max", max, envir = env)
    log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
        # the log density function
        phis <- exp(phis)
        mu <- mu_fun(eta)
        Z <- function (lambda, nu, sumTo) {
            out <- lambda
            j <- seq(1, sumTo)
            log_lambda <- log(lambda)
            nu_log_factorial <- nu * cumsum(log(j))
            for (i in seq_along(out)) {
                out[i] <- 1 + sum(exp(j * log_lambda[i] - nu_log_factorial))
            }
            out
        }
        out <- y * log(mu) - phis * lgamma(y + 1) - log(Z(mu, phis, .max))
        attr(out, "mu_y") <- mu
        out
    }
    score_eta_fun <- function (y, mu, phis, eta_zi) {
        phi <- exp(phis)
        Y <- function (lambda, nu, sumTo) {
            out <- lambda
            j <- seq(1, sumTo)
            log_lambda <- log(lambda)
            log_j <- log(j)
            nu_log_factorial <- nu * cumsum(log_j)
            for (i in seq_along(out)) {
                F1 <- j * log_lambda[i] - nu_log_factorial
                num <- sum(exp(log_j + F1))
                den <- 1 + sum(exp(F1))
                out[i] <- num / den
            }
            out
        }
        y - Y(mu, phi, .max)
    }
    score_phis_fun <- function (y, mu, phis, eta_zi) {
        phi <- exp(phis)
        W <- function (lambda, nu, sumTo) {
            out <- lambda
            j <- seq(1, sumTo)
            log_lambda <- log(lambda)
            log_factorial <- cumsum(log(j))
            log_log_factorial <- log(log_factorial)
            nu_log_factorial <- nu * log_factorial
            for (i in seq_along(out)) {
                num <- sum(exp(j * log_lambda[i] + log_log_factorial - nu_log_factorial))
                den <- 1 + sum(exp(j * log_lambda[i] - nu_log_factorial))
                out[i] <- num / den
            }
            out
        }
        (- lgamma(y + 1) + W(mu, phi, .max)) * phi
    }
    simulate <- function (n, mu, phis, eta_zi) {
        phi <- exp(phis)
    }
    structure(list(family = "Conway Maxwell Poisson", link = stats$name, 
                   linkfun = stats$linkfun, linkinv = stats$linkinv, log_dens = log_dens,
                   score_eta_fun = score_eta_fun, score_phis_fun = score_phis_fun),
              class = "family")
}

find_lambda <- function (mu, nu, sumTo = 100) {
    j <- seq(1, sumTo)
    nu_log_factorial <- nu * cumsum(log(j))
    f <- function (lambda, mu) {
        fact <- exp(j * log(lambda) - nu_log_factorial)
        sum(c(-mu, (j - mu) * fact))
    }
    out <- mu
    init_lambda <- (mu + (nu - 1) / (2 * nu))^nu
    for (i in seq_along(mu)) {
        int <- c(max(1e-06, init_lambda[i] - 10), min(sumTo, init_lambda[i] + 10))
        test <- try(uniroot(f, interval = int, mu = mu[i])$root, silent = TRUE)
        if (inherits(test, "try-error")) {
            test <- try(uniroot(f, interval = c(1e-06, sumTo), mu = mu[i])$root, 
                        silent = TRUE)
        }
        if (inherits(test, "try-error")) {
            stop("it was not possible to find lambda parameter of the ", 
                 "Conway Maxwell Poisson distribution;\nre-fit the model using ",
                 "\n\n\tmixed_model(..., family = compoisson(max = XXX))\n\n",
                 "where 'XXX' is a big enough count.")
        }
        out[i] <- test
    }
    out
}


