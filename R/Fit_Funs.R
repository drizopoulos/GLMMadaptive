logLik_mixed <- function (thetas, id, y, N, X, Z, offset, phis, Ztb, GH, canonical,
                          user_defined, Xty, log_dens, mu_fun, var_fun, mu.eta_fun,
                          score_eta_fun, score_phis_fun, list_thetas, diag_D) {
    thetas <- relist(thetas, skeleton = list_thetas)
    betas <- thetas$betas
    phis <- thetas$phis
    D <- if (diag_D) diag(exp(thetas$D), length(thetas$D)) else chol_transf(thetas$D)
    ##
    b <- GH$b
    Ztb <- GH$Ztb
    wGH <- GH$wGH
    dets <- GH$dets
    ##
    eta_y <- as.vector(X %*% betas) + Ztb
    if (!is.null(offset))
        eta_y <- eta_y + offset
    log_p_yb <- rowsum(log_dens(y, eta_y, mu_fun, phis), id, reorder = FALSE)
    log_p_b <- matrix(dmvnorm(b, rep(0, ncol(Z)), D, TRUE),
                      nrow(log_p_yb), ncol(log_p_yb), byrow = TRUE)
    p_yb <- exp(log_p_yb + log_p_b)
    if (any(zero_ind <- p_yb == 0.0)) {
        p_yb[zero_ind] <- 1e-30
    }
    p_y <- c(p_yb %*% wGH) * dets
    - sum(log(p_y), na.rm = TRUE)
}

score_mixed <- function (thetas, id, y, N, X, Z, offset, phis, Ztb, GH, canonical,
                         user_defined, Xty, log_dens, mu_fun, var_fun, mu.eta_fun,
                         score_eta_fun, score_phis_fun, list_thetas, diag_D) {
    thetas <- relist(thetas, skeleton = list_thetas)
    betas <- thetas$betas
    phis <- thetas$phis
    D <- if (diag_D) diag(exp(thetas$D), length(thetas$D)) else chol_transf(thetas$D)
    ##
    b <- GH$b
    b2 <- GH$b2
    Ztb <- GH$Ztb
    wGH <- GH$wGH
    ##
    eta_y <- as.vector(X %*% betas) + Ztb
    ncz <- ncol(Z)
    if (!is.null(offset))
        eta_y <- eta_y + offset
    log_Lik <- log_dens(y, eta_y, mu_fun, phis)
    log_p_yb <- rowsum(log_Lik, id, reorder = FALSE)
    log_p_b <- matrix(dmvnorm(b, rep(0, ncz), D, TRUE), nrow(log_p_yb),
                      ncol(log_p_yb), byrow = TRUE)
    p_yb <- exp(log_p_yb + log_p_b)
    if (any(zero_ind <- p_yb == 0.0)) {
        p_yb[zero_ind] <- 1e-30
    }
    p_y <- c(p_yb %*% wGH)
    p_by <- p_yb / p_y
    t_p_by <- t(p_by)
    n <- length(p_y)
    post_b <- apply(b, 2, function (b_k) colSums(t_p_by * matrix(b_k, ncol(Ztb), n) * wGH))
    post_b2 <- apply(b2, 2, function (b_k) colSums(t_p_by * matrix(b_k, ncol(Ztb), n) * wGH))
    post_vb <- post_b2 - if (ncz > 1) t(apply(post_b, 1, function (x) x %o% x)) else
        as.matrix(apply(post_b, 1, function (x) x %o% x))
    ###
    mu_y <- if (!is.null(attr(log_Lik, "mu_y"))) attr(log_Lik, "mu_y") else mu_fun(eta_y)
    score.betas <- if (user_defined) {
        if (!is.null(score_eta_fun)) {
            z <- score_eta_fun(y, mu_y, phis)
            ncx <- ncol(X)
            sc <- numeric(ncx)
            for (l in seq_len(ncx)) {
                cc <- rowsum(X[, l] * z, id, reorder = FALSE)
                sc[l] <- sum(c((cc * p_by) %*% wGH))
            }
            - sc
        } else {
            l1 <- log_dens(y, eta_y + 1e-04, mu_fun, phis)
            l2 <- log_dens(y, eta_y - 1e-04, mu_fun, phis)
            z <- (l1 - l2) / (2 * 1e-04)
            ncx <- ncol(X)
            sc <- numeric(ncx)
            for (l in seq_len(ncx)) {
                cc <- rowsum(X[, l] * z, id, reorder = FALSE)
                sc[l] <- sum(c((cc * p_by) %*% wGH))
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
                cc <- rowsum(X[, l] * mu_y, id, reorder = FALSE)
                sc[l] <- sum(c((cc * p_by) %*% wGH))
            }
            - Xty + sc
        } else {
            var <- var_fun(mu_y)
            deriv <- mu.eta_fun(eta_y)
            z <- if (!is.null(N)) (y[, 1] - N * mu_y) * deriv / var else (y - mu_y) * deriv / var
            ncx <- ncol(X)
            sc <- numeric(ncx)
            for (l in seq_len(ncx)) {
                cc <- rowsum(X[, l] * z, id, reorder = FALSE)
                sc[l] <- sum(c((cc * p_by) %*% wGH))
            }
            - sc
        }
    }
    ###
    score.phis <- if (all(!is.na(phis))) {
        if (is.null(score_phis_fun)) {
            n_phis <- length(phis)
            sc <- numeric(n_phis)
            for (i in seq_len(n_phis)) {
                phis1 <- phis2 <- phis
                phis1[i] <- phis[i] + 1e-03
                phis2[i] <- phis[i] - 1e-03
                l1 <- log_dens(y, eta_y, mu_fun, phis1)
                l2 <- log_dens(y, eta_y, mu_fun, phis2)
                z <- (l1 - l2) / (phis1[i] - phis2[i])
                sc[i] <- sum(c((rowsum(z, id, reorder = FALSE) * p_by) %*% wGH), na.rm = TRUE)
            }
            - sc
        } else {
            z <- score_phis_fun(y, mu_y, phis)
            -sum(c((rowsum(z, id, reorder = FALSE) * p_by) %*% wGH), na.rm = TRUE)
        }
    }
    ###
    score.D <- if (diag_D) {
        D <- diag(D)
        svD <- 1/D
        svD2 <- svD^2
        cS.postVB <- colSums(as.matrix(post_vb), na.rm = TRUE)
        dim(cS.postVB) <- c(ncz, ncz)
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
            dim(D.mat) <- c(ncz, ncz)
            out[i] <- sum(D2[i, ] * cS.postVB, na.rm = TRUE) +
                sum((post_b %*% D.mat) * post_b, na.rm = TRUE)
        }
        J <- jacobian2(attr(D, "L"), ncz)
        drop(0.5 * (n * D1 - out) %*% J)
    }
    ###
    c(score.betas, score.D, score.phis)
}

score_betas <- function (betas, y, N, X, id, offset, phis, Ztb, p_by, wGH, canonical,
                         user_defined, Xty, log_dens, mu_fun, var_fun, mu.eta_fun,
                         score_eta_fun, score_phis_fun) {
    eta_y <- as.vector(X %*% betas) + Ztb
    if (!is.null(offset))
        eta_y <- eta_y + offset
    mu_y <- mu_fun(eta_y)
    ncx <- ncol(X)
    sc <- numeric(ncx)
    if (user_defined) {
        if (!is.null(score_eta_fun)) {
            z <- score_eta_fun(y, mu_y, phis)
            for (l in seq_len(ncx)) {
                cc <- rowsum(X[, l] * z, id, reorder = FALSE)
                sc[l] <- sum(c((cc * p_by) %*% wGH))
            }
            - sc
        } else {
            l1 <- log_dens(y, eta_y + 1e-05, mu_fun, phis)
            l2 <- log_dens(y, eta_y - 1e-05, mu_fun, phis)
            z <- (l1 - l2) / (2 * 1e-05)
            for (l in seq_len(ncx)) {
                cc <- rowsum(X[, l] * z, id, reorder = FALSE)
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
}

score_phis <- function (phis, y, X, betas, Ztb, offset, id, p_by,
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
            l1 <- log_dens(y, eta_y, mu_fun, phis1)
            l2 <- log_dens(y, eta_y, mu_fun, phis2)
            z <- (l1 - l2) / (phis1[i] - phis2[i])
            sc[i] <- sum(c((rowsum(z, id, reorder = FALSE) * p_by) %*% wGH), na.rm = TRUE)
        }
        - sc
    } else {
        mu_y <- mu_fun(eta_y)
        z <- score_phis_fun(y, mu_y, phis)
        -sum(c((rowsum(z, id, reorder = FALSE) * p_by) %*% wGH), na.rm = TRUE)
    }
}

binomial_log_dens = function (y, eta, mu_fun, phis) {
    mu_y <- mu_fun(eta)
    out <- if (NCOL(y) == 2) {
        dbinom(y[, 1], y[, 1] + y[, 2], mu_y, TRUE)
    } else {
        dbinom(y, 1, mu_y, TRUE)
    }
    attr(out, "mu_y") <- mu_y
    out
}

poisson_log_dens = function (y, eta, mu_fun, phis) {
    mu_y <- mu_fun(eta)
    out <- dpois(y, mu_y, TRUE)
    attr(out, "mu_y") <- mu_y
    out
}

negative.binomial_log_dens <- function (y, eta, mu_fun, phis) {
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

negative.binomial <- function (theta = stop("'theta' must be specified"), link = "log") {
    linktemp <- substitute(link)
    if (!is.character(linktemp))
        linktemp <- deparse(linktemp)
    if (linktemp %in% c("log", "identity", "sqrt"))
        stats <- make.link(linktemp)
    else if (is.character(link)) {
        stats <- make.link(link)
        linktemp <- link
    } else {
        if (inherits(link, "link-glm")) {
            stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        } else {
            stop(gettextf("\"%s\" link not available for negative binomial family; available links are \"identity\", \"log\" and \"sqrt\"",
                           linktemp))
        }
    }
    .Theta <- theta
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", theta, envir = env)
    variance <- function (mu) mu + mu^2/.Theta
    validmu <- function (mu) all(mu > 0)
    dev.resids <- function (y, mu, wt) {
        2 * wt * (y * log(pmax(1, y)/mu) - (y + .Theta) * log((y + .Theta)/ (mu + .Theta)))
    }
    aic <- function(y, n, mu, wt, dev) {
        term <- (y + .Theta) * log(mu + .Theta) - y * log(mu) +
            lgamma(y + 1) - .Theta * log(.Theta) + lgamma(.Theta) - lgamma(.Theta+y)
        2 * sum(term * wt)
    }
    initialize <- expression({
        if (any(y < 0))
            stop("negative values not allowed for the negative binomial family")
        n <- rep(1, nobs)
        mustart <- y + (y == 0)/6
    })
    environment(variance) <- environment(validmu) <- environment(dev.resids) <- environment(aic) <- env
    score_phis_fun <- function (y, mu, phis) {
        phis <- exp(phis)
        mu_phis <- mu + phis
        comp1 <- digamma(y + phis) - digamma(phis)
        comp2 <- log(phis) + 1 - log(mu_phis) - phis / (mu_phis)
        comp3 <- - y / (mu_phis)
        comp1 + comp2 + comp3
    }
    famname <- "negative binomial"
    structure(list(family = famname, link = linktemp, linkfun = stats$linkfun,
                   linkinv = stats$linkinv, variance = variance, mu.eta = stats$mu.eta,
                   dev.resids = dev.resids, initialize = initialize, validmu = validmu,
                   aic = aic, valideta = stats$valideta, score_phis_fun = score_phis_fun),
              class = "family")
}

