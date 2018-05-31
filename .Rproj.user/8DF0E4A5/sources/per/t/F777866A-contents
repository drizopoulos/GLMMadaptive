print.MixMod <- function (x, digits = max(4, getOption("digits") - 4), ...) {
    cat("\nCall:\n", printCall(x$call), "\n\n", sep = "")
    cat("\nModel:")
    user_defined <- is.null(x$family)
    cat("\n family:", if (user_defined) "user-defined" else x$family$family)
    cat("\n link:", if (user_defined) "user-defined" else x$family$link, "\n")
    cat("\nRandom effects covariance matrix:\n")
    D <- x$D
    ncz <- nrow(D)
    diag.D <- all(abs(D[lower.tri(D)]) < sqrt(.Machine$double.eps))
    sds <- sqrt(diag(D))
    if (ncz > 1) {
        if (diag.D) {
            dat <- data.frame("StdDev" = round(sds, digits), row.names = rownames(D))
        } else {
            corrs <- cov2cor(D)
            corrs[upper.tri(corrs, TRUE)] <- 0
            mat <- round(cbind(sds, corrs[, -ncz]), digits)
            mat <- apply(mat, 2, sprintf, fmt = "% .4f")
            mat[mat == mat[1, 2]] <- ""
            mat[1, -1] <- abbreviate(colnames(mat)[-1], 6)
            colnames(mat) <- c(colnames(mat)[1], rep("", ncz - 1))
            dat <- data.frame(mat, check.rows = FALSE, check.names = FALSE)
            names(dat) <- c("StdDev", "Corr", if (ncz > 2) rep(" ", ncz - 2) else NULL)
            row.names(dat) <- dimnames(D)[[1]]
        }
    } else {
        dat <- data.frame("StdDev" = sds, row.names = rownames(D),
                          check.rows = FALSE, check.names = FALSE)
    }
    print(dat)
    cat("\nFixed effects:\n")
    print(x$coefficients)
    if (!is.null(x$phis)) {
        cat("\nphi parameters:\n", x$phis)
    }
    cat("\nlog-Lik:", x$logLik)
    cat("\n\n")
    invisible(x)

}

vcov.MixMod <- function (object, ...) {
    solve(object$Hessian)
}

logLik.MixMod <- function (object, ...) {
    out <- object$logLik
    attr(out, "df") <- nrow(object$Hessian)
    attr(out, "nobs") <- length(unique(object$id))
    class(out) <- "logLik"
    out
}

coef.MixMod <- function (object, ...) {
    betas <- fixef(object)
    b <- ranef(object)
    out <- matrix(betas, nrow = nrow(b), ncol = length(betas), byrow = TRUE)
    colnames(out) <- names(betas)
    rownames(out) <- rownames(b)
    out[, colnames(b)] <- out[, colnames(b)] + b
    out
}

fixef <- function (object, ...) UseMethod("fixef")

fixef.MixMod <- function(object, ...) {
    object$coefficients
}

ranef <- function (object, ...) UseMethod("ranef")

ranef.MixMod <- function(object, ...) {
    object$post_modes
}

summary.MixMod <- function (object, ...) {
    betas <- fixef(object)
    n_betas <- length(betas)
    V <- vcov(object)
    var_betas <- V[seq_len(n_betas), seq_len(n_betas)]
    ses <- sqrt(diag(var_betas))
    D <- object$D
    n_D <- length(D[lower.tri(D, TRUE)])
    coef_table <- cbind("Value" = betas, "Std.Err" = ses, "z-value" = betas / ses,
                        "p-value" = 2 * pnorm(abs(betas / ses), lower.tail = FALSE))
    out <- list(coef_table = coef_table, D = D, logLik = logLik(object),
                AIC = AIC(object), BIC = BIC(object), call = object$call,
                N = length(object$id))
    if (!is.null(object$phis)) {
        phis <- object$phis
        var_phis <- as.matrix(V[-seq_len(n_betas + n_D), -seq_len(n_betas + n_D)])
        out$phis_table <- cbind("Value" = phis, "Std.Err" = sqrt(diag(var_phis)))
    }
    out$control <- object$control
    out$family <- object$family
    out$converged <- object$converged
    class(out) <- 'summary.MixMod'
    out
}


print.summary.MixMod <- function (x, digits = max(4, getOption("digits") - 4), ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    cat("Data Descriptives:")
    cat("\nNumber of Observations:", x$N)
    cat("\nNumber of Groups:", attr(x$logLik, 'n'), "\n")
    cat("\nModel:")
    user_defined <- is.null(x$family)
    cat("\n family:", if (user_defined) "user-defined" else x$family$family)
    cat("\n link:", if (user_defined) "user-defined" else x$family$link, "\n")
    cat("\nFit statistics:\n")
    model.sum <- data.frame(log.Lik = x$logLik, AIC = x$AIC, BIC = x$BIC, row.names = "")
    print(model.sum)
    cat("\nRandom effects covariance matrix:\n")
    D <- x$D
    ncz <- nrow(D)
    diag.D <- all(abs(D[lower.tri(D)]) < sqrt(.Machine$double.eps))
    sds <- sqrt(diag(D))
    if (ncz > 1) {
        if (diag.D) {
            dat <- data.frame("StdDev" = round(sds, digits), row.names = rownames(D))
        } else {
            corrs <- cov2cor(D)
            corrs[upper.tri(corrs, TRUE)] <- 0
            mat <- round(cbind(sds, corrs[, -ncz]), digits)
            mat <- apply(mat, 2, sprintf, fmt = "% .4f")
            mat[mat == mat[1, 2]] <- ""
            mat[1, -1] <- abbreviate(colnames(mat)[-1], 6)
            colnames(mat) <- c(colnames(mat)[1], rep("", ncz - 1))
            dat <- data.frame(mat, check.rows = FALSE, check.names = FALSE)
            names(dat) <- c("StdDev", "Corr", if (ncz > 2) rep(" ", ncz - 2) else NULL)
            row.names(dat) <- dimnames(D)[[1]]
        }
    } else {
        dat <- data.frame("StdDev" = sds, row.names = rownames(D),
                          check.rows = FALSE, check.names = FALSE)
    }
    print(dat)
    cat("\nFixed effects:\n")
    coef_table <- as.data.frame(x$coef_table)
    coef_table[1:3] <- lapply(coef_table[1:3], round, digits = digits)
    coef_table[["p-value"]] <- format.pval(coef_table[["p-value"]], eps = 1e-04)
    print(coef_table)
    if (!is.null(x$phis_table)) {
        cat("\nphi parameters:\n")
        phis_table <- as.data.frame(x$phis_table)
        phis_table[] <- lapply(phis_table, round, digits = digits)
        print(phis_table)
    }
    cat("\nIntegration:")
    cat("\nmethod: adaptive Gauss-Hermite quadrature rule")
    cat("\nquadrature points:", x$control$nAGQ)
    cat("\n\nOptimization:")
    cat("\nmethod: hybrid EM and quasi-Newton")
    cat("\nconverged:", x$converged, "\n")
    invisible(x)
}

coef.summary.MixMod <- function (object, ...) {
    object$coef_table
}

confint.MixMod <- function (object, parm, level = 0.95, ...) {
    betas <- fixef(object)
    n_betas <- length(betas)
    V <- vcov(object)
    ses_betas <- sqrt(diag(V[seq_len(n_betas), seq_len(n_betas)]))
    out <- cbind(betas + qnorm((1 - level) / 2) * ses_betas,
                 betas + qnorm((1 + level) / 2) * ses_betas)
    colnames(out) <- paste(round(100 * c((1 - level) / 2, (1 + level) / 2), 1), "%")
    out
}

anova.MixMod <- function (object, object2, test = TRUE, L = NULL, ...) {
    if (missing(object2) && is.null(L))
        stop("either argument 'object2' or argument 'L' needs to be specified.\n")
    if (!missing(object2)) {
        if (!object$converged)
            warning("it seems that '", deparse(substitute(object)),
                    "' has not converged.\n")
        if (!object2$converged)
            warning("it seems that '", deparse(substitute(object2)),
                    "' has not converged.")
        L0 <- logLik(object)
        L1 <- logLik(object2)
        nb0 <- attr(L0, "df")
        nb1 <- attr(L1, "df")
        df <- nb1 - nb0
        if (test && df < 0)
            stop("'object' should be nested in 'object2'.\n")
        out <- list(nam0 = deparse(substitute(object)), L0 = L0,
                    aic0 = AIC(object), bic0 = BIC(object),
                    nam1 = deparse(substitute(object2)), L1 = L1, aic1 = AIC(object2),
                    bic1 = AIC(object2), df = df, test = test)
        if (test) {
            LRT <- - 2 * (L0 - L1)
            attributes(LRT) <- NULL
            if (LRT < 0)
                warning("either the two models are not nested or the model ",
                        "represented by 'object2' fell on a local maxima.\n")
            out$LRT <- LRT
            out$p.value <- pchisq(LRT, df, lower.tail = FALSE)
        }
    } else {
        betas <- fixef(object)
        n_betas <- length(betas)
        if (!is.numeric(L) || ncol(L) != n_betas)
            stop("L matrix not of appropriate dimensions. ",
                 "It should have ", n_betas, " columns.\n")
        colnames(L) <- abbreviate(names(betas), 6)
        rownames(L) <- rep("", nrow(L))
        V <- vcov(object)
        var_betas <- V[seq_len(n_betas), seq_len(n_betas)]
        Lbetas <- c(L %*% betas)
        LVtL <- L %*% tcrossprod(var_betas, L)
        stat <- c(crossprod(Lbetas, solve(LVtL, Lbetas)))
        pval <- pchisq(stat, nrow(L), lower.tail = FALSE)
        res <- data.frame(Chisq = stat, df = nrow(L),
                          "Pr(>|Chi|)" = pval, check.names = FALSE, row.names = " ")
        out <- list(aovTab.L = res, L = L)
    }
    class(out) <- "aov.MixMod"
    out
}

print.aov.MixMod <- function (x, ...) {
    if (is.null(x$L0)) {
        f <- function (dat) {
            dat[] <- lapply(dat, function (x)
                round(unlist(x), 4))
            dat$'Pr(>|Chi|)' <- format.pval(dat$'Pr(>|Chi|)', eps = 1e-04)
            dat
        }
        cat("\nMarginal Wald Tests Table\n")
        if (!is.null(x$aovTab.L)) {
            cat("\nUser-defined contrasts matrix:\n")
            print(x$L)
            cat("\n")
            print(f(x$aovTab.L))
        }
        cat("\n")
    } else {
        dat <- if (x$test) {
            p.val <- round(x$p.value, 4)
            p.val <- if (p.val < 0.0001) "<0.0001" else p.val
            data.frame(AIC = round(c(x$aic0, x$aic1), 2),
                       BIC = round(c(x$bic0, x$bic1), 2),
                       log.Lik = round(c(x$L0, x$L1), 2),
                       LRT = c(" ", round(x$LRT, 2)), df = c("", x$df),
                       p.value = c("", p.val), row.names = c(x$nam0, x$nam1))
        } else {
            data.frame(AIC = round(c(x$aic0, x$aic1), 2),
                       BIC = round(c(x$bic0, x$bic1), 2),
                       log.Lik = round(c(x$L0, x$L1), 2), df = c("", x$df),
                       row.names = c(x$nam0, x$nam1))
        }
        cat("\n")
        print(dat)
        cat("\n")
    }
    invisible(x)
}

fitted.MixMod <- function (object, type = c("mean_subject", "subject_specific", "marginal"),
                           link_fun = NULL, ...) {
    type <- match.arg(type)
    X <- model.matrix(object$Terms$termsX, object$model_frames$mfX)
    eta <- if (type == "mean_subject") {
        betas <- fixef(object)
        eta <- c(X %*% betas)
    } else if (type == "subject_specific") {
        betas <- fixef(object)
        b <- ranef(object)
        Z <- model.matrix(object$Terms$termsZ, object$model_frames$mfZ)
        id <- match(object$id, unique(object$id))
        eta <- c(X %*% betas) + rowSums(Z * b[id, , drop = FALSE])
    } else {
        betas <- marginal_coefs(object, link_fun = link_fun)$betas
        eta <- c(X %*% betas)
    }
    out <- object$Funs$mu_fun(eta)
    names(out) <- rownames(X)
    out
}

residuals.MixMod <- function (object, type = c("mean_subject", "subject_specific",
                                               "marginal"), link_fun = NULL, ...) {
    type <- match.arg(type)
    fits <- fitted(object, type = type, link_fun = link_fun)
    y <- model.response(object$model_frames$mfX)
    y - fits
}

marginal_coefs <- function (object, ...) UseMethod("marginal_coefs")

marginal_coefs.MixMod <- function (object, std_errors = FALSE, link_fun = NULL, 
                                   M = 3000, K = 100,
                                   seed = 1, cores = max(parallel::detectCores() - 1, 1), 
                                   ...) {
    X <- model.matrix(object$Terms$termsX, object$model_frames$mfX)
    Z <- model.matrix(object$Terms$termsZ, object$model_frames$mfZ)
    betas <- fixef(object)
    D <- object$D
    compute_marg_coefs <- function (object, X, betas, Z, D, M, link_fun, seed) {
        mu_fun <- object$Funs$mu_fun
        if (is.null(link_fun)) {
            link_fun <- object$family$linkfun
        }
        if (is.null(link_fun)) {
            stop("you must specify the 'link_fun' argument.\n")
        }
        Xbetas <- c(X %*% betas)
        nRE <- ncol(Z)
        n <- nrow(X)
        eS <- eigen(D, symmetric = TRUE)
        ev <- eS$values
        V <- eS$vectors %*% diag(sqrt(pmax(ev, 0)), nRE)
        marg_inv_mu <- numeric(n)
        for (i in seq_len(n)) {
            set.seed(seed + i)
            b <- t(V %*% t(matrix(rnorm(M * nRE), M, nRE)))
            Zb <- c(Z[i, , drop = FALSE] %*% t(b))
            mu <- mu_fun(Xbetas[i] + Zb)
            marg_inv_mu[i] <- link_fun(mean(mu))
        }
        rm(list = ".Random.seed", envir = globalenv())
        c(solve(crossprod(X), crossprod(X, marg_inv_mu)))
        res <- c(solve(crossprod(X), crossprod(X, marg_inv_mu)))
        names(res) <- names(betas)
        res
    }
    out <- list(betas = compute_marg_coefs(object, X, betas, Z, D, M, link_fun, seed))
    if (std_errors) {
        blocks <- split(seq_len(K), rep(seq_len(cores), each = ceiling(K / cores),
                                        length.out = K))

        list_thetas <- list(betas = betas, D = chol_transf(D))
        tht <- unlist(as.relistable(list_thetas))
        V <- vcov(object)
        cluster_compute_marg_coefs <- function (block, tht, list_thetas, V, XX, Z, M,
                                                compute_marg_coefs, chol_transf,
                                                object, link_fun, seed) {
            n_block <- length(block)
            m_betas <- matrix(0.0, n_block, length(list_thetas[["betas"]]))
            for (b in seq_along(block)) {
                new_tht <- relist(MASS::mvrnorm(1, tht, V), skeleton = list_thetas)
                new_betas <- new_tht$betas
                new_D <- chol_transf(new_tht$D)
                m_betas[b, ] <- compute_marg_coefs(object, XX, new_betas, Z, new_D, M,
                                                   link_fun, seed = seed + block[b])
            }
            m_betas
        }

        cl <- makeCluster(cores)
        res <- parLapply(cl, blocks, cluster_compute_marg_coefs, tht = tht,
                         list_thetas = list_thetas, V = V, XX = X, Z = Z, M = M,
                         object = object, compute_marg_coefs = compute_marg_coefs,
                         chol_transf = chol_transf, link_fun = link_fun, seed = seed)
        stopCluster(cl)
        out$var_betas <- var(do.call("rbind", res))
        dimnames(out$var_betas) <- list(names(out$betas), names(out$betas))
        ses <- sqrt(diag(out$var_betas))
        coef_table <- cbind("Value" = out$betas, "Std.Err" = ses,
                            "z-value" = out$betas / ses,
                            "p-value" = 2 * pnorm(abs(out$betas / ses), lower.tail = FALSE))
        out$coef_table <- coef_table
    }
    class(out) <- "m_coefs"
    out
}

print.m_coefs <- function (x, digits = max(4, getOption("digits") - 4), ...) {
    if (is.null(x$coef_table)) {
        print(round(x$betas, digits = digits))
    } else {
        coef_table <- as.data.frame(x$coef_table)
        coef_table[1:3] <- lapply(coef_table[1:3], round, digits = digits)
        coef_table[["p-value"]] <- format.pval(coef_table[["p-value"]], eps = 1e-04)
        print(coef_table)
        cat("\n")
    }
    invisible(x)
}

effectPlotData <- function (object, newdata, level, ...) UseMethod("effectPlotData")

effectPlotData.MixMod <- function (object, newdata, level = 0.95, marginal = FALSE, ...) {
    if (marginal) {
        mcoefs <- marginal_coefs(object, std_errors = TRUE, ...)
        betas <- mcoefs$betas
        var_betas <- mcoefs$var_betas
    } else {
        betas <- fixef(object)
        n_betas <- length(betas)
        V <- vcov(object)
        var_betas <- V[seq_len(n_betas), seq_len(n_betas)]
    }
    termsX <- delete.response(object$Terms$termsX)
    mfX <- model.frame(termsX, newdata, 
                       xlev = .getXlevels(termsX, object$model_frames$mfX))
    X <- model.matrix(termsX, mfX)
    pred <- c(X %*% betas)
    ses <- sqrt(diag(X %*% var_betas %*% t(X)))
    newdata$pred <- pred
    newdata$low <- pred + qnorm((1 - level) / 2) * ses
    newdata$upp <- pred + qnorm((1 + level) / 2) * ses
    newdata
}

predict.MixMod <- function (object, newdata, type = c("link", "response"),
                            level = c("mean_subject", "subject_specific", "marginal"),
                            se.fit = FALSE, ...) {
    type <- match.arg(type)
    level <- match.arg(level)
    termsX <- delete.response(object$Terms$termsX)
    mfX <- model.frame(termsX, newdata, 
                       xlev = .getXlevels(termsX, object$model_frames$mfX))
    X <- model.matrix(termsX, mfX)
    if (level %in% c("mean_subject", "marginal")) {
        if (level == "mean_subject") {
            betas <- fixef(object)
            n_betas <- length(betas)
            V <- vcov(object)
            var_betas <- V[seq_len(n_betas), seq_len(n_betas)]
            pred <- if (type == "link") c(X %*% betas) else object$family$linkinv(c(X %*% betas))
            se.fit <- if (se.fit) sqrt(diag(X %*% var_betas %*% t(X)))
        } else {
            mcoefs <- marginal_coefs(object, std_errors = TRUE, ...)
            betas <- mcoefs$betas
            var_betas <- mcoefs$var_betas
            pred <- if (type == "link") c(X %*% betas) else object$family$linkinv(c(X %*% betas))
            se.fit <- if (se.fit) sqrt(diag(X %*% var_betas %*% t(X)))
        }
    } else {
        y <- model.response(model.frame(object$Terms$termsX, newdata))
        if (is.factor(y)) {
            if (family$family == "binomial")
                y <- as.numeric(y != levels(y)[1L])
            else
                stop("the response variable should not be a factor.\n")
        }
        offset <- model.offset(mfX)
        termsZ <- delete.response(object$Terms$termsZ)
        mfZ <- model.frame(termsZ, newdata, 
                           xlev = .getXlevels(termsZ, object$model_frames$mfZ))
        Z <- model.matrix(termsZ, mfZ)
        id_nam <- object$id_name
        id <- newdata[[id_nam]]
        id <- match(id, unique(id))
        id_unq <- unique(id)
        y_lis <- if (is.matrix(y)) lapply(id_unq, function (i) y[id == i, ]) else split(y, id)
        X_lis <- lapply(id_unq, function (i) X[id == i, , drop = FALSE])
        Z_lis <- lapply(id_unq, function (i) Z[id == i, , drop = FALSE])
        offset_lis <- if (!is.null(offset)) split(offset, id)
        Zty_lis <- lapply(mapply(crossprod, Z_lis, y_lis, SIMPLIFY = FALSE), drop)
        Xty <- drop(crossprod(X, y))
        log_dens <- object$Funs$log_dens
        mu_fun <- object$Funs$mu_fun
        var_fun <- object$Funs$var_fun
        mu.eta_fun <- object$Funs$mu.eta_fun
        score_eta_fun <- object$Funs$score_eta_fun
        score_phis_fun <- object$Funs$score_phis_fun
        family <- object$family
        canonical <- !is.null(family$family) &&
            ((family$family == "binomial" && family$link == "logit") ||
                 (family$family == "poisson" && family$link == "log"))
        known_families <- c("binomial", "poisson", "negative binomial")
        user_defined <- !family$family %in% known_families
        numer_deriv <- if (object$control$numeric_deriv == "fd") fd else cd
        numer_deriv_vec <- if (object$control$numeric_deriv == "fd") fd_vec else cd_vec
        start <- matrix(0.0, length(y_lis), ncol(Z))
        betas <- fixef(object)
        invD <- solve(object$D)
        phis <- object$phis
        post_modes <- find_modes(start, y_lis, X_lis, Z_lis, offset_lis, 
                                 betas, invD, phis, canonical, user_defined, Zty_lis, 
                                 log_dens, mu_fun, var_fun, mu.eta_fun,
                                 score_eta_fun, score_phis_fun)$post_modes
        eta <- c(X %*% betas) + rowSums(Z * post_modes[id, , drop = FALSE])
        pred <- if (type == "link") eta else object$family$linkinv(eta)
    }
    if (se.fit) {
        list(pred = pred, se.fit = se.fit)
    } else {
        pred
    }
}



