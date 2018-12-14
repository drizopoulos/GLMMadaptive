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
    if (!is.null(x$gammas)) {
        cat("\nZero-part coefficients:\n")
        print(x$gammas)
    }
    if (!is.null(x$phis)) {
        if (x$family$family %in% c("negative binomial", "zero-inflated negative binomial")) {
            cat("\ndispersion parameter:\n", exp(x$phis), "\n")
        } else if (x$family$family %in% c("hurdle log-normal")) {
            cat("\nResidual std. dev.:\n", exp(x$phis), "\n")
        } else {
            cat("\nphi parameters:\n", x$phis, "\n")
        }
    }
    cat("\nlog-Lik:", x$logLik)
    cat("\n\n")
    invisible(x)
}

vcov.MixMod <- function (object, parm = c("all", "fixed-effects", "var-cov","extra", 
                                          "zero_part"), sandwich = FALSE, ...) {
    parm <- match.arg(parm)
    V <- solve(object$Hessian)
    if (sandwich) {
        meat <- object$score_vect_contributions
        ind <- !names(meat) %in% "score.D" & !sapply(meat, is.null)
        meat[ind] <- lapply(meat[ind], rowsum, group = object$id, reorder = FALSE)
        meat <- do.call('cbind', meat)
        meat <- Reduce("+", lapply(split(meat, row(meat)), function (x) x %o% x))
        V <- V %*% meat %*% V
    }
    if (parm == "all") {
        return(V)
    }
    if (parm == "fixed-effects") {
        n_betas <- length(object$coefficients)
        return(V[seq_len(n_betas), seq_len(n_betas), drop = FALSE])
    }
    if (parm == "var-cov") {
        D <- object$D
        diag_D <- ncol(D) > 1 && all(abs(D[lower.tri(D)]) < sqrt(.Machine$double.eps))
        include <- if (diag_D) {
            unconstr_D <- log(diag(D))
            n_betas <- length(object$coefficients)
            seq(n_betas + 1, n_betas + length(unconstr_D))
        } else {
            unconstr_D <- chol_transf(D)
            n_betas <- length(object$coefficients)
            seq(n_betas + 1, n_betas + length(unconstr_D))
        }
        return(V[include, include, drop = FALSE])
    }
    if (parm == "extra") {
        if (is.null(object$phis)) {
            stop("the model behind 'object' contains no extra (phis) parameters.\n")
        } else {
            ind_phis <- grep("phi_", colnames(V), fixed = TRUE)
            return(V[ind_phis, ind_phis, drop = FALSE])
        }
    }
    if (parm == "zero_part") {
        if (is.null(object$gammas)) {
            stop("the fitted model does not have an extra zero part.")
        } else {
            gammas <- object$gammas
            ind_gammas <- grep("zi_", colnames(V), fixed = TRUE)
            return(V[ind_gammas, ind_gammas, drop = FALSE])
        }
    }
}

logLik.MixMod <- function (object, ...) {
    out <- object$logLik
    attr(out, "df") <- nrow(object$Hessian)
    attr(out, "nobs") <- length(unique(object$id))
    class(out) <- "logLik"
    out
}

coef.MixMod <- function (object, sub_model = c("main", "zero_part"), ...) {
    sub_model <- match.arg(sub_model)
    b <- ranef(object)
    RE_zi <- grep("zi_", colnames(b), fixed = TRUE)
    if (sub_model == "main") {
        betas <- fixef(object, sub_model = "main")
        if (length(RE_zi)) 
            b <- b[, -RE_zi, drop = FALSE]
        out <- matrix(betas, nrow = nrow(b), ncol = length(betas), byrow = TRUE)
        colnames(out) <- names(betas)
        rownames(out) <- rownames(b)
        out[, colnames(b)] <- out[, colnames(b)] + b
        out
    } else {
        gammas <- fixef(object, sub_model = "zero_part")
        if (length(RE_zi)) {
            b <- b[, RE_zi, drop = FALSE]
            colnames(b) <- gsub("zi_", "", colnames(b), fixed = TRUE)
            out <- matrix(gammas, nrow = nrow(b), ncol = length(gammas), byrow = TRUE)
            colnames(out) <- names(gammas)
            rownames(out) <- rownames(b)
            out[, colnames(b)] <- out[, colnames(b)] + b
            out
        } else {
            gammas
        }
        
    }
}

fixef.MixMod <- function(object, sub_model = c("main", "zero_part"), ...) {
    sub_model <- match.arg(sub_model)
    if (sub_model == "main") {
        object$coefficients
    } else {
        if (!is.null(object$gammas)) 
            object$gammas
        else
            stop("the fitted model does not have an extra zero-part.")
    }
}

ranef.MixMod <- function(object, post_vars = FALSE, ...) {
    out <- object$post_modes
    if (post_vars)
        attr(out, "post_vars") <- object$post_vars
    out
}

summary.MixMod <- function (object, sandwich = FALSE, ...) {
    betas <- fixef(object)
    n_betas <- length(betas)
    V <- vcov(object, sandwich = sandwich)
    var_betas <- V[seq_len(n_betas), seq_len(n_betas), drop = FALSE]
    ses <- sqrt(diag(var_betas))
    D <- object$D
    n_D <- length(D[lower.tri(D, TRUE)])
    coef_table <- cbind("Estimate" = betas, "Std.Err" = ses, "z-value" = betas / ses,
                        "p-value" = 2 * pnorm(abs(betas / ses), lower.tail = FALSE))
    if (!is.null(object$gammas)) {
        gammas <- object$gammas
        ind_gammas <- grep("zi_", colnames(V), fixed = TRUE)
        ses <- sqrt(diag(V[ind_gammas, ind_gammas, drop = FALSE]))
        coef_table_zi <- cbind("Estimate" = gammas, "Std.Err" = ses, "z-value" = gammas / ses,
                              "p-value" = 2 * pnorm(abs(gammas / ses), lower.tail = FALSE))
    }
    out <- list(coef_table = coef_table, 
                coef_table_zi = if (!is.null(object$gammas)) coef_table_zi,D = D, 
                logLik = logLik(object),
                AIC = AIC(object), BIC = BIC(object), call = object$call,
                N = length(object$id))
    if (!is.null(object$phis)) {
        phis <- object$phis
        ind_phis <- grep("phi_", colnames(V), fixed = TRUE)
        var_phis <- V[ind_phis, ind_phis, drop = FALSE]
        out$phis_table <- cbind("Estimate" = phis, "Std.Err" = sqrt(diag(var_phis)))
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
    coef_table <- as.data.frame(x[["coef_table"]])
    coef_table[1:3] <- lapply(coef_table[1:3], round, digits = digits)
    coef_table[["p-value"]] <- format.pval(coef_table[["p-value"]], eps = 1e-04)
    print(coef_table)
    if (!is.null(x[["coef_table_zi"]])) {
        cat("\nZero-part coefficients:\n")
        coef_table <- as.data.frame(x[["coef_table_zi"]])
        coef_table[1:3] <- lapply(coef_table[1:3], round, digits = digits)
        coef_table[["p-value"]] <- format.pval(coef_table[["p-value"]], eps = 1e-04)
        print(coef_table)
    }
    if (!is.null(x$phis_table)) {
        if (NB <- x$family$family %in% c("negative binomial", "zero-inflated negative binomial",
                                         "hurdle negative binomial")) {
            cat("\nlog(dispersion) parameter:\n")
        } else if (NB <- x$family$family %in% c("hurdle log-normal")) {
            cat("\nlog(residual std. dev.):\n")
        } else {
            cat("\nphi parameters:\n")
        }
        phis_table <- as.data.frame(x$phis_table)
        if (NB) 
            row.names(phis_table) <- " "
        phis_table[] <- lapply(phis_table, round, digits = digits)
        print(phis_table)
    }
    cat("\nIntegration:")
    cat("\nmethod: adaptive Gauss-Hermite quadrature rule")
    cat("\nquadrature points:", x$control$nAGQ)
    cat("\n\nOptimization:")
    methd <- if (x$control$iter_EM == 0) "quasi-Newton" 
    else if (isTRUE(attr(x$converged, "during_EM"))) "EM" 
    else "hybrid EM and quasi-Newton"
    cat("\nmethod:", methd)
    cat("\nconverged:", as.logical(x$converged), "\n")
    invisible(x)
}

coef.summary.MixMod <- function (object, ...) {
    object$coef_table
}

confint.MixMod <- function (object, parm = c("fixed-effects", "var-cov","extra", 
                                             "zero_part"), 
                            level = 0.95, sandwich = FALSE, ...) {
    parm <- match.arg(parm)
    V <- vcov(object, sandwich = sandwich)
    if (parm == "fixed-effects") {
        betas <- fixef(object)
        n_betas <- length(betas)
        ses_betas <- sqrt(diag(V[seq_len(n_betas), seq_len(n_betas), drop = FALSE]))
        out <- cbind(betas + qnorm((1 - level) / 2) * ses_betas, betas,
                     betas + qnorm((1 + level) / 2) * ses_betas)
    } else if (parm == "var-cov") {
        D <- object$D
        diag_D <- ncol(D) > 1 && all(abs(D[lower.tri(D)]) < sqrt(.Machine$double.eps))
        if (diag_D) {
            unconstr_D <- log(diag(D))
            n_betas <- length(object$coefficients)
            include <- seq(n_betas + 1, n_betas + length(unconstr_D))
            ses_unconstr_D <- sqrt(diag(V[include, include, drop = FALSE]))
            out <- cbind(unconstr_D + qnorm((1 - level) / 2) * ses_unconstr_D,
                         unconstr_D,
                         unconstr_D + qnorm((1 + level) / 2) * ses_unconstr_D)
            out <- exp(out)
            rownames(out) <- paste0("var.", rownames(D))
        } else {
            unconstr_D <- chol_transf(D)
            n_betas <- length(object$coefficients)
            include <- seq(n_betas + 1, n_betas + length(unconstr_D))
            ses_unconstr_D <- sqrt(diag(V[include, include, drop = FALSE]))
            out <- cbind(unconstr_D + qnorm((1 - level) / 2) * ses_unconstr_D,
                         unconstr_D,
                         unconstr_D + qnorm((1 + level) / 2) * ses_unconstr_D)
            ind <- lower.tri(D, TRUE)
            out[, 1] <- chol_transf(out[, 1])[ind]
            out[, 2] <- chol_transf(out[, 2])[ind]
            nams <- rownames(D)
            rownames(out) <- apply(which(ind, arr.ind = TRUE), 1, function (k) {
                if (k[1L] == k[2L]) paste0("var.", nams[k[1]]) else {
                    paste0("cov.", abbreviate(nams[k[2]], 5), "_", abbreviate(nams[k[1]], 5)) 
                }
            })
        }
    } else if (parm == "extra") {
        if (is.null(object$phis)) {
            stop("the model behind 'object' contains no extra (phis) parameters.\n")
        } else {
            phis <- object$phis
            ind_phis <- grep("phi_", colnames(V), fixed = TRUE)
            ses_phis <- sqrt(diag(V[ind_phis, ind_phis, drop = FALSE]))
            out <- cbind(phis + qnorm((1 - level) / 2) * ses_phis, phis,
                         phis + qnorm((1 + level) / 2) * ses_phis)
            if (object$family$family %in% c("negative binomial", 
                                            "zero-inflated negative binomial"))
                out <- exp(out)
        }
    } else {
        if (is.null(object$gammas)) {
            stop("the fitted model does not have an extra zero part.")
        } else {
            gammas <- object$gammas
            ind_gammas <- grep("zi_", colnames(V), fixed = TRUE)
            ses_gammas <- sqrt(diag(V[ind_gammas, ind_gammas, drop = FALSE]))
            out <- cbind(gammas + qnorm((1 - level) / 2) * ses_gammas, gammas,
                         gammas + qnorm((1 + level) / 2) * ses_gammas)
        }
    }
    colnames(out) <- c(paste(round(100 * c((1 - level) / 2, 
                                         (1 + level) / 2), 1), "%"), "Estimate")[c(1,3,2)]
    out
}

anova.MixMod <- function (object, object2, test = TRUE, L = NULL, 
                          sandwich = FALSE, ...) {
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
        if (L0 > L1) {
            L0 <- logLik(object2)
            L1 <- logLik(object)
        }
        nb0 <- attr(L0, "df")
        nb1 <- attr(L1, "df")
        df <- nb1 - nb0
        if (test && df == 0) {
            test <- FALSE
            warning("the two objects represent models with the same number of parameters;",
                    " argument 'test' is set to FALSE.")
        }
        fam <- object$family$family
        fam2 <- object2$family$family
        if (test &&  (fam != fam2 && ((fam != "poisson" & fam2 != "negative binomial") &&
                      (fam != "zero-inflated poisson" & fam2 != "zero-inflated negative binomial")))) {
            warning("it seems that the two objects represent model with different families;",
                    " are the models nested? If not, you should set 'test' to FALSE.")
        }
        out <- list(nam0 = deparse(substitute(object)), L0 = L0,
                    aic0 = AIC(object), bic0 = BIC(object),
                    nam1 = deparse(substitute(object2)), L1 = L1, aic1 = AIC(object2),
                    bic1 = AIC(object2), df = df, test = test)
        if (test) {
            LRT <- - 2 * (L0 - L1)
            attributes(LRT) <- NULL
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
        V <- vcov(object, sandwich = sandwich)
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
        RE_zi <- grep("zi_", colnames(b), fixed = TRUE)
        if (length(RE_zi))
            b <- b[, -RE_zi, drop = FALSE]
        Z <- model.matrix(object$Terms$termsZ, object$model_frames$mfZ)
        id <- match(object$id, unique(object$id))
        eta <- c(X %*% betas) + rowSums(Z * b[id, , drop = FALSE])
    } else {
        betas <- marginal_coefs(object, link_fun = link_fun)$betas
        eta <- c(X %*% betas)
    }
    if (!is.null(object$offset))
        eta <- eta + object$offset
    mu <- object$Funs$mu_fun(eta)
    if (!is.null(object$gammas)) {
        X_zi <- model.matrix(object$Terms$termsX_zi, object$model_frames$mfX_zi)
        offset_zi <- model.offset(object$model_frames$mfX_zi)
        gammas <- fixef(object, "zero_part")
        eta_zi <- c(X_zi %*% gammas)
        if (type == "subject_specific" && !is.null(object$Terms$termsZ_zi)) {
            b <- ranef(object)
            RE_zi <- grep("zi_", colnames(b), fixed = TRUE)
            if (length(RE_zi))
                b <- b[, RE_zi, drop = FALSE]
            Z_zi <- model.matrix(object$Terms$termsZ_zi, object$model_frames$mfZ_zi)
            id <- match(object$id, unique(object$id))
            eta_zi <- eta_zi + rowSums(Z_zi * b[id, , drop = FALSE])
        }
        if (!is.null(offset_zi))
            eta_zi <- eta_zi + offset_zi
        plogis(eta_zi, lower.tail = FALSE) * mu
    }
    names(mu) <- rownames(X)
    mu
}

residuals.MixMod <- function (object, type = c("mean_subject", "subject_specific",
                                               "marginal"), link_fun = NULL, 
                              tasnf_y = function (x) x, ...) {
    type <- match.arg(type)
    fits <- fitted(object, type = type, link_fun = link_fun)
    y <- model.response(object$model_frames$mfX)
    tasnf_y(y) - fits
}

marginal_coefs <- function (object, ...) UseMethod("marginal_coefs")

marginal_coefs.MixMod <- function (object, std_errors = FALSE, link_fun = NULL, 
                                   M = 3000, K = 100,
                                   seed = 1, cores = max(parallel::detectCores() - 1, 1), 
                                   sandwich = FALSE, ...) {
    if (!is.null(object$gammas)) {
        stop("marginal_coefs() is not yet implemented for models with an extra zero-part.")
    }
    X <- model.matrix(object$Terms$termsX, object$model_frames$mfX)
    Z <- model.matrix(object$Terms$termsZ, object$model_frames$mfZ)
    betas <- fixef(object)
    D <- object$D
    compute_marg_coefs <- function (object, X, betas, Z, D, M, link_fun, seed) {
        if (!exists(".Random.seed", envir = .GlobalEnv)) 
            runif(1)
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
        on.exit(assign(".Random.seed", RNGstate, envir = .GlobalEnv))
        mu_fun <- object$Funs$mu_fun
        if (is.null(link_fun)) {
            link_fun <- object$family$linkfun
        }
        if (is.null(link_fun)) {
            stop("you must specify the 'link_fun' argument.\n")
        }
        Xbetas <- c(X %*% betas)
        if (!is.null(object$offset)) {
            Xbetas <- Xbetas + object$offset
        }
        id <- match(object$id, unique(object$id))
        nRE <- ncol(Z)
        N <- nrow(X)
        n <- length(unique(id))
        eS <- eigen(D, symmetric = TRUE)
        ev <- eS$values
        V <- eS$vectors %*% diag(sqrt(pmax(ev, 0)), nRE)
        marg_inv_mu <- numeric(N)
        for (i in seq_len(n)) {
            set.seed(seed + i)
            id_i <- id == i
            b <- V %*% t(matrix(rnorm(M * nRE), M, nRE))
            Zb <- Z[id_i, , drop = FALSE] %*% b 
            mu <- mu_fun(Xbetas[id_i] + Zb)
            marg_inv_mu[id_i] <- link_fun(rowMeans(mu))
        }
        res <- c(solve(crossprod(X), crossprod(X, marg_inv_mu)))
        names(res) <- names(betas)
        res
    }
    out <- list(betas = compute_marg_coefs(object, X, betas, Z, D, M, link_fun, seed))
    if (std_errors) {
        blocks <- split(seq_len(K), rep(seq_len(cores), each = ceiling(K / cores),
                                        length.out = K))
        D <- object$D
        diag_D <- ncol(D) > 1 && all(abs(D[lower.tri(D)]) < sqrt(.Machine$double.eps))
        list_thetas <- list(betas = betas, D = if (diag_D) log(diag(D)) else chol_transf(D))
        tht <- unlist(as.relistable(list_thetas))
        V <- vcov(object, sandwich = sandwich)
        cluster_compute_marg_coefs <- function (block, tht, list_thetas, V, XX, Z, M,
                                                compute_marg_coefs, chol_transf,
                                                object, link_fun, seed) {
            if (!exists(".Random.seed", envir = .GlobalEnv)) 
                runif(1)
            RNGstate <- get(".Random.seed", envir = .GlobalEnv)
            on.exit(assign(".Random.seed", RNGstate, envir = .GlobalEnv))
            n_block <- length(block)
            m_betas <- matrix(0.0, n_block, length(list_thetas[["betas"]]))
            for (b in seq_along(block)) {
                seed. <- seed + block[b]
                set.seed(seed.)
                new_tht <- relist(MASS::mvrnorm(1, tht, V), skeleton = list_thetas)
                new_betas <- new_tht$betas
                new_D <- if (diag_D) diag(exp(new_tht$D), length(new_tht$D)) else chol_transf(new_tht$D)
                m_betas[b, ] <- compute_marg_coefs(object, XX, new_betas, Z, new_D, M,
                                                   link_fun, seed = seed.)
            }
            m_betas
        }
        cl <- parallel::makeCluster(cores)
        res <- parallel::parLapply(cl, blocks, cluster_compute_marg_coefs, tht = tht,
                                   list_thetas = list_thetas, V = V, XX = X, Z = Z, M = M,
                                   object = object, compute_marg_coefs = compute_marg_coefs,
                                   chol_transf = chol_transf, link_fun = link_fun, seed = seed)
        parallel::stopCluster(cl)
        out$var_betas <- var(do.call("rbind", res))
        dimnames(out$var_betas) <- list(names(out$betas), names(out$betas))
        ses <- sqrt(diag(out$var_betas))
        coef_table <- cbind("Estimate" = out$betas, "Std.Err" = ses,
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

effectPlotData.MixMod <- function (object, newdata, level = 0.95, marginal = FALSE, 
                                   K = 200, seed = 1, sandwich = FALSE, ...) {
    termsX <- delete.response(object$Terms$termsX)
    mfX <- model.frame(termsX, newdata, 
                       xlev = .getXlevels(termsX, object$model_frames$mfX))
    X <- model.matrix(termsX, mfX)
    if (is.null(object$gammas)) {
        if (marginal) {
            mcoefs <- marginal_coefs(object, std_errors = TRUE, ...)
            betas <- mcoefs$betas
            var_betas <- mcoefs$var_betas
        } else {
            betas <- fixef(object)
            n_betas <- length(betas)
            V <- vcov(object, sandwich = sandwich)
            var_betas <- V[seq_len(n_betas), seq_len(n_betas)]
        }
        pred <- c(X %*% betas)
        ses <- sqrt(diag(X %*% var_betas %*% t(X)))
        newdata$pred <- pred
        newdata$low <- pred + qnorm((1 - level) / 2) * ses
        newdata$upp <- pred + qnorm((1 + level) / 2) * ses
    } else {
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
            runif(1)
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
        mu_fun <- object$Funs$mu_fun
        betas <- fixef(object)
        gammas <- fixef(object, sub_model = "zero_part")
        termsX_zi <- object$Terms$termsX_zi
        mfX_zi <- model.frame(termsX_zi, newdata, 
                              xlev = .getXlevels(termsX_zi, object$model_frames$mfX_zi))
        X_zi <- model.matrix(termsX_zi, mfX_zi)
        list_thetas <- list(betas = betas, gammas = gammas)
        tht <- unlist(as.relistable(list_thetas))
        V <- vcov(object, sandwich = sandwich)
        ind <- c(seq_len(length(betas)), grep("zi_", colnames(V), fixed = TRUE))
        V <- V[ind, ind, drop = FALSE]
        new_tht <- MASS::mvrnorm(K, tht, V)
        if (marginal) {
            stop("the 'marginal = TRUE' option of effectPlotData() is not yet ", 
                 "implemented for models with an extra zero-part.")
            termsZ <- delete.response(object$Terms$termsZ)
            mfZ <- model.frame(termsZ, newdata, 
                               xlev = .getXlevels(termsZ, object$model_frames$mfZ))
            Z <- model.matrix(termsZ, mfZ)
            if (!is.null(object$Terms$termsZ_zi)) {
                termsZ_zi <- object$Terms$termsZ_zi
                mfZ_zi <- model.frame(termsZ_zi, newdata, 
                                      xlev = .getXlevels(termsZ_zi, object$model_frames$mfZ_zi))
                Z_zi <- model.matrix(termsZ_zi, mfZ_zi)
            }
        } else {
            eta_y <- c(X %*% betas)
            eta_zi <- c(X_zi %*% gammas)
            pred <- mu_fun(eta_y) / (1 + exp(eta_zi))
            Preds <- matrix(0.0, length(pred), K)
            for (k in seq_len(K)) {
                thetas_k <- relist(new_tht[k, ], skeleton = list_thetas)
                betas_k <- thetas_k$betas
                gammas_k <- thetas_k$gammas
                eta_y <- c(X %*% betas_k)
                eta_zi <- c(X_zi %*% gammas_k)
                Preds[, k] <- mu_fun(eta_y) / (1 + exp(eta_zi))
            }
            newdata$pred <- pred
            Qs <- apply(Preds, 1, quantile, probs = c((1 - level) / 2, (1 + level) / 2))
            newdata$low <- Qs[1, ]
            newdata$upp <- Qs[2, ]
        }
    }
    newdata
}

create_lists <- function (object, newdata) {
    if (!inherits(object, "MixMod")) {
        stop("only works for 'MixMod' objects.")
    }
    termsX <- delete.response(object$Terms$termsX)
    mfX <- model.frame(termsX, newdata, 
                       xlev = .getXlevels(termsX, object$model_frames$mfX))
    y <- model.response(model.frame(object$Terms$termsX, newdata))
    if (is.factor(y)) {
        if (family$family == "binomial")
            y <- as.numeric(y != levels(y)[1L])
        else
            stop("the response variable should not be a factor.\n")
    }
    offset <- model.offset(mfX)
    X <- model.matrix(termsX, mfX)
    termsZ <- delete.response(object$Terms$termsZ)
    mfZ <- model.frame(termsZ, newdata, 
                       xlev = .getXlevels(termsZ, object$model_frames$mfZ))
    Z <- model.matrix(termsZ, mfZ)
    if (!is.null(object$gammas)) {
        termsX_zi <- delete.response(object$Terms$termsX_zi)
        mfX_zi <- model.frame(termsX_zi, newdata, 
                              xlev = .getXlevels(termsX_zi, object$model_frames$mfX_zi))
        offset_zi <- model.offset(mfX_zi)
        X_zi <- model.matrix(termsX_zi, mfX_zi)
        if (!is.null(object$Terms$termsZ_zi)) {
            termsZ_zi <- delete.response(object$Terms$termsZ_zi)
            mfZ_zi <- model.frame(termsZ_zi, newdata, 
                                  xlev = .getXlevels(termsZ_zi, object$model_frames$mfZ_zi))
            Z_zi <- model.matrix(termsZ_zi, mfZ_zi)
        } else {
            Z_zi <- NULL
        }
    } else {
        X_zi <- offset_zi <- Z_zi <- NULL
    }
    id_nam <- object$id_name
    id <- newdata[[id_nam]]
    id <- match(id, unique(id))
    id_unq <- unique(id)
    y_lis <- if (is.matrix(y)) lapply(id_unq, function (i) y[id == i, ]) else split(y, id)
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
    X_zi_lis <- if (!is.null(X_zi)) lapply(id_unq, function (i) X_zi[id == i, , drop = FALSE])
    Z_zi_lis <- if (!is.null(Z_zi)) lapply(id_unq, function (i) Z_zi[id == i, , drop = FALSE])
    offset_zi_lis <- if (!is.null(offset_zi)) split(offset_zi, id)
    log_dens <- object$Funs$log_dens
    mu_fun <- object$Funs$mu_fun
    var_fun <- object$Funs$var_fun
    mu.eta_fun <- object$Funs$mu.eta_fun
    score_eta_fun <- object$Funs$score_eta_fun
    score_phis_fun <- object$Funs$score_phis_fun
    score_eta_zi_fun <- object$Funs$score_eta_zi_fun
    family <- object$family
    canonical <- !is.null(family$family) &&
        ((family$family == "binomial" && family$link == "logit") ||
             (family$family == "poisson" && family$link == "log"))
    known_families <- c("binomial", "poisson")
    user_defined <- !family$family %in% known_families
    numer_deriv <- if (object$control$numeric_deriv == "fd") fd else cd
    numer_deriv_vec <- if (object$control$numeric_deriv == "fd") fd_vec else cd_vec
    nRE <- if (!is.null(Z_zi)) ncol(Z) + ncol(Z_zi) else ncol(Z)
    start <- matrix(0.0, length(y_lis), nRE)
    betas <- fixef(object)
    invD <- solve(object$D)
    phis <- object$phis
    gammas <- object$gammas
    list(y_lis = y_lis, N_lis = N_lis, X_lis = X_lis, Z_lis = Z_lis, Z = Z, Z_zi = Z_zi,
         X_zi_lis = X_zi_lis, Z_zi_lis = Z_zi_lis, offset_zi_lis = offset_zi_lis,
         id = id, id_nam = id_nam, offset_lis = offset_lis, betas = betas, invD = invD, 
         phis = phis, gammas = gammas, start = start, canonical = canonical, 
         user_defined = user_defined, Zty_lis = Zty_lis, log_dens = log_dens, 
         mu_fun = mu_fun, var_fun = var_fun, 
         mu.eta_fun = mu.eta_fun, score_eta_fun = score_eta_fun, termsZ = termsZ,
         score_phis_fun = score_phis_fun, score_eta_zi_fun = score_eta_zi_fun)
}

predict.MixMod <- function (object, newdata, newdata2 = NULL, 
                            type_pred = c("response", "link"),
                            type = c("mean_subject", "subject_specific", "marginal", "zero_part"),
                            se.fit = FALSE, M = 300, df = 10, scale = 0.3, level = 0.95, 
                            seed = 1, return_newdata = FALSE, sandwich = FALSE, ...) {
    type_pred <- match.arg(type_pred)
    type <- match.arg(type)
    if (!is.null(object$gammas) && type != "zero_part" && type_pred == "link") {
        warning("for model with an extra zero-part only predictions at the level of the ",
                "response variable are returned;\n'type_pred' is set to 'response'.")
        type_pred <- "response"
    }
    termsX <- delete.response(object$Terms$termsX)
    mfX <- model.frame(termsX, newdata, 
                       xlev = .getXlevels(termsX, object$model_frames$mfX))
    X <- model.matrix(termsX, mfX)
    offset <- model.offset(mfX)
    if (!is.null(object$gammas)) {
        termsX_zi <- delete.response(object$Terms$termsX_zi)
        mfX_zi <- model.frame(termsX_zi, newdata, 
                              xlev = .getXlevels(termsX_zi, object$model_frames$mfX_zi))
        X_zi <- model.matrix(termsX_zi, mfX_zi)
        offset_zi <- model.offset(mfX_zi)
        gammas <- fixef(object, sub_model = "zero_part")
        eta_zi <- c(X_zi %*% gammas)
        if (!is.null(offset_zi)) {
            eta_zi <- eta_zi + offset_zi
        }
    }
    if (type %in% c("mean_subject", "marginal")) {
        if (type == "mean_subject") {
            betas <- fixef(object)
            var_betas <- vcov(object, parm = "fixed-effects", sandwich = sandwich)
            eta_y <- c(X %*% betas)
            if (!is.null(offset)) {
                eta_y <- eta_y + offset
            }
            pred <- if (type_pred == "link") eta_y else object$family$linkinv(eta_y)
            if (!is.null(object$gammas)) {
                pred <- plogis(eta_zi, lower.tail = FALSE) * pred
            }
            names(pred) <- row.names(newdata)
            se_fit <- if (se.fit && is.null(object$gammas)) sqrt(diag(X %*% var_betas %*% t(X)))
        } else {
            if (!is.null(object$gammas)) {
                stop("the predict() method is not yet implemented for models with an extra zero-part.")
            }
            mcoefs <- marginal_coefs(object, std_errors = TRUE, ...)
            betas <- fixef(object)
            var_betas <- mcoefs$var_betas
            pred <- if (type_pred == "link") c(X %*% betas) else object$family$linkinv(c(X %*% betas))
            names(pred) <- row.names(newdata)
            se_fit <- if (se.fit) sqrt(diag(X %*% var_betas %*% t(X)))
        }
    } else if (type == "zero_part") {
        if (is.null(object$gammas))
            stop("the fitted model does not have an extra zero part.")
        pred <- if (type_pred == "link") eta_zi else plogis(eta_zi)
        names(pred) <- row.names(newdata)
        var_gammas <- vcov(object, parm = "zero_part", sandwich = sandwich)
        se_fit <- if (se.fit) sqrt(diag(X_zi %*% var_gammas %*% t(X_zi)))
    } else {
        Lists <- create_lists(object, newdata)
        id <- Lists[["id"]]
        betas <- Lists[["betas"]]
        gammas <- Lists[["gammas"]]
        Z <- Lists[["Z"]]
        ncz <- ncol(Z)
        Z_zi <- Lists[["Z_zi"]]
        EBs <- find_modes(b = Lists$start, y_lis = Lists[["y_lis"]], 
                          N_lis = Lists[["N_lis"]], X_lis = Lists[["X_lis"]], 
                          Z_lis = Lists[["Z_lis"]], offset_lis = Lists[["offset_lis"]], 
                          X_zi_lis = Lists[["X_zi_lis"]], Z_zi_lis = Lists[["Z_zi_lis"]], 
                          offset_zi_lis = Lists[["offset_zi_lis"]], 
                          betas = betas, invD = Lists[["invD"]], phis = Lists[["phis"]], 
                          gammas = gammas, canonical = Lists[["canonical"]], 
                          user_defined = Lists[["user_defined"]], 
                          Zty_lis = Lists[["Zty_lis"]], log_dens = Lists[["log_dens"]], 
                          mu_fun = Lists[["mu_fun"]], var_fun = Lists[["var_fun"]], 
                          mu.eta_fun = Lists[["mu.eta_fun"]], 
                          score_eta_fun = Lists[["score_eta_fun"]], 
                          score_phis_fun = Lists[["score_phis_fun"]], 
                          score_eta_zi_fun = Lists[["score_eta_zi_fun"]])
        eta <- c(X %*% betas) + rowSums(Z * EBs$post_modes[id, seq_len(ncz), drop = FALSE])
        if (!is.null(offset))
            eta <- eta + offset
        if (!is.null(object$Terms$termsZ_zi)) {
            eta_zi <- eta_zi + rowSums(Z_zi * EBs$post_modes[id, -seq_len(ncz), drop = FALSE])
        }
        pred <- if (type_pred == "link") eta else object$family$linkinv(eta)
        if (!is.null(object$gammas)) {
            pred <- plogis(eta_zi, lower.tail = FALSE) * pred
            attr(pred, "zi_probs") <- plogis(eta_zi)
        }
        names(pred) <- row.names(newdata)
        if (se.fit) {
            if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
                runif(1)
            R.seed <- get(".Random.seed", envir = .GlobalEnv)
            set.seed(seed)
            RNGstate <- structure(seed, kind = as.list(RNGkind()))
            on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
            log_post_b <- function (b_i, y_i, X_i, Z_i, offset_i, X_zi_i, Z_zi_i, offset_zi_i,
                                    betas, invD, phis, gammas, log_dens, mu_fun) {
                ncz <- ncol(Z_i)
                eta_y <- as.vector(X_i %*% betas + Z_i %*% b_i[seq_len(ncz)])
                if (!is.null(offset_i))
                    eta_y <- eta_y + offset_i
                if (!is.null(X_zi_i)) {
                    eta_zi <- as.vector(X_zi_i %*% gammas)
                    if (!is.null(Z_zi_i))
                        eta_zi <- eta_zi + as.vector(Z_zi_i %*% b_i[-seq_len(ncz)])
                    if (!is.null(offset_zi_i))
                        eta_zi <- eta_zi + offset_zi_i
                }
                sum(log_dens(y_i, eta_y, mu_fun, phis, eta_zi), na.rm = TRUE) -
                    c(0.5 * crossprod(b_i, invD) %*% b_i)
            }
            log_dens <- object$Funs$log_dens
            mu_fun <- object$Funs$mu_fun
            calc_alpha <- function (log_post_new, log_post_old, log_prop_new, 
                                    log_prop_old) {
                min(exp(log_post_new + log_prop_old - log_post_old - log_prop_new), 1)
            }
            phis <- object$phis
            D <- object$D
            diag_D <- all(abs(D[lower.tri(D)]) < sqrt(.Machine$double.eps))
            list_thetas <- list(betas = betas, D = if (diag_D) log(diag(D)) else chol_transf(D))
            if (!is.null(phis)) {
                list_thetas <- c(list_thetas, list(phis = phis))
            }
            if (!is.null(gammas)) {
                list_thetas <- c(list_thetas, list(gammas = gammas))
            }
            tht <- unlist(as.relistable(list_thetas))
            V <- vcov(object, sandwich = sandwich)
            tht_new <- MASS::mvrnorm(M, tht, V)
            row_split_ind <- row(EBs$post_modes)
            mu <- split(EBs$post_modes, row_split_ind)
            Sigma <- lapply(EBs$post_hessians, solve)
            scale <- rep(scale, length.out = length(Sigma))
            Sigma <- mapply("*", scale, Sigma, SIMPLIFY = FALSE)
            EBs_proposed <- mapply(rmvt, mu = mu, Sigma = Sigma, SIMPLIFY = FALSE,
                                   MoreArgs = list(n = M, df = df))
            dmvt_proposed <- mapply(dmvt, x = EBs_proposed, mu = mu, Sigma = Sigma,
                                    MoreArgs = list(df = df, log = TRUE, prop = FALSE),
                                    SIMPLIFY = FALSE)
            b_current <- mu
            dmvt_current <- mapply(dmvt, x = mu, mu = mu, Sigma = Sigma, SIMPLIFY = FALSE,
                                   MoreArgs = list(df = df, log = TRUE, prop = FALSE))
            y_lis <- Lists[["y_lis"]]
            X_lis <- Lists[["X_lis"]]
            Z_lis <- Lists[["Z_lis"]]
            offset_lis <- Lists[["offset_lis"]]
            if (is.null(offset_lis))
                offset_lis <- rep(list(NULL), length(y_lis))
            X_zi_lis <- Lists[["X_zi_lis"]]
            if (is.null(X_zi_lis))
                X_zi_lis <- rep(list(NULL), length(y_lis))
            Z_zi_lis <- Lists[["Z_zi_lis"]]
            if (is.null(Z_zi_lis))
                Z_zi_lis <- rep(list(NULL), length(y_lis))
            offset_zi_lis <- Lists[["offset_zi_lis"]]
            if (is.null(offset_zi_lis))
                offset_zi_lis <- rep(list(NULL), length(y_lis))
            n <- length(pred)
            Preds <- matrix(0.0, n, M)
            b <- vector("list", M)
            if (!is.null(object$gammas)) {
                zi_probs <- matrix(0.0, n, M)
            }
            success_rate <- matrix(FALSE, M, length(y_lis))
            for (m in seq_len(M)) {
                # Extract simulared new parameter values
                new_pars <- relist(tht_new[m, ], skeleton = list_thetas)
                betas_m <- new_pars$betas
                phis_m <- new_pars$phis
                gammas_m <- new_pars$gammas
                D_m <- if (diag_D) diag(exp(new_pars$D), length(new_pars$D)) else chol_transf(new_pars$D)
                invD_m <- solve(D)
                # Simulate new EBs
                log_post_b_current <- mapply(log_post_b, b_i = b_current, y_i = y_lis, 
                                             X_i = X_lis, Z_i = Z_lis, offset_i = offset_lis,
                                             X_zi_i = X_zi_lis, Z_zi_i = Z_zi_lis, 
                                             offset_zi_i = offset_zi_lis,
                                             MoreArgs = list(betas = betas_m, invD = invD_m, 
                                                             phis = phis_m, gammas = gammas_m,
                                                             log_dens = log_dens, 
                                                             mu_fun = mu_fun),
                                             SIMPLIFY = FALSE)
                b_new <- lapply(EBs_proposed, function (x, m) x[m, ], m = m)
                log_post_b_new <- mapply(log_post_b, b_i = b_new, y_i = y_lis, X_i = X_lis, 
                                         Z_i = Z_lis, offset_i = offset_lis, 
                                         X_zi_i = X_zi_lis, Z_zi_i = Z_zi_lis, 
                                         offset_zi_i = offset_zi_lis,
                                         MoreArgs = list(betas = betas_m, invD = invD_m, 
                                                         phis = phis_m, gammas = gammas_m,
                                                         log_dens = log_dens, 
                                                         mu_fun = mu_fun),
                                         SIMPLIFY = FALSE)
                alphas <- mapply(calc_alpha, log_post_b_new, log_post_b_current, dmvt_current, 
                                 lapply(dmvt_proposed, "[", m))
                keep_ind <- runif(length(alphas)) <= alphas
                if (any(keep_ind)) {
                    b_current[keep_ind] <- b_new[keep_ind]
                    dmvt_current[keep_ind] <- lapply(dmvt_proposed, "[", m)[keep_ind]
                }
                success_rate[m, ] <- keep_ind
                # Calculate Predictions
                b[[m]] <- do.call("rbind", b_current)
                eta <- c(X %*% betas_m) + rowSums(Z * b[[m]][id, seq_len(ncz), drop = FALSE])
                if (!is.null(offset))
                    eta <- eta + offset
                Preds[, m] <- if (type_pred == "link") eta else object$family$linkinv(eta)
                if (!is.null(object$gammas)) {
                    eta_zi <- as.vector(X_zi %*% gammas_m)
                    if (!is.null(object$Terms$termsZ_zi)) {
                        eta_zi <- eta_zi + rowSums(Z_zi * b[[m]][id, -seq_len(ncz), drop = FALSE])
                    }
                    if (!is.null(offset_zi))
                        eta_zi <- eta_zi + offset_zi
                    Preds[, m] <- plogis(eta_zi, lower.tail = FALSE) * Preds[, m]
                    zi_probs[, m] <- plogis(eta_zi)
                }
            }
            se_fit <- apply(Preds, 1, sd, na.rm = TRUE)
            Qs <- apply(Preds, 1, quantile, 
                        probs = c((1 - level) / 2, (1 + level) / 2))
            low <- Qs[1, ]
            upp <- Qs[2, ]
            if (!is.null(gammas)) {
                Qs_zi <- apply(zi_probs, 1, quantile, 
                            probs = c((1 - level) / 2, (1 + level) / 2))
                
                attr(low, "zi_probs") <- Qs_zi[1, ]
                attr(upp, "zi_probs") <- Qs_zi[2, ]
            }
            names(se_fit) <- names(low) <- names(upp) <- names(pred)
        }
        if (!is.null(newdata2)) {
            id_nam <- Lists[["id_nam"]]
            id2 <- newdata2[[id_nam]]
            id2 <- match(id2, unique(newdata[[id_nam]]))
            mfX2 <- model.frame(termsX, newdata2, 
                                xlev = .getXlevels(termsX, object$model_frames$mfX))
            X2 <- model.matrix(termsX, mfX2)
            offset2 <- model.offset(mfX2)
            termsZ <- Lists[["termsZ"]]
            mfZ2 <- model.frame(termsZ, newdata2, 
                                xlev = .getXlevels(termsZ, object$model_frames$mfZ))
            Z2 <- model.matrix(termsZ, mfZ2)
            ncz <- ncol(Z2)
            eta2 <- c(X2 %*% betas) + rowSums(Z2 * EBs$post_modes[id2, seq_len(ncz), drop = FALSE])
            if (!is.null(offset2)) {
                eta2 <- eta2 + offset2
            }
            pred2 <- if (type_pred == "link") eta2 else object$family$linkinv(eta2)
            if (!is.null(object$gammas)) {
                mfX2_zi <- model.frame(termsX_zi, newdata2, 
                                    xlev = .getXlevels(termsX_zi, object$model_frames$mfX_zi))
                X2_zi <- model.matrix(termsX_zi, mfX2_zi)
                offset2_zi <- model.offset(mfX2_zi)
                eta2_zi <- c(X2_zi %*% gammas)
                if (!is.null(object$Terms$termsZ_zi)) {
                    termsZ_zi <- object$Terms$termsZ_zi
                    mfZ2_zi <- model.frame(termsZ_zi, newdata2, 
                                           xlev = .getXlevels(termsZ_zi, object$model_frames$mfZ_zi))
                    Z2_zi <- model.matrix(termsZ_zi, mfZ2_zi)
                    eta2_zi <- eta2_zi + rowSums(Z2_zi * EBs$post_modes[id2, -seq_len(ncz), drop = FALSE])
                }
                if (!is.null(offset2_zi))
                    eta2_zi <- eta2_zi + offset2_zi
                pred2 <- plogis(eta2_zi, lower.tail = FALSE) * pred2
                attr(pred2, "zi_probs") <- plogis(eta2_zi)
            }
            names(pred2) <- row.names(newdata2)
            if (se.fit) {
                Preds2 <- matrix(0.0, length(pred2), M)
                if (!is.null(gammas)) {
                    zi_probs2 <- matrix(0.0, length(pred2), M)
                }
                for (m in seq_len(M)) {
                    new_pars <- relist(tht_new[m, ], skeleton = list_thetas)
                    betas_m <- new_pars$betas
                    b_m <- b[[m]]
                    gammas_m <- new_pars$gammas
                    eta2 <- c(X2 %*% betas_m) + rowSums(Z2 * b_m[id2, seq_len(ncz), drop = FALSE])
                    if (!is.null(offset2))
                        eta2 <- eta2 + offset2
                    Preds2[, m] <- if (type_pred == "link") eta2 else object$family$linkinv(eta2)
                    if (!is.null(gammas)) {
                        eta2_zi <- c(X2_zi %*% gammas_m)
                        if (!is.null( object$Terms$termsZ_zi))
                            eta2_zi <- eta2_zi + rowSums(Z2_zi * b_m[id2, -seq_len(ncz), drop = FALSE])
                        if (!is.null(offset2_zi))
                            eta2_zi <- eta2_zi + offset2_zi
                        Preds2[, m] <- plogis(eta2_zi, lower.tail = FALSE) * Preds2[, m]
                        zi_probs2[, m] <- plogis(eta2_zi)
                    }
                }
                se_fit2 <- apply(Preds2, 1, sd, na.rm = TRUE)
                Qs2 <- apply(Preds2, 1, quantile, probs = c((1 - level) / 2, (1 + level) / 2))
                low2 <- Qs2[1, ]
                upp2 <- Qs2[2, ]
                if (!is.null(gammas)) {
                    Qs2_zi <- apply(zi_probs2, 1, quantile, 
                                   probs = c((1 - level) / 2, (1 + level) / 2))
                    
                    attr(low2, "zi_probs") <- Qs2_zi[1, ]
                    attr(upp2, "zi_probs") <- Qs2_zi[2, ]
                }
                names(se_fit2) <- names(low2) <- names(upp2) <- names(pred2)
            }
        }
    }
    if (return_newdata) {
        na_exclude <- attr(mfX, "na.action")
        if (!is.null(object$gammas)) { 
            na_exclude <- union(attr(mfX_zi, "na.action"), na_exclude)
        }
        if (!is.null(na_exclude))
            newdata <- newdata[-na_exclude, ]
        newdata$pred <- pred
        if (se.fit && type == "subject_specific") {
            newdata$low <- low
            newdata$upp <- upp
        }
        if (!is.null(newdata2)) {
            na_exclude2 <- attr(mfX2, "na.action")
            if (!is.null(object$gammas)) { 
                na_exclude2 <- union(attr(mfX2_zi, "na.action"), na_exclude)
            }
            if (!is.null(na_exclude))
                newdata2 <- newdata2[-na_exclude2, ]
            newdata2$pred <- pred2
            if (se.fit && type == "subject_specific") {
                newdata2$low <- low2
                newdata2$upp <- upp2
            }
            return(list(newdata = newdata, newdata2 = newdata2))
        } else {
            return(newdata)
        }
    } else {
        if (se.fit) {
            if (is.null(newdata2)) {
                if (type == "subject_specific") 
                    list(pred = pred, se.fit = se_fit, low = low, upp = upp,
                         success_rate = colMeans(success_rate))
                else 
                    list(pred = pred, se.fit = se_fit)
            } else {
                if (type == "subject_specific")
                    list(pred = pred, pred2 = pred2, se.fit = se_fit, se.fit2 = se_fit2,
                         low = low, upp = upp, low2 = low2, upp2 = upp2)
                else
                    list(pred = pred, pred2 = pred2, se.fit = se_fit, se.fit2 = se_fit2)
            }
        } else {
            if (is.null(newdata2)) pred else list(pred = pred, pred2 = pred2)
        }
    }
}

simulate.MixMod <- function (object, nsim = 1, seed = NULL, 
                             type = c("subject_specific", "mean_subject"),
                             acount_MLEs_var = FALSE, sim_fun = NULL, 
                             sandwich = FALSE, ...) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) 
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    type <- match.arg(type)
    if (is.null(sim_fun)) {
        if (object$family$family == "binomial") {
            N <- if ((y <- NCOL(model.response(object$model_frames$mfX))) == 2) 
                y[, 1] + y[, 2] else 1
            .N <- N
            env <- new.env(parent = .GlobalEnv)
            assign(".N", N, envir = env)
            sim_fun <- function (n, mu, phis, eta_zi) {
                rbinom(n = n, size = .N, prob = mu)
            }
            environment(sim_fun) <- env
        } else if (object$family$family == "poisson") {
            sim_fun <- function (n, mu, phis, eta_zi) {
                rpois(n = n, lambda = mu)
            }
        } else if (object$family$family == "negative binomial") {
            sim_fun <- function (n, mu, phis, eta_zi) {
                rnbinom(n = n, size = exp(phis), mu = mu)
            }
        } else if (object$family$family == "zero-inflated poisson") {
            sim_fun <- function (n, mu, phis, eta_zi) {
                out <- rpois(n = n, lambda = mu)
                out[as.logical(rbinom(n, 1, plogis(eta_zi)))] <- 0
                out
            }
        } else if (object$family$family == "zero-inflated negative binomial") {
            sim_fun <- function (n, mu, phis, eta_zi) {
                out <- rnbinom(n = n, size = exp(phis), mu = mu)
                out[as.logical(rbinom(n, 1, plogis(eta_zi)))] <- 0
                out
            }
        } else if (!is.null(object$family$simulate) && is.function(object$family$simulate)) {
            sim_fun <- object$family$simulate
        } else {
            stop("'sim_fun()' needs to be specified; check the help page.")
        }
    }
    id <- object$id
    id <- match(id, unique(id))
    n <- length(unique(id))
    X <- model.matrix(object$Terms$termsX, object$model_frames$mfX)
    Z <- model.matrix(object$Terms$termsZ, object$model_frames$mfZ)
    offset <- model.offset(object$model_frames$mfX)
    if (has_X_Zi <- !is.null(object$Terms$termsX_zi)) {
        X_zi <- model.matrix(object$Terms$termsX_zi, object$model_frames$mfX_zi)
        offset_zi <- model.offset(object$model_frames$mfX_zi)
    }
    if (has_Z_Zi <- !is.null(object$Terms$termsZ_zi)) {
        Z_zi <- model.matrix(object$Terms$termsZ_zi, object$model_frames$mfZ_zi)
    }
    betas <- fixef(object)
    D <- object$D
    gammas <- object$gammas
    phis <- object$phis
    diag_D <- all(abs(D[lower.tri(D)]) < sqrt(.Machine$double.eps))
    nRE <- ncol(D)
    ind <- vector("logical", nRE)
    ind[grep("zi_", colnames(D), fixed = TRUE, invert = TRUE)] <- TRUE
    if (acount_MLEs_var) {
        list_thetas <- list(betas = betas, 
                            D = if (diag_D) log(diag(D)) else chol_transf(D))
        if (!is.null(phis)) {
            list_thetas <- c(list_thetas, list(phis = phis))
        }
        if (!is.null(gammas)) {
            list_thetas <- c(list_thetas, list(gammas = gammas))
        }
        tht <- unlist(as.relistable(list_thetas))
        new_thetas <- MASS::mvrnorm(nsim, tht, vcov(object, sandwich = sandwich))
    }
    out <- matrix(0.0, nrow(X), nsim)
    for (i in seq_len(nsim)) {
        if (acount_MLEs_var) {
            new_thetas_i <- relist(new_thetas[i, ], skeleton = list_thetas)
            betas <- new_thetas_i$betas
            phis <- new_thetas_i$phis
            gammas <- new_thetas_i$gammas
            D <- if (diag_D) diag(exp(new_thetas_i$D), length(new_thetas_i$D)) 
            else chol_transf(new_thetas_i$D)
        }
        b_i <- MASS::mvrnorm(n, rep(0, nRE), D)
        if (type == "mean_subject")
            b_i <- b_i * 0
        eta_y <- c(X %*% betas) + rowSums(Z * b_i[id, ind, drop = FALSE])
        if (!is.null(offset))
            eta_y <- eta_y + offset
        mu <- object$Funs$mu_fun(eta_y)
        if (has_X_Zi)
            eta_zi <- c(X_zi %*% gammas)
        if (has_Z_Zi)
            eta_zi <- eta_zi + rowSums(Z_zi * b_i[id, !ind, drop = FALSE])
        if (has_X_Zi && !is.null(offset_zi))
            eta_zi <- eta_zi + offset_zi
        out[, i] <- sim_fun(nrow(X), mu, phis, eta_zi)
    }
    out
}

model.matrix.MixMod <- function (object, type = c("fixed", "random"), ...) {
    type <- match.arg(type)
    if (type == "fixed") {
        model.matrix(object$Terms$termsX, object$model_frames$mfX)
    } else {
        model.matrix(object$Terms$termsZ, object$model_frames$mfZ)
    }
}

model.frame.MixMod <- function (formula, type = c("fixed", "random"), ...) {
    type <- match.arg(type)
    if (type == "fixed") {
        formula$model_frames$mfX
    } else {
        formula$model_frames$mfZ
    }
}

terms.MixMod <- function (x, type = c("fixed", "random"), ...) {
    type <- match.arg(type)
    if (type == "fixed") {
        x$Terms$termsX
    } else {
        x$Terms$termsZ
    }
}

recover_data.MixMod <- function (object, ...) {
    fcall <- object$call
    emmeans::recover_data(fcall, delete.response(terms(object)), object$na.action, ...)
}

emm_basis.MixMod <- function (object, trms, xlev, grid, ...) { 
    m <- model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X <- model.matrix(trms, m, contrasts.arg = object$contrasts) 
    bhat <- fixef(object) 
    V <- vcov(object, parm = "fixed-effects")
    nbasis <- matrix(NA) 
    dfargs <- list(df = Inf)
    dffun <- function (k, dfargs) dfargs$df
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, dfargs = dfargs)
}

scoring_rules <- function (object, newdata, newdata2 = NULL, max_count = 2000, 
                           return_newdata = FALSE) {
    termsX <- object$Terms$termsX
    ND <- if (is.null(newdata2)) newdata else newdata2
    y <- model.response(model.frame(termsX, data = ND,
                                    xlev = .getXlevels(termsX, object$model_frames$mfX)))
    if (is.null(y)) {
        stop("the outcome variable is not in 'newdata' and/or 'newdata2'.")
    }
    n <- length(y)
    if (object$family$family == "binomial" && NCOL(y) == 2) {
        N <- max_count <- y[, 1] + y[, 2]
        y <- y[, 1]
    } else if (object$family$family == "binomial" && NCOL(y) == 1) {
        N <- max_count <- rep(1, n)
    } else {
        N <- NULL
    }
    max_count <- rep(max_count, length.out = n)
    prob_fun <- if (object$family$family == "binomial") {
        function (x, mean, pis, N) dbinom(x, size = N, prob = mean)
    } else if (object$family$family == "poisson") {
        function (x, mean, pis, N) dpois(x, lambda = mean)
    } else if (object$family$family == "negative binomial") {
        function (x, mean, pis, N) dnbinom(x, mu = mean, size = exp(object$phis))
    } else if (object$family$family == "zero-inflated poisson") {
        function (x, mean, pis, N) {
            ind0 <- x == 0
            out <- (1 - pis) * dpois(x, lambda = mean)
            out[ind0] <- pis + out[ind0]
            out
        }
    } else if (object$family$family == "zero-inflated negative binomial") {
        function (x, mean, pis, N) {
            ind0 <- x == 0
            out <- (1 - pis) * dnbinom(x, mu = mean, size = exp(object$phis))
            out[ind0] <- pis + out[ind0]
            out
        }
    } else if (object$family$family == "hurdle poisson") {
        function (x, mean, pis, N) {
            ind0 <- x == 0
            trunc_zero <- dpois(x, lambda = mean) / 
                ppois(0, lambda = mean, lower.tail = FALSE)
            out <- (1 - pis) * trunc_zero
            out[ind0] <- pis
            out
        }
    } else if (object$family$family == "hurdle negative binomial") {
        function (x, mean, pis, N) {
            ind0 <- x == 0
            trunc_zero <- dnbinom(x, mu = mean, size = exp(object$phis)) / 
                pnbinom(0, mu = mean, size = exp(object$phis), lower.tail = FALSE)
            out <- (1 - pis) * trunc_zero
            out[ind0] <- pis
            out
        }
    }
    max_count_seq <- lapply(max_count, seq, from = 0)
    pred <- predict(object, newdata = newdata, newdata2 = newdata2, 
                    type = "subject_specific")
    pred_zi <- if (!is.null(object$gammas)) attr(pred, "zi_probs")
    if (!is.null(newdata2)) {
        pred <- pred$pred2
        pred_zi <- attr(pred, "zi_probs")
    }
    logarithmic <- quadratic <- spherical <- numeric(n)
    for (i in seq_len(n)) {
        p_y <- prob_fun(y[i], mean = pred[i], pis = pred_zi[i], N[i])
        quadrat_p <- sum(prob_fun(max_count_seq[[i]], mean = pred[i], 
                                  pis = pred_zi[i], N[i])^2)
        logarithmic[i] <- log(p_y)
        quadratic[i] <- 2 * p_y + quadrat_p
        spherical[i] <- p_y / sqrt(quadrat_p)
    }
    result <- data.frame(logarithmic = logarithmic, quadratic = quadratic, 
                         spherical = spherical)
    if (return_newdata) cbind(ND, result) else result
}


