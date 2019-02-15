find_modes <- function (b, y_lis, N_lis, X_lis, Z_lis, offset_lis, X_zi_lis, Z_zi_lis, 
                        offset_zi_lis, betas, invD, phis, gammas, 
                        canonical, user_defined, Zty_lis, log_dens, mu_fun, var_fun, 
                        mu.eta_fun, score_eta_fun, score_phis_fun, score_eta_zi_fun) {
    log_post_b <- function (b_i, y_i, N_i, X_i, Z_i, offset_i, X_zi_i, Z_zi_i, offset_zi_i, 
                            betas, invD, phis, gammas, canonical,
                            user_defined, Zty_i, log_dens, mu_fun, var_fun, mu.eta_fun,
                            score_eta_fun, score_phis_fun, score_eta_zi_fun) {
        ind_Z <- seq_len(ncol(Z_i))
        eta_y <- as.vector(X_i %*% betas + Z_i %*% b_i[ind_Z])
        if (!is.null(offset_i))
            eta_y <- eta_y + offset_i
        eta_zi <- if (!is.null(X_zi_i)) as.vector(X_zi_i %*% gammas)
        if (!is.null(Z_zi_i))
            eta_zi <- eta_zi + as.vector(Z_zi_i %*% b_i[-ind_Z])
        if (!is.null(offset_zi_i))
            eta_zi <- eta_zi + offset_zi_i
        - sum(log_dens(y_i, eta_y, mu_fun, phis, eta_zi), na.rm = TRUE) +
                c(0.5 * crossprod(b_i, invD) %*% b_i)
    }
    score_log_post_b <- function (b_i, y_i, N_i, X_i, Z_i, offset_i,X_zi_i, Z_zi_i, offset_zi_i,
                                  betas, invD, phis, gammas,
                                  canonical, user_defined, Zty_i, log_dens, mu_fun,
                                  var_fun, mu.eta_fun, score_eta_fun, score_phis_fun,
                                  score_eta_zi_fun) {
        eta_y <- as.vector(X_i %*% betas + Z_i %*% b_i[seq_len(ncol(Z_i))])
        if (!is.null(offset_i))
            eta_y <- eta_y + offset_i
        eta_zi <- if (!is.null(X_zi_i)) as.vector(X_zi_i %*% gammas)
        if (!is.null(Z_zi_i))
            eta_zi <- eta_zi + as.vector(Z_zi_i %*% b_i[-seq_len(ncol(Z_i))])
        if (!is.null(offset_zi_i))
            eta_zi <- eta_zi + offset_zi_i
        mu_y <- mu_fun(eta_y)
        log_dens_part <- if (user_defined) {
            out <- if (!is.null(score_eta_fun)) {
                - crossprod(Z_i, score_eta_fun(y_i, mu_y, phis, eta_zi))
            } else {
                l1 <- log_dens(y_i, eta_y + 1e-04, mu_fun, phis, eta_zi)
                l2 <- log_dens(y_i, eta_y - 1e-04, mu_fun, phis, eta_zi)
                - crossprod(Z_i, (l1 - l2) / (2 * 1e-04))
            }
            if (!is.null(Z_zi_i)) {
                out <- if (!is.null(score_eta_zi_fun)) {
                    c(out, - crossprod(Z_zi_i, score_eta_zi_fun(y_i, mu_y, phis, eta_zi)))
                } else {
                    l1 <- log_dens(y_i, eta_y, mu_fun, phis, eta_zi + 1e-04)
                    l2 <- log_dens(y_i, eta_y, mu_fun, phis, eta_zi - 1e-04)
                    c(out, - crossprod(Z_zi_i, (l1 - l2) / (2 * 1e-04)))
                }
            }
            out
        } else {
            if (canonical) {
                if (!is.null(N_i))- Zty_i + crossprod(Z_i, N_i * mu_y) else 
                    - Zty_i + crossprod(Z_i, mu_y)
            } else {
                var <- var_fun(mu_y)
                deriv <- mu.eta_fun(eta_y)
                if (!is.null(N_i)) - crossprod(Z_i, (y_i[, 1] - N_i * mu_y) * deriv / var) else
                - crossprod(Z_i, (y_i - mu_y) * deriv / var)
            }
        }
        drop(log_dens_part + invD %*% b_i)
    }
    n <- length(y_lis)
    post_modes <- b
    post_hessians <- vector("list", n)
    for (i in seq_len(n)) {
        y_i <- y_lis[[i]]
        N_i <- if (!is.null(N_lis)) N_lis[[i]]
        X_i <- X_lis[[i]]
        Z_i <- Z_lis[[i]]
        offset_i <- if (!is.null(offset_lis)) offset_lis[[i]]
        Zty_i <- Zty_lis[[i]]
        X_zi_i <- if (!is.null(X_zi_lis)) X_zi_lis[[i]]
        Z_zi_i <- if (!is.null(Z_zi_lis)) Z_zi_lis[[i]]
        offset_zi_i <- if (!is.null(offset_zi_lis)) offset_zi_lis[[i]]
        b_i <- b[i, ]
        opt_i <- optim(par = b_i, fn = log_post_b, gr = score_log_post_b, method = "BFGS",
                       y_i = y_i, N_i = N_i, X_i = X_i, Z_i = Z_i, offset_i = offset_i,
                       X_zi_i = X_zi_i, Z_zi_i = Z_zi_i, offset_zi_i = offset_zi_i,
                       betas = betas, invD = invD, phis = phis, gammas = gammas, 
                       canonical = canonical,
                       user_defined = user_defined, Zty_i = Zty_i, log_dens = log_dens,
                       mu_fun = mu_fun, var_fun = var_fun, mu.eta_fun = mu.eta_fun,
                       score_eta_fun = score_eta_fun, score_phis_fun = score_phis_fun,
                       score_eta_zi_fun = score_eta_zi_fun)
        post_modes[i, ] <- opt_i$par
        post_hessians[[i]] <- cd_vec(post_modes[i, ], score_log_post_b,
                                     y_i = y_i, N_i = N_i, X_i = X_i, Z_i = Z_i, 
                                     offset_i = offset_i, X_zi_i = X_zi_i, Z_zi_i = Z_zi_i, 
                                     offset_zi_i = offset_zi_i, betas = betas, invD = invD, 
                                     phis = phis, gammas = gammas, canonical = canonical, 
                                     user_defined = user_defined,
                                     Zty_i = Zty_i, log_dens = log_dens, mu_fun = mu_fun,
                                     var_fun = var_fun, mu.eta_fun = mu.eta_fun,
                                     score_eta_fun = score_eta_fun,
                                     score_phis_fun = score_phis_fun, 
                                     score_eta_zi_fun = score_eta_zi_fun)
    }
    list(post_modes = post_modes, post_hessians = post_hessians)
}

GHfun <- function (b, y_lis, N_lis, X_lis, Z_lis, offset_lis, X_zi_lis, Z_zi_lis, offset_zi_lis,
                   betas, inv_D, phis, gammas, k, q,
                   canonical, user_defined, Zty_lis, log_dens, mu_fun, var_fun, mu.eta_fun,
                   score_eta_fun, score_phis_fun, score_eta_zi_fun) {
    GH <- gauher(k)
    aGH <- find_modes(b, y_lis, N_lis, X_lis, Z_lis, offset_lis, X_zi_lis, Z_zi_lis, 
                      offset_zi_lis, betas, inv_D, phis, gammas,
                      canonical, user_defined, Zty_lis, log_dens, mu_fun, var_fun, 
                      mu.eta_fun, score_eta_fun, score_phis_fun, score_eta_zi_fun)
    modes <- aGH$post_modes
    chol_hessians <- lapply(aGH$post_hessian, chol)
    b <- as.matrix(expand.grid(lapply(seq_len(q), function (k, u) u$x, u = GH)))
    n <- nrow(modes)
    b_new <- vector("list", n)
    dets <- numeric(n)
    for (i in seq_len(n)) {
        b_new[[i]] <- t(sqrt(2) * solve(chol_hessians[[i]], t(b)) + modes[i, ])
        dets[i] <- 1 / determinant.matrix(chol_hessians[[i]], logarithm = FALSE)$modulus
    }
    wGH <- as.matrix(expand.grid(lapply(seq_len(q), function (k, u) u$w, u = GH)))
    wGH <- 2^(q/2) * apply(wGH, 1, prod) * exp(rowSums(b * b))
    b2 <- lapply(b_new, function (b) if (q == 1) b * b else
        t(apply(b, 1, function (x) x %o% x)))
    ind_Z <- seq_len(ncol(Z_lis[[1]]))
    Ztb <- do.call('rbind', mapply(function (z, b) z %*% t(b[, ind_Z, drop = FALSE]), 
                                   Z_lis, b_new, SIMPLIFY = FALSE))
    Z_zitb <- if (!is.null(Z_zi_lis[[1]])) {
        do.call('rbind', mapply(function (z, b) z %*% t(b[, -ind_Z, drop = FALSE]), 
                                Z_zi_lis, b_new, SIMPLIFY = FALSE))  
    } 
    list(b = do.call('rbind', b_new), b2 = do.call('rbind', b2), Ztb = Ztb, Z_zitb = Z_zitb,
               wGH = wGH, dets = dets, post_modes = modes, 
         post_vars = lapply(aGH$post_hessian, solve))
}

chol_transf <- function (x) {
    if (any(is.na(x) | !is.finite(x)))
        stop("NA or infinite values in 'x'.\n")
    if (is.matrix(x)) {
        k <- nrow(x)
        U <- chol(x)
        U[cbind(1:k, 1:k)] <- log(U[cbind(1:k, 1:k)])
        U[upper.tri(U, TRUE)]
    } else {
        nx <- length(x)
        k <- round((-1 + sqrt(1 + 8 * nx))/2)
        mat <- matrix(0, k, k)
        mat[upper.tri(mat, TRUE)] <- x
        mat[cbind(1:k, 1:k)] <- exp(mat[cbind(1:k, 1:k)])
        res <- crossprod(mat)
        attr(res, "L") <- t(mat)[lower.tri(mat, TRUE)]
        res
    }
}

deriv_D <- function (D) {
    ncz <- nrow(D)
    ind <- which(lower.tri(D, TRUE), arr.ind = TRUE)
    dimnames(ind) <- NULL
    nind <- nrow(ind)
    svD <- solve(D)
    lapply(seq_len(nind), function (x, ind) {
        mat <- matrix(0, ncz, ncz)
        ii <- ind[x, , drop = FALSE]
        mat[ii[1], ii[2]] <- mat[ii[2], ii[1]] <- 1
        mat
    }, ind = ind[, 2:1, drop = FALSE])
}

jacobian2 <- function (L, ncz) {
    ind <- which(lower.tri(matrix(0, ncz, ncz), TRUE), arr.ind = TRUE)
    dimnames(ind) <- NULL
    nind <- nrow(ind)
    id <- 1:nind
    rind <- which(ind[, 1] == ind[, 2])
    lind <- vector("list", length(rind))
    for (i in seq_along(rind)) {
        tt <- matrix(0, ncz - i + 1, ncz - i + 1)
        tt[lower.tri(tt, TRUE)] <- seq(rind[i], nind)
        tt <- tt + t(tt)
        diag(tt) <- diag(tt)/2
        lind[[i]] <- tt
    }
    out <- matrix(0, nind, nind)
    for (g in 1:ncz) {
        gind <- id[g == ind[, 2]]
        vals <- L[gind]
        for (j in gind) {
            k <- which(j == gind)
            out[cbind(lind[[g]][k, ], j)] <- if (j %in% rind) vals[1] * vals else vals
        }
    }
    out[rind, ] <- 2 * out[rind, ]
    col.ind <- matrix(0, ncz, ncz)
    col.ind[lower.tri(col.ind, TRUE)] <- seq(1, length(L))
    col.ind <- t(col.ind)
    out[, col.ind[upper.tri(col.ind, TRUE)]]
}

fd <- function (x, f, ..., eps = .Machine$double.eps^0.25) {
    n <- length(x)
    res <- numeric(n)
    ex <- eps * (abs(x) + eps)
    f0 <- f(x, ...)
    for (i in seq_len(n)) {
        x1 <- x
        x1[i] <- x[i] + ex[i]
        diff.f <- c(f(x1, ...) - f0)
        diff.x <- x1[i] - x[i]
        res[i] <- diff.f / diff.x
    }
    res
}

fd_vec <- function (x, f, ..., eps = .Machine$double.eps^0.25) {
    n <- length(x)
    res <- matrix(0, n, n)
    ex <- pmax(abs(x), 1)
    f0 <- f(x, ...)
    for (i in 1:n) {
        x1 <- x
        x1[i] <- x[i] + eps * ex[i]
        diff.f <- c(f(x1, ...) - f0)
        diff.x <- x1[i] - x[i]
        res[, i] <- diff.f / diff.x
    }
    0.5 * (res + t(res))
}

cd <- function (x, f, ..., eps = 0.001) {
    n <- length(x)
    res <- numeric(n)
    ex <- pmax(abs(x), 1)
    for (i in seq_len(n)) {
        x1 <- x2 <- x
        x1[i] <- x[i] + eps * ex[i]
        x2[i] <- x[i] - eps * ex[i]
        diff.f <- c(f(x1, ...) - f(x2, ...))
        diff.x <- x1[i] - x2[i]
        res[i] <- diff.f / diff.x
    }
    res
}

cd_vec <- function (x, f, ..., eps = 0.001) {
    n <- length(x)
    res <- matrix(0, n, n)
    ex <- pmax(abs(x), 1)
    for (i in seq_len(n)) {
        x1 <- x2 <- x
        x1[i] <- x[i] + eps * ex[i]
        x2[i] <- x[i] - eps * ex[i]
        diff.f <- c(f(x1, ...) - f(x2, ...))
        diff.x <- x1[i] - x2[i]
        res[, i] <- diff.f / diff.x
    }
    0.5 * (res + t(res))
}

dmvnorm <- function (x, mu, Sigma, log = FALSE)  {
    if (!is.matrix(x))
        x <- rbind(x)
    p <- length(mu)
    if (p == 1) {
        dnorm(x, mu, sqrt(Sigma), log = log)
    } else {
        t1 <- length(mu) == length(Sigma)
        t2 <- all(abs(Sigma[lower.tri(Sigma)]) < sqrt(.Machine$double.eps))
        if (t1 || t2) {
            if (!t1)
                Sigma <- diag(Sigma)
            nx <- nrow(x)
            ff <- rowSums(dnorm(x, rep(mu, each = nx),
                                sd = rep(sqrt(Sigma), each = nx), log = TRUE))
            if (log) ff else exp(ff)
        } else {
            ed <- eigen(Sigma, symmetric = TRUE)
            ev <- ed$values
            evec <- ed$vectors
            if (!all(ev >= -1e-06 * abs(ev[1])))
                stop("'Sigma' is not positive definite")
            ss <- x - rep(mu, each = nrow(x))
            inv.Sigma <- evec %*% (t(evec)/ev)
            quad <- 0.5 * rowSums((ss %*% inv.Sigma) * ss)
            fact <- - 0.5 * (p * log(2 * pi) + sum(log(ev)))
            if (log) as.vector(fact - quad) else as.vector(exp(fact - quad))
        }
    }
}

unattr <- function (x) {
    if (is_mat <- is.matrix(x)) {
        d <- dim(x)
    }
    attributes(x) <- NULL
    if (is_mat) {
        dim(x) <- d
    }
    x
}

gauher <- function (n) {
    m <- trunc((n + 1) / 2)
    x <- w <- rep(-1, n)
    for (i in seq_len(m)) {
        z <- if (i == 1) {
            sqrt(2 * n + 1) - 1.85575 * (2 * n + 1)^(-0.16667)
        } else if (i == 2) {
            z - 1.14 * n^0.426/z
        } else if (i == 3) {
            1.86 * z - 0.86 * x[1]
        } else if (i == 4) {
            1.91 * z - 0.91 * x[2]
        } else {
            2 * z - x[i - 2]
        }
        for (its in seq_len(10)) {
            p1 <- 0.751125544464943
            p2 <- 0
            for (j in seq_len(n)) {
                p3 <- p2
                p2 <- p1
                p1 <- z * sqrt(2 / j) * p2 - sqrt((j - 1) / j) * p3
            }
            pp <- sqrt(2 * n) * p2
            z1 <- z
            z <- z1 - p1/pp
            if (abs(z - z1) <= 3e-14)
                break
        }
        x[i] <- z
        x[n + 1 - i] <- -z
        w[i] <- 2 / (pp * pp)
        w[n + 1 - i] <- w[i]
    }
    list(x = x, w = w)
}

nearPD <- function (M, eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08,
                    maxits = 100) {
    if (!(is.numeric(M) && is.matrix(M) && identical(M, t(M))))
        stop("Input matrix M must be square and symmetric.\n")
    inorm <- function(x) max(rowSums(abs(x)))
    n <- ncol(M)
    U <- matrix(0.0, n, n)
    X <- M
    iter <- 0
    converged <- FALSE
    while (iter < maxits && !converged) {
        Y <- X
        T <- Y - U
        e <- eigen(Y, symmetric = TRUE)
        Q <- e$vectors
        d <- e$values
        D <- if (length(d) > 1) diag(d) else as.matrix(d)
        p <- (d > eig.tol * d[1])
        QQ <- Q[, p, drop = FALSE]
        X <- QQ %*% D[p, p, drop = FALSE] %*% t(QQ)
        U <- X - T
        X <- (X + t(X)) / 2
        conv <- inorm(Y - X)/inorm(Y)
        iter <- iter + 1
        converged <- conv <= conv.tol
    }
    X <- (X + t(X)) / 2
    e <- eigen(X, symmetric = TRUE)
    d <- e$values
    Eps <- posd.tol * abs(d[1L])
    if (d[n] < Eps) {
        d[d < Eps] <- Eps
        Q <- e$vectors
        o.diag <- diag(X)
        X <- Q %*% (d * t(Q))
        D <- sqrt(pmax(Eps, o.diag) / diag(X))
        X[] <- D * X * rep(D, each = n)
    }
    (X + t(X)) / 2
}

getRE_Formula <- function (form) {
    if (!(inherits(form, "formula"))) {
        stop("formula(object) must return a formula")
    }
    form <- form[[length(form)]]
    if (length(form) == 3 && (form[[1]] == as.name("|") || form[[1]] == as.name("||"))) {
        form <- form[[2]]
    }
    eval(substitute(~form))
}

getID_Formula <- function (form) {
    if (is.list(form)) {
        nams <- names(form)
        as.formula(paste0("~", nams[1L], "/", nams[2L]))
    } else {
        form <- form[[length(form)]]
        asOneSidedFormula(form[[3]])
    }
}

printCall <- function (call) {
    d <- deparse(call)
    if (length(d) <= 3) {
        paste(d, sep = "\n", collapse = "\n")
    } else {
        d <- d[1:3]
        d[3] <- paste0(d[3], "...")
        paste(d, sep = "\n", collapse = "\n")
    }
}

dgt <- function (x, mu = 0, sigma = 1, df = stop("no df argument."), log = FALSE) {
    if (log) {
        dt(x = (x - mu) / sigma, df = df, log = TRUE) - log(sigma)
    } else {
        dt(x = (x - mu) / sigma, df = df) / sigma
    }
}

dmvt <- function (x, mu, Sigma = NULL, invSigma = NULL, df, log = TRUE, prop = TRUE) {
    if (!is.numeric(x)) 
        stop("'x' must be a numeric matrix or vector")
    if (!is.matrix(x)) 
        x <- rbind(x)
    p <- length(mu)
    if (is.null(Sigma) && is.null(invSigma))
        stop("'Sigma' or 'invSigma' must be given.")
    if (!is.null(Sigma)) {
        if (is.list(Sigma)) {
            ev <- Sigma$values
            evec <- Sigma$vectors
        } else {
            ed <- eigen(Sigma, symmetric = TRUE)
            ev <- ed$values
            evec <- ed$vectors
        }
        if (!all(ev >= -1e-06 * abs(ev[1]))) 
            stop("'Sigma' is not positive definite")
        invSigma <- evec %*% (t(evec)/ev)
        if (!prop)
            logdetSigma <- sum(log(ev))
    } else {
        if (!prop)
            logdetSigma <- c(-determinant(invSigma)$modulus)
    }
    ss <- x - rep(mu, each = nrow(x))
    quad <- rowSums((ss %*% invSigma) * ss)/df
    if (!prop)
        fact <- lgamma((df + p)/2) - lgamma(df/2) - 
        0.5 * (p * (log(pi) + log(df)) + logdetSigma)
    if (log) {
        if (!prop) as.vector(fact - 0.5 * (df + p) * log(1 + quad)) else 
            as.vector(- 0.5 * (df + p) * log(1 + quad))
    } else {
        if (!prop) as.vector(exp(fact) * ((1 + quad)^(-(df + p)/2))) else 
            as.vector(((1 + quad)^(-(df + p)/2)))
    }
}

rmvt <- function (n, mu, Sigma, df) {
    p <- length(mu)
    if (is.list(Sigma)) {
        ev <- Sigma$values
        evec <- Sigma$vectors
    } else {
        ed <- eigen(Sigma, symmetric = TRUE)
        ev <- ed$values
        evec <- ed$vectors
    }
    X <- drop(mu) + tcrossprod(evec * rep(sqrt(pmax(ev, 0)), each = p), 
                               matrix(rnorm(n * p), n)) / rep(sqrt(rchisq(n, df)/df), each = p)
    if (n == 1L) drop(X) else t.default(X)
}

register_s3_method <- function (pkg, generic, class) {
    fun <- get(paste0(generic, ".", class), envir = parent.frame())
    if (isNamespaceLoaded(pkg))
        registerS3method(generic, class, fun, envir = asNamespace(pkg))
    # Also ensure registration is done if pkg is loaded later:
    setHook(
        packageEvent(pkg, "onLoad"),
        function (...)
            registerS3method(generic, class, fun, envir = asNamespace(pkg))
    )
}

.onLoad <- function (libname, pkgname) {
    if (requireNamespace("emmeans", quietly = TRUE)) {
        register_s3_method("emmeans", "recover_data", "MixMod")
        register_s3_method("emmeans", "emm_basis", "MixMod")
    }
    if (requireNamespace("effects", quietly = TRUE)) {
        register_s3_method("effects", "Effect", "MixMod")
    }
}

constructor_form_random <- function (formula, data) {
    groups <- all.vars(getID_Formula(formula))
    ngroups <- length(groups)
    formula <- if (!is.list(formula)) {
        form_random <- vector("list", ngroups)
        names(form_random) <- groups
        form_random[] <- lapply(form_random, function (x) getRE_Formula(formula))
    } else formula
    if (ngroups > 1) {
        nesting <- function (form, group_name) {
            terms_form <- attr(terms(form), "term.labels")
            if (length(terms_form)) {
                interaction_terms <- paste0(group_name, ":", terms_form, collapse = " + ")
                as.formula(paste0("~ 0 + ", group_name, " + ", interaction_terms))
            } else {
                as.formula(paste0("~ 0 + ", group_name))
            }
        }
        formula[-1] <- mapply(nesting, formula[-1], groups[-1], SIMPLIFY = FALSE)
    }
    formula
}

constructor_Z <- function (termsZ_i, mfZ_i, id) {
    n <- length(unique(id))
    Zmats <- vector("list", n)
    for (i in seq_len(n)) {
        mf <- model.frame(termsZ_i, mfZ_i[id == i, , drop = FALSE],
                          drop.unused.levels = TRUE)
        mm <- model.matrix(termsZ_i, mf)
        assign <- attr(mm, "assign")
        Zmats[[i]] <- mm[, c(t(sapply(unique(assign), function (x) which(assign == x)))), 
                         drop = FALSE]
    }
    do.call("rbind", Zmats)
}

cr_setup <- function (y, direction = c("forward", "backward")) {
    direction <- match.arg(direction)
    yname <- as.character(substitute("y"))
    if (!is.factor(y)) {
        y <- factor(y)
    }
    ylevels <- levels(y)
    ncoefs <- length(ylevels) - 1
    if (ncoefs < 2) {
        stop("it seems that variable ", yname, " has two levels; use a mixed effects ",
             "logistic regression instead.\n")
    }
    y <- as.numeric(unclass(y) - 1)
    if (direction == "forward") {
        reps <- ifelse(is.na(y), 1, ifelse(y < ncoefs - 1, y + 1, ncoefs))
        subs <- rep(seq_along(y), reps)
        cuts <- vector("list", ncoefs + 2)
        cuts[[1]] <- NA
        for (j in seq(0, ncoefs)) {
            cuts[[j + 2]] <- seq(0, if (j < ncoefs - 1) j else ncoefs - 1)
        } 
        cuts <- unlist(cuts[ifelse(is.na(y), 1, y + 2)], use.names = FALSE)
        labels <- c("all", paste0(yname, ">=", ylevels[2:ncoefs]))
        y <- rep(y, reps)
        Y <- as.numeric(y == cuts)
    } else {
        reps <- ifelse(is.na(y), 1, ifelse(y > ncoefs - 3, ncoefs - (y - 1), ncoefs))
        subs <- rep(seq_along(y), reps)
        cuts <- vector("list", ncoefs + 2)
        cuts[[ncoefs + 2]] <- NA
        for (j in seq(ncoefs, 0)) {
            cuts[[j + 1]] <- seq(0, ncoefs - if (j > ncoefs - 3) j else 1)
        } 
        cuts <- unlist(cuts[ifelse(is.na(y), 1, y + 1)], use.names = FALSE)
        labels <- c("all", paste0(yname, "<=", ylevels[ncoefs:2]))
        y <- rep(y, reps)
        Y <- as.numeric(y == (ncoefs - cuts))
    }
    cohort <- factor(cuts, levels = seq(0, ncoefs - 1), labels = labels)
    list(y = Y, cohort = cohort, subs = subs, reps = reps)
}

cr_marg_probs <- function (eta, direction = c("forward", "backward")) {
    direction <- match.arg(direction)
    ncoefs <- ncol(eta)
    if (direction == "forward") {
        cumsum_1_minus_p <- t(apply(plogis(eta[, -ncoefs], log.p = TRUE, 
                                           lower.tail = FALSE), 1, cumsum))
        probs <- exp(plogis(eta, log.p = TRUE) + cbind(0, cumsum_1_minus_p))
        cbind(probs, 1 - rowSums(probs))
    } else {
        cumsum_1_minus_p <- t(apply(plogis(eta[, seq(ncoefs, 2)], log.p = TRUE, 
                                           lower.tail = FALSE), 1, cumsum))
        probs <- exp(plogis(eta, log.p = TRUE) + 
                         cbind(cumsum_1_minus_p[, seq(ncoefs - 1, 1)], 0))
        cbind(1 - rowSums(probs), probs)
    }
}