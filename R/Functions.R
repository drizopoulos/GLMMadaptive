find_modes <- function (b, y_lis, N_lis, X_lis, Z_lis, offset_lis, betas, invD, phis, 
                        canonical, user_defined, Zty_lis, log_dens, mu_fun, var_fun, 
                        mu.eta_fun, score_eta_fun, score_phis_fun) {
    log_post_b <- function (b_i, y_i, N_i, X_i, Z_i, offset_i, betas, invD, phis, canonical,
                            user_defined, Zty_i, log_dens, mu_fun, var_fun, mu.eta_fun,
                            score_eta_fun, score_phis_fun) {
        eta_y <- as.vector(X_i %*% betas + Z_i %*% b_i)
        if (!is.null(offset_i))
            eta_y <- eta_y + offset_i
        - sum(log_dens(y_i, eta_y, mu_fun, phis), na.rm = TRUE) +
                c(0.5 * crossprod(b_i, invD) %*% b_i)
    }
    score_log_post_b <- function (b_i, y_i, N_i, X_i, Z_i, offset_i, betas, invD, phis,
                                  canonical, user_defined, Zty_i, log_dens, mu_fun,
                                  var_fun, mu.eta_fun, score_eta_fun, score_phis_fun) {
        eta_y <- c(X_i %*% betas + Z_i %*% b_i)
        if (!is.null(offset_i))
            eta_y <- eta_y + offset_i
        mu_y <- mu_fun(eta_y)
        log_dens_part <- if (user_defined) {
            if (!is.null(score_eta_fun)) {
                - crossprod(Z_i, score_eta_fun(y_i, mu_y, phis))
            } else {
                l1 <- log_dens(y_i, eta_y + 1e-04, mu_fun, phis)
                l2 <- log_dens(y_i, eta_y - 1e-04, mu_fun, phis)
                (l1 - l2) / (2 * 1e-04)
                - crossprod(Z_i, (l1 - l2) / (2 * 1e-04))
            }
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
        b_i <- b[i, , drop = FALSE]
        opt_i <- optim(par = b_i, fn = log_post_b, gr = score_log_post_b, method = "BFGS",
                       y_i = y_i, N_i = N_i, X_i = X_i, Z_i = Z_i, offset_i = offset_i,
                       betas = betas, invD = invD, phis = phis, canonical = canonical,
                       user_defined = user_defined, Zty_i = Zty_i, log_dens = log_dens,
                       mu_fun = mu_fun, var_fun = var_fun, mu.eta_fun = mu.eta_fun,
                       score_eta_fun = score_eta_fun, score_phis_fun = score_phis_fun)
        post_modes[i, ] <- opt_i$par
        post_hessians[[i]] <- cd_vec(post_modes[i, ], score_log_post_b,
                                     y_i = y_i, N_i = N_i, X_i = X_i, Z_i = Z_i, 
                                     offset_i = offset_i, betas = betas, invD = invD, 
                                     phis = phis, canonical = canonical, 
                                     user_defined = user_defined,
                                     Zty_i = Zty_i, log_dens = log_dens, mu_fun = mu_fun,
                                     var_fun = var_fun, mu.eta_fun = mu.eta_fun,
                                     score_eta_fun = score_eta_fun,
                                     score_phis_fun = score_phis_fun)
    }
    list(post_modes = post_modes, post_hessians = post_hessians)
}

GHfun <- function (b, y_lis, N_lis, X_lis, Z_lis, offset_lis, betas, inv_D, phis, k, q,
                   canonical, user_defined, Zty_lis, log_dens, mu_fun, var_fun, mu.eta_fun,
                   score_eta_fun, score_phis_fun) {
    GH <- gauher(k)
    aGH <- find_modes(b, y_lis, N_lis, X_lis, Z_lis, offset_lis, betas, inv_D, phis, 
                      canonical, user_defined, Zty_lis, log_dens, mu_fun, var_fun, 
                      mu.eta_fun, score_eta_fun, score_phis_fun)
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
    Ztb <- do.call('rbind', mapply(function (z, b) z %*% t(b), Z_lis, b_new, SIMPLIFY = FALSE))
    list(b = do.call('rbind', b_new), b2 = do.call('rbind', b2), Ztb = Ztb,
               wGH = wGH, dets = dets, post_modes = modes)
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
    form <- form[[length(form)]]
    asOneSidedFormula(form[[3]])
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

