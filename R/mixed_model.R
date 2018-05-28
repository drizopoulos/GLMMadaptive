mixed_model <- function (fixed, random, data, family = NULL, na.action = na.exclude,
                         n_phis = NULL, initial_values = NULL, control = list(), ...) {
    call <- match.call()
    # set family
    if (is.null(family)) {
        stop("argument 'family' needs to be specified; check the online help file.\n")
    } else {
        if (is.character(family))
            family <- get(family, mode = "function", envir = parent.frame())
        if (is.function(family))
            family <- family()
        if (is.null(family$family)) {
            print(family)
            stop("'family' not recognized.\n")
        }
        if (family$family == "gaussian")
            stop("use function lme() from package 'nlme' or function lmer() from ",
                 "package 'lme4'.\n")
    }
    known_families <- c("binomial", "poisson", "negative binomial")
    # extract response vector, design matrices, offset
    mfX <- model.frame(terms(fixed, data = data), data = data, na.action = na.action)
    na_exclude <- attr(mfX, "na.action")
    termsX <- terms(mfX)
    y <- model.response(mfX)
    if (is.factor(y)) {
        if (family$family == "binomial")
            y <- as.numeric(y != levels(y)[1L])
        else
            stop("the response variable should not be a factor.\n")
    }
    X <- model.matrix(termsX, mfX)
    offset <- model.offset(mfX)
    form_random <- getRE_Formula(random)
    mfZ <- model.frame(terms(form_random, data = data), data = data)
    if (!is.null(na_exclude))
        mfZ <- mfZ[-na_exclude, ]
    termsZ <- terms(mfZ)
    Z <- model.matrix(termsZ, mfZ)
    id_orig <- model.frame(terms(getID_Formula(random)), data)[[1L]]
    id <- match(id_orig, unique(id_orig))
    ###########################
    # control settings
    con <- list(iter_EM = 30, iter_qN_outer = 15, iter_qN = 10, iter_qN_incr = 10,
                optim_method = "BFGS", parscale_betas = 0.1, parscale_D = 0.01,
                parscale_phis = 0.01, tol1 = 1e-03, tol2 = 1e-04, tol3 = 1e-07,
                numeric_deriv = "fd", nAGQ = 11, update_GH_every = 10, verbose = FALSE)
    control <- c(control, list(...))
    namC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% namC]) > 0)
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    ###########################
    # initial values
    inits <- if (family$family %in% known_families || inherits(initial_values, 'family')) {
        betas <- if (family$family %in% known_families) {
            glm.fit(X, y, family = family)$coefficients
        } else {
            glm.fit(X, y, family = initial_values)$coefficients
        }
        list(betas = betas * sqrt(1.346), D = diag(ncol(Z)))
    } else {
        list(betas = rep(0, ncol(X)), D = diag(ncol(Z)))
    }
    if (!is.null(initial_values) && is.list(initial_values) &&
        !inherits(initial_values, 'family')) {
        lngths <- lapply(inits[(nams.initial_values <- names(initial_values))], length)
        if (!isTRUE(all.equal(lngths, lapply(initial_values, length)))) {
            warning("'initial_values' is not a list with elements of appropriate ",
                    "length; default initial_values are used instead.\n")
        } else {
            inits[nams.initial_values] <- initial_values
        }
    }
    ##########################
    # Functions
    Funs <- list(
        mu_fun = family$linkinv,
        var_fun = family$variance,
        mu.eta_fun = family$mu.eta
    )
    if (family$family %in% known_families && is.null(family$log_den)) {
        Funs$log_dens <- switch(family$family,
               'binomial' = binomial_log_dens,
               'poisson' = poisson_log_dens,
               'negative binomial' = negative.binomial_log_dens)
    } else if (!family$family %in% known_families && !is.null(family$log_den)) {
        Funs$log_dens <- family$log_dens
    } else {
        stop("'log_dens' component of the 'family' argument is NULL with no default.\n")
    }
    if (!is.function(Funs$log_dens)) {
        stop("'log_dens' component of the 'family' argument must be a function.\n")
    }
    if (!is.function(Funs$mu_fun)) {
        stop("'linkinv' component of the 'family' argument must be a function.\n")
    }
    if (!is.function(Funs$mu_fun)) {
        stop("'linkinv' component of the 'family' argument must be a function.\n")
    }
    if (!is.null(family$score_eta_fun) && is.function(family$score_eta_fun)) {
        Funs$score_eta_fun <- score_eta_fun
    }
    if (!is.null(family$score_phis_fun) && is.function(family$score_phis_fun)) {
        Funs$score_phis_fun <- score_phis_fun
    }
    has_phis <- inherits(try(Funs$log_dens(y, 0, Funs$mu_fun, phis = NULL), TRUE),
                         "try-error")
    if (has_phis) {
        if (is.null(n_phis)) {
            stop("argument 'n_phis' needs to be specified.\n")
        }
        inits$phis <- rep(0.0, n_phis)
    }
    ###############
    # Fit the model
    out <- mixed_fit(y, X, Z, id, offset, family, inits, Funs, con)
    # fix names
    names(out$coefficients) <- colnames(X)
    dimnames(out$D) <- list(colnames(Z), colnames(Z))
    out$id <- id_orig
    dimnames(out$post_modes) <- list(unique(id_orig), colnames(Z))
    out$Terms <- list(termsX = termsX, termsZ = termsZ)
    out$model_frames <- list(mfX = mfX, mfZ = mfZ)
    out$control <- con
    out$Funs <- Funs
    out$family <- family
    out$call <- call
    class(out) <- "MixMod"
    out
}
