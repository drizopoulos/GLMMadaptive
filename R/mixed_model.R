mixed_model <- function (fixed, random, data, family, na.action = na.exclude,
                         zi_fixed = NULL, zi_random = NULL, 
                         penalized = FALSE, n_phis = NULL, initial_values = NULL, 
                         control = list(), ...) {
    call <- match.call()
    # set family
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        stop("'family' not recognized.\n")
    }
    if (family$family == "gaussian")
        stop("use function lme() from package 'nlme' or function lmer() from ",
             "package 'lme4'.\n")
    if (length(grep("Negative Binomial", family$family))) {
        stop("Because the namespace of the MASS package seems also to be loaded\n",
             "  use 'family = GLMMadaptive::negative.binomial(xx)' with 'xx' ", 
             "denoting a value for the\n  'theta' parameter of the family.")
    }
    if (family$family %in% c("zero-inflated poisson", "zero-inflated negative binomial",
                             "hurdle poisson", "hurdle negative binomial", "hurdle beta") && 
        is.null(zi_fixed)) {
        stop("you have defined a family with an extra zero-part;\nat least argument ",
             "'zi_fixed' needs to be defined, and potentially also argument 'zi_random'.")
    }
    if (family$family %in% c("binomial", "poisson", "negative binomial", "beta") && 
        !is.null(zi_fixed)) {
        stop("\nyou have defined a family object *without* an extra zero-part but\n ",
             "you have also specified the 'zi_fixed' argument; use instead a family\n ", 
             "object with an extra zero-part.")
    }
    known_families <- c("binomial", "poisson", "negative binomial")
    data <- as.data.frame(data) # in case 'data' is a tibble
    groups <- unique(c(all.vars(getID_Formula(random)), 
                       if (!is.null(zi_random)) all.vars(getID_Formula(zi_random))))
    data[groups] <- lapply(data[groups], function (x) if (!is.factor(x)) factor(x) else x)
    # construct model frames
    # fixed effects
    mfX <- model.frame(terms(fixed, data = data), data = data, na.action = na.pass)
    termsX <- terms(mfX)
    # random effects
    form_random <- constructor_form_random(random, data)
    mfZ <- lapply(form_random, function (form, data) 
        model.frame(terms(form, data = data), data = data, na.action = na.pass), 
        data = data)
    termsZ <- lapply(mfZ, terms)
    # fixed effects ZI part
    if (!is.null(zi_fixed)) {
        mfX_zi <- model.frame(terms(zi_fixed, data = data), data = data, 
                              na.action = na.pass)
        termsX_zi <- terms(mfX_zi)
    } else {
        mfX_zi <- termsX_zi <- NULL
    }
    # random effects ZI part
    if (!is.null(zi_random)) {
        form_zi_random <- constructor_form_random(zi_random, data)
        mfZ_zi <- lapply(form_zi_random, function (form, data) 
            model.frame(terms(form, data = data), data = data, na.action = na.pass),
            data = data)
        termsZ_zi <- lapply(mfZ_zi, terms)
    } else {
        mfZ_zi <- termsZ_zi <- NULL
    }
    # delete missing data if na.action is na.fail or na.omit
    chr_na_action <- as.character(body(na.action))[2]
    if (chr_na_action == "na.exclude" || chr_na_action == "na.omit") {
        complete_cases <- cbind(complete.cases(mfX), sapply(mfZ, complete.cases))
        if (!is.null(zi_fixed))
            complete_cases <- cbind(complete_cases, complete.cases(mfX_zi))
        if (!is.null(zi_random))
            complete_cases <- cbind(complete_cases, sapply(mfZ_zi, complete.cases))
        keep <- apply(complete_cases, 1, all)
        mfX <- mfX[keep, , drop = FALSE]
        mfZ[] <- lapply(mfZ, function (mf) mf[keep, , drop = FALSE])
        if (!is.null(zi_fixed))
            mfX_zi <- mfX_zi[keep, , drop = FALSE]
        if (!is.null(zi_random))
            mfZ_zi[] <- lapply(mfZ_zi, function (mf) mf[keep, , drop = FALSE])
    }
    # id variable
    id_nam <- all.vars(getID_Formula(random))
    id_orig <- model.frame(terms(getID_Formula(random), data = data), 
                           data = data, na.action = na.pass)
    if (chr_na_action == "na.exclude" || chr_na_action == "na.omit")
        id_orig <- id_orig[keep, , drop = FALSE]
    id <- match(id_orig[[1L]], unique(id_orig[[1L]]))
    # model response
    y <- model.response(mfX)
    if (is.factor(y)) {
        if (family$family == "binomial")
            y <- as.numeric(y != levels(y)[1L])
        else
            stop("the response variable should not be a factor.\n")
    }
    # construct model matrices
    # fixed effects
    offset <- model.offset(mfX)
    X <- model.matrix(termsX, mfX)
    # random effects
    Z <- mapply(constructor_Z, termsZ, mfZ, MoreArgs = list(id = id), SIMPLIFY = FALSE)
    Z <- do.call("cbind", Z)
    # fixed effects ZI part
    if (!is.null(zi_fixed)) {
        offset_zi <- model.offset(mfX_zi)
        X_zi <- model.matrix(termsX_zi, mfX_zi)
    } else {
        offset_zi <- X_zi <- NULL
    }
    # random effects ZI part
    if (!is.null(zi_random)) {
        id_nam_zi <- all.vars(getID_Formula(zi_random))
        if (id_nam_zi[1L] != id_nam[1L]) {
            stop("the 'random' and 'zi_random' formulas for the random effects ",
                 "should have the same grouping variable.")
        }
        Z_zi <- mapply(constructor_Z, termsZ_zi, mfZ_zi, MoreArgs = list (id = id), 
                       SIMPLIFY = FALSE)
        Z_zi <- do.call("cbind", Z_zi)
    } else {
        Z_zi <- NULL
    }
    nRE <- ncol(Z) + if (!is.null(Z_zi)) ncol(Z_zi) else 0
    ###########################
    # control settings
    con <- list(iter_EM = 30, iter_qN_outer = 15, iter_qN = 10, iter_qN_incr = 10,
                optimizer = "optim", optim_method = "BFGS", parscale_betas = 0.1, 
                parscale_D = 0.01, parscale_phis = 0.01, parscale_gammas = 0.01, 
                tol1 = 1e-04, tol2 = 1e-05, tol3 = 1e-08, numeric_deriv = "fd", 
                nAGQ = if (nRE < 3) 11 else 7, update_GH_every = 10, max_coef_value = 12, 
                max_phis_value = exp(10), verbose = FALSE)
    control <- c(control, list(...))
    namC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% namC]) > 0)
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    ###########################
    # initial values
    diag_D <- (random[[length(random)]])[[1]] == as.name("||")
    inits <- if (family$family %in% known_families || (is.list(initial_values) &&
                                                       inherits(initial_values$betas, 'family'))) {
        betas <- if (family$family %in% known_families) {
            if (family$family == "negative binomial")
                glm.fit(X, y, family = poisson(), offset = offset)$coefficients
            else 
                glm.fit(X, y, family = family, offset = offset)$coefficients
        } else {
            glm.fit(X, y, family = initial_values$betas, offset = offset)$coefficients
        }
        list(betas = betas * sqrt(1.346), D = if (diag_D) rep(1, nRE) else diag(nRE))
    } else {
        list(betas = rep(0, ncol(X)), D = if (diag_D) rep(1, nRE) else diag(nRE))
    }
    if (!is.null(zi_fixed)) {
        inits <- c(inits, 
                   list(gammas = glm.fit(X_zi, 
                                         as.numeric(if (NCOL(y) == 2) y[, 1] == 0 else y == 0), 
                                         family = binomial(),
                                         offset = offset_zi)$coefficients))
        if (family$family %in% c("zero-inflated poisson", "zero-inflated negative binomial",
                                 "hurdle poisson", "hurdle negative binomial", 
                                 "Conway Maxwell Poisson"))
            inits$betas <- glm.fit(X, y, family = poisson())$coefficients
    }
    ##########################
    # penalized
    penalized <- if (is.logical(penalized) && !penalized) {
        list(penalized = penalized)
    } else if (is.logical(penalized) && penalized) {
        list(penalized = penalized, pen_mu = 0, pen_sigma = 1, pen_df = 3)
    } else if (is.list(penalized)) {
        if (!all(names(penalized) %in% c("pen_mu", "pen_sigma", "pen_df")))
            stop("when argument 'penalized' is a list it needs to have the components ",
                 "'pen_mu', 'pen_sigma' and 'pen_df'.\n")
        c(list(penalized = TRUE), penalized)
    } else {
        stop("argument 'penalized' must be a logical or a list.\n")
    }
    ##########################
    # Functions
    Funs <- list(
        mu_fun = family$linkinv,
        var_fun = family$variance,
        mu.eta_fun = family$mu.eta
    )
    if (family$family %in% known_families && is.null(family$log_dens)) {
        Funs$log_dens <- switch(family$family,
                                'binomial' = binomial_log_dens,
                                'poisson' = poisson_log_dens,
                                'negative binomial' = negative.binomial_log_dens)
    } else if (family$family %in% known_families && !is.null(family$log_dens)) {
        Funs$log_dens <- family$log_dens
    } else if (!family$family %in% known_families && !is.null(family$log_dens)) {
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
        Funs$score_eta_fun <- family$score_eta_fun
    }
    if (!is.null(family$score_eta_zi_fun) && is.function(family$score_eta_zi_fun)) {
        Funs$score_eta_zi_fun <- family$score_eta_zi_fun
    }
    if (!is.null(family$score_phis_fun) && is.function(family$score_phis_fun)) {
        Funs$score_phis_fun <- family$score_phis_fun
    }
    has_phis <- inherits(try(Funs$log_dens(y, rep(0, length(y)), Funs$mu_fun, 
                                           phis = NULL, rep(0, length(y))), TRUE),
                         "try-error")
    if (has_phis) {
        if (family$family %in% c("negative binomial", "zero-inflated negative binomial",
                                 "hurdle negative binomial", "hurdle log-normal",
                                 "beta", "hurdle beta", "Conway Maxwell Poisson")) {
            n_phis <- 1
        } else if (is.null(n_phis)) {
            stop("argument 'n_phis' needs to be specified.\n")
        }
        inits$phis <- rep(0.0, n_phis)
    }
    if (!is.null(initial_values) && is.list(initial_values) &&
        !inherits(initial_values$betas, 'family')) {
        lngths <- lapply(inits[(nams.initial_values <- names(initial_values))], length)
        if (!isTRUE(all.equal(lngths, lapply(initial_values, length)))) {
            warning("'initial_values' is not a list with elements of appropriate ",
                    "length; default initial_values are used instead.\n")
        } else {
            inits[nams.initial_values] <- initial_values
        }
    }
    ###############
    # Fit the model
    out <- mixed_fit(y, X, Z, X_zi, Z_zi, id, offset, offset_zi, family, inits, Funs, 
                     con, penalized)
    # check whether Hessian is positive definite at convergence
    H <- out$Hessian
    if (any(is.na(H) | !is.finite(H))) {
        warning("infinite or missing values in Hessian at convergence.\n")
    } else {
        ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
        if (!all(ev >= -1e-06 * abs(ev[1L]))) 
            warning("Hessian matrix at convergence is not positive definite; ", 
                    "unstable solution.\n")
    }
    # fix names
    names(out$coefficients) <- colnames(X)
    RE_nams <- c(colnames(Z), if (!is.null(Z_zi)) paste0("zi_", colnames(Z_zi)))
    dimnames(out$D) <- list(RE_nams, RE_nams)
    if (!is.null(out$phis))
        names(out$phis) <- paste0("phi_", seq_along(out$phis))
    if (!is.null(out$gammas))
        names(out$gammas) <- colnames(X_zi)
    all_nams <- if (diag_D) {
        nams_D <- paste0("D_", seq_len(nRE), seq_len(nRE))
        c(names(out$coefficients), nams_D, names(out$phis), 
          if (!is.null(out$gammas)) paste0("zi_", names(out$gammas)))
    } else {
        nams_D <- paste0("D_", apply(which(upper.tri(out$D, TRUE), arr.ind = TRUE), 1, 
                                     paste0, collapse = ""))
        c(names(out$coefficients), nams_D, names(out$phis), 
          if (!is.null(out$gammas)) paste0("zi_", names(out$gammas)))
    }
    dimnames(out$Hessian) <- list(all_nams, all_nams)
    out$id <- id_orig
    out$id_name <- id_nam 
    out$offset <- offset
    dimnames(out$post_modes) <- list(unique(id_orig[[1L]]), RE_nams)
    names(out$post_vars) <- unique(id_orig)
    out$post_vars[] <- lapply(out$post_vars, function (v) {
        dimnames(v) <- list(RE_nams, RE_nams)
        v
    })
    out$Terms <- list(termsX = termsX, termsZ = termsZ, termsX_zi = termsX_zi, 
                      termsZ_zi = termsZ_zi)
    out$model_frames <- list(mfX = mfX, mfZ = mfZ, mfX_zi = mfX_zi, mfZ_zi = mfZ_zi)
    out$control <- con
    out$Funs <- Funs
    out$family <- family
    out$penalized <- penalized
    out$na.action <- na.action
    out$contrasts <- attr(X, "contrasts")
    out$call <- call
    class(out) <- "MixMod"
    out
}
