\name{mixed_model}
\alias{mixed_model}

\title{
Generalized Linear Mixed Effects Models
}

\description{
Fits generalized linear mixed effects models under maximum likelihood using adaptive 
Gaussian quadrature.
}

\usage{
mixed_model(fixed, random, data, family, weights = NULL,
  na.action = na.exclude, zi_fixed = NULL, zi_random = NULL, 
  penalized = FALSE, n_phis = NULL, initial_values = NULL, 
  control = list(), \dots)
}

\arguments{
\item{fixed}{a formula for the fixed-effects part of the model, including the outcome.}
\item{random}{a formula for the random-effects part of the model. This should only contain
    the right-hand side part, e.g., \code{~ time | id}, where \code{time} is a variable,
    and \code{id} the grouping factor. When the symbol \code{||} is used in the definition
    of this argument (instead of \code{|}), then the covariance matrix of the random 
    effects is assumed to be diagonal.}
\item{data}{a data.frame containing the variables required in \code{fixed} and 
    \code{random}.}
\item{family}{a \code{\link{family}} object specifying the type of the repeatedly measured response 
    variable, e.g., \code{binomial()} or \code{poisson()}. The function also allows for 
    user-defined family objects, but with specific extra components; see the example
    is \code{\link{negative.binomial}} for more details. Contrary to the standard 
    practice in model fitting R functions with a \code{family} argument 
    (e.g., \code{\link{glm}}) in which the default
    family is \code{gaussian()}, in \code{mixed_model()} no default is provided. If the
    users wish to fit a mixed model for a Gaussian outcome, this could be done with 
    function \code{lme()} from the \bold{nlme} package or function \code{lmer()} from the
    \bold{lme4} package.}
\item{weights}{a numeric vector of weights. These are simple multipliers on the 
  log-likelihood contributions of each group/cluster, i.e., we presume that there are multiple 
  replicates of each group/cluster denoted by the weights. The length of `weights` need to
  be equal to the number of independent groups/clusters in the data.}
\item{na.action}{what to do with missing values in \code{data}.}
\item{zi_fixed, zi_random}{formulas for the fixed and random effects of the zero inflated 
    part.}
\item{penalized}{logical or a list. If logical and equal to \code{FALSE}, then no penalty
    is used. If logical and equal to \code{TRUE}, for the fixed effects a Student's-t 
    penalty/prior with mean 0, scale equal to 1 and 3 degrees of freedom is used. If a
    list, then it is expected to have the components \code{pen_mu}, \code{pen_sigma} and
    \code{pen_df}, denoting the mean, scale and degrees of freedom of the Student's-t 
    penalty/prior for the fixed effects.}
\item{n_phis}{a numeric scalar; in case the family corresponds to a distribution that has extra 
    (dispersion/shape) parameters, you need to specify how many extra parameters are needed.}
\item{initial_values}{a list of initial values. This can have up to three components, 
    namely,
    \describe{
        \item{betas}{a numeric vector of fixed effects. This can also be 
        \code{\link{family}} object. In this case initial values for the fixed effects 
        will be calculated by using \code{\link{glm}} to the data ignoring the 
        correlations in the repeated measurements. For example, for a negative binomial
        response outcome, we could set \code{betas = poisson()}.}
        \item{D}{a numeric matrix denoting the covariance matrix of the random effects.}
        \item{phis}{a numeric vector for the extra (dispersion/shape) parameters.}
    }
    }
\item{control}{a list with the following components:
    \describe{
        \item{iter_EM}{numeric scalar denoting the number of EM iterations; default is 30.}
        \item{iter_qN_outer}{numeric scalar denoting the number of outer iterations during 
            the quasi-Newton phase; default is 15. In each outer iteration the locations 
            of the quadrature points are updated.}
        \item{iter_qN}{numeric scalar denoting the starting number of iterations for the 
            quasi-Newton; default is 10.}
        \item{iter_qN_incr}{numeric scalar denoting the increment in \code{iter_qN} for
            each outer iteration; default is 10.}
        \item{optimizer}{character string denoting the optimizer to be used; available 
            options are \code{"optim"} (default), \code{"nlminb"} and 
            \code{"optimParallel"}, the last option implemented in the 
            \bold{optimParallel} package.}
        \item{optim_method}{character string denoting the type of \code{\link{optim}} 
            algorithm to be used when \code{optimizer = "optim"}; default is the BFGS 
            algorithm.}
        \item{parscale_betas}{the control argument \code{parscale} of \code{\link{optim}} for
            the fixed-effects; default is 0.1.}
        \item{parscale_D}{the control argument \code{parscale} of \code{\link{optim}} for
            the unique element of the covariance matrix of the random effects; default is 
            0.01.}
        \item{parscale_phis}{the control argument \code{parscale} of \code{\link{optim}} for
            the extra (dispersion/shape) parameters; default is 0.01.}
        \item{tol1, tol2, tol3}{numeric scalars controlling tolerances for declaring 
            convergence; \code{tol1} and \code{tol2} are for checking convergence in 
            successive parameter values; \code{tol3} is similar to \code{reltop} of 
            \code{\link{optim}}; default values are \code{1e-03}, \code{1e-04}, and 
            \code{1e-08}, respectively.}
        \item{numeric_deriv}{character string denoting the type of numerical derivatives
            to be used. Options are \code{"fd"} for forward differences, and \code{cd}
            for central difference; default is \code{"fd"}.}
        \item{nAGQ}{numeric scalar denoting the number of quadrature points; default is 
            11 when the number of random effects is one or two, and 7 otherwise.}
        \item{update_GH_every}{numeric scalar denoting every how many iterations to 
            update the quadrature points during the EM-phase; default is 10.}
        \item{max_coef_value}{numeric scalar denoting the maximum allowable value for the
            fixed effects coefficients during the optimization; default is 10.}
        \item{max_phis_value}{numeric scalar denoting the maximum allowable value for the
            shape/dispersion parameter of the negative binomial distribution during the 
            optimization; default is \code{exp(10)}.}
        \item{verbose}{logical; print information during the optimization phase; default 
            is \code{FALSE}.}
    }
    }
\item{\dots}{arguments passed to \code{control}.}
}

\details{
\bold{General:} The \code{mixed_model()} function fits mixed effects models in which the
integrals over the random effects in the definition of the marginal log-likelihood cannot
be solved analytically and need to be approximated. The function works under the 
assumption of normally distributed random effects with mean zero and variance-covariance
matrix \eqn{D}. These integrals are approximated numerically using an adaptive 
Gauss-Hermite quadrature rule. Using the control argument \code{nAGQ}, the user can 
specify the number of quadrature points used in the approximation.

\bold{User-defined family:} The user can define its own family object; for an example, 
see the help page of \code{\link{negative.binomial}}.

\bold{Optimization:} A hybrid approach is used, starting with \code{iter_EM} iterations
and unless convergence was achieved it continuous with a direct optimization of the 
log-likelihood using function \code{\link{optim}} and the algorithm specified by 
\code{optim_method}. For stability and speed, the derivative of the log-likelihood with 
respect to the parameters are internally programmed.
}

\value{
An object of class \code{"MixMod"} with components:
\item{coefficients}{a numeric vector with the estimated fixed effects.}
\item{phis}{a numeric vector with the estimated extra parameters.}
\item{D}{a numeric matrix denoting the estimated covariance matrix of the random effects.}
\item{post_modes}{a numeric matrix with the empirical Bayes estimates of the random effects.}
\item{post_vars}{a list of numeric matrices with the posterior variances of the random effects.}
\item{logLik}{a numeric scalar denoting the log-likelihood value at the end of the 
    optimization procedure.}
\item{Hessian}{a numeric matrix denoting the Hessian matrix at the end of the optimization 
    procedure.}
\item{converged}{a logical indicating whether convergence was attained.}
\item{data}{a copy of the \code{data} argument.}
\item{id}{a copy of the grouping variable from \code{data}.}
\item{id_name}{a character string with the name of the grouping variable.}
\item{Terms}{a list with two terms components, \code{termsX} derived from the \code{fixed}
    argument, and \code{termsZ} derived from the \code{random} argument.}
\item{model_frames}{a list with two model.frame components, \code{mfX} derived from the 
    \code{fixed} argument, and \code{mfZ} derived from the \code{random} argument.}
\item{control}{a copy of the (user-specific) \code{control} argument.}
\item{Funs}{a list of functions used in the optimization procedure.}
\item{family}{a copy of the \code{family} argument.}
\item{call}{the matched call.}
}

\author{
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}

\seealso{
\code{\link{methods.MixMod}},
\code{\link{effectPlotData}},
\code{\link{marginal_coefs}}
}


\examples{
# simulate some data
set.seed(123L)
n <- 200
K <- 15
t.max <- 25

betas <- c(-2.13, -0.25, 0.24, -0.05)
D <- matrix(0, 2, 2)
D[1:2, 1:2] <- c(0.48, -0.08, -0.08, 0.18)

times <- c(replicate(n, c(0, sort(runif(K-1, 0, t.max)))))
group <- sample(rep(0:1, each = n/2))
DF <- data.frame(year = times, group = factor(rep(group, each = K)))
X <- model.matrix(~ group * year, data = DF)
Z <- model.matrix(~ year, data = DF)

b <- cbind(rnorm(n, sd = sqrt(D[1, 1])), rnorm(n, sd = sqrt(D[2, 2])))
id <- rep(1:n, each = K)
eta.y <- as.vector(X \%*\% betas + rowSums(Z * b[id, ]))
DF$y <- rbinom(n * K, 1, plogis(eta.y))
DF$id <- factor(id)

################################################

fm1 <- mixed_model(fixed = y ~ year * group, random = ~ 1 | id, data = DF,
                   family = binomial())

# fixed effects
fixef(fm1)

# random effects
head(ranef(fm1))

# detailed output
summary(fm1)

# fitted values for the 'mean subject', i.e., with
# random effects values equal to 0
head(fitted(fm1, type = "mean_subject"))

# fitted values for the conditioning on the estimated random effects
head(fitted(fm1, type = "subject_specific"))

##############
\donttest{
fm2 <- mixed_model(fixed = y ~ year, random = ~ 1 | id, data = DF,
                   family = binomial())

# likelihood ratio test between the two models
anova(fm2, fm1)

# the same hypothesis but with a Wald test
anova(fm1, L = rbind(c(0, 0, 1, 0), c(0, 0, 0, 1)))

##############

# An effects plot for the mean subject (i.e., with random effects equal to 0)
nDF <- with(DF, expand.grid(year = seq(min(year), max(year), length.out = 15),
    group = levels(group)))
    
plot_data <- effectPlotData(fm2, nDF)

require("lattice")
xyplot(pred + low + upp ~ year | group, data = plot_data,
       type = "l", lty = c(1, 2, 2), col = c(2, 1, 1), lwd = 2,
       xlab = "Follow-up time", ylab = "log odds")

expit <- function (x) exp(x) / (1 + exp(x))
xyplot(expit(pred) + expit(low) + expit(upp) ~ year | group, data = plot_data,
       type = "l", lty = c(1, 2, 2), col = c(2, 1, 1), lwd = 2,
       xlab = "Follow-up time", ylab = "Probabilities")

# An effects plots for the marginal probabilities
plot_data_m <- effectPlotData(fm2, nDF, marginal = TRUE, cores = 1L)

expit <- function (x) exp(x) / (1 + exp(x))
xyplot(expit(pred) + expit(low) + expit(upp) ~ year | group, data = plot_data_m,
       type = "l", lty = c(1, 2, 2), col = c(2, 1, 1), lwd = 2,
       xlab = "Follow-up time", ylab = "Probabilities")

##############

# include random slopes
fm1_slp <- update(fm1, random = ~ year | id)

# increase the number of quadrature points to 15
fm1_slp_q15 <- update(fm1_slp, nAGQ = 15)

# a diagonal covariance matrix for the random effects
fm1_slp_diag <- update(fm1, random = ~ year || id)

anova(fm1_slp_diag, fm1_slp)
}
}