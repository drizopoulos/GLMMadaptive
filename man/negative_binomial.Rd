\name{negative.binomial}
\alias{negative.binomial}

\title{
Family function for Negative Binomial GLMs
}

\description{
Specifies the information required to fit a Negative Binomial generalized linear model,
with known \code{theta} parameter, using \code{glm()}.
}

\usage{
negative.binomial(theta = stop("'theta' must be specified"), link = "log")
}

\arguments{
\item{theta}{The known value of the additional parameter \code{theta}.}
\item{link}{The link function, as a character string, name or one-element
    character vector specifying one of \code{log}, \code{sqrt}
    or \code{identity}, or an object of class
\code{"\link[=family]{link-glm}"}.}
}

\references{
Venables, W. N. and Ripley, B. D. (1999) \emph{Modern Applied Statistics with S-PLUS.} 
Third Edition. Springer.
}

\author{
The authors of the \bold{MASS} package, with a small addition from 
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}

\examples{
# simulate some data
set.seed(101)
dd <- expand.grid(f1 = factor(1:3), f2 = LETTERS[1:2], g = 1:30, rep = 1:15,
                  KEEP.OUT.ATTRS = FALSE)
mu <- 5*(-4 + with(dd, as.integer(f1) + 4 * as.numeric(f2)))
dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)

gm <-  mixed_model(fixed = y ~ f1 * f2, random = ~ 1 | g, data = dd, 
                   family = negative.binomial(0.5))

summary(gm)


my_negBinom <- function (link = "log") {
    stats <- make.link(link)
    log_dens <- function (y, eta, mu_fun, phis) {
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
    score_eta_fun <- function (y, mu, phis) {
        # the derivative of the log density w.r.t. mu
        phis <- exp(phis)
        mu_phis <- mu + phis
        comp2 <- - phis / mu_phis
        comp3 <- y / mu - y / mu_phis
        # the derivative of mu w.r.t. eta (this depends on the chosen link function)
        mu.eta <- mu
        (comp2 + comp3) * mu.eta
    }
    score_phis_fun <- function (y, mu, phis) {
        # the derivative of the log density w.r.t. phis
        phis <- exp(phis)
        mu_phis <- mu + phis
        comp1 <- digamma(y + phis) - digamma(phis)
        comp2 <- log(phis) + 1 - log(mu_phis) - phis / mu_phis
        comp3 <- - y / mu_phis
        comp1 + comp2 + comp3
    }
    structure(list(family = "user Neg Binom", link = stats$name, linkfun = stats$linkfun,
                   linkinv = stats$linkinv, log_dens = log_dens,
                   score_eta_fun = score_eta_fun,
                   score_phis_fun = score_phis_fun),
              class = "family")
}

fm <-  mixed_model(fixed = y ~ f1 * f2, random = ~ 1 | g, data = dd, 
                   family = my_negBinom(), n_phis = 1, 
                   initial_values = list("betas" = poisson()))

summary(fm)
}
