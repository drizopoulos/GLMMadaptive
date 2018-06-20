\name{zi.poisson}
\alias{zi.poisson}

\title{
Family function for Negative Binomial GLMs
}

\description{
Specifies the information required to fit a Negative Binomial generalized linear model,
with known \code{theta} parameter, using \code{glm()}.
}

\usage{
zi.poisson(link = "log")
}

\arguments{
\item{link}{The link function, as a character string or name; currently only \code{log} 
is allowed.}
}

\examples{
# simulate some data from a negative binomial model
set.seed(102)
dd <- expand.grid(f1 = factor(1:3), f2 = LETTERS[1:2], g = 1:30, rep = 1:15,
                  KEEP.OUT.ATTRS = FALSE)
mu <- 5*(-4 + with(dd, as.integer(f1) + 4 * as.numeric(f2)))
dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)

# Fit a zero-inflated Poisson model, with only fixed effects in the 
# zero-inflated part
fm1 <- mixed_model(fixed = y ~ f1 * f2, random = ~ 1 | g, data = dd, 
                  family = zi.poisson(), zi_fixed = ~ 1)

summary(fm1)

\donttest{
# We extend the previous model allowing also for a random intercept in the
# zero-inflated part
fm2 <- mixed_model(fixed = y ~ f1 * f2, random = ~ 1 | g, data = dd, 
                  family = zi.poisson(), zi_fixed = ~ 1, zi_random = ~ 1 | g)

# We do a likelihood ratio test between the two models
anova(fm1, fm2)
}
}
