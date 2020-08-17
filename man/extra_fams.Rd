\name{Extra Family Objects}
\alias{students.t}
\alias{beta.fam}
\alias{beta.binomial}
\alias{zi.poisson}
\alias{zi.negative.binomial}
\alias{hurdle.poisson}
\alias{hurdle.negative.binomial}
\alias{hurdle.lognormal}
\alias{hurdle.beta.fam}
\alias{unit.lindley}
\alias{Gamma.fam}


\title{
Family functions for Student's-t, Beta, Zero-Inflated and Hurdle Poisson and Negative Binomial,
Hurdle Log-Normal, and Hurdle Beta Mixed Models}

\description{
Specifies the information required to fit a Beta, zero-inflated and hurdle Poisson, 
zero-inflated and hurdle Negative Binomial, a hurdle normal and a hurdle Beta 
mixed-effects model, using 
\code{mixed_model()}.
}

\usage{
students.t(df = stop("'df' must be specified"), link = "identity")
beta.fam()
zi.poisson()
zi.negative.binomial()
hurdle.poisson()
hurdle.negative.binomial()
hurdle.lognormal()
hurdle.beta.fam()
unit.lindley()
beta.binomial(link = "logit")
Gamma.fam()
}

\arguments{
\item{link}{name of the link function.}
\item{df}{the degrees of freedom of the Student's t distribution.}
}

\note{Currently only the log-link is implemented for the Poisson, negative binomial 
and Gamma models, the logit link for the beta and hurdle beta models and the identity link for the 
log-normal model.}

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

#############################################################################
#############################################################################

# The same as above but with a negative binomial model
gm1 <- mixed_model(fixed = y ~ f1 * f2, random = ~ 1 | g, data = dd, 
                  family = zi.negative.binomial(), zi_fixed = ~ 1)

summary(gm1)

# We do a likelihood ratio test between the Poisson and negative binomial models
anova(fm1, gm1)
}
}
