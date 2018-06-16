\name{MixMod Methods}
\alias{methods.MixMod}
\alias{coef}
\alias{coef.MixMod}
\alias{fixef}
\alias{fixef.MixMod}
\alias{ranef}
\alias{ranef.MixMod}
\alias{confint}
\alias{confint.MixMod}
\alias{anova}
\alias{anova.MixMod}
\alias{fitted}
\alias{fitted.MixMod}
\alias{residuals}
\alias{residuals.MixMod}
\alias{predict}
\alias{predict.MixMod}

\title{
Various Methods for Standard Generics
}

\description{
Methods for object of class \code{"MixMod"} for standard generic functions.
}

\usage{

coef(object, \dots)

\method{coef}{MixMod}(object, \dots)

fixef(object, \dots)

\method{fixef}{MixMod}(object, \dots)

ranef(object, \dots)

\method{ranef}{MixMod}(object, post_vars = FALSE, \dots)

confint(object, parm, level = 0.95, \dots)

\method{confint}{MixMod}(object, 
  parm = c("fixed-effects", "var-cov","extra"), 
  level = 0.95, \dots)

anova(object, \dots)

\method{anova}{MixMod}(object, object2, test = TRUE, L = NULL, \dots)

fitted(object, \dots)

\method{fitted}{MixMod}(object, 
  type = c("mean_subject", "subject_specific", "marginal"),
  link_fun = NULL, \dots)

residuals(object, \dots)

\method{residuals}{MixMod}(object, 
  type = c("mean_subject", "subject_specific", "marginal"), 
  link_fun = NULL, \dots)
  
predict(object, \dots)

\method{predict}{MixMod}(object, newdata, newdata2 = NULL, 
  type = c("link", "response"),
  level = c("mean_subject", "subject_specific", "marginal"),
  se.fit = FALSE, M = 200, df = 10, scale = 0.5, CI_level = 0.05, 
  seed = 1, \dots)
}

\arguments{
\item{object, object2}{objects inheriting from class \code{"MixMod"}. When \code{object2}
  is also provided, then the model behind \code{object} must be nested within the model 
  behind \code{object2}.}
\item{post_vars}{logical; if \code{TRUE} the posterior variances of the random effects
  are returned as an extra attribute of the numeric matrix produced by \code{ranef()}.}
\item{parm}{character string; for which type of parameters to calculate confidence 
  intervals. Option \code{"var-cov"} corresponds to the variance-covariance matrix of the 
  random effects. Option \code{extra} corresponds to extra (shape/dispersion) 
  parameters in the distribution of the outcome (e.g., the \eqn{\theta} parameter in
  the negative binomial family).}
\item{level, CI_level}{numeric scalar between 0 and 1 denoting the level of the confidence interval.}
\item{test}{logical; should a p-value be calculated.}
\item{L}{a numeric matrix representing a contrasts matrix. This is only used when in 
  \code{anova()} only \code{object} is provided, and it can only be specified for the fixed 
  effects. When \code{L} is used, a Wald test is performed.}
\item{type}{character string indicating the type of fitted values / residuals to calculate.
    Option \code{"mean_subject"} corresponds to only using the fixed-effects part; option
    \code{"subject_specific"} corresponds to using both the fixed- and random-effects 
    parts; option \code{"marginal"} is based in multiplying the fixed effects design 
    matrix with the marginal coefficients obtained by \code{\link{marginal_coefs}}.}
\item{link_fun}{the \code{link_fun} of \code{\link{marginal_coefs}}.}
\item{newdata, newdata2}{a data frame based on which predictions are to be calculated.
    \code{newdata2} is only relevant when \code{level = "subject_specific"}; see 
    \bold{Details} for more information.}
\item{se.fit}{logical, if \code{TRUE} standard errors of predictions are returned.}
\item{M}{numeric scalar denoting the number of Monte Carlo samples; see 
    \bold{Details} for more information.}
\item{df}{numeric scalar denoting the degrees of freedom for the Student's t proposal
    distribution; see \bold{Details} for more information.}
\item{scale}{numeric scalar or vector denoting the scaling applied to the subject-specific
    covariance matrices of the random effects; see \bold{Details} for more information.}
\item{seed}{numerical scalar giving the seed to used in the Monte Carlo scheme.}
\item{\dots}{further arguments; currently none is used.}
}

\details{
 ***.
}


\value{
The estimated fixed and random effects, coefficients (this is similar as in package 
\bold{nlme}), confidence intervals fitted values (on the scale on the response) and 
residuals.
}

\author{
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}

\seealso{
\code{\link{mixed_model}},
\code{\link{marginal_coefs}}
}

\examples{
\donttest{
# simulate some data
set.seed(123L)
n <- 500
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

fm1 <- mixed_model(fixed = y ~ year + group, random = ~ year | id, data = DF,
                   family = binomial())

head(coef(fm1))
fixef(fm1)
head(ranef(fm1))


confint(fm1)
confint(fm1, "var-cov")

head(fitted(fm1, "subject_specific"))
head(residuals(fm1, "marginal"))

fm2 <- mixed_model(fixed = y ~ year * group, random = ~ year | id, data = DF,
                   family = binomial())

# likelihood ratio test between fm1 and fm2
anova(fm1, fm2)

# the same but with a Wald test
anova(fm2, L = rbind(c(0, 0, 0, 1)))
}
}
