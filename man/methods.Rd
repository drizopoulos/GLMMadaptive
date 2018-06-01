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
\alias{fitted}
\alias{fitted.MixMod}
\alias{residuals}
\alias{residuals.MixMod}

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

\method{ranef}{MixMod}(object, \dots)

confint(object, parm, level = 0.95, \dots)

\method{confint}{MixMod}(object, 
  parm = c("fixed-effects", "var-cov","extra"), 
  level = 0.95, \dots)

fitted(object, \dots)

\method{fitted}{MixMod}(object, 
  type = c("mean_subject", "subject_specific", "marginal"),
  link_fun = NULL, \dots)

residuals(object, \dots)

\method{residuals}{MixMod}(object, 
  type = c("mean_subject", "subject_specific", "marginal"), 
  link_fun = NULL, \dots)
}

\arguments{
\item{object}{an object inheriting from class \code{"MixMod"}.}
\item{parm}{character string; for which type of parameters to calculate confidence 
  intervals. Option \code{"var-cov"} corresponds to the variance-covariance matrix of the 
  random effects. Option \code{extra} corresponds to extra (shape/dispersion) 
  parameters in the distribution of the outcome (e.g., the \eqn{\theta} parameter in
  the negative binomial family).}
\item{level}{numeric scalar between 0 and 1 denoting the level of the confidence interval.}
\item{type}{character string indicating the type of fitted values / residuals to calculate.
    Option \code{"mean_subject"} corresponds to only using the fixed-effects part; option
    \code{"subject_specific"} corresponds to using both the fixed- and random-effects 
    parts; option \code{"marginal"} is based in multiplying the fixed effects design 
    matrix with the marginal coefficients obtained by \code{\link{marginal_coefs}}.}
\item{link_fun}{the \code{link_fun} of \code{\link{marginal_coefs}}.}
\item{\dots}{further arguments; currently none is used.}
}

\value{
The estimated fixed and random effects, coefficients (this is similar as in package 
\bold{nlme}), confidence intervals fitted values (on the scase on the response) and 
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
\dontrun{
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
confint(fm1, "var")

head(fitted(fm1, "subject_specific"))
head(residuals(fm1, "marginal"))
}
}
