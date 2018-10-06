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
\alias{simulate}
\alias{simulate.MixMod}
\alias{terms}
\alias{terms.MixMod}
\alias{model.matrix}
\alias{model.matrix.MixMod}
\alias{model.frame}
\alias{model.frame.MixMod}

\title{
Various Methods for Standard Generics
}

\description{
Methods for object of class \code{"MixMod"} for standard generic functions.
}

\usage{

coef(object, \dots)

\method{coef}{MixMod}(object, sub_model = c("main", "zero_part"), 
    \dots)

fixef(object, \dots)

\method{fixef}{MixMod}(object, sub_model = c("main", "zero_part"), \dots)

ranef(object, \dots)

\method{ranef}{MixMod}(object, post_vars = FALSE, \dots)

confint(object, parm, level = 0.95, \dots)

\method{confint}{MixMod}(object, 
  parm = c("fixed-effects", "var-cov","extra", "zero_part"), 
  level = 0.95, sandwich = FALSE, \dots)

anova(object, \dots)

\method{anova}{MixMod}(object, object2, test = TRUE, 
  L = NULL, sandwich = FALSE, \dots)

fitted(object, \dots)

\method{fitted}{MixMod}(object, 
  type = c("mean_subject", "subject_specific", "marginal"),
  link_fun = NULL, \dots)

residuals(object, \dots)

\method{residuals}{MixMod}(object, 
  type = c("mean_subject", "subject_specific", "marginal"), 
  link_fun = NULL, tasnf_y = function (x) x, \dots)
  
predict(object, \dots)

\method{predict}{MixMod}(object, newdata, newdata2 = NULL, 
    type_pred = c("response", "link"),
    type = c("mean_subject", "subject_specific", "marginal", "zero_part"),
    se.fit = FALSE, M = 300, df = 10, scale = 0.3, level = 0.95, 
    seed = 1, sandwich = FALSE, \dots)
    
simulate(object, nsim = 1, seed = NULL, \dots)

\method{simulate}{MixMod}(object, nsim = 1, seed = NULL, 
    acount_MLEs_var = FALSE, sim_fun = NULL, 
    sandwich = FALSE, \dots)
    
terms(x, \dots)

\method{terms}{MixMod}(x, type = c("fixed", "random"), \dots)

model.frame(formula, \dots)

\method{model.frame}{MixMod}(formula, type = c("fixed", "random"), \dots)

model.matrix(object, \dots)

\method{model.matrix}{MixMod}(object, type = c("fixed", "random"), \dots)
}

\arguments{
\item{object, object2, x, formula}{objects inheriting from class \code{"MixMod"}. When \code{object2}
  is also provided, then the model behind \code{object} must be nested within the model 
  behind \code{object2}.}
\item{sub_model}{character string indicating for which sub-model to extract the estimated
  coefficients; it is only relevant for zero-inflated models.}
\item{post_vars}{logical; if \code{TRUE} the posterior variances of the random effects
  are returned as an extra attribute of the numeric matrix produced by \code{ranef()}.}
\item{parm}{character string; for which type of parameters to calculate confidence 
  intervals. Option \code{"var-cov"} corresponds to the variance-covariance matrix of the 
  random effects. Option \code{extra} corresponds to extra (shape/dispersion) 
  parameters in the distribution of the outcome (e.g., the \eqn{\theta} parameter in
  the negative binomial family). Option \code{zero_inflated} corresponds to the 
  coefficients of the zero-inflated sub-model.}
\item{level}{numeric scalar between 0 and 1 denoting the level of the confidence interval.}
\item{test}{logical; should a p-value be calculated.}
\item{L}{a numeric matrix representing a contrasts matrix. This is only used when in 
  \code{anova()} only \code{object} is provided, and it can only be specified for the fixed 
  effects. When \code{L} is used, a Wald test is performed.}
\item{sandwich}{logical; if \code{TRUE} the sandwich estimator is used in the calculation
  of standard errors.}
\item{type}{character string indicating the type of fitted values / residuals / predictions
    to calculate. Option \code{"mean_subject"} corresponds to only using the fixed-effects
    part; option \code{"subject_specific"} corresponds to using both the fixed- and 
    random-effects parts; option \code{"marginal"} is based in multiplying the fixed 
    effects design matrix with the marginal coefficients obtained by 
    \code{\link{marginal_coefs}}.}
\item{link_fun}{the \code{link_fun} of \code{\link{marginal_coefs}}.}
\item{tasnf_y}{a function to transform the grouped / repeated measurements outcome before
    calculating the residuals; for example, relevant in two-part models for semi-continuous
    data, in which it is assumed that the log outcome follows a normal distribution.}
\item{newdata, newdata2}{a data frame based on which predictions are to be calculated.
    \code{newdata2} is only relevant when \code{level = "subject_specific"}; see 
    \bold{Details} for more information.}
\item{type_pred}{character string indicating at which scale to calculate predictions. 
    Options are \code{"link"} indicating to calculate predictions at the link function / 
    linear predictor scale, and \code{"response"} indicating to calculate predictions at
    the scale of the response variable.}
\item{se.fit}{logical, if \code{TRUE} standard errors of predictions are returned.}
\item{M}{numeric scalar denoting the number of Monte Carlo samples; see 
    \bold{Details} for more information.}
\item{df}{numeric scalar denoting the degrees of freedom for the Student's t proposal
    distribution; see \bold{Details} for more information.}
\item{scale}{numeric scalar or vector denoting the scaling applied to the subject-specific
    covariance matrices of the random effects; see \bold{Details} for more information.}
\item{seed}{numerical scalar giving the seed to be used in the Monte Carlo scheme.}
\item{nsim}{numeric scalar giving the number of times to simulate the response variable.}
\item{acount_MLEs_var}{logical; if \code{TRUE} it accounts for the variability of the 
    maximum likelihood estimates (MLEs) by simulating a new value for the parameters from 
    a multivariate normal distribution with mean the MLEs and covariance matrix the 
    covariance matrix of the MLEs.}
\item{sim_fun}{a function based on which to simulate the response variable. This is 
    relevant for non-standard models. The \code{simulate()} function also tries to extract 
    this function from the \code{family} component of \code{object}. The function should 
    have the following four arguments: \code{n} a numeric scalar denoting the number of
    observations to simulate, \code{mu} a numeric vector of means, \code{phis} a numeric
    vector of extra dispersion/scale parameters, and \code{eta_zi} a numeric vector for 
    the zero-part of the model, if this is relevant.}
\item{\dots}{further arguments; currently none is used.}
}

\details{
 In generic terms, we assume that the mean of the outcome \eqn{y_i} (\eqn{i = 1, ..., n} 
 denotes the subjects) conditional on the random effects is given by the equation: 
 \deqn{g{E(y_i | b_i)} = \eta_i = X_i \beta + Z_i b_i,}
 where \eqn{g(.)} denotes the link function, \eqn{b_i} the vector of random effects, 
 \eqn{\beta} the vector of fixed effects, and \eqn{X_i} and \eqn{Z_i} the design matrices
 for the fixed and random effects, respectively.
 
 Argument \code{type_pred} of \code{predict()} specifies whether predictions will be 
 calculated in the link / linear predictor scale, i.e., \eqn{\eta_i} or in the response 
 scale, i.e., \eqn{g{E(y_i | b_i)}}.
 
 When \code{type = "mean_subject"}, predictions are calculated using only the fixed 
 effects, i.e., the \eqn{X_i \beta} part, where \eqn{X_i} is evaluated in \code{newdata}. 
 This corresponds to predictions for the 'mean' subjects, i.e., subjects who have 
 random effects value equal to 0. Note, that in the case of nonlinear link functions this
 does not correspond to the averaged over the population predictions (i.e., marginal 
 predictions).
 
 When \code{type = "marginal"}, predictions are calculated using only the fixed 
 effects, i.e., the \eqn{X_i \beta} part, where \eqn{X_i} is evaluated in \code{newdata},
 but with \eqn{\beta} coefficients the marginalized coefficients obtain from 
 \code{\link{marginal_coefs}}. These predictions will be marginal / population averaged 
 predictions.
 
 When \code{type = "subject_specific"}, predictions are calculated using both the fixed- 
 and random-effects parts, i.e., \eqn{X_i \beta + Z_i b_i}, where \eqn{X_i} and \eqn{Z_i} 
 are evaluated in \code{newdata}. Estimates for the random effects of each subject are 
 obtained as modes from the posterior distribution \eqn{[b_i | y_i; \theta]} evaluated in 
 \code{newdata} and with \eqn{theta} (denoting the parameters of the model, fixed effects
 and variance components) replaced by their maximum likelihood estimates. 
 
 When \code{se.fit = TRUE} and \code{type = "subject_specific"}, standard errors and 
 confidence intervals for the subject-specific predictions are obtained by a Monte Carlo 
 scheme entailing three steps repeated \code{M} times, namely
 \describe{
    \item{Step I}{Account for the variability of maximum likelihood estimates (MLES) by 
        simulating a new value \eqn{\theta^*} for the parameters \eqn{\theta} from a 
        multivariate normal distribution with mean the MLEs and covariance matrix the 
        covariance matrix of the MLEs.}
    
    \item{Step II}{Account for the variability in the random effects estimates by 
    simulating a new value \eqn{b_i^*} for the random effects \eqn{b_i} from the posterior
    distribution \eqn{[b_i | y_i; \theta^*]}. Because the posterior distribution does not
    have a closed-form, a Metropolis-Hastings algorithm is used to sample the new value 
    \eqn{b_i^*} using as proposal distribution a multivariate Student's-t distribution
    with degrees of freedom \code{df}, centered at the mode of the posterior distribution
    \eqn{[b_i | y_i; \theta]} with \eqn{\theta} the MLEs, and scale matrix the inverse
    Hessian matrix of the log density of \eqn{[b_i | y_i; \theta]} evaluated at the modes,
    but multiplied by \code{scale}. The \code{scale} and \code{df} parameters can be used 
    to adjust the acceptance rate.}
    
    \item{Step III}{The predictions are calculated using \eqn{X_i \beta^* + Z_i b_i^*}.}
 }
 
 Argument \code{newdata2} can be used to calculate dynamic subject-specific predictions. 
 I.e., using the observed responses \eqn{y_i} in \code{newdata}, estimates of the random
 effects of each subject are obtained. For the same subjects we want to obtain predictions
 in new covariates settings for which no response data are yet available. For example, 
 in a longitudinal study, for a subject we have responses up to a follow-up \eqn{t} 
 (\code{newdata}) and we want the prediction at \eqn{t + \Delta t} (\code{newdata2}).
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
