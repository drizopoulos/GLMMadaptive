\name{marginal_coefs}
\alias{marginal_coefs}
\alias{marginal_coefs.MixMod}

\title{
Marginal Coefficients from Generalized Linear Mixed Models
}

\description{
Calculates marginal coefficients and their standard errors from fitted generalized linear 
mixed models.
}

\usage{
marginal_coefs(object, \dots)

\method{marginal_coefs}{MixMod}(object, std_errors = FALSE, 
   link_fun = NULL, M = 5000, K = 100, seed = 1, 
   cores = max(parallel::detectCores() - 1, 1), 
   sandwich = FALSE, \dots)
}

\arguments{
\item{object}{an object inheriting from class \code{"MixMod"}.}
\item{std_errors}{logical indicating whether standard errors are to be computed.}
\item{link_fun}{a function transforming the mean of the repeated measurements outcome to the
    linear predictor scale. Typically, this derived from the \code{family} argument of
    \code{\link{mixed_model}}.}
\item{M}{numeric scalar denoting the number of Monte Carlo samples.}
\item{K}{numeric scalar denoting the number of samples from the sampling distribution of
    the maximum likelihood estimates.}
\item{seed}{integer denoting the seed for the random number generation.}
\item{cores}{integer giving the number of cores to use; applicable only when 
    \code{std_errors = TRUE}.}
\item{sandwich}{logical; if \code{TRUE} robust/sandwich standard errors are used in the 
    calculations.}
\item{\dots}{extra arguments; currently none is used.}
}

\details{
It uses the approach of Hedeker et al. (2017) to calculate marginal coefficients from
mixed models with nonlinear link functions. The marginal probabilities are calculated 
using Monte Carlo integration over the random effects with \code{M} samples, by sampling 
from the estimated prior distribution, i.e., a multivariate normal distribution with mean 
0 and covariance matrix \eqn{\hat{D}}, where \eqn{\hat{D}} denotes the estimated 
covariance matrix of the random effects.

To calculate the standard errors, the Monte Carlo integration procedure is repeated 
\code{K} times, where each time instead of the maximum likelihood estimates of the fixed
effects and the covariance matrix of the random effects, a realization is used from the
sampling distribution of the maximum likelihood estimates. To speed-up this process, 
package \bold{parallel} is used.
}

\value{
A list of class \code{"m_coefs"} with components \code{betas} the marginal coefficients, 
and when \code{std_errors = TRUE}, the extra components \code{var_betas} the estimated 
covariance matrix of the marginal coefficients, and \code{coef_table} a numeric matrix
with the estimated marginal coefficients, their standard errors and corresponding 
p-values using the normal approximation.
}

\references{
Hedeker, D., du Toit, S. H., Demirtas, H. and Gibbons, R. D. (2018), A note on 
marginalization of regression parameters from mixed models of binary outcomes. 
\emph{Biometrics} \bold{74}, 354--361. doi:10.1111/biom.12707
}

\author{
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}

\seealso{
\code{\link{mixed_model}}
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

fm1 <- mixed_model(fixed = y ~ year * group, random = ~ 1 | id, data = DF,
                   family = binomial())

fixef(fm1)                   
marginal_coefs(fm1)
marginal_coefs(fm1, std_errors = TRUE)
}
}