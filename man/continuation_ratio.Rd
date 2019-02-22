\name{Continuation Ratio Set-Up}
\alias{cr_setup}
\alias{cr_marg_probs}

\title{
Functions to Set-Up Data for a Continuation Ratio Mixed Model
}

\description{
Data set-up and calculation of marginal probabilities from a continuation ratio model
}

\usage{
cr_setup(y, direction = c("forward", "backward"))

cr_marg_probs(eta, direction = c("forward", "backward"))
}

\arguments{
  \item{y}{a numeric vector denoting the ordinal response variable.}
  \item{direction}{character string specifying the direction of the continuation ratio 
    model; \code{"forward"} corresponds to a discrete hazard function.}
  \item{eta}{a numeric matrix of the linear predictor, with columns corresponding to the 
    different levels of the ordinal response.}
}

\author{
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}

Frank Harrell
}

\note{Function \code{cr_setup()} is based on the \code{cr.setup()} function from package 
\bold{rms}.}

\examples{
n <- 300 # number of subjects
K <- 8 # number of measurements per subject
t_max <- 15 # maximum follow-up time

# we constuct a data frame with the design: 
# everyone has a baseline measurment, and then measurements at random follow-up times
DF <- data.frame(id = rep(seq_len(n), each = K),
                 time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

# design matrices for the fixed and random effects
X <- model.matrix(~ sex * time, data = DF)[, -1]
Z <- model.matrix(~ 1, data = DF)

thrs <- c(-1.5, 0, 0.9) # thresholds for the different ordinal categories
betas <- c(-0.25, 0.24, -0.05) # fixed effects coefficients
D11 <- 0.48 # variance of random intercepts
D22 <- 0.1 # variance of random slopes

# we simulate random effects
b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)))[, 1, drop = FALSE]
# linear predictor
eta_y <- drop(X \%*\% betas + rowSums(Z * b[DF$id, , drop = FALSE]))
# linear predictor for each category
eta_y <- outer(eta_y, thrs, "+")
# marginal probabilities per category
mprobs <- cr_marg_probs(eta_y)
# we simulate ordinal longitudinal data
DF$y <- unname(apply(mprobs, 1, sample, x = ncol(mprobs), size = 1, replace = TRUE))

# If you want to simulate from the backward formulation of the CR model, you need to
# change `eta_y <- outer(eta_y, thrs, "+")` to `eta_y <- outer(eta_y, rev(thrs), "+")`,
# and `mprobs <- cr_marg_probs(eta_y)` to `mprobs <- cr_marg_probs(eta_y, "backward")`

#################################################

# prepare the data
cr_vals <- cr_setup(DF$y)
cr_data <- DF[cr_vals$subs, ]
cr_data$y_new <- cr_vals$y
cr_data$cohort <- cr_vals$cohort
fm <- mixed_model(y_new ~ cohort + sex * time, random = ~ 1 | id, 
                  data = cr_data, family = binomial())

summary(fm)

# If you want to fit the CR model under backward formulation, you need to change
# `cr_vals <- cr_setup(DF$y)` to `cr_vals <- cr_setup(DF$y, "backward")`
}
