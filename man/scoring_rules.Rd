\name{scoring_rules}
\alias{scoring_rules}

\title{
Proper Scoring Rules for Categorical Data
}

\description{
Calculates the logarithmic, quadratic/Brier and spherical based on a fitted mixed model
for categorical data.
}

\usage{
scoring_rules(object, newdata, newdata2 = NULL, max_count = 2000, 
    return_newdata = FALSE)
}

\arguments{
\item{object}{an object inheriting from class \code{"MixMod"}.}
\item{newdata}{a data.frame based on which to estimate the random effect and calculate 
    predictions. It should contain the response variable.}
\item{newdata2}{a data.frame based on which to estimate the random effect and calculate 
    predictions. It should contain the response variable.}
\item{max_count}{numeric scalar denoting the maximum count up to which to calculate 
    probabilities; this is relevant for count response data.}
\item{return_newdata}{logical; if \code{TRUE} the values of the scoring rules are 
    ruturned as extra columns of the \code{newdata} or \code{newdata2} data.frame.}
}

\value{
A data.frame with (extra) columns the values of the logarithmic, quadratic and spherical
scoring rules calculated based on the fitted model and the observed responses in 
\code{newdata} or \code{newdata2}.
}

\references{
Carvalho, A. (2016). An overview of applications of proper scoring rules. 
\emph{Decision Analysis} \bold{13}, 223--242. doi:10.1287/deca.2016.0337
}

\author{
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}

\seealso{
\code{\link{mixed_model}},
\code{\link{predict.MixMod}}
}


\examples{
\donttest{
NA
NA
NA
}
}