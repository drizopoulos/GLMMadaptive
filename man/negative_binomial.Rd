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
