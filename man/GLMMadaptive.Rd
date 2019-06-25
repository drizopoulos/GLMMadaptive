\name{GLMMadaptive}
\alias{GLMMadaptive-package}
\alias{GLMMadaptive}
\docType{package}

\title{
Generalized Linear Mixed Models using Adaptive Gaussian Quadrature
}

\description{
This package fits generalized linear mixed models for a single grouping factor under 
maximum likelihood approximating the integrals over the random effects with an adaptive 
Gaussian quadrature rule.
}

\details{
\tabular{ll}{
Package: \tab GLMMadaptive\cr
Type: \tab Package\cr
Version: \tab 0.6-0\cr
Date: \tab 2019-06-25\cr
License: \tab GPL (>=3)\cr
}

This package fits mixed effects models for grouped / repeated measurements data for which
the integral over the random effects in the definition of the marginal likelihood cannot
be solved analytically. The package approximates these integrals using the adaptive 
Gauss-Hermite quadrature rule.

Multiple random effects terms can be included for the grouping factor (e.g., random 
intercepts, random linear slopes, random quadratic slopes), but currently only a single
grouping factor is allowed.

The package also offers several utility functions that can extract useful information from
fitted mixed effects models. The most important of those are included in the 
\bold{See also} Section below.
}

\author{
Dimitris Rizopoulos

Maintainer: Dimitris Rizopoulos <d.rizopoulos@erasmusmc.nl>
}

\references{
Fitzmaurice, G., Laird, N. and Ware J. (2011). \emph{Applied Longitudinal Analysis}, 
2nd Ed. New York: John Wiley & Sons. 

Molenberghs, G. and Verbeke, G. (2005). \emph{Models for Discrete Longitudinal Data}. 
New York: Springer-Verlag. 
}

\keyword{package}
\keyword{multivariate}

\seealso{
    \code{\link{mixed_model}},
    \code{\link{methods.MixMod}},
    \code{\link{effectPlotData}},
    \code{\link{marginal_coefs}}
}

