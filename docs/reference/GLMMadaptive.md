# Generalized Linear Mixed Models using Adaptive Gaussian Quadrature

This package fits generalized linear mixed models for a single grouping
factor under maximum likelihood approximating the integrals over the
random effects with an adaptive Gaussian quadrature rule.

## Details

|          |              |
|----------|--------------|
| Package: | GLMMadaptive |
| Type:    | Package      |
| Version: | 0.9-9        |
| Date:    | 2026-11-17   |
| License: | GPL (\>=3)   |

This package fits mixed effects models for grouped / repeated
measurements data for which the integral over the random effects in the
definition of the marginal likelihood cannot be solved analytically. The
package approximates these integrals using the adaptive Gauss-Hermite
quadrature rule.

Multiple random effects terms can be included for the grouping factor
(e.g., random intercepts, random linear slopes, random quadratic
slopes), but currently only a single grouping factor is allowed.

The package also offers several utility functions that can extract
useful information from fitted mixed effects models. The most important
of those are included in the **See also** Section below.

## Author

Dimitris Rizopoulos

Maintainer: Dimitris Rizopoulos \<d.rizopoulos@erasmusmc.nl\>

## References

Fitzmaurice, G., Laird, N. and Ware J. (2011). *Applied Longitudinal
Analysis*, 2nd Ed. New York: John Wiley & Sons.

Molenberghs, G. and Verbeke, G. (2005). *Models for Discrete
Longitudinal Data*. New York: Springer-Verlag.

## See also

[`mixed_model`](https://drizopoulos.github.io/GLMMadaptive/reference/mixed_model.md),
[`methods.MixMod`](https://drizopoulos.github.io/GLMMadaptive/reference/methods.md),
[`effectPlotData`](https://drizopoulos.github.io/GLMMadaptive/reference/effectPlotData.md),
[`marginal_coefs`](https://drizopoulos.github.io/GLMMadaptive/reference/marginal_coefs.md)
