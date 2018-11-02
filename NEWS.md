# GLMMadaptive 0.4.0

## General
* The `predict()` method now works for zero-inflated and hurdle models.

* Hurdle Beta mixed effects models are now available using the `hurdle.beta.fam` family
object.

* The new function `scoring_rules()` calculates proper scoring rule for subject-specific
predictions from mixed models for categorical data.

* Added support for the **emmeans** package.

# GLMMadaptive 0.3.0

## General
* Hurdle Poisson and negative binomial models are now implemented using the family objects `hurdle.poisson` and `hurdle.negative.binomial`, respectively.

* added S3 methods for the `terms()`, `model.frame()` and `model.matrix()` generics in order to
work with the **multcomp** package.

* A new vignette illustrating multiple comparisons with the **multcomp** package.

* Methods `vcov()`, `summary()`, `confint()`, `anova()`, `marginal_coefs()`, `effectPlotData()`, `predict()`, and `simulate()` have gained the logical argument `sandwich` to invoke the use of robust/sandwich standard errors in the calculations.

# GLMMadaptive 0.2.0

## General
* Zero-inflated Poisson and negative binomial models are now implemented using the family objects `zi.poisson()` and `zi.negative.binomial()`, respectively. In addition, taking into advantage of the fact that users can specify their own log density functions for the outcome, two-part / hurdle model can also be implemented. 

* A new vignette illustrates how the zero-inflated models can be fitted.

* The `predict()` method is now fully available. It calculates predictions, and standard errors for models returned by `mixed_model()` at three levels:
     + "mean subject": only the fixed effects part corresponding to predictions for the average subject (but not population averaged predictions in case of nonlinear link functions).
     
     + "marginal": predictions using the marginalized coefficients that correspond to population averaged predictions.
     
     + "subject specific": predictions at the subject level. These can be also calculated for subjects not originally in the dataset (i.e., estimates of the random effects are internally obtained).

* The `simulate()` method is available to simulate data from fitted mixed models. This can be used for instance to perform replication / posterior predictive checks.

# GLMMadaptive 0.1.6

## General
* Added vignettes.

# GLMMadaptive 0.1.3

## General
* First version of the package.

