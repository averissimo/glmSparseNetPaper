
> Support R-:package: for paper Veríssimo et al. (2020)

This package contains all methods and workflow necessary to reproduce
the results in the paper using the vignettes in the folder with the same
name:

  - `Analysis.Rmd`: Choosing an alpha value and a target number of
    variables, it compares between Network-based models and classical
    ElasticNet approach. It also performs a stability study on the
    selected variables.
  - `OptimalVariableSize.Rmd`: Performs cross-validation on multiple
    values of `alpha` and shows Kaplan-Meier estimators and Concordance
    C-Index values.
  - `Normalization_tests.Rmd`: Test different methods to normalize the
    data and then performs cross-validation to help decide on which to
    use.
