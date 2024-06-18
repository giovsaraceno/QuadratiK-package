QuadratiK (Next Release) 
=========================

### NEW FEATURES

  * Added NEWS.md file
  * Added badge for rOpenSci review status


QuadratiK 1.1.1 (2024-06-05)
=========================

### NEW FEATURES

  * Improve computation of variance for k-sample test
  * Add `breast_cancer` and `wine` data sets
  
### MINOR IMPROVEMENTS

  * `print()` and `summary()` methods of the `kb.test` object for the two and k-sample tests
  
### BUG FIXES

  * Computation of test statistics for normality test

### DOCUMENTATION FIXES

  * Correction of typos
  

QuadratiK 1.1.0 (2024-05-14)
=========================

### NEW FEATURES

  * Add GitHub links to the DESCRIPTION file
  * Add unit tests  
  * Add CONTRIBUTING.md
  * Add srr standards for the categories: General, Clustering and Probability Distributions
  * Add methods for the clustering object `pkbc`: `show()`, `summary()`, `estract_stats()`, `plot()` and `predict()`
  * Add vignette `generate_rpkb.Rmd`
  * Add CODE_OF_CONDUCT.md
  * Create the GitHub page using pkgdown
  * Add V-statistic for Normality test, with corresponding Critical Value and Degrees of Freedom with the function `DOF_norm()`
  * Update `kb.test` S4 class slots, `kb.test()` function, unit-test codes, `summary()` and `print()` method for including the V-statistic
  * Add Trace statistic for two and k-sample tests
  * Add computation of the Variance for two and k-sample tests

### MINOR IMPROVEMENTS

  * Addition of GitHub workflows for R-CMD-check, test-coverage and pkgcheck
  * Update codes to check arguments of function `pk.test()`
  * Update codes to check arguments of function `kb.test()` and `select_h()`
  * Modify indentation and assignments (using '<-') following `goodpractice` standards 
  * Update references in README.md file

### BUG FIXES

  * Add codemeta.json to .Rbuildignore
  * Add CONTRIBUTING.md to .Rbuildignore
  * Fix warning from internal function `compare_qq()`
  * Add `doc` folder to .gitignore
  
  * Modify returned test statistics multiplying them by the sample size
  * Remove "Nonparam" centering option for the Normality test
  * Critical value of the V-statistic divided by the variance in the function `pk.test()`

### DEPRECATED AND DEFUNCT

  * `summary_stat()` and `validation()` substituted by `estract_stats()` and `pkbc_validation()`

### DOCUMENTATION FIXES

  * C++ functions are set as internal functions
  * Correction of typos
  * Remove Date from DESCRIPTION
  * Update References in the DESCRIPTION file


QuadratiK 1.0.0 (2024-02-23)
=========================

### NEW FEATURES

  * released to CRAN
