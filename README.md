<!-- badges: start -->

| Usage  | Release | Development |
|------------------|------------------------|------------------------------|
| [![CRAN_Downloads_Badge](https://cranlogs.r-pkg.org/badges/grand-total/QuadratiK)](https://cran.r-project.org/package=QuadratiK) [![License: GPL v3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://docs.ropensci.org/QuadratiK/LICENSE.html)|[![arXiv](https://img.shields.io/badge/doi-arXiv:2402.02290v2-green.svg)](https://doi.org/10.48550/arXiv.2402.02290) [![CRAN version](https://www.r-pkg.org/badges/version/QuadratiK)](https://CRAN.R-project.org/package=QuadratiK) [![GitHub version](https://img.shields.io/badge/devel%20version-1.1.3-blue.svg)](https://github.com/ropensci/QuadratiK) | [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) [![Status at rOpenSci Software Peer Review](https://badges.ropensci.org/632_status.svg)](https://github.com/ropensci/software-review/issues/632)  [![R-CMD-check](https://github.com/ropensci/QuadratiK/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ropensci/QuadratiK/actions/workflows/R-CMD-check.yaml) [![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable) |             |

<!-- badges: end -->

# Collection of Methods Constructed using the Kernel-Based Quadratic Distances

`QuadratiK` provides the first implementation, in R and Python, of a comprehensive set of goodness-of-fit tests and a clustering technique for $d$-dimensional spherical data $d \ge 2$ using kernel-based quadratic distances. It includes:

-   **Goodness-of-Fit Tests**: The software implements one, two, and *k*-sample tests for goodness of fit, offering an efficient and mathematically sound way to assess the fit of probability distributions. Our tests are particularly useful for large, high dimensional data sets where the assessment of fit of probability models is of interest. Specifically, we offer tests for multivariate normality, as well as two- and *k*-sample tests, where testing equality of two or more distributions is of interest, that is $H_0: F_1 = F_2$ and $H_0: F_1 = \ldots = F_k$ respectively. The proposed tests perform well in terms of level and power for contiguous alternatives, heavy tailed distributions and in higher dimensions.\
Expanded capabilities include supporting tests for uniformity on the *d*-dimensional Sphere based on the Poisson kernel, exhibiting excellent results especially in the case of multimodal distributions.

-   **Poisson kernel-based distribution (PKBD)**: the package offers functions for computing the density value and for generating random samples from a PKBD. The Poisson kernel-based densities are based on the normalized Poisson kernel and are defined on the $d$-dimensional unit sphere. Given a vector $\mu \in \mathcal{S}^{d-1}$, and a parameter $\rho$ such that $0 < \rho < 1$, the probability density function of a $d$-variate Poisson kernel-based density is defined by: $$f(\mathbf{x}|\rho, \mathbf{\mu}) = \frac{1-\rho^2}{\omega_d ||\mathbf{x} - \rho \mathbf{\mu}||^d},$$ where $\mu$ is a vector orienting the center of the distribution, $\rho$ is a parameter to control the concentration of the distribution around the vector $\mu$ and it is related to the variance of the distribution. Furthermore, $\omega_d = 2\pi^{d/2} [\Gamma(d/2)]^{-1}$ is the surface area of the unit sphere in $\mathbb{R}^d$ (see Golzy and Markatou, 2020).

-   **Clustering Algorithm for Spherical Data**: the package incorporates a unique clustering algorithm specifically tailored for $d$-dimensional spherical data and it is especially useful in the presence of noise in the data and the presence of non-negligible overlap between clusters. This algorithm leverages a mixture of Poisson kernel-based densities on the $d$-dimensional Sphere, enabling effective clustering of spherical data or data that has been spherically transformed. 

-   **Additional Features**: Alongside these functionalities, the software includes additional graphical functions, aiding users in validating and representing the cluster results as well as enhancing the interpretability and usability of the analysis.

For an introduction to the usage of `QuadratiK` see the vignette [Introduction to the QuadratiK Package](https://docs.ropensci.org/QuadratiK/articles/Introduction.html).

## Installation

You can install the version published on CRAN of `QuadratiK`

``` r
install.packages("QuadratiK")
```

or the development version from GitHub

``` r
library(devtools)
install_github('ropensci/QuadratiK')
# or via the rOpenSci organization repository
install.packages("QuadratiK", repos = "https://ropensci.r-universe.dev")
```

The `QuadratiK` package is also available in Python on PyPI <https://pypi.org/project/QuadratiK/> and also as a Dashboard application. Usage instruction for the Dashboard can be found at <https://quadratik.readthedocs.io/en/latest/user_guide/dashboard_application_usage.html>.

## Authors

Giovanni Saraceno, Marianthi Markatou, Raktim Mukhopadhyay, Mojgan Golzy\
Maintainer: Giovanni Saraceno \<[giovanni.saraceno\@unipd.it](mailto:giovanni.saraceno@unipd.it)\>


## Citation

If you use this package in your research or work, please cite it as follows:

Saraceno, G., Markatou, M., Mukhopadhyay, R. and Golzy, M. (2024). QuadratiK: Collection of Methods Constructed using Kernel-Based Quadratic Distances. <https://cran.r-project.org/package=QuadratiK>, <https://github.com/ropensci/QuadratiK>, <https://docs.ropensci.org/QuadratiK/>.

```         
@Manual{saraceno2024QuadratiK,  
   title = {QuadratiK: Collection of Methods Constructed using Kernel-Based 
            Quadratic Distances},  
   author = {Giovanni Saraceno and Marianthi Markatou and Raktim Mukhopadhyay 
             and Mojgan Golzy},  
   year = {2024},  
   note = {<https://cran.r-project.org/package=QuadratiK>,
            <https://github.com/ropensci/QuadratiK>,
            <https://docs.ropensci.org/QuadratiK/>},  
}
```

and the associated paper:

Saraceno, G., Markatou, M., Mukhopadhyay, R. and Golzy, M. (2024). Goodness-of-Fit and Clustering of Spherical Data: the QuadratiK package in R and Python. arXiv preprint [arXiv:2402.02290v2](https://arxiv.org/abs/2402.02290).

```         
@misc{saraceno2024package, 
      title={Goodness-of-Fit and Clustering of Spherical Data: the QuadratiK 
             package in R and Python},  
      author={Giovanni Saraceno and Marianthi Markatou and Raktim Mukhopadhyay 
              and Mojgan Golzy}, 
      year={2024}, 
      eprint={2402.02290}, 
      archivePrefix={arXiv}, 
      primaryClass={stat.CO}, url={https://arxiv.org/abs/2402.02290}
}
```

## References

-   Ding, Y., Markatou, M. and Saraceno, G. (2023). “Poisson Kernel-Based Tests for Uniformity on the d-Dimensional Sphere.” Statistica Sinica. doi: [10.5705/ss.202022.0347](https://doi.org/10.5705/ss.202022.0347).

-   Golzy, M. & Markatou, M. (2020) Poisson Kernel-Based Clustering on the Sphere: Convergence Properties, Identifiability, and a Method of Sampling, Journal of Computational and Graphical Statistics, 29:4, 758-770, DOI: [10.1080/10618600.2020.1740713](https://doi.org/10.1080/10618600.2020.1740713).

-   Markatou, M. and Saraceno, G. (2024). “A Unified Framework for Multivariate Two- and k-Sample Kernel-based Quadratic Distance Goodness-of-Fit Tests.” [arXiv:2407.16374](https://doi.org/10.48550/arXiv.2407.16374)

## Details

The work has been supported by Kaleida Health Foundation and National Science Foundation.

Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
