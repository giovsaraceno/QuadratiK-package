[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/QuadratiK)](https://cran.r-project.org/package=QuadratiK)

# Collection of Methods Constructed using the Kernel-Based Quadratic Distances

\code{QuadratiK} provides the first implementation, in R and Python, of a comprehensive set of goodness-of-fit tests and a clustering technique for spherical data using kernel-based quadratic distances. This framework aims to bridge the gap between the statistical and machine learning literature. It includes:

- **Goodness-of-Fit Tests**: The software implements one, two, and k-sample tests for goodness of fit, offering an efficient and mathematically sound way to assess the fit of probability distributions. Expanded capabilities include supporting tests for uniformity on the *d*-dimensional Sphere based on Poisson kernel densities.

- **Clustering Algorithm for Spherical Data**: the package incorporates a unique clustering algorithm specifically tailored for spherical data. This algorithm leverages a mixture of Poisson kernel-based densities on the Sphere, enabling effective clustering of spherical data or data that has been spherically transformed. The package also provides the functions for density evaluation and random sampling from the Poisson kernel-based distribution.

- **Additional Features**: Alongside these functionalities, the software includes additional graphical functions, aiding users in validating and representing the cluster results as well as enhancing the interpretability and usability of the analysis.

## Details

The work has been supported by Kaleida Health Foundation, National Science Foundation and Department of Biostatistics, University at Buffalo.

## Authors

Giovanni Saraceno, Marianthi Markatou, Raktim Mukhopadhyay, Mojgan Golzy  
Email: [gsaracen@buffalo.edu](mailto:gsaracen@buffalo.edu)

## References

- Saraceno Giovanni, Markatou Marianthi, Mukhopadhyay Raktim, Golzy Mojgan 
(2024). Goodness-of-Fit and Clustering of Spherical Data: the QuadratiK 
package in R and Python. arXiv preprint [arXiv:2402.02290](https://arxiv.org/abs/2402.02290).

- Ding Yuxin, Markatou Marianthi, Saraceno Giovanni (2023). “Poisson Kernel-Based Tests for Uniformity on the d-Dimensional Sphere.” Statistica Sinica. doi: [10.5705/ss.202022.0347](https://doi.org/10.5705/ss.202022.0347).

- Mojgan Golzy & Marianthi Markatou (2020) Poisson Kernel-Based Clustering on the Sphere: Convergence Properties, Identifiability, and a Method of Sampling, Journal of Computational and Graphical Statistics, 29:4, 758-770, DOI: [10.1080/10618600.2020.1740713](https://doi.org/10.1080/10618600.2020.1740713).

- Markatou M, Saraceno G, Chen Y (2023). “Two- and k-Sample Tests Based on Quadratic Distances.” Manuscript, Department of Biostatistics, University at Buffalo.
