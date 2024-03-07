\name{QuadratiK-package}
\alias{QuadratiK-package}
\docType{package}
\title{
 Collection of methods constructed using kernel-based quadratic distances
}
\description{
  The QuadratiK package is implemented in R and Python, providing a comprehensive set of goodness-of-fit tests and a clustering technique using kernel-based quadratic distances. This framework aims to bridge the gap between the statistical and machine learning literature. It includes: 
  \itemize{
      \item \bold{Goodness-of-Fit Tests}: The software implements one, two, and k-sample tests for goodness of fit, offering an efficient and mathematically sound way to assess the fit of probability distributions. Expanded capabilities include supporting tests for uniformity on the \eqn{d}-dimensional sphere based on Poisson kernel densities.
      \item \bold{Clustering Algorithm for Spherical Data}: the package incorporates a unique clustering algorithm specifically tailored for spherical data. This algorithm leverages a mixture of Poisson kernel-based densities on the sphere, enabling effective clustering of spherical data or data that has been spherically transformed. The package also provides the functions for density evaluation and random sampling from the Poisson kernel-based distribution.
      \item \bold{Additional Features}: Alongside these functionalities, the software includes additional graphical functions, aiding users in validating and representing the cluster results as well as enhancing the interpretability and usability of the analysis. 
   }
}
\details{ 
    The work has been supported by Kaleida Health Foundation, National Science Foundation and Department of Biostatistics, Univerisity at Buffalo. 
}
\author{
Giovanni Saraceno, Marianthi Markatou, Raktim Mukhopadhyay, Mojgan Golzy

Mantainer: Giovanni Saraceno <gsaracen@buffalo.edu>
}
\references{
Saraceno, G., Markatou, M., Mukhopadhyay, R., Golzy, M. (2023). Goodness-of-fit and clustering of spherical data: The QuadratiK package in R and Python. Technical Report, Department of Biostatistics, University at Buffalo.

Ding, Y., Markatou, M., Saraceno, G. (2023). “Poisson Kernel-Based Tests for Uniformity on the d-Dimensional Sphere.” Statistica Sinica. doi: doi:10.5705/ss.202022.0347.

Golzy, M., Markatou, M. (2020) Poisson Kernel-Based Clustering on the Sphere: Convergence Properties, Identifiability, and a Method of Sampling, Journal of Computational and Graphical Statistics, 29:4, 758-770, DOI: 10.1080/10618600.2020.1740713.

Markatou, M., Saraceno, G., Chen, Y. (2023). “A unified framework for multivariate two-sample and k-sample kernel-based quadratic distance goodness-of-fit tests” Manuscript, (Department of Biostatistics, University at Buffalo).
}