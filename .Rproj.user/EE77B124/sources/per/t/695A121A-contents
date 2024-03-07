#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>
using namespace Eigen;
using namespace Rcpp;
using namespace std;
// [[Rcpp::depends(RcppEigen)]]

//'
//' Compute the Gaussian kernel matrix between two samples
//'
//' @param x_mat A matrix containing the observations of X
//' @param y_mat A matrix containing the observations of Y
//' @param H Covariance matrix of the Gaussian kernel
//'
//' @return Matrix of computed Gaussian kernel values 
//'
//' @useDynLib QuadratiK
//' @rdname computeKernelMatrix
//' @keywords internal
//' 
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd computeKernelMatrix(const Eigen::MatrixXd& x_mat, const Eigen::MatrixXd& y_mat, const Eigen::MatrixXd& H)
{
   int k = x_mat.cols();
   
   Eigen::MatrixXd H_inv_sqrt_x, H_inv_sqrt_y, x_mat_star, y_mat_star, x_mat_rowsum, y_mat_rowsum, Qd;
   Eigen::MatrixXd H_inv = H.inverse();
   double Hdet = pow(H.determinant(), -0.5);
   double Pi_const = pow(2 * M_PI, -0.5 * k);
   
   Eigen::ArrayXd Diag_H_inv = H_inv.diagonal();
   
   H_inv_sqrt_x = Diag_H_inv.sqrt().replicate(1, x_mat.rows()).transpose();
   H_inv_sqrt_y = Diag_H_inv.sqrt().replicate(1, y_mat.rows()).transpose();
   x_mat_star = x_mat.cwiseProduct(H_inv_sqrt_x);
   y_mat_star = y_mat.cwiseProduct(H_inv_sqrt_y);
   x_mat_rowsum = (x_mat_star.cwiseProduct(x_mat_star)).rowwise().sum();
   y_mat_rowsum = (y_mat_star.cwiseProduct(y_mat_star)).rowwise().sum();
   Qd = x_mat_rowsum.replicate(1, y_mat.rows()) - 2 * x_mat_star * (y_mat_star.transpose()) + (y_mat_rowsum.replicate(1, x_mat.rows())).transpose();
   
   Eigen::MatrixXd kmat_zz = Pi_const * Hdet * ((-0.5 * Qd).array().exp());
   
   return kmat_zz;
}
//' Compute the Poisson kernel matrix between observations in a sample.
//'
//' @param x_mat A matrix containing the observations.
//' @param rho concentration parameter of the Poisson kernel.
//'
//' @return Matrix of computed Poisson kernel values
//'
//' @useDynLib QuadratiK
//' @rdname computePoissonMatrix
//' @keywords internal
//' 
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd computePoissonMatrix(const Eigen::MatrixXd& x_mat, double rho)
{
   int n_x = x_mat.rows();
   int k = x_mat.cols();
   
   Eigen::MatrixXd x_mat_star = x_mat * x_mat.transpose();
   //Eigen::MatrixXd pmat = (1 - pow(rho, 2)) / pow((1 + pow(rho, 2) - 2 * rho * x_mat_star.array()).array(), k/2.0).array() - MatrixXd::Ones(n_x, n_x);
   
   Eigen::ArrayXXd pmat_array = (1 - pow(rho, 2)) / (1 + pow(rho, 2) - 2 * rho * x_mat_star.array()).pow(k/2.0) - Eigen::ArrayXXd::Ones(n_x, n_x);
   Eigen::MatrixXd pmat = pmat_array.matrix();
   
   return pmat;
}
//' Non-parametric centered kernel
//' 
//' Compute the values of the non-parametric centered kernel given the values of the 
//' kernel matrix. 
//'
//' @param kmat_zz A matrix containing the values of the kernel functions for each pair of observations
//' @param n_z The number of total observations
//'
//' @return Matrix of centered kernel
//'
//'
//'
//' @useDynLib QuadratiK
//' @rdname NonparamCentering
//' @keywords internal
//' 
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd NonparamCentering(const Eigen::MatrixXd& kmat_zz, int n_z)
{
   
   Eigen::MatrixXd k_center = kmat_zz.array() -
      kmat_zz.rowwise().sum().replicate(1, n_z).array() / n_z -
      kmat_zz.colwise().sum().replicate(n_z, 1).array() / n_z +
      kmat_zz.sum()/ (n_z * (n_z-1));
   
   return k_center;
}
//' Parametric centered kernel
//' 
//' Compute the Gaussian kernel centered with respect to a Normal distribution
//' with mean vector mu and covariance Sigma
//'
//' @param kmat_zz A matrix containing the values of the kernel functions for each pair of observations
//' @param z_mat Matrix of observations used for computing the kernel matrix
//' @param H Covariance matrix of the Normal kernel
//' @param mu_hat Mean of centering of Normal distribution
//' @param Sigma_hat Covariance of centering Normal distribution
//'
//' @return Matrix of centered kernel
//'
//' @useDynLib QuadratiK
//' @rdname ParamCentering
//' @keywords internal
//' 
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd ParamCentering(const Eigen::MatrixXd& kmat_zz, const Eigen::MatrixXd& z_mat,
                               const Eigen::MatrixXd& H, const Eigen::MatrixXd& mu_hat,
                               const Eigen::MatrixXd& Sigma_hat)
{
   int n_z = z_mat.rows();
   
   Eigen::MatrixXd k_center = kmat_zz - computeKernelMatrix(z_mat, mu_hat, H + Sigma_hat).replicate(1, n_z) -
      computeKernelMatrix(mu_hat, z_mat, H + Sigma_hat).replicate(n_z, 1) +
      computeKernelMatrix(mu_hat, mu_hat, H + 2*Sigma_hat)(0,0)*Eigen::MatrixXd::Constant(n_z, n_z, 1.0);
   
   return k_center;
}
//' Compute kernel-based quadratic distance two-sample test with Normal kernel
//'
//' @param x_mat A matrix containing observations from the first sample
//' @param y_mat A matrix containing observations from the second sample
//' @param h The bandwidth parameter for the kernel function
//' @param centeringType String indicating the method used for centering the normal kernel
//' @param mu_hat Mean vector for the reference distribution 
//' @param Sigma_hat Covariance matrix of the reference distribution
//'
//' @return A scalar value representing the test statistic
//'
//' @details 
//' \code{mu_hat} and \code{Sigma_hat} need to be provided even if they are used 
//' only in case of "Param" \code{centeringType}.
//' 
//'
//' @useDynLib QuadratiK
//' @rdname stat2sample
//' @keywords internal
//' 
//' @noRd
// [[Rcpp::export]]
double stat2sample(Eigen::MatrixXd& x_mat, Eigen::MatrixXd& y_mat, double h,
                   const Eigen::MatrixXd& mu_hat, const Eigen::MatrixXd& Sigma_hat, const std::string& centeringType = "Nonparam")
{
   int n_x = x_mat.rows();
   int k = x_mat.cols();
   int n_y = y_mat.rows();
   int n_z = n_x + n_y;
   
   
   Eigen::MatrixXd z_mat(n_z, k);
   z_mat << x_mat, y_mat;
   
   Eigen::MatrixXd H = pow(h, 2) * Eigen::MatrixXd::Identity(k, k);
   Eigen::MatrixXd kmat_zz = computeKernelMatrix(z_mat,z_mat, H);
   kmat_zz.diagonal().setZero();
   
   Eigen::MatrixXd k_center;
   
   if(centeringType == "Nonparam")
   {
      k_center = NonparamCentering(kmat_zz, n_z);
   }
   else if(centeringType == "Param")
   {
      k_center = ParamCentering( kmat_zz,  z_mat, H, mu_hat, Sigma_hat);
   }
   else
   {
      // Handle error for unknown centeringType or use a default centering
      stop("Unknown centering type.");
   }
   
   k_center.diagonal().setZero();
   
   double Test_NonPar = (k_center.block(0, 0, n_x, n_x).sum() / (n_x * (n_x - 1))) -
      2 * (k_center.block(0, n_x, n_x, n_y).sum() / (n_x * n_y)) +
      (k_center.block(n_x, n_x, n_y, n_y).sum() / (n_y * (n_y - 1)));
   return Test_NonPar;
   
}
//' Compute kernel-based quadratic distance test for Normality
//'
//' @param x_mat A matrix containing the observations.
//' @param h The bandwidth parameter for the Gaussian kernel function.
//' @param mu_hat Mean vector for the reference distribution (if available)
//' @param Sigma_hat Covariance matrix of the reference distribution (if available)
//' @param centeringType String indicating the method used for centering the normal kernel ('Param' or 'Nonparam').
//'
//' @return A scalar value representing the test statistic.
//'
//' @useDynLib QuadratiK
//' @rdname kbNormTest
//' @keywords internal
//' 
//' @noRd
// [[Rcpp::export]]
double kbNormTest(Eigen::MatrixXd x_mat, double h, const Eigen::MatrixXd& mu_hat, const Eigen::MatrixXd Sigma_hat, const std::string& centeringType = "Param")
{
   int n_x = x_mat.rows();
   int k = x_mat.cols();
   
   
   Eigen::MatrixXd H = pow(h, 2) * Eigen::MatrixXd::Identity(k, k);
   Eigen::MatrixXd kmat_zz = computeKernelMatrix(x_mat, x_mat, H);
   kmat_zz.diagonal().setZero();
   
   Eigen::MatrixXd k_center;
   
   if (centeringType == "Nonparam")
   {
      k_center = NonparamCentering(kmat_zz, n_x);
   }
   else if (centeringType == "Param")
   {
      k_center = ParamCentering( kmat_zz,  x_mat, H, mu_hat, Sigma_hat);
   }
   else
   {
      stop("Unknown centering type.");
   }
   
   k_center.diagonal().setZero();
   
   // Compute your normality test statistic here
   double Test_Normality = k_center.sum() / (n_x * (n_x - 1));
   return Test_Normality;
   
}
//' Poisson kernel-based test for Uniformity on the Sphere
//' 
//' Compute the Poisson kernel-based test for Uniformity given a sample of observations on the Sphere
//'
//' @param x_mat A matrix containing the observations.
//' @param rho Concentration parameter of the Poisson kernel.
//'
//' @return Vector with the values of the U-statistic and V-statistic
//'
//' @useDynLib QuadratiK
//' @rdname statPoissonUnif
//' @keywords internal
//' 
//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd statPoissonUnif(const Eigen::MatrixXd& x_mat, double rho)
{
   int n_x = x_mat.rows();
   
   Eigen::MatrixXd pmat = computePoissonMatrix(x_mat, rho);
   //Eigen::MatrixXd pmat = (1 - pow(rho, 2)) / pow((1 + pow(rho, 2) - 2 * rho * x_mat_star.array()).array(), k/2.0).array() - MatrixXd::Ones(n_x, n_x);
   
   
   double Vn = pmat.sum() / n_x;
   
   // For the lower triangular sum, you'd typically loop through the matrix
   double lower_tri_sum = 0.0;
   for (int i = 0; i < n_x; ++i) {
      for (int j = 0; j < i; ++j) {  // Notice j < i for lower triangular
         lower_tri_sum += pmat(i, j);
      }
   }
   double Un = 2 * lower_tri_sum / (n_x * (n_x - 1));
   
   //double Var_Un = (2 / (n_x * (n_x - 1))) * ((1 + pow(rho, 2)) / pow(1 - pow(rho, 2), k - 1) - 1);
   
   //Un = Un / sqrt(Var_Un);
   
   Eigen::VectorXd results(2);
   results(0) = Un;
   results(1) = Vn;
   
   return results;
}
//' Kernel-based quadratic distance k-sample tests
//' 
//' Compute the kernel-based quadratic distance k-sample tests with the Normal kernel and bandwidth parameter h.
//'
//' @param x A matrix containing observations from the pooled sample, from the k samples.
//' @param y A vector containing observations' memberships to the k samples.
//' @param h The bandwidth parameter for the kernel function.
//' @param sizes Vector with sample sizes of the considered samples
//' @param cum_size Vector indicating the cumulative sizes, adding one sample at the time.
//'
//' @return A vector containing the two k-sample test statistics
//'
//' @useDynLib QuadratiK
//' @rdname stat_ksample_cpp
//' @keywords internal
//' 
//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd stat_ksample_cpp(const Eigen::MatrixXd& x, const Eigen::VectorXd& y,
                                 double h, const Eigen::VectorXd& sizes, const Eigen::VectorXd& cum_size) {
   int n = x.rows();
   int d = x.cols();
   
   int K = sizes.size();
   
   // Matrix H initialization
   Eigen::MatrixXd H = pow(h, 2) * Eigen::MatrixXd::Identity(d, d);
   
   Eigen::MatrixXd kmat_zz = computeKernelMatrix(x,x, H);
   kmat_zz.diagonal().setZero();
   
   Eigen::MatrixXd k_center = NonparamCentering(kmat_zz, n);
   
   
   double TraceK = 0;
   double Tn = 0;
   for (int l = 1; l <= K; ++l) {
      for (int r = l; r <= K; ++r) {
         if (l == r) {
            
            double sumVal=0;
            for (int i = cum_size[l - 1]; i < cum_size[l]; ++i) {
               for (int j = cum_size[r - 1]; j < cum_size[r]; ++j) {
                  if (i != j) {
                     sumVal += k_center(i, j);
                  }
               }
            }
            TraceK += sumVal / (sizes[l-1] * (sizes[l-1] - 1));
         } else {
            double meanVal = 0;
            for (int i = cum_size[l - 1]; i < cum_size[l]; ++i) {
               for (int j = cum_size[r - 1]; j < cum_size[r]; ++j) {
                  meanVal += k_center(i, j);
               }
            }
            Tn -= 2 * meanVal / (sizes[l-1] * sizes[r-1]);
         }
      }
   }
   
   Eigen::VectorXd result(2);
   result(0) = (K - 1) * TraceK + Tn;
   result(1) = TraceK + Tn / (K - 1);
   
   return result;
}

