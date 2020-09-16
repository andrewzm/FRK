#include <TMB.hpp>
#include <cmath>
#include "FRK-init.h"


// Transform correlation parameter to 
template <class Type>
Type transform_minus_one_to_one(Type x){
  double epsilon1{10e-8};
  double epsilon2{2 / (2 - epsilon1) - 1 + 10e-8};
  return Type(2)/(Type(1) + exp(-Type(2) * x) + epsilon2) - Type(1) + epsilon1;
}

// Forward declarations of functions defined after objective_function template
// (function definitions are after the objective_function template)

// Determine whether the link is canonical for this distribution
bool isCanonicalLink(std::string response, std::string link);

// Constructs the (lower) Cholesky factor of an AR1 precision matrix
template<class Type>
Eigen::SparseMatrix<Type> choleskyAR1(Type sigma2, Type rho, int n);

// Logarithm of the diagonal entries of a sparse matrix, and compute their sum.
template<class Type> 
Type diagLogSum(Eigen::SparseMatrix<Type> mat);

// Computes the mean given the latent process (i.e., applies the inverse-link function)
template<class Type>
Type inverseLinkFunction(Type Y_Z, std::string link);

// Canonical parameter, lambda, as a function of the mean, mu
template<class Type> 
Type canonicalParameter(Type mu_Z, Type k_Z, std::string response); 

// Computes the cumulant function
template<class Type>
Type cumulantFunction(Type x, Type k_Z, std::string response, std::string parameterisation);




template<class Type>
Type objective_function<Type>::operator() ()
{
  // typedef's:
  typedef Eigen::SparseMatrix<Type> SpMat;    
  typedef Eigen::Triplet<Type> T;             
  
  // ---- 0. Data and Parameters ----
  
  DATA_VECTOR(Z);             // Vector of observations
  int m = Z.size();           // Sample size
  DATA_MATRIX(X_O);           // Design matrix of fixed effects
  DATA_SPARSE_MATRIX(S_O);    // Design matrix for observed basis function random weights
  DATA_SPARSE_MATRIX(C_O);    // Incidence matrix, mapping the BAUs to the observations
  int mstar = C_O.cols();     // The number of observed BAUs
  DATA_STRING(K_type);        // Indicates the desired model formulation of eta prior (K or Q)
  DATA_STRING(response);      // String specifying the response distribution
  DATA_STRING(link);          // String specifying the link function
  DATA_VECTOR(k_BAU);         // Known size parameter at the BAU level (only relevant for negative-binomial and binomial)
  DATA_VECTOR(k_Z);           // Known size parameter at the DATA support level (only relevant for negative-binomial and binomial)
  DATA_INTEGER(temporal);     // Boolean indicating whether we are in space-time or not (1 if true, 0 if false)
  
  DATA_INTEGER(r_t);         // Total number of temporal basis functions
  DATA_IVECTOR(r_si);        // Vector containing the number of spatial basis functions at each resolution
  int nres = r_si.size();    // Number of resolutions (of spatial basis functions) 
  int r_s  = r_si.sum();     // Total number of spatial basis functions
  int r = r_s * r_t;
  DATA_IVECTOR(spatial_BAU_id);
  
  DATA_VECTOR(beta);          // Tapering parameters (only relevant for block-exponential formulation)
  DATA_IVECTOR(row_indices);
  DATA_IVECTOR(col_indices);
  DATA_VECTOR(x);             
  DATA_IVECTOR(nnz);          // Integer vector indicating the number of non-zeros at each resolution of K_tap or Q
  DATA_IVECTOR(n_r);          // Integer vector indicating the number of rows at each resolution (applicable only if K-type == separable)
  DATA_IVECTOR(n_c);          // Integer vector indicating the number of columns at each resolution (applicable only if K-type == separable)

  DATA_SCALAR(sigma2e);  // measurement error for Gaussian data model (fixed)
  
  DATA_VECTOR(sigma2fs_hat);   // estimate of sigma2fs (the fine-scale variance component)
  DATA_INTEGER(fix_sigma2fs);  // Flag indicating whether we should fix sigma2fs or not (1 if true, 0 if false)
  DATA_INTEGER(include_fs);    // Flag indicating whether fine-scale variation is included in the model (1 if true, 0 if false)
  
  // Fixed effects and variance components relating directly to data
  PARAMETER_VECTOR(alpha);
  PARAMETER(logphi);          Type phi      = exp(logphi);
  
  // Variance components relating to eta
  PARAMETER_VECTOR(logsigma2);  vector<Type> sigma2 = exp(logsigma2);
  PARAMETER_VECTOR(logtau);     // The transformation we apply to tau depends on which K_type is used
  PARAMETER(logsigma2_t);       Type sigma2_t       = exp(logsigma2_t);
  PARAMETER(frho_t);            Type rho_t          = transform_minus_one_to_one(frho_t);
  
  
  // Fine-scale variation variance parameter
  // If we are not estimating sigma2fs, fix it to the estimate. Otherwise, 
  // treat it as a parameter.
  PARAMETER_VECTOR(logsigma2fs);   vector<Type> sigma2fs = exp(logsigma2fs);
  DATA_INTEGER(fs_by_spatial_BAU);
  if (fix_sigma2fs)
    sigma2fs = sigma2fs_hat;
  
  
  // Latent random effects (will be integrated out).
  // Write it this way so that we have the option to exclude xi_O from within R.
  PARAMETER_VECTOR(random_effects);
  vector<Type> eta = random_effects.head(r);
  vector<Type> xi_O(mstar);
  
  if (include_fs) {
    xi_O = random_effects.tail(mstar);
  } else {
    xi_O.fill(0.0);
  }
  
  // Small, positive constant used to avoid division and logarithm of zero:
  Type epsilon = 10.0e-8;
  
  
  // ---- 1. Construct prior covariance matrix K or precision matrix Q  ---- //
  
  Type logdetQ_inv{0}; 
  Type quadform_eta{0}; 
  
  
  // 1.1. Temporal precision matrix (if relevant) 
  
  // NB: this code assumes the temporal basis functions are at one resolution only. 
  
  // In C++ variables declared in if-statements or for-loops cannot be accessed 
  // outside the scope of the if-statement or for-loop in which they were 
  // declared. Hence we need to define a matrix variable here, even though this 
  // variable will not be used if temporal == 0. This matrix, J, will store 
  // the various products of L_sk, H_k, and L_t needed for the quadratic form. 
  matrix<Type> J; 
  if (temporal) {
    // Construct the temporal Cholesky factor, L_t:
    Eigen::SparseMatrix<Type> L_t = choleskyAR1(sigma2_t, rho_t, r_t);
    
    // Log-determinant:
    logdetQ_inv -= 2.0 * r_s * diagLogSum(L_t);
    
    // Quadratic form:
    J = eta;
    J.resize(r_s, r_t);
    J *= L_t; 
  }
  
  
  // 1.2. Spatial variance/precision matrix 
  
  int start_x{0};    // Keep track of starting point in x (the vector of non-zero coefficients)
  int start_eta{0};  // Keep track of starting point in eta
  Type coef{0};      // Variable to store the current element (coefficient) of the matrix
  
  for (int k = 0; k < nres; k++) { // For each resolution of spatial basis functions
    
    // Separable is very different, so treat it as a special case
    if (K_type == "separable") {
      
      // NB: "separable" has only been implemented for spatial
      
      // Variance components associated with row and column direction
      // The first half of the sigma2 vector corresponds to row direction variances
      // and second half the column variances. Same for rho.
      vector<Type> sigma2_r = sigma2.head(nres);
      vector<Type> sigma2_c = sigma2.tail(nres);
      // The correlation parameter in an AR1 must be between -1 and 1. 
      vector<Type> tau      = transform_minus_one_to_one(logtau);
      vector<Type> rho_r    = tau.head(nres); 
      vector<Type> rho_c    = tau.tail(nres);
      
      // First construct M_r and M_c
      Eigen::SparseMatrix<Type> M_r = choleskyAR1(sigma2_r[k], rho_r[k], n_r[k]);
      Eigen::SparseMatrix<Type> M_c = choleskyAR1(sigma2_c[k], rho_c[k], n_c[k]);
      
      // Log-determinant
      logdetQ_inv += -2.0 * n_c[k] * diagLogSum(M_r);
      logdetQ_inv += -2.0 * n_r[k] * diagLogSum(M_c);
      
      
      // Quadratic form
      if (temporal) {
        // First compute L_{s, k}' H_k
        // FIXME: HAVENT IMPLEMENTED SPATIO-TEMPORAL WITH SEPARABLE YET
        // J.block(start_eta, 0, r_si[k], r_t) = Mk * J.block(start_eta, 0, r_si[k], r_t);
        
      } else {
        matrix<Type> Hi = eta.segment(start_eta, r_si[k]);
        Hi.resize(n_c[k], n_r[k]);
        matrix<Type> vi = M_c.transpose() * Hi * M_r;
        vi.resize(r_si[k], 1); // vec() operator
        quadform_eta += (vi.array() * vi.array()).sum();
      }
      
    } else {
      
      // For all K_type != "separable", we only require tau > 0, rather than -1 < tau < 1.
      vector<Type> tau = exp(logtau);
      
      // Construct kth block as a sparse matrix: use triplet list (row, column, value)
      std::vector<T> tripletList;    // Create a vector of triplet lists, called 'tripletList'
      tripletList.reserve(nnz[k]);   // Reserve number of non-zeros in the matrix
      
      
      // Compute the matrix coefficients and store them in the triplet list.
      // FIXME: I COULD DO THIS IN A FUNCTION
      
      if (K_type == "neighbour") {
        for (int j = start_x; j < start_x + nnz[k]; j++) {  
          (row_indices[j] == col_indices[j]) ? coef = tau[k] * (x[j] + sigma2[k]) : coef = tau[k] * x[j];
          tripletList.push_back(T(row_indices[j] - start_eta, col_indices[j] - start_eta, coef));
        }
      } else if (K_type == "block-exponential") {
        for (int j = start_x; j < start_x + nnz[k]; j++) {  
          coef = sigma2[k] * exp( -x[j] / tau[k] ) * pow( 1.0 - x[j] / beta[k], 2.0) * ( 1.0 + x[j] / (2.0 * beta[k]));
          tripletList.push_back(T(row_indices[j] - start_eta, col_indices[j] - start_eta, coef));
        }
      }
      
      
      // Convert triplet list of non-zero entries to a true SparseMatrix.
      // NB: mat is "Kk" the variance matrix if K_type == "block-exponential",
      // and the precision matrix "Qk" for all other formulations.
      SpMat mat(r_si[k], r_si[k]);
      mat.setFromTriplets(tripletList.begin(), tripletList.end());
      
      // Compute the (upper) Cholesky factor of mat
      Eigen::SimplicialLLT< SpMat, Eigen::Upper > llt;
      llt.compute(mat);
      SpMat Uk = llt.matrixU();
      
      // Log-determinant
      if (K_type == "block-exponential") {
        logdetQ_inv += 2.0 * r_t * diagLogSum(Uk);
      } else {
        logdetQ_inv += -2.0 * r_t * diagLogSum(Uk);
      }
      
      // P (the permutation matrix)
      Eigen::PermutationMatrix<Eigen::Dynamic> P = llt.permutationP();
      
      // Construct the matrix Mk such that Qk = Mk' Mk.
      SpMat Mk(r_si[k], r_si[k]);
      if (K_type != "block-exponential") { // We don't explicitly construct Mk with block-exponential
        Mk = Uk * P; 
      } 
      
      // Quadratic form
      if (temporal) {
        if (K_type == "block-exponential") {
          J.block(start_eta, 0, r_si[k], r_t) = (Uk.transpose()).template triangularView<Eigen::Lower>().solve(P * J.block(start_eta, 0, r_si[k], r_t));
        } else {
          J.block(start_eta, 0, r_si[k], r_t) = Mk * J.block(start_eta, 0, r_si[k], r_t);
        }
      } else {
        vector<Type> vk(r_si[k]);
        if (K_type == "block-exponential") {
          vk = ((Uk.transpose()).template triangularView<Eigen::Lower>().solve(P * eta.segment(start_eta, r_si[k]).matrix())).array();
        } else {
          vk = Mk * eta.segment(start_eta, r_si[k]).matrix();
        }
        quadform_eta += (vk * vk).sum();
      }
    }
    
    start_eta += r_si[k];
    start_x   += nnz[k];
  }
  
  
  // 1.3. Quadratic form in the case of space-time 
  if (temporal) { 
    J.resize(r_s * r_t, 1); // apply the vec operator
    quadform_eta += (J.array() * J.array()).sum();
  }
  
  
  // ---- 2. Construct ln[eta|Q] and ln[xi_O|sigma2fs].  ---- //
  
  Type ld_eta =  -0.5 * r * log(2.0 * M_PI) - 0.5 * logdetQ_inv - 0.5 * quadform_eta;
 
  Type ld_xi_O{0}; 
  if (include_fs) {
    
    ld_xi_O += -0.5 * mstar * log(2.0 * M_PI); // constant term in the log-density of xi_O
    
    if (fs_by_spatial_BAU) {
      
      vector<Type> sigma2fs_long(mstar); 
      for (int i = 0; i < mstar; i++) {
        sigma2fs_long[i] = sigma2fs[spatial_BAU_id[i]]; 
      }
      Type quadform_xi_O = ((xi_O * xi_O) / sigma2fs_long).sum();
      ld_xi_O += -0.5 * (sigma2fs_long).log().sum() - 0.5 * quadform_xi_O;

    } else {
      Type quadform_xi_O = pow(sigma2fs[0], -1.0) * (xi_O * xi_O).sum();
      ld_xi_O += - 0.5 * mstar * log(sigma2fs[0]) - 0.5 * quadform_xi_O;
    }
  }
    


  // ---- 2a. Prior for sigma^2_\xi, the fine-scale variance component ---- //
  
  // We assume an inverse-gamma prior distribution for the fine-scale variance, 
  // wherein the mean and variance of the prior is equal to sigmafs_hat.
  // Type a = sigma2fs_hat + 2.0;
  // Type b = pow(sigma2fs_hat, 2.0) + 1.0;
  // Type ld_sigma_xi = -(a + 1) * sigma2fs - b / sigma2fs;
  
  
  // ---- 3. Construct ln[Z|Y_Z]  ---- //
  
  // 3.1 Construct Y, the latent process at the BAU level, and mu 
  // (NB: mu is the probability parameter if a logit, probit, or cloglog link is used)
  vector<Type> Y_O  = X_O * alpha + S_O * eta + xi_O;
  vector<Type> mu_O = inverseLinkFunction(Y_O, link);
  
  // If response is negative-binomial or binomial, link the probability parameter to the mean:
  if (link == "logit" || link == "probit" || link == "cloglog") {
    if (response == "negative-binomial") {
      mu_O = k_BAU * (1.0 / (mu_O + epsilon) - 1);
    } else if (response == "binomial") {
      mu_O *= k_BAU; 
    }
  } 
  
  // Compute the mean over the observed data supports:
  vector<Type> mu_Z = C_O * mu_O;
  
  // Compute the canonical parameter and cumulant function using the mean
  vector<Type> lambda = canonicalParameter(mu_Z, k_Z, response); 
  vector<Type> blambda = cumulantFunction(mu_Z, k_Z, response, "mu");

  
  // 4.3. Construct a(phi) and c(Z, phi).
  // NB: computation of C(Z, phi) for one-parameter exponential family is done within R
  Type aphi{1.0}; // initialise to 1.0, only change for two-parameter exponential families
  vector<Type> cZphi(m);
  
  if (response == "poisson" || response == "negative-binomial" ||
      response == "binomial" || response == "bernoulli") 
    phi = 1.0;
  
  
  if (response == "gaussian") {
    phi = sigma2e;
    aphi = phi;
    cZphi = -0.5 * (Z * Z / phi + log(2.0 * M_PI * phi));
  } else if (response == "gamma") {
    aphi    =   -phi;
    cZphi   =   log(Z/phi)/phi - log(Z) - lgamma(1.0/phi);
  } else if (response == "inverse-gaussian") {
    aphi    =   - 2.0 * phi;
    cZphi   =   - 0.5 / (phi * Z) - 0.5 * log(2.0 * M_PI * phi * Z * Z * Z);
  } else if (response == "poisson") {
    for (int i = 0; i < m; i++) {
      cZphi[i] = -lfactorial(Z[i]);
    }
  } else if (response == "negative-binomial") {
    for (int i = 0; i < m; i++) {
      cZphi[i] = lfactorial(Z[i] + k_Z[i] - 1.0) - lfactorial(Z[i]) - lfactorial(k_Z[i] - 1.0);
    }
  } else if (response == "binomial") {
    for (int i = 0; i < m; i++) {
      cZphi[i] = lfactorial(k_Z[i]) - lfactorial(Z[i]) - lfactorial(k_Z[i] - Z[i]);
    }
  } else if (response == "bernoulli") {
    cZphi = 0.0;
  }
  

  // 4.4. Construct ln[Z|Y_Z]
  Type ld_Z  =  ((Z * lambda - blambda)/aphi).sum() + cZphi.sum();
  
  
  // -------- 5. Define Objective function -------- //
  
  // ln[Z, eta, xi_O | ...] =  ln[Z|Y_Z] + ln[eta|K] + ln[xi_O|sigma2fs]
  // Specify the negative joint log-likelihood function,
  // as R optimisation routines minimise by default.
  
  // Type nld = -(ld_Z  + ld_xi_O + ld_eta + ld_sigma_xi);
  Type nld = -(ld_Z  + ld_xi_O + ld_eta);
  
  
  return nld;
}


// Logarithm of the diagonal entries of a sparse matrix, and compute their sum.
template<class Type> 
Type diagLogSum(Eigen::SparseMatrix<Type> mat){
  Type x = mat.diagonal().array().log().sum();
  return x;
}


// Function to construct the (lower) Cholesky factor of an AR1 precision matrix
template<class Type>
Eigen::SparseMatrix<Type> choleskyAR1(Type sigma2, Type rho, int n){
  
  double epsilon{10e-8};
  sigma2 += epsilon;
  typedef Eigen::Triplet<Type> T;   // typedef: Triplet lists called 'T'
  
  std::vector< T > tripletList;
  tripletList.reserve(2 * n - 1);
  // NB: Must have -1 < rho < 1
  Type common_term = 1 / sqrt(sigma2 * (1 - rho * rho));
  
  // Diagonal entries (except last diagonal entry), lower diagonal entries
  for (int j = 0; j < (n - 1); j++) {
    tripletList.push_back(T(j, j, 1));
    tripletList.push_back(T(j + 1, j, -rho * common_term));
  }
  // Final diagonal entry
  tripletList.push_back(T(n - 1, n - 1, 1 / sqrt(sigma2)));
  
  // Convert triplet list of non-zero entries to a true SparseMatrix.
  Eigen::SparseMatrix<Type> L(n, n);
  L.setFromTriplets(tripletList.begin(), tripletList.end());
  
  return L;
}

// Check whether the link is canonical
bool isCanonicalLink(std::string response, std::string link) {
  bool canonical_link{false};
  if ((response == "gaussian"         && link == "identity") ||
      (response == "gamma"            && link == "inverse") ||
      (response == "inverse-gaussian" && link == "inverse-squared") ||
      (response == "poisson"          && link == "log") ||
      (response == "bernoulli"        && link == "logit")){
    canonical_link = true;
  }
  return canonical_link;
}

// Computes the mean (i.e., applies the inverse-link function)
template<class Type>
Type inverseLinkFunction(Type Y_Z, std::string link) {
  double epsilon{10e-8};
  Type mu_Z;
  if       (link == "identity"){         mu_Z = Y_Z;
  }else if (link == "inverse"){          mu_Z = 1.0 / Y_Z;
  }else if (link == "inverse-squared"){  mu_Z = 1.0 / sqrt(Y_Z);
  }else if (link == "log"){              mu_Z = exp(Y_Z) + epsilon;
  }else if (link == "square-root"){      mu_Z = Y_Z * Y_Z + epsilon;
  }else if (link == "logit"){            mu_Z = 1.0 / (1.0 + exp(-1.0 * Y_Z));
  }else if (link == "probit"){           mu_Z = pnorm(Y_Z);
  }else if (link == "cloglog"){          mu_Z = 1.0 - exp(-exp(Y_Z));
  }
  return mu_Z;
}

// Canonical parameter, lambda, as a function of the mean, mu
template<class Type>
Type canonicalParameter(Type mu_Z, Type k_Z, std::string response) {
  double epsilon{10e-8};
  Type lambda;
  if (response == "gaussian") { 
    lambda = mu_Z;
  } else if (response == "gamma") {
    lambda = 1.0 / (mu_Z + epsilon);
  } else if (response == "inverse-gaussian") {
    lambda = 1.0 / (mu_Z * mu_Z + epsilon);
  } else if (response == "poisson") {
    lambda  =   log(mu_Z + epsilon);
  } else if (response == "negative-binomial") {
    lambda = -log(1.0 + k_Z/(mu_Z + epsilon));
  } else if (response == "binomial") {
    lambda = log((mu_Z + epsilon) / (k_Z - mu_Z + epsilon));
  } else if (response == "bernoulli") {
    lambda = log((mu_Z + epsilon) / (1.0 - mu_Z + epsilon));
  }
  return lambda;
}


// Computes the cumulant function
template<class Type>
Type cumulantFunction(Type x, Type k_Z, std::string response, std::string parameterisation) {
  
  // NB: if the parameterisation is lambda, the canonical parameter should be passed in for x.
  // if the parameterisation is lambda, the mean should be passed in for x.
  double epsilon{10e-8};
  Type b;
  if (parameterisation == "lambda") {
    if       (response == "bernoulli") {        b =  log(1.0 + exp(x));
    }else if (response == "gaussian") {         b = (x * x) / 2.0;
    }else if (response == "gamma"){             b =  log(x);
    }else if (response == "inverse-gaussian"){  b =  2.0 * sqrt(x);
    }else if (response == "poisson"){           b =  exp(x);
    }
  } 
  
  if (parameterisation == "mu") {
    if       (response == "bernoulli"){           b =  -log(1.0 - x + epsilon);
    }else if (response == "gaussian"){            b = (x * x) / 2.0;
    }else if (response == "gamma"){               b =  -log(x + epsilon);
    }else if (response == "inverse-gaussian"){    b =  2.0 / (x + epsilon);
    }else if (response == "poisson"){             b =  x;
    }else if (response == "negative-binomial"){   b = k_Z * log(1 + x / k_Z);
    }else if (response == "binomial")             b = -k_Z * log(1.0 - (x - epsilon) / k_Z);
  }
  
  return b;
}


// PARAMETER_VECTOR(eta) defines eta as an ARRAY object in Eigen i.e. vector<type>
// is equivalent to an array.  - be careful,  as matrices/vectors cannot be mixed 
// with arrays! 
// PARAMETER_MATRIX() casts as matrix<type>, which is indeed a matrix.
