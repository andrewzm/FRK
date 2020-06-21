#include <TMB.hpp>
#include <cmath>
#include "FRK-init.h"

// forward decleration of functions defined after objective_function template:
bool isCanonicalLink(std::string response, std::string link);

// FIXME: 
//  - Forward declerations for templates
//  - epsilon as an argument rather than created in the function


// Logarithm of the diagonal entries of a sparse matrix, and compute their sum.
template<class Type>
Type diaglnSum(Eigen::SparseMatrix<Type> mat){
  Type x = mat.diagonal().array().log().sum();
  return x;
}


// Computes the mean (i.e., applies the inverse-link function)
template<class Type>
Type inverseLinkFunction(Type Y_O, std::string link) {
  double epsilon{10e-8};
  Type mu_O;
  if (link == "identity")         mu_O = Y_O;
  if (link == "inverse")          mu_O = 1.0 / Y_O;
  if (link == "inverse-squared")  mu_O = 1.0 / sqrt(Y_O);
  if (link == "log")              mu_O = exp(Y_O) + epsilon;
  if (link == "square-root")      mu_O = Y_O * Y_O + epsilon;
  if (link == "logit")            mu_O = 1.0 / (1.0 + exp(-1.0 * Y_O));
  if (link == "probit")           mu_O = pnorm(Y_O);
  if (link == "cloglog")          mu_O = 1.0 - exp(-exp(Y_O));
  return mu_O;
}

// Canonical parameter, lambda, as a function of the mean, mu
template<class Type>
Type canonicalParameter(Type mu_O, Type k_Z, std::string response) {
  double epsilon{10e-8};
  Type lambda;
  if (response == "gaussian") { 
    lambda = mu_O;
  } else if (response == "gamma") {
    lambda = 1.0 / (mu_O + epsilon);
  } else if (response == "inverse-gaussian") {
    lambda = 1.0 / (mu_O * mu_O + epsilon);
  } else if (response == "poisson") {
    lambda  =   log(mu_O + epsilon);
  } else if (response == "negative-binomial") {
    lambda = -log(1.0 + k_Z/(mu_O + epsilon));
  } else if (response == "binomial") {
    lambda = log((mu_O + epsilon) / (k_Z - mu_O + epsilon));
  } else if (response == "bernoulli") {
    lambda = log((mu_O + epsilon) / (1.0 - mu_O + epsilon));
  }
  return lambda;
}

// FIXME: Change this function to if-elseif statements (more efficient, as it can break the chain)
template<class Type>
Type cumulantFunction(Type x, Type k_Z, std::string response, std::string parameterisation) {
  
  // NB: if the parameterisation is lambda, the canonical parameter should be passed in for x.
  // if the parameterisation is lambda, the mean should be passed in for x.
  double epsilon{10e-8};
  Type b;
  if (parameterisation == "lambda") {
    if (response == "gaussian")           b = (x * x) / 2.0;
    if (response == "gamma")              b =  log(x);
    if (response == "inverse-gaussian")   b =  2.0 * sqrt(x);
    if (response == "poisson")            b =  exp(x);
    if (response == "bernoulli")          b =  log(1.0 + exp(x));
  } else if (parameterisation == "mu") {
    if (response == "gaussian") b = (x * x) / 2.0;
    if (response == "gamma") b =  -log(x + epsilon);
    if (response == "inverse-gaussian") b =  2.0 / (x + epsilon);
    if (response == "poisson") b =  x;
    if (response == "negative-binomial") b = k_Z * log(1 + x / k_Z);
    if (response == "binomial") b = -k_Z * log(1.0 - (x - epsilon) / k_Z);
    if (response == "bernoulli") b =  -log(1.0 - x + epsilon);
  }
  return b;
}



template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // To do:
  // - I want to change sigma2 to kappa
  // - Change the latticeKrig formulation to match how I describe it in the report
  // - Add spatio-temporal separable K_type
  // - Change r_si to r_sk

  
  // typedef's:
  typedef Eigen::SparseMatrix<Type> SpMat;    // Sparse matrices called 'SpMat'
  typedef Eigen::Triplet<Type> T;             // Triplet lists called 'T'
  
  // ---- 0. Data and Parameters ----
  
  DATA_VECTOR(Z);             // Vector of observations
  int m    = Z.size();        // Sample size
  DATA_MATRIX(X);             // Design matrix of fixed effects
  DATA_SPARSE_MATRIX(S);      // Design matrix for basis function random weights
  DATA_SCALAR(sigma2e);       // Measurement error estimate (relevant to Gaussian case only)
  DATA_STRING(K_type);        // Indicates the desired model formulation of eta prior (K or Q)
  DATA_STRING(response);      // String specifying the response distribution
  DATA_STRING(link);          // String specifying the link function
  DATA_VECTOR(k_Z);           // "Known-constant" parameter (only relevant for negative-binomial and binomial)
  DATA_INTEGER(temporal);     // Boolean indicating whether we are in space-time or not
  
  DATA_INTEGER(r_t);         // Total number of temporal basis functions
  DATA_IVECTOR(r_si);        // Vector containing the number of spatial basis functions at each resolution
  int nres = r_si.size();    // Number of resolutions (of spatial basis functions) 
  int r_s  = r_si.sum();     // Total number of spatial basis functions
  int r = r_s * r_t;
  
  DATA_VECTOR(alpha);         // Tapering parameters (only relevant for block-exponential formulation)
  DATA_IVECTOR(row_indices);
  DATA_IVECTOR(col_indices);
  DATA_VECTOR(x);             
  DATA_IVECTOR(nnz);          // Integer vector indicating the number of non-zeros at each resolution of K_tap or Q
  DATA_IVECTOR(n_r);          // Integer vector indicating the number of rows at each resolution (applicable only if K-type == separable)
  DATA_IVECTOR(n_c);          // Integer vector indicating the number of columns at each resolution (applicable only if K-type == separable)
  
  // Parameters/basis function variance components/latent random effects 
  PARAMETER_VECTOR(beta);
  PARAMETER(logsigma2xi);     Type sigma2xi = exp(logsigma2xi);
  PARAMETER(logphi);          Type phi      = exp(logphi);
  
  // Variance components
  PARAMETER_VECTOR(logsigma2);  vector<Type> sigma2 = exp(logsigma2);
  PARAMETER_VECTOR(logtau);     vector<Type> tau    = exp(logtau);
  PARAMETER(logsigma2_t);       Type sigma2_t       = exp(logsigma2_t);
  PARAMETER(logrho_t);          Type rho_t          = exp(logrho_t);
  PARAMETER_VECTOR(logdelta);   vector<Type> delta  = exp(logdelta);
  
  // Latent random effects (will be integrated out)
  PARAMETER_VECTOR(eta);
  PARAMETER_VECTOR(xi_O);
  
  // Small, positive constant used to avoid division and logarithm of zero:
  Type epsilon = 10.0e-8;
  
  // ---- 1. Construct prior covariance matrix K or precision matrix Q  ---- //
  
  Type logdetQ_inv = 0; 
  Type quadform_eta = 0; 
  
  
  // ---- Temporal precision matrix (if relevant) ----
  
  // NB: this code assumes the temporal basis functions are at one resolution only. 
  
  // In C++ variables declared in if-statements or for-loops cannot be accessed 
  // outside the scope of the if-statement or for-loop in which they were 
  // declared. Hence we need to define a matrix variable here, even though this 
  // variable will not be used if temporal == 0. This matrix, J, will store 
  // the various products of L_sk, H_k, and L_t needed for the quadratic form. 
  matrix<Type> J; 
  
  if (temporal == 1) {
    // Construct the temporal Cholesky factor, L_t.
    std::vector<T> tripletList_L_t;
    tripletList_L_t.reserve(2 * r_t - 1);
    Type common_term = 1.0 / sqrt(sigma2_t * (1.0 - rho_t * rho_t));
    
    // Diagonal entries (except last diagonal entry), lower diagonal entries
    for (int j = 0; j < (r_t - 1); j++) {
      tripletList_L_t.push_back(T(j, j, 1));
      tripletList_L_t.push_back(T(j + 1, j, -rho_t * common_term));
    }
    // Final diagonal entry
    tripletList_L_t.push_back(T(r_t - 1, r_t - 1, 1.0 / sqrt(sigma2_t)));
    
    // Convert triplet list of non-zero entries to a true SparseMatrix.
    SpMat L_t(r_t, r_t);
    L_t.setFromTriplets(tripletList_L_t.begin(), tripletList_L_t.end());
    
    // Log-determinant
    logdetQ_inv -= 2.0 * r_s * L_t.diagonal().array().log().sum();
    
    // Quadratic form
    J = eta;
    J.resize(r_s, r_t);
    J *= L_t; 
  }
  
  
  // ---- Spatial variance/precision matrix ----
  
  int start_x = 0;    // Keep track of starting point in x (the vector of non-zero coefficients)
  int start_eta = 0;  // Keep track of starting point in eta
  Type coef = 0;      // Variable to store the current element (coefficient) of the matrix
  
  for (int k = 0; k < nres; k++) { // For each resolution of spatial basis functions
    
    // Separable is very different, so treat it as a special case
    if (K_type == "separable") {
      
      // Variance components associated with row and column direction
      // The first half of the sigma2 vector corresponds to row direction variances
      // and second half the column variances. Same for rho.
      // FIX: Check that tail() gives the correct order (I want n-3, n-2, n-1, n, NOT n, n-1, n-2, n-3)
      vector<Type> sigma2_r = sigma2.head(nres);
      vector<Type> sigma2_c = sigma2.tail(nres);
      vector<Type> rho_r    = tau.head(nres);
      vector<Type> rho_c    = tau.tail(nres);
      
      // First construct M_r and M_c
      // Look into this for possible better alternative: https://eigen.tuxfamily.org/dox/group__TutorialAdvancedInitialization.html
      std::vector<T> tripletList_M_r;
      std::vector<T> tripletList_M_c;
      tripletList_M_r.reserve(2 * n_r[k] - 1);
      tripletList_M_c.reserve(2 * n_c[k] - 1);
      Type common_term_c = 1.0 / sqrt(sigma2_c[k] * (1.0 - rho_c[k] * rho_c[k]));
      Type common_term_r = 1.0 / sqrt(sigma2_r[k] * (1.0 - rho_r[k] * rho_r[k]));
      
      // Diagonal entries (except last diagonal entry) AND lower diagonal entries
      for (int j = 0; j < (n_r[k] - 1); j++) {
        tripletList_M_r.push_back(T(j, j, 1));
        tripletList_M_r.push_back(T(j + 1, j, -rho_r[k] * common_term_r));
      }
      for (int j = 0; j < (n_c[k] - 1); j++) {
        tripletList_M_c.push_back(T(j, j, 1));
        tripletList_M_c.push_back(T(j + 1, j, -rho_c[k] * common_term_c));
      }
      // Last diagonal entry
      tripletList_M_r.push_back(T(n_r[k] - 1, n_r[k] - 1, 1.0 / sqrt(sigma2_r[k])));
      tripletList_M_c.push_back(T(n_c[k] - 1, n_c[k] - 1, 1.0 / sqrt(sigma2_c[k])));
      
      // Convert triplet list of non-zero entries to a true SparseMatrix.
      SpMat M_r(n_r[k], n_r[k]);
      SpMat M_c(n_c[k], n_c[k]);
      M_r.setFromTriplets(tripletList_M_r.begin(), tripletList_M_r.end());
      M_c.setFromTriplets(tripletList_M_c.begin(), tripletList_M_c.end());
      
      // Log-determinant
      logdetQ_inv += -2.0 * n_c[k] * M_r.diagonal().array().log().sum();
      logdetQ_inv += -2.0 * n_r[k] * M_c.diagonal().array().log().sum();
      
      // Quadratic form
      matrix<Type> Hi = eta.segment(start_eta, r_si[k]);
      Hi.resize(n_c[k], n_r[k]);
      matrix<Type> vi = M_c.transpose() * Hi * M_r;
      vi.resize(r_si[k], 1); // vec() operator
      quadform_eta += (vi.array() * vi.array()).sum();
      
      
      
    } else {
      
      // Construct kth block as a sparse matrix: use triplet list (row, column, value)
      std::vector<T> tripletList;    // Create a vector of triplet lists, called 'tripletList'
      tripletList.reserve(nnz[k]);   // Reserve number of non-zeros in the matrix
      
      // Vector to store the row sums
      vector<Type> rowSums(r_si[k]);
      rowSums.fill(0);
      
      // make a quantity which is between 0 and 1 (for the correlation parameters)
      // vector<Type> rho = 1 / (1 + exp(-1 * logtau));
      
      // Compute the matrix coefficients and store them in the triplet list.
      if (K_type == "precision_exp") {
        
        for (int j = start_x; j < start_x + nnz[k]; j++){ // For each non-zero entry within resolution k
          if (col_indices[j] != row_indices[j]) {
            coef = -sigma2[k] * exp( -x[j] / tau[k] ) * pow(1.0 - x[j] / alpha[k], 2.0) * ( 1.0 + x[j] / (2.0 * alpha[k]));
            tripletList.push_back(T(row_indices[j] - start_eta, col_indices[j] - start_eta, coef));
            rowSums[row_indices[j] - start_eta] += coef;
          }
        }
        // Add the diagonal elements (these depend on the row sums)
        for (int j = 0; j < r_si[k]; j++) {
          tripletList.push_back(T(j, j, delta[k] - rowSums[j]));
        }
      }
      
      if (K_type == "neighbour") {
        for (int j = start_x; j < start_x + nnz[k]; j++){  // For each non-zero entry within resolution k
          if (row_indices[j] == col_indices[j]) {
            coef = tau[k] * (x[j] + sigma2[k]);
          } else {
            coef = tau[k] * x[j]; // The "neighbour matrix" in R is responsible for setting the weights
          }
          tripletList.push_back(T(row_indices[j] - start_eta, col_indices[j] - start_eta, coef));
        }
      }
      
      bool rhoInB{true};
      
      if (K_type == "latticekrig" && rhoInB) {
        for (int j = start_x; j < start_x + nnz[k]; j++){  // For each non-zero entry within resolution k
          if (row_indices[j] == col_indices[j]) {
            coef = (x[j] + sigma2[k]) / sqrt(tau[k] + epsilon);
          } else {
            coef = x[j] / sqrt(tau[k] + epsilon); // The "neighbour matrix" in R is responsible for setting the weights
          }
          tripletList.push_back(T(row_indices[j] - start_eta, col_indices[j] - start_eta, coef));
        }
      }
      
      if (K_type == "latticekrig" && !rhoInB) {
        for (int j = start_x; j < start_x + nnz[k]; j++){  // For each non-zero entry within resolution k
          if (row_indices[j] == col_indices[j]) {
            coef = x[j] + sigma2[k];
          } else {
            coef = x[j];
          }
          tripletList.push_back(T(row_indices[j] - start_eta, col_indices[j] - start_eta, coef));
        }
      }
      
      if (K_type == "block-exponential") {
        
        for (int j = start_x; j < start_x + nnz[k]; j++){   // For each non-zero entry within resolution i
          coef = sigma2[k] * exp( -x[j] / tau[k] ) * pow( 1.0 - x[j] / alpha[k], 2.0) * ( 1.0 + x[j] / (2.0 * alpha[k]));
          tripletList.push_back(T(row_indices[j] - start_eta, col_indices[j] - start_eta, coef));
        }
      }
      
      // Convert triplet list of non-zero entries to a true SparseMatrix.
      // NB: this is "Kk" the variance matrix if K_type == "block-exponential",
      // the "Bk" matrix if K_type == "latticekrig" (and rhoInB == false),
      // and the precision matrix "Qk" for all other formulations.
      SpMat mat(r_si[k], r_si[k]);
      mat.setFromTriplets(tripletList.begin(), tripletList.end());
      
      bool constructQ{false};
      if (K_type == "latticekrig" && !constructQ) {
        rhoInB ? mat = mat * mat :  mat = mat * mat / (tau[k] + epsilon);
      }
      
      
      // Compute the (upper) Cholesky factor of mat
      Eigen::SimplicialLLT< SpMat, Eigen::Upper > llt;
      llt.compute(mat);
      SpMat Uk = llt.matrixU();
      
      // Log-determinant
      if (K_type == "latticekrig" && !constructQ && rhoInB) {
        logdetQ_inv += -4.0 * r_t * Uk.diagonal().array().log().sum();
      } else if (K_type == "latticekrig" && !constructQ && !rhoInB) {
        logdetQ_inv += r_si[k] * log(tau[k] + epsilon) - 4.0 * r_t * Uk.diagonal().array().log().sum();
      } else if (K_type == "block-exponential") {
        logdetQ_inv += 2.0 * r_t * Uk.diagonal().array().log().sum();
      } else {
        logdetQ_inv += -2.0 * r_t * Uk.diagonal().array().log().sum();
      }
      
      // P (the permutation matrix)
      Eigen::PermutationMatrix<Eigen::Dynamic> P = llt.permutationP();
      
      // Construct the matrix Mk such that Qk = Mk' Mk.
      SpMat Mk(r_si[k], r_si[k]);
      if (K_type == "latticekrig" && !constructQ && rhoInB) {
        Mk = mat;
      } else if (K_type == "latticekrig" && !constructQ && !rhoInB) {
        Mk = mat / sqrt(tau[k] + epsilon);
      } else if (K_type == "block-exponential") {
        // Don't construct Mk explicitly with block-exponential
      } else {
        Mk = Uk * P;
      }
      
      // Quadratic form
      if (temporal == 1) {
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
  
  
  // ---- Quadratic form in the case of space-time ----
  if (temporal == 1) { 
    J.resize(r_s * r_t, 1); // apply the vec operator
    quadform_eta += (J.array() * J.array()).sum();
  }
  
  
  // -------- 3. Construct ln[eta|Q] and ln[xi|sigma2xi].  -------- //
  
  Type quadform_xi = pow(sigma2xi, -1.0) * (xi_O * xi_O).sum();

  Type ld_eta =  -0.5 * r * log(2.0 * M_PI) - 0.5 * logdetQ_inv - 0.5 * quadform_eta;
  Type ld_xi  =  -0.5 * m * log(2.0 * M_PI) - 0.5 * m * log(sigma2xi) - 0.5 * quadform_xi;
  
  
  // -------- 4. Construct ln[Z|Y_O]  -------- //
  
  // 4.1. Construct Y_O, the latent spatial process at observed locations
  vector<Type> Y_O  = X * beta + S * eta + xi_O;
  
  // 4.2 Compute the canonical parameter and cumulant function
  vector<Type> lambda(m);
  vector<Type> blambda(m);
  
  if (isCanonicalLink(response, link)){
    // Compute canonical parameter (simply equal to Y when g is canonical)
    lambda = Y_O;
    
    // Compute the cumulant function using the canonical parameter
    blambda = cumulantFunction(lambda, k_Z, response, "lambda");
    
  } else {
    
    // Compute the mean (or probability parameter if a logit, probit, or cloglog link is used)
    vector<Type> mu_O(m);
    mu_O = inverseLinkFunction(Y_O, link);
    
    // If response is negative-binomial or binomial, link the probability parameter to the mean:
    if (link == "logit" || link == "probit" || link == "cloglog") {
      if (response == "negative-binomial"){
        mu_O = k_Z * (1.0 / (mu_O + epsilon) - 1);
      } else if (response == "binomial") {
        mu_O *= k_Z; 
      }
    } 
    
    // Compute the canonical parameter and cumulant function using the mean
    lambda = canonicalParameter(mu_O, k_Z, response); 
    blambda = cumulantFunction(mu_O, k_Z, response, "mu");
  }
  
  // 4.3. Construct a(phi) and c(Z, phi).
  // NB: computation of C(Z, phi) for one-parameter exponential family is done within R
  Type aphi{1.0}; // initialise to 1.0, only change for two-parameter exponential families
  vector<Type> cZphi(m);
  
  if (response == "poisson" || response == "negative-binomial" ||
      response == "binomial" || response == "bernoulli") {
    phi = 1.0;
    cZphi = 0.0;
  } else if (response == "gaussian") {
    phi = sigma2e;
    aphi = phi;
    cZphi = -0.5 * (Z * Z / phi + log(2.0 * M_PI * phi));
  } else if (response == "gamma") {
    aphi    =   -phi;
    cZphi   =   log(Z/phi)/phi - log(Z) - lgamma(1.0/phi);
  } else if (response == "inverse-gaussian") {
    aphi    =   - 2.0 * phi;
    cZphi   =   - 0.5 / (phi * Z) - 0.5 * log(2.0 * M_PI * phi * Z * Z * Z);
  } 

  // 4.4. Construct ln[Z|Y_O]
  Type ld_Z  =  ((Z * lambda - blambda)/aphi).sum() + cZphi.sum();
  
  
  // -------- 5. Define Objective function -------- //
  
  // ln[Z, eta, xi | ...] =  ln[Z|Y_O] + ln[eta|K] + ln[xi|sigma2xi]
  // Specify the negative joint log-likelihood function,
  // as R optimisation routines minimise by default.
  
  Type nld = -(ld_Z  + ld_xi + ld_eta);
  
  return nld;
}


// Some notes:

// PARAMETER_VECTOR(eta) defines eta as an ARRAY object in Eigen i.e. vector<type>
// is equivalent to an array.  - be careful,  as matrices/vectors cannot be mixed 
// with arrays! 
// PARAMETER_MATRIX() casts as matrix<type>, which is indeed a matrix.

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


// // Function to construct the (lower) Cholesky factor of an AR1 precision matrix
// template<class Type>
// Type choleskyAR1(Type sigma2, Type rho, int n){
// 
//   // typedef's:
//   typedef Eigen::SparseMatrix<Type> SpMat;    // Sparse matrices called 'SpMat'
//   typedef Eigen::Triplet<Type> T;             // Triplet lists called 'T'
// 
//   std::vector< T > tripletList;
//   tripletList.reserve(2 * n - 1);
//   Type common_term = 1 / sqrt(sigma2 * (1 - rho * rho));
// 
//   // Diagonal entries (except last diagonal entry), lower diagonal entries
//   for (int j = 0; j < (n - 1); j++) {
//     tripletList.push_back(T(j, j, 1));
//     tripletList.push_back(T(j + 1, j, -rho * common_term));
//   }
//   // Final diagonal entry
//   tripletList.push_back(T(n - 1, n - 1, 1 / sqrt(sigma2)));
// 
//   // Convert triplet list of non-zero entries to a true SparseMatrix.
//   SpMat L(n, n);
//   L.setFromTriplets(tripletList.begin(), tripletList.end());
// 
//   return L;
// }



