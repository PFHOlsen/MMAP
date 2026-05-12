#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Epsilon for uniformization
const double uni_epsilon = 0.0001;

//// Matrix functions //////////////////////////////////////////////////////////

// Matrix exponential function and related functions are from 
// https://github.com/martinbladt/matrixdist_1.0

double inf_norm(arma::mat A) {
  double value{0.0};
  
  for (int i{0}; i < A.n_rows; ++i) {
    double row_sum{0.0};
    for (int j{0}; j < A.n_cols; ++j) {
      row_sum += fabs(A(i,j));
    }
    value = std::max(value, row_sum);
  }
  return value;
}

// [[Rcpp::export]]
arma::mat matrix_exponential(arma::mat A) {
  const int q{6};
  
  arma::mat expm(A.n_rows,A.n_cols);
  
  double a_norm{inf_norm(A)};
  
  int ee{static_cast<int>(log2(a_norm)) + 1};
  
  int s{std::max(0, ee + 1)};
  
  double t{1.0 / pow(2.0, s)};
  
  arma::mat a2 = A * t;
  arma::mat x = a2;
  
  double c{0.5};
  
  expm.eye(size(A));
  expm = expm + (a2 * c);
  
  arma::mat d;
  d.eye(size(A));
  d = (d + a2 * (-c));
  
  int p{1};
  
  for (int k{2}; k <= q; ++k) {
    c = c * static_cast<double>(q - k + 1) / static_cast<double>(k * (2 * q - k + 1));
    x = (a2 * x);
    expm = (x * c) + expm;
    if (p) {
      d = (x * c) + d;
    }
    else {
      d = (x * (-c)) + d;
    }
    p = !p;
  }
  expm = inv(d) * expm;
  for (int k = 1; k <= s; ++k) {
    expm = expm * expm;
  }
  return(expm);
}


// Assign submatrix to i,j'th submatrix-location in matrix A
void assign_submatrix(arma::mat &A, 
                      const arma::mat& submatrix, 
                      int i, 
                      int j) {
  
  int block_size = submatrix.n_cols;  
  
  int row_start = (i - 1) * block_size;
  int col_start = (j - 1) * block_size;
  
  A.submat(row_start, col_start,
           row_start + block_size - 1,
           col_start + block_size - 1) = submatrix;
}

// Extract extract submatrix with dimension blocksize in i,j't location of 
// matrix A
arma::mat extract_submatrix(arma::mat A, 
                            int i, 
                            int j, 
                            int block_size) {
  
  int row_start = (i - 1) * block_size;
  int col_start = (j - 1) * block_size;
  
  return A.submat(row_start, col_start,
                  row_start + block_size - 1,
                  col_start + block_size - 1);
}

double max_diagonal(const arma::mat & A) {
  double maximum{A(0,0)};
  for (int i{0}; i < A.n_rows; ++i) {
    if (A(i,i) > maximum) {
      maximum = A(i,i);
    }
  }
  return maximum;
}

void vector_of_matrices(std::vector<arma::mat> & vect, const arma::mat & S, double a, int vect_size) {
  arma::mat I;
  I.eye(size(S));
  
  arma::mat P = I + S * (1.0 / a);
  
  vect.push_back(I);
  
  for (int k{1}; k <= vect_size; ++k) {
    vect.push_back((P * (1.0 / k) ) * vect[k - 1]);
  }
}


arma::mat m_exp_sum(double x, int n, const std::vector<arma::mat> & pow_vector, double a) {
  arma::mat res_mat = pow_vector[0];
  
  for (int i{1}; i <= n; ++i) {
    res_mat = res_mat + pow_vector[i] * exp(i * std::log(a * x));
  }
  res_mat = res_mat * exp(-a * x);
  
  return res_mat;
}

void pow2_matrix(int n , arma::mat & A) {
  arma::mat aux_mat(size(A));
  
  for (int i{1}; i <= n; ++i) {
    aux_mat = A * A;
    A = aux_mat;
  }
}

int find_n(double h, double lambda) {
  int n{0};
  double cum_prob{0.0};
  
  do {
    cum_prob += R::dpois(n, lambda, false);
    ++n;
  } while (cum_prob < 1.0 - h);
  
  return (n - 1);
}

// Compute vector of matrix exponentials by scaling and squaring with uniformization
std::vector<arma::mat> expA_vec(const std::vector<double>& x, const arma::mat& A){
  int n = x.size();
  int p = A.n_cols;
  
  double a = max_diagonal(A * (-1.0));
  
  int m{find_n(uni_epsilon, 1.0)};
  
  std::vector<arma::mat> expA(n, arma::mat(p, p));
  arma::mat expAx(p,p);
  
  std::vector<arma::mat> aux_vect;
  vector_of_matrices(aux_vect, A, a, m);
  
  double x_i;
  
  for (int i = 0; i < n; ++i){
    x_i = x[i];
    
    if (x_i * a <= 1.0) {
      expAx = m_exp_sum(x_i, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x_i) / std::log(2.0);
      ++n;
      
      expAx = m_exp_sum(x_i / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expAx);
    }
    expA[i] = expAx;
  }
  return expA;
}

// Compute Van Loan matrix of matrix with blocks A, B, C
arma::mat VanLoan(double t, 
                  const arma::mat& A, 
                  const arma::mat& B, 
                  const arma::mat& C) {
  
  int p = A.n_rows;

  arma::mat VL_full(p+p, p+p, arma::fill::zeros);
  assign_submatrix(VL_full, A, 1, 1);
  assign_submatrix(VL_full, B, 1, 2);
  assign_submatrix(VL_full, C, 2, 2);
  
  double a = max_diagonal(VL_full * (-1.0));
  
  int m{find_n(uni_epsilon, 1.0)};
  
  arma::mat expVL;
  
  std::vector<arma::mat> aux_vect;
  vector_of_matrices(aux_vect, VL_full, a, m);
  
  if (t * a <= 1.0) {
    expVL = m_exp_sum(t, m, aux_vect, a);
  }
  else {
    int n{};
    n = std::log(t * a) / std::log(2.0);
    ++n;
    
    expVL = m_exp_sum(t / pow(2.0, n), m, aux_vect, a);
    
    pow2_matrix(n, expVL);
  }
  
  arma::mat VL = extract_submatrix(expVL, 1, 2, p);
  
  return VL;
}

// Fast Van Loan with precomputed matrix A with p x p blocks
arma::mat VanLoan_fast(double t,
                       const arma::mat& A,
                       int p) {
  
  double a = max_diagonal(A * (-1.0));
  
  int m{find_n(uni_epsilon, 1.0)};
  
  arma::mat expVL;
  
  std::vector<arma::mat> aux_vect;
  vector_of_matrices(aux_vect, A, a, m);
  
  if (t * a <= 1.0) {
    expVL = m_exp_sum(t, m, aux_vect, a);
  }
  else {
    int n{};
    n = std::log(t * a) / std::log(2.0);
    ++n;
    
    expVL = m_exp_sum(t / pow(2.0, n), m, aux_vect, a);
    
    pow2_matrix(n, expVL);
  }
  
  arma::mat VL = extract_submatrix(expVL, 1, 2, p);
  
  return VL;
}

//// Robust EM algorithm for MAPs //////////////////////////////////////////////

// Get position of element (i,j) in lexicographically 
// ordered state space of dimension p
int lex_idx(int i, int j, int p){ 
  return i * p + j;
}

// Initial vector starting in state i
arma::rowvec start_in_i(int i, arma::rowvec alpha){
  int p = alpha.size();
  arma::rowvec pi_i_vec(p, arma::fill::zeros);
  pi_i_vec[i] = alpha[i];
  return pi_i_vec;
  
}

// e_i * e_j' used for E-step in Breuer (2002) 
arma::mat e_ij(int i, int j, int p){
  arma::mat e_ij(p, p, arma::fill::zeros);
  e_ij(i,j) = 1.0;
  return e_ij;
}

// Used for E-step in Breuer (2002) 
double c_nij(const std::vector<double>& x, int n, int p, const arma::mat& A, const arma::mat& alpha_matrix, const arma::mat& eta_matrix){
  arma::mat VL = VanLoan_fast(x[n], A, p);
  const arma::rowvec& alpha_row = alpha_matrix.row(n);
  const arma::colvec& eta_col = eta_matrix.col(n);
  double val = 0.0;
  for (int k = 0; k < p; ++k) {
    for (int l = 0; l < p; ++l) {
      val += alpha_row[k] * VL(k, l) * eta_col[l];
    }
  }
  
  return val;
}

// Used for E-step in Breuer (2002) 
double l_nij(const std::vector<double>& x, int n, int i, int j,
             const std::vector<arma::mat>& expC, const arma::mat& D,
             const arma::mat& alpha_matrix, const arma::mat& eta_matrix) {
  double d_ij = D(i,j);
  double val = 0.0;
  const arma::rowvec& alpha_row = alpha_matrix.row(n);
  const arma::colvec& eta_col = eta_matrix.col(n + 1);
  const arma::mat& Cn     = expC[n];
  const arma::mat& Cnplus = expC[n+1];
  
  // Avoid allocating vectors: extract scalar directly
  double left = 0.0;
  double right = 0.0;
  for (int k = 0; k < Cn.n_cols; ++k)
    left += alpha_row[k] * Cn(k, i);
  for (int k = 0; k < Cnplus.n_cols; ++k)
    right += Cnplus(j, k) * eta_col[k];
  
  val = left * d_ij * right;
  
  return val;
}


// Weights for D matrices for improved numerical stability
// [[Rcpp::export]]
std::vector<double> D_weights(const std::vector<double>& x,
                              const arma::rowvec& alpha,
                              const arma::mat& C,
                              const arma::mat& D){
  
  int N = x.size();
  int p = alpha.size();
  
  std::vector<double> weights(N);
  
  arma::rowvec start_vec = alpha;
  arma::rowvec vec(p);
  double vec_sum;
  arma::mat M(p,p);
  
  std::vector<arma::mat> expC = expA_vec(x, C);
  
  for (int i = 0; i < N; ++i) {
    M = expC[i] * D;
    vec = start_vec * M;
    vec_sum = sum(vec);
    weights[i] = vec_sum;
    start_vec = vec / vec_sum;
  }
  
  return weights;
  
}

// Weights for D matrices with precomputed matrix exponentials
std::vector<double> D_weights_EM(const std::vector<double>& x,
                                 const arma::rowvec& alpha,
                                 const std::vector<arma::mat>& expC,
                                 const arma::mat& D){
  
  int N = x.size();
  int p = alpha.size();
  
  std::vector<double> weights(N);
  
  arma::rowvec start_vec = alpha;
  arma::rowvec vec(p);
  double vec_sum;
  arma::mat M(p,p);
  
  for (int i = 0; i < N; ++i) {
    M = expC[i] * D;
    vec = start_vec * M;
    vec_sum = sum(vec);
    weights[i] = vec_sum;
    start_vec = vec / vec_sum;
  }
  
  return weights;
  
}

// Joint density of arrivals x in MAP(alpha, C, D)
// Direct implementation - highly numerically unstable
double MAP_joint_density(const std::vector<double>& x,
                         const arma::rowvec& alpha,
                         const arma::mat& C,
                         const arma::mat& D){
  
  int n = x.size();
  int p = alpha.size();
  
  arma::vec e(p);
  e.ones();
  
  arma::mat trans_matrix = arma::eye(p,p);
  
  std::vector<arma::mat> expC = expA_vec(x, C);
  
  for (int i = 0; i < n; ++i){
    
    trans_matrix = trans_matrix * expC[i] * D;
  }
  
  arma::mat val_mat = alpha * trans_matrix * e;
  return arma::as_scalar(val_mat);
}

// Joint density of arrivals x in MAP(alpha, C, D)
// Robust implementation - numerically unstable
double RMAP_joint_density(const std::vector<double>& x,
                          const arma::rowvec& alpha,
                          const std::vector<arma::mat>& expC,
                          const arma::mat& D,
                          const std::vector<double>& D_w // D_weights
                          ){
  
  int n = x.size();
  int p = alpha.size();
  
  arma::vec e(p);
  e.ones();
  
  arma::mat trans_matrix = arma::eye(p,p);
  for (int i = 0; i < n; ++i){
    trans_matrix = trans_matrix * expC[i] * D * (1.0/ D_w[i]);
  }
  
  arma::mat val_mat = alpha * trans_matrix * e;
  
  return arma::as_scalar(val_mat);
}

// Robust alpha matrix for E-step
arma::mat Ralpha_matrix(const std::vector<double>& x,
                        const arma::rowvec& alpha,
                        const std::vector<arma::mat>& expC,
                        const arma::mat& D,
                        const std::vector<double>& D_w // D_weights
                        ){
  int N = x.size();
  int p = alpha.size();
  
  arma::mat alpha_matrix(N+1, p, arma::fill::zeros);
  alpha_matrix.row(0) = alpha;
  
  for (int n = 1; n < N+1; ++n) {
    alpha_matrix.row(n) = alpha_matrix.row(n - 1) * expC[n-1] * D * (1.0/ D_w[n-1]);
  }
  
  return alpha_matrix;
}

// Robust eta matrix for E-step
arma::mat Reta_matrix(const std::vector<double>& x,
                      const std::vector<arma::mat>& expC,
                      const arma::mat& D,
                      const std::vector<double>& D_w) {
  int N = x.size();
  int p = D.n_cols;
  
  arma::mat eta_matrix(p, N, arma::fill::zeros);
  
  eta_matrix.col(N - 1) = arma::sum(D / D_w[N - 1], 1);
  
  for (int n = N - 2; n >= 0; --n) {
    eta_matrix.col(n) = D * (1.0 / D_w[n]) * expC[n+1] * eta_matrix.col(n + 1);
  }
  
  return eta_matrix;
}

// Robust expectation of sufficient statistic B_i in Breuer (2002)
double REB_i(const std::vector<double> x, int i, 
             const arma::rowvec& alpha, const std::vector<arma::mat>& expC, 
             const arma::mat& D, const std::vector<double>& D_w){
  double val = 0.0;
  if(alpha[i] != 0.0){
    val = RMAP_joint_density(x, start_in_i(i, alpha), expC, D, D_w);
  }
  return val;
}

// Robust expectation of sufficient statistic Z_i in Breuer (2002)
double REZ_i(const std::vector<double>& x, int i, int p,
             const std::vector<arma::mat>& VL_vec,
             const arma::mat& Ralpha_matrix, const arma::mat& Reta_matrix) {
  int N = x.size();
  arma::mat VL_ii = VL_vec[lex_idx(i,i,p)];
  double val = 0.0;
  for (int n = 0; n < N; ++n) {  // n from 0 to N
    val += c_nij(x, n, p, VL_ii, Ralpha_matrix, Reta_matrix);
  }
  return val;
}

// Robust expectation of sufficient statistic Y_i in Breuer (2002)
double REY_i(const std::vector<double>& x, int i, int p,
             const std::vector<arma::mat>& VL_vec,
             const arma::mat& Ralpha_matrix, const arma::mat& Reta_matrix) {
  int N = x.size();
  arma::mat VL_ii = VL_vec[lex_idx(i,i,p)];
  double val = 0.0;
  for (int n = 0; n < N-1; ++n) {  // n from 0 to N-1
    val += c_nij(x, n, p, VL_ii, Ralpha_matrix, Reta_matrix);
  }
  return val;
}

// Robust expectation of sufficient statistic N_ij in Breuer (2002)
double REN_ij(const std::vector<double>& x, int i, int j, 
              int p,
              const arma::mat& C,
              const std::vector<arma::mat>& VL_vec,
              const arma::mat& Ralpha_matrix, 
              const arma::mat& Reta_matrix, const std::vector<double>& D_w) {
  int N = x.size();
  arma::mat VL_ij = VL_vec[lex_idx(i,j,p)];
  double c_ij = C(i,j);
  double val = 0.0;
  if(c_ij != 0.0){
    for (int n = 0; n < N; ++n) {
      val += c_ij * c_nij(x, n, p, VL_ij, Ralpha_matrix, Reta_matrix);
    }
  }
  return val;
}

// Robust expectation of sufficient statistic L_ij in Breuer (2002)
double REL_ij(const std::vector<double>&  x, int i, int j, const std::vector<arma::mat>& expC, const arma::mat& D, const arma::mat& Ralpha_matrix, const arma::mat& Reta_matrix, const std::vector<double>& D_w) {
  int N = x.size();
  double val = 0.0;
  double d_ij = D(i,j);
  if(d_ij != 0.0) {
    for (int n = 0; n < N-1; ++n) {  // n from 0 to N-1
      val += l_nij(x, n, i, j, expC, D, Ralpha_matrix, Reta_matrix) * (1.0/D_w[n]);
    }
  }
  return val;
}

// EM algorithm for MAPs
// [[Rcpp::export]]
void EM_MAP(const std::vector<double> x,
            arma::rowvec &alpha,
            arma::mat &C,
            arma::mat &D,
            bool benchmark = false, // Set true to see runtime for each part of E-step
            double zero_eps = 0.0000000001 // Values below zero_eps are set to zero
            ){
  
  R_CheckUserInterrupt();
  
  int N = x.size();
  int p = alpha.size();
  
  // Van Loan matrices
  arma::mat VL_base(2*p, 2*p, arma::fill::zeros);
  assign_submatrix(VL_base, C, 1, 1);
  assign_submatrix(VL_base, C, 2, 2);
  
  std::vector<arma::mat> VL_vec(p*p, arma::mat(2*p, 2*p));
  
  // Compute VanLoan matrices
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) {
      int idx = lex_idx(i, j, p);
      VL_vec[idx] = VL_base;
      assign_submatrix(VL_vec[idx], e_ij(i,j,p), 1, 2);
    }  
  }
  
  // Compute matrix exponentials
  std::vector<arma::mat> expC = expA_vec(x, C);
  
  arma::rowvec alpha_update(p, arma::fill::zeros);
  arma::mat C_update(p,p, arma::fill::zeros);
  arma::mat D_update(p,p, arma::fill::zeros);
  
  std::vector<double> D_w = D_weights_EM(x,alpha, expC, D);
  
  // Compute auxiliary matrices
  arma::mat Ralpha_mat = Ralpha_matrix(x, alpha, expC, D, D_w);
  arma::mat Reta_mat = Reta_matrix(x, expC, D, D_w);
  
  // Vectors / matrices to store estimates for sufficient statistics
  arma::rowvec REB(p, arma::fill::zeros);
  arma::rowvec REZ(p, arma::fill::zeros);
  arma::rowvec REY(p, arma::fill::zeros);
  
  arma::mat REN(p,p, arma::fill::zeros);
  arma::mat REL(p,p, arma::fill::zeros);
  
  // E step //
  
  if(benchmark){Rcpp::Rcout << "Bi" << std::endl;}
  auto start_BI = std::chrono::high_resolution_clock::now();
  
  // Bi
  for (int i = 0; i < p; ++i) {
    REB[i] = REB_i(x, i, alpha, expC, D, D_w);
  }
  
  auto end_BI = std::chrono::high_resolution_clock::now();
  double dur_BI = std::chrono::duration<double, std::milli>(end_BI - start_BI).count();
  if(benchmark){Rcpp::Rcout << "Bi done in " << dur_BI << " ms" << std::endl;}
  
  if(benchmark){Rcpp::Rcout << "Yi" << std::endl;}
  auto start_YI = std::chrono::high_resolution_clock::now();
  
  // Yi
  for (int i = 0; i < p; ++i) {
    REY[i] = REY_i(x, i, p, VL_vec , Ralpha_mat, Reta_mat);
  }
  
  auto end_YI = std::chrono::high_resolution_clock::now();
  double dur_YI = std::chrono::duration<double, std::milli>(end_YI - start_YI).count();
  if(benchmark){Rcpp::Rcout << "Yi done in " << dur_YI << " ms" << std::endl;}
  
  if(benchmark){Rcpp::Rcout << "Zi" << std::endl;}
  auto start_ZI = std::chrono::high_resolution_clock::now();
  
  
  // Zi # Reusing Yi and adding N'th term
  for (int i = 0; i < p; ++i) {
    REZ[i] = REY[i] + c_nij(x, N-1, p, VL_vec[lex_idx(i,i,p)], Ralpha_mat, Reta_mat);
  }
  
  auto end_ZI = std::chrono::high_resolution_clock::now();
  double dur_ZI = std::chrono::duration<double, std::milli>(end_ZI - start_ZI).count();
  if(benchmark){Rcpp::Rcout << "Zi done in " << dur_ZI << " ms" << std::endl;}
  
  if(benchmark){Rcpp::Rcout << "Nij" << std::endl;}
  auto start_NIJ = std::chrono::high_resolution_clock::now();
  
  // Nij
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) {
      if(i != j){REN(i,j) = REN_ij(x, i, j, p, C, VL_vec, Ralpha_mat, Reta_mat, D_w);}
    }
  }
  
  auto end_NIJ = std::chrono::high_resolution_clock::now();
  double dur_NIJ = std::chrono::duration<double, std::milli>(end_NIJ - start_NIJ).count();
  if(benchmark){Rcpp::Rcout << "Nij done in " << dur_NIJ << " ms" << std::endl;}
  
  if(benchmark){Rcpp::Rcout << "Lij" << std::endl;}
  auto start_LIJ = std::chrono::high_resolution_clock::now();
  
  // Lij
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) {
      REL(i,j) = REL_ij(x, i, j, expC, D, Ralpha_mat, Reta_mat, D_w);
    }
  }
  
  auto end_LIJ = std::chrono::high_resolution_clock::now();
  double dur_LIJ = std::chrono::duration<double, std::milli>(end_LIJ - start_LIJ).count();
  if(benchmark){Rcpp::Rcout << "Lij done in " << dur_LIJ << " ms" << std::endl;}
  
  // M step //
  
  alpha_update = REB;
  
  // Update off-diagonal elements
  for (int i = 0; i < p; ++i) {
    double REZi = REZ[i];
    double REYi = REY[i];
    for (int j = 0; j < p; ++j) {
      if(i != j){
        C_update(i,j) = REN(i,j) / REZi;
      }
      D_update(i,j) = REL(i,j) / REYi;
    }
  }
  
  // Set (very!) small elements to zero
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) {
      if(i != j){
        if(abs(C_update(i,j))< zero_eps){C_update(i,j) = 0.0;}
        if(abs(D_update(i,j))<zero_eps){D_update(i,j) = 0.0;}
      }
    }
  }
  
  // Update diagonal elements
  arma::colvec C_rowsum = arma::sum(C_update, 1);
  arma::colvec D_rowsum = arma::sum(D_update, 1);
  for (int i = 0; i < C_update.n_rows; ++i) {
    C_update(i, i) = -(C_rowsum[i] + D_rowsum[i]);
  }
  
  // Update parameters //
  alpha = alpha_update;
  C = C_update;
  D = D_update;

}

//// Robust EM algorithm for Marked MAPs ///////////////////////////////////////

// Converts a List of matrices to std::vector<arma::mat>
std::vector<arma::mat> list_to_vector_of_matrices(const Rcpp::List& mat_list) {
  std::vector<arma::mat> vector_of_matrices;
  vector_of_matrices.reserve(mat_list.size());
  for (int i = 0; i < mat_list.size(); ++i) {
    vector_of_matrices.push_back(as<arma::mat>(mat_list[i]));
  }
  return vector_of_matrices;
}

// Converts a List of vectors to std::vector<std::vector<double>>
std::vector<std::vector<double>> list_to_vector_of_vectors(const Rcpp::List& vec_list) {
  std::vector<std::vector<double>> vector_of_vectors;
  vector_of_vectors.reserve(vec_list.size());
  for (int i = 0; i < vec_list.size(); ++i) {
    vector_of_vectors.push_back(as<std::vector<double>>(vec_list[i]));
  }
  return vector_of_vectors;
}


// Construct C matrix for k'th marginal MAP
// [[Rcpp::export]]
arma::mat C_k(int k,
              const arma::mat C,
              const arma::mat D,
              const Rcpp::List R_list){
  
  int n = R_list.size(); // Number of marks
  std::vector<arma::mat> R_mats = list_to_vector_of_matrices(R_list);
  arma::mat M = C; // Marginal C matrix
  for (int i = 0; i < n; ++i) {
    if(i != k){
      M = M + R_mats[i] % D;
    }
  }
  return M;
}

// Construct D matrix for k'th marginal MAP
// [[Rcpp::export]]
arma::mat D_k(int k,
              const arma::mat D,
              const Rcpp::List R_list){
  std::vector<arma::mat> R_mats = list_to_vector_of_matrices(R_list);
  return R_mats[k] % D;
}

// [[Rcpp::export]]
Rcpp::List EM_MMAP(const arma::rowvec alpha,
                   const arma::mat C,
                   const arma::mat D,
                   const Rcpp::List x_list,
                   Rcpp::List R_list){
  
  R_CheckUserInterrupt();
  std::vector<std::vector<double>> x_vecs = list_to_vector_of_vectors(x_list);

  int p = alpha.size();
  int n = R_list.size();
  
  auto D_hat = [&alpha, &p](std::vector<double> x_auto, arma::mat C_auto, arma::mat D_auto) {
    // Van Loan matrices
    arma::mat VL_base(2*p, 2*p, arma::fill::zeros);
    assign_submatrix(VL_base, C_auto, 1, 1);
    assign_submatrix(VL_base, C_auto, 2, 2);
    
    std::vector<arma::mat> VL_vec(p*p, arma::mat(2*p, 2*p));
    
    // Compute VanLoan matrices
    for (int i = 0; i < p; ++i) {
      for (int j = 0; j < p; ++j) {
        int idx = lex_idx(i, j, p);
        VL_vec[idx] = VL_base;
        assign_submatrix(VL_vec[idx], e_ij(i,j,p), 1, 2);
      }
    }
    
    // Compute matrix exponentials
    std::vector<arma::mat> expC = expA_vec(x_auto, C_auto);
    
    arma::rowvec alpha_update(p, arma::fill::zeros);
    arma::mat C_update(p,p, arma::fill::zeros);
    arma::mat D_update(p,p, arma::fill::zeros);
    
    std::vector<double> D_w = D_weights_EM(x_auto, alpha, expC, D_auto);
    
    // Compute auxiliary matrices
    arma::mat Ralpha_mat = Ralpha_matrix(x_auto, alpha, expC, D_auto, D_w);
    arma::mat Reta_mat = Reta_matrix(x_auto, expC, D_auto, D_w);
    
    // Vectors / matrices to store estimates for sufficient statistics
    
    arma::rowvec REY(p, arma::fill::zeros);
    arma::mat REL(p,p, arma::fill::zeros);
    
    // E step //
    
    // Yi
    for (int i = 0; i < p; ++i) {
      // REY[i] = REY_i(x, i, C, Ralpha_mat, Reta_mat);
      REY[i] = REY_i(x_auto, i, p, VL_vec , Ralpha_mat, Reta_mat);
    }
    
    // Lij
    for (int i = 0; i < p; ++i) {
      for (int j = 0; j < p; ++j) {
        REL(i,j) = REL_ij(x_auto, i, j, expC, D_auto, Ralpha_mat, Reta_mat, D_w);
      }
    }
    // M step //
    // Update off-diagonal elements
    for (int i = 0; i < p; ++i) {
      double REYi = REY[i];
      for (int j = 0; j < p; ++j) {
        D_update(i,j) = REL(i,j) / REYi;
      }
    }
    return D_update;
  };
  
  std::vector<arma::mat> C_k_vec(n, arma::mat(p, p));
  std::vector<arma::mat> D_k_vec(n, arma::mat(p, p)); 
  
  for (int k = 0; k < n; ++k) {
    C_k_vec[k] = C_k(k, C, D, R_list);
    D_k_vec[k] = D_k(k, D, R_list);
  }
  
  std::vector<arma::mat> D_k_hat_vec(n, arma::mat(p, p));
  for (int k = 0; k < n; ++k) {
    D_k_hat_vec[k] = D_hat(x_vecs[k], C_k_vec[k], D_k_vec[k]);
  }
  
  arma::mat D_hat_full(p,p, arma::fill::zeros);
  for (int k = 0; k < n; ++k) {
    D_hat_full = D_hat_full + D_k_hat_vec[k];
  }
  
  for (int k = 0; k < n; ++k) {
    arma::mat D_k_hat_temp = D_k_hat_vec[k];
    arma::mat R_k(p,p, arma::fill::zeros);
    for (int i = 0; i < p; ++i) {
      for (int j = 0; j < p; ++j) {
        if(D_hat_full(i,j) == 0){continue;}
        R_k(i,j) = D_k_hat_temp(i,j) / D_hat_full(i,j);
      }
    }
    arma::mat R_k_copy = R_k;  // make a deep copy
    R_list[k] = R_k_copy;
  }
  
  return R_list;
}
