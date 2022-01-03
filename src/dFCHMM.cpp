#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
arma::vec colSums_cpp(const arma::mat& A){
  arma::vec result(A.n_cols); result.fill(0);
  for(int i = 0; i < A.n_cols; i++){
    result(i) = accu(A.col(i));
  }
  return(result);
}

// [[Rcpp::export]]
arma::vec rowSums_cpp(const arma::mat& A){
  arma::vec result(A.n_rows); result.fill(0);
  for(int i = 0; i < A.n_rows; i++){
    result(i) = accu(A.row(i));
  }
  return(result);
}

// [[Rcpp::export]]
arma::vec logvec_c(const arma::vec& x){
  arma::vec result(x.n_elem); result.fill(0);
  for(int i = 0; i < x.n_elem; i++){
    result(i) = log(x(i));
  }
  return(result);
}

// [[Rcpp::export]]
unsigned long factorial_cpp(unsigned n){
  return (n == 0) ? 1 : n*factorial_cpp(n-1);
}

// [[Rcpp::export]]
arma::mat covmat_c(arma::mat Xt, arma::vec mu, arma::vec u){
  int P = Xt.n_cols;
  int TIME = Xt.n_rows;
  arma::mat covest(P,P); covest.fill(0);
  for(int time = 0; time < TIME; time++){
    arma::vec Xtmu = Xt.row(time).t() - mu;
    covest = covest + u(time)*Xtmu*Xtmu.t();
  }
  covest = covest/accu(u);
  return(covest);
}


// [[Rcpp::export]]
arma::mat standardize_rows_cpp(arma::mat A){
  arma::mat Anew = A;
  for(int i = 0; i < A.n_rows; i++){
    arma::rowvec A_i = A.row(i);
    Anew.row(i) = A_i/accu(A_i.t());
  }
  return(Anew);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec dmvnrm_arma(const arma::mat& x,  
                      const arma::rowvec& mean,  
                      const arma::mat& sigma, 
                      bool logd = false) { 
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
  }  
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::export]]
arma::mat eye_cpp(int P){
  arma::mat Id(P,P); Id.fill(0);
  for(int i = 0; i < P; i++){
    Id(i,i) = 1;
  }
  return(Id);
}


// [[Rcpp::export]]
double logsumexp_cpp(arma::vec X){
  double maxvec = X.max();
  return(maxvec + log(accu(exp(X-maxvec))));
}

// [[Rcpp::export]]
List forward_backward_cpp(const arma::mat& Xt,
                         const arma::mat& B,
                         const arma::mat& mu,
                         const arma::cube& Sigma,
                         const arma::vec& delta){
  // B is transition probability matrix
  // mu is KxP matrix containing means, each row being the mean vector for that state
  // Sigma is PxPxK array, each slice being that state's covariance matrix
  int numStates = mu.n_rows;
  int TIME = Xt.n_rows;
  
  double fb_verysmallvalue = pow(4.940656,-142);
  
  // Get stationary distribution
  arma::mat ones(numStates,numStates); ones.fill(1);
  
  // Compute pdf at each time point
  arma::mat allprobs(TIME,numStates); allprobs.fill(fb_verysmallvalue);
  for(int j = 0; j < numStates; j++){
      allprobs.col(j) = dmvnrm_arma(Xt,mu.row(j),Sigma.slice(j),FALSE);  
  }
  
  // Check if valid probability
  for(int time = 0; time < TIME; time++){
    if(all(allprobs.row(time) <= 1e-70)){
      allprobs.row(time).fill(1);
    }
    for(int j = 0; j < numStates; j++){
      if(allprobs(time,j) < 0){
        allprobs(time,j) = -999;
        cout << "probability < 0" << endl;
      }
      if(!std::isfinite(allprobs(time,j))){
        allprobs(time,j) = -999;
        cout << "probability NaN" << endl;
      }
    }
  }
  
  // Forward variables
  arma::mat lalpha(numStates,TIME); lalpha.fill(fb_verysmallvalue);
  arma::vec foo = delta%allprobs.row(0).t(); // element-wise multiplication
  double sumfoo = accu(foo) + fb_verysmallvalue;
  foo = foo/sumfoo;
  double lscale = std::log(sumfoo);
  lalpha.col(0) = logvec_c(foo) + lscale;
  for(int time = 1; time < TIME; time++){
    arma::mat tempB(B.n_rows,B.n_cols); tempB.fill(0);
    for(int i = 0; i < B.n_rows; i++){
      tempB.row(i) = (B.row(i)%allprobs.row(time));
    }
    foo = (foo.t()*tempB).t();
    sumfoo = accu(foo) + fb_verysmallvalue;
    lscale = lscale + std::log(sumfoo);
    foo = foo/accu(foo);
    lalpha.col(time) = logvec_c(foo) + lscale;
  }
  double llk_alpha = lscale;
  
  // Backward variables
  arma::mat lbeta(numStates,TIME); lbeta.fill(fb_verysmallvalue);
  arma::vec foo2(numStates); foo2.fill((double)(1/(double)numStates));
  lscale = std::log(numStates);
  for(int time = TIME-2; time >= 0; time--){
    foo2 = B*(allprobs.row(time+1)%foo2.t()).t();
    lbeta.col(time) = logvec_c(foo2) + lscale;
    double sumfoo2 = accu(foo2) + fb_verysmallvalue;
    foo2 = foo2/sumfoo2;
    lscale = lscale + std::log(sumfoo2);
  }
  double llk_beta = lscale;
  return(List::create(_["lalpha"] = lalpha, 
                      _["lbeta"] = lbeta,
                      _["llk_alpha"] = llk_alpha,
                      _["llk_beta"] = llk_beta));
}

// [[Rcpp::export]]
List hmm_cpp(const arma::cube& Xt, 
             int num_states, 
             const arma::cube& Sigma_init, 
             int maxiter = 10000, 
             double tol = 1e-4, 
             bool verbose = 1){
  int TIME = Xt.n_rows;
  int N = Xt.n_slices;
  int P = Xt.n_cols;
  
  // Initialize mean vector and TPM
  arma::mat mu(num_states,P); mu.fill(0);
  arma::cube A(num_states,num_states,N); A.fill((1-.8)/(num_states-1));
  for(int n = 0; n < N; n++){
    for(int i = 0; i < num_states; i++){
      A(i,i,n) = .8;
    }  
  }
  
  arma::cube A_hat_next = A; 
  arma::mat mu_hat_next = mu;
  arma::cube Sigma_hat = Sigma_init;
  arma::cube Sigma_hat_next = Sigma_init;
  arma::mat llk(N,maxiter); llk.fill(-999);
  arma::cube u_all(TIME,num_states,N);
  arma::mat delta(N,num_states); delta.fill((double)1/((double)num_states));
  arma::mat delta_next = delta;
  
  for(int iter = 0; iter < maxiter; iter++){
    if(verbose == 1){
      cout << "Iteration " << iter+1 << " of " << maxiter << "..." << endl;
    }
    
    // Reset parameters
    A_hat_next.fill(0);
    mu_hat_next.fill(0);
    Sigma_hat_next.fill(0);
    
    for(int n = 0; n < N; n++){
      arma::mat Xtn = Xt.slice(n);
      
      // Forward-backward algorithm
      List fb = forward_backward_cpp(Xtn,A.slice(n),mu,Sigma_hat,delta.row(n).t()); 
      arma::mat lalpha = fb[0];
      arma::mat lbeta = fb[1];
      arma::mat lprobs(TIME,num_states); lprobs.fill(0);
      for(int j = 0; j < num_states; j++){
        lprobs.col(j) = dmvnrm_arma(Xtn,mu.row(j),Sigma_hat.slice(j),TRUE);
      }
      
      // Get transition probabilities
      llk(n,iter) = fb[2];
      arma::mat Atemp_mat(num_states, num_states); Atemp_mat.fill(0);
      
      for(int j = 0; j < num_states; j++){
        for(int k = 0; k < num_states; k++){
          arma::vec lalpha_temp = lalpha.row(j).subvec(0,TIME-2).t();
          arma::vec lbeta_temp = lbeta.row(k).subvec(1,TIME-1).t();
          arma::vec lprobs_temp = lprobs.col(k).subvec(1,TIME-1);
          arma::vec logAjk_temp(TIME-1); logAjk_temp.fill(std::log(A(j,k,n)));
          Atemp_mat(j,k) = accu(exp(lalpha_temp+lbeta_temp+lprobs_temp + logAjk_temp - llk(n,iter)));
        }
      }
      
      // Standardize each row
      Atemp_mat = standardize_rows_cpp(Atemp_mat);
      A_hat_next.slice(n) = Atemp_mat;
      
      // Stationary probabilities
      for(int j = 0; j < num_states; j++){
        delta_next(n,j) = std::exp(lalpha(j,0) + lbeta(j,0) - llk(n,iter));
        if(delta_next(n,j)<=0){
          delta_next(n,j) = pow(4.940656,-142);
        }
      }
      delta_next.row(n) = delta_next.row(n) / accu(delta_next.row(n));
      
      // Get emission distribution parameters
      arma::mat u_mat(TIME,num_states); u_mat.fill(0);
      for(int j = 0; j < num_states; j++){
        arma::rowvec u = exp(lalpha.row(j) + lbeta.row(j) - llk(n,iter));
        u_mat.col(j) = u.t();
        for(int p = 0; p < P; p++){
          mu_hat_next(j,p) += accu(u.t()%Xtn.col(p))/accu(u.t())/N; // % is element-wise multiplication
        }
        Sigma_hat_next.slice(j) += covmat_c(Xtn,mu_hat_next.row(j).t(),u.t())/N;
      }
      u_all.slice(n) = u_mat;
    }
    // Determine convergence
    if(iter > 0){
      if(verbose == 1){
        cout << "Relative Change in likelihood = " << accu(llk.col(iter)-llk.col(iter-1))/fabs(accu(llk.col(iter-1))) << endl;
      }
      if(fabs(accu(llk.col(iter)-llk.col(iter-1))/fabs(accu(llk.col(iter-1)))) < tol){
        cout << "Converged after " << iter << " iterations" << endl;
        return(List::create(_["mu.hat"] = mu_hat_next,
                            _["Sigma.hat"] = Sigma_hat_next,
                            _["A.hat"] = A_hat_next,
                            _["delta.hat"] = delta_next,
                            _["llk"] = llk.submat(0,0,N-1,iter),
                            _["u.all"] = u_all));
      }
    }
    
    // Update parameters
    mu = mu_hat_next;
    Sigma_hat = Sigma_hat_next;
    A = A_hat_next;
    delta = delta_next;
  }
  cout << "No convergence after " << maxiter << " iterations" << endl;
  return(List::create(_["mu.hat"] = mu_hat_next,
                      _["Sigma.hat"] = Sigma_hat_next,
                      _["A.hat"] = A_hat_next,
                      _["delta.hat"] = delta_next,
                      _["llk"] = llk,
                      _["u.all"] = u_all));
}
