#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
// [[Rcpp::export]]
mat costfunction(mat a, mat b,int k = 2){
  mat cost;
  cost.zeros(a.n_rows,b.n_rows);
  for(int i=0;i<a.n_rows;i++){
    for(int j=0;j<b.n_rows;j++){
      cost(i,j) = norm(a.row(i)-b.row(j),k);
    }
  }
  return cost;
}

