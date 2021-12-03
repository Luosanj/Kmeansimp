
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector acCtr(NumericMatrix exdata, NumericMatrix centroids, int i){
  int size = exdata.nrow();
  NumericVector exdist (size);
  if (i == 2) {
    for(int l = 1; l <= size; ++l){
      exdist[l - 1] = sum(pow(centroids(0, _) - exdata(l - 1, _), 2)) ;
    }
  }else{
    for(int l = 1; l <= size; ++l){
      for(int k = 1; k <= (i - 1); ++k){
        exdist[l - 1] += sum(pow(centroids(k - 1, _) - exdata(l - 1, _), 2)) ;
      }
    }
  }
  return exdist ;
}