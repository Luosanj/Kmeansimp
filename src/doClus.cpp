
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector doClus(NumericMatrix dataset, NumericMatrix centroids){
  int size = dataset.nrow();
  NumericVector curCluster (size);
  int k = centroids.nrow();
  for (int i = 1; i <= size; ++i){
    NumericVector dist (k);
    for (int j = 1; j <= k; ++j){
      NumericVector distone = dataset(i - 1, _ ) - centroids(j - 1, _ ) ;
      dist[j - 1] = sum(pow(distone, 2)) ;
    }
    curCluster[i - 1] = which_min(dist) + 1 ;
  }
  return curCluster ;
}