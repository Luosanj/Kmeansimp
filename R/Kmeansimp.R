#' Improvement Kmeans
#'
#'@description
#'This package is an improved version of Lloyd kmeans clustering
#'Standard sampling method is to randomly select k points( k cluster number)
#'as initial centroids. Yet this way usually converges too slow.
#'The improved version has 3 steps:
#'Step 1: From all objects randomly select one.
#'Step 2: select next initial centroids from n-objects in such a way that the
#'Euclidean distance of that object is maximum from other selected initial centroids.
#'Step 3: repeat step 2 until we get k initial centroids.
#'
#' @param data The dataset needed to be clustered(type: matrix or dataframe)
#' @param n Number of clusters
#' @param iter Maximum iteration times for updating clustering result
#'
#' @useDynLib Kmeansimp 
#' @importFrom Rcpp sourceCpp
#' @return A list including all important information of clustering
#' @export
#'
#' @examples
#' library(Kmeansimp)
#' X = matrix(runif(100, 0, 100), 25, 4)
#' Kmeansimp(X, 3, iter = 15)
#' clu = Kmeansimp(X, 2)
#' clu$Clustering_Vector
Kmeansimp = function(data, n, iter = 10){
  # Kmeans algorithm is 'Lloyd'
  if(!is.matrix(data) && !is.data.frame(data)){
    stop('Input dataset must be a dataframe or matrix !')
  }
  if(!(n == round(n)) || n < 1L){
    stop('Number of clusters must be a positive integer!')
  }
  size = nrow(data)
  if(n >= size){
    stop('Number of clusters must be above 0 and below sample size !')
  }
  if(n == 1L){
    warning('1 cluster is the dataset itself!')
  }
  
  dataset = as.matrix(data)
  maxiter = iter
  k = n
  oriCluster = c()
  curCluster = c()
  
  # ********Innovation : generate centroids********
  centroids = matrix(nrow = k, ncol = ncol(dataset))
  centroids[1,] = dataset[sample(nrow(dataset), 1),]
  exdata = dataset
  sizex = nrow(exdata)
  i = 2
  while(i <= k){
    exdist = acCtr(exdata, centroids, i)
    dindex = which.max(exdist)
    centroids[i, ] = exdata[dindex, ]
    exdata = exdata[-dindex,, drop = FALSE]
    i = i + 1
  }
  
  # do first clustering
  curCluster = doClus(dataset, centroids)
  
  # recluster until converge or maximum iterations
  flag = FALSE
  numiter = 2L
  while(numiter <= maxiter && !flag){
    oriCluster = curCluster
    # update centroids
    for(j in c(1L : k)){
      clu_index = which(oriCluster == j)
      points = dataset[clu_index,, drop = FALSE]
      centroids[j, ] = colMeans(points)
    }
    curCluster = doClus(dataset, centroids)
    if(identical(oriCluster, curCluster)){
      flag = TRUE
    }
    numiter = numiter + 1L
  }
  if(!flag){
    warning('Clustering did not converge, more iterations are suggested!')
  }
  
  # count number of each cluster
  sizeVec = c()
  for(j in c(1L : k)){
    sizeVec[j] = sum(curCluster == j)
  }
  
  # clusters within sum of squares and proportion of total sum of squares
  distSS = c(rep(0L, k))
  for(j in c(1L : k)){
    index = which(curCluster == j)
    distSS[j] = sum(sweep(dataset[index, ], 2L, centroids[j, ])^2L)
  }
  centroidAll = colMeans(dataset)
  totalSS = sum(sweep(dataset, 2L, colMeans(dataset))^2L)
  prop = sum(distSS) / totalSS
  
  # output list of different information
  result = list()
  result[[1]] = sizeVec
  result[[2]] = centroids
  result[[3]] = curCluster
  result[[4]] = distSS
  result[[5]] = round(1 - prop, 2L)
  names(result) = c('each_cluster_size',
                    'Cluster_means',
                    'Clustering_Vector:',
                    'Within_cluster_sum_of_squares(withinSS)',
                    'proportion of betweenSS with totalSS')
  return(result)
}











