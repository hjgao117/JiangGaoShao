
################  Auxiliary Functions ################ 

compute_d_multi = function(y, metric="euc"){
  # Compute the distance matrix for a multivariate time series
  # Input: 
  #   y: n by p matrix, a multivariate time series of size n and dimension p
  #   metric: the metric to be used
  # Output:
  #   d: n by n distance matrix, the (i,j)-component denotes the distance d(yi,yj) between yi and yj 
  
  y = cbind(y)
  if (metric == "euc"){
    # Euclidean metric
    d = as.matrix(dist(y, method="euclidean"))
  }
  return(d)
}


compute_d_func = function(y, metric="L2"){
  # Compute the distance matrix for a functional time series
  # Input: 
  #   y: n by p matrix, a functional time series of size n, each observation is evaluated on p equally-spaced points between [0,1]
  #   metric: the metric to be used
  # Output:
  #   d: n by n distance matrix, the (i,j)-component denotes the distance d(yi,yj) between yi and yj 
  
  p = ncol(y)
  if (metric == "L2"){
    d = as.matrix(dist(y, method="euclidean")) / sqrt(p)
  }
  return(d)
}


compute_d_cov = function(y, metric="euc"){
  # Compute the distance matrix for a covariance matrix time series
  # Input: 
  #   y: array with dimension (n,p,p), a covariance matrix time series of size n and dimension p by p
  #   metric: the metric to be used
  # Output:
  #   d: n by n distance matrix, the (i,j)-component denotes the distance d(yi,yj) between yi and yj 
  
  if (metric == "euc"){
    # Euclidean metric
    euc = function(A,B) sqrt( sum((A-B)^2) )
    d = matrix(0,n,n)
    for (i in 1:(n-1)){
      for (j in (i+1):n)
        d[i,j] = euc(y[i,,], y[j,,])
    }
    d[lower.tri(d)] = t(d)[lower.tri(d)]
  }
  
  if (metric == "log-euc"){
    # log-Euclidean metric
    euc = function(A,B) sqrt( sum((A-B)^2) )
    y.log = array(t(apply(y, 1, 
                          function(x){
                            temp = eigen(x)
                            temp$vectors %*% diag(log(temp$values)) %*% t(temp$vectors)
                          })), 
                  dim=dim(y))
    d = matrix(0,n,n)
    for (i in 1:(n-1)){
      for (j in (i+1):n)
        d[i,j] = euc(y.log[i,,], y.log[j,,])
    }
    d[lower.tri(d)] = t(d)[lower.tri(d)]
  }
  
  if (metric == "chol"){
    # Cholesky metric
    euc = function(A,B) sqrt( sum((A-B)^2) )
    y.chol = array(t(apply(y, 1, function(x) t(chol(x)))), dim=dim(y))
    d = matrix(0,n,n)
    for (i in 1:(n-1)){
      for (j in (i+1):n)
        d[i,j] = euc(y.chol[i,,], y.chol[j,,])
    }
    d[lower.tri(d)] = t(d)[lower.tri(d)]
  }
  
  if (metric == "riemann"){
    # Riemann metric
    riemann = function(A,B){
      e.A = eigen(A)
      inv.sqrt.A = e.A$vectors %*% diag((sqrt(e.A$values))^{-1}) %*% t(e.A$vectors)
      temp = inv.sqrt.A %*% B %*% inv.sqrt.A
      e.temp = eigen(temp)
      return( sqrt( sum( ( e.temp$vectors %*% diag(log(e.temp$values)) %*% t(e.temp$vectors) )^2 ) ) )
    } 
    d = matrix(0,n,n)
    for (i in 1:(n-1)){
      for (j in (i+1):n)
        d[i,j] = riemann(y[i,,], y[j,,])
    }
    d[lower.tri(d)] = t(d)[lower.tri(d)]
  }
  
  return(d)
}


compute_d_dist = function(y, metric="W1"){
  # Compute the distance matrix for a univariate distributional time series
  # Input: 
  #   y: n by p matrix, each row denotes a quantile/distribution function equally sampled from p points between [0,1]
  #   metric: the metric to be used
  # Output:
  #   d: n by n distance matrix, the (i,j)-component denotes the distance d(yi,yj) between yi and yj 
  
  library(seewave)
  
  if (metric == "W1"){
    # W1 metric, requires y to be quantile functions
    d = matrix(0,n,n)
    for (i in 1:(n-1)){
      for (j in (i+1):n)
        d[i,j] = mean(abs( y[i,]-y[j,] ))
    }
    d[lower.tri(d)] = t(d)[lower.tri(d)]
  }
  
  if (metric == "W2"){
    # W2 metric, requires y to be quantile functions
    d = matrix(0,n,n)
    for (i in 1:(n-1)){
      for (j in (i+1):n)
        d[i,j] = sqrt(mean( (y[i,]-y[j,])^2 ))
    }
    d[lower.tri(d)] = t(d)[lower.tri(d)]
  }
  
  if (metric == "KS"){
    # Kolmogorov-Smirnov distance, requires y to be distribution functions
    p = 200
    x = seq(0,1,length.out=p+1)[-1]
    
    d = matrix(0,n,n)
    for (i in 1:(n-1)){
      for (j in (i+1):n)
        d[i,j] = ks.dist(cbind(x,y[i,]),cbind(x,y[j,]))$D
    }
    d[lower.tri(d)] = t(d)[lower.tri(d)]
  }
  
  if (metric == "KL"){
    # Kullback-Leibler distance, requires y to be distribution functions
    p = 200
    x = seq(0,1,length.out=p+1)[-1]
    
    d = matrix(0,n,n)
    for (i in 1:(n-1)){
      for (j in (i+1):n)
        d[i,j] = kl.dist(cbind(x,y[i,]),cbind(x,y[j,]))$D
    }
    d[lower.tri(d)] = t(d)[lower.tri(d)]
  }
  
  if (metric == "IS"){
    # Itakuro-Saito distance, requires y to be distribution functions
    p = 200
    x = seq(0,1,length.out=p+1)[-1]
    
    d = matrix(0,n,n)
    for (i in 1:(n-1)){
      for (j in (i+1):n)
        d[i,j] = itakura.dist(cbind(x,y[i,]),cbind(x,y[j,]))$D
    }
    d[lower.tri(d)] = t(d)[lower.tri(d)]
  }
  
  if (metric == "LS"){
    # Log-spectral distance, requires y to be distribution functions
    p = 200
    x = seq(0,1,length.out=p+1)[-1]
    
    d = matrix(0,n,n)
    for (i in 1:(n-1)){
      for (j in (i+1):n)
        d[i,j] = logspec.dist(cbind(x,y[i,]),cbind(x,y[j,]))
    }
    d[lower.tri(d)] = t(d)[lower.tri(d)]
  }
  
  return(d)
}

