
################  Example of multivariate time series ################ 

source("auxiliary.R")
source("main.R")

# Generate ARCH(2) model in Appendix A.2
generate_ARCH = function(n,Sigma){
  eps = MASS::mvrnorm(n=n+50,mu=c(0,0),Sigma=Sigma)
  y = h = matrix(0,nrow=n+50,ncol=2)
  h[1,] = c(0.003,0.005)
  y[1,] = sqrt(h[1,]) * eps[1,]
  
  for (t in 2:(n+50)){
    h[t,] = c(0.003,0.005) +
      rbind(c(0.2,0.1),c(0.1,0.3)) %*% (y[t-1,]^2) +
      rbind(c(0.4,0.05),c(0.05,0.5)) %*% h[t-1,]
    y[t,] = sqrt(h[t,]) * eps[t,]
  }
  return(y[-c(1:50),])
}


# Example 
set.seed(1)

# Generate time series y
n = 100
Sigma = matrix(0.4,2,2)+diag(rep(0.6,2))
y = generate_ARCH(n,Sigma)

# Compute distance matrix d
d = compute_d_multi(y, "euc")

# Proposed test based on permutation
res.permt = prop_test(d, type="permt", B=300, alpha=0.05)
res.permt

# Proposed test based on bootstrap
res.boot = prop_test(d, type="boot", B=300, alpha=0.05)
res.boot
