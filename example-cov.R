
################  Example of covariance matrix time series ################ 

source("auxiliary.R")
source("main.R")

library(MASS)
library(mvtnorm)

# CAW model in Section 5.2
generate_Wishart = function(n,df,Sigma){
  y = aperm(rWishart(n,df,Sigma), c(3,1,2))
  return(y)
}

generate_ConditionalAuto = function(n,p=3,df=10,rho=1){
  Omega = diag(1, p)
  A = diag(0.7, p) * rho
  B = diag(0.5, p) * rho
  
  eps = generate_Wishart(n+50,df,diag(rep(1,p))) / df
  y = Sigma = array(0, dim=c(n+50,p,p))
  
  Sigma[1,,] = Omega
  temp = t(chol(Sigma[1,,]))
  y[1,,] = temp %*% eps[1,,] %*% t(temp)
  for (i in 2:(n+50)){
    Sigma[i,,] = Omega + A %*% y[i-1,,] %*% t(A) + B %*% Sigma[i-1,,] %*% t(B)
    temp = t(chol(Sigma[i,,]))
    y[i,,] = temp %*% eps[i,,] %*% t(temp)
  }
  return(y[-c(1:50),,])
}


# Example 
set.seed(3)

# Generate time series y
n = 100
p = 2
y = generate_ConditionalAuto(n,p,df=10,rho=0.75)

# Compute distance matrix d
d = compute_d_cov(y, "euc")

# Proposed test based on permutation
res.permt = prop_test(d, type="permt", B=300, alpha=0.05)
res.permt

# Proposed test based on bootstrap
res.boot = prop_test(d, type="boot", B=300, alpha=0.05)
res.boot
