
################  Example of distributional matrix time series ################ 

source("auxiliary.R")
source("main.R")

library(MASS)
library(transport)
library(seewave)
library(fdapace)

# ATM(0) model in Section 5.3
map2boundary = function(Logfit, sup) {
  if (!is.unsorted(Logfit + sup))
    return (1)
  eta0 = 1
  eta1 = 0.5
  while (abs(eta0 - eta1) > 1e-6 | is.unsorted(eta1 * Logfit + sup)) {
    if (is.unsorted(eta1 * Logfit + sup)){
      tmp = eta1
      eta1 = eta1 - (eta0 - eta1) / 2
      eta0 = tmp
    } else {
      eta1 = (eta0 + eta1) / 2
    }
  }
  return(eta1)
}

generate_eps = function(n,p){
  x = seq(0,1,length.out=p+1)
  grid = x[-c(1, p+1)]
  
  xi = runif(n, min = -1, max = 1)
  f = splinefun(c(0, 0.33, 0.66, 1), c(0, 0.2, 0.8, 1), method = "hyman")
  Tf = vector("list", length = n)
  ITf = vector("list", length = n)
  
  for (i in c(1:n)){
    fxi = splinefun((1-xi[i])*f(x)/2+(1+xi[i])*x/2, (1+xi[i])*f(x)/2+(1-xi[i])*x/2, method = "hyman")
    temp = fxi(x)
    if(is.unsorted(temp)){
      temp = x + (temp-x)*map2boundary(temp-x, x)
      temp = fdapace::Lwls1D(bw=0.05, kernel_type = "epan", xin=x, yin=temp, xout=x)
    }
    Tf[[i]] = splinefun(x, temp, method = "hyman")
    ITf[[i]] = splinefun(temp, x, method = "hyman")
  }
  
  Q = t(sapply(c(1:n), function(i) {return(Tf[[i]](x, 0))})) # nxp quantile function
  return(Q[,-1])
}


# Example 
set.seed(4)

# Generate time series y
n = 100
p = 200
y = generate_eps(n,p)

# Compute distance matrix d
d = compute_d_dist(y, "W1")

# Proposed test based on permutation
res.permt = prop_test(d, type="permt", B=300, alpha=0.05)
res.permt

# Proposed test based on bootstrap
res.boot = prop_test(d, type="boot", B=300, alpha=0.05)
res.boot
