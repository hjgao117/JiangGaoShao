
################  Example of functional time series ################ 

source("auxiliary.R")
source("main.R")

library(MASS)
library(fda)

# Generate BM model in Section 5.1
generate_BM = function(n,t){
  sigma = sqrt(1/t)
  x = cbind(0,
            matrix(rnorm(n*t,mean=0,sd=sigma), 
                   nrow=n, ncol=t))
  eps = t(apply(x,1,cumsum))
  return(eps)
}


# Example 
set.seed(2)

# Generate time series y
n = 100
t = 1000
p = 1000
pt.y = t(generate_BM(n,t))
basis = create.fourier.basis(rangeval=c(0,1), nbasis=20)
fd.y = smooth.basis(seq(from=0,to=1,length.out=t+1), pt.y, basis)$fd
data = eval.fd(fd.y,(1:p)/p)
y = t(data)

# Compute distance matrix d
d = compute_d_func(y, "L2")

# Proposed test based on permutation
res.permt = prop_test(d, type="permt", B=300, alpha=0.05)
res.permt

# Proposed test based on bootstrap
res.boot = prop_test(d, type="boot", B=300, alpha=0.05)
res.boot
