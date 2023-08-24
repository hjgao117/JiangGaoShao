
################  Proposed Method ################ 

compute_AutoDist = function(d){
  # Input:
  #   d: n by n distance matrix, the (i,j)-component denotes the distance d(yi,yj) between yi and yj
  # Output:
  #   Ax, Ay: U-centered versions of the auto-distance matrix
  
  n = dim(d)[1]
  num.k = n - 4
  Ax = Ay = array(NA,dim=c(n,n,num.k))
  diag(d) = 0
  for (k in 1:num.k){
    dx = d[(k+1):n,(k+1):n]
    dy = d[1:(n-k),1:(n-k)]
    Ax[1:(n-k),1:(n-k),k] =  dx - outer(rowSums(dx)/(n-k-2), colSums(dx)/(n-k-2), "+") + sum(dx) / (n-k-1) / (n-k-2)
    Ay[1:(n-k),1:(n-k),k] =  dy - outer(rowSums(dy)/(n-k-2), colSums(dy)/(n-k-2), "+") + sum(dy) / (n-k-1) / (n-k-2)
    diag(Ax[1:(n-k),1:(n-k),k]) = diag(Ay[1:(n-k),1:(n-k),k]) = 0
  }
  
  return(list(Ax = Ax,
              Ay = Ay))
}


prop_test = function(d, type="permt", B=500, alpha=0.05){
  # Input:    
  #   d: n by n distance matrix, the (i,j)-component denotes the distance d(yi,yj) between yi and yj 
  #   type: permt or boot
  #   B: number of replications for permutation/bootstrap
  #   alpha: significance level
  # Output:
  #   the test statistics of CM and KS
  #   the empirical quantiles
  #   the test decision
  #   the computing time
  
  n = dim(d)[1]
  num.k = n-4
  zeta = seq(from=0,to=pi,by=0.01)
  
  res = numeric(7)
  names(res) = c(paste(rep(c("stat","quantile","decision"),2),
                       rep(c("CM","KS"),each=3),
                       sep="-"),
                 "time")
  rep.CM = rep.KS = numeric(B)
  
  # Perform the test based on permutation
  if (type == "permt"){
    st = Sys.time()
    
    # Compute the U-centered auto-distance matrices
    A = compute_AutoDist(d)
    Ax = A$Ax
    Ay = A$Ay
    
    # Compute the test statistic
    Sn = rbind(apply(Ax*Ay, 3, sum, na.rm=T) / (n-(1:num.k)-3) / ((1:num.k)*pi)) %*% sin(outer(1:num.k, zeta, "*"))
    stat.CM = mean(Sn^2) * pi 
    stat.KS = max(abs(Sn))
    
    # Compute the empirical quantile via permutation
    for (b in 1:B){
      order_b = sample(1:n,n)
      A.b = compute_AutoDist(d[order_b,order_b])
      Ax.b = A.b$Ax
      Ay.b = A.b$Ay
      Sn.b = rbind(apply(Ax.b*Ay.b, 3, sum, na.rm=T) / (n-(1:num.k)-3) / ((1:num.k)*pi)) %*% sin(outer(1:num.k, zeta, "*"))
      rep.CM[b] = mean(Sn.b^2) * pi 
      rep.KS[b] = max(abs(Sn.b))
    }
    qt.CM = quantile(rep.CM, 1-alpha)
    qt.KS = quantile(rep.KS, 1-alpha)
    dc.CM = stat.CM > qt.CM
    dc.KS = stat.KS > qt.KS
    
    ed = Sys.time()
    res[1:7] = c(stat.CM, qt.CM, dc.CM, 
                 stat.KS, qt.KS, dc.KS,
                 round(as.numeric(difftime(ed,st),units="secs"),2))
  }
  
  # Perform the test based on bootstrap
  if (type == "boot"){
    st = Sys.time()
    
    # Compute the U-centered auto-distance matrices
    A = compute_AutoDist(d)
    Ax = A$Ax
    Ay = A$Ay
    
    # Compute the test statistic
    Sn = rbind(apply(Ax*Ay, 3, sum, na.rm=T) / (n-(1:num.k)-3) / ((1:num.k)*pi)) %*% sin(outer(1:num.k, zeta, "*"))
    stat.CM = mean(Sn^2) * pi 
    stat.KS = max(abs(Sn))
    
    # Compute the empirical quantile via permutation
    temp = c(1,-1)
    prob = c(1/2,1/2)
    for (b in 1:B){
      aux.b = numeric(num.k)
      for (k in 1:num.k){
        w = sample(temp, n-k, replace=TRUE, prob=prob)
        aux.b[k] = sum(Ax[1:(n-k),1:(n-k),k] * Ay[1:(n-k),1:(n-k),k] * outer(w,w,"*"), na.rm=T)  / (n-k-3)
      }
      Sn.b = rbind(aux.b / ((1:num.k)*pi)) %*% sin(outer(1:num.k, zeta, "*"))
      rep.CM[b] = mean(Sn.b^2) * pi 
      rep.KS[b] = max(abs(Sn.b))
    }
    qt.CM = quantile(rep.CM, 1-alpha)
    qt.KS = quantile(rep.KS, 1-alpha)
    dc.CM = stat.CM > qt.CM
    dc.KS = stat.KS > qt.KS
    
    ed = Sys.time()
    res[1:7] = c(stat.CM, qt.CM, dc.CM, stat.KS, qt.KS, dc.KS,
                 round(as.numeric(difftime(ed,st),units="secs"),2))
  }
  
  return(res)
}
