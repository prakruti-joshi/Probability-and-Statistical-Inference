##original jackknife

Jackknife<-function( v1, statfunc=sd ) {
 
  n1 <- length(v1)
  jackvec <- NULL
  
  
  mu0 <- statfunc(v1) 
  
  for(i in 1:n1)
  {
    mua<-statfunc(v1[-i])
    jackvec<-c(jackvec, n1*(mu0)-(n1-1)*mua)
    
  }
  jackbias<-mean(jackvec)-mu0 
  jacksd<-sd(jackvec)
  
  #list(mu0=mu0,jackbias=jackbias,jacksd=jacksd, mean0 = mean0, jack_mean = jack_mean) 
  list(mu0=mu0,jackbias=jackbias,jacksd=jacksd) 

}

##original jackknife

#mean(pseudo) ±(t0.975,n−1)* root(var(pseudo)/n)
#mean(pseudo) + qt(0.975,length(x)-1)*sqrt(var(pseudo)/length(x))
#mean(pseudo) - qt(0.975,length(x)-1)*sqrt(var(pseudo)/length(x))
#bootstrap:
#mean0-(sd0/sqrt(n0))*qnorm(1-alpha/2) 

Jackknife_confidence_interval <- function( v1, statfunc=sd, alpha = 0.05 ) {
  
  n1 <- length(v1)
  jackvec <- NULL
  
  
  mu0 <- statfunc(v1) 
  
  
  for(i in 1:n1)
  {
    mua<-statfunc(v1[-i])
    jackvec<-c(jackvec, n1*(mu0)-(n1-1)*mua)
    
  }
  jackbias<-mean(jackvec)-mu0 
  jacksd<-sd(jackvec)
  
  lq<-quantile(jackvec,alpha/2) 
  uq<-quantile(jackvec,1-alpha/2)
  
  #ADD the other two confidence intervals.
  #incorporate into the bootstrap confidence interval (what algebra supports this?) and output result 
  LB<-mean(v1)-(sd(v1)/sqrt(n1))*uq
  UB<-mean(v1)-(sd(v1)/sqrt(n1))*lq
  
  NLB<-mean(v1)-(sd(v1)/sqrt(n1))*qnorm(1-alpha/2) 
  NUB<-mean(v1)+(sd(v1)/sqrt(n1))*qnorm(1-alpha/2)
  
  #list(mu0=mu0,jackbias=jackbias,jacksd=jacksd, mean0 = mean0, jack_mean = jack_mean) 
  list(mu0=mu0,jackbias=jackbias,jacksd=jacksd, jack.confidence.interval=c(LB,UB), normal.confidence.interval=c(NLB,NUB)) 
  
}

Jackknife_confidence_interval_trial <-function( v1, statfunc=sd, alpha = 0.05 ) {
  
  n1 <- length(v1)
  jackvec <- NULL
  jackbiasvec <- NULL
  
  mean0 <- statfunc(v1) 
  sd0 <- sd(v1)
  
  
  for(i in 1:n1)
  {
    meanb <- statfunc(v1[-i])
    sdb<-sqrt(var(v1[-i]))
    
    #jackvec<-c(jackvec, n1*(mu0)-(n1-1)*mua)
    
    jackvec<-c(jackvec,(meanb-mean0)/(sdb/sqrt(n1-1))) 
    jackbiasvec<-c(jackbiasvec,meanb-mean0)
    
  }
  jackbias<-mean(jackbiasvec)
  #jacksd<-sd(jackvec)
  
  lq<-quantile(jackvec, alpha/2) 
  uq<-quantile(jackvec, 1-alpha/2)
  
  LB<-mean0-(sd0/sqrt(n1))*uq
  UB<-mean0-(sd0/sqrt(n1))*lq
  
  ##NORMAL 
  NLB<-mean0-(sd0/sqrt(n1))*qnorm(1-alpha/2) 
  NUB<-mean0+(sd0/sqrt(n1))*qnorm(1-alpha/2)
  
  list(pivotal.confidence.interval = c(LB,UB), normal.confidence.interval=c(NLB,NUB)) 
}

v1 = c(1:1000)
Jackknife(v1)
Jackknife_confidence_interval_trial(v1,mean)

## Q1 ##
