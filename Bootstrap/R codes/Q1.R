# Jacknife Estimator

Jackknife_confidence_interval <- function( v1, statfunc=sd, alpha = 0.05 ) {
  
  n1 <- length(v1)
  jackvec <- NULL
  mu0 <- statfunc(v1)   # Original value from the data
  for(i in 1:n1)
  {
    mua<-statfunc(v1[-i])                     # Function value by removing the ith data 
    jackvec<-c(jackvec, n1*(mu0)-(n1-1)*mua)  # Jacknife bias corrected term
    
  }
  jackbias<-mean(jackvec)-mu0  # Bias calculation 
  jacksd<-sd(jackvec)          # Standard deviation for jacknife estimation

  # Jacknife confidence interval 
  NLB<-mean(v1)-(sd(v1)/sqrt(n1))*qt(1-alpha/2, df = n1-1)  # Lower bound
  NUB<-mean(v1)+(sd(v1)/sqrt(n1))*qt(1-alpha/2, df = n1-1)  # Upper bound
  
  # Jacknife results
  list(mu0=mu0,jackestimate = mean(jackvec), jackbias=jackbias,jacksd=jacksd, normal.confidence.interval=c(NLB,NUB)) 
}

my.bootstrapci.ml<-function(vec0, nboot=10000, alpha=0.1){
  
  # Sample size, mean and standard deviation from the original data
  n0 <- length(vec0) 
  mean0 <- mean(vec0) 
  sd0 <- sqrt(var(vec0))
  
  # Vector to store the location of the bootstrap studentized deviation vector 
  bootvec<-NULL
  bootbiasvec<-NULL
  bootvec_estimate <- NULL
  jacckbiasvec<-NULL
  jacksdvec<-NULL
  
  # Bootstrap simulation nboot times using a for loop 
  for( i in 1:nboot){
    vecb<-sample(vec0,replace=T) # samples from original data with replacement
    meanb<-mean(vecb)
    sdb<-sqrt(var(vecb))
    
    # Note: since resampling full vector we can use n0 for sample size of vecb
    
    bootvec<-c(bootvec,(meanb-mean0)/(sdb/sqrt(n0)))   # Bias Correction by subtracting 1st moment and dividing by 2d moment
    bootbiasvec<-c(bootbiasvec, meanb-mean0 )          # Bias calculation
    bootvec_estimate <- c(bootvec_estimate, meanb)     # Parameter (mean) estimation
  }
  hist(bootvec)
  bootsd <- sd(bootvec)  
  bootbias<- mean(bootbiasvec)
  
  # Calculate lower and upper quantile of the bootstrap distribution
  
  lq<-quantile(bootvec,alpha/2)    # Lower quantile
  uq<-quantile(bootvec,1-alpha/2)  # Upper quantile
  
  LB<-mean0-(sd0/sqrt(n0))*uq      # Lower bound
  UB<-mean0-(sd0/sqrt(n0))*lq      # Upper bound
 
  # Normal confidence interval
  NLB<-mean0-(sd0/sqrt(n0))*qnorm(1-alpha/2) 
  NUB<-mean0+(sd0/sqrt(n0))*qnorm(1-alpha/2)
  
  # Bootstrap results
  list( bootestimate = mean(bootvec_estimate), bootbias = bootbias, bootsd = bootsd, bootstrap.confidence.interval=c(LB,UB),normal.confidence.interval=c(NLB,NUB)) 
}

# Test by generating data
v1 = rnorm(500, mean = 10, sd = 2)
Jackknife_confidence_interval(v1,mean)
my.bootstrapci.ml(v1, nboot = 1000,alpha=0.05)
