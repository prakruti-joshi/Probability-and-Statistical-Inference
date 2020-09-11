## Original bootstap 

my.bootstrapci.ml<-function(vec0,nboot=10000,alpha=0.1)
{
  #extract sample size, mean and standard deviation from the original data
  n0 <- length(vec0) 
  mean0 <- mean(vec0) 
  sd0 <- sqrt(var(vec0))
  
  # create a vector to store the location of the bootstrap studentized deviation vector 
  bootvec<-NULL
  bootbiasvec<-NULL
  jacckbiasvec<-NULL
  jacksdvec<-NULL
  bootvec_trial <- NULL       ## addition
  
  #create the bootstrap distribution using a for loop
  for( i in 1:nboot){
      vecb<-sample(vec0,replace=T) #create mean and standard deviation to studentize
      meanb<-mean(vecb)
      sdb<-sqrt(var(vecb))
      #note since resampling full vector we can use n0 for sample size of vecb
      bootvec<-c(bootvec,(meanb-mean0)/(sdb/sqrt(n0))) 
      bootbiasvec<-c(bootbiasvec,meanb-mean0)
      bootvec_trial <- c(bootvec_trial, meanb)  ## addition
  }
  
  #Calculate lower and upper quantile of the bootstrap distribution
  bootbias<-mean(bootbiasvec)
  
  mean_estimate <- mean(bootvec_trial)  ## addition
  
  lq<-quantile(bootvec,alpha/2) 
  uq<-quantile(bootvec,1-alpha/2)
  
  #ADD the other two confidence intervals.
  #incorporate into the bootstrap confidence interval (what algebra supports this?) and output result 
  LB<-mean0-(sd0/sqrt(n0))*uq
  UB<-mean0-(sd0/sqrt(n0))*lq
  
  #print(lq)
  #print(uq)
  
  #since I have the mean and standard deviation calculate the normal confidence interval here as well 
  NLB<-mean0-(sd0/sqrt(n0))*qnorm(1-alpha/2) 
  NUB<-mean0+(sd0/sqrt(n0))*qnorm(1-alpha/2)
  
  list( bootstrap.confidence.interval=c(LB,UB),normal.confidence.interval=c(NLB,NUB)) 
}

v1 = c(1:1000)
my.bootstrapci.ml(v1, nboot = 10000, alpha = 0.05)

