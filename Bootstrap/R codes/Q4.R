# Jackknife estimation
Jackknife_confidence_interval <- function( v1, statfunc=sd, alpha = 0.05 ) {
  
  n1 <- length(v1)
  jackvec <- NULL
  mu0 <- statfunc(v1)  # Original parameter/statfunc value from the data
  for(i in 1:n1)
  {
    mua<-statfunc(v1[-i])                    # Function value by removing the ith data
    jackvec<-c(jackvec, n1*(mu0)-(n1-1)*mua) # Jacknife bias corrected term
    
  }
  jackbias<-mean(jackvec)-mu0 # Bias of the estimator
  jacksd<-sd(jackvec)         # Standard deviation for jacknife estimation
  
  # Confidence interval
  NLB<-sd(v1)-(sd(v1)/sqrt(n1))*qt(1-alpha/2, df = n1-1) 
  NUB<-sd(v1)+(sd(v1)/sqrt(n1))*qt(1-alpha/2, df = n1-1)
  
  # Results of jackknife
  list(jackestimate = mean(jackvec), jackbias=jackbias,jacksd=jacksd, normal.confidence.interval=c(NLB,NUB)) 
}

## Boostrap estimation
my.bootstrapci.ml<-function(vec0, nboot=10000, alpha=0.1){
  
  # Sample size, mean and standard deviation from the original data
  n0 <- length(vec0) 
  mean0 <- mean(vec0) 
  sd0 <- sqrt(var(vec0))
  
  # create a vector to store the location of the bootstrap studentized deviation vector 
  bootvec<-NULL
  bootbiasvec<-NULL
  jacckbiasvec<-NULL
  jacksdvec<-NULL
  
  #create the bootstrap distribution using a for loop
  for( i in 1:nboot){
    vecb<-sample(vec0,replace=T)   # samples from original data with replacement
    meanb<-mean(vecb)
    sdb<-sqrt(var(vecb))
    
    #note since resampling full vector we can use n0 for sample size of vecb
    bootvec<-c(bootvec,(meanb-mean0)/(sdb/sqrt(n0)))   # Bias Correction by subtracting 1st moment and dividing by 2d moment
    bootbiasvec<-c(bootbiasvec,meanb-mean0)
  }
  hist(bootvec)
  bootsd <- sd(bootvec)   # Boostrap estimation of standard deviation
  bootbias<-mean(bootbiasvec)
  
  # Normal CI 
  NLB<-mean0-(sd0/sqrt(n0))*qnorm(1-alpha/2) 
  NUB<-mean0+(sd0/sqrt(n0))*qnorm(1-alpha/2)
  
  # T-dist CT
  TLB <- mean0-(sd0/sqrt(n0))*qt(1-alpha/2, df = n0) 
  TUB<- mean0+(sd0/sqrt(n0))*qt(1-alpha/2, df = n0)
  
  list(bootestimate = mean(bootvec),   bootbias = bootbias, bootsd = bootsd, bootstrap.confidence.interval=c(LB,UB),normal.confidence.interval=c(NLB,NUB), t.normal.confidence.interval = c(TLB, TUB) ) 
}

## estimation of SD:
my.bootstrapci.sd<-function(vec0, nboot=10000, alpha=0.1){
  #extract sample size, mean and standard deviation from the original data
  n0 <- length(vec0) 
  mean0 <- mean(vec0) 
  sd0 <- sqrt(var(vec0))   # StatFunction 
  
  # create a vector to store the location of the bootstrap studentized deviation vector 
  bootvec<-NULL
  bootbiasvec<-NULL
  jacckbiasvec<-NULL
  jacksdvec<-NULL
  
  #create the bootstrap distribution using a for loop
  for( i in 1:nboot){
    vecb<-sample(vec0,replace=T) #create mean and standard deviation to studentize
    meanb<-mean(vecb)
    sdb<-sqrt(var(vecb))
    #note since resampling full vector we can use n0 for sample size of vecb
    bootvec<-c(bootvec, sdb/sqrt(n0)) 
    bootbiasvec<-c( bootbiasvec, sdb-sd0 )
  }
  bootsd <- sd(bootvec)
  bootbias<-mean(bootbiasvec)
  
  list(bootestimate = sqrt(n0) * mean(bootvec),   bootbias = bootbias, bootsd = bootsd, bootstrap.confidence.interval=c(LB,UB) ) 
}

## simulation function
Sim.func<-function(mu.val=3,n=30,nsim=1000)
{
  #create coverage indicator vectors for bootstrap and normal 
  boot.estimate<-NULL
  boot.bias<-NULL
  boot.sd <- NULL
  cvec.boot<-NULL
  jack.sd <- NULL
  cvec.jack <- NULL
  
  #calculate real mean
  mulnorm <- 1 #taking standard normal
  
  # run simulation
  for(i in 1:nsim){
    if((i/10)==floor(i/10)){ 
      print(i)
    }
    #sample the simulation vector 
    vec.sample <- rnorm(n)
    #bootstrap it
    boot.list<-my.bootstrapci.sd(vec.sample,alpha = 0.05) 
    boot.estimate<- c(boot.estimate, boot.list$bootestimate) 
    boot.bias<-c(boot.bias, boot.list$bootbias)
    boot.sd <- c(boot.sd, boot.list$bootsd)
    boot.conf<-boot.list$bootstrap.confidence.interval 
    
    jack.list <- Jackknife_confidence_interval(vec.sample, sd)
    jack.sd <- c(jack.sd, jack.list$jacksd)
    jack.conf <- jack.list$normal.confidence.interval
    
    cvec.boot<-c(cvec.boot, (boot.conf[1]<mulnorm) * (boot.conf[2]>mulnorm) )  
    cvec.jack <- c(cvec.jack,  (jack.conf[1]<mulnorm) * (jack.conf[2]>mulnorm) )
    
  }
  list(avg.boot.estimate = (sum(boot.estimate)/nsim) , avg.boot.bias = (sum(boot.bias)/nsim), avg.boot.sd = (sum(boot.sd)/nsim),
       avg.jack.sd = (sum(jack.sd)/nsim), boot.coverage=(sum(cvec.boot)/nsim), jack.coverage=(sum(cvec.jack)/nsim), boot.sd = boot.sd, jack.sd = jack.sd) 
}

out <- Sim.func(nsim=1000)
hist(out$ boot.sd)
hist(out$ jack.sd)

# Boostrap and jackknife for estimation of standard deviation of normal distribution
v1 <- rnorm(1000)
Jackknife_confidence_interval(v1,sd)
my.bootstrapci.sd(v1, nboot = 1000,alpha=0.05)
