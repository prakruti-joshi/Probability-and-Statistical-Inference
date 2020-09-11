library(tidyverse)
library(dplyr)

# Jacknife Estimator
Jackknife_confidence_interval <- function( v1, statfunc=sd, alpha = 0.05 ) {
  
  n1 <- length(v1)
  jackvec <- NULL
  mu0 <- statfunc(v1)  # Original value from the data
  
  for(i in 1:n1)
  {
    mua<-statfunc(v1[-i])                    # Function value by removing the ith data
    jackvec<-c(jackvec, n1*(mu0)-(n1-1)*mua) # Jacknife bias corrected term
    
  }
  jackbias<-mean(jackvec)-mu0   # Bias calculation
  jacksd<-sd(jackvec)           # Standard deviation for jacknife estimation
  
  # Confidence interval
  NLB<-mean(v1)-(sd(v1)/sqrt(n1))*qt(1-alpha/2, df = n1-1) 
  NUB<-mean(v1)+(sd(v1)/sqrt(n1))*qt(1-alpha/2, df = n1-1)
  
  # Jackknife result 
  list(mu0=mu0,jackbias=jackbias,jacksd=jacksd, normal.confidence.interval=c(NLB,NUB)) 
}

# Bootstrap Estimation
my.bootstrapci.ml<-function(vec0, nboot=10000, alpha=0.1){
  
  # Sample size, mean and standard deviation from the original data
  
  n0 <- length(vec0) 
  mean0 <- mean(vec0) 
  sd0 <- sqrt(var(vec0))
  
  # create a vector to store the location of the bootstrap studentized deviation vector 
  bootvec<-NULL
  bootbiasvec<-NULL
  
  #create the bootstrap distribution using a for loop
  for( i in 1:nboot){
    vecb<-sample(vec0,replace=T) # samples from original data with replacement
    meanb<-mean(vecb)
    sdb<-sqrt(var(vecb))
    
    #note since resampling full vector we can use n0 for sample size of vecb
    bootvec<-c(bootvec,(meanb-mean0)/(sdb/sqrt(n0))) # Bias Correction by subtracting 1st moment and dividing by 2d moment
    bootbiasvec<-c(bootbiasvec,meanb-mean0)
  }
 
  bootsd <- sd(bootvec)   # Bootstrap estimation of parameter
  bootbias<-mean(bootbiasvec)
  
  # Calculate lower and upper quantile of the bootstrap distribution
  
  lq<-quantile(bootvec,alpha/2) 
  uq<-quantile(bootvec,1-alpha/2)
  
  LB<-mean0-(sd0/sqrt(n0))*uq
  UB<-mean0-(sd0/sqrt(n0))*lq

  # Normal confidence interval 
  NLB<-mean0-(sd0/sqrt(n0))*qnorm(1-alpha/2) 
  NUB<-mean0+(sd0/sqrt(n0))*qnorm(1-alpha/2)
  
  # t- distribution confidence interval 
  TLB <- mean0-(sd0/sqrt(n0))*qt(1-alpha/2, df = n0) 
  TUB<- mean0+(sd0/sqrt(n0))*qt(1-alpha/2, df = n0)
  
  # Bootstrap result
  list( bootbias = bootbias, bootsd = bootsd, bootstrap.confidence.interval=c(LB,UB),normal.confidence.interval=c(NLB,NUB), t.confidence.interval = c(TLB, TUB) ) 
}

Sim.func<-function(mu.val=3,n=30,nsim=1000)
{
  #create coverage indicator vectors for bootstrap and normal 
  cvec.boot<-NULL
  cvec.norm<-NULL
  cvec.t <- NULL
  
  #calculate real mean
  mulnorm<-(exp(mu.val+1/2)) #taking standard normal
  
  #run simulation
  
  for(i in 1:nsim){
    if((i/10)==floor(i/10)){ 
      print(i)
      #let me know computer hasnt died
    }
    
    #sample the simulation vector 
    vec.sample<-rlnorm(n,mu.val)
    
    #bootstrap 
    boot.list<-my.bootstrapci.ml(vec.sample)           # Running bootstrap function and storing results
    boot.conf<-boot.list$bootstrap.confidence.interval # Bootstrp confidence interval
    norm.conf<-boot.list$normal.confidence.interval    # Normal confidence interval
    t.conf <- boot.list$t.confidence.interval          # T-dist confidence interval
    
    # Checking if estimated value lies in the confidence intervals:
    
    ## bootstap CI
    cvec.boot<-c(cvec.boot, (boot.conf[1]<mulnorm) * (boot.conf[2]>mulnorm) )  
    ## Normal CI
      # Note: no need to call Jackknife since normal of both will be the same
    cvec.norm<-c(cvec.norm,(norm.conf[1]<mulnorm)*(norm.conf[2]>mulnorm)) 
    ## t CI
    cvec.t<-c(cvec.t,(t.conf[1]<mulnorm)*(t.conf[2]>mulnorm)) 
    
  }
  # calculate and output coverage probability estimates 
  list(boot.coverage=(sum(cvec.boot)/nsim), norm.coverage=(sum(cvec.norm)/nsim), t.coverage=(sum(cvec.t)/nsim)) 
}

# Fucntion call
Sim.func(n=30, nsim=1000)

############ Question:3 Comparison of coverage rates for different confidence intervals ############


# dataframe to store results
df <-  data.frame(n_sample = numeric(), boot.cov = numeric(), norm.cov = numeric(), t.cov = numeric()) 
# Simulating for different sample sizes:
sample_size = c(10,30,100)
nsim = 1000
for (n_sample in sample_size)
{
  coverage.list <- Sim.func(mu.val=3,n=n_sample,nsim=nsim)
  boot.cov<-coverage.list$boot.coverage
  norm.cov<-coverage.list$norm.coverage
  t.cov <- coverage.list$t.coverage
  
  l <- list(n_sample = n_sample,boot.cov = boot.cov, norm.cov = norm.cov, t.cov=t.cov )
  df <- rbind(df, l)
}

# For plotting the coverage rate results for different confidence intervals

df %>%
  gather(key = "variable", value = "value", -n_sample) %>%
  ggplot(aes(n_sample, value, col = variable)) +
  geom_line( position=position_dodge(0.5)) +
  geom_point(position=position_dodge(0.5), size=2) +
  scale_colour_manual(name="Confidence Interval", labels = c("Bootstrap","CLT based", "Jackknife Normal \n Approximation (t distibution)"), values =  c( "#ee5253", "#2e86de", "#10ac84")) +
  labs(y = "Coverage Rate",
       x = "Sample size", 
       subtitle = paste0('#simulations = ',nsim) ) + 
  theme_bw()

