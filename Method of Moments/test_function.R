# Function to test the estimated parameters for each distribution. 
#100 samples are generated for a particular distribution.

testFunction <- function(dist)
{
  num_samples <- 100 # Number of samples generated
  
  if(dist=="bernoulli")
  {
    data <- rbinom(num_samples,1,0.6) # p = 0.6
    p <- myfunction(dist,data)
    return(p)
  }
  else if(dist=="binomial")
  {
    data <- rbinom(num_samples,10,0.6) # n = 10, p = 0.6
    theta <- myfunction(dist,data)
    return(theta)
  }
  else if(dist=="geometric")
  {
    data <- rgeom(num_samples,0.6) # p = 0.6
    p <- myfunction(dist,data)
    return(p)
  }
  else if(dist=="poisson")
  {
    data <- rpois(num_samples,0.02) # lambda = 0.02
    lamda <- myfunction(dist,data)
    return(lamda)
  }
  else if(dist=="uniform")
  {
    data <- runif(num_samples,0,5) # a = 0, b = 5
    theta <- myfunction(dist,data)
    return(theta)
  }
  else if(dist=="normal")
  {
    data <- rnorm(num_samples,0,1) # mean = 0, variance = 1
    theta <- myfunction(dist,data)
    return(theta)
  }
  else if(dist=="exponential")
  {
    data <- rexp(num_samples,6) # beta = 6
    b <- myfunction(dist,data)
    return(b)
  }
  else if(dist=="gamma")
  {
    data <- rgamma(num_samples,2,scale = 0.4) # alpha = 2, beta = 0.4
    theta <- myfunction(dist,data)
    return(theta)
  }
  else if(dist=="beta")
  {
    data <- rbeta(num_samples, 3, 5) # alpha = 3, beta = 5
    theta <- myfunction(dist,data)
    return(theta)
  }
  else if(dist=="Chi-square")
  {
    data <- rchisq(num_samples,3) # degree of freedom (p) = 3
    p <- myfunction(dist,data)
    return(p)
  }
  else
  {
    return("Wrong distribution")
  }
}

# Test by passing distribution as parameter
res <- testFunction("Chi-square")
