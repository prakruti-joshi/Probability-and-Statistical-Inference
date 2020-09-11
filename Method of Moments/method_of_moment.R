# Function to estimate paramters of a distribution from sample data using method of moments estimator.

myfunction<-function(dist,data)
{
  sum_total <- sum(data,na.rm=TRUE)
  n <- length(data)
  alpha_1 <- sum_total/n
  
  square_total <- sum((data-alpha_1)^2,na.rm=TRUE)
  alpha_2 <-square_total/n
  
  if(dist=="point_mass")
  {
    a <- alpha_1
    return(a)
  }
  else if(dist=="bernoulli")
  {
    p <- alpha_1
    return(p)
  }
  else if(dist=="binomial")
  {
    m <- alpha_1/(1-(alpha_2/alpha_1))
    p <- 1-(alpha_2/alpha_1)
    theta <- c(m,p)
    return(theta)
  }
  else if(dist=="geometric")
  {
    p <- (1/(alpha_1))-1
    return(p)
  }
  else if(dist=="poisson")
  {
    lamda <- alpha_1
    return(lamda)
  }
  else if(dist=="uniform")
  {
    a <- alpha_1 - sqrt(3*alpha_2)
    b <- alpha_1 + sqrt(3*alpha_2)
    theta <- c(a,b)
    return(theta)
  }
  else if(dist=="normal")
  {
    mean <- alpha_1
    variance <- alpha_2
    theta <- c(mean,variance)
    return(theta)
  }
  else if(dist=="exponential")
  {
    b <- 1/alpha_1
    return(b)
  }
  else if(dist=="gamma")
  {
    a <- (alpha_1^2)/(alpha_2)
    b <- (alpha_2)/alpha_1
    theta <- c(a,b)
    return(theta)
  }
  else if(dist=="beta")
  {
    a <- (alpha_1^2)*(1-alpha_1)/alpha_2 - alpha_1
    b <- (((alpha_1/alpha_2)*(1-alpha_1))-1)*(1-alpha_1)
    theta <-c(a,b)
    return(theta)
  }
  else if(dist=="Chi-square")
  {
    p <- alpha_1
    return(p)
  }
  else
  {
    return("Wrong parameters")
  }
}