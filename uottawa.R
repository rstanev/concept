#########################################################################
# author: 	roger stanev 
# date: 	may 15, 2015
# version:	v 0.3
##########################################################################
# a)

library(LearnBayes) 

# Function to compute log of the posterior density of hyper 
# parameters mu and log thau
logpost <- function(theta, data){ 
  y <- data[, 1] 
  sigma <- data[, 2] 
  mu <- theta[1] 
  tau <- exp(theta[2]) 
  val <- sum(dnorm(y, mean = mu, sd = (sigma^2 + tau^2)^0.5, log = T)) + log(tau) 
  return(val) 
}

# Initializing data
data <- data.frame(y = c(28, 8, -3, 7, -1, 1, 18, 12), 
    sigma = c(15, 10, 16, 11, 9, 11, 10, 18)) 

fit <- laplace(logpost, c(2, 3), data) 
mycontour(logpost, c(-20, 40, -7, 4), data) 
r <- gibbs(logpost, start = c(10, 2), 1000, c(8, 1.4), data) 
points(r[[1]])

# b)

J <- 8 
y <- data[, 1] 
sigma <- data[, 2] 
probint <- NULL 

for(j in 1:J) { 
  tau <- exp(r[[1]][, 2]) 
  mu <- r[[1]][, 1] 
  thetaj <- (y[j]/sigma[j]^2+mu/tau^2)/(1/sigma[j]^2+1/tau^2) 
  Vj <- 1/(1/sigma[j]^2+1/tau^2) 
  probint <- cbind(probint, rnorm(1000, thetaj, Vj^0.5)) 
}

apply(probint, 2, median) 
apply(probint, 2, mean) 
apply(probint, 2, sd) 

# question 3 

# part (a) 
logpost <- function(theta, data){ 
  y <- data[, 1] 
  sigma <- data[, 2] 
  mu <- theta[1] 
  tau <- exp(theta[2]) 
  val <- sum(dnorm(y, mean = mu, sd = (sigma^2 + tau^2)^0.5, log = T)) + log(tau) 
  return(val) 
} 

data <- data.frame(y = c(28, 8, -3, 7, -1, 1, 18, 12), 
                   sigma = c(15, 10, 16, 11, 9, 11, 10, 18)) 

fit <- laplace(logpost, c(2, 3), data) 
mycontour(logpost, c(-20, 40, -7, 4), data) 
r <- gibbs(logpost, start = c(10, 2), 1000, c(8, 1.4), data) 
points(r[[1]]) 

J <- 8 
y <- data[, 1] 
sigma <- data[, 2] 
probint <- NULL 
thetaj <- NULL 

for(j in 1:J) { 
  tau <- exp(r[[1]][, 2]) 
  mu <- r[[1]][, 1] 
  thetaj <- cbind(thetaj, (y[j]/sigma[j]^2+mu/tau^2)/(1/sigma[j]^2+1/tau^2)) 
  Vj <- 1/(1/sigma[j]^2+1/tau^2) 
  probint <- cbind(probint, rnorm(10000, thetaj[,j], Vj^0.5)) 
} 
apply(probint, 2, median) 
apply(probint, 2, mean) 
apply(probint, 2, sd) 
B <- NULL 
for (j in 1:J) { 
  B <- cbind(B, tau^(-2)/( tau^(-2)+sigma[j]^(-2))) 
} 
B.e <- apply(B, 2, mean) 
sort(B.e) 
order(B.e) 

# part (b) 
sum(apply(thetaj, 1, which.max) ==1)/1000
