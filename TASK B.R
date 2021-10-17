# ASSIGNMENT 3
# APPLIED BAYESIAN DATA ANALYSIS
# AUTHOR: MARC VILA
# DATE: 15 OCTOBER 2021

# TASK B

library("tidyverse")
library("dplyr")
library("Rlab")

library("stats")
library("ggplot2")

source("BernGrid.R")
source("DBDA2E-utilities.R")

# We will start by reproducing the "inferring bias in a coin" example in the book (Ch 5 & 6).
# As pointed out by the author, this is a toy example but the same model applies to many 
# different situations that you may be faced with (see e.g. the analogies given by the 
# author on page 73 and 124 in the book).

### 1. The probability of each outcome is given by the Bernoulli 
### distribution Eq. (6.1).

### 2. It is a discrete distribution of the outcome y=0 or y=1 for a fixed theta.
### Plot the outcome probabilities given  theta = 0.25 and theta = 0.5 and using a bar or a stem plot.

  # Sample size
N <- 1000
theta1 = 0.25
theta2 = 0.5

  # Generate random variables using rbern( ) function in R, setting the seed for reproducibility
set.seed(98999)  
random_values <- rbern(N, prob = 0.25)  

  # plot of randomly drawn density
hist(random_values,breaks = 10, ylim = c(0,1000), col = "blue", border = "pink", main = "")

### 3. Given a certain outcome, say y=1, we can plot the likelihood of our parameter theta. 
### Note that when evaluating (6.1) with respect to theta for a fixed outcome Y=1, 
### it is called a likelihood function and it is not a probability distribution 
### since it will not integrate to 1. For y=0 and y=1 plot the likelihood function for theta E [0,1]

  # Bernoulli probability density function
Bernoulli_pdf <- function (gamma, theta) {
  bpdf <- (theta^gamma)*((1.0-theta)^(1.0-gamma))  
}

  # Generate sequence of theta values
theta <- seq(0, 1, length.out = 1000)

  # 2 gamma (y) values
y1 <- 1
y0 <- 0

  # Generating both likelihood functions
likelihood1 <- Bernoulli_pdf(y1,theta)
likelihood0 <- Bernoulli_pdf(y0,theta)

  # Plot of both functions
plot(x=theta, y=likelihood1,col = "blue", border = "pink")
lines(x=theta, y=likelihood0, col = "red")

plot( theta, likelihood0, type="p", col="red" )
par(new=TRUE)
plot( x=theta, likelihood1, type="p", col="green" )

### 4. Now we will assume independence, row one in Eq. (6.2), and that the coin flips 
### are independent of each other.
  
## a. Implement 2 functions:

  # i) The first used in Eq. 6.1

likelihood <- Bernoulli_pdf(y1,theta)

  # ii) The other function is just the product of the first function over the data 
  #     as the first row in Eq. (6.2)

likelihood <- function (gamma, theta) {
  product_likelihood <- prod(Bernoulli_pdf(gamma, theta))
}

  # Create vectors of random sample of n flips (heads=1, tails=0):

set.seed(12345) # Set seed
n <- c(10,100,1000,10000) # Vector with different number of flips.
pHeads = 0.25 # Probability of heads

  flipSequence1 = sample( x=c(0,1), prob=c(1-pHeads,pHeads), size=n[1], replace=TRUE)
  flipSequence2 = sample( x=c(0,1), prob=c(1-pHeads,pHeads), size=n[2], replace=TRUE)
  flipSequence3 = sample( x=c(0,1), prob=c(1-pHeads,pHeads), size=n[3], replace=TRUE)
  flipSequence4 = sample( x=c(0,1), prob=c(1-pHeads,pHeads), size=n[4], replace=TRUE)

  # Likelihood as a product calculation with theta = 0.5
  
lk1 <- likelihood(flipSequence1,0.5)
lk2 <- likelihood(flipSequence2,0.5)
lk3 <- likelihood(flipSequence3,0.5)
lk4 <- likelihood(flipSequence4,0.5)

print(paste("Likelihood for n = ", n[1], " flips is ",lk1))
print(paste("Likelihood for n = ", n[2], " flips is ",lk2))
print(paste("Likelihood for n = ", n[3], " flips is ",lk3))
print(paste("Likelihood for n = ", n[4], " flips is ",lk4))

## b. Implement the logarithm versions of these two functions,
## i.e. the log-likelihood and the log-pdf, using the trick of
## logarithm and summation instead. 

  # pdf log function and sum (equivalent to the product of the non log pdf distribution)

log_Bernoulli_pdf <- function (gamma, theta) {
  log_bpdf <- gamma*log(theta) + ((1.0-gamma)*log(1.0-theta))
}

log_likelihood <- function (gamma, theta) {
  sum_likelihood <- sum(log_Bernoulli_pdf(gamma, theta))
}

### Can you evaluate the log-likelihood for larger  without problems of under- or overflow?

log_lk1 <- log_likelihood(flipSequence1,0.5)
log_lk2 <- log_likelihood(flipSequence2,0.5)
log_lk3 <- log_likelihood(flipSequence3,0.5)
log_lk4 <- log_likelihood(flipSequence4,0.5)

print(paste("Log-Likelihood for n = ", n[1], " flips is ",log_lk1))
print(paste("Log-Likelihood for n = ", n[2], " flips is ",log_lk2))
print(paste("Log-Likelihood for n = ", n[3], " flips is ",log_lk3))
print(paste("Log-Likelihood for n = ", n[4], " flips is ",log_lk4))

## c. At the final step, if absolutely necessary, one can exponentiate
## the log-likelihood function to obtain the likelihood value. However, 
## for this particular example it will still be problematic as 
## exponentiating large negative values will result in zero.

print(paste("Log-Likelihood for n = ", n[1], " flips is ",exp(log_lk1)))
print(paste("Log-Likelihood for n = ", n[2], " flips is ",exp(log_lk2)))
print(paste("Log-Likelihood for n = ", n[3], " flips is ",exp(log_lk3)))
print(paste("Log-Likelihood for n = ", n[4], " flips is ",exp(log_lk4)))

## d. Plot the likelihood function (use either the likelihood or preferable
## exp(log_likelihood) with respect to theta E [1,0] given:
## y=[1], y=[1,1] and y=[1,1,0,1]

y <- list(c(1), c(1,1), c(1,1,0,1)) # List of vectors of different sequences of y
l <- seq(0.0, 0.0, length.out=length(theta)) # Initializing a vector of the same length as seq of theta values

# Calculating th likelihood values

  for(x in y){
    for (i in 1:length(theta)){
      l[i] <- likelihood(x,theta[i])
    }
    plot(x=theta, y=l, type="l", col="blue")
    par(new=TRUE)
  }
dev.off()

### 5. Convert the likelihood of theta to a probability density of theta
### given the observed values of y.

## a. Apply Eq. (5.7) Bayes rule

  # i. beta distribution with a=1 and b=1
y_beta <- dbeta(theta, shape1 = 1, shape2 = 1)

  # ii. Use posterior distribution in Eq. (6.8) and make new plots of data
  # sets of y. Implement the log_posterior and at the very utmost expression
  # take the exponential

z1 <-1
z2 <- 2
z3 <- 3
N1 <- 1
N2 <- 2
N3 <- 4
a <- 1
b <- 1

# The posterior distribution is a beta(theta|z+a,N-z+b) distribution

posterior_dbf1 <- dbeta (theta, shape1=z1+a, shape2=N1-z1+b)
posterior_dbf2 <- dbeta (theta, shape1=z2+a, shape2=N2-z2+b)
posterior_dbf3 <- dbeta (theta, shape1=z3+a, shape2=N3-z3+b)

# Plot the three distributions

plot(x=theta, y=posterior_dbf1, type="l", col="blue")
lines(x=theta, y=posterior_dbf2, type="l", col="red")
lines(x=theta, y=posterior_dbf3, type="l", col="green")

# Plot the Log of the three distributions

plot(x=theta, y=log(posterior_dbf1),type="l", col="blue")
lines(x=theta, y=log(posterior_dbf2), type="l", col="red")
lines(x=theta, y=log(posterior_dbf3),type="l", col="green")

### Questions 1, 2 and 3

  ## iii. Assume now more informative prior information and reproduce
  ## Figure 6.4 in the book. Note that the last three figures in 
  ## column 3 you have already managed to work out from the above tasks. 


# As done in the book, with book's code and libraries

source("DBDA2E-utilities.R") # Load definitions of graphics functions etc.
source("BernBeta.R") # Load the definition of the BernBeta function

# Specify the prior:
t = 0.5 # Specify the prior mode.
n = 200 # Specify the effective prior sample size.
a = t*(n-2) + 1 # Convert to beta shape parameter a.
b = (1-t)*(n-2) + 1 # Convert to beta shape parameter b.
Prior = c(a,b) # Specify Prior as vector with the two shape parameters

N = 20 # The total number of flips.
z = 17 # The number of heads.

Data = c(rep(0,N-z),rep(1,z)) # Convert N and z into vector of 0's and 1's.

openGraph(width=5,height=7)

posterior = BernBeta( priorBetaAB=Prior, Data=Data , plotType="Bars" ,
                      showCentTend="Mode" , showHDI=TRUE , showpD=FALSE )
saveGraph(file="BernBetaExample",type="jpg")
dev.off()

  ## My code

# Specify the prior (it will change depending on the example)
t = 0.5 # Specify the prior mode.
n = 200 # Specify the effective prior sample size.
a = t*(n-2) + 1 # Convert to beta shape parameter a.
b = (1-t)*(n-2) + 1 # Convert to beta shape parameter b.

N = 20 # The total number of flips.
z = 17 # The number of heads.

# The prior is a beta distribution with a and b parameters, performed with dbeta function
prior_Beta <- dbeta(theta,shape1=a, shape2=b)

# Likelihood is a Bernoulli likelihood distribution
likelihood_Bernoulli <- (theta^z)*((1.0-theta)^(N-z))

# Posterior is a beta distribution, like prior but with modified shape (depending on N and z)
posterior_Beta <- dbeta (theta, shape1=z+a, shape2=N-z+b)

  ## 3 plots, prior, likelihood and posterior

par(mfrow=c(3,1))
plot(theta, prior_Beta, type="l", col="red")
plot(theta, likelihood_Bernoulli, type="l", col="blue")
plot(theta, posterior_Beta, type="l", col="green")
dev.off

### b. Read summary in section 6.5 carefully and discuss




