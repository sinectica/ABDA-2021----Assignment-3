# ASSIGNMENT 3
# APPLIED BAYESIAN DATA ANALYSIS
# AUTHOR: MARC VILA
# DATE: 13 OCTOBER 2021

# TASK A

library("ggplot2")
library("dplyr")
library("stats")

### 1.Use Eq. (4.9) graphically in Table 4.1 to calculate the conditional probability 
### p(blue|red U blond). 

table41 <- cbind(c(0.11,0.03,0.03,0.01,0.18),c(0.20,0.14,0.09,0.05,0.48),c(0.04,0.03,0.02,0.03,0.12),
                c(0.01,0.16,0.02,0.03,0.21), c(0.37,0.36,0.16,0.11,1.00))

Eyes <- c("Black","Brunette","Red","Blond","MarginalEyeColor")
Hair <- c("Brown","Blue","Hazel","Green","MarginalHairColor")

# Create the data frame
df <- expand.grid(Hair, Eyes)
df$value <- c(0.11,0.03,0.03,0.01,0.18,0.20,0.14,0.09,0.05,0.48,0.04,0.03,0.02,0.03,0.12,
                  0.01,0.16,0.02,0.03,0.21,0.37,0.36,0.16,0.11,1.00)

#Plot the Data
g <- ggplot(df, aes(Var1, Var2)) + geom_point(aes(size = value), colour = "green") + theme_bw() + xlab("") + ylab("")
g + scale_size_continuous(range=c(10,30)) + geom_text(aes(label = value))

# Probability p(blue|red U blond)

pblueredUblond <- (table41[15]+table41[20])*(table41[12]+table41[17])/(table41[15]+table41[20])

print(paste("Numerator is",(table41[15]+table41[20])*(table41[12]+table41[17])))
print(paste("Denominator is",(table41[15]+table41[20])))
print(paste("Probability p(blue|red U blond)",pblueredUblond))

### 2. Solve exercise 5.1 using a function that takes an input argument that represents the person's 
### sequential test results and returns the posterior probability of having the disease 
### and not having the disease

P <- vector(mode="list", length=6)
names(P) <- c(":(", ":)", "+|:(", "+|:)", "-|:(", "-|:)")
P[[1]] <- 0.001; P[[2]] <- 0.999; P[[3]] <- 0.99; P[[4]]<-0.05; P[[5]]<-1-0.99; P[[6]]<-1-0.05

P$'+' = P$`+|:)`*P$`:)` + P$`+|:(`*P$`:(`
P$':(|+' = P$'+|:('*P$':(' / P$'+'

P$'-' = P$'-|:)'*P$':)' + P$'-|:('*P$':('
P$':(|-' = P$'-|:('*P$':(' / P$'-'

post <- function (T) {
  
# Initialize prior and likelihood vectors 
prior <- c(P$':(', P$':)')
likelihood = c(1.0,1.0)
posterior = as.double(likelihood)*prior/(sum(as.double(likelihood)*prior)) # Bayes formula with no data, so posterior = prior

T <- unlist(strsplit(T, split="")) # Splitting string T into a vector to be able to work with


# The algorithm works sequentially, starting from the prior probabilities (i.e. without any test results) and then moving 
# on sequentially with access to new test results.

for (i in T) {

  likelihood <- c( P[paste(i,'|:(',sep = "")], P[paste(i,'|:)',sep = "")]) # get likelihood for result in t given sick or healthy
  posterior <- as.double(likelihood)*prior/(sum(as.double(likelihood)*prior)) # posterior calculus
  prior = posterior # update prior for the next experiment
}

names(posterior) <- c(':(|T', ':)|T')
print(posterior)

}

