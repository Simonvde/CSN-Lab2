setwd("~/Desktop/Complex-and-Social-Networks/Lab2-Degree-Distribution")
library(igraph)

require("stats4") # for MLE
require("VGAM") # for the Riemann-zeta function

#install.packages("stats4")
#install.packages("VGAM")
library("VGAM")
library("stats4")

# Set a seed so random results will be the same when file is executed again.
set.seed(1)

degree_sequence = read.table("./data/English_in-degree_sequence.txt",
                             header = FALSE)

degree_sequence=read.table("./samples_from_discrete_distributions/data/sample_of_zeta_with_parameter_2.5.txt",header = FALSE)


# Intro -------------------------------------------------------------------


dim(degree_sequence)[1]
nrow(degree_sequence)
sum(degree_sequence)
sum(degree_sequence)/dim(degree_sequence)[1]

degree_spectrum <- table(degree_sequence)
degree_spectrum

barplot(degree_spectrum, main = "English",
        xlab = "degree", ylab = "number of vertices")
barplot(degree_spectrum, main = "English",
        xlab = "degree", ylab = "number of vertices",log="y")



# Zeta_Distribution -------------------------------------------------------

x <- degree_sequence$V1

N <- length(x)
M <- sum(x)
Mp <- sum(log(x))
C <- sum(sapply(x,function(j) sum(log(2:j))))
H <- function(k_max,gamma) sum((1:k_max)^(-gamma))

minus_log_likelihood_zeta <- function(gamma) {
  N * log(zeta(gamma)) + gamma * Mp
}

mle_zeta <- mle(minus_log_likelihood_zeta,
                start = list(gamma = 2),
                method = "L-BFGS-B",
                lower = c(1.0000001))

get_AIC <- function(m2logL,K,N) {
  m2logL + 2*K*N/(N-K-1) # AIC with a correction for sample size
}

summary(mle_zeta)
#AIC(mle_zeta)
#confint(mle_zeta)

attributes(summary(mle_zeta))
attributes(summary(mle_zeta))$coef[1]

m2logL <- attributes(summary(mle_zeta))$m2logL

get_AIC(attributes(summary(mle_zeta))$m2logL, 1, length(x))


minus_log_likelihood_RT_zeta <- function(k_max,gamma) {
  gamma*Mp+N*log(H(k_max,gamma))
}

# Geometric distribution --------------------------------------------------

#Invert formula!
minus_log_likelihood_geom <- function(q) {
  -(sum(x)-length(x))*log(1-q)-length(x)*log(q)
}


q0 <- length(x)/sum(x)

mle_geom <- mle(minus_log_likelihood_geom,
                start = list(q = N/M),
                method = "L-BFGS-B",
                lower = c(0.001),
                upper = c(.999))

summary(mle_geom)



minus_log_likelihood_poisson <- function(lambda) {
  -M*log(lambda)+N*(lambda+log(1-exp(-lambda)))+C
}

mle_poisson <- mle(minus_log_likelihood_poisson,
                start = list(lambda = M/N),
                method = "L-BFGS-B",
                lower = c(0.00001))


# Test on sample data -----------------------------------------------------
get_AIC <- function(m2logL,K,N) {
  print(c(m2logL,2*K*N/(N-K-1)))
  m2logL + 2*K*N/(N-K-1) # AIC with a correction for sample size
}

x <- degree_sequence$V1
AICs <- function(x){
  
  N <- length(x)
  M <- sum(x)
  Mp <- sum(log(x))
  C <- sum(sapply(x,function(j) sum(log(2:j))))
  H <- function(k_max,gamma) sum((1:k_max)^(-gamma))
  
  minus_log_likelihood_zeta <- function(gamma) {
    N * log(zeta(gamma)) + gamma * Mp
  }
  minus_log_likelihood_RT_zeta <- function(k_max,gamma) {
    gamma*Mp+N*log(H(k_max,gamma))
  }
  minus_log_likelihood_geom <- function(q) {
    -(sum(x)-length(x))*log(1-q)-length(x)*log(q)
  }
  minus_log_likelihood_poisson <- function(lambda) {
    -M*log(lambda)+N*(lambda+log(1-exp(-lambda)))+C
  }
  mle_zeta <- mle(minus_log_likelihood_zeta,
                  start = list(gamma = 2),
                  method = "L-BFGS-B",
                  lower = c(1.0000001))
  mle_zeta_G2 <- mle(minus_log_likelihood_zeta,
                  fixed = list(gamma = 2),
                  method = "L-BFGS-B",
                  lower = c(1.0000001))
  mle_RT_zeta <- mle(minus_log_likelihood_RT_zeta,
                  start = list(k_max=max(x),gamma = 2),
                  method = "CG")
  mle_geom <- mle(minus_log_likelihood_geom,
                  start = list(q = N/M),
                  method = "L-BFGS-B",
                  lower = c(0.001),
                  upper = c(.999))
  mle_poisson <- mle(minus_log_likelihood_poisson,
                     start = list(lambda = M/N),
                     method = "L-BFGS-B",
                     lower = c(0.00001))
  att_z <- attributes(summary(mle_zeta))
  att_z_G2 <- attributes(summary(mle_zeta_G2))
  att_rt_z <- attributes(summary(mle_RT_zeta))
  att_geom <- attributes(summary(mle_geom))
  att_pois <- attributes(summary(mle_poisson))
  aics <- c(get_AIC(att_z$m2logL, 1, N),
            get_AIC(att_z_G2$m2logL, 0, N),
            get_AIC(att_rt_z$m2logL, 2, N),
            get_AIC(att_geom$m2logL, 1, N),
            get_AIC(att_pois$m2logL, 1, N))
  coefficients <- c(att_z$coef[1],
                    att_rt_z$coef[1,1],
                    att_rt_z$coef[2,1],
                    att_geom$coef[1],
                    att_pois$coef[1])
  return (aics-min(aics))
}

AICs(degree_sequence$V1)


# Calculating RSS ---------------------------------------------------------

#RT_Zeta
mle_rt_zeta <- mle(minus_log_likelihood_RT_zeta,
                start = list(k_max=max(N),gamma = 2),
                method = "CG")
att_rt_z <- attributes(summary(mle_rt_zeta))
k_max <- att_rt_z$coef[1,1]
gamma_rt <- att_rt_z$coef[2,1]
p_RT_zeta <- function(k) k^(-gamma_rt)/(H(k_max,gamma_rt))

#Zeta
mle_zeta <- mle(minus_log_likelihood_zeta,
                   start = list(gamma = 2),
                   method = "L-BFGS-B",
                   lower=c(1.001))
att_z <- attributes(summary(mle_zeta))
gamma <- att_z$coef[1,1]
p_zeta <- function(k) k^(-gamma)/zeta(gamma)

#Poisson
mle_poisson <- mle(minus_log_likelihood_poisson,
                                start = list(lambda = M/N),
                                method = "L-BFGS-B",
                                lower = c(0.00001))
att_poisson <- attributes(summary(mle_poisson))
lambda <- att_poisson$coef[1,1]
p_poisson <- function(k) lambda^k*exp(lambda)/factorial(k)


s <- sum(degree_sequence)




rowNames <- sort(unique(unlist(degree_sequence)))
degree_spectrum <- table(degree_sequence)
deg_spec_prob <- function(k) degree_spectrum[[k]]/s

RSS <- function(max,FUN){
  arr <- numeric()
  sum <- 0
  i <- 1
  counter=0
  for (n in 1:max){
    if(n<=max(rowNames) && rowNames[i]==n){
      sum <- sum+(deg_spec_prob(which(rowNames==n,arr.ind=TRUE))-FUN(n))^2
      i <- i+1
    }
    else{
      sum <-  sum+(FUN(n))^2
    }
    arr[n] <- sum
  }

  print(sum)
  plot(arr)
}

RSS(3000,p_zeta)
RSS(3000,p_RT_zeta)
RSS(150,p_poisson) #k! is to large to compute values >170. But we already have converge for RSS at 50.
