setwd("~/Desktop/CSN-Labs/Lab2/")
library(igraph)

library("VGAM")
library("stats4")

# Set a seed so random results will be the same when file is executed again.
set.seed(1)

degree_sequence <- read.table("./data/Catalan_in-degree_sequence.txt",
                              header = FALSE)

# Test on sample data -----------------------------------------------------
get_AIC <- function(m2logL,K,N) {
  #print(c(m2logL,2*K*N/(N-K-1)))
  m2logL + 2*K*N/(N-K-1) # AIC with a correction for sample size
}

x <- degree_sequence$V1
AICs <- function(x,graph_=F){
  
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
    -(M-N)*log(1-q)-N*log(q)
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
  if(graph_){
    return (coefficients)
  }
  return (aics-min(aics))
}

AICs(degree_sequence$V1)



# Calculating RSS ---------------------------------------------------------

RSS_full <- function(x,zeta_rt=F){
  pars<-AICs(x,T)
  gamma<-pars[1]
  p_zeta <- function(k) k^(-gamma)/zeta(gamma)
  k_max<-pars[2]
  gamma_rt<-pars[3]
  H <- function(k_max,gamma) sum((1:k_max)^(-gamma))
  p_RT_zeta <- function(k) k^(-gamma_rt)/(H(k_max,gamma_rt))
  
  s <- length(x)
  
  rowNames <- sort(unique(unlist(x)))
  degree_spectrum <- table(x)
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

    return (sum)
  }
  if(zeta_rt) return(RSS(max(x),p_RT_zeta))
  else return(RSS(max(x),p_zeta))
}
RSS_full(x)
RSS(3000,p_zeta)
RSS(3000,p_RT_zeta)
RSS(150,p_poisson) #k! is to large to compute values >170. But we already have converge for RSS at 50.


# Output ------------------------------------------------------------------

source = read.table("list.txt", 
                    header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
                    as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
)

aic_diff <- matrix(nrow=length(source$language),ncol=5)
rss <- matrix(nrow=length(source$language),ncol=2)

write_RSS <- function(language,file){
  degree_sequence = read.table(file, header = FALSE)
  return(c(RSS_full(degree_sequence$V1),RSS_full(degree_sequence$V1,T)))
}

write_AICs <- function(language,file) {
  degree_sequence = read.table(file, header = FALSE)
  return(AICs(degree_sequence$V1))
}

for (i in 1:nrow(source)) aic_diff[i,] <- write_AICs(source$language[i], source$file[i])
for (i in 1:nrow(source)) rss[i,] <- write_RSS(source$language[i], source$file[i])

dimnames(aic_diff) <- list(source$language, c("zeta","zeta_2","RT_zeta","geom","poisson"))
aic_diff
library(xtable)
xtable(aic_diff)

dimnames(rss) <- list(source$language, c("zeta","RT_zeta"))
rss
library(xtable)
xtable(rss,display=rep("e",3))


# Graphical ---------------------------------------------------------------

#Zeta


plot(table(x)/length(x),xlim=c(-14,14))
curve(p_zeta,add=T,col='red')

#All models
pars<-AICs(x,T)
gamma<-pars[1]
p_zeta <- function(k) k^(-gamma)/zeta(gamma)
k_max<-pars[2]
gamma_rt<-pars[3]
H <- function(k_max,gamma) sum((1:k_max)^(-gamma))
p_RT_zeta <- function(k) k^(-gamma_rt)/(H(k_max,gamma_rt))
lambda<-pars[5]
p_poisson <- function(k) lambda^k*exp(-lambda)/(factorial(k)*(1-exp(-lambda)))
prob<-pars[4]
p_geom <- function(k) prob*(1-prob)^(k-1)

plot(table(x)/length(x),xlim=c(0,27))
curve(p_zeta,add=T,col='red',lty=2)
curve(p_poisson,add = T,col='blue')
points(1:27,p_geom(1:27),add = T,col='green')
curve(p_RT_zeta,add = T,col='yellow',lty=2,lwd=3)

curve(p_RT_zeta,add=T,col='purple')

plot(table(x)/length(x),xlim=c(100,200))

p_zeta(100)
p_RT_zeta(10000)
curve(dgeom)
