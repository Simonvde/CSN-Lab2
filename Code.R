setwd("~/GitHub/CSN-Lab2")

library(igraph)

library("VGAM")
library("stats4")

# Set a seed so random results will be the same when file is executed again.
set.seed(1)

degree_sequence <- read.table("./data/Catalan_in-degree_sequence.txt",
                              header = FALSE)

# Akaike no_Altmann -----------------------------------------------------
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
  minus_log_likelihood_geom_corrected <- function(q) {
    x_1=x-1
    N <- length(x_1)
    M <- sum(x_1)
    -(M-N)*log(1-q)-N*log(q)
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
  mle_geom_corrected <- mle(minus_log_likelihood_geom_corrected,
                  start = list(q = N/M),
                  method = "L-BFGS-B",
                  lower = c(0.001),
                  upper = c(.999))
  att_z <- attributes(summary(mle_zeta))
  att_z_G2 <- attributes(summary(mle_zeta_G2))
  att_rt_z <- attributes(summary(mle_RT_zeta))
  att_geom <- attributes(summary(mle_geom))
  att_pois <- attributes(summary(mle_poisson))
  att_geom_corrected <- attributes(summary(mle_geom_corrected))
  aics <- c(get_AIC(att_z$m2logL, 1, N),
            get_AIC(att_z_G2$m2logL, 0, N),
            get_AIC(att_rt_z$m2logL, 2, N),
            get_AIC(att_geom$m2logL, 1, N),
            get_AIC(att_pois$m2logL, 1, N),
            get_AIC(att_geom_corrected$m2logL, 1, N)
            )
  coefficients <- c(att_z$coef[1],
                    att_rt_z$coef[1,1],
                    att_rt_z$coef[2,1],
                    att_geom$coef[1],
                    att_pois$coef[1],
                    att_geom_corrected$coef[1])
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

# Output ------------------------------------------------------------------

source = read.table("list.txt", 
                    header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
                    as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
)

aic_diff <- matrix(nrow=length(source$language),ncol=6)
rss <- matrix(nrow=length(source$language),ncol=2)

write_RSS <- function(language,file){
  degree_sequence = read.table(file, header = FALSE)
  return(c(RSS_full(degree_sequence$V1),RSS_full(degree_sequence$V1,T)))
}

write_AICs <- function(language,file) {
  degree_sequence = read.table(file, header = FALSE)
  return(AICs(degree_sequence$V1))
}

plots <- function(language, file){
  x=read.table(file, header = FALSE)$V1
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
  p_geom <- function(k) dgeom(k,prob)
  prob_1 <- pars[6]
  p_geom_corrected <- function(k) dgeom(k-1,prob_1)
  
  plot(table(x)/length(x),xlim=c(0,27),ylim=c(0,max((table(x)/length(x))[1],p_geom_corrected(1))),ylab='',xlab='Degree')
  curve(p_zeta,xlim=c(1,27),add=T,col='red',lty=2)
  curve(p_poisson,xlim=c(1,27),add = T,col='blue')
  points(0:27,p_geom(0:27),col='green',lwd=5,lty=10,pch=15)
  points(0:27,p_geom_corrected(0:27),col='brown',lwd=5,lty=10,pch=15)
  curve(p_RT_zeta,xlim=c(1,27),add = T,col='yellow',lty=2,lwd=3)
  legend('topright',cex=0.8,legend = c('Observe degree frequency','Zeta','Right Trunc. Zeta ','Geometric','Geometric from 1','Poisson'),col = c('black','red','yellow','green','brown','blue'),lty = c(1,2,2,NA,NA,1),pch = c(NA,NA,NA,15,15,NA))
  title(language)
  }

for (i in 1:nrow(source)) aic_diff[i,] <- write_AICs(source$language[i], source$file[i])
for (i in 1:nrow(source)) rss[i,] <- write_RSS(source$language[i], source$file[i])


png(file='images//General_1.png')
par(mfrow=c(2,1))
for (i in 1:2) plots(source$language[i], source$file[i])
dev.off()
png(file='images//General_2.png')
par(mfrow=c(2,1))
for (i in 3:4) plots(source$language[i], source$file[i])
dev.off()
png(file='images//General_3.png')
par(mfrow=c(2,1))
for (i in 5:6) plots(source$language[i], source$file[i])
dev.off()

png(file='images//General_4.png')
par(mfrow=c(2,1))
for (i in 7:8) plots(source$language[i], source$file[i])
dev.off()
png(file='images//General_5.png')
par(mfrow=c(2,1),)
?par
for (i in 9:10) plots(source$language[i], source$file[i])
dev.off()

for (i in 1:nrow(source)) plots(source$language[i], source$file[i])

dimnames(aic_diff) <- list(source$language, c("zeta","zeta_2","RT_zeta","geom","poisson",'geom_corrected'))
aic_diff
library(xtable)
xtable(aic_diff)

dimnames(rss) <- list(source$language, c("zeta","RT_zeta"))
rss
library(xtable)
xtable(rss,display=rep("e",3))






# Akaike & Altmann -----------------------------------------------------

#Altmann Distriution
altm_func <- function(gamma,delta,k){
  N=length(k)
  return(((sum((1:N)^(-gamma)*exp(-delta*(1:N))))^(-1))*k^(-gamma)*exp(-delta*k))
}

altm_func=Vectorize(altm_func,vectorize.args = c('delta','gamma'))

AICs_Alt <- function(x,graph_=F){
  
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
  minus_log_likelihood_geom_corrected <- function(q) {
    x_1=x-1
    N <- length(x_1)
    M <- sum(x_1)
    -(M-N)*log(1-q)-N*log(q)
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
  mle_geom_corrected <- mle(minus_log_likelihood_geom_corrected,
                            start = list(q = N/M),
                            method = "L-BFGS-B",
                            lower = c(0.001),
                            upper = c(.999))
  
  #Altmann 
  
  mle_log_altm <- function(par){
    -sum(-N*log(sum((1:N)^(-par[1])*exp(-(1:N)*par[2]))) + sum(log(x^(-par[1]))) + sum(-x*par[2]) ) 
  }
  
  gamma_delta=optim(c(1,1),mle_log_altm,lower = c(0,0), upper = c(1,1),method = 'L-BFGS-B')$par
  
  att_z <- attributes(summary(mle_zeta))
  att_z_G2 <- attributes(summary(mle_zeta_G2))
  att_rt_z <- attributes(summary(mle_RT_zeta))
  att_geom <- attributes(summary(mle_geom))
  att_pois <- attributes(summary(mle_poisson))
  att_geom_corrected <- attributes(summary(mle_geom_corrected))
  
  aics <- c(get_AIC(att_z$m2logL, 1, N),
            get_AIC(att_z_G2$m2logL, 0, N),
            get_AIC(att_rt_z$m2logL, 2, N),
            get_AIC(att_geom$m2logL, 1, N),
            get_AIC(att_pois$m2logL, 1, N),
            get_AIC(att_geom_corrected$m2logL, 1, N),
            get_AIC(-2*mle_log_altm(gamma_delta),2,N)
  )
  coefficients <- c(att_z$coef[1],
                    att_rt_z$coef[1,1],
                    att_rt_z$coef[2,1],
                    att_geom$coef[1],
                    att_pois$coef[1],
                    att_geom_corrected$coef[1],
                    gamma_delta)
  
  if(graph_){
    return (coefficients)
  }
  
  return (aics-min(aics))
}

aic_diff <- matrix(nrow=length(source$language),ncol=7)
dimnames(aic_diff) <- list(source$language, c("zeta","zeta_2","RT_zeta","geom","poisson",'geom_corrected','Altmann'))


# Output ------------------------------------------------------------------


plots <- function(language, file){
  x=read.table(file, header = FALSE)$V1
  pars<-AICs_Alt(x,T)
  gamma_delta <- pars[7:8]
  altm_func_c <- function(k) {altm_func(gamma_delta[1],gamma_delta[2],k)}
  plot(table(x)/length(x),xlim=c(0,27),ylim=c(0,(table(x)/length(x))[1]),ylab='',xlab='Degree')
  curve(altm_func_c,xlim=c(1,27),add=T,col='red',lty=2)
  title(paste(language, 'Altmann'))
  
}


write_AICs <- function(language,file) {
  degree_sequence = read.table(file, header = FALSE)
  return(AICs_Alt(degree_sequence$V1))
}

for (i in 1:nrow(source)) aic_diff[i,] <- write_AICs(source$language[i], source$file[i])

png(file='images//Altman_1.png')
par(mfrow=c(2,1),lheight=15)
for (i in 1:2) plots(source$language[i], source$file[i])
dev.off()
png(file='images//Altman_2.png')
par(mfrow=c(2,1))
for (i in 3:4) plots(source$language[i], source$file[i])
dev.off()
png(file='images//Altman_3.png')
par(mfrow=c(2,1))
for (i in 5:6) plots(source$language[i], source$file[i])
dev.off()

png(file='images//Altman_4.png')
par(mfrow=c(2,1))
for (i in 7:8) plots(source$language[i], source$file[i])
dev.off()
png(file='images//Altman_5.png')
par(mfrow=c(2,1))
for (i in 9:10) plots(source$language[i], source$file[i])
dev.off()

xtable(aic_diff)
