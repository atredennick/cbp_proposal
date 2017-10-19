##  Two Species Annual Plant Population Model: simulates the dynamics of
##  the seedbanks of two competing annual plants. The plants coxist by the
##  storage effect. This is just to show two time series (one with coexistence,
##  one with exclusion) and to show that reducing environmental variability
##  causes exclusion and a reduction in ecosystem stability, here measured by
##  total community abundance.
##
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##  Date created: November 11, 2016


rm(list=ls(all.names = TRUE))



####
#### Parameters ----------------------------------------------------------------
####

sigE   <- c(0.05,0.5)
rho    <- c(0.5,-0.8,-0.5)
s      <- c(0.5, 0.5)
alpha  <- c(1,1)
lambda <- c(101,99)
nTime  <- 200



####
####  Libraries ----------------------------------------------------------------
####
library(mvtnorm)
library(plyr)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(zoo)



####
####  Annual Plant Model Functions ---------------------------------------------
####

### Update seedbank function
updateN <- function(g, s, alpha, lambda, lastN){
  newN <- numeric(2)
  for(i in 1:2){
    newN[i] <- lastN[i]*s[i]*(1-g[i]) + ((lambda[i]*g[i]*lastN[i]) / (1 + (alpha[i]*g[i]*lastN[i] + alpha[-i]*g[-i]*lastN[-i])))
    if(newN[i] < 1) { newN[i] <- 0 }
  }
  return(newN)
}

### Get germination fractions function
getG <- function(sigE, rho, nTime){
  varcov <- matrix(c(sigE, rho*sigE, rho*sigE, sigE), 2, 2)
  e      <- rmvnorm(n = nTime, mean = c(0,0), sigma = varcov)
  g      <- exp(e) / (1+exp(e))
  return(g)
}



####
#### Simulations ---------------------------------------------------------------
####

png("../figures/sim_example.png", width = 5, height=3, units="in", res = 120)

par(mgp=c(2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1), las=1)
layout(matrix(c(1,2,3,4,4,4), 3, 2, byrow = FALSE))

### Low variability // competitive exclusion
gSeries <- getG(sigE[1], rho[2], nTime)
N       <- matrix(data = NA, nrow = nTime, ncol = 2)
N[1,]   <- c(100,20)
for(t in 2:nTime){
  N[t,] <- updateN(g=gSeries[t,], s, alpha, lambda, lastN=N[t-1,])
}
matplot(c(1:nTime), N, type="l", lwd=1, bty="n", xlab="Time", col="black", 
        ylab="N", main="A) Competitive Exclusion", ylim=c(0,150))
totpop1 <- rowSums(N)

### High variability // coexistence
gSeries <- getG(sigE[2], rho[2], nTime)
N       <- matrix(data = NA, nrow = nTime, ncol = 2)
N[1,]   <- c(100,0)
for(t in 2:nTime){
  N[t,] <- updateN(g=gSeries[t,], s, alpha, lambda, lastN=N[t-1,])
}
matplot(c(1:nTime), N, type="l", lwd=1, bty="n", xlab="Time", 
        col="blue",ylab="N", main="B) Coexistence Possible")
totpop2 <- rowSums(N)


### High variability // coexistence
gSeries <- getG(sigE[2], rho[2], nTime)
N       <- matrix(data = NA, nrow = nTime, ncol = 2)
N[1,]   <- c(100,100)
for(t in 2:nTime){
  N[t,] <- updateN(g=gSeries[t,], s, alpha, lambda, lastN=N[t-1,])
}
matplot(c(1:nTime), N, type="l", lwd=1, bty="n", xlab="Time", 
        col="red", ylab="N", main="C) Species Coexistence", ylim=c(0,150))
totpop3 <- rowSums(N)


cv1 <- sd(totpop1)/mean(totpop1)
cv2 <- sd(totpop2)/mean(totpop2)
cv3 <- sd(totpop3)/mean(totpop3)
barplot(c(cv1, cv2, cv3), ylab="CV of Community Abundance", xlab="Environmental Variance",
        ylim=c(0,0.1), col=c("grey45","blue","red"))
text(x=c(0.7, 1.9, 3.1), y=c(cv1, cv2, cv3)+0.005, labels = LETTERS[1:3])
axis(1, at=c(0.7, 1.9, 3.1), labels = c(expression(paste(sigma[E]^2, " = ", 0.05)), 
                                        expression(paste(sigma[E]^2, " = ", 0.5)),
                                        expression(paste(sigma[E]^2, " = ", 0.5))))

dev.off()

### Plot the rolling CV
# mycv <- function(x) {sd(x) / mean(x)}
# rolling_cv <- data.frame(cv1 = rollapply(totpop1, width=10, FUN=mycv, fill=NA),
#                          cv2 = rollapply(totpop2, width=10, FUN=mycv, fill=NA),
#                          iteration = 1:length(totpop1))
# matplot(rolling_cv[,c(1:2)], type="l", col=c("grey45","steelblue"), bty="n")



