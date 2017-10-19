##  R script to fit an annual plant model to the winter annual plant
##  dataset from Larry Veneble, Tucsan, AZ

rm(list=ls(all.names = TRUE))

VERBOSE <- FALSE


####
####  Load Libraries -----------------------------------------------------------
####
library(ggplot2)
library(ggthemes)
library(plyr)
library(reshape2)
library(MASS)



####
####  Read In Winter Annual Datasets -------------------------------------------
####
plant_dat         <- read.csv("../data/species_x_year.csv", na.strings = ".")
plant_dat$species <- tolower(plant_dat$species)
dominant_species  <- c("erla", "erci", "erte", "evmu", "mobe", "pehe", "pere", 
                       "plpa", "plin", "scba", "stmi", "vuoc")
dominant_data     <- subset(plant_dat, species %in% dominant_species)
fitting_dat       <- subset(dominant_data, species=="erla" & !is.na(germ.fraction))



####
####  Fit a distribution
####
logit<-function(x) log(x/(1-x))
antilogit <- function(x) { exp(x) / (1 + exp(x) ) }
x <- logit(fitting_dat$germ.fraction)
x[which(x==Inf)] <- 36
x[which(x==-Inf)] <- -36
norm_params <- fitdistr(x, "normal")

germ_sim <- rnorm(100000, norm_params$estimate[1], norm_params$estimate[2])

png("../figures/germ_distr_example.png", width = 8, height = 3, units = "in", res = 120)
par(mfrow=c(1,2), mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1), las=1)
plot(fitting_dat$year, fitting_dat$germ.fraction, type="l", frame.plot = FALSE, 
     xlim=c(1990,2015), xlab="Year", ylab="Germination Fraction")
points(fitting_dat$year, fitting_dat$germ.fraction, pch=19)

hist(antilogit(germ_sim), breaks = 20, col="grey", border = "white", 
     freq = F, main="", xlab="Germination Rate")

sds <- seq(0.5,5,by = 1)
count=1
for(do_sd in sds){
  germ_sim <- rnorm(100000, norm_params$estimate[1], do_sd)
  lines(density(antilogit(germ_sim), adjust=2), col=count, lwd=2)
  # hist(antilogit(germ_sim), add=TRUE, col=NA, border=count, freq = F)
  count <- count+1
}
dev.off()