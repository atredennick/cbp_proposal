####################################################################
##  wintann_pop_model.R: simulates long term population dynamics  ##
##  of species in a Sonoran Desert winter annual plant community. ##
##  The model is a generic annual plant model with seedbank.      ##
####################################################################

# Clear the workspace
rm(list=ls(all.names = TRUE))

# Load some libraries
library(plyr)
library(reshape2)
library(ggthemes)
library(ggplot2)

do_species <- "erla"

####
####  ANNUAL PLANT POPULATION MODEL -- DENSITY-INDEPENDENT
####
# From Gremer and Venable 2014 (Ecology Letters)
project_popmod <- function(seed_surv_old, seed_surv_new, seeds_per_seedling, 
                           germination_frac, nongerm_seed_surv)
{
  pgr <- seeds_per_seedling*germination_frac*seed_surv_new + nongerm_seed_surv*seed_surv_old*(1-germination_frac)
}



####
####  READ IN ANNUAL PLANT DATA
####
species_year_data <- read.csv("../data/species_x_year.csv")
single_spp_data <- subset(species_year_data, species==do_species)



####
####  SIMULATE LONG-TERM PER CAPITA GROWTH RATE
####
single_spp_data <- subset(single_spp_data, germ.fraction!="." & lxbx!=".")
seeds_per_seedling_vec <- as.numeric(as.character(single_spp_data$lxbx))
seed_surv_new <- 0.153 # from Gremer and Veneble 2014
seed_surv_old <- 0.56 # from Gremer and Veneble 2014
nongerm_seed_surv <- 0.823 # from Gremer and Veneble 2014
germination_frac_vec <- as.numeric(as.character(single_spp_data$germ.fraction))
iters <- 1000
pgr <- numeric(iters)
germ_fracs <- seq(0,1,0.01)
geom_mean <- numeric(length(germ_fracs))
arith_mean <- numeric(length(germ_fracs))
count <- 1
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
random_years <- sample(c(1:length(seeds_per_seedling_vec)), iters, replace = TRUE)
for(germ in germ_fracs){
  for(t in 1:iters){
    randyr <- random_years[t]
    pgr[t] <- project_popmod(seed_surv_old = seed_surv_old, 
                             seed_surv_new = seed_surv_new, 
                             seeds_per_seedling = seeds_per_seedling_vec[randyr],
                             germination_frac = germ,
                             nongerm_seed_surv = nongerm_seed_surv)
  }
  geom_mean[count] <- gm_mean(pgr)
  arith_mean[count] <- log(mean(pgr[pgr>0]))
  count <- count+1
}



####
####  PLOT SOME RESULTS
####
par(mfrow=c(1,2))
plot(germ_fracs, geom_mean, type="l", main=do_species)
plot(germ_fracs, arith_mean, type="l", main=do_species)
