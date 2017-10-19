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
library(rjags)
library(coda)
# library(devtools)
# install_github("atredennick/ecoforecastR") # get latest version
library(ecoforecastR)



####
####  Read In Winter Annual Datasets -------------------------------------------
####
plant_dat         <- read.csv("../data/species_x_year.csv", na.strings = ".")
plant_dat$species <- tolower(plant_dat$species)
dominant_species  <- c("erla", "erci", "erte", "evmu", "mobe", "pehe", "pere", 
                       "plpa", "plin", "scba", "stmi", "vuoc")
dominant_data     <- subset(plant_dat, species %in% dominant_species)
fitting_dat       <- subset(dominant_data, species=="erla")
fitting_dat       <- subset(fitting_dat, is.na(seeds_per_m2)==FALSE)


####
####  JAGS State-Space Model ---------------------------------------------------
####
ann_pop_model <- "
model{

  #### Variance Priors
  tau_proc ~ dgamma(0.0001, 0.0001)
  sigma_proc <- 1/sqrt(tau_proc)
  tau_germ ~ dgamma(0.0001, 0.0001)
  sigma_germ <- 1/sqrt(tau_germ)
  sigma.o ~ dunif(1,1000)
  
  #### Fixed Effects Priors
  surv ~ dunif(0,1)
  # lambda ~ dunif(0,100)
  alpha ~ dunif(0,10)
  for(t in 1:yrs){
    germ[t] ~ dunif(0,1)
    lambda[t] ~ dunif(0,100)
  }
  # germ ~ dunif(0,1)
  
  #### Initial Conditions
  N0 ~ dunif(1,10)
  Nmed[1] <- log(max(1, (surv*(1-germ[1])*N0 + ((lambda[1]*germ[1]*N0) / (1+alpha*germ[1]*N0)) )))
  N[1] ~ dlnorm(Nmed[1], tau_proc)

  ####  Process Model
  for(t in 2:yrs){
    Nmed[t] <- log(max(1, (surv*(1-germ[t])*N[t-1] + ((lambda[t]*germ[t]*N[t-1]) / (1+alpha*germ[t]*N[t-1]))  )))
    N[t] ~ dlnorm(Nmed[t], tau_proc)
  }
  
  ####  Data Model
  var.o <- sigma.o*sigma.o
  for(t in 1:yrs){
    shape[t] <- N[t]*N[t]/var.o
    rate[t] <- N[t]/var.o
    lambda2[t] ~ dgamma(shape[t], rate[t])
    Nobs[t] ~ dpois(lambda2[t])
  }

}"



####
####  Fit Annual Plant State-Space Model ---------------------------------------
####
##  Prepare data list
mydat         <- list(Nobs = round(fitting_dat$seeds_per_m2)+1, 
                      yrs = nrow(fitting_dat))
out_variables <- c("germ","surv","lambda","sigma_proc","N","alpha")

##  Send to JAGS
mc3     <- jags.model(file=textConnection(ann_pop_model), data=mydat, n.chains=3)
mc3_out <- coda.samples(model=mc3, variable.names=out_variables, n.iter=10000)

## Split output
out          <- list(params=NULL, predict=NULL, model=ann_pop_model, data=mydat)
mfit         <- as.matrix(mc3_out,chains=TRUE)
pred.cols    <- union(grep("N[",colnames(mfit),fixed=TRUE),grep("Nmed[",colnames(mfit),fixed=TRUE))
chain.col    <- which(colnames(mfit)=="CHAIN")
out$predict  <- mat2mcmc.list(mfit[,c(chain.col,pred.cols)])
out$params   <- mat2mcmc.list(mfit[,-pred.cols])
fitted_model <- out

## Collate predictions
predictions        <- rbind(fitted_model$predict[[1]],
                            fitted_model$predict[[2]],
                            fitted_model$predict[[3]])
median_predictions <- apply(predictions, MARGIN = 2, FUN = "median")
upper_predictions  <- apply(predictions, MARGIN = 2, FUN = function(x){quantile(x, probs = 0.975)})
lower_predictions  <- apply(predictions, MARGIN = 2, FUN = function(x){quantile(x, probs = 0.025)})
prediction_df      <- data.frame(year = fitting_dat$year,
                                 observation = fitting_dat$seeds_per_m2,
                                 median_prediction = median_predictions,
                                 upper_prediction = upper_predictions,
                                 lower_prediction = lower_predictions)
if(VERBOSE){ 
  plot(fitted_model$params)
  gelman.diag(mc3_out)
  heidel.diag(mc3_out)
  }



####
####  Plot The Calibration Data And Predictions --------------------------------
####
pred_color <- "dodgerblue"
ggplot(prediction_df, aes(x=year))+
  geom_ribbon(aes(ymax=upper_prediction, ymin=lower_prediction), fill=pred_color, color=NA, alpha=0.2)+
  geom_line(aes(y=median_prediction), color=pred_color)+
  geom_point(aes(y=observation), size=0.5)+
  ylab("Number of Seeds")+
  xlab("Year")+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(linetype="dotted", color="grey65"))
