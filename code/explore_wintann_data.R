#################################################################
##  explore_wintann_data.R: R script to explore Sonoran Desert ##
##  winter annual demographic dataset                          ##
#################################################################

# Clear the workspace
rm(list=ls(all.names = TRUE))

# Load some libraries
library(plyr)
library(reshape2)
library(ggthemes)
library(ggplot2)



####
####  LOAD SEED BANK DATA
####
seed_bank <- read.csv("../data/Seed_Bank.csv")
seed_bank <- subset(seed_bank, viable.seeds!=-99)
seed_bank$species <- tolower(seed_bank$species)

# Summarise data by species and year
avg_seed_bank <- ddply(seed_bank, .(year, species), summarise,
                       viable_seeds = mean(viable.seeds))
ggplot(avg_seed_bank, aes(x=year, y=viable_seeds, color=species))+
  geom_line()


# Plot time series for abundant species from Gremer and Venable
gremer_species <- c("erla", "erci", "erte", "evmu", "mobe", "pehe", "pere", 
                    "plpa", "plin", "scba", "stmi", "vuoc")
abundant_seeds <- avg_seed_bank[which(avg_seed_bank$species %in% gremer_species),]
ggplot(abundant_seeds, aes(x=year, y=viable_seeds, color=species))+
  geom_line(alpha=0.5)+
  geom_point()



####
####  LOAD INDIVIDUAL CENSUS DATA
####
census_data <- read.csv("../data/census_data.csv")
census_data$species <- tolower(census_data$species)

# Summarise data to get average abundance per year
abund_per_plot <- ddply(census_data, .(Year, plot.habitat.replicate, species), summarise,
                        abundance = length(which(seeds!=-99)))
avg_abundance <- ddply(abund_per_plot, .(Year, species), summarise,
                       mean_abundance = mean(abundance))

# Plot abundant species only
abundant_abundance <- avg_abundance[which(avg_abundance$species %in% gremer_species),]
total_abundance <- ddply(abundant_abundance, .(Year), summarise,
                         abundance = sum(mean_abundance))
cv_abundance <- sd(total_abundance$abundance) / mean(total_abundance$abundance)
ggplot(abundant_abundance, aes(x=Year, y=mean_abundance, color=species))+
  geom_line(alpha=0.5)+
  geom_point()+
  geom_line(data=total_abundance, aes(x=Year, y=abundance), color="black", size=1, alpha=0.5)+
  geom_point(data=total_abundance, aes(x=Year, y=abundance), color="black", size=3)+
  scale_x_continuous(breaks=unique(abundant_abundance$Year))+
  ylab("Average Seedling Abundance")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(color=FALSE)


