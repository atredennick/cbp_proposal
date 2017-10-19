##  Conceptual plots for proposal
library(ggplot2)

boots <- 10
samps <- 10000
dens_df <- list()
sdvec <- seq(1,2,length.out = boots)

for(i in 1:boots){
  x <- rnorm(samps, mean = 1, sd = sdvec[i])
  xdf <- data.frame(boot=i, dens=x)
  dens_df <- rbind(dens_df, xdf)
}

ggplot(dens_df, aes(x=dens, color=as.factor(boot)))+
  geom_line(stat = "density", adjust=3)+
  theme_bw()+
  guides(color=FALSE)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = boots
cols = gg_color_hue(n)
ggplot(subset(dens_df, boot==boots), aes(x=dens))+
  geom_line(stat = "density", adjust=3, color=cols[boots])+
  theme_bw()+
  scale_y_continuous(limits=c(0,0.38))+
  guides(color=FALSE)





