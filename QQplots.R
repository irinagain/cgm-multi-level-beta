# Investigate the beta model fit via QQ plots
#############################################################################################
# Author: Irina Gaynanova

# ATTENTION: results in manucscript are based on real CGM data, however the code will work with provided synthetic data

# Load real CGM data (how the code was used in the paper)
load("CGMdataReal.Rdata")
data = CGMdataReal
minval = minvalReal
maxval = maxvalReal
subjectIDs = idsReal

# In absence of real CGM data, the code below will work with provided synthetic data 
# load("CGMdataSynthetic.Rdata")
# data = simul_beta_real
# minval = minvalSythetic
# maxval = maxvalSynthetic
# subjectIDs = idsSynthetic


# Path to save figures to
figure.path = getwd()


library(tidyverse)
library(refund)
library(reshape2)

# Fit multi-level Beta model
source("BetaEstimatingFunctions.R")
out_mfbeta <- estimateBetaParameters(data, minval, maxval)

# Calculate quantiles for QQ plot
nS = length(subjectIDs)
nT = ncol(data[[1]])
plotdataAll <- c()
for (i in 1:nS){
  for (t in 1:nT){
    yit = data[[i]][ , t]
    yit = yit[!is.na(yit)] # remove NA
    nights = length(yit) # number of nights available for subject i at time t
    if (nights > 1){
      percs = ((1:nights) - 0.5)/nights # get percentiles
      quantiles <- minval[i] + (maxval[i] - minval[i]) * qbeta(percs, shape1 = out_mfbeta$alphaM[i, t], shape2 = out_mfbeta$betaM[i, t]) # get theoretical beta quantiles
      quantNormal <- qnorm(percs, mean = out_mfbeta$means_smooth[i, t], sd = out_mfbeta$sd_smooth[i, t]) #get theoretical normal quantiles
      plotdata <- data.frame(Subject = subjectIDs[i], Time = round((t * 5)/60, 2), Sample = sort(yit), Theoretical = quantiles, TheoreticalNorm = quantNormal) # create joint data frame
      plotdataAll <- rbind(plotdataAll, plotdata)
    }
  }
}

# Set some of the plot attributes
library(ggplot2);
theme_update(plot.title    = element_text(size = 22, hjust = 0.5), plot.subtitle = element_text(size = 20, hjust = 0.5), axis.title.x  = element_text(size = 20, hjust = 0.5), axis.title.y  = element_text(size = 20, vjust = 0.5), plot.margin   = unit(c(1, 1, 1, 1), "cm"))


size_text = 16 # text size for plots

# Global quantile plot, everything vs everything
###########################################################
global_qqplot <- ggplot(data = plotdataAll, aes(x = Theoretical, y = Sample)) + geom_point(size = 2, colour = 'blue', alpha = 0.5, shape = 16) + geom_abline(intercept = 0, slope = 1, linetype = 'dashed') + 
  ggtitle('Beta model, global') + 
  xlab('Theoretical Quantiles') + 
  ylab('Sample Quantiles') + ylim(c(0,425)) + xlim(c(0,425)) + geom_rect(aes(xmin = 0, xmax = 100, ymin = 0, ymax = 100), color = "black", fill = "transparent")


png(file = paste(figure.path, "/Global_beta_qqplot.png", sep=""), width = 800, height = 800, res = 150)
print(global_qqplot)
dev.off()

# Same for normal distribution

global_qqplotN <- ggplot(data = plotdataAll, aes(x = TheoreticalNorm, y = Sample)) +
  geom_point(size = 2, colour = 'red', alpha = 0.5, shape = 17) + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') + 
  ggtitle('Normal model, global') + 
  xlab('Theoretical Quantiles') + 
  ylab('Sample Quantiles') + ylim(c(0,425)) + xlim(c(0,425)) + geom_rect(aes(xmin = 0, xmax = 100, ymin = 0, ymax = 100), color = "black", fill = "transparent")


png(file = paste(figure.path, "/Global_normal_qqplot.png", sep=""), width = 800, height = 800, res = 150)
print(global_qqplotN)
dev.off()

# Subject-specific quantile plots for selected subjects
###########################################################
library(tidyverse)

# Create plots for subjects 70142 and 70186 (real ids)
for (i in c(18,28)){
  # Make sure the quantiles are symmetric, but different range for each subject
  if (i == 18){
    ranges = c(125, 400)
    xmin = 125
    ymin = 125
    xmax = 200
    ymax = 200
  }else{
    ranges = c(25, 275)
    xmin = 25
    ymin = 25
    xmax = 100
    ymax = 100
  }
  
  png(file = paste(figure.path, "/Subject", subjectIDs[i],
                   "_beta_normal_qqplot.png", sep = ""), width = 800, height = 800, res = 150)
  subject_qqplot <- plotdataAll%>%
    filter(Subject == subjectIDs[i])%>%
    ggplot(aes(x = Theoretical, y = Sample)) +
    geom_point(size = 2, colour = 'blue', alpha = 0.5, shape = 16) + 
    geom_point(aes(x = TheoreticalNorm, y = Sample), size = 2, col = "red", alpha = 0.5, shape = 17) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') + 
    ggtitle(paste('Subject ', subjectIDs[i], ", all times", sep = "")) + 
    xlab('Theoretical Quantiles') + 
    ylab('Sample Quantiles') + xlim(ranges) + ylim(ranges)+ geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), color = "black", fill = "transparent")
  print(subject_qqplot)
  dev.off()
}


# Time-specific quantile plots
###########################################################

# Create plots for t = 1 and 2 hours from sleep onset
for (t in c(12,24)){
  png(file = paste(figure.path, "/Time", round((t * 5)/60),
                  "_beta_normal_qqplot.png", sep = ""), width = 800, height = 800, res = 150)
  time_qqplot <- plotdataAll%>%
    filter(Time == round((t * 5)/60, 2))%>%
    ggplot(aes(x = Theoretical, y = Sample)) +
    geom_point(size = 2, colour = 'blue', alpha = 0.5, shape = 16) + 
    geom_point(aes(x = TheoreticalNorm, y = Sample), size = 2, col = "red", alpha = 0.5, shape = 17) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') + 
    ggtitle(paste('Time t = ', round((t * 5)/60), "hour(s), all subjects", sep = " ")) + 
    xlab('Theoretical Quantiles') + 
    ylab('Sample Quantiles')+ ylim(c(0,400)) + xlim(c(0,400))+ geom_rect(aes(xmin = 0, xmax = 100, ymin = 0, ymax = 100), color = "black", fill = "transparent")
  print(time_qqplot)
  dev.off()
}

