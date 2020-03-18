# Assess model validity in estimating quantiles based on simulated data
# ###########################################################
# Author: Irina Gaynanova

# ATTENTION: requires real CGM data for valid comparison

# Load real CGM data (how the code was used in the paper)
load("CGMdataReal.Rdata")
data = CGMdataReal
minval = minvalReal
maxval = maxvalReal
subjectIDs = idsReal

## Load the simulated CGM data - each subject has 14 curves
#############################################################################
load("CGMdataSynthetic.Rdata")

# Supporting libraries
library(tidyverse)
library(refund)

# Specify path to save figures, working directory by default
figure.path = getwd()

# Fit Beta model on real data, those were the true parameters underlying synthetic data generator
source('BetaEstimatingFunctions.R')
out_mfbeta <- estimateBetaParameters(data, minval, maxval)



# If in the selected data actually see something smaller or larger, need to adjust, minval/maxval are misspecified compared to the true minvalnew/maxvalnew
minvalnew2 = minval
maxvalnew2 = maxval
n = length(simul_beta_real) # Total number of subjects
for (i in 1:n){
  # Adjist min/max in case see something smaller
  minvalnew2[i] = min(minval[i], min(simul_beta_real[[i]]))
  maxvalnew2[i] = max(maxval[i], max(simul_beta_real[[i]]))
}

################################################
# Apply on beta-simulated data  
################################################
# Fit the model on simulated data
modelfit <- estimateBetaParameters(simul_beta_real, minvalnew2, maxvalnew2)

# For each subject, construct percentiles and compare with true data generating mechanism on multiple quantiles
quantiles <- c(0.025, 0.05, 0.1, 0.9, 0.95, 0.975)
nq <- length(quantiles)
maxH = 7 # Set the maximal hour
tdexcom = round(maxH * 60 / 5) # Length of time grid in 5 minute intervals from the estimated sleep onset

perc_norm <- array(0, dim = c(nq, n, tdexcom))
perc_beta <- perc_norm
perc_true <- perc_norm

skewness <- matrix(0, n, tdexcom)
for (i in 1:n){
  # True percentiles for synthetic data
  for (q in 1:nq){
    perc_true[q,i,] = minvalSynthetic[i] + (maxvalSynthetic[i]-minvalSynthetic[i])*qbeta(quantiles[q], out_mfbeta$alphaM[i,], out_mfbeta$betaM[i,])
  }
    
  # Percentiles estimated using normal distribution  
  meani = modelfit$means_smooth[i,]
  sdi = modelfit$sd_smooth[i,]
  sdi[sdi <= 0] = min(sdi[sdi > 0])
  for (q in 1:nq){
    perc_norm[q,i,] = qnorm(quantiles[q], mean = meani, sd = sdi)
  }
  
  # Percentiles estimated using beta distribution 
  for (q in 1:nq){
    perc_beta[q,i,] = minvalnew2[i] + (maxvalnew2[i]-minvalnew2[i])*qbeta(quantiles[q], modelfit$alphaM[i,], modelfit$betaM[i,])
  }
  
  # Calculate true skewness
  skewness[i, ] = 2 * (out_mfbeta$betaM[i, ] - out_mfbeta$alphaM[i, ]) * sqrt(out_mfbeta$alphaM[i, ] + out_mfbeta$betaM[i, ] + 1) / ((out_mfbeta$alphaM[i, ] + out_mfbeta$betaM[i, ] + 2) * sqrt(out_mfbeta$alphaM[i, ] * out_mfbeta$betaM[i, ]))
}

# Skewed subjects with positive skew
###########################################
skew_thresh = 0.6
averages = rowMeans(skewness)
# Select subjects with positive skew
selecteds = averages > skew_thresh
sum(selecteds) # 23 subjects

# Pretty picture for 5th and 10th percentiles
###########################################
iq = which(quantiles == 0.05)
# Comparison with normal on skewed subjects
summary(c(perc_norm[iq, selecteds,] - perc_true[iq, selecteds,])) #median = -9.5
IQR(c(perc_norm[iq, selecteds,] - perc_true[iq, selecteds,])) # IQR = 8.5

# Comparison with beta on skewed subjects
summary(c(perc_beta[iq, selecteds,] - perc_true[iq, selecteds,])) #median = 0.6
IQR(c(perc_beta[iq, selecteds,] - perc_true[iq, selecteds,])) #IQR = 3.0

data_first <- data.frame(Beta = c(perc_beta[iq, selecteds,] - perc_true[iq, selecteds,]), Normal = c(perc_norm[iq, selecteds,] - perc_true[iq, selecteds,]), percentile = rep("5%", sum(selecteds)*tdexcom), type = rep(paste("Beta distribution, skewness > ", skew_thresh, sep=""), sum(selecteds)*tdexcom))

# Comparison with normal on other subjects
summary(c(perc_norm[iq, !selecteds,] - perc_true[iq, !selecteds,])) #median = -3.1
IQR(c(perc_norm[iq, !selecteds,] - perc_true[iq, !selecteds,])) #IQR = 4.3

# Comparison with beta on other subjects
summary(c(perc_beta[iq, !selecteds,] - perc_true[iq, !selecteds,])) #median = 0.4
IQR(c(perc_beta[iq, !selecteds,] - perc_true[iq, !selecteds,])) #IQR = 2.4

data_first <- rbind(data_first, data.frame(Beta = c(perc_beta[iq, !selecteds,] - perc_true[iq, !selecteds,]), Normal = c(perc_norm[iq, !selecteds,] - perc_true[iq, !selecteds,]), percentile = rep("5%", (n-sum(selecteds))*tdexcom), type = rep(paste("Beta distribution, skewness <= ", skew_thresh, sep=""), (n-sum(selecteds))*tdexcom)))

# 0.1 quantile
################################
iq = which(quantiles == 0.1)

# Comparison with normal on skewed subjects
summary(c(perc_norm[iq, selecteds,] - perc_true[iq, selecteds,])) #median = -3.0
IQR(c(perc_norm[iq, selecteds,] - perc_true[iq, selecteds,])) #IQR = 3.9

# Comparison with beta on skewed subjects
summary(c(perc_beta[iq, selecteds,] - perc_true[iq, selecteds,])) #median = 0.1
IQR(c(perc_beta[iq, selecteds,] - perc_true[iq, selecteds,])) #IQR = 3.1

# Comparison with normal on other subjects
summary(c(perc_norm[iq, !selecteds,] - perc_true[iq, !selecteds,])) #median = -0.5
IQR(c(perc_norm[iq, !selecteds,] - perc_true[iq, !selecteds,])) #IQR = 2.0

# Comparison with beta on other subjects
summary(c(perc_beta[iq, !selecteds,] - perc_true[iq, !selecteds,])) # median = -0.3
IQR(c(perc_beta[iq, !selecteds,] - perc_true[iq, !selecteds,])) #IQR = 2.0

data_first <- rbind(data_first, data.frame(Beta = c(perc_beta[iq, selecteds,] - perc_true[iq, selecteds,]), Normal = c(perc_norm[iq, selecteds,] - perc_true[iq, selecteds,]), percentile = rep("10%", sum(selecteds)*tdexcom), type = rep(paste("Beta distribution, skewness > ", skew_thresh, sep=""), sum(selecteds)*tdexcom)))

data_first <- rbind(data_first, data.frame(Beta = c(perc_beta[iq, !selecteds,] - perc_true[iq, !selecteds,]), Normal = c(perc_norm[iq, !selecteds,] - perc_true[iq, !selecteds,]), percentile = rep("10%", (n-sum(selecteds))*tdexcom), type = rep(paste("Beta distribution, skewness <= ", skew_thresh, sep=""), (n-sum(selecteds))*tdexcom)))



library(reshape)
data <- melt(data_first, variable_name = "Model")

library(ggplot2)

p <- ggplot(data, aes(x = Model, y = value, fill = Model)) + geom_boxplot() + facet_grid(percentile ~ type) + ylab("Difference in glucose values [mg / dL]") + guides(fill = FALSE)
print(p)

################################################
# Apply on normal-simulated data  
################################################
# If in the selected data actually see something smaller or larger, need to adjust, minval/maxval are misspecified compared to the true minvalnew/maxvalnew
minvalnew2 = minval
maxvalnew2 = maxval

for (i in 1:n){
  # Adjist min/max in case see something smaller
  minvalnew2[i] = min(minval[i], min(simul_normal_real[[i]]))
  maxvalnew2[i] = max(maxval[i], max(simul_normal_real[[i]]))
}

# Fit the model on normal simulated data 
modelfit <- estimateBetaParameters(simul_normal_real, minvalnew2, maxvalnew2)

# For each subject, construct percentiles and compare with true data generating mechanism
perc_norm <- array(0, dim = c(nq, n, tdexcom))
perc_beta <- perc_norm
perc_true <- perc_norm

for (i in 1:n){
  # True percentiles for simulated normal data
  for (q in 1:nq){
    perc_true[q,i,] = qnorm(quantiles[q], out_mfbeta$pca_means$Yhat[i,], out_mfbeta$pca_sds$Yhat[i,])
  }
  
  # Percentiles estimated using normal distribution  
  meani = modelfit$means_smooth[i,]
  sdi = modelfit$sd_smooth[i,]
  sdi[sdi <= 0] = min(sdi[sdi > 0])
  for (q in 1:nq){
    perc_norm[q,i,] = qnorm(quantiles[q], mean = meani, sd = sdi)
  }
  
  # Percentiles estimated using beta distribution  
  for (q in 1:nq){
    perc_beta[q,i,] = minvalnew2[i] + (maxvalnew2[i]-minvalnew2[i])*qbeta(quantiles[q], modelfit$alphaM[i,], modelfit$betaM[i,])
  }
  
}


# Create data frame for underlying normal distribution
###########################################
iq = which(quantiles == 0.05)

# Comparison with normal
summary(c(perc_norm[iq, ,] - perc_true[iq, ,])) #median = 0
IQR(c(perc_norm[iq, ,] - perc_true[iq, ,])) # almost 0

# Comparison with beta
summary(c(perc_beta[iq, ,] - perc_true[iq, ,])) #median = 2.1
IQR(c(perc_beta[iq, ,] - perc_true[iq, ,])) # IQR = 3

data_second <- data.frame(Beta = c(perc_beta[iq, ,] - perc_true[iq, ,]), Normal = c(perc_norm[iq, ,] - perc_true[iq, ,]), percentile = rep("5%", n*tdexcom), type = rep("Normal distribution", n*tdexcom))

iq = which(quantiles == 0.1)
data_second <- rbind(data_second, data.frame(Beta = c(perc_beta[iq, ,] - perc_true[iq, ,]), Normal = c(perc_norm[iq, ,] - perc_true[iq, ,]), percentile = rep("10%", n*tdexcom), type = rep("Normal distribution", n*tdexcom)))

library(reshape)
data2 <- melt(data_second, variable_name = "Model")

library(ggplot2)

p <- ggplot(data2, aes(x = Model, y = value, fill = Model)) + geom_boxplot() + facet_grid(percentile ~ type) + ylab("Difference in glucose values [mg / dL]") + guides(fill = FALSE)
print(p)

# Combine everything together for the comparative picture of 5% and 10% quantiles
data_all <- rbind(data, data2)

pdf(paste(figure.path, "/Quantiles_5_10.pdf", sep = ""), onefile = TRUE, width = 14, height = 8)
p <- ggplot(data_all, aes(x = Model, y = value, fill = Model)) + geom_boxplot() + facet_grid(percentile ~ type) + ylab("Difference in glucose values [mg / dL]") + guides(fill = FALSE) + scale_y_continuous(breaks=c(-40,-20,0,20))+ theme(text = element_text(size = 18))
print(p)
dev.off()

