# Fit the multi-level functional Beta model to CGM data
##############################################################
# Author: Irina Gaynanova
# Use both original and perturbed values of min/max; compare the difference

rm(list=ls())

# Supporting libraries
library(tidyverse)
library(refund)

# Specify path to save results
results.path = '~/Documents/Research/Glucometers/NightTrajectories/Data_rda/'

# Load real CGM data (how the code was used in the paper)
load("CGMdataReal.Rdata")
data = CGMdataReal
minval = minvalReal
maxval = maxvalReal
subjectIDs = idsReal

# In absence of real CGM data, the code below will work with provided synthetic data (except for cases where real subject IDs are used)
# load("CGMdataSynthetic.Rdata")
# data = simul_beta_real
# minval = minvalSythetic
# maxval = maxvalSynthetic
# subjectIDs = idsSynthetic

# Estimate marginal parameters for multi-level functional Beta model, use mean/sd for normal data generation; alpha/beta/m/M for beta data generation
source("BetaEstimatingFunctions.R")
# Original model fit
out_mfbeta1 <- estimateBetaParameters(data, minval, maxval)
# Model fit with 5% perturbation on the values of min and max
out_mfbeta2 <- estimateBetaParameters(data, 0.95*minval, 1.05*maxval)


# Some exploration
# Quantile comparison between original and misspecified model
percentiles = seq(0.01, 0.99, by = 0.01)
n = nrow(out_mfbeta1$alphaM)
nt = ncol(out_mfbeta1$alphaM)
q1 <- array(0, dim = c(length(percentiles), nt, n))
q2 <- q1
q3 <- q1
for (i in 1:n){
  for (t in 1:nt){
    q1[, t, i] = minval[i] + (maxval[i] - minval[i]) * qbeta(percentiles, out_mfbeta1$alphaM[i, t], out_mfbeta1$betaM[i, t])
    q2[, t, i] = 0.95*minval[i] + (1.05*maxval[i] - 0.95*minval[i]) * qbeta(percentiles, out_mfbeta2$alphaM[i, t], out_mfbeta2$betaM[i, t])
    q3[, t, i] = qnorm(percentiles, mean = out_mfbeta1$pca_means$Yhat[i, t], sd = out_mfbeta1$pca_sds$Yhat[i, t])
  }
}

# Need to compare resulting quantiles: 
# q1 - Beta with m, M;
# q2 - Beta with perturbed m, M;
# q3 - Normal marginal

# Average difference by subject
q12_diff = rep(0, n)
q13_diff = rep(0, n)
for (i in 1:n){
  q12_diff[i] = mean(abs(q1[, , i] - q2[, , i]))
  q13_diff[i] = mean(abs(q1[, , i] - q3[, , i]))
}
summary(q12_diff)
# median 0.2, mean 0.3, IQR 0.12, 0.35, min 0.17, max 2.8
summary(q13_diff)
# median 2.4, mean 2.8, IQR 1.12 - 3.22, max 15, min 0.24

# Look at particular quantiles
#####################################################
# Average behaviour across all quantiles
q12_diff = rep(0, 99)
q13_diff = rep(0, 99)
q12_diff_max = rep(0, 99)
q13_diff_max = rep(0, 99)
for (q in 1:99){
  q12_diff[q] = mean(abs(q1[q, , ] - q2[q, , ])/q1[q, , ])
  q13_diff[q] = mean(abs(q1[q, , ] - q3[q, , ])/q1[q, , ])
  q12_diff_max[q] = max(abs(q1[q, , ] - q2[q, , ])/q1[q, , ])
  q13_diff_max[q] = max(abs(q1[q, , ] - q3[q, , ])/q1[q, , ])
}

np = length(percentiles)
quantile_diff <- data.frame(percentiles = rep(percentiles, 4), type = c(rep("Average", 2*np), rep("Maximal", 2*np)), method = rep(c(rep("Perturbed Beta", np), rep("Normal", np)), 2), difference = c(q12_diff, q13_diff, q12_diff_max, q13_diff_max))

library(ggplot2)

p = ggplot(quantile_diff, aes(x = percentiles, y = round(difference*100, 2), col = method, shape = method)) + geom_point() + facet_grid(type~., scales = "free_y") + ylab("Percent difference") + xlab("Quantile level") + geom_hline(yintercept = 5, col = "black") + xlim(c(0.01, 0.99))+ theme(text = element_text(size = 16)) 

figure.path = getwd()

pdf(paste(figure.path, "/mM_Stability.pdf", sep = ""), onefile = F, width = 10, height = 4)
print(p)
dev.off()