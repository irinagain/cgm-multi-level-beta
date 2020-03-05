# Investigate the performance of naive methods for estimation of marginal quantile bands, specifically try parametric Beta and nonparametric density estimation approach
##############################################################
# Author: Irina Gaynanova


rm(list=ls())

# Supporting libraries
library(tidyverse)
library(refund)


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

n = length(subjectIDs)
ngrid = ncol(data[[1]]) #length of time grid
time_grid_dexcom = c(0, 1:(ngrid-1) * 5 /60)

# Estimate marginal parameters for multi-level functional Beta model, use mean/sd for normal data generation; alpha/beta/m/M for beta data generation
source("BetaEstimatingFunctions.R")
# Original model fit
out_mfbeta <- estimateBetaParameters(data, minval, maxval)

# NOT RUN: below only works with subjects ids from real data

# ID numbers for 6 selected subjects
selectedID <- c(70058, 70134, 70139, 70142, 70153, 70186)

##########################################################################################
# Figure with 95% confidence bands for each of 6 subjects and corresponding score values
##########################################################################################

# Extract smoothed alphas
selectedAlphas <- as.data.frame(t(out_mfbeta$alphaM[subjectIDs %in% selectedID,]))
selectedAlphas$t <- time_grid_dexcom
selectedAlphas <- gather(selectedAlphas, key = idnew, value = alpha, -7)

# Extract smoothed betas
selectedBetas <- as.data.frame(t(out_mfbeta$betaM[subjectIDs %in% selectedID,]))
selectedBetas$t <- time_grid_dexcom
selectedBetas <- gather(selectedBetas, key = idnew, value = beta, -7)

# Merge
selectedAlphas$beta <- selectedBetas$beta


datagl = c()
nid = which(subjectIDs %in% selectedID)
for (i in nid){
  subject_data <- reshape2::melt(data[[i]], value.name = "gl")
  subject_data[,2] <- subject_data[, 2] * 5/60
  colnames(subject_data) <- c("sleepday", "sleeptime", "gl")
  subject_data$id <- subjectIDs[i]
  datagl <- rbind(datagl, subject_data)
}

datagl$idnew <- as.factor(as.character(datagl$id))

selectedAlphas$idnew <- as.factor(selectedAlphas$idnew)
levels(selectedAlphas$idnew) <- levels(datagl$idnew)
# Get minimal and maximal values
order <- match(selectedAlphas$idnew, subjectIDs)
selectedAlphas$minval = minval[order]
selectedAlphas$maxval = maxval[order]

figure.path = getwd()

text_size = 24

p = ggplot(datagl,aes(x = as.numeric(sleeptime), y = gl)) + geom_line(aes(group=sleepday)) + geom_line(aes(x = t, y = minval + (maxval-minval)*qbeta(1-0.05/2, alpha, beta), col = idnew), data = selectedAlphas, size = 1.2)+ geom_line(aes(x = t, y = minval + (maxval-minval)*qbeta(0.05/2, alpha, beta), col = idnew), data = selectedAlphas, size = 1.2)+ geom_line(aes(x = t, y = minval + (maxval-minval)*qbeta(0.5, alpha, beta), col = idnew), data = selectedAlphas, size = 1.2)+ facet_wrap(~as.character(idnew), ncol = 2)+ xlim(c(0,7))+ xlab("Hours of sleep") + ylab("Glucose  [mg / dL]") + ylim(c(50, 400))+ theme(legend.position="none")+ theme(text = element_text(size = text_size))
print(p)

# Altnernative 1 - same min/max, but separate Beta on each time point
####################################################################

medians = matrix(NA, n, ngrid)
lowerQ = matrix(NA, n, ngrid)
upperQ = matrix(NA, n, ngrid)

for (i in 1:n){
  for (t in 1:ngrid){
    x = as.numeric(data[[i]][,t])
    x = x[!is.na(x)]
    if (length(x) > 1){
      out = getAlphaBeta(mean(x), sd(x), lower = minval[i], upper = maxval[i])
      medians[i, t] = minval[i] + (maxval[i]-minval[i])*qbeta(0.5, out$alpha, out$beta)
      lowerQ[i, t] = minval[i] + (maxval[i]-minval[i])*qbeta(0.05/2, out$alpha, out$beta)
      upperQ[i, t] = minval[i] + (maxval[i]-minval[i])*qbeta(1-0.05/2, out$alpha, out$beta)
    }
  }
}

BetaPointFit <- data.frame(median = as.vector(medians), upperQ = as.vector(upperQ), lowerQ = as.vector(lowerQ), id = rep(subjectIDs, ngrid), t = rep(time_grid_dexcom + 0.08333333, each = n))

BetaPointFitSelect <- BetaPointFit[BetaPointFit$id %in% selectedID, ]
BetaPointFitSelect$idnew <-  as.factor(as.character(BetaPointFitSelect$id))

pdf(file = paste(figure.path, "/RepresentativeSBeta_Col_pointwise.pdf", sep = ""), onefile = F, width = 16, height = 12)
p = ggplot(datagl,aes(x = as.numeric(sleeptime), y = gl)) + geom_line(aes(group=sleepday)) + geom_line(aes(x = t, y = upperQ, col = idnew), data = BetaPointFitSelect, size = 1.2)+ geom_line(aes(x = t, y = lowerQ, col = idnew), data = BetaPointFitSelect, size = 1.2)+ geom_line(aes(x = t, y = median, col = idnew), data = BetaPointFitSelect, size = 1.2)+ facet_wrap(~as.character(idnew), ncol = 2)+ xlim(c(0,7))+ xlab("Hours of sleep") + ylab("Glucose  [mg / dL]") + ylim(c(40, 410))+ theme(legend.position="none")+ theme(text = element_text(size = text_size))
print(p)
dev.off()

# Altnernative 2 - nonparametric approach (kernel density estimator)
####################################################################
library(KernSmooth)

# Function to return pth quantie from bkde density estimate (stored in out)
quantile_bkde <- function(out, p){
  # Get the grid of points
  xx <- out$x
  yy <- out$y
  nn <- length(xx)
  # Approximate cdf using trapezoid rule for integrals
  Fx <- cumsum(yy * c(0, diff(xx)))
  # Renormalize to 1 just in case
  Fx <- Fx/Fx[nn]

  # Do inverse cdf transformation for quantiles
  ii <- min(which(Fx >= p))
  if (!is.na(ii) && ii >= 1 && ii <= nn){ 
    qs <- xx[ii]
  }else{
    qs <- NA
  }
  return(qs)
}

medians = matrix(NA, n, ngrid)
lowerQ = matrix(NA, n, ngrid)
upperQ = matrix(NA, n, ngrid)

nid = which(subjectIDs %in% selectedID)
for (i in 1:n){
  for (t in 1:ngrid){
    x = as.numeric(data[[i]][,t])
    x = x[!is.na(x)]
    if (length(x) > 1){
      # Kernel density estimation
      out = bkde(x, range.x = c(40, 500)) 
      # Find quantiles
      medians[i, t] = quantile_bkde(out, p = 0.5)
      lowerQ[i, t] = quantile_bkde(out, p = 0.05/2)
      upperQ[i, t] = quantile_bkde(out, p = 1 - 0.05/2)
    }
  }
}

KernelPointFit <- data.frame(median = as.vector(medians), upperQ = as.vector(upperQ), lowerQ = as.vector(lowerQ), id = rep(subjectIDs, ngrid), t = rep(time_grid_dexcom + 0.08333333, each = n))

KernelPointFitSelect <- KernelPointFit[KernelPointFit$id %in% selectedID, ]
KernelPointFitSelect$idnew <-  as.factor(as.character(KernelPointFitSelect$id))

pdf(file = paste(figure.path, "/RepresentativeSKernel_Col_pointwise.pdf", sep = ""), onefile = F, width = 16, height = 12)
p = ggplot(datagl,aes(x = as.numeric(sleeptime), y = gl)) + geom_line(aes(group=sleepday)) + geom_line(aes(x = t, y = upperQ, col = idnew), data = KernelPointFitSelect, size = 1.2)+ geom_line(aes(x = t, y = lowerQ, col = idnew), data = KernelPointFitSelect, size = 1.2)+ geom_line(aes(x = t, y = median, col = idnew), data = KernelPointFitSelect, size = 1.2)+ facet_wrap(~as.character(idnew), ncol = 2)+ xlim(c(0,7))+ xlab("Hours of sleep") + ylab("Glucose  [mg / dL]") + ylim(c(40, 410))+ theme(legend.position="none")+ theme(text = element_text(size = text_size))
print(p)
dev.off()