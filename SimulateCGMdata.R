# Create synthetic CGM data by using multi-level fPCA model combined with inverse-cdf transformation to mimick marginal subject-specific distributions
##############################################################
# Author: Irina Gaynanova
# Date of last update: June 10th, 2019

# NOT RUN: requires real CGM data, however the code will work with provided synthetic data

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


maxH = 7 # Maximal hour
tdexcom = round(maxH * 60 / 5) # Length of time grid in 5 minute intervals from the estimated sleep onset
n = length(minval) # Number of subjects

# Supporting libraries
library(tidyverse)
library(refund)

# Create a vector of curves and ids
Y <- c()
id <- c()
for (i in 1:n){
  # Extract estimated sleep trajectories for subject i
  raw_dexcom <- data[[i]]
  # Combine all curves together
  Y <- rbind(Y, raw_dexcom)
  # Combine all ids together
  id <- c(id, rep(subjectIDs[i], nrow(raw_dexcom)))
}

# Apply multi-level fPCA on all the data (TAKES TIME TO RUN, approx 5 min)
out_mfpca <- mfpca.sc(Y = Y, id = id, twoway = F)

# Calculate subject-level means and sds for level 2 scores
scores2 <- data.frame(scores = out_mfpca$scores$level2, subjects = id)

scoresSummaries <- scores2%>%
  group_by(subjects)%>%
  summarize(m1 = mean(scores.1), m2 = mean(scores.2), m3 = mean(scores.3), m4 = mean(scores.4), m5 = mean(scores.5), m6 = mean(scores.6), sd1 = sd(scores.1), sd2 = sd(scores.2), sd3 = sd(scores.3), sd4 = sd(scores.4), sd5 = sd(scores.5), sd6 = sd(scores.6))

# Estimate marginal parameters for multi-level functional Beta model, use mean/sd for normal data generation; alpha/beta/m/M for beta data generation
source("BetaEstimatingFunctions.R")
out_mfbeta <- estimateBetaParameters(data, minval, maxval)


# Create synthetic data with 14 curves per subject following either Normal or Beta marginal distributions
##################################################
simul_beta_real <- list()
simul_normal_real <- list()

nj = 14
nlevel1 = length(out_mfpca$evalues$level1)
nlevel2 = length(out_mfpca$evalues$level2) # number of scores needed for each night

# Modify min and max to check robustness of beta model
set.seed(23984)
minnew = pmax(minval + round(rnorm(n, mean = -20, sd = 10)), 40)
minnew = pmin(minnew, minval)
maxnew = pmin(maxval + round(rnorm(n, mean = 20, sd = 10)), 450)
maxnew = pmax(maxnew, maxval)


# Simulate n subjects with nj curves
set.seed(537349)
for (j in 1:n){
  # Draw level 1 scores
  ksi_means <- matrix(0, nlevel1, 1)
  for (i in 1:nlevel1){
    ksi_means[i,] <- rnorm(1, mean = 0, sd = sqrt(out_mfpca$evalues$level1[i]))
  }
  
  # Draw level 2 scores using a specific number of nights
  zeta_sds <- matrix(0, nlevel2, nj)
  for (i in 1:nlevel2){
    zeta_sds[i, ] <- rnorm(nj, mean = as.numeric(scoresSummaries[j,i + 1]), sd = as.numeric(scoresSummaries[j,i + 7])) 
  }
  
  # Generate nj curves strictly following mfPCA model
  subject_curves <- matrix(0, nj, tdexcom)
  for (i in 1:nj){
    subject_curves[i, ] = out_mfpca$mu + as.vector(out_mfpca$efunctions[[1]]%*%ksi_means) + as.vector(out_mfpca$efunctions[[2]]%*%zeta_sds[,i,drop = F])
  }
  
  # Add noise
  noise <- matrix(rnorm(nj*tdexcom, sd = sqrt(out_mfpca$sigma2)), nj, tdexcom)
  # For more realistic noise, set 2/3 of it to zero
  noise[sample(1:(nj*tdexcom), round(2*nj*tdexcom/3))] <- 0
  subject_curves <- subject_curves + noise
  
  # Generate uniforms by applying normal cdfs
  unif <- subject_curves
  subject_beta_curves <- subject_curves
  subject_normal_curves <- subject_curves
  means <- colMeans(subject_curves)
  sds <- apply(subject_curves, 2, sd)
  centered <- (subject_curves - matrix(means, nj, tdexcom, byrow = T))
  scaled <- centered %*% diag(1/sds)
  unif <- pnorm(scaled)
  
  # Apply inverse CDF to uniforms following either marginal normal or marginal beta distribution
  for (d in 1:nj){
    # Marginal normal model
    subject_normal_curves[d, ] <- qnorm(unif[d,], mean = out_mfbeta$means_smooth[j,], sd = out_mfbeta$sd_smooth[j,])
    # Marginal beta model
    out1 <- qbeta(unif[d,], shape1 = out_mfbeta$alphaM[j,], shape2 = out_mfbeta$betaM[j,])
    subject_beta_curves[d,] <- (maxnew[j] - minnew[j])*out1 + minnew[j] 
  }
  
  simul_beta_real[[j]] <- subject_beta_curves
  simul_normal_real[[j]] <- subject_normal_curves
}

# Save the generated synthetic data
# Min and max values for synthetic data, only apply to beta simulation
minvalSynthetic = minnew
maxvalSynthetic = maxnew
idsSynthetic = c(1:n)
save(simul_beta_real, simul_normal_real, minvalSynthetic, maxvalSynthetic, idsSynthetic, file = "CGMdataSynthetic.Rdata")



# Create figure of simulated data vs real data for 3 selected subjects
#####################################################################################

# NOT RUN: below only works with subjects ids from real data

# Path to save figures, use working directory by default
figure.path = getwd()

# Selected  subject ids
selectedIDs <- c(70058, 70134, 70186)
nid <- which(subjectIDs %in% selectedIDs)

# Code for figure
pdf(paste(figure.path, "Simulated_vs_real_examples.pdf", sep = ""), onefile = TRUE, width = 14, height = 8)

nj = 14
plotdata <- c()
for (j in nid){
  raw = data[[j]]
  Y = c(c(simul_normal_real[[j]]),c(simul_beta_real[[j]]), c(raw))
  plotdata <- rbind(plotdata, data.frame(Y = Y, day = c(rep(1:nj, tdexcom), rep(1:nj, tdexcom), rep(1:nrow(raw), tdexcom)) , time = c(rep(0:(tdexcom-1), each = nj),rep(0:(tdexcom-1), each = nj),rep(0:(tdexcom-1), each = nrow(raw))) * 5/60, type = c(rep("Normal", tdexcom*nj), rep("Beta", tdexcom*nj), rep("Real", tdexcom*nrow(raw))), subject = rep(subjectIDs[j], length(Y))))
}

p = ggplot(plotdata, aes(x = time, y = Y, group = as.factor(day))) + geom_line() + ylim(c(0,400)) + facet_grid(subject~factor(type,levels=c("Normal", "Beta", "Real"))) 
print(p)

dev.off()
