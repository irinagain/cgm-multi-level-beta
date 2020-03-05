# R code for creating figures for the paper on "Modeling CGM trajectories during sleep"
########################################################################################
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

library(tidyverse)
library(refund)
library(reshape2)

# Specify path to save figures, by default use working directory
figure.path = getwd()
text_size = 20 # font size for plotting

#######################################################################################
# Figure that visualizes the ranges of [m_i, M_i] sorted by the value of M_i
#######################################################################################

# Create data frame for plotting
ranges <- data.frame(id = subjectIDs, min = minval, max = maxval)
orderM <- order(maxval)
ranges$idnew <- orderM
ranges$idnew[orderM] <- c(1:length(orderM))

pdf(paste(figure.path, "/mMplot.pdf", sep = ""), onefile = F, width = 16, height = 4)
p = ranges %>% gather(key = "type", value = "value", -idnew, -id)%>% ggplot(aes(x = idnew, y = value, group = idnew)) + geom_line() + geom_point() + xlab("Subject") + ylab("Glucose [mg/dL]")+ theme(text = element_text(size = text_size))
print(p)
dev.off()

# Fit multi-level functional Beta model
#######################################################################################
source("BetaEstimatingFunctions.R")
out_mfbeta <- estimateBetaParameters(data, minval, maxval)

##########################################################################################
# Figure with (1) pointwise versus smoothed means for 6 selected subjects 
#             (2) smoothed means for all subjects with 6 highlighted in blue
##########################################################################################

# Time grid
grid_length = ncol(data[[1]])
time_grid_dexcom = c(0, 1:(grid_length-1) *5 /60)

# NOT RUN: below only works with subjects ids from real data

# ID numbers for 6 selected subjects
selectedID <- c(70058, 70134, 70139, 70142, 70153, 70186)


# Extract smoothed means for 6 selected subjects
selectedMeans <- as.data.frame(t(out_mfbeta$pca_means$Yhat[subjectIDs %in% selectedID,]))
selectedMeans$t <- time_grid_dexcom
selectedMeans <- gather(selectedMeans, key = idnew, value = smean, -7)

# Extract raw pointwise means
raw_means <- t(sapply(data, function(x) colMeans(x, na.rm = T)))
selectedRawMeans <- as.data.frame(t(raw_means[subjectIDs %in% selectedID,]))
selectedRawMeans$t <- time_grid_dexcom
selectedRawMeans <- gather(selectedRawMeans, key = idnew, value = rawmean, -7)

# Pointwise versus smoothed for 6 selected subjects
pdf(paste(figure.path, "/PointwiseVsSmoothedMeans.pdf", sep = ""), onefile = F, width = 6, height = 6)
p = ggplot(selectedMeans,aes(x = t, y = smean)) + geom_line(aes(group=idnew), col = "blue")  + geom_line(aes(x = t, y = rawmean, group=idnew), data = selectedRawMeans, col = "red")+ xlab("Hours of sleep") + ylab("Glucose  [mg / dL]") + ylim(c(50, 300))+ theme(text = element_text(size = 16))
print(p)
dev.off()  

# Consider all smoothed means
rownames(out_mfbeta$pca_means$Yhat) = subjectIDs
smoothed_means <- melt(out_mfbeta$pca_means$Yhat)
colnames(smoothed_means) <- c("id", "sleeptime", "gl")
smoothed_means$sleeptime <- (smoothed_means$sleeptime - 1) * 5/60

# Create a figure for all smoothed means
pdf(paste(figure.path, "/AllSmoothedMeans.pdf", sep = ""), onefile = F, width = 6, height = 6)
p = ggplot(smoothed_means,aes(x = sleeptime, y = gl)) + geom_line(aes(group=id), col = "grey")  + geom_line(aes(group=id), data = smoothed_means[smoothed_means$id %in% selectedID,], col = "blue")+ xlab("Hours of sleep") + ylab("Glucose  [mg / dL]") + ylim(c(50, 300))+ theme(text = element_text(size = 16))
print(p)
dev.off()  



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
  subject_data <- melt(data[[i]], value.name = "gl")
  subject_data[,2] <- subject_data[, 2] * 5/60
  colnames(subject_data) <- c("sleepday", "sleeptime", "gl")
  subject_data$id <- subjectIDs[i]
  datagl <- rbind(datagl, subject_data)
}

datagl$idnew <- as.factor(as.character(datagl$id))

# Plot of representative subjects without any bands
pdf(paste(figure.path, "/RepresentativeS.pdf", sep = ""), onefile = F, width = 16, height = 12)
p = ggplot(datagl,aes(x = as.numeric(sleeptime), y = gl)) + geom_line(aes(group=sleepday)) + facet_wrap(~as.character(idnew), ncol = 2)+ xlim(c(0,7))+ xlab("Hours of sleep") + ylab("Glucose  [mg / dL]") + ylim(c(50, 400))+ theme(legend.position="none") + theme(text = element_text(size = text_size))
print(p)
dev.off()


# Extract all the scores
scores_means = out_mfbeta$pca_means$scores
scores_sds = out_mfbeta$pca_sds$scores
# Combine all the scores and the subject
dexcom = data.frame(id = subjectIDs, scores_sd1 = -scores_sds[,1], scores_sd2 = scores_sds[,2], scores_sd3 = scores_sds[,3], scores1 = scores_means[,1],scores2 = scores_means[,2], scores3 = scores_means[,3])
dexcom_selected <- dexcom[dexcom$id %in% selectedID,]

# Change labels to show the values of the scores
subject_names <- c(
  "70058"=paste("70058, mean scores (", round(dexcom_selected$scores1[1],1),", ",round(dexcom_selected$scores2[1],1),"), ","sd scores (",round(dexcom_selected$scores_sd1[1],1),", ",round(dexcom_selected$scores_sd2[1],1),")", sep=""),
  '70134'=paste("70134, mean scores (", round(dexcom_selected$scores1[2],1),", ",round(dexcom_selected$scores2[2],1),"), ","sd scores (",round(dexcom_selected$scores_sd1[2],1),", ",round(dexcom_selected$scores_sd2[2],1),")", sep=""),
  '70139'=paste("70139, mean scores (", round(dexcom_selected$scores1[3],1),", ",round(dexcom_selected$scores2[3],1),"), ","sd scores (",round(dexcom_selected$scores_sd1[3],1),", ",round(dexcom_selected$scores_sd2[3],1),")", sep=""),
  '70142'=paste("70142, mean scores (", round(dexcom_selected$scores1[4],1),", ",round(dexcom_selected$scores2[4],1),"), ","sd scores (",round(dexcom_selected$scores_sd1[4],1),", ",round(dexcom_selected$scores_sd2[4],1),")", sep=""),
  '70153'=paste("70153, mean scores (", round(dexcom_selected$scores1[5],1),", ",round(dexcom_selected$scores2[5],1),"), ","sd scores (",round(dexcom_selected$scores_sd1[5],1),", ",round(dexcom_selected$scores_sd2[5],1),")", sep=""),
  '70186'=paste("70186, mean scores (", round(dexcom_selected$scores1[6],1),", ",round(dexcom_selected$scores2[6],1),"), ","sd scores (",round(dexcom_selected$scores_sd1[6],1),", ",round(dexcom_selected$scores_sd2[6],1),")", sep="")
)

selectedAlphas$idnew <- as.factor(selectedAlphas$idnew)
levels(selectedAlphas$idnew) <- levels(datagl$idnew)
# Get minimal and maximal values
order <- match(selectedAlphas$idnew, subjectIDs)
selectedAlphas$minval = minval[order]
selectedAlphas$maxval = maxval[order]

# Plot of representative subjects with quantile bands and scores from PCA
pdf(paste(figure.path, "/RepresentativeSBeta_Col.pdf", sep = ""), onefile = F, width = 16, height = 12)
p = ggplot(datagl,aes(x = as.numeric(sleeptime), y = gl)) + geom_line(aes(group=sleepday)) + geom_line(aes(x = t, y = minval + (maxval-minval)*qbeta(1-0.05/2, alpha, beta), col = idnew), data = selectedAlphas, size = 1.2)+ geom_line(aes(x = t, y = minval + (maxval-minval)*qbeta(0.05/2, alpha, beta), col = idnew), data = selectedAlphas, size = 1.2)+ geom_line(aes(x = t, y = minval + (maxval-minval)*qbeta(0.5, alpha, beta), col = idnew), data = selectedAlphas, size = 1.2)+ facet_wrap(~as.character(idnew), ncol = 2, labeller = as_labeller(subject_names))+ xlim(c(0,7))+ xlab("Hours of sleep") + ylab("Glucose  [mg / dL]") + ylim(c(50, 400))+ theme(legend.position="none")+ theme(text = element_text(size = text_size))
print(p)
dev.off()


# Fit multi-level functional Beta model with 3 components for mean and 6 for sd
#######################################################################################
source("BetaEstimatingFunctions.R")
out_mfbeta <- estimateBetaParameters(data, minval, maxval, npc_sd = 6)

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
  subject_data <- melt(data[[i]], value.name = "gl")
  subject_data[,2] <- subject_data[, 2] * 5/60
  colnames(subject_data) <- c("sleepday", "sleeptime", "gl")
  subject_data$id <- subjectIDs[i]
  datagl <- rbind(datagl, subject_data)
}

datagl$idnew <- as.factor(as.character(datagl$id))


# Extract all the scores
scores_means = out_mfbeta$pca_means$scores
scores_sds = out_mfbeta$pca_sds$scores
# Combine all the scores and the subject
# Use -scores for sd 1st component to match orientation for the mean for easer interpretation
dexcom = data.frame(id = subjectIDs, scores_sd1 = -scores_sds[,1], scores_sd2 = scores_sds[,2], scores_sd3 = scores_sds[,3], scores1 = scores_means[,1],scores2 = scores_means[,2], scores3 = scores_means[,3])
dexcom_selected <- dexcom[dexcom$id %in% selectedID,]

# Change labels to show the values of the scores
subject_names <- c(
  "70058"=paste("70058, mean scores (", round(dexcom_selected$scores1[1],1),", ",round(dexcom_selected$scores2[1],1),"), ","sd scores (",round(dexcom_selected$scores_sd1[1],1),", ",round(dexcom_selected$scores_sd2[1],1),")", sep=""),
  '70134'=paste("70134, mean scores (", round(dexcom_selected$scores1[2],1),", ",round(dexcom_selected$scores2[2],1),"), ","sd scores (",round(dexcom_selected$scores_sd1[2],1),", ",round(dexcom_selected$scores_sd2[2],1),")", sep=""),
  '70139'=paste("70139, mean scores (", round(dexcom_selected$scores1[3],1),", ",round(dexcom_selected$scores2[3],1),"), ","sd scores (",round(dexcom_selected$scores_sd1[3],1),", ",round(dexcom_selected$scores_sd2[3],1),")", sep=""),
  '70142'=paste("70142, mean scores (", round(dexcom_selected$scores1[4],1),", ",round(dexcom_selected$scores2[4],1),"), ","sd scores (",round(dexcom_selected$scores_sd1[4],1),", ",round(dexcom_selected$scores_sd2[4],1),")", sep=""),
  '70153'=paste("70153, mean scores (", round(dexcom_selected$scores1[5],1),", ",round(dexcom_selected$scores2[5],1),"), ","sd scores (",round(dexcom_selected$scores_sd1[5],1),", ",round(dexcom_selected$scores_sd2[5],1),")", sep=""),
  '70186'=paste("70186, mean scores (", round(dexcom_selected$scores1[6],1),", ",round(dexcom_selected$scores2[6],1),"), ","sd scores (",round(dexcom_selected$scores_sd1[6],1),", ",round(dexcom_selected$scores_sd2[6],1),")", sep="")
)

selectedAlphas$idnew <- as.factor(selectedAlphas$idnew)
levels(selectedAlphas$idnew) <- levels(datagl$idnew)
# Get minimal and maximal values
order <- match(selectedAlphas$idnew, subjectIDs)
selectedAlphas$minval = minval[order]
selectedAlphas$maxval = maxval[order]

# Plot of representative subjects with quantile bands and scores from PCA
pdf(paste(figure.path, "/RepresentativeSBeta_Col_6compsd.pdf", sep = ""), onefile = F, width = 16, height = 12)
p = ggplot(datagl,aes(x = as.numeric(sleeptime), y = gl)) + geom_line(aes(group=sleepday)) + geom_line(aes(x = t, y = minval + (maxval-minval)*qbeta(1-0.05/2, alpha, beta), col = idnew), data = selectedAlphas, size = 1.2)+ geom_line(aes(x = t, y = minval + (maxval-minval)*qbeta(0.05/2, alpha, beta), col = idnew), data = selectedAlphas, size = 1.2)+ geom_line(aes(x = t, y = minval + (maxval-minval)*qbeta(0.5, alpha, beta), col = idnew), data = selectedAlphas, size = 1.2)+ facet_wrap(~as.character(idnew), ncol = 2, labeller = as_labeller(subject_names))+ xlim(c(0,7))+ xlab("Hours of sleep") + ylab("Glucose  [mg / dL]") + ylim(c(50, 400))+ theme(legend.position="none")+ theme(text = element_text(size = text_size))
print(p)
dev.off()