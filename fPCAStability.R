# Investigate the stability of principal eigenfunctions/scores with respect to subjects selection
#############################################################################################
# Author: Irina Gaynanova
# Date of last update: June 10th, 2019

# ATTENTION: results in manucscript are based on real CGM data, however the code will work with provided synthetic data (exception: use of real subject IDs for coloring of some of the figures)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

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

# Supporting libraries
library(tidyverse)
library(refund)

# Specify path to save figures, by default use working directory
figure.path = getwd()

# Estimate marginal parameters for multi-level functional Beta model, use mean/sd for normal data generation; alpha/beta/m/M for beta data generation
source("BetaEstimatingFunctions.R")
out_mfbeta <- estimateBetaParameters(data, minval, maxval)


##########################################################################################
# Figures with overall mean and two eigenfunctions for time-varying mean and standard deviation processes
##########################################################################################
ngrid = ncol(data[[1]])
time_grid_dexcom = c(0, 1:(ngrid-1) *5 /60)


# Plot of means and eigenfunctions for mean process
#########################################
# Look at all subjects mean
mean_fPCA_all <- data.frame(time = time_grid_dexcom, mu = out_mfbeta$pca_means$mu)
mean_fPCA_all$e1 <- out_mfbeta$pca_means$efunctions[,1]*sqrt(ngrid)
mean_fPCA_all$e2 <- out_mfbeta$pca_means$efunctions[,2]*sqrt(ngrid)
mean_fPCA_all$e3 <- out_mfbeta$pca_means$efunctions[,3]*sqrt(ngrid)

p1 <- mean_fPCA_all %>% 
  ggplot(aes(x = time, y = mu)) + geom_line() + xlab("Hours of sleep") + ylab("Glucose  [mg / dL]") + ggtitle("Overall glucose mean")
p1


# Helper functions to create facet labels
function_names <- list(
  'e1'="Mean eigenfunction 1",
  'e2'="Mean eigenfunction 2"
)

elabeller <- function(variable, value){
  return(function_names[value])
}

# Mean eigenfunctions plot
p2 <- mean_fPCA_all %>%
  dplyr::select(-mu, -e3)%>%
  gather(key = "type", value = "value", -time)%>%
  ggplot(aes(x = time, y = value)) + facet_grid(~type, labeller = elabeller) + geom_line()+ xlab("Hours of sleep") + ylab("")
p2

pdf(paste(figure.path, "MeanPCA.pdf", sep = ""), onefile = F, width = 12, height = 4)
multiplot(p1,p2, cols = 2)
dev.off()

# Plot of means and eigenfunctions for sd process
#########################################
# Create sd data frame
sd_fPCA_all <- data.frame(time = time_grid_dexcom, mu = out_mfbeta$pca_sds$mu)
sd_fPCA_all$e1 <- out_mfbeta$pca_sds$efunctions[,1]*sqrt(ngrid)
sd_fPCA_all$e2 <- out_mfbeta$pca_sds$efunctions[,2]*sqrt(ngrid)
sd_fPCA_all$e3 <- out_mfbeta$pca_sds$efunctions[,3]*sqrt(ngrid)

# Overall mean plot for sd process
p1 <- sd_fPCA_all %>% 
  dplyr::select(time, mu) %>% 
  ggplot(aes(x = time, y = mu)) + geom_line() + xlab("Hours of sleep") + ylab("[mg / dL]") + ggtitle("Overall glucose standard deviation")
p1

# Helper functions to create facet labels
function_names <- list(
  'e1'="Sd eigenfunction 1",
  'e2'="Sd eigenfunction 2"
)

elabeller <- function(variable, value){
  return(function_names[value])
}

# Sd eigenfunctions plot
p2 <- sd_fPCA_all %>%
  dplyr::select(-mu, -e3)%>%
  gather(key = "type", value = "value", -time)%>%
  ggplot(aes(x = time, y = value)) + facet_grid(~type, labeller = elabeller) + geom_line()+ xlab("Hours of sleep") + ylab("")
p2

pdf(paste(figure.path, "SdPCA.pdf", sep = ""), onefile = F, width = 12, height = 4)
multiplot(p1,p2, cols = 2)
dev.off()

###########################################################################################
# Create a figure for pairwise scores and calculate pairwise correlations
###########################################################################################

# Extract all the scores
scores_means = out_mfbeta$pca_means$scores
scores_sds = out_mfbeta$pca_sds$scores

# Combine all the scores and the subject
dexcom = data.frame(id = subjectIDs, scores_sd1 = scores_sds[,1], scores_sd2 = scores_sds[,2], scores_sd3 = scores_sds[,3], scores1 = scores_means[,1],scores2 = scores_means[,2], scores3 = scores_means[,3])

# Calculate correlations between pairs of scores
cor(dexcom$scores1, dexcom$scores_sd1) # 0.61
cor(dexcom$scores2, dexcom$scores_sd2) # 0.47
cor(dexcom$scores3, dexcom$scores_sd3) # 0.44

# Figures for representative 6 subjects selected
# NOT RUN: selected ID numbers for 6 selected subjects from real CGM data only
selectedID <- c(70058, 70134, 70139, 70142, 70153, 70186)
dexcom_selected <- dexcom[dexcom$id %in% selectedID,]
dexcom_selected$id <- as.factor(as.character(dexcom_selected$id))


# Create a figure of pairwise scores, use different color for each subject. Only do the first 2 scores.
pdf(paste(figure.path, "PairwiseScores_col.pdf", sep = ""), onefile = TRUE, width = 12, height = 5)
p1 = dexcom %>% ggplot(aes(scores1, scores_sd1)) + geom_point(alpha = 0.3, size = 2) + xlab("Mean score 1") + ylab("Sd score 1") + geom_point(aes(scores1, scores_sd1, col = id), data = dexcom_selected, size = 3)+ theme(legend.position="none")
p2 = dexcom %>% ggplot(aes(scores2, scores_sd2)) + geom_point(alpha = 0.3, size = 2) + xlab("Mean score 2") + ylab("Sd score 2") + geom_point(aes(scores2, scores_sd2, col = id), data = dexcom_selected, size = 3)+ theme(legend.position="none")
multiplot(p1,p2, cols = 2)
dev.off()


########################################################################################
# Evaluate stability by 
# (i) randomly taking a subset of subjects out; 
# (ii) recalculating smoothed means/sds and eigenfunctions;
# (iii) plot eigenfunctions across replications with overall eigenfunction overlayed
########################################################################################
# Number of splits
nrep = 100
# How many subjects to take out?
n = length(subjectIDs)
percent = 0.8
nselect = round(n * percent) 
# Number of conponents to extract
ncomp = 3
# Store first two mean eigenfunctions + overall mean

overall_mean <- matrix(0, nrep, ngrid)
eigenfunction_mean_1 <- matrix(0, nrep, ngrid)
eigenfunction_mean_2 <- matrix(0, nrep, ngrid)
# Store first two sd eigenfunctions + overall sd
overall_sd <- matrix(0, nrep, ngrid)
eigenfunction_sd_1 <- matrix(0, nrep, ngrid)
eigenfunction_sd_2 <- matrix(0, nrep, ngrid)
# Start the loop 
set.seed(23946)
for (i in 1:nrep){
  # Choose nselect subjects out of n
  select_id = sample(1:n, nselect)
  # Apply estimation procedure on remaining subject
  out_mfbeta_subsample <- estimateBetaParameters(data[select_id], minval[select_id], maxval[select_id])
  
  # Get overall mean and eigenfunctions for mean process
  overall_mean[i, ] <- out_mfbeta_subsample$pca_means$mu
  eigenfunction_mean_1[i, ] <- out_mfbeta_subsample$pca_means$efunctions[,1]*sqrt(ngrid)
  eigenfunction_mean_2[i, ] <- out_mfbeta_subsample$pca_means$efunctions[,2]*sqrt(ngrid)
  
  # Get overall mean and eigenfunctions for sd process
  overall_sd[i, ] <- out_mfbeta_subsample$pca_sds$mu
  eigenfunction_sd_1[i, ] <- out_mfbeta_subsample$pca_sds$efunctions[,1]*sqrt(ngrid)
  eigenfunction_sd_2[i, ] <- out_mfbeta_subsample$pca_sds$efunctions[,2]*sqrt(ngrid)
}

########################################################################################
# Make plots of stability
#######################################################################################
# Create overall mean data frame
mean_fPCA <- melt(t(overall_mean))
colnames(mean_fPCA) <- c("time", "nrep", "mu")
mean_fPCA$time <- rep(time_grid_dexcom, nrep)


# Compare all subjects with subjects removed plot
p1 = mean_fPCA %>% ggplot(aes(x = time, y = mu, group = nrep)) + geom_line(alpha = 0.8) + xlab("Hours of sleep") + ylab("Glucose  [mg / dL]") + ggtitle("Overall glucose mean") + geom_line(aes(x = time, y = mu), data = mean_fPCA_all, col = "red", lwd = 1.5)
p1



# Mean eigenfunctions plot, separately for each eigenfunction
########################################
# Eigenfunction 1
########################################
# Flip the sign in accordance with overall eigenfunction
# Calculate inner-product
innerprod <- eigenfunction_mean_1 %*% mean_fPCA_all$e1
# For negative ones, flip the sign
eigenfunction_mean_1[innerprod < 0, ] = - eigenfunction_mean_1[innerprod < 0, ]


# Create a data frame
e1_fPCA <- melt(t(eigenfunction_mean_1))
colnames(e1_fPCA) <- c("time", "nrep", "e1")
e1_fPCA$time <- rep(time_grid_dexcom, nrep)

p2 <- e1_fPCA %>%
  ggplot(aes(x = time, y = e1, group = nrep)) + geom_line(alpha = 0.8)+ xlab("Hours of sleep") + ylab("") + geom_line(aes(x = time, y = e1), data = mean_fPCA_all, col = "red", lwd = 1.5)+ ylim(c(-1.7, 1.7))+ ggtitle("Mean eigenfunction 1")
p2

########################################
# Eigenfunction 2
########################################
# Flip the sign in accordance with overall eigenfunction
# Calculate inner-product
innerprod <- eigenfunction_mean_2 %*% mean_fPCA_all$e2
# For negative ones, flip the sign
eigenfunction_mean_2[innerprod < 0, ] = - eigenfunction_mean_2[innerprod < 0, ]


# Create a data frame
e2_fPCA <- melt(t(eigenfunction_mean_2))
colnames(e2_fPCA) <- c("time", "nrep", "e2")
e2_fPCA$time <- rep(time_grid_dexcom, nrep)

p3 <- e2_fPCA %>%
  ggplot(aes(x = time, y = e2, group = nrep)) + geom_line(alpha = 0.8)+ xlab("Hours of sleep") + ylab("") + geom_line(aes(x = time, y = e2), data = mean_fPCA_all, col = "red", lwd = 1.5)+ ylim(c(-1.7, 1.7))+ ggtitle("Mean eigenfunction 2")
p3

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

pdf(paste(figure.path, "MeanPCA_stability.pdf", sep = ""), onefile = F, width = 12, height = 4)
multiplot(p1,p2,p3, cols = 3)
dev.off()

# ########################################
# # Repeat for standard deviation process
# ########################################

# Melt sd data frame for overall sd mean
sd_fPCA <- melt(t(overall_sd))
colnames(sd_fPCA) <- c("time", "nrep", "mu")
sd_fPCA$time <- rep(time_grid_dexcom, nrep)

# Compare all subjects with subjects removed plot
p1 = sd_fPCA %>% ggplot(aes(x = time, y = mu, group = nrep)) + geom_line(alpha = 0.8) + xlab("Hours of sleep") + ylab("Glucose  [mg / dL]") + ggtitle("Overall glucose sd") + geom_line(aes(x = time, y = mu), data = sd_fPCA_all, col = "red", lwd = 1.5)
p1

# sd eigenfunctions plot, separately for each eigenfunction
########################################
# Sd Eigenfunction 1
########################################
# Flip the sign in accordance with overall eigenfunction
# Calculate inner-product
innerprod <- eigenfunction_sd_1 %*% sd_fPCA_all$e1
# For negative ones, flip the sign
eigenfunction_sd_1[innerprod < 0, ] = - eigenfunction_sd_1[innerprod < 0, ]


# Create a data frame
e1_fPCA <- melt(t(eigenfunction_sd_1))
colnames(e1_fPCA) <- c("time", "nrep", "e1")
e1_fPCA$time <- rep(time_grid_dexcom, nrep)

p2 <- e1_fPCA %>%
  ggplot(aes(x = time, y = e1, group = nrep)) + geom_line(alpha = 0.8)+ xlab("Hours of sleep") + ylab("") + geom_line(aes(x = time, y = e1), data = sd_fPCA_all, col = "red", lwd = 1.5)+ ylim(c(-1.9, 1.6))+ ggtitle("Sd eigenfunction 1")
p2

########################################
# Eigenfunction 2
########################################
# Flip the sign in accordance with overall eigenfunction
# Calculate inner-product
innerprod <- eigenfunction_sd_2 %*% sd_fPCA_all$e2
# For negative ones, flip the sign
eigenfunction_sd_2[innerprod < 0, ] = - eigenfunction_sd_2[innerprod < 0, ]


# Create a data frame
e2_fPCA <- melt(t(eigenfunction_sd_2))
colnames(e2_fPCA) <- c("time", "nrep", "e2")
e2_fPCA$time <- rep(time_grid_dexcom, nrep)

p3 <- e2_fPCA %>%
  ggplot(aes(x = time, y = e2, group = nrep)) + geom_line(alpha = 0.8)+ xlab("Hours of sleep") + ylab("") + geom_line(aes(x = time, y = e2), data = sd_fPCA_all, col = "red", lwd = 1.5)+ ylim(c(-1.9, 1.6))+ ggtitle("Sd eigenfunction 2")
p3


pdf(paste(figure.path, "SdPCA_Stability.pdf", sep = ""), onefile = F, width = 12, height = 4)
multiplot(p1, p2, p3, cols = 3)
dev.off()
