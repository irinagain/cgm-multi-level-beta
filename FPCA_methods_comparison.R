# Comparison of different FPCA methods
################################################
# Author: Irina Gaynanova

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


figure.path = getwd()

train = data
require(refund)
n = length(train)
maxH = 7 # Set the maximal hour
tdexcom = round(maxH * 60 / 5)
ngrid = tdexcom

# Calculate pointwise means
##############################################
raw_means = matrix(NA, n, tdexcom)
nsleep = rep(NA, n)
for (i in 1:n){
  # Extract estimated sleep trajectories for subject i
  raw_dexcom <- train[[i]][, 1:tdexcom]
  # Calculate pointwise means
  raw_means[i,] <- colMeans(raw_dexcom, na.rm = T)
  # Calculate number of sleep periods
  nsleep[i] = nrow(train[[i]])
}
summary(nsleep)

# Weighted methods
######################################################
# Make explicit centering using smoothing splines
argvals <- (1:tdexcom)/tdexcom
meanX <- colMeans(raw_means, na.rm = TRUE)
meanXsmooth <- smooth.spline(argvals, meanX, all.knots = TRUE)$y
Ycentered <- t(t(raw_means) - meanXsmooth)

#Alternative centering to be consistent with SC method
nbasis = 15
d.vec = rep(argvals, each = nrow(raw_means))
gam0 = mgcv::gam(as.vector(raw_means) ~ s(d.vec, bs = "cr", k = nbasis))
mu = predict(gam0, newdata = data.frame(d.vec = argvals))
mu = as.vector(mu)
Ycentered2 <- t(t(raw_means) - mu)

# Construct the weights as nI/sum(nI) * sqrt(I) (I - total numver of subjects, probably should be n)
weights = sqrt(nsleep * n / sum(nsleep)) # from 0.62 to 1.2, the larger nsleep, the larger is the weight

Ywcentered = diag(weights) %*% Ycentered
Ywcentered2 = diag(weights) %*% Ycentered2

pca_means_w <- fpca.face(Y = Ywcentered, center = F, argvals=(1:tdexcom)/tdexcom, var = T, knots = 5) # 3 components are selected to explain 99% of variance; error variance is estimated at 0.97
pca_means_w$npc #3
pca_means_w$sigma2 #3.7

pca_means_default_w <- fpca.face(Y = Ywcentered, center = F, argvals=(1:tdexcom)/tdexcom, var = T) # 3 components are selected to explain 99% of variance; error variance is estimated at 0.97
pca_means_default_w$npc #3
pca_means_default_w$sigma2 #0.86
# Very very similar to unweighted version, still 3 components, very little difference in scores

pca_means_default2_w <- fpca.sc(Y = Ywcentered2, center = F, argvals=(1:tdexcom)/tdexcom, var = T)
pca_means_default2_w$npc # 3
pca_means_default2_w$sigma2 # 11.9


# Unweighted methods
######################################################
pca_means_default <- fpca.face(Y = raw_means, center = T, argvals=(1:tdexcom)/tdexcom, var = T) # 3 components are selected to explain 99% of variance; error variance is estimated at 0.97
pca_means_default$npc
pca_means_default$sigma2

# Calculates values differently, appear to be a magnitude smaller than in face. Difference appear to be the magnitude of grid length (84). This centers using gam rather than smooth.spline
pca_means_default2 <- fpca.sc(Y = raw_means, center = T, argvals=(1:tdexcom)/tdexcom, var = T)
pca_means_default2$npc # 3
pca_means_default2$sigma2 # 13

# Between the two methods, the means and function look pretty much the same, although the second one is more smooth. I wonder how that compares to when I use less knots below. Very very similar
pca_means_default3 <- fpca.sc(Y = raw_means, center = T, argvals=(1:tdexcom)/tdexcom, var = T, useSymm = T)
pca_means_default3$npc # 4
pca_means_default3$sigma2 # 7.4

ncomp = 3
pca_means <- fpca.face(Y = raw_means, center = T, argvals=(1:tdexcom)/tdexcom, knots=5, npc = ncomp, var = T) # error variance is estimated at 4


# Compare fits on all eigenfunctions between FACE (5 knots vs 35 knots) and SC
###############################################################################################################
time_grid_dexcom = c(0, 1:(tdexcom-1) * 5 /60)

fpca_face5 <- data.frame(time = time_grid_dexcom, mu = pca_means$mu, method = rep("FACE, 5 knots", tdexcom))
fpca_face5$e1 <- pca_means$efunctions[,1]*sqrt(ngrid)
fpca_face5$e2 <- pca_means$efunctions[,2]*sqrt(ngrid)
fpca_face5$e3 <- pca_means$efunctions[,3]*sqrt(ngrid)


fpca_face5w <- data.frame(time = time_grid_dexcom, mu = meanXsmooth, method = rep("FACEw, 5 knots", tdexcom))
fpca_face5w$e1 <- pca_means_w$efunctions[,1]*sqrt(ngrid)
fpca_face5w$e2 <- pca_means_w$efunctions[,2]*sqrt(ngrid)
fpca_face5w$e3 <- pca_means_w$efunctions[,3]*sqrt(ngrid)

fpca_face35 <- data.frame(time = time_grid_dexcom, mu = pca_means_default$mu, method = rep("FACE, 35 knots", tdexcom))
fpca_face35$e1 <- -pca_means_default$efunctions[,1]*sqrt(ngrid)
fpca_face35$e2 <- -pca_means_default$efunctions[,2]*sqrt(ngrid)
fpca_face35$e3 <- pca_means_default$efunctions[,3]*sqrt(ngrid)

fpca_face35w <- data.frame(time = time_grid_dexcom, mu = meanXsmooth, method = rep("FACEw, 35 knots", tdexcom))
fpca_face35w$e1 <- -pca_means_default_w$efunctions[,1]*sqrt(ngrid)
fpca_face35w$e2 <- -pca_means_default_w$efunctions[,2]*sqrt(ngrid)
fpca_face35w$e3 <- pca_means_default_w$efunctions[,3]*sqrt(ngrid)

fpca_sc10 <- data.frame(time = time_grid_dexcom, mu = pca_means_default2$mu, method = rep("SC, 10 basis", tdexcom))
fpca_sc10$e1 <- -pca_means_default2$efunctions[,1]
fpca_sc10$e2 <- pca_means_default2$efunctions[,2]
fpca_sc10$e3 <- -pca_means_default2$efunctions[,3]

fpca_sc10w <- data.frame(time = time_grid_dexcom, mu = mu, method = rep("SCw, 10 basis", tdexcom))
fpca_sc10w$e1 <- -pca_means_default2_w$efunctions[,1]
fpca_sc10w$e2 <- pca_means_default2_w$efunctions[,2]
fpca_sc10w$e3 <- -pca_means_default2_w$efunctions[,3]

mean_fPCA_all <- rbind(fpca_face5, fpca_face5w, fpca_face35, fpca_face35w, fpca_sc10, fpca_sc10w)

p1 <- mean_fPCA_all %>% 
  ggplot(aes(x = time, y = mu, col = as.factor(method))) + geom_line() + xlab("Hours of sleep") + ylab("Glucose  [mg / dL]") + ggtitle("Overall glucose mean") + theme(legend.position = "none") 
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
  gather(key = "type", value = "value", -c(time, method))%>%
  ggplot(aes(x = time, y = value, col = as.factor(method))) + facet_grid(~type, labeller = elabeller) + geom_line()+ xlab("Hours of sleep") + ylab("")+ labs(col = "Method")
p2

pdf(paste(figure.path, "/MeanPCA_methods_comparison.pdf", sep = ""), onefile = F, width = 12, height = 4)
multiplot(p1, p2, layout = matrix(c(1,2,2), ncol=3, byrow=TRUE))
dev.off()

# Conclusion: for the mean eigenfunctions, the shape is so simple that it doesn't matter how many knots, and there is very high agreement with SC and FACE.



# Compare fits on all fitted means between FACE (5 knots vs 35 knots) and SC
##########################################################################################
rownames(pca_means$Yhat) = subjectIDs
smoothed_means <- reshape2::melt(pca_means$Yhat)
colnames(smoothed_means) <- c("id", "sleeptime", "gl")
smoothed_means$sleeptime <- (smoothed_means$sleeptime - 1) * 5/60
smoothed_means$method <- "FACE, 5 knots"

rownames(pca_means_default$Yhat) = subjectIDs
smoothed_means2 <- reshape2::melt(pca_means_default$Yhat)
colnames(smoothed_means2) <- c("id", "sleeptime", "gl")
smoothed_means2$sleeptime <- (smoothed_means2$sleeptime - 1) * 5/60
smoothed_means2$method <- "FACE, 35 knots"

rownames(pca_means_default2$Yhat) = subjectIDs
smoothed_means3 <- reshape2::melt(pca_means_default2$Yhat)
colnames(smoothed_means3) <- c("id", "sleeptime", "gl")
smoothed_means3$sleeptime <- (smoothed_means3$sleeptime - 1) * 5/60
smoothed_means3$method <- "SC, 10 basis"

rownames(raw_means) = subjectIDs
smoothed_means4 <- reshape2::melt(raw_means)
colnames(smoothed_means4) <- c("id", "sleeptime", "gl")
smoothed_means4$sleeptime <- (smoothed_means4$sleeptime - 1) * 5/60
smoothed_means4$method <- "Raw"

rownames(pca_means_w$Yhat) = subjectIDs
pca_means_w$Yhat = diag(1/weights) %*% pca_means_w$Yhat + matrix(meanXsmooth, nrow(pca_means_w$Yhat), ncol(pca_means_w$Yhat), byrow = T)
smoothed_means5 <- reshape2::melt(pca_means_w$Yhat)
colnames(smoothed_means5) <- c("id", "sleeptime", "gl")
smoothed_means5$sleeptime <- (smoothed_means5$sleeptime - 1) * 5/60
smoothed_means5$method <- "FACEw, 5 knots"


smoothed_means_all <- rbind(smoothed_means, smoothed_means2, smoothed_means3, smoothed_means5)

p = ggplot(smoothed_means_all,aes(x = sleeptime, y = gl, group = interaction(id, method), col = method)) + geom_line() + xlab("Hours of sleep") + ylab("Glucose  [mg / dL]") + ylim(c(50, 300))


# Create a figure for all smoothed means - all look same essentially, but 35 knots overfit

pdf(paste(figure.path, "/AllSmoothedMeansComparison.pdf", sep = ""), onefile = F, width = 6, height = 10)
print(p)
dev.off()  

###############################################################################################################
# Comparison for sd processes
###############################################################################################################
dexcom_sds = matrix(NA, n, tdexcom)
dexcom_sds_point = matrix(NA, n, tdexcom)
rownames(dexcom_sds) = subjectIDs
rownames(dexcom_sds_point) = rownames(dexcom_sds)
for (i in 1:n){
  # Extract estimated sleep trajectories for subject i
  raw_dexcom <- data[[i]]
  # Subtract smoothed mean from raw data
  subtract_mean <- raw_dexcom - matrix(pca_means$Yhat[i,], nrow = nrow(raw_dexcom), ncol = ncol(raw_dexcom), byrow = T)
  # Estimate st.dev on the residuals
  dexcom_sds[i, ] <- apply(subtract_mean, 2, function(x) sqrt(sum(x^2, na.rm = T)/(sum(!is.na(x))-1))) 
  # Also estimate standard normal pointwise mean
  dexcom_sds_point[i, ] <- apply(raw_dexcom, 2, sd, na.rm = T)
}
dexcom_sds[dexcom_sds == Inf] = NA


########################################################
# Smooth pointwise sds
########################################################
ncomp = 3
pca_sds <- fpca.face(Y = dexcom_sds, center = T, argvals=(1:tdexcom)/tdexcom, knots = 10, npc = ncomp) # 3

pca_sds_default <- fpca.face(Y = dexcom_sds, center = T, argvals=(1:tdexcom)/tdexcom, knots = 35, npc = ncomp) # still 3

pca_sds_default2 <- fpca.face(Y = dexcom_sds, center = T, argvals=(1:tdexcom)/tdexcom, knots = 5, pve = 0.9999999) # 9 components with 35 knots, 8 components with 20 knots, going below 6 gives less than 1%, first 3 explain 93%, with 4th it's 95%

pca_sds_default3 <- fpca.face(Y = dexcom_sds, center = T, argvals=(1:tdexcom)/tdexcom, knots = 5) # 6 components with 5 or 7 knots, 7 components with 10 knots (7 and 8 are less than 0.5%, although there are extra 2.4% in 4th, and 1.3% in 5th)

### Pointwise
pca_sds_point <- fpca.face(Y = dexcom_sds_point, center = T, argvals=(1:tdexcom)/tdexcom, knots = 10, npc = ncomp)

pca_sds_point_default <- fpca.face(Y = dexcom_sds_point, center = T, argvals=(1:tdexcom)/tdexcom, knots = 35, npc = ncomp)

pca_sds_point_default2 <- fpca.face(Y = dexcom_sds_point, center = T, argvals=(1:tdexcom)/tdexcom, knots = 35) #9

pca_sds_point_default3 <- fpca.face(Y = dexcom_sds_point, center = T, argvals=(1:tdexcom)/tdexcom, knots = 5) #6

# Compare fits on all eigenfunctions between FACE (5 knots vs 35 knots) 
############################################################################################################
time_grid_dexcom = c(0, 1:(tdexcom-1) * 5 /60)
ngrid = length(time_grid_dexcom)

fpca_face5 <- data.frame(time = time_grid_dexcom, mu = pca_sds$mu, method = rep("FACE, 5 knots", tdexcom))
fpca_face5$e1 <- pca_sds$efunctions[,1]*sqrt(ngrid)
fpca_face5$e2 <- pca_sds$efunctions[,2]*sqrt(ngrid)
fpca_face5$e3 <- pca_sds$efunctions[,3]*sqrt(ngrid)

fpca_face35 <- data.frame(time = time_grid_dexcom, mu = pca_sds_default$mu, method = rep("FACE, 35 knots", tdexcom))
fpca_face35$e1 <- -pca_sds_default$efunctions[,1]*sqrt(ngrid)
fpca_face35$e2 <- pca_sds_default$efunctions[,2]*sqrt(ngrid)
fpca_face35$e3 <- -pca_sds_default$efunctions[,3]*sqrt(ngrid)

fpca_face5_point <- data.frame(time = time_grid_dexcom, mu = pca_sds_point$mu, method = rep("FACE, 5 knots, point", tdexcom))
fpca_face5_point$e1 <- pca_sds_point$efunctions[,1]*sqrt(ngrid)
fpca_face5_point$e2 <- pca_sds_point$efunctions[,2]*sqrt(ngrid)
fpca_face5_point$e3 <- -pca_sds_point$efunctions[,3]*sqrt(ngrid)

fpca_face35_point <- data.frame(time = time_grid_dexcom, mu = pca_sds_point_default$mu, method = rep("FACE, 35 knots, point", tdexcom))
fpca_face35_point$e1 <- -pca_sds_point_default$efunctions[,1]*sqrt(ngrid)
fpca_face35_point$e2 <- -pca_sds_point_default$efunctions[,2]*sqrt(ngrid)
fpca_face35_point$e3 <- -pca_sds_point_default$efunctions[,3]*sqrt(ngrid)


sd_fPCA_all <- rbind(fpca_face5, fpca_face35, fpca_face5_point, fpca_face35_point)

library(tidyverse)

# Overall glucose variability is slightly different, but not by much.
p1 <- sd_fPCA_all %>% 
  ggplot(aes(x = time, y = mu, col = as.factor(method))) + geom_line() + xlab("Hours of sleep") + ylab("Glucose  [mg / dL]") + ggtitle("Overall glucose variability mean") 
p1


# Helper functions to create facet labels
function_names <- list(
  'e1'="Sd eigenfunction 1",
  'e2'="Sd eigenfunction 2",
  'e3'="Sd eigenfunction 3",
  'e4'="Sd eigenfunction 4"
)

elabeller <- function(variable, value){
  return(function_names[value])
}

figure.path = getwd()
text_size = 16

# Sd eigenfunctions plot
p2 <- sd_fPCA_all %>%
  dplyr::select(-mu)%>%
  gather(key = "type", value = "value", -c(time, method))%>%
  ggplot(aes(x = time, y = value, col = as.factor(method))) + facet_grid(~type, labeller = elabeller) + geom_line()+ xlab("Hours of sleep") + ylab("")+ labs(col = "Method") + theme(text = element_text(size = text_size))
p2


pdf(paste(figure.path, "/SdPCA_methods_comparison.pdf", sep = ""), onefile = F, width = 12, height = 4)
print(p2)
dev.off()

# Conclusion: for the sd eigenfunctions, the shape is so simple that it doesn't matter how many knots, and there if very high agreement between smooth and pointwise


# What if we only use first 6 hours rather than first 7 hours of sleep?
#########################################################################
# Grid length 84 fot 7 hours, only 72 for the first 6 hours
ngrid = 72

# Unweighted mean smoothing as default in the paper
ncomp = 3
pca_means_6h <- fpca.face(Y = raw_means[, 1:ngrid], center = T, argvals=(1:ngrid)/ngrid, knots = 5, npc = ncomp, var = T) # error variance is estimated at 4

# Construct sd estimates
dexcom_sds = matrix(NA, n, ngrid)
rownames(dexcom_sds) = subjectIDs
for (i in 1:n){
  # Extract estimated sleep trajectories for subject i
  raw_dexcom <- data[[i]]
  # Subtract smoothed mean from raw data
  subtract_mean <- raw_dexcom[, 1:ngrid] - matrix(pca_means_6h$Yhat[i,], nrow = nrow(raw_dexcom), ncol = ngrid, byrow = T)
  # Estimate st.dev on the residuals
  dexcom_sds[i, ] <- apply(subtract_mean, 2, function(x) sqrt(sum(x^2, na.rm = T)/(sum(!is.na(x))-1))) # changed this part
}
dexcom_sds[dexcom_sds == Inf] = NA

# Smooth pointwise sds
ncomp = 3
pca_sds_6h <- fpca.face(Y = dexcom_sds, center = T, argvals=(1:ngrid)/ngrid, knots = 5, npc = ncomp) 


# Create new figure of means and components (how different it is from 7 hours?)
time_grid_dexcom = c(0, 1:(ngrid-1) *5 /60)

# Look at all subjects mean
mean_fPCA_all <- data.frame(time = time_grid_dexcom, mu = pca_means_6h$mu)
mean_fPCA_all$e1 <- pca_means_6h$efunctions[,1]*sqrt(ngrid)
mean_fPCA_all$e2 <- -pca_means_6h$efunctions[,2]*sqrt(ngrid)
mean_fPCA_all$e3 <- pca_means_6h$efunctions[,3]*sqrt(ngrid)

text_size = 16

p1 <- mean_fPCA_all %>% 
  ggplot(aes(x = time, y = mu)) + geom_line() + xlab("Hours of sleep") + ylab("Glucose  [mg / dL]") + ggtitle("Overall glucose mean") + theme(text = element_text(size = text_size))
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
  ggplot(aes(x = time, y = value)) + facet_grid(~type, labeller = elabeller) + geom_line()+ xlab("Hours of sleep") + ylab("")+ theme(text = element_text(size = text_size))
p2


# Plot of means and eigenfunctions for sd process
#########################################
# Create sd data frame
sd_fPCA_all <- data.frame(time = time_grid_dexcom, mu = pca_sds_6h$mu)
# Use -eigenfunction for sd 1st component to match orientation for the mean for easer interpretation
sd_fPCA_all$e1 <- pca_sds_6h$efunctions[,1]*sqrt(ngrid)
sd_fPCA_all$e2 <- -pca_sds_6h$efunctions[,2]*sqrt(ngrid)
sd_fPCA_all$e3 <- pca_sds_6h$efunctions[,3]*sqrt(ngrid)

# Overall mean plot for sd process
p1 <- sd_fPCA_all %>% 
  dplyr::select(time, mu) %>% 
  ggplot(aes(x = time, y = mu)) + geom_line() + xlab("Hours of sleep") + ylab("Glucose [mg / dL]") + ggtitle("Overall glucose standard deviation")+ theme(text = element_text(size = text_size))
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
  ggplot(aes(x = time, y = value)) + facet_grid(~type, labeller = elabeller) + geom_line()+ xlab("Hours of sleep") + ylab("")+ theme(text = element_text(size = text_size))
p2

# Conclusion: using 6 hours gives results that are identical to the (first 6 hours) of 7 hours results