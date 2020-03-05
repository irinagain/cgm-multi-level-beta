# Author: Irina Gaynanova

###########################################################################
# Supporting function for fitting Beta distribution using method of moments
###########################################################################
# Apply method of moments to find alpha and beta for 4 parameter beta distribution using 
# xbar - mean
# sd - standard deviation
# [lower, upper] - support of the distribution
getAlphaBeta <- function(xbar, sd, lower, upper){
  if (is.na(xbar)|is.na(sd)){
    return(list(alpha = NA, beta = NA))
  }
  ybar = (xbar - lower)/(upper - lower)
  yvar = sd^2/(upper-lower)^2
  if (yvar < ybar*(1-ybar)){
    tmp = (ybar*(1-ybar)/yvar - 1)
    alpha = ybar*tmp
    beta = (1-ybar)*tmp
    return(list(alpha = alpha, beta = beta))
  }else{
    return(list(alpha = NA, beta = NA))
  }
}

###########################################################################
# Estimate parameters of multi-level functional Beta model
###########################################################################
# train - list of length n (number of subjects), each elements is a nj (number of sleep periods) by tdexcom (number of points in 7 hour sleep period at which data are collected) matrix of sleep trajectories for corresponding subject
# minval - length n vector with values of minimal m_i, minimal observed glucose value for subject i
# maxval - length n vector with values of maximal M_i, maximal observed glucose value for subejct i
# npc_mean - number of components to use for mean FPCA, default 3
# npc_sd - number of components to use for sd FPCA, default 3
estimateBetaParameters <- function(train, minval, maxval, npc_mean = 3, npc_sd = 3){
  # Package refund is used for functional PCA
  require(refund)
  n = length(train)
  maxH = 7 # Set the maximal hour
  tdexcom = round(maxH * 60 / 5)
  
  # Calculate pointwise means
  ##############################################
  raw_means = matrix(NA, n, tdexcom)
  for (i in 1:n){
    # Extract estimated sleep trajectories for subject i
    raw_dexcom <- train[[i]][, 1:tdexcom]
    # Calculate pointwise means
    raw_means[i,] <- colMeans(raw_dexcom, na.rm = T)
  }

  # Smooth pointwise means
  ##############################################
  ncomp = npc_mean
  pca_means <- fpca.face(Y = raw_means, center = T, argvals=(1:tdexcom)/tdexcom, knots=5, npc = ncomp)
  
  # Calculate pointwise sds adjusting for smoothed mean
  ########################################################
  dexcom_sds = matrix(NA, n, tdexcom)
  for (i in 1:n){
    # Extract estimated sleep trajectories for subject i
    raw_dexcom <- train[[i]][, 1:tdexcom]
    # Subtract smoothed mean from raw data
    subtract_mean <- raw_dexcom - matrix(pca_means$Yhat[i,], nrow = nrow(raw_dexcom), ncol = ncol(raw_dexcom), byrow = T)
    # Estimate st.dev on the residuals
    dexcom_sds[i, ] <- apply(subtract_mean, 2, function(x) sqrt(sum(x^2, na.rm = T)/(sum(!is.na(x))-1)))
  }
  dexcom_sds[dexcom_sds == Inf] = NA # Adjust for possibility of 1 curve for some time points which would put zero in the denominator; should really be NA there, gets corrected after smoothing
  
  # Smooth pointwise sds
  ########################################################
  ncomp = npc_sd
  pca_sds <- fpca.face(Y = dexcom_sds, center = T, argvals=(1:tdexcom)/tdexcom, knots=5, npc = ncomp)
  
  # Perform correction for possibly negative sds due to smoothing - doesn't happen on our data but just in case
  ########################################################
  indneg = which(pca_sds$Yhat < 0, arr.ind = T)
  if (length(indneg) > 0){
    pca_sds$Yhat[indneg] = min(pca_sds$Yhat[pca_sds$Yhat>0])
  }
  
  # Calculate moment estimates for alpha and beta
  ######################################################################################
  alphaM <- matrix(0, n, tdexcom)
  betaM <- matrix(0, n, tdexcom)
  for (i in 1:n){
    # Find alpha_i(t) and beta_i(t) using method of moments
    for (t in 1:tdexcom){
      out <- getAlphaBeta(pca_means$Yhat[i,t], pca_sds$Yhat[i,t], lower = minval[i], upper = maxval[i])
      alphaM[i,t] = out$alpha
      betaM[i,t] = out$beta
    }
  }
  # Return output
  ######################################################################################
  return(list(alphaM = alphaM, betaM = betaM, means_smooth = pca_means$Yhat, sd_smooth = pca_sds$Yhat, pca_means = pca_means, pca_sds = pca_sds))
}



