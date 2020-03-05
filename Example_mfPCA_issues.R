# Fit the multi-level FPCA model to CGM data and check level 1 fits
##############################################################
# Author: Irina Gaynanova

# Load real CGM data (how the code was used in the paper)
load("CGMdataReal.Rdata")
data = CGMdataReal
minval = minvalReal
maxval = maxvalReal
subjectIDs = idsReal

## ATTENTION: the issue with level 1 eigenfunctions could not be directly reproduced on synthetic data right now, but is consistently reproduced on real data using refund version 0.1-21 (in 0.1-20 the issue was even more pronounced)

# In absence of real CGM data, the code below will work with provided synthetic data (except for cases where real subject IDs are used).
# load("CGMdataSynthetic.Rdata")
# data = simul_beta_real
# minval = minvalSythetic
# maxval = maxvalSynthetic
# subjectIDs = idsSynthetic


figure.path = getwd()

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

maxH = 7 # Maximal hour
tdexcom = round(maxH * 60 / 5) # Length of time grid in 5 minute intervals from the estimated sleep onset
n = length(minval) # Number of subjects

# Supporting libraries
library(tidyverse)
library(refund) # version 0.1-21, 2019-12-06 from CRAN 

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
# 2 components at level 1, 7 components at level 2

# Calculate posterior means and st.dev. on subject level and overall for scores
scores2 <- data.frame(scores = out_mfpca$scores$level2, subjects = id)

scoresSummaries <- scores2%>%
  group_by(subjects)%>%
  summarize(m1 = mean(scores.1), m2 = mean(scores.2), m3 = mean(scores.3), m4 = mean(scores.4), m5 = mean(scores.5), m6 = mean(scores.6), sd1 = sd(scores.1), sd2 = sd(scores.2), sd3 = sd(scores.3), sd4 = sd(scores.4), sd5 = sd(scores.5), sd6 = sd(scores.6))

### ATTENTION: Below uses real subjects' ids, will not run on synthetic data

subject = 70324
data <- data.frame(Y = c(Y[id == subject, ]), Yhat = c(out_mfpca$Yhat[id == subject,]), day = rep(1:nrow(Y[id == subject, ]), ncol(Y)), time = rep(0:(ncol(Y)-1), each = sum(id == subject)) * 5/60, subjectY = c(out_mfpca$Yhat.subject[id == subject,]))
p1 = ggplot(data, aes(x = time, y = Y, group = as.factor(day))) + geom_line() + geom_line(aes(x = time, y = subjectY, group = as.factor(day)), data = data, col = 'red', size = 2)+ ylim(c(50,400))+ggtitle(subject) + ylab("Glucose [mg/dL]") + xlab("Hours of sleep")
print(p1)


subject = 90021
data <- data.frame(Y = c(Y[id == subject, ]), Yhat = c(out_mfpca$Yhat[id == subject,]), day = rep(1:nrow(Y[id == subject, ]), ncol(Y)), time = rep(0:(ncol(Y)-1), each = sum(id == subject)) * 5/60, subjectY = c(out_mfpca$Yhat.subject[id == subject,]))
p2 = ggplot(data, aes(x = time, y = Y, group = as.factor(day))) + geom_line() + geom_line(aes(x = time, y = subjectY, group = as.factor(day)), data = data, col = 'red', size = 2)+ ylim(c(50,400))+ggtitle(subject) + ylab("Glucose [mg/dL]") + xlab("Hours of sleep")
print(p2)

subject = 70227
data <- data.frame(Y = c(Y[id == subject, ]), Yhat = c(out_mfpca$Yhat[id == subject,]), day = rep(1:nrow(Y[id == subject, ]), ncol(Y)), time = rep(0:(ncol(Y)-1), each = sum(id == subject)) * 5/60, subjectY = c(out_mfpca$Yhat.subject[id == subject,]))
p3 = ggplot(data, aes(x = time, y = Y, group = as.factor(day))) + geom_line() + geom_line(aes(x = time, y = subjectY, group = as.factor(day)), data = data, col = 'red', size = 2)+ ylim(c(50,400))+ggtitle(subject) + ylab("Glucose [mg/dL]") + xlab("Hours of sleep")
print(p3)

pdf(file = paste(figure.path, "/MFPCAfail_examples.pdf", sep=""), width = 10, height = 4)
multiplot(p1, p2, p3, cols = 3)
dev.off()