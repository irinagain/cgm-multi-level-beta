# R code for creating figures for the paper on "Modeling CGM trajectories during sleep"
# Creates figures with scores and calculates association with A1C
########################################################################################
# Author: Irina Gaynanova

# ATTENTION: this code relies on actual A1C values for real CGM data, so can only be run with real data


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


# Specify path to save figures
figure.path = getwd()
text_size = 20


# Load real CGM data (how the code was used in the paper)
load("CGMdataReal.Rdata")
data = CGMdataReal
minval = minvalReal
maxval = maxvalReal
subjectIDs = idsReal

# Estimate marginal parameters for multi-level functional Beta model, use mean/sd for normal data generation; alpha/beta/m/M for beta data generation
source("BetaEstimatingFunctions.R")
out_mfbeta <- estimateBetaParameters(data, minval, maxval)

# Extract all the scores
scores_means = out_mfbeta$pca_means$scores
scores_sds = out_mfbeta$pca_sds$scores

# Combine all the scores and the subject
# Use -scores for sd 1st component to match orientation for the mean for easer interpretation
dexcom = data.frame(id = subjectIDs, scores_sd1 = -scores_sds[,1], scores_sd2 = scores_sds[,2], scores_sd3 = scores_sds[,3], scores1 = scores_means[,1],scores2 = scores_means[,2], scores3 = scores_means[,3])


########################################################
# Load the A1C information (NOT RUN - requires real A1C data)
########################################################
library(readxl)
newa1c.path = '../../CGMS/UpdatedA1C-2018-04-06/Punjabi_Hypnos_2018_01_23.xlsx'
newa1c = read_excel(newa1c.path, sheet = "results")


###########################################################################################
# Create a figure of first scores versus A1C plus calculate correlation and R2 values for linear model
####################################################################################

# Match a1c by subject - do not have a1c values for all subjects
a1c_match <- match(dexcom$id, newa1c$`sample ID`)
a1c = newa1c$`Glycated Hemoglobin, HbA1c (%)`[a1c_match]
dexcom$a1c = a1c

# NOT RUN: selected ID numbers for 6 selected subjects from real CGM data only
selectedID <- c(70058, 70134, 70139, 70142, 70153, 70186)
dexcom_selected <- dexcom[dexcom$id %in% selectedID,]
dexcom_selected$id <- as.factor(as.character(dexcom_selected$id))

# Scores with color for each subject
pdf(paste(figure.path, "/FirstScoresVsA1c.pdf", sep = ""), onefile = TRUE, width = 16, height = 4)
p1 = dexcom%>%ggplot(aes(x = scores1, y = a1c)) + geom_point(alpha = 0.3, size = 2) + xlab("Mean score 1") + ylab(expression(paste("HbA"["1c"]," %"))) + geom_point(aes(scores1, a1c, col = id), data = dexcom_selected, size = 3)+ theme(legend.position="none")+ theme(text = element_text(size = text_size))
p2 = dexcom%>%ggplot(aes(x = scores_sd1, y = a1c)) + geom_point(alpha = 0.3, size = 2) + xlab("Sd score 1") + ylab(expression(paste("HbA"["1c"]," %"))) + geom_point(aes(scores_sd1, a1c, col = id), data = dexcom_selected, size = 3)+ theme(legend.position="none")+ theme(text = element_text(size = text_size))
multiplot(p1,p2, cols = 2)
dev.off()

# Calculate correlations between scores and a1c
cor(dexcom$a1c, dexcom$scores1, use = "complete.obs") # 0.79
cor(dexcom$a1c, dexcom$scores_sd1, use = "complete.obs") # 0.60

# Calculate correlations with all other scores and a1c -> all below 0.19
cor(dexcom$a1c, dexcom$scores2, use = "complete.obs") # -0.19
cor(dexcom$a1c, dexcom$scores3, use = "complete.obs") # 0.10
cor(dexcom$a1c, dexcom$scores_sd2, use = "complete.obs") #-0.18
cor(dexcom$a1c, dexcom$scores_sd3, use = "complete.obs") # 0.02


# Calculate R2 values from predicting A1C using linear model on scores
# (1) Only the first scores from mean and sd
lm1 <- lm(a1c ~ scores1 + scores_sd1, data = dexcom)
summary(lm1)$r.squared # 0.64
# (2) All three scores from mean and sd
lm2 <- lm(a1c ~ scores1 + scores2+ scores3 + scores_sd1 + scores_sd2 + scores_sd3, data = dexcom)
summary(lm2)$r.squared # 0.70
summary(lm2) # first 2 scores significant for mean, first score significant for st.dev.

cor(dexcom$scores1, dexcom$scores_sd1) # 0.63
cor(dexcom$scores2, dexcom$scores_sd2) # 0.46