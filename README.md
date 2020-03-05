# Analysis of CGM night trajectories using multi-level functional Beta model

Synthetic continuous glucose monitoring (CGM) data and corresponding R code to fit multi-level functional Beta model.

  * **BetaEstimatingFunctions.R** - functions for fitting multi-level functional Beta model
  
  * **CGMdataSynthetic.Rdata** - synthetic CGM data. Contains list of n=124 subjects with 14 night trajectories per subject over 5 min grid in 7 hour sleep period. Also contains minimal and maximal glucose values for each subject. The data have been generated based on real CGM data described in the paper following **SimulateCGMdata.R** script
  
  * **fPCAStability.R** - code for visual display of functional PCA results for mean and standard deviation processes, as well as for evaluating stability of corresponding eigenfunctions via subsampling
  
  * **Code_for_figures.R** - code for creating figures of subject-specific night trajectories, corresponding quantile bands, and minimal/maximal glucose values
  
  * **Code_for_A1C_analyses.R** - code for the analysis of association between fPCA scores and A1C as described in the paper together with code for creating corresponding figures
  
  * **ComparePercentiles.R** - code for comparing percentiles of functional Normal and Beta models based on simulated data
  
  * **QQplots.R** - code for creating QQplots at a global, time-specific and subject-specific levels
  
  * **NaiveBands.R** - code for estimation of quantiles based on separate subject/time approaches (parametric Beta and nonparametric kernel density estimators)

  * **Beta_model_fit_perturb.R** - comparison of estimated quantiles between Beta model with estimated values of min/max, and Beta model with perturbed values of min/max
  
  * **Example_mfPCA_issues.R** - code for fitting multi-level fPCA model on CGM data with code for figure that highlight issues with level 1 fits
  
  * **FPCA_methods_comparison.R** - code for comparing different FPCA methods in terms of estimated mean and eigenfunctions for subjects' means and standard deviations