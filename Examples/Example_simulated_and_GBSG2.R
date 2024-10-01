#'########################################################
#' September 2024,  code by Diana Shamsutdinova
#'
#' Illustrative example of using the survcompare package
#' a) linear data
#' b) non-linear and interaction terms data
#' c) GBSG2 data https://rdrr.io/cran/pec/man/GBSG2.html
#'    German Breast Cancer Study Group, Schumacher et al. (1994)
#'
#'#######################################################

### Aims of survcompare() method:
### 1) Check if there are non-linear and interaction terms in the data
### 2) Quantify their contribution to the models' performance
### 3) Help researchers to make an informed decision on whether
###     the benefit of using a flexible but less transparent machine learning method is high enough,
###     or the classical (or regularized) Cox model should be preferred.
### Method:
### 1) Using repeated cross-validation, validate the performance of the Cox model (or Cox Lasso)
### 2) Validate the performance of the Survival Random Forest(SRF) (ensembled with the Cox model) which captures non-linearity and interactions automatically.
### 3) Perform statistical test to compare the models
###
### Why ensemble and not just SRF?
### -> The ensemble of Cox and SRF takes the predictions of Cox model and adds to the list of
### predictors to train SRF. This way, we make sure that linearity is captured by SRF at
### least as good as the Cox model, and hence the marginal outperformance can be attributed to
### the qualities of SRF that Cox does not have, i.e. data complexity.
#'##############################################################

#### a) Application to the simulated data with simple linear dependencies ####

# simulate data
mydata_1 <- simulate_linear(200)
predictors <- names(mydata_1)[1:4]
compare_models_1 <- survcompare(mydata_1, predictors, randomseed = 100)

# # Internally validated test performance of CoxPH and SRF_ensemble over 2
# repeated 3 fold cross-validations (inner k = 3 ). Mean performance:
# T C_score  AUCROC Calib_slope  sec
# CoxPH        8.043  0.7320  0.7265      0.8733 0.33
# SRF_ensemble 8.043  0.6988  0.7007      1.0866 8.62
# Diff         0.000 -0.0332 -0.0258      0.2133 8.29
# pvalue         NaN  0.9864  0.9828      0.1238  NaN
#
# Median performance:
#   T C_score  AUCROC Calib_slope  sec
# CoxPH        8.043  0.7295  0.7309      0.7962 0.33
# SRF_ensemble 8.043  0.7169  0.7226      0.9778 8.62
# Diff         0.000 -0.0126 -0.0082      0.1816 8.29
# pvalue         NaN  0.9864  0.9828      0.1238  NaN
#
# SRF_ensemble has NOT outperformed CoxPH with the mean c-index difference of -0.0332.
# The difference is not statistically significant with the p-value = 0.9864.
# The data may NOT contain considerable non-linear or cross-term dependencies
# that could be captured by SRF_ensemble.
# Mean C-score:
#   CoxPH  0.732(95CI=0.721-0.743;SD=0.0163)
# SRF_ensemble 0.6988(95CI=0.6898-0.7079;SD=0.0135)
# Mean AUCROC:
#   CoxPH  0.7265(95CI=0.7119-0.741;SD=0.0217)
# SRF_ensemble 0.7007(95CI=0.6869-0.7146;SD=0.0206)

# NOTE: the same results can be obtained by cross-validating CoxPH and SRF separately,
# and then using survcompare2() function;
# outer_cv, inner_cv, repeat_cv, and randomseed should be the same.
# cv1<- survcox_cv(mydata_1, predictors, randomseed = 100)
# cv2<- survsrfens_cv(mydata_1, predictors, randomseed = 100)
# survcompare2(cv1, cv2)


compare_models_1$main_stats_pooled
#                       mean      sd  95CILow 95CIHigh
# C_score_CoxPH        0.7320 0.01633  0.7210   0.7430
# C_score_SRF_ensemble 0.6898 0.02260  0.6746   0.7050
# AUCROC_CoxPH         0.7265 0.02165  0.7119   0.7410
# AUCROC_SRF_ensemble  0.6919 0.03111  0.6710   0.7128

# More information, including Brier Scores, Calibration slope and alphas can be seen in the
# compare_models_1$results_mean. These are averaged performance metrics on the test sets.
compare_models_1$results_mean
#                 T  C_score   AUCROC      BS BS_scaled Calib_slope Calib_alpha  sec
# CoxPH        8.043  0.73202  0.72649 0.14823   0.12633      0.8733    -0.14200 0.29
# SRF_ensemble 8.043  0.68978  0.69192 0.16037   0.05228      1.0610     0.02838 8.65
# Diff         0.000 -0.04224 -0.03457 0.01214  -0.07405      0.1877     0.17037 8.36
# pvalue         NaN  0.99883  0.98781      NA   0.91152      0.1481     0.01964  NaN


#### b) Application to the simulated data with complex dependencies ####

# simulate data using simsurv_crossterms()
mydata_2 <- simulate_crossterms(200)
mypredictors <- names(mydata_2)[1:4]
compare_models_2 <- survcompare(mydata_2, mypredictors, randomseed = 101)
compare_models_2
# Internally validated test performance of CoxPH and SRF_ensemble over 2 repeated
# 3 fold cross-validations (inner k = 3 ). Mean performance:
#   T C_score AUCROC Calib_slope  sec
# CoxPH        8.315  0.5970 0.5716      0.6427 0.31
# SRF_ensemble 8.315  0.6243 0.6130      0.4244 8.44
# Diff         0.000  0.0273 0.0414     -0.2183 8.13
# pvalue         NaN  0.0523 0.0328      0.8257  NaN
#
# Median performance:
#   T C_score AUCROC Calib_slope  sec
# CoxPH        8.315  0.6101 0.5983      0.7172 0.31
# SRF_ensemble 8.315  0.6159 0.6145      0.3881 8.44
# Diff         0.000  0.0058 0.0162     -0.3291 8.13
# pvalue         NaN  0.0523 0.0328      0.8257  NaN
#
# SRF_ensemble has NOT outperformed CoxPH with the mean c-index difference of 0.0273.
# The difference is not statistically significant with the p-value = 0.0523.
# The data may NOT contain considerable non-linear or cross-term dependencies
# that could be captured by SRF_ensemble.
# Mean C-score:
#   CoxPH  0.597(95CI=0.5691-0.625;SD=0.0416)
# SRF_ensemble 0.6243(95CI=0.6222-0.6264;SD=0.0031)
# Mean AUCROC:
#   CoxPH  0.5716(95CI=0.5339-0.6094;SD=0.0562)
# SRF_ensemble 0.613(95CI=0.6044-0.6216;SD=0.0128)


compare_models_2$main_stats_pooled
#                                 mean      sd 95CILow 95CIHigh
# C_score_CoxPH                  0.5970 0.04160  0.5691   0.6250
# C_score_Survival Random Forest 0.6310 0.01541  0.6207   0.6414
# AUCROC_CoxPH                   0.5858 0.05611  0.5481   0.6235
# AUCROC_Survival Random Forest  0.6540 0.00758  0.6489   0.6591

##################################################################
# c) Application to the GBSG2 data https://rdrr.io/cran/pec/man/GBSG2.html
#  German Breast Cancer Study Group, Schumacher et al. (1994)

library(pec) #for GBSG2 data
data("GBSG2")
# re-format the data - hot-coding binary variables horTh and menostat,
# also assuming that tgrade is ordinary variable 1<2<3
gbsg_data = GBSG2
gbsg_data$horTh = ifelse(gbsg_data$horTh== "no", 0, 1)
gbsg_data$menostat = ifelse(gbsg_data$menostat == "Post", 1, 0)
gbsg_data$tgrade = ifelse(gbsg_data$tgrade == "I", 1, ifelse(gbsg_data$tgrade=="II", 2, 3))
for (i in 1:dim(gbsg_data)[2]) print(c(names(gbsg_data)[i], class(gbsg_data[,i])))
params = c("age", "horTh", "menostat", "tsize", "tgrade", "pnodes", "progrec", "estrec")
gbsg_data$event = gbsg_data$cens
gbsg_data$time = gbsg_data$time/365 #convert into years

# final data for the analyses
gbsg_data = gbsg_data[c("time", "event", params)]
# choose the time horizon for predictions
quantile(gbsg_data[gbsg_data$event==1, "time"],0.95) #4.96
quantile(gbsg_data[gbsg_data$event==0, "time"],0.95) #6.34

compare_models_gbsg <-
  survcompare(gbsg_data,params,fixed_time = 5,repeat_cv = 5, randomseed= 102)

# [1] "Cross-validating CoxPH using 5 repeat(s), 3 outer, 3 inner loops)."
# [1] "Repeated CV 1 / 5"
# |===================================================================| 100%
# [1] "Repeated CV 2 / 5"
# |===================================================================| 100%
# [1] "Repeated CV 3 / 5"
# |===================================================================| 100%
# [1] "Repeated CV 4 / 5"
# |===================================================================| 100%
# [1] "Repeated CV 5 / 5"
# |===================================================================| 100%
# Time difference of 1.347 secs
# [1] "Cross-validating Survival Random Forest using 5 repeat(s), 3 outer, 3 inner loops).
# For SRF inner CV is not used if oob = TRUE (default)"
# [1] "Repeated CV 1 / 5"
# |===================================================================| 100%
# [1] "Repeated CV 2 / 5"
# |===================================================================| 100%
# [1] "Repeated CV 3 / 5"
# |===================================================================| 100%
# [1] "Repeated CV 4 / 5"
# |===================================================================| 100%
# [1] "Repeated CV 5 / 5"
# |===================================================================| 100%
# Time difference of 57.9 secs
#
# Internally validated test performance of CoxPH and Survival Random Forest
# over 5 repeated 3 fold cross-validations (inner k = 3 ). Mean performance:
#                         T C_score  AUCROC Calib_slope   sec
# CoxPH                    5  0.6737  0.7222      0.7779  1.35
# Survival Random Forest   5  0.6898  0.7065      1.2132 57.90
# Diff                     0  0.0161 -0.0157      0.4353 56.55
# pvalue                 NaN  0.0014  0.9709      0.0009   NaN
#
# Median performance:
#                          T C_score  AUCROC Calib_slope   sec
# CoxPH                    5  0.6761  0.7444      0.7792  1.35
# Survival Random Forest   5  0.6909  0.7107      1.1435 57.90
# Diff                     0  0.0148 -0.0336      0.3643 56.55
# pvalue                 NaN  0.0014  0.9709      0.0009   NaN

# Survival Random Forest has outperformed CoxPHby 0.0161 in C-index.
# The difference is statistically significant with the p-value 0.0014**.
# The supplied data may contain non-linear or cross-term dependencies,
# better captured by Survival Random Forest.
# Mean C-score:
#   CoxPH  0.6737(95CI=0.6686-0.6787;SD=0.0043)
# Survival Random Forest 0.6898(95CI=0.6795-0.6969;SD=0.0074)
# Mean AUCROC:
#   CoxPH  0.7222(95CI=0.7046-0.7372;SD=0.0138)
# Survival Random Forest 0.7065(95CI=0.6828-0.7205;SD=0.0158)


round(compare_models_gbsg$main_stats_pooled,5)
#                                  mean      sd 95CILow 95CIHigh
# C_score_CoxPH                  0.6737 0.00434  0.6686   0.6787
# C_score_Survival Random Forest 0.6898 0.00744  0.6794   0.6969
# AUCROC_CoxPH                   0.7222 0.01379  0.7046   0.7371
# AUCROC_Survival Random Forest  0.7065 0.01580  0.6828   0.7205
