#########################################################
#' September 2023,  code by Diana Shamsutdinova
#' 
#' Illustrative example of using the survcompare package
#' a) linear data
#' b) non-linear and interaction terms data
#' c) GBSG2 data https://rdrr.io/cran/pec/man/GBSG2.html  
#'    German Breast Cancer Study Group, Schumacher et al. (1994) 
#' 
########################################################

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


###############################################################

# a) Application to the simulated data with simple linear dependencies

# simulate data 
mydata_1 <- simulate_linear(500)
predictors <- names(mydata_1)[1:4]
compare_models_1 <- survcompare(mydata_1, predictors, 10,
                              outer_cv = 3, inner_cv = 3, repeat_cv = 1, 
                              useCoxLasso = TRUE)

summary(compare_models_1)
# Internally validated test performance of CoxLasso   and Survival Random Forest ensemble:
#                 T C_score  AUCROC     BS BS_scaled Calib_slope
# CoxLasso      10  0.7684  0.7562 0.1133    0.3045      1.1366
# SRF_Ensemble  10  0.7657  0.7489 0.1186    0.2726      1.4142
# Diff           0 -0.0027 -0.0073 0.0052   -0.0319      0.2775
# pvalue       NaN  0.5334  0.6005 0.6239    0.9887      0.2409
#               Calib_alpha  sec
# CoxLasso          0.2148 0.34
# SRF_Ensemble      0.1627 3.03
# Diff             -0.0521 2.69
# pvalue            0.7563  NaN
# 
# Survival Random Forest ensemble has NOT outperformed CoxLasso  with mean c-index difference of-0.0027.
# The difference is not statistically significant, p-value = 0.5334. The data may NOT contain considerable non-linear or cross-term dependencies, better captured by the Survival Random Forest.
# C-score: 
#   CoxLasso    0.7684(95CI=0.7075-0.823;SD=0.0611)
#   SRF_Ensemble 0.7657(95CI=0.7229-0.8217;SD=0.0537)
# AUCROC:
#  CoxLasso   0.7562(95CI=0.6938-0.8617;SD=0.0987)
#  SRF_Ensemble 0.7489(95CI=0.6508-0.8798;SD=0.1251)


compare_models_1$main_stats
#                          mean         sd   95CILow  95CIHigh
# CoxPH_______C_score 0.7387869 0.04301220 0.6700542 0.8055633
# SRFensemble_C_score 0.7135155 0.04300221 0.6417006 0.7851860
# CoxPH________AUCROC 0.7842373 0.04593391 0.7053326 0.8534232
# SRFensemble__AUCROC 0.7555365 0.04596786 0.6655202 0.8234735


# b) Application to the simulated data with complex dependencies

# simulate data using simsurv_crossterms()
mydata_2 <- simulate_crossterms(500)
mypredictors <- names(mydata_2)[1:4]
compare_models_2 <- survcompare(mydata_2, mypredictors, 10,
                              outer_cv = 3, inner_cv = 3, repeat_cv = 5, 
                              useCoxLasso = FALSE)

# [1] "Cross-validating Cox-PH ( 5 repeat(s), 3 loops)"
# |===============================================================| 100%
# [1] "Cross-validating Survival Random Forest - Cox model ensemble ( 5 repeat(s), 3 outer, 3 inner loops)"
# |===============================================================| 100%
# Time difference of 22.68073 secs
# 
# Internally validated test performance of CoxPH     and Survival Random Forest ensemble:
#                T C_score AUCROC      BS BS_scaled Calib_slope Calib_alpha   sec
# CoxPH         10  0.6484 0.6502  0.1306    0.1714      0.7985      0.2242  1.18
# SRF_Ensemble  10  0.7479 0.7601  0.1061    0.3282      0.7336      0.2361 22.68
# Diff           0  0.0995 0.1098 -0.0245    0.1568     -0.0648      0.0119 21.50
# pvalue       NaN  0.0000 0.0000  0.0104    0.0000      0.7684      0.3294   NaN
# 
# Survival Random Forest ensemble has outperformed CoxPH    by 0.0995 in C-index.
# The difference is statistically significant with the p-value = 0***.
# The supplied data may contain non-linear or cross-term dependencies, 
# better captured by the Survival Random Forest.
# C-score: 
#   CoxPH      0.6484(95CI=0.5491-0.7104;SD=0.0476)
# SRF_Ensemble 0.7479(95CI=0.6933-0.8316;SD=0.0433)
# AUCROC:
#   CoxPH      0.6502(95CI=0.5428-0.7078;SD=0.0484)
# SRF_Ensemble 0.7601(95CI=0.6986-0.8381;SD=0.045)
                    

compare_models_2$main_stats
#                           mean         sd   95CILow  95CIHigh
# C_score_CoxPH        0.6483860 0.04764071 0.5491227 0.7103612
# C_score_SRF_Ensemble 0.7478547 0.04332483 0.6933300 0.8315516
# AUCROC_CoxPH         0.6502224 0.04842263 0.5428135 0.7077804
# AUCROC_SRF_Ensemble  0.7600702 0.04495999 0.6985573 0.8381287


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
final_time = 5

compare_models_gbsg <-
  survcompare(
    gbsg_data,
    params,
    predict_t = final_time,
    outer_cv = 3,
    inner_cv = 3,
    repeat_cv = 5,
    useCoxLasso = FALSE
  )

# [1] "Cross-validating Cox-PH ( 5 repeat(s), 3 loops)"
# |=============================================================| 100%
# [1] "Cross-validating Survival Random Forest - Cox model ensemble ( 5 repeat(s), 3 outer, 3 inner loops)"
# |=============================================================| 100%
# Time difference of 56.00061 secs
# 
# Internally validated test performance of CoxPH and Survival Random Forest ensemble:
#                T C_score  AUCROC      BS BS_scaled Calib_slope
# CoxPH          5  0.6773  0.7266  0.2435    0.0631      1.1869
# SRF_Ensemble   5  0.6938  0.7230  0.2188    0.1584      1.3555
# Diff           0  0.0165 -0.0036 -0.0247    0.0953      0.1686
# pvalue       NaN  0.0001  0.7043  0.0819    0.0236      0.0596
# Calib_alpha   sec
# CoxPH             0.2393  1.52
# SRF_Ensemble      0.8103 56.00
# Diff              0.5709 54.48
# pvalue            0.0067   NaN
# 
# Survival Random Forest ensemble has outperformed CoxPH by 0.0165*** in C-index.
# The difference is statistically significant, p-value = 9.1e-05. 
# The supplied data may contain non-linear or cross-term dependencies, 
# better captured by the Survival Random Forest.
# C-score: 
#   CoxPH      0.6773(95CI=0.6392-0.7207;SD=0.0254)
# SRF_Ensemble 0.6938(95CI=0.6617-0.7475;SD=0.0282)
# AUCROC:
#   CoxPH      0.7266(95CI=0.6529-0.8073;SD=0.0501)
# SRF_Ensemble 0.723(95CI=0.6599-0.8206;SD=0.0515)

round(compare_models_gbsg$main_stats,5)
#                         mean      sd 95CILow 95CIHigh
# C_score_CoxPH        0.67727 0.02538 0.63917  0.72065
# C_score_SRF_Ensemble 0.69381 0.02821 0.66174  0.74751
# AUCROC_CoxPH         0.72660 0.05008 0.65287  0.80733
# AUCROC_SRF_Ensemble  0.72299 0.05149 0.65988  0.82055

