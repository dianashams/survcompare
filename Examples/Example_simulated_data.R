#########################################################
#' September 2023,  code by Diana Shamsutdinova
#' 
#' Illustrative example of using the survcompare package
#' 
########################################################

### Aims: 
### 1) Check if there are non-linear and interaction terms in the data
### 2) Quantify their contribution to the models' performance
### Methods:
### 1) Compare the 

# 1) simulate observations with survival outcome depending 
# linearly on log-hazards using package's function simsurv_linear()

mydata_1 <- simsurv_linear(200)
predictors <- names(mydata_1)[1:4]
compare_models_1 <- survcompare(mydata_1, predictors, predict_t = 10,
                              outer_cv = 3, inner_cv = 3, repeat_cv = 1, 
                              useCoxLasso = TRUE)

summary.survcompare(compare_models_1, 1)
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


# 2) simulate 500 observations with survival outcome depending 
# on interaction and non-linear terms (using package's function simsurv_crossterms)

mydata_2 <- simsurv_crossterms(500)
compare_models_2 <- survcompare(mydata, names(mydata)[1:4], predict_t = 10,
                              outer_cv = 3, inner_cv = 3, repeat_cv = 5, 
                              useCoxLasso = FALSE)


# [1] "Cross-validating Cox-PH ( 5 repeat(s), 3 loops)"
# |=============================================================| 100%
# [1] "Cross-validating Survival Random Forest - Cox PH Ensemble ( 5 repeat(s), 3 outer, 3 inner loops)"
# |=============================================================| 100%
# Time difference of 23.88578 secs
# 
# Internally validated test performance of Cox-PH, Cox-Lasso and Ensemble 1 (Cox-Survival Random Forest):
#               T  C_score   AUCROC        BS BS_scaled Calib_slope
# CoxPH         10 0.628820 0.620980  0.131362  0.153630    0.809373
# SRF_ensemble  10 0.724878 0.730674  0.112147  0.276923    0.823634
# Diff           0 0.096058 0.109694 -0.019215  0.123293    0.014261
# pvalue       NaN 0.000018 0.000001  0.039339  0.000002    0.452790
# Calib_alpha   sec
# CoxPH           0.202257  1.40
# SRF_ensemble    0.201479 23.89
# Diff           -0.000777 22.49
# pvalue          0.508386   NaN
# 
# Survival Random Forest ensemble has outperformed Cox-PH by 0.0961 in C-index.
# The difference is statistically significant with p-value = 1.8e-05***.
# The supplied data may contain non-linear or cross-term dependencies 
# better captured by survival random forest.
# C-score: 
#   CoxPH  0.6288(95CI=0.5219-0.7209;SD=0.0615)
# Ensemble1 0.7249(95CI=0.6881-0.7814;SD=0.0311)
# AUCROC:
#   CoxPH  0.621(95CI=0.501-0.7263;SD=0.0656)
# Ensemble1 0.7307(95CI=0.6705-0.7896;SD=0.038)
#  

compare_models_2$main_stats
#                         mean         sd   95CILow  95CIHigh
# CoxPH_______C_score 0.6288199 0.06148075 0.5219370 0.7208713
# SRFensemble_C_score 0.7248778 0.03114461 0.6881374 0.7813728
# CoxPH________AUCROC 0.6209804 0.06564835 0.5010091 0.7263161
# SRFensemble__AUCROC 0.7306740 0.03795796 0.6705293 0.7896202


##################
library(pec) #for GBSG2 data
data("GBSG2")
names(gbsg_data)
dim(GBSG2)
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
quantile(gbsg_data[gbsg_data$event==1, "time"],0.95) #4.96
quantile(gbsg_data[gbsg_data$event==0, "time"],0.95) #6.34
final_time = 5