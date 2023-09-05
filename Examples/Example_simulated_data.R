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
# [1] "Cross-validating Survival Random Forest - Cox model ensemble ( 5 repeat(s), 3 outer, 3 inner loops)"
# |=============================================================| 100%
# Time difference of 21.4098 secs
# 
# Internally validated test performance of CoxPH     and Survival Random Forest ensemble:
#                T C_score  AUCROC     BS BS_scaled Calib_slope
# CoxPH         10  0.7392  0.7822 0.0978    0.3514      0.9453
# SRF_Ensemble  10  0.7217  0.7640 0.1027    0.3199      0.9806
# Diff           0 -0.0175 -0.0182 0.0049   -0.0316      0.0353
# pvalue       NaN  0.9571  0.9612 0.7068    0.9915      0.3044
#              Calib_alpha   sec
# CoxPH             0.2229  1.22
# SRF_Ensemble      0.2153 21.41
# Diff             -0.0077 20.19
# pvalue            0.6798   NaN
# 
# Survival Random Forest ensemble has NOT outperformed CoxPH   with mean c-index difference of-0.0175.
# The difference is not statistically significant, p-value = 0.9571. The data may NOT contain considerable non-linear or cross-term dependencies, better captured by the Survival Random Forest.
# C-score: 
#   CoxPH      0.7392(95CI=0.6958-0.7933;SD=0.0338)
# SRF_Ensemble 0.7217(95CI=0.6368-0.7833;SD=0.0446)
# AUCROC:
#   CoxPH      0.7822(95CI=0.7349-0.8316;SD=0.0334)
# SRF_Ensemble 0.764(95CI=0.6886-0.8281;SD=0.0433)

compare_models_2$main_stats
#                         mean         sd   95CILow  95CIHigh
# C_score_CoxPH        0.7391683 0.03378723 0.6958052 0.7933356
# C_score_SRF_Ensemble 0.7216899 0.04459629 0.6368015 0.7833033
# AUCROC_CoxPH         0.7822150 0.03336371 0.7349474 0.8316363
# AUCROC_SRF_Ensemble  0.7640463 0.04331564 0.6885705 0.8280955


##################
library(pec) #for GBSG2 data
data("GBSG2")
#names(gbsg_data)
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

compare_models_gbsg <- survcompare(gbsg_data, params, predict_t = final_time,
                                outer_cv = 3, inner_cv = 3, repeat_cv = 5, 
                                useCoxLasso = FALSE)

# [1] "Cross-validating Cox-PH ( 5 repeat(s), 3 loops)"
# |=============================================================| 100%
# [1] "Cross-validating Survival Random Forest - Cox model ensemble ( 5 repeat(s), 3 outer, 3 inner loops)"
# |=============================================================| 100%
# Time difference of 56.00061 secs
# 
# Internally validated test performance of CoxPH     and Survival Random Forest ensemble:
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

