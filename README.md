# R package survcompare
internally validates and compares performance of the Cox Proportionate Hazards and Survival Random Forest 

## How to use this package

You can install the package from its github directory by running '''devtools::install_github("dianashams/ensemblesurv")''' command. 
The files in the "Example/" illustrate application of the main survcompare::survcompare() function to the simulated and GBSG2 data sets (https://rdrr.io/cran/pec/man/GBSG2.html), and compute internally-validated performance metrics (using repeated nested cross-validation) along with the results of the statistical testing of whether Survival Random Forest outperforms the Cox Proportionate Hazard model, or Cox Lasso model depending on the user's input.  

#Example:
mydata_2 <- simsurv_crossterms(500)
compare_models_2 <- survcompare(mydata, names(mydata)[1:4], predict_t = 10,
                              outer_cv = 3, inner_cv = 3, repeat_cv = 5, 
                              useCoxLasso = FALSE)

/' [1] "Cross-validating Cox-PH ( 5 repeat(s), 3 loops)"
/' # |=============================================================| 100%
/'# [1] "Cross-validating Survival Random Forest - Cox model ensemble ( 5 repeat(s), 3 outer, 3 inner loops)"
/'# |=============================================================| 100%
/'# Time difference of 21.4098 secs
/'# 
/'# Internally validated test performance of CoxPH     and Survival Random Forest ensemble:
/'#                T C_score  AUCROC     BS BS_scaled Calib_slope
/'# CoxPH         10  0.7392  0.7822 0.0978    0.3514      0.9453
/'# SRF_Ensemble  10  0.7217  0.7640 0.1027    0.3199      0.9806
/'# Diff           0 -0.0175 -0.0182 0.0049   -0.0316      0.0353
/'# pvalue       NaN  0.9571  0.9612 0.7068    0.9915      0.3044
/'#              Calib_alpha   sec
/'# CoxPH             0.2229  1.22
/'# SRF_Ensemble      0.2153 21.41
/'# Diff             -0.0077 20.19
/'# pvalue            0.6798   NaN
/'# 
/'# Survival Random Forest ensemble has NOT outperformed CoxPH   with mean c-index difference of-0.0175.
/'# The difference is not statistically significant, p-value = 0.9571. The data may NOT contain considerable non-linear or cross-term dependencies, better /'captured by the Survival Random Forest.
/'# C-score: 
/'#   CoxPH      0.7392(95CI=0.6958-0.7933;SD=0.0338)
/'# SRF_Ensemble 0.7217(95CI=0.6368-0.7833;SD=0.0446)
/'# AUCROC:
/'#   CoxPH      0.7822(95CI=0.7349-0.8316;SD=0.0334)
/'# SRF_Ensemble 0.764(95CI=0.6886-0.8281;SD=0.0433)


**All package functions assume the data is a data frame,  "time" and "event" columns define survival outcome; most functions require a list of column names that correspond to model's predictors**

## If you use this project's code, please cite:

Shamsutdinova, D., Stamate, D., Roberts, A., & Stahl, D. (2022). Combining Cox Model and Tree-Based Algorithms to Boost Performance and Preserve Interpretability for Health Outcomes. In IFIP International Conference on Artificial Intelligence Applications and Innovations (pp. 170-181). Springer, Cham.

#### Support or Contact
diana.shamsutdinova@kcl.ac.uk
