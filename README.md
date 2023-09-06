# R package "survcompare": 
### Internal validation and comparison of the Cox Proportionate Hazards model and Survival Random Forest performances using repeated cross-validation.

#### The package aims to help researchers to make an informed decision on whether the benefit of using a flexible but less transparent machine learning method is high enough, or the classical (or regularized) Cox model should be preferred.

### Method: 
The package checks whether there are considerable non-linear and interaction terms in the data, and quantifies their contributions to the models' performance. Using repeated nested cross-validation, the package
  * validates Cox Proportionate Hazards model (or Cox Lasso depending on the user's input)
  * validates Survival Random Forest(SRF) ensembled with the baseline Cox model
  * performs statistical testing of whether the Survival Random Forest ensemble outperforms the Cox model

The performance metrics include
 * discrimination measures: Harrell's concordance index, time-dependent AUCROC,
 * calibration measures: calibration slope, calibration alpha
 * overall fit: Brier score, Scaled Brier score 

#### FAQ: Why the ensemble and not just SRF? 
The ensemble of Cox and SRF takes the predictions of the Cox model and adds to the list of predictors to train SRF. This way, we make sure that linearity is captured by SRF at  least as good as the Cox model, and hence the marginal outperformance can be attributed to the qualities of SRF that Cox does not have, i.e. data complexity.

which captures non-linearity and interactions automatically is validated.

### How to use the package 
You can install the package from its github directory by running the `devtools::install_github("dianashams/survcompare")` command. 

**The main function is `survcompare(data, predictors)`. The data should be a data frame, with "time" and "event" columns defining the survival outcome. A list of column names corresponding to the predictors to be used should also be supplied.**

The files in the "Example/" folder illustrate `survcompare`'s  application to the simulated and GBSG2  (https://rdrr.io/cran/pec/man/GBSG2.html) datasets. The outputs contain  internally-validated performance metrics along with the results of the statistical testing of whether Survival Random Forest outperforms the Cox Proportionate Hazard model (or Cox Lasso).  

### Example:

```R
mydata <- simsurv_crossterms(500)
mypredictors <- names(mydata)[1:4]
compare_models <- survcompare(mydata, mypredictors, predict_t = 10, inner_cv = 3, repeat_cv = 5)

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

compare_models$main_stats

#                           mean         sd   95CILow  95CIHigh
# C_score_CoxPH        0.6483860 0.04764071 0.5491227 0.7103612
# C_score_SRF_Ensemble 0.7478547 0.04332483 0.6933300 0.8315516
# AUCROC_CoxPH         0.6502224 0.04842263 0.5428135 0.7077804
# AUCROC_SRF_Ensemble  0.7600702 0.04495999 0.6985573 0.8381287
```

### If you use this project's code, please cite:

Shamsutdinova, D., Stamate, D., Roberts, A., & Stahl, D. (2022). Combining Cox Model and Tree-Based Algorithms to Boost Performance and Preserve Interpretability for Health Outcomes. In IFIP International Conference on Artificial Intelligence Applications and Innovations (pp. 170-181). Springer, Cham.

### Support or Contact
diana.shamsutdinova@kcl.ac.uk

### Links and references: 
[1] Ishwaran, H., Lauer, M.S., Blackstone, E.H., Lu, M.: randomForestSRC: Random Survival Forests Vignette (2021)

[2] Cox, D. R. (1972). Regression models and life‐tables. Journal of the Royal Statistical Society: Series B (Methodological), 34(2), 187-202.

[3] Shamsutdinova, D., Stamate, D., Roberts, A., & Stahl, D. (2022). Combining Cox Model and Tree-Based Algorithms to Boost Performance and Preserve Interpretability for Health Outcomes. In IFIP International Conference on Artificial Intelligence Applications and Innovations (pp. 170-181). Springer, Cham.

[4] Amunategui, M.: Data Exploration & Machine Learning, Hands-on. https://amunategui.github.io/survival-ensembles/index.html

