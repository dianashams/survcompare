# R package "survcompare": 

### Internal validation and comparison of the Cox Proportionate Hazards model (CoxPH) and Survival Random Forest (SRF) performances using repeated cross-validation.

#### The package aims to help researchers to make an informed decision on whether the benefit of using a flexible but less transparent machine learning method is high enough, or the classical (or regularized) Cox model should be preferred.

![image](https://github.com/dianashams/ensemble-methods-for-survival-analysis/blob/gh-pages/survcompare_cartoon.png)

### Method: 
The package checks whether there are considerable non-linear and interaction terms in the data, and quantifies their contributions to the models' performance. Using repeated nested cross-validation, the package:
  * Validates Cox Proportionate Hazards model, or Cox Lasso depending on the user's input. Uses 'survival::coxph()' [1,2] or 'glmnet::glmnet(..., family="cox")'[3] functions.
  * Validates Survival Random Forest ensembled with the baseline Cox model. Uses 'randomForestSRC::rfsrc()' [4] function.
  * Performs statistical testing of whether the Survival Random Forest ensemble outperforms the Cox model.

The performance metrics include [5]:
 * Discrimination measures: Harrell's concordance index, time-dependent AUCROC.
 * Calibration measures: calibration slope, calibration alpha.
 * Overall fit: Brier score, Scaled Brier score. 

### Getting started 
You can install the package from its github directory by running the `devtools::install_github("dianashams/survcompare")` command. The main function to use is `survcompare(data, predictors)`. The data should be in a form of a data frame, with "time" and "event" columns defining the survival outcome. A list of column names corresponding to the predictors to be used should also be supplied.

#### FAQ1: Why these (CoxPH and SRF) models? 
CoxPH model is a widely used survival model proved to be robust and easy to interprete. It assumes linear dependency of the log-hazards on the predictors; in its classical form, the effect of predictors is time-invariant which underlies the proportionality assumption. This  means that models' estimates are the averaged over time effects of predictors on the instant chances of survival. 

SRF is a machine learning algorithm that recursively splits the data into the sub-samples with similar survival profiles. It can deal with non-proportionate hazards and automatically captures non-linear and interaction terms, on top of the linear dependencies that CoxPH handles. However, it tends to overfit, and interpretation of radom forests' predictions is not straightforward especially for the survival data.

Given these qualities, SRF vs CoxPH's comparison is indicative of compex data dependencies, and quantifies the cost of using a simpler CoxPH model versus more flexible alternatives.

#### FAQ2: Why the ensemble and not just SRF? 
The ensemble of Cox and SRF takes the predictions of the Cox model and adds to the list of predictors to train SRF. This way, we make sure that linearity is captured by SRF at least as good as in the Cox model, and hence the marginal outperformance of the ensemble over the Cox model can be fully attributed to the qualities of SRF that Cox does not have, that is, data complexity.

#### FAQ3: How do I use the results? 
First, try to run sufficient number of repetitions (repeat_cv), at least 5, ideally 20-50 depending on the data heterogeneity and size.
There are two possible outcomes: "Survival Random Forest ensemble has outperformed CoxPH by ... in C-index", or "Survival Random Forest ensemble has NOT outperformed CoxPH". 
  * If there is **no outperformance**, this result can justify the employment of CoxPH model and indicate a negligible advantage of using a more flexible model such as Survival Random Forest.
  * In the case of **outperformance**, a researcher can 1) decide to go for a more complex model, 2) look for the interaction and non-linear terms that could be added to the CoxPH and re-run the test again, or 3) consider still using the CoxPH model if the difference is not large in the context of the performed task, or not enough to sacrifice model interpretability.

### Example:
The files in the "Example/" folder illustrate `survcompare`'s  application to the simulated and GBSG2  (https://rdrr.io/cran/pec/man/GBSG2.html) datasets. The outputs contain  internally-validated performance metrics along with the results of the statistical testing of whether Survival Random Forest outperforms the Cox Proportionate Hazard model (or Cox Lasso).  
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

### If you use the package or its code, please cite:

Shamsutdinova, D., Stamate, D., Roberts, A., & Stahl, D. (2022). Combining Cox Model and Tree-Based Algorithms to Boost Performance and Preserve Interpretability for Health Outcomes. In IFIP International Conference on Artificial Intelligence Applications and Innovations (pp. 170-181). Springer, Cham.

### Support or Contact
diana.shamsutdinova@kcl.ac.uk

### Links and references: 
##### References:

[1] Cox, D. R. (1972). Regression models and life‐tables. Journal of the Royal Statistical Society: Series B (Methodological), 34(2), 187-202.

[2] Therneau T (2023). *A Package for Survival Analysis in R*. R package version 3.5-7, <https://CRAN.R-project.org/package=survival>.

[3] Simon N, Friedman J, Tibshirani R, Hastie T (2011). "Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent." *Journal of Statistical Software, 39(5)*, 1--13. <doi:10.18637/jss.v039.i05>.

[4] Ishwaran H, Kogalur U (2023). *Fast Unified Random Forests for Survival, Regression, and Classification (RF-SRC).* R package version 3.2.2, <https://cran.r-project.org/package=randomForestSRC>

[5] Steyerberg EW, Vergouwe Y. (2014). Towards better clinical prediction models: seven steps for development and an ABCD for validation. *European heart journal, 35(29)*, 1925-1931 <https://doi.org/10.1093/eurheartj/ehu207>

[6] Shamsutdinova, D., Stamate, D., Roberts, A., & Stahl, D. (2022, June). Combining Cox Model and Tree-Based Algorithms to Boost Performance and Preserve Interpretability for Health Outcomes. In IFIP International Conference on Artificial Intelligence Applications and Innovations (pp. 170-181). Cham: Springer International Publishing. <https://www.researchgate.net/profile/Diana-Shamsutdinova-3/publication/361355831_Combining_Cox_Model_and_Tree-Based_Algorithms_to_Boost_Performance_and_Preserve_Interpretability_for_Health_Outcomes/links/650c47dcd5293c106ccb7043/Combining-Cox-Model-and-Tree-Based-Algorithms-to-Boost-Performance-and-Preserve-Interpretability-for-Health-Outcomes.pdf>
