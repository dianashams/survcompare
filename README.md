# R package "survcompare": 

**Comparing linear (Cox Proportionate Hazards model) and non-linear (Survival Random Forest) survival models to quantify predictive value of non-linear and interaction terms in the data.**

![image](https://github.com/dianashams/ensemble-methods-for-survival-analysis/blob/gh-pages/survcompare_cartoon.png)
### Method: 
The primary goal is to assist researchers in making informed decisions regarding whether they should choose a flexible yet less transparent machine learning approach or a traditional linear method. This is achieved by examining the presence of non-linear and interaction terms within the data and quantifying their impact on the models' performance. 

The package fits, tunes, and internaly validates the baseline Cox Proportionate Hazards model (CoxPH), and compares its performance to one of the machine learning alternatives (Survival Random Forest (SRF), or its ensemble with CoxPH). A comparison with the deep learning model DeepHit is available through the package's GitHub version. 
The package 
   * Fits and tunes the underlying models using cross-validated random hyperparameter search 
   * Performs statistical testing to compare the outperformance of the SRF (or SRF-CoxPH ensemble) over the baseline Cox-PH model
     
In more details, the baseline models are: 
  * Cox Proportionate Hazards model. The underlying model is 'survival::coxph()' [1,2]. 
  * Cox-Lasso regularized version of the CoxpH. The underlying model is 'glmnet::glmnet(..., family="cox")'[3].

The alternatives are:  
  * Survival Random Forest (SRF) model, the underlying  model is 'randomForestSRC::rfsrc()' [4]. 
  * Sequential ensemble of the SRF with the baseline Cox model [6]. The ensemble takes CoxPH predictions and supplied them to  Survival Random Forest as an additional predictor, see more explanations below.
  * Stacked ensemble of the CoxPH and Survival Random Forest, $`\lambda CoxPH + (1-\lambda) SRF`$. Lambda parameter is tuned within the package and the value shows what share of SRF predictions can improve the baseline CoxPH performance. $`\lambda = 0`$ means only CoxPH is used, $`\lambda = 1`$  means the model only relies on Survival Random Forest.
  * GitHit only: deep learning model DeepHit 'survivalmodels::deephit()', as well as its sequential and stacked ensembles with CoxPH or Cox-Lasso. 

The performance metrics include [5]:
 * Discrimination measures: Harrell's concordance index, time-dependent AUCROC.
 * Calibration measures: calibration slope, calibration alpha.
 * Overall fit: Brier score, Scaled Brier score. 

NB: Sequential ensemble is the first ensemble method described in https://dianashams.github.io/ensemble-methods-for-survival-analysis/ as published in Shamsutdinova, Stamate, Roberts, & Stahl (2022, June) [6]. 

### Getting started 
You can install the package from CRAN as `install.packages("survcompare")`, or from its github directory by running the `devtools::install_github("dianashams/survcompare")` command. The main function to use is `survcompare(data, predictors)`. The data should be in a form of a data frame, with "time" and "event" columns defining the survival outcome. A list of column names corresponding to the predictors to be used should also be supplied.

#### FAQ1: Why these (CoxPH and SRF) models? 
CoxPH model is a widely used survival model proved to be robust and easy to interprete. It assumes linear dependency of the log-hazards on the predictors; in its classical form, the effect of predictors is time-invariant which underlies the proportionality assumption. This  means that models' estimates are the averaged over time effects of predictors on the instant chances of survival. 

SRF is a machine learning algorithm that recursively splits the data into the sub-samples with similar survival profiles. It can deal with non-proportionate hazards and automatically captures non-linear and interaction terms, on top of the linear dependencies that CoxPH handles. However, it can overfit, especially in smaller datasets, often seen in clinical data, and interpretation of random forests' predictions is not straightforward especially for the survival data.

Given these qualities, SRF vs CoxPH's comparison is indicative of compex data dependencies, and quantifies the cost of using a simpler CoxPH model versus more flexible alternatives.

#### FAQ2: Why the ensemble and not just SRF? 
First, you can use the package to compare the performances of the CoxPH and SRF themselves. 

Second, the ensembles aim to single out the predictive value of the non-linearities and other data relationships that could not be captured by the baseline models. In both ensembles, the final models has a direct access to the predictions of the baseline CoxPH, and hence, the outperformance can be fully attributed to such complex relationships. 

For example, the sequential ensemble of Cox and SRF takes the predictions of the Cox model and adds to the list of predictors to train SRF. This way, we make sure that linearity is captured by SRF at least as good as in the Cox model, and hence the marginal outperformance of the ensemble over the Cox model can be fully attributed to the qualities of SRF that Cox does not have, that is, data complexity.

#### FAQ3: How do I interpret and use the results? 
First, try to run sufficient number of repetitions (repeat_cv), at least 5, ideally 20-50 depending on the data heterogeneity and size.
There are two possible outcomes: "Survival Random Forest ensemble has outperformed CoxPH by ... in C-index", or "Survival Random Forest ensemble has NOT outperformed CoxPH". 
  * If there is **no outperformance**, this result can justify the employment of CoxPH model and indicate a negligible advantage of using a more flexible model such as Survival Random Forest.
  * In the case of **outperformance**, a researcher can 1) decide to go for a more complex model, 2) look for the interaction and non-linear terms that could be added to the CoxPH and re-run the test again, or 3) consider still using the CoxPH model if the difference is not large in the context of the performed task, or not enough to sacrifice model interpretability.

### Example:
```R
mydata <- simulate_crossterms()
mypredictors <- names(mydata)[1:4]
compare_models <- survcompare(mydata, mypredictors, fixed_time = 9)

# [1] "Cross-validating CoxPH using 2 repeat(s), 3 outer, 3 inner loops)."
# [1] "Repeated CV 1 / 2"
# |====================================================================================| 100%
# [1] "Repeated CV 2 / 2"
# |====================================================================================| 100%
# Time difference of 0.4649661 secs
# [1] "Cross-validating Survival Random Forest using 2 repeat(s), 3 outer, 3 inner loops)."
# [1] "Repeated CV 1 / 2"
# |====================================================================================| 100%
# [1] "Repeated CV 2 / 2"
# |====================================================================================| 100%
# Time difference of 9.842596 secs
# Internally validated test performance of CoxPH and Survival Random Forest over 2 repeated 3 fold cross-validations (inner k = 3 ). Mean performance:
#   T C_score AUCROC Calib_slope  sec
# CoxPH                    9  0.6774 0.7096      0.8407 0.46
# Survival Random Forest   9  0.6974 0.7277      1.0101 9.84
# Diff                     0  0.0200 0.0182      0.1695 9.38
# pvalue                 NaN  0.0602 0.1177      0.1564  NaN
# 
# Median performance:
#   T C_score AUCROC Calib_slope  sec
# CoxPH                    9  0.6623 0.7054      0.7447 0.46
# Survival Random Forest   9  0.7042 0.7514      1.0164 9.84
# Diff                     0  0.0419 0.0460      0.2717 9.38
# pvalue                 NaN  0.0602 0.1177      0.1564  NaN
# 
# Survival Random Forest has NOT outperformed CoxPH with the mean c-index difference of 0.02.
# The difference is not statistically significant with the p-value = 0.0602. 
# The data may NOT contain considerable non-linear or cross-term dependencies
# that could be captured by Survival Random Forest.
# Mean C-score: 
#   CoxPH  0.6774(95CI=0.6737-0.6811;SD=0.0055)
# Survival Random Forest 0.6974(95CI=0.6791-0.7157;SD=0.0272)
# Mean AUCROC:
#   CoxPH  0.7096(95CI=0.6949-0.7242;SD=0.0218)
# Survival Random Forest 0.7277(95CI=0.6917-0.7638;SD=0.0536)

round(compare_models$main_stats_pooled,4)

#                                  mean     sd 95CILow 95CIHigh
# C_score_CoxPH                  0.6774 0.0055  0.6737   0.6811
# C_score_Survival Random Forest 0.6974 0.0272  0.6791   0.7157
# AUCROC_CoxPH                   0.7096 0.0218  0.6949   0.7242
# AUCROC_Survival Random Forest  0.7277 0.0536  0.6917   0.7638

# -------------  Stacked ensemble: -------------
cvstack <- survsrfstack_cv(mydata2, mypredictors2, randomseed = 100, repeat_cv = 3)
# get lambdas:
unlist(cvstack$bestparams$lambda)
#[1] 0.99 1.00 1.00 0.78 0.98 0.55 0.44 0.96 0.60
# mean lambda
mean(unlist(cvstack$bestparams$lambda)) #0.811 - the meta-learner mostly relies on SRF 

# Compare stacked ensemble performance to the basline CoxLasso using survcompare2() function:
cv1 <- survcox_cv(mydata2, mypredictors2, randomseed = 100, repeat_cv = 3, useCoxLasso = TRUE)
compare2 <- survcompare2(cv1, cvstack)

# Internally validated test performance of CoxLasso and Stacked_SRF_CoxPH over 
# 3 repeated 3 fold cross-validations (inner k = 3 ). Mean performance:
#                     T C_score AUCROC Calib_slope   sec
# CoxLasso            9  0.6377 0.6395      1.4429  0.69
# Stacked_SRF_CoxPH   9  0.7720 0.8119      1.0922 14.09
# Diff                0  0.1343 0.1723     -0.3508 13.40
# pvalue            NaN  0.0000 0.0000      0.8327   NaN
# 
# Stacked_SRF_CoxPH has outperformed CoxLassoby 0.1343 in C-index.
# The difference is statistically significant with the p-value 1.27e-06***.
# The supplied data may contain non-linear or cross-term dependencies, 
# better captured by Stacked_SRF_CoxPH.
# Mean C-score: 
#   CoxLasso  0.6377(95CI=0.6312-0.6479;SD=0.0096)
# Stacked_SRF_CoxPH 0.772(95CI=0.7543-0.7885;SD=0.018)
# Mean AUCROC:
#   CoxLasso  0.6395(95CI=0.6325-0.6492;SD=0.0091)
# Stacked_SRF_CoxPH 0.8119(95CI=0.7885-0.8296;SD=0.0224)

```
NB: More examples are located in the "Example/" folder. 

### If you use the package or its code, please cite:
Shamsutdinova, D., Stamate, D., Roberts, A., & Stahl, D. (2022). Combining Cox Model and Tree-Based Algorithms to Boost Performance and Preserve Interpretability for Health Outcomes. In IFIP International Conference on Artificial Intelligence Applications and Innovations (pp. 170-181). Springer, Cham.

### Support or Contact
If you have any comments, suggestions, corrections or enchancements, kindly submit an issue on the
<https://github.com/dianashams/survcompare/issues> or email to diana.shamsutdinova.github@gmail.com.

### Disclaimer
This R package is offered free and without warranty of any kind, either expressed or implied. The package authors will not be held liable to you for any damage arising out of the use, modification or inability to use this program. This R package can be used, redistributed and/or modified freely for non-commercial purposes subject to the original source being properly cited. Licensed under GPL-3.

The authors received financial support by the National Institute for Health Research (NIHR) Biomedical Research Centre at South London and Maudsley NHS Foundation Trust and King’s College London. The views expressed are those of the author(s) and not necessarily those of the NHS, the NIHR or the Department of Health.

### Links and references: 
##### References:

[1] Cox, D. R. (1972). Regression models and life‐tables. Journal of the Royal Statistical Society: Series B (Methodological), 34(2), 187-202.

[2] Therneau T (2023). *A Package for Survival Analysis in R*. R package version 3.5-7, <https://CRAN.R-project.org/package=survival>.

[3] Simon N, Friedman J, Tibshirani R, Hastie T (2011). "Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent." *Journal of Statistical Software, 39(5)*, 1--13. <doi:10.18637/jss.v039.i05>.

[4] Ishwaran H, Kogalur U (2023). *Fast Unified Random Forests for Survival, Regression, and Classification (RF-SRC).* R package version 3.2.2, <https://cran.r-project.org/package=randomForestSRC>

[5] Steyerberg EW, Vergouwe Y. (2014). Towards better clinical prediction models: seven steps for development and an ABCD for validation. *European heart journal, 35(29)*, 1925-1931 <https://doi.org/10.1093/eurheartj/ehu207>

[6] Shamsutdinova, D., Stamate, D., Roberts, A., & Stahl, D. (2022, June). Combining Cox Model and Tree-Based Algorithms to Boost Performance and Preserve Interpretability for Health Outcomes. In IFIP International Conference on Artificial Intelligence Applications and Innovations (pp. 170-181). Cham: Springer International Publishing.
