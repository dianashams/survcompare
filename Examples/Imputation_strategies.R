
#'########################################################
#' July 2025,  code by Diana Shamsutdinova
#'
#' Illustrative example of using the survcompare package
#' to cross-validate survival models in the presence of missing data
#' using different imputation strategies defined by parameter "impute":
#' 
#'   - no imputation (impute = 0): no imputation to be performed, 
#'                                  the function is aborted if missing values 
#'                                  are in predictors, event or time.
#'                                  
#'   - proper imputation (impute = 1) : imputation by missForest is performed 
#'                                      in a proper way: during cross-validation, 
#'                                      missForest imputer is trained on the 
#'                                      train dataset, and used to impute both
#'                                       the train and the test. 
#'                                       
#'   - fast imputation (impute = 2): imputation by missForest in a fast way,
#'                                      where the entire data is imputed, 
#'                                      then nested cross-validation is performed.
#'                                      There is some leakage of information 
#'                                      from the test set into the train 
#'                                      set as imputed values in the train 
#'                                      set used all the data including test.
#'                                      
#'   - complete cases analysis (impute = 3): only instances (rows) with no missing 
#'                                    data in predict.factors, event, and time 
#'                                    columns are used.
#'   

#'#######################################################

# Simulate data 
n=1000
d0 = simulate_nonlinear(n, randomseed = 10000)

# Copy data and introduce some missing values 
d1 = d0
set.seed(111)
d1$age[sample(1:n, 220)] <- NA
d1$bmi[sample(1:n, 170)] <- NA
d1$hyp[sample(1:n, 250)] <- NA
mean(complete.cases(d1))*100 #48.7% are complete cases

# Cross-validate CoxPH with different imputation strategies

# impute = 0 for no imputation. This should stop with a warning (s)
s = survcox_cv(d1, names(d1)[1:4], repeat_cv = 5, impute = 0)

# impute = 0 for the original data without missing values (s0)
# impute = 1 for fast imputation (first impute, then cross-validate) (s1)
# impute = 2 for proper imputation (s2) (impute train and test for each split in external loop of cross-validation)
# impute = 3 for complete cases analysis (s3)
s0 = survcox_cv(d0, names(d1)[1:4], repeat_cv = 5, impute = 0, randomseed = 100)
s1 = survcox_cv(d1, names(d1)[1:4], repeat_cv = 5, impute = 1, randomseed = 100)
s2 = survcox_cv(d1, names(d1)[1:4], repeat_cv = 5, impute = 2, randomseed = 100)
s3 = survcox_cv(d1, names(d1)[1:4], repeat_cv = 5, impute = 3, randomseed = 100)

r1 = rbind(no_missing = s0$testaverage,
      impute_proper = s1$testaverage,
      complete_cases = s3$testaverage,
      impute_fast = s2$testaverage
      )
r1 
#                     T    AUCROC        BS  BS_scaled   C_score Calib_slope
# no_missing     8.1553 0.6597002 0.1631174 0.06983837 0.6385180   0.9712738
# impute_proper  8.1553 0.6348548 0.1668051 0.04871387 0.6133009   1.0062060
# complete_cases 8.1553 0.6307287 0.1646967 0.03777169 0.6109208   0.7905322
# impute_fast    8.1553 0.6910390 0.1602510 0.08609562 0.6619527   1.0143664

s0 = survcox_cv(d0, names(d1)[1:4], repeat_cv = 5, impute = 0, useCoxLasso = TRUE, randomseed = 100)
s1 = survcox_cv(d1, names(d1)[1:4], repeat_cv = 5, impute = 1, useCoxLasso = TRUE, randomseed = 100)
s2 = survcox_cv(d1, names(d1)[1:4], repeat_cv = 5, impute = 2, useCoxLasso = TRUE, randomseed = 100)
s3 = survcox_cv(d1, names(d1)[1:4], repeat_cv = 5, impute = 3, useCoxLasso = TRUE, randomseed = 100)

r2 = rbind(no_missing = s0$testaverage,
      impute_proper = s1$testaverage,
      complete_cases = s3$testaverage,
      impute_fast = s2$testaverage
)
r2 

#                     T    AUCROC        BS  BS_scaled   C_score Calib_slope
# no_missing     8.1553 0.6627159 0.1635507 0.06734287 0.6404243    1.166076
# impute_proper  8.1553 0.6414736 0.1670026 0.04758242 0.6190577    1.277912
# complete_cases 8.1553 0.6364252 0.1649744 0.03629009 0.6151604    1.057965
# impute_fast    8.1553 0.6953897 0.1599165 0.08801780 0.6665221    1.148341

s0 = survsrf_cv(d0, names(d1)[1:4], repeat_cv = 5, impute = 0, randomseed = 100)
s1 = survsrf_cv(d1, names(d1)[1:4], repeat_cv = 5, impute = 1, randomseed = 100)
s2 = survsrf_cv(d1, names(d1)[1:4], repeat_cv = 5, impute = 2, randomseed = 100)
s3 = survsrf_cv(d1, names(d1)[1:4], repeat_cv = 5, impute = 3, randomseed = 100)

r3 = rbind(no_missing = s0$testaverage,
      impute_proper = s1$testaverage,
      complete_cases = s3$testaverage,
      impute_fast = s2$testaverage
)
r3 

#                     T    AUCROC        BS  BS_scaled   C_score Calib_slope
# no_missing     8.1553 0.7559426 0.1460732 0.16684962 0.7233609   0.9442337
# impute_proper  8.1553 0.7208506 0.1552944 0.11449474 0.6875546   0.9433139
# complete_cases 8.1553 0.7381347 0.1547121 0.09579353 0.6965497   0.8112521
# impute_fast    8.1553 0.8110991 0.1369203 0.21831277 0.7708002   0.9436761

s0 = survsrfens_cv(d0, names(d1)[1:4], repeat_cv = 5, impute = 0, randomseed = 100)
s1 = survsrfens_cv(d1, names(d1)[1:4], repeat_cv = 5, impute = 1, randomseed = 100)
s2 = survsrfens_cv(d1, names(d1)[1:4], repeat_cv = 5, impute = 2, randomseed = 100)
s3 = survsrfens_cv(d1, names(d1)[1:4], repeat_cv = 5, impute = 3, randomseed = 100)

r4 = rbind(no_missing = s0$testaverage,
           impute_proper = s1$testaverage,
           complete_cases = s3$testaverage,
           impute_fast = s2$testaverage
)
r4 

#                     T    AUCROC        BS BS_scaled   C_score Calib_slope
# no_missing     8.1553 0.7417882 0.1526104 0.1290705 0.7082885   0.8755592
# impute_proper  8.1553 0.7245468 0.1545608 0.1185132 0.6908927   0.9864300
# complete_cases 8.1553 0.7385567 0.1533414 0.1040091 0.6967883   0.8041718
# impute_fast    8.1553 0.8002766 0.1401648 0.1989690 0.7598287   1.2498746






