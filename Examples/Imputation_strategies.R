devtools::load_all("~/Documents/GitHub/survcompare")

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
s0 = survcox_cv(d0, names(d1)[1:4], repeat_cv = 5, impute = 0)
s1 = survcox_cv(d1, names(d1)[1:4], repeat_cv = 5, impute = 1)
s2 = survcox_cv(d1, names(d1)[1:4], repeat_cv = 5, impute = 2)
s3 = survcox_cv(d1, names(d1)[1:4], repeat_cv = 5, impute = 3)

rbind(no_missing = s0$testaverage,
      complete_cases = s3$testaverage,
      impute_fast = s2$testaverage,
      impute_proper = s1$testaverage
      )

#                     T    AUCROC        BS  BS_scaled   C_score Calib_slope
# no_missing     8.1553 0.6566159 0.1638284 0.06604365 0.6339907   0.9098623
# complete_cases 8.1678 0.6424947 0.1657198 0.04537381 0.6213333   0.8307362
# impute_fast    8.7278 0.6381525 0.1645724 0.04701499 0.6152537   0.9828823
# impute_proper  8.7278 0.6370180 0.1639297 0.05237755 0.6144680   0.9698049
