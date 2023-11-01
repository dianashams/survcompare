library(devtools) 
load_all()
roxygen2::roxygenize()
load_all()

###### For Vignette 1 ######
library(survcompare)
mydata <- simulate_linear(200)
mypredictors <- names(mydata)[1:4]
compare_models_lasso <-
  survcompare(mydata, mypredictors, useCoxLasso = TRUE)
summary(compare_models_lasso)

mydata2 <- simulate_crossterms()
mypredictors2 <- names(mydata)[1:4]
compare_models2 <-
  survcompare(mydata2, mypredictors2, repeat_cv = 5)
summary(compare_models2)
round(compare_models_lasso$test$CoxPH, 4)
round(compare_models_lasso$test$SRF_ensemble, 4)
round(compare_models_lasso$main_stats, 4)
round(compare_models_lasso$results_mean, 4)

# GBSG #
#data(package = "survival")
#gbsg
psych::describe(gbsg)
names(gbsg)
mygbsg <- gbsg
mygbsg$time <- gbsg$rfstime / 365
mygbsg$event <- gbsg$status
myfactors <-
  c("age", "meno", "size", "grade", "nodes", "pgr", "er", "hormon")
mygbsg[1:15, myfactors]
sum(is.na(mygbsg[myfactors]))
survcompare(mygbsg, myfactors)

coxm <- survcox_train(mygbsg, myfactors)
t <- mygbsg[1:10, "time"]
mine = rep(NA, 10)
for (i in 1:10) {
  mine[i] = survcox_predict(coxm, mygbsg[i,], fixed_time = t[i])
}
s = 1-predict(coxm, mygbsg, type = "survival")[1:10]
data.frame(mine, "survival" = s, "time" = mygbsg[1:10, "time"], "diff" = mine - s)
#          mine   survival diff
# 1  0.48747623 0.48747623    0
# 2  0.25834711 0.25834711    0
# 3  0.69051210 0.69051210    0

pec <- 1-pec::predictSurvProb(coxm, mygbsg[1:10,],6)
mine6 <- as.vector(survcox_predict(coxm, mygbsg[1:10,], fixed_time = 6))
data.frame(mine6, pec, "time" = mygbsg[1:10, "time"], "diff" = round(mine6 - pec,10))
#        mine6       pec      time diff
# 1  0.5730218 0.5730218 5.0356164    0
# 2  0.9415808 0.9415808 1.1041096    0
# 3  0.8202664 0.8202664 4.3917808    0
# 4  0.6912045 0.6912045 0.4849315    0

###### For Vignette 2 ######

mytime = 5
coxm <- survcox_train(mygbsg, myfactors, mytime, useCoxLasso = TRUE)

# apparent validation: 
p5<- survcox_predict(coxm, mygbsg, mytime) #event probabilities
app_validity <- surv_validate(p5, mytime, mygbsg, mygbsg)
round(app_validity, 4)  
#   T AUCROC     BS BS_scaled C_score Calib_slope Calib_alpha
# 1 5 0.7287 0.2156    0.1725  0.6861      1.6904       0.726

#sub-group validation
#let's validate the model for a subgroup of people aged 65 and over
older_group <- mygbsg[mygbsg$age>=65, ] # 92 people 

#compute 5 year probability of event 
cox_prob5 <- survcox_predict(coxm, older_group, mytime)

# compute performance statistics (validate)
group_val <- 
  surv_validate(cox_prob5, mytime, mygbsg, older_group)

round(group_val, 4)
#   T AUCROC     BS BS_scaled C_score Calib_slope Calib_alpha
# 1 5 0.7884 0.2046    0.2197  0.7346       1.505      0.4373

mytime = 5
coxm <- survcox_train(mygbsg, myfactors, mytime, useCoxLasso = TRUE)

# apparent validation: 
p5<- survcox_predict(coxm, mygbsg, mytime) #event probabilities
app_validity <- surv_validate(p5, mytime, mygbsg, mygbsg)
round(app_validity, 4)  
#   T AUCROC     BS BS_scaled C_score Calib_slope Calib_alpha
# 1 5 0.7287 0.2156    0.1725  0.6861      1.6904       0.726

#sub-group validation
#let's validate the model for a subgroup of people aged 65 and over
older_group <- mygbsg[mygbsg$age>=65, ] # 92 people 

#compute 5 year probability of event 
cox_prob5 <- survcox_predict(coxm, older_group, mytime)

# compute performance statistics (validate)
group_val <- 
  surv_validate(cox_prob5, mytime, mygbsg, older_group)

round(group_val, 4)
#   T AUCROC     BS BS_scaled C_score Calib_slope Calib_alpha
# 1 5 0.7884 0.2046    0.2197  0.7346       1.505      0.4373

####  vignette 2 SRF ####

srfmodel1 <- 
  survsrf_train(mygbsg, myfactors, mytime, inner_cv = 3,randomseed = 42)
survcox_predict()
srfmodel1$tuning
# $mtry
# [1] 2 3 4 5
# 
# $nodesize
# [1] 15 20 25 30 35 40 45 50
# 
# $nodedepth
# [1]  5 25

srfmodel1$beststats
#    mtry nodesize nodedepth time   AUCROC       BS BS_scaled
# V1    2       25        25    5 0.725787 0.215723 0.1721251
#      C_score Calib_alpha Calib_slope
# V1 0.6929181   0.7828155    1.271496

#try to tune the model further or customize SRF tuning
# mtry was chosen 2 out of 2-5 range, so we will leave it as 2 
# nodesize 25 from 15-50 range, so we can leave it as such
#nodedepth was the higher number from 5 and 25, 
# we will add even higher values and values between 5 and 25

my_srf_tuning<- 
  list("mtry"= c(2,3), "nodesize" = c(15,25), 
       "nodedepth" = c(10, 15, 20,25, 30))
srfmodel2<- 
  survsrf_train(mygbsg, myfactors, mytime, srf_tuning = my_srf_tuning, 
                inner_cv = 3,randomseed = 42)
#Now we can compare the results of the default and customised SRF tuning: 
rbind(default_tune = srfmodel1$beststats, 
      custom_tune = srfmodel2$beststats)
#=> seem largely similar.

# SRF (apparent) validation for older_group
# Note that we measured performance for the data that 
# was used for training, so this is still an apparent performance

srf_p6<- survsrf_predict(srfmodel2, older_group,mytime)
group_val_srf<- surv_validate(srf_p6,mytime,mygbsg, older_group)
round(group_val_srf,4)

rbind("CoxPH_for_oldergroup" = round(group_val,4),
      "SRF_for_oldergroup" = round(group_val_srf,4))


#######################
coxcv <- survcox_cv(mygbsg, myfactors, fixed_time = 5)
summary(coxcv)

# Cross-validation results
# Call:
#   survcox_cv(df = mygbsg, predict.factors = myfactors, fixed_time = 5)
# 
# The stats are computed from the  6  data splits.
#             test.mean test.sd test.95CIHigh test.95CILow
# T              5.0000  0.0000        5.0000       5.0000
# AUCROC         0.7368  0.0351        0.6950       0.6841
# BS             0.2142  0.0287        0.1805       0.1697
# BS_scaled      0.1829  0.0596        0.1229       0.1227
# C_score        0.6817  0.0164        0.6640       0.6571
# Calib_slope    1.2789  0.1514        1.1291       1.1058
# Calib_alpha    0.7679  0.1737        0.5854       0.5563
#             train.mean train.sd train.95CIHigh train.95CILow
# T               5.0000   0.0000         5.0000        5.0000
# AUCROC          0.7421   0.0162         0.7248        0.7241
# BS              0.2105   0.0065         0.2033        0.2011
# BS_scaled       0.1911   0.0229         0.1684        0.1624
# C_score         0.6894   0.0045         0.6844        0.6841
# Calib_slope     1.3699   0.1557         1.2184        1.1782
# Calib_alpha     0.7584   0.0720         0.6796        0.6561


######################


coxcvlin <- survcox_cv(mydata, names(mydata)[1:4],10 )
summary(coxcvlin)

# Cross-validation results
# Call:
#   survcox_cv(df = mydata, predict.factors = names(mydata)[1:4], 
#              fixed_time = 10)
# 
# The stats are computed from the  6  data splits.
# test.mean test.sd test.95CIHigh test.95CILow
# T             10.0000  0.0000       10.0000      10.0000
# AUCROC         0.8418  0.0713        0.7749       0.7351
# BS             0.0923  0.0360        0.0638       0.0630
# BS_scaled      0.4552  0.1653        0.2754       0.1973
# C_score        0.8003  0.0593        0.7317       0.7046
# Calib_slope    1.2403  0.7339        0.7522       0.6591
# Calib_alpha    0.1372  0.6245       -0.3996      -0.4511
# train.mean train.sd train.95CIHigh train.95CILow
# T              10.0000   0.0000        10.0000       10.0000
# AUCROC          0.8545   0.0284         0.8241        0.8112
# BS              0.0855   0.0084         0.0760        0.0759
# BS_scaled       0.4906   0.0386         0.4538        0.4337
# C_score         0.8133   0.0215         0.7905        0.7894
# Calib_slope     1.1480   0.0652         1.0829        1.0803
# Calib_alpha     0.1528   0.0454         0.1044        0.0958


############ TESTS ######## 

usethis::use_test("name") 

#######   Vignettes' theory #######
# first vignette:
#Creates a vignettes/ directory.
#Adds the necessary dependencies to DESCRIPTION, i.e. adds knitr to the VignetteBuilder field and adds both knitr and rmarkdown to Suggests.
#Drafts a vignette, vignettes/my-vignette.Rmd.Adds some patterns to .gitignore to ensure that files created as a side effect of previewing your vignettes are kept out of source control (we’ll say more about this later).


usethis::use_vignette("my-vignette")

# Other vignettes:
# You also call use_vignette() to create your second and all subsequent vignettes; it will just skip any setup that’s already been done.

use_vignette("Cross_Validation_using_survcompare")

# Checking vignette: 
# An option is to use devtools::build_rmd("vignettes/my-vignette.Rmd") 
# to render the vignette. 
# This builds your vignette against a (temporarily installed) development 
# version of your package.

##########################

type <- tolower(type)
if(any(type %in% c("lp","hazard","cumhazard","survival") == FALSE)){
  stop("type can only be \"lp\", \"hazard\", \"cumhazard\" or/and \"survival\" \n") 
}


