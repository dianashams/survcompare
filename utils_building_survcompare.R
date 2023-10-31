library(devtools) 
load_all()
roxygen2::roxygenize()
load_all()


library(survcompare)
mydata <- simulate_linear()

mypredictors <- names(mydata)[1:4]

compare_models_lasso <-
  survcompare(mydata, mypredictors, useCoxLasso = TRUE)

summary(compare_models_lasso)

mydata2 <- simulate_crossterms()

mypredictors2 <- names(mydata)[1:4]

compare_models2 <-
  survcompare(mydata2, mypredictors2, repeat_cv = 5)

summary(compare_models2)

data(package = "survival")
gbsg
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

round(compare_models_lasso$test$CoxPH, 4)
round(compare_models_lasso$test$SRF_ensemble, 4)

round(compare_models_lasso$main_stats, 4)
round(compare_models_lasso$results_mean, 4)

survcompare(sdlkjfgh, dklsjfgh)
coxph(dklsjgh, dflsjhg)

coxm <- survcox_train(mygbsg, myfactors)
m2<- rms::cph(Surv(time, event)~age+meno+size+grade+nodes+pgr+er+hormon, data = mygbsg)

t <- mygbsg[1:10, "time"]
mine = rep(NA, 10)
for (i in 1:10) {
  mine[i] = survcox_predict(coxm, mygbsg[i,], fixed_time = t[i])
}
s = 1-predict(coxm, mygbsg, type = "survival")[1:10]
pec = 1-pec::predictSurvProb(coxm, mygbsg[1:10,],t[3])
data.frame(mine,s,pec, "time" = mygbsg[1:10, "time"])

coxcv <- survcox_cv(mygbsg, myfactors)

riskRegression::predictCox
coxm <- survcox_cv()


############ TESTS ######## 

usethis::use_test("name") 

#######   Vignettes #######
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
