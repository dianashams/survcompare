#' Cross-validates and compares Cox Proportionate Hazards and Survival Random Forest models
#' @description
#' The function performs a repeated nested cross-validation for
#' 1) Cox-PH (survival package, survival::coxph) or Cox-Lasso (glmnet package, glmnet::cox.fit)
#' 2) Survival Random Forest (randomForestSRC::rfsrc), or its ensemble with the Cox model (if use_ensemble =TRUE)
#'
#' The same random seed for the train/test splits are used for all models to aid fair comparison;
#' and the performance metrics are computed for the tree models including Harrel's c-index,
#' time-dependent AUC-ROC, time-dependent Brier Score, and calibration slope.
#' The statistical significance of the performance differences between Cox-PH and Cox-SRF Ensemble is tested and reported.
#'
#' The function is designed to help with the model selection by quantifying the loss of predictive
#' performance (if any) if Cox-PH is used instead of a more complex model such as SRF
#' which can capture non-linear and interaction terms, as well as non-proportionate hazards.
#' The difference in performance of the Ensembled Cox and SRF and the baseline Cox-PH
#' can be viewed as quantification of the non-linear and cross-terms contribution to
#' the predictive power of the supplied predictors.
#'
#' The function is a wrapper for survcompare2(), for comparison of the CoxPH and SRF models, and
#' an alternative way to do the same analysis is to run survcox_cv() and survsrf_cv(), then using survcompare2()
#'
#' @importFrom survival Surv
#' @importFrom timeROC timeROC
#' @importFrom stats as.formula
#' @importFrom stats quantile
#' @importFrom stats runif
#' @importFrom stats coef
#' @importFrom stats binomial
#' @importFrom stats poisson
#' @importFrom stats predict
#' @importFrom stats sd
#' @importFrom stats approxfun
#' @importFrom stats median
#' @importFrom randomForestSRC rfsrc
#' @importFrom survival coxph
#' @param df_train training data, a data frame with "time" and "event" columns to define the survival outcome
#' @param predict_factors list of column names to be used as predictors
#' @param fixed_time prediction time of interest. If NULL, 0.90th quantile of event times is used
#' @param randomseed random seed for replication
#' @param useCoxLasso TRUE / FALSE, for whether to use regularized version of the Cox model, FALSE is default
#' @param outer_cv k in k-fold CV
#' @param inner_cv k in k-fold CV for internal CV to tune survival random forest hyper-parameters
#' @param tuningparams list of tuning parameters for random forest: 1) NULL for using a default tuning grid, or 2) a list("mtry"=c(...), "nodedepth" = c(...), "nodesize" = c(...))
#' @param return_models TRUE/FALSE to return the trained models; default is FALSE, only performance is returned
#' @param repeat_cv if NULL, runs once, otherwise repeats several times with different random split for CV, reports average of all
#' @param ml this is currently for Survival Random Forest only ("SRF")
#' @param use_ensemble TRUE/FALSE for whether to train SRF on its own, apart from the CoxPH->SRF ensemble. Default is FALSE as there is not much information in SRF itself compared to the ensembled version.
#' @param max_grid_size number of random grid searches for model tuning
#' @param suppresswarn TRUE/FALSE, TRUE by default
#' @return outcome - cross-validation results for CoxPH, SRF, and an object containing the comparison results
#' @examples
#' \dontshow{rfcores_old <- options()$rf.cores; options(rf.cores=1)}
#' df <-simulate_nonlinear(100)
#' predictors <- names(df)[1:4]
#' srf_params <- list("mtry" = c(2), "nodedepth"=c(25), "nodesize" =c(15))
#' mysurvcomp <- survcompare(df, predictors, tuningparams = srf_params, max_grid_size = 1)
#' summary(mysurvcomp)
#' \dontshow{options(rf.cores=rfcores_old)}
#' @export
survcompare <- function(df_train,
                           predict_factors,
                           fixed_time = NaN,
                           randomseed = NaN,
                           useCoxLasso = FALSE,
                           outer_cv = 3,
                           inner_cv = 3,
                           tuningparams = list(),
                           return_models = FALSE,
                           repeat_cv = 2,
                           ml = "SRF",
                           use_ensemble = FALSE,
                           max_grid_size = 10,
                        suppresswarn = TRUE) {

  Call <- match.call()
  inputs <- list(
    df_train,
    predict_factors,
    fixed_time,
    outer_cv,
    inner_cv,
    repeat_cv,
    randomseed,
    return_models,
    tuningparams,
    useCoxLasso
  )

  inputclass <-
    list(
      df_train = "data.frame",
      predict_factors = "character",
      fixed_time = "numeric",
      outer_cv = "numeric",
      inner_cv = "numeric",
      repeat_cv = "numeric",
      randomseed = "numeric",
      return_models = "logical",
      tuningparams = "list",
      useCoxLasso = "logical"
    )
  cp <- check_call(inputs, inputclass, Call)
  if (cp$anyerror)
    stop (paste(cp$msg[cp$msg != ""], sep = ""))

  if (sum(is.na(df_train[c("time", "event", predict_factors)])) > 0) {
    stop("Missing data can not be handled. Please impute first.")
  }

  if (is.nan(randomseed)) {
    randomseed <- round(stats::runif(1) * 1e9, 0) + 1
  }
  if (is.nan(fixed_time)) {
    fixed_time <- quantile(df_train[df_train$event == 1, "time"], 0.9, na.rm = TRUE)
  }
  #SRF_ensemble or DeepHit_ensemble
  ensemble_name <- paste(ml, "ensemble", sep="_")

  if (suppresswarn){ user_warn <-options()$warn; options(warn=-1)}

  # CoxPH
  cv1 <- survcox_cv(
    df = df_train,
    predict.factors = predict_factors,
    fixed_time = fixed_time ,
    outer_cv = outer_cv,
    randomseed = randomseed,
    useCoxLasso = useCoxLasso,
    return_models = return_models,
    repeat_cv = repeat_cv
  )

  if(use_ensemble){
    cv2 <- survsrfens_cv(
      df = df_train,
      predict.factors = predict_factors,
      fixed_time = fixed_time,
      outer_cv = outer_cv,
      inner_cv = inner_cv,
      randomseed = randomseed,
      return_models = return_models,
      tuningparams = tuningparams,
      useCoxLasso = useCoxLasso,
      repeat_cv = repeat_cv,
      max_grid_size = max_grid_size
    )
  } else{
    cv2 <- survsrf_cv(
      df = df_train,
      predict.factors = predict_factors,
      fixed_time = fixed_time,
      outer_cv = outer_cv,
      inner_cv = inner_cv,
      repeat_cv = repeat_cv,
      randomseed = randomseed,
      return_models = return_models,
      tuningparams = tuningparams,
      max_grid_size = max_grid_size)
  }

  if (suppresswarn){ options(warn=user_warn)}

  output<- survcompare2(cv1, cv2)
  output$cv1 <- cv1
  output$cv2 <- cv2

  summary.survcompare(output)

  return(output)
}
