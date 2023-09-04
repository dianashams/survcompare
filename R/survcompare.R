
#' Cross-validates predictive performance of the baseline Cox-PH model and Survival random forest
#'
#' @description
#' The following models are evaluated:
#' 1) Cox-PH (survival package, survival::coxph)
#' 2) Cox-Lasso (glmnet package, glmnet::cox.fit)
#' 3) Survival random forest, SRF (randomForestSRC::rfsrc)
#' 4) Ensemble of Cox-PH or Cox-Lasso and SRF
#'
#' The same random seed for the train/test splits are used for all models.
#' The comparison can aid model selection and quantify the loss of predictive
#' performance (if any) if Cox-PH is used instead of a more complex model such as SRF
#' The difference in performance of the Ensembled Cox and SRF and the baseline Cox-PH
#' can be viewed as quantification of the non-linear and cross-terms contribution to
#' the predictive power of the supplied predictors
#'
#' @param df_train data frame
#' @param predict_factors list of predictors to use
#' @param predict_t prediction time of interest. If NULL, 90% quantile of event times is used
#' @param randomseed random seed for replication
#' @param useCoxLasso TRUE / FALSE, FALSE by default
#' @param outer_cv k in k-fold CV
#' @param inner_cv k in k-fold CV for internal tuning of survival random forest
#' @param srf_tuning list of tuning parameters for srf, e.g. list("mtry"=c(3,6,16), "ndepth" = c(3,6,19))
#' @param return_models TRUE/FALSE, if returns all the trained models; default is FALSE, only performance is returned
#' @param repeat_cv if NULL, runs once, otherwise repeats several times, reports average of all
#' @param train_srf TRUE/FALSE, default is FALSE, as Ensemble 1 by definition is very similar
#' @return outcome = list(data frame with performance results, fitted Cox models, fitted SRF)
#' @export
survcompare <- function(df_train,
                        predict_factors,
                        predict_t = NULL,
                        randomseed = NULL,
                        useCoxLasso = FALSE,
                        outer_cv = 5,
                        inner_cv = 3,
                        srf_tuning = NaN,
                        return_models = FALSE,
                        repeat_cv = 5,
                        train_srf= FALSE) {

  stopifnot(expr = {
    is.data.frame(df_train)
    predict_factors %in% colnames(df_train)
  })
  if (is.null(randomseed)) {
    randomseed <- round(stats::runif(1) * 1e9, 0) + 1
  }
  if (is.null(predict_t)) {
    predict_t <- quantile(df_train[df_train$event == 1, "time"], 0.9)
  }
  # __________________________________________________
  # cross-validating
  cox_cv <- survcox_cv(
    df = df_train,
    predict.factors = predict_factors,
    fixed_time=predict_t ,
    cv_number = outer_cv,
    randomseed = randomseed,
    useCoxLasso = useCoxLasso,
    return_models = return_models,
    repeat_cv = repeat_cv
  )
  if (train_srf){
    srf_cv <- survrf_cv(
      df = df_train,
      predict.factors = predict_factors,
      fixed_time= predict_t,
      cv_number = outer_cv,
      inner_cv = inner_cv,
      randomseed = randomseed,
      return_models = return_models,
      srf_tuning = srf_tuning,
      repeat_cv = repeat_cv
    )
  }
  ens1_cv <- survensemble1_cv(
    df = df_train,
    predict.factors = predict_factors,
    fixed_time=predict_t,
    cv_number = outer_cv,
    inner_cv = inner_cv,
    randomseed = randomseed,
    return_models = return_models,
    srf_tuning = srf_tuning,
    useCoxLasso = useCoxLasso,
    repeat_cv = repeat_cv
  )
  # __________________________________________________
  # gathering results for the output: test&train performance
  results_mean <- as.data.frame(
    rbind(
      "CoxPH" = cox_cv$testaverage,
      #      "SRF" = srf_cv$testaverage,
      "SRF_ensemble" = ens1_cv$testaverage
    )
  )
  results_mean$sec <-
    round(as.numeric(c(
      cox_cv$time,
      # srf_cv$time,
      ens1_cv$time)), 2)
  results_mean <-
    results_mean[c("T",
                   "C_score",
                   "AUCROC",
                   "BS",
                   "BS_scaled",
                   "Calib_slope",
                   "Calib_alpha",
                   "sec")]
  results_mean_train <- as.data.frame(
    rbind(
      "CoxPH" = cox_cv$trainaverage,
      # "SRF" = srf_cv$trainaverage,
      "SRF_ensemble" = ens1_cv$trainaverage
    ))
  results_mean_train$sec <- results_mean$sec
  results_mean_train <-
    results_mean_train[c("T",
                         "C_score",
                         "AUCROC",
                         "BS",
                         "BS_scaled",
                         "Calib_slope",
                         "Calib_alpha",
                         "sec")]
  # _____________________________________________________
  # testing outperformance of Ensemble 1 over CoxPH and CoxLasso
  t_coxph <-
    difftest(ens1_cv$test,
             cox_cv$test,
             dim(df_train)[1],
             length(predict_factors))

  t_coxph_train <-
    difftest(ens1_cv$train,
             cox_cv$train,
             dim(df_train)[1],
             length(predict_factors))

  # adding results line for the differences with Cox-PH
  results_mean["Diff", ]=results_mean[2, ]-results_mean[1,]
  results_mean_train["Diff", ]=results_mean[2, ]-results_mean[1,]
  results_mean["pvalue", ]= c(t_coxph[3,], NaN) #NaN for "sec" column
  results_mean_train["pvalue", ]= c(t_coxph_train[3,], NaN)
  # __________________________________________________

  stats_ci <- function(x, col = "C_score") {
    temp <- x[, col]
    c(
      "mean" = mean(temp, na.rm = 1),
      "sd" = sd(temp, na.rm = 1),
      "95CILow" = unname(quantile(temp, 0.025)),
      "95CIHigh" = unname(quantile(temp, 0.975))
    )
  }
  auc_c_stats <-
    rbind(
      "CoxPH_______C_score" = stats_ci(cox_cv$test,  "C_score"),
      "SRFensemble_C_score" = stats_ci(ens1_cv$test, "C_score"),
      "CoxPH________AUCROC" = stats_ci(cox_cv$test,  "AUCROC"),
      "SRFensemble__AUCROC" = stats_ci(ens1_cv$test, "AUCROC")
    )

  # output
  output <- list()
  output$results_mean <- results_mean
  output$results_mean_train <- results_mean_train
  #output$results_apparent <- results_apparent
  output$return_models <- list(
    "CoxPH" = cox_cv$tuned_cv_models,
    #"SRF" = srf_cv$tuned_cv_models,
    "SRF_ensemble" = ens1_cv$tuned_cv_models
  )
  output$test <- list(
    "CoxPH" = cox_cv$test,
    # "SRF" = srf_cv$test,
    "SRF_ensemble" = ens1_cv$test
  )
  output$train <- list(
    "CoxPH" = cox_cv$train,
    # "SRF" =  srf_cv$train,
    "SRF_ensemble" = ens1_cv$train
  )
  output$difftest <- t_coxph
  output$main_stats <-auc_c_stats
  output$randomseed <- randomseed
  summary.ensemblesurv(output)
  return(output)
}

#' Testing statistical significance of the ensembled model outperformance
#' over the baseline model (Model 1 vs Model 0)
#'
#' @description
#' Testing statistical significance of the outperformance of
#' Model 1 over Model 0
#' T-test is used for AUC and C-score, and Fisher test for Brier Scores
#' Models are given with the results of the survval function applied
#' to several cross-validation or bootstrap iterations, one per row
#' It is assumed that the same training data is used for Model 1 and Model 2
#' in each row
#'
#' @param res1 data frame
#' @param res0 list of predictors to use
#' @param sample_n prediction time of interest. If NULL, 90% quantile of event times is used
#' @param param_n random seed for replication
#' @return matrix: mean, std and p-values for AUC, BS, BS_scaled, C_score, Calib_slope, Calib_alpha
difftest <- function(res1, res0, sample_n, param_n) {
  m <- apply(res1 - res0, FUN = mean, 2, na.rm = 1)
  std <- apply(res1 - res0, FUN = stats::sd, 2, na.rm = 1)
  tpval <-
    function(x) {
      return(stats::t.test(x, alternative = "greater")$p.value)
    }
  #Fisher test for diff in R2 of 2 models
  #https://sites.duke.edu/bossbackup/files/2013/02/FTestTutorial.pdf
  pv_bs <-
    1 - stats::pf(
      mean(res1$BS, na.rm = 1) / mean(res0$BS, na.rm = 1),
      sample_n - param_n,
      sample_n - param_n - 1,
      lower.tail = FALSE
    )
  pvalue <-  apply(res1 - res0, FUN = tpval, 2)
  res <- rbind(m, std, pvalue)
  res["pvalue", "BS"] = pv_bs
  return(
    res[,c("T","C_score","AUCROC","BS","BS_scaled",
           "Calib_slope","Calib_alpha")]
  )
}

