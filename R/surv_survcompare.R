#' Cross-validates and compares Cox Proportionate Hazards and Survival Random Forest models
#' @description
#' The function performs a repeated nested cross-validation for
#' 1) Cox-PH (survival package, survival::coxph) or Cox-Lasso (glmnet package, glmnet::cox.fit)
#' 2) Ensemble of the Cox model and Survival Random Forest (randomForestSRC::rfsrc)
#' 3) Survival Random Forest on its own, if train_srf = TRUE
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
#' @importFrom randomForestSRC rfsrc
#' @importFrom survival coxph
#' @param df_train training data, a data frame with "time" and "event" columns to define the survival outcome
#' @param predict_factors list of column names to be used as predictors
#' @param predict_time prediction time of interest. If NULL, 0.90th quantile of event times is used
#' @param randomseed random seed for replication
#' @param useCoxLasso TRUE / FALSE, for whether to use regularized version of the Cox model, FALSE is default
#' @param outer_cv k in k-fold CV
#' @param inner_cv k in k-fold CV for internal CV to tune survival random forest hyper-parameters
#' @param srf_tuning list of tuning parameters for random forest: 1) NULL for using a default tuning grid, or 2) a list("mtry"=c(...), "nodedepth" = c(...), "nodesize" = c(...))
#' @param return_models TRUE/FALSE to return the trained models; default is FALSE, only performance is returned
#' @param repeat_cv if NULL, runs once, otherwise repeats several times with different random split for CV, reports average of all
#' @param train_srf TRUE/FALSE for whether to train SRF on its own, apart from the CoxPH->SRF ensemble. Default is FALSE as there is not much information in SRF itself compared to the ensembled version.
#' @return outcome = list(data frame with performance results, fitted Cox models, fitted SRF)
#' @examples
#' \dontshow{rfcores_old <- options()$rf.cores; options(rf.cores=1)}
#' df <-simulate_nonlinear(100)
#' srf_params <- list("mtry" = c(2), "nodedepth"=c(25), "nodesize" =c(15))
#' mysurvcomp <- survcompare(df, names(df)[1:4], srf_tuning = srf_params, outer_cv = 2, inner_cv =2)
#' summary(mysurvcomp)
#' \dontshow{options(rf.cores=rfcores_old)}
#' @export
survcompare <- function(df_train,
                        predict_factors,
                        predict_time = NULL,
                        randomseed = NULL,
                        useCoxLasso = FALSE,
                        outer_cv = 3,
                        inner_cv = 3,
                        srf_tuning = list(),
                        return_models = FALSE,
                        repeat_cv = 2,
                        train_srf = FALSE) {

  Call <- match.call()
  inputs <- list(
    df_train,
    predict_factors,
    predict_time,
    outer_cv,
    inner_cv,
    repeat_cv,
    randomseed,
    return_models,
    srf_tuning,
    useCoxLasso
  )

  inputclass <-
    list(
      df_train = "data.frame",
      predict_factors = "character",
      predict_time = "numeric",
      outer_cv = "numeric",
      inner_cv = "numeric",
      repeat_cv = "numeric",
      randomseed = "numeric",
      return_models = "logical",
      srf_tuning = "list",
      useCoxLasso = "logical"
    )
  cp <- check_call(inputs, inputclass, Call)
  if (cp$anyerror)
    stop (paste(cp$msg[cp$msg != ""], sep = ""))

  if (sum(is.na(df_train[c("time", "event", predict_factors)])) > 0) {
    stop("Missing data can not be handled. Please impute first.")
  }

  if (is.null(randomseed)) {
    randomseed <- round(stats::runif(1) * 1e9, 0) + 1
  }
  if (is.null(predict_time)) {
    predict_time <- quantile(df_train[df_train$event == 1, "time"], 0.9)
  }

  # cross-validation
  cox_cv <- survcox_cv(
    df = df_train,
    predict.factors = predict_factors,
    fixed_time = predict_time ,
    outer_cv = outer_cv,
    randomseed = randomseed,
    useCoxLasso = useCoxLasso,
    return_models = return_models,
    repeat_cv = repeat_cv
  )
  if (train_srf) {
    srf_cv <- survsrf_cv(
      df = df_train,
      predict.factors = predict_factors,
      fixed_time = predict_time,
      outer_cv = outer_cv,
      inner_cv = inner_cv,
      randomseed = randomseed,
      return_models = return_models,
      srf_tuning = srf_tuning,
      repeat_cv = repeat_cv
    )
  }
  ens1_cv <- survensemble_cv(
    df = df_train,
    predict.factors = predict_factors,
    fixed_time = predict_time,
    outer_cv = outer_cv,
    inner_cv = inner_cv,
    randomseed = randomseed,
    return_models = return_models,
    srf_tuning = srf_tuning,
    useCoxLasso = useCoxLasso,
    repeat_cv = repeat_cv
  )

  # gathering the output: test&train performance

  stats_ci <- function(x, col = "C_score") {
    temp <- x[, col]
    c(
      "mean" = mean(temp, na.rm = 1),
      "sd" = sd(temp, na.rm = 1),
      "95CILow" = unname(quantile(temp, 0.025)),
      "95CIHigh" = unname(quantile(temp, 0.975))
    )
  }

  if (train_srf) {
    modelnames <-
      c(ifelse(useCoxLasso, "CoxLasso", "CoxPH"),
        "SRF_Ensemble",
        "SRF")
    results_mean <-
      as.data.frame(rbind(
        cox_cv$testaverage,
        ens1_cv$testaverage,
        srf_cv$testaverage
      ))
    results_mean_train <-
      as.data.frame(rbind(
        cox_cv$trainaverage,
        ens1_cv$trainaverage,
        srf_cv$trainaverage
      ))
    results_mean$sec <-
      round(as.numeric(c(
        cox_cv$time, ens1_cv$time, srf_cv$time
      )), 2)
    results_mean_train$sec <- results_mean$sec
    auc_c_stats <- as.data.frame(rbind(
      stats_ci(cox_cv$test,  "C_score"),
      stats_ci(ens1_cv$test, "C_score"),
      stats_ci(srf_cv$test, "C_score"),
      stats_ci(cox_cv$test,  "AUCROC"),
      stats_ci(ens1_cv$test, "AUCROC"),
      stats_ci(srf_cv$test, "AUCROC"),
    ))
  } else{
    modelnames <-
      c(ifelse(useCoxLasso, "CoxLasso", "CoxPH"), "SRF_Ensemble")
    results_mean <-
      as.data.frame(rbind(cox_cv$testaverage, ens1_cv$testaverage))
    results_mean_train <-
      as.data.frame(rbind(cox_cv$trainaverage, ens1_cv$trainaverage))
    results_mean$sec <-
      round(as.numeric(c(cox_cv$time, ens1_cv$time)), 2)
    results_mean_train$sec <- results_mean$sec
    auc_c_stats <- as.data.frame(rbind(
      stats_ci(cox_cv$test,  "C_score"),
      stats_ci(ens1_cv$test, "C_score"),
      stats_ci(cox_cv$test,  "AUCROC"),
      stats_ci(ens1_cv$test, "AUCROC")
    ))
  }

  col_order <- c("T",
                 "C_score",
                 "AUCROC",
                 "BS",
                 "BS_scaled",
                 "Calib_slope",
                 "Calib_alpha",
                 "sec")
  row.names(results_mean_train) <- modelnames
  row.names(results_mean) <- modelnames
  row.names(auc_c_stats) <-
    c(paste("C_score", modelnames, sep = "_"),
      paste("AUCROC", modelnames, sep = "_"))
  results_mean <- results_mean[col_order]
  results_mean_train <- results_mean_train[col_order]

  # testing outperformance of the SRF ensemble over the Cox model
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
  results_mean["Diff",] = results_mean[2,] - results_mean[1, ]
  results_mean_train["Diff",] = results_mean[2,] - results_mean[1, ]
  results_mean["pvalue",] = c(t_coxph[3, ], NaN) #NaN for "sec" column
  results_mean_train["pvalue",] = c(t_coxph_train[3, ], NaN)
  # __________________________________________________

  # output
  output <- list()
  output$results_mean <- results_mean
  output$results_mean_train <- results_mean_train
  output$return_models <- list("CoxPH" = cox_cv$tuned_cv_models,
                               "SRF_ensemble" = ens1_cv$tuned_cv_models)
  output$test <- list("CoxPH" = cox_cv$test,
                      "SRF_ensemble" = ens1_cv$test)
  output$train <- list("CoxPH" = cox_cv$train,
                       "SRF_ensemble" = ens1_cv$train)
  output$difftest <- t_coxph
  output$main_stats <- auc_c_stats
  output$randomseed <- randomseed
  output$useCoxLasso <- useCoxLasso
  class(output) <- "survcompare"
  summary.survcompare(output)
  return(output)
}


# Testing statistical significance of the ensembled model outperformance
# over the baseline model (Model 1 vs Model 0)
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
  return(res[, c("T",
                 "C_score",
                 "AUCROC",
                 "BS",
                 "BS_scaled",
                 "Calib_slope",
                 "Calib_alpha")])
}
