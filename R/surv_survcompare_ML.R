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
#' @param train_ml TRUE/FALSE for whether to train DeepHit on its own, apart from the CoxPH->DeepHit ensemble. Default is FALSE as there is not much information in DeespSurv itself compared to the ensembled version.
#' @param ml "DeepHit" or "SRF", "DeepHit" by default
#' @param max_grid_size max number of hyperperam grid searches (random)
#' @return outcome = list(data frame with performance results, fitted Cox models, fitted DeespSurv)
#' @examples
#' df <-simulate_nonlinear(100)
#'  tuningparams <- list("dropout" = 0.2,"learning_rate" = 0.1,
#'  "num_nodes" = c(32, 32),"batch_size" = 100,"epochs" = 10)
#' mysurvcomp <-
#'    survcompare_ml(
#'    df,
#'    names(df)[1:4],
#'    tuningparams = tuningparams,
#'    ml = "DeepSurv"
#'    )
#' summary(mysurvcomp)
#' @export
survcompare_ml <- function(df_train,
                           predict_factors,
                           fixed_time = NULL,
                           randomseed = NULL,
                           useCoxLasso = FALSE,
                           outer_cv = 3,
                           inner_cv = 3,
                           tuningparams = list(),
                           return_models = FALSE,
                           repeat_cv = 2,
                           train_ml = FALSE,
                           ml = "DeepHit",
                           max_grid_size = 10) {

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

  if (is.null(randomseed)) {
    randomseed <- round(stats::runif(1) * 1e9, 0) + 1
  }
  if (is.null(fixed_time)) {
    fixed_time <- quantile(df_train[df_train$event == 1, "time"], 0.9, na.rm = TRUE)
  }
  #SRF_ensemble or DeepHit_ensemble
  ensemble_name <- paste(ml, "ensemble", sep="_")

  # cross-validation
  cox_cv <- survcox_cv(
    df = df_train,
    predict.factors = predict_factors,
    fixed_time = fixed_time ,
    outer_cv = outer_cv,
    randomseed = randomseed,
    useCoxLasso = useCoxLasso,
    return_models = return_models,
    repeat_cv = repeat_cv
  )
  if (train_ml) {
    if (ml == "DeepHit"){
      ml_cv <- deephit_cv(
        df = df_train,
        predict.factors = predict_factors,
        fixed_time = fixed_time,
        outer_cv = outer_cv,
        inner_cv = inner_cv,
        randomseed = randomseed,
        return_models = return_models,
        tuningparams = tuningparams,
        repeat_cv = repeat_cv,
        max_grid_size = max_grid_size
      )} else if (ml=="DeepSurv") {
        ml_cv <- deepsurv_cv(
          df = df_train,
          predict.factors = predict_factors,
          fixed_time = fixed_time,
          outer_cv = outer_cv,
          inner_cv = inner_cv,
          randomseed = randomseed,
          return_models = return_models,
          tuningparams = tuningparams,
          repeat_cv = repeat_cv
          )
        }else{  # ml == "SRF"
          ml_cv <- survsrf_cv(
            df = df_train,
            predict.factors = predict_factors,
            fixed_time = fixed_time,
            outer_cv = outer_cv,
            inner_cv = inner_cv,
            randomseed = randomseed,
            return_models = return_models,
            tuningparams = tuningparams,
            repeat_cv = repeat_cv
          )
        }
    }
  if (ml == "DeepHit"){
  ens1_cv <- ens_deephit_cv(
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
  } else if (ml=="DeepSurv") {
    ens1_cv <- ens_deepsurv_cv(
      df = df_train,
      predict.factors = predict_factors,
      fixed_time = fixed_time,
      outer_cv = outer_cv,
      inner_cv = inner_cv,
      randomseed = randomseed,
      return_models = return_models,
      tuningparams = tuningparams,
      useCoxLasso = useCoxLasso,
      repeat_cv = repeat_cv
    )
  }else{ # ml == "SRF"
    ens1_cv <- survensemble_cv(
      df = df_train,
      predict.factors = predict_factors,
      fixed_time = fixed_time,
      outer_cv = outer_cv,
      inner_cv = inner_cv,
      randomseed = randomseed,
      return_models = return_models,
      tuningparams = tuningparams,
      useCoxLasso = useCoxLasso,
      repeat_cv = repeat_cv
    )
    }

  # gathering the output: test&train performance

  stats_ci <- function(x, col = "C_score") {
    temp <- x[, col]
    c(
      "mean" = mean(temp, na.rm = 1),
      "sd" = sd(temp, na.rm = 1),
      "95CILow" = unname(quantile(temp, 0.025, na.rm= TRUE)),
      "95CIHigh" = unname(quantile(temp, 0.975, na.rm = TRUE))
    )
  }
  # combined results, if ML was trained (train_ml = TRUE), and then if only the Ensemble (FALSE)
  if (train_ml) {
    modelnames <-
      c(ifelse(useCoxLasso, "CoxLasso", "CoxPH"),ensemble_name,ml)
    results_mean <-
      as.data.frame(rbind(cox_cv$testaverage,ens1_cv$testaverage,ml_cv$testaverage))
    results_mean_train <-
      as.data.frame(rbind(cox_cv$trainaverage,ens1_cv$trainaverage,ml_cv$trainaverage))
    results_median <-
      as.data.frame(rbind(cox_cv$testmedian,ens1_cv$testmedian,ml_cv$testmedian))
    results_mean$sec <-
      round(as.numeric(c(cox_cv$time, ens1_cv$time, ml_cv$time)), 2)
    results_mean_train$sec <- results_mean$sec
    results_median$sec <- results_mean$sec
    auc_c_stats <- as.data.frame(rbind(
      stats_ci(cox_cv$test,  "C_score"),
      stats_ci(ens1_cv$test, "C_score"),
      stats_ci(ml_cv$test, "C_score"),
      stats_ci(cox_cv$test,  "AUCROC"),
      stats_ci(ens1_cv$test, "AUCROC"),
      stats_ci(ml_cv$test, "AUCROC")
    ))
    if(repeat_cv == 1){auc_c_stats_pooled= auc_c_stats
    }else{
    auc_c_stats_pooled <- as.data.frame(rbind(
      stats_ci(cox_cv$test_pooled,  "C_score"),
      stats_ci(ens1_cv$test_pooled, "C_score"),
      stats_ci(ml_cv$test_pooled, "C_score"),
      stats_ci(cox_cv$test_pooled,  "AUCROC"),
      stats_ci(ens1_cv$test_pooled, "AUCROC"),
      stats_ci(ml_cv$test_pooled, "AUCROC")
    ))
    }
  } else{ #train_ml = FALSE
    modelnames <-
      c(ifelse(useCoxLasso, "CoxLasso", "CoxPH"), ensemble_name)
    results_mean <-
      as.data.frame(rbind(cox_cv$testaverage, ens1_cv$testaverage))
    results_mean_train <-
      as.data.frame(rbind(cox_cv$trainaverage, ens1_cv$trainaverage))
    results_mean$sec <-
      round(as.numeric(c(cox_cv$time, ens1_cv$time)), 2)
    results_median <-
      as.data.frame(rbind(cox_cv$testmedian, ens1_cv$testmedian))
    results_median$sec <-
      round(as.numeric(c(cox_cv$time, ens1_cv$time)), 2)
    results_mean_train$sec <- results_mean$sec
    auc_c_stats <- as.data.frame(rbind(
      stats_ci(cox_cv$test,  "C_score"),
      stats_ci(ens1_cv$test, "C_score"),
      stats_ci(cox_cv$test,  "AUCROC"),
      stats_ci(ens1_cv$test, "AUCROC")
    ))
    if(repeat_cv == 1){auc_c_stats_pooled= auc_c_stats
    }else{
      auc_c_stats_pooled=
        as.data.frame(rbind(
          stats_ci(cox_cv$test_pooled,  "C_score"),
          stats_ci(ens1_cv$test_pooled, "C_score"),
          stats_ci(cox_cv$test_pooled,  "AUCROC"),
          stats_ci(ens1_cv$test_pooled, "AUCROC")
        ))
    }
  }

  # results_mean and results_mean_train row and col names

  col_order <- c("T","C_score","AUCROC","BS","BS_scaled","Calib_slope","Calib_alpha","sec")
  row.names(results_mean_train) <- modelnames
  row.names(results_mean) <- modelnames
  row.names(results_median) <- modelnames
  row.names(auc_c_stats) <-c(paste("C_score", modelnames, sep = "_"),paste("AUCROC", modelnames, sep = "_"))
  row.names(auc_c_stats_pooled) <-c(paste("C_score", modelnames, sep = "_"),paste("AUCROC", modelnames, sep = "_"))

  results_mean <- results_mean[col_order]
  results_median <- results_median[col_order]
  results_mean_train <- results_mean_train[col_order]

  # testing outperformance of the DeespSurv ensemble over the Cox model
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

  results_median["Diff",]=results_median[2,] - results_median[1, ]
  results_median["pvalue", ] =results_mean["pvalue",]

  # __________________________________________________
  # output
  output <- list()
  output$results_mean <- results_mean
  output$results_median <- results_median
  output$results_mean_train <- results_mean_train
  output$return_models <- list("CoxPH" = cox_cv$tuned_cv_models,
                               "ml" = ifelse(train_ml, ml_cv$tuned_cv_models,NaN),
                               "ensemble" = ens1_cv$tuned_cv_models)
  output$bestparams <- list("ensemble" = ens1_cv$bestparams,
                            "ml" = ifelse(train_ml, ml_cv$bestparams, NaN))
  output$test <- list("CoxPH" = cox_cv$test,
                      "ensemble" = ens1_cv$test)
  output$train <- list("CoxPH" = cox_cv$train,
                       "ensemble" = ens1_cv$train)
  output$difftest <- t_coxph
  output$main_stats <- auc_c_stats
  output$main_stats_pooled <- auc_c_stats_pooled
  output$randomseed <- randomseed
  output$useCoxLasso <- useCoxLasso
  output$model_name_base <- cox_cv$model_name
  output$model_name <- ens1_cv$model_name
  output$cv <- c(outer_cv, inner_cv,repeat_cv)
  class(output) <- "survcompare"
  summary.survcompare(output)
  return(output)
}
