#' @description
#' Compares two cross-validated models using surv____cv functions of this package.
#' Usually, this is a Cox Proportionate Hazards Model (or Cox LASSO), and Survival Random Forest,
#' or DeepHit (if installed from GitHub, not in CRAN version). Please see examples below.
#'
#' The same random seed, number of repetitions (repeat_cv), outer and inner folds numbers should have
#' been used for cross-validation objects such that those can be compared using survcompare2().
#' This ensures that the same data splits were used, and hence model performance on the same train/test splits are compared.
#' Harrel's c-index,time-dependent AUC-ROC, time-dependent Brier Score, and calibration slopes are reported.
#' The statistical significance of the performance differences is performed based on the C-index.
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
#' @param base an object of type "survensemble_cv", for example, outcomes of survcox_cv, survsrf_cv, survsrfens_cv, survsrfstack_cv
#' @param alternative an object of type "survensemble_cv", to compare to "base"
#' @return outcome = list(data frame with performance results, fitted Cox models, fitted DeespSurv)
#' @examples
#' df <-simulate_nonlinear(100)
#' params <- names(df)[1:4]
#' cv1 <- survcox_cv(df, params, randomseed = 42)
#' cv2 <- survsrf_cv(df, params, randomseed = 42)
#' survcompare2(cv1, cv2)
#' @export
survcompare2 <- function(base, alternative) {

  Call <- match.call()
  name1 <- base$model_name
  name2 <- alternative$model_name
  modelnames <- c(name1, name2)
  inputs <- list(base, alternative)
  inputclass <-
    list(base = "survensemble_cv", alternative = "survensemble_cv")
  cp <- check_call(inputs, inputclass, Call)

  if (cp$anyerror)
    stop (paste("Error: Objects class is not 'survensemble_cv'"))
  if (sum(base$cv[1:3] != alternative$cv[1:3]) >0)
    stop (paste(
      "Error: cross-validation parameters differ (base$cv != alternative$cv)."
    ))
  if (base$randomseed != alternative$randomseed)
    stop (paste(
      "Error: different randomseeds (base$randomseed!= alternative$randomseed)."
    ))

  results_mean <-
    as.data.frame(rbind(base$testaverage, alternative$testaverage))
  results_mean_train <-
    as.data.frame(rbind(base$trainaverage, alternative$trainaverage))
  results_mean$sec <- round(as.numeric(c(base$time, alternative$time)), 2)
  results_median <-
    as.data.frame(rbind(base$testmedian, alternative$testmedian))
  results_median$sec <- round(as.numeric(c(base$time, alternative$time)), 2)
  results_mean_train$sec <- results_mean$sec
  stats_ci <- function(x, col = "C_score") {
    temp <- x[, col]
    c(
      "mean" = mean(temp, na.rm = TRUE),
      "sd" = sd(temp, na.rm = TRUE),
      "95CILow" = unname(quantile(temp, 0.025,na.rm = TRUE)),
      "95CIHigh" = unname(quantile(temp, 0.975, na.rm=TRUE))
    )
  }
  auc_c_stats <- as.data.frame(rbind(
    stats_ci(base$test,  "C_score"),
    stats_ci(alternative$test, "C_score"),
    stats_ci(base$test,  "AUCROC"),
    stats_ci(alternative$test, "AUCROC")
  ))

  repeat_cv <- base$cv[3]
  if (repeat_cv == 1) {
    auc_c_stats_pooled = auc_c_stats
  } else{
    auc_c_stats_pooled =
      as.data.frame(rbind(
        stats_ci(base$test_pooled,  "C_score"),
        stats_ci(alternative$test_pooled, "C_score"),
        stats_ci(base$test_pooled,  "AUCROC"),
        stats_ci(alternative$test_pooled, "AUCROC")
      ))
  }
  # results_mean and results_mean_train row and col names
  row.names(results_mean_train) <- modelnames
  row.names(results_mean) <- modelnames
  row.names(results_median) <- modelnames
  row.names(auc_c_stats) <-
    c(paste("C_score", modelnames, sep = "_"),
      paste("AUCROC", modelnames, sep = "_"))
  row.names(auc_c_stats_pooled) <-
    c(paste("C_score", modelnames, sep = "_"),
      paste("AUCROC", modelnames, sep = "_"))

  col_order <-  c("T", "C_score", "AUCROC", "BS", "BS_scaled",
                  "Calib_slope", "Calib_alpha", "sec")
  results_mean <- results_mean[col_order]
  results_median <- results_median[col_order]
  results_mean_train <- results_mean_train[col_order]

  # testing outperformance
  t_coxph <- difftest(alternative$test, base$test, 1000, 25)
  t_coxph_train <- difftest(alternative$train, base$train, 1000, 25)

  # adding results' row for the differences with Cox-PH
  results_mean["Diff", ] = results_mean[2, ] - results_mean[1,]
  results_mean_train["Diff", ] = results_mean[2, ] - results_mean[1,]
  results_mean["pvalue", ] = c(t_coxph[3,], NaN) #NaN for "sec" column
  results_mean_train["pvalue", ] = c(t_coxph_train[3,], NaN)

  results_median["Diff", ] = results_median[2, ] - results_median[1,]
  # we test significant difference in means, so p-value is from the "means" table
  results_median["pvalue",] = results_mean["pvalue", ]

  # output
  output <- list()
  output$results_mean <- results_mean
  output$results_median <- results_median
  output$results_mean_train <- results_mean_train
  output$return_models <- list(name1 = NaN,name2 = NaN)
  output$bestparams <- list(name1 = NaN,name2 = alternative$bestparams)
  output$test <- list(name1 = base$test, name2 = alternative$test)
  output$train <- list(name1 = base$train, name2 = alternative$train)
  output$difftest <- t_coxph
  output$main_stats <- auc_c_stats
  output$main_stats_pooled <- auc_c_stats_pooled
  output$randomseed <- base$randomseed
  output$useCoxLasso <- FALSE
  output$model_name_base <- name1
  output$model_name <- name2
  output$cv <- base$cv

  class(output) <- "survcompare"
  return(output)
}

# Testing statistical significance of the ensembled model outperformance
# over the baseline model (Model 1 vs Model 0)
difftest <- function(res1, res0, sample_n, param_n) {
  m <- apply(res1 - res0, FUN = mean, 2, na.rm = 1)
  std <- apply(res1 - res0, FUN = stats::sd, 2, na.rm = 1)
  tpval <-
    function(x) {
      if (class(try(stats::t.test(x, alternative = "greater")$p.value)
      )  ==  "try-error")
        return(NaN)
      return(stats::t.test(x, alternative = "greater")$p.value)
    }
  #Fisher test for diff in R2 of 2 models
  # https://sites.duke.edu/bossbackup/files/2013/02/FTestTutorial.pdf
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
