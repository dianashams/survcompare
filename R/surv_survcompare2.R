#' survcompare2() compares if model evaluated in the second object (cv2) outperforms the first one (cv1)
#' it assumes that the same CV folds were used in each CV (if the repeated CV),
#' i.e. the same randomseed were used while performing cv1 and cv2
#'
#' @param cv1 an object of type "survensemble_cv", an outcome of survcv1, survsrf_cv, deephit_cv and other "_cv" functions of 'survcompare' packge
#' @param cv2 an object of type "survensemble_cv", the outcome of survcv1, survsrf_cv, deephit_cv and other "_cv" functions of 'survcompare' packge
#' @return outcome = list(data frame with performance results, fitted Cox models, fitted DeespSurv)
#' @examples
#' df <-simulate_nonlinear(100)
#' params <- names(df)[1:4]
#' cv1 <- survcv1(df, params, randomseed = 42)
#' cv1 <- survsrf_cv(df, params, randomseed = 42)
#' survcompare2(cv1, cv2)
#' @export
survcompare2 <- function(cv1, cv2) {
  name1 <- cv1$model_name
  name2 <- cv2$model_name

  Call <- match.call()
  inputs <- list(cv1, cv2)
  inputclass <-
    list(cv1 = "survensemble_cv", cv2 = "survensemble_cv")
  cp <- check_call(inputs, inputclass, Call)

  if (cp$anyerror)
    stop (paste(cp$msg[cp$msg != ""], sep = ""))

  if (sum(cv1$cv != cv2$cv) >0)
    stop (paste(
      "Error: cross-validation parameters differ (cv1$cv != cv2$cv)."
    ))
  if (cv1$randomseed != cv2$randomseed)
    stop (paste(
      "Error: different randomseeds (cv1$randomseed!= cv2$randomseed)."
    ))

  # gathering the output: test&train performance
  stats_ci <- function(x, col = "C_score") {
    temp <- x[, col]
    c(
      "mean" = mean(temp, na.rm = TRUE),
      "sd" = sd(temp, na.rm = TRUE),
      "95CILow" = unname(quantile(temp, 0.025,na.rm = TRUE)),
      "95CIHigh" = unname(quantile(temp, 0.975, na.rm=TRUE))
    )
  }

  repeat_cv <- cv1$cv[1]
  # combined results, if ML was trained (train_ml = TRUE), and then if only the Ensemble (FALSE)

  modelnames <- c(name1, name2)
  results_mean <-
    as.data.frame(rbind(cv1$testaverage, cv2$testaverage))
  results_mean_train <-
    as.data.frame(rbind(cv1$trainaverage, cv2$trainaverage))
  results_mean$sec <- round(as.numeric(c(cv1$time, cv2$time)), 2)
  results_median <-
    as.data.frame(rbind(cv1$testmedian, cv2$testmedian))
  results_median$sec <- round(as.numeric(c(cv1$time, cv2$time)), 2)
  results_mean_train$sec <- results_mean$sec
  auc_c_stats <- as.data.frame(rbind(
    stats_ci(cv1$test,  "C_score"),
    stats_ci(cv2$test, "C_score"),
    stats_ci(cv1$test,  "AUCROC"),
    stats_ci(cv2$test, "AUCROC")
  ))
  if (repeat_cv == 1) {
    auc_c_stats_pooled = auc_c_stats
  } else{
    auc_c_stats_pooled =
      as.data.frame(rbind(
        stats_ci(cv1$test_pooled,  "C_score"),
        stats_ci(cv2$test_pooled, "C_score"),
        stats_ci(cv1$test_pooled,  "AUCROC"),
        stats_ci(cv2$test_pooled, "AUCROC")
      ))
  }
  # results_mean and results_mean_train row and col names
  col_order <-  c("T", "C_score", "AUCROC", "BS", "BS_scaled",
                  "Calib_slope", "Calib_alpha", "sec")
  row.names(results_mean_train) <- modelnames
  row.names(results_mean) <- modelnames
  row.names(results_median) <- modelnames
  row.names(auc_c_stats) <-
    c(paste("C_score", modelnames, sep = "_"),
      paste("AUCROC", modelnames, sep = "_"))
  row.names(auc_c_stats_pooled) <-
    c(paste("C_score", modelnames, sep = "_"),
      paste("AUCROC", modelnames, sep = "_"))

  results_mean <- results_mean[col_order]
  results_median <- results_median[col_order]
  results_mean_train <- results_mean_train[col_order]

  # testing outperformance
  t_coxph <- difftest(cv2$test, cv1$test, 1000, 25)
  t_coxph_train <- difftest(cv2$train, cv1$train, 1000, 25)

  # adding results line for the differences with Cox-PH
  results_mean["Diff", ] = results_mean[2, ] - results_mean[1,]
  results_mean_train["Diff", ] = results_mean[2, ] - results_mean[1,]
  results_mean["pvalue", ] = c(t_coxph[3,], NaN) #NaN for "sec" column
  results_mean_train["pvalue", ] = c(t_coxph_train[3,], NaN)

  results_median["Diff", ] = results_median[2, ] - results_median[1,]
  results_median["pvalue",] = results_mean["pvalue", ]

  # output
  output <- list()
  output$results_mean <- results_mean
  output$results_median <- results_median
  output$results_mean_train <- results_mean_train
  output$return_models <- list(name1 = NaN,
                               name2 = NaN)
  output$bestparams <- list(name1 = NaN,
                            name2 = cv2$bestparams)
  output$test <- list(name1 = cv1$test, name2 = cv2$test)
  output$train <- list(name1 = cv1$train, name2 = cv2$train)
  output$difftest <- t_coxph
  output$main_stats <- auc_c_stats
  output$main_stats_pooled <- auc_c_stats_pooled
  output$randomseed <- cv1$randomseed
  output$useCoxLasso <- FALSE
  output$model_name_base <- name1
  output$model_name <- name2
  output$cv <- cv1$cv
  class(output) <- "survcompare"
  #summary.survcompare(output)
  return(output)
}
