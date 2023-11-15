surv_CV <-
  function(df,
           predict.factors,
           fixed_time = NaN,
           outer_cv = 3,
           inner_cv = 3,
           repeat_cv = 2,
           randomseed = NULL,
           return_models = FALSE,
           train_function,
           predict_function,
           model_args = list(),
           predict_args = list(),
           model_name = "a model") {
    time_0 <- Sys.time()
    if (is.null(randomseed)) {
      randomseed <- round(stats::runif(1) * 1e9, 0)
    }

    if (!any(
      is.data.frame(df),
      predict.factors %in% colnames(df),
      c("time", "event") %in% colnames(df)
    )) {
      stop("Data should be a data frame, predictors should correspond to the columns.")
    }

    if (sum(is.na(df[c("time", "event", predict.factors)])) > 0) {
      stop("Missing data can not be handled. Please impute first.")
    }

    if (sum(is.nan(fixed_time)) > 0) {
      fixed_time <- round(quantile(df[df$event == 1, "time"], 0.9), 1)
    }
    predict.factors <- eligible_params(predict.factors, df)

    if (length(predict.factors) == 0) {
      print("No eligible params")
      return(NULL)
    }
    Call <- match.call()

    #defining number of repeated cv
    if (is.null(repeat_cv)) {
      repeat_cv = 1
    }
    if (is.numeric(repeat_cv) & repeat_cv > 1) {
      repeat_cv = round(repeat_cv, 0)
    } else{
      repeat_cv = 1
    }

    note_srf = ifelse(
      model_name == "Survival Random Forest",
      "",
      "For SRF inner CV is not used if oob = TRUE (default)."
    )

    print(
      paste(
        "Cross-validating ",
        model_name,
        " using ",
        repeat_cv,
        " repeat(s), ",
        outer_cv,
        " outer, ",
        inner_cv,
        " inner loops).",
        note_srf,
        sep = ""
      )
    )

    modelstats_train <- list()
    modelstats_test <- list()
    models_for_each_cv <- list()
    #progress bar
    pb <- utils::txtProgressBar(0, outer_cv * repeat_cv, style = 3)
    utils::setTxtProgressBar(pb, outer_cv * repeat_cv / 50)

    # cross-validation loop
    for (rep_cv in 1:repeat_cv) {
      set.seed(randomseed + rep_cv)
      if (rep_cv != 1) {
        df <- df[sample(1:nrow(df)), ]
      }
      cv_folds <-
        caret::createFolds(df$event, k = outer_cv, list = FALSE)

      for (cv_iteration in 1:outer_cv) {
        df_train_cv <- df[cv_folds != cv_iteration,]
        df_test_cv <- df[cv_folds == cv_iteration,]
        # tune the model using train_function
        trained_model <-
          do.call(train_function,
                  append(list(df_train_cv, predict.factors), model_args))
        # compute event probability predictions for given times
        y_predict_test <-
          do.call(predict_function,
                  append(
                    list(trained_model, df_test_cv, fixed_time),
                    predict_args
                  ))
        y_predict_train <-
          do.call(predict_function,
                  append(
                    list(trained_model, df_train_cv, fixed_time),
                    predict_args
                  ))
        modelstats_test[[cv_iteration + (rep_cv - 1) * outer_cv]] <-
          surv_validate(y_predict_test,
                        fixed_time,
                        df_train_cv,
                        df_test_cv,
                        weighted = 1)
        modelstats_train[[cv_iteration + (rep_cv - 1) * outer_cv]] <-
          surv_validate(y_predict_train,
                        fixed_time,
                        df_train_cv,
                        df_train_cv,
                        weighted = 1)
        if (return_models) {
          models_for_each_cv[[cv_iteration + (rep_cv - 1) * outer_cv]] <-
            trained_model
        }
        utils::setTxtProgressBar(pb, cv_iteration + (rep_cv - 1) * outer_cv)
      }
    }

    df_modelstats_test <- data.frame(modelstats_test[[1]])
    df_modelstats_train <- data.frame(modelstats_train[[1]])

    for (i in 2:(outer_cv * repeat_cv)) {
      df_modelstats_test[i, ] <- modelstats_test[[i]]
      df_modelstats_train[i, ] <- modelstats_train[[i]]
    }
    row.names(df_modelstats_train) <- 1:(outer_cv * repeat_cv)
    row.names(df_modelstats_test) <- 1:(outer_cv * repeat_cv)

    utils::setTxtProgressBar(pb, outer_cv * repeat_cv)
    close(pb)
    df_modelstats_test$test <- 1
    df_modelstats_train$test <- 0

    #summary for printing and summary(obj)
    stats_summary <- function(x) {
      #remove the last column ("test" 1 or 0)
      x = x[, 1:(dim(x)[2] - 1)]
      as.data.frame(round(
        cbind(
          "mean" = apply(x, 2, mean, na.rm = 1),
          "sd" = apply(x, 2, sd, na.rm = 1),
          "95CILow" = apply(x, 2, quantile, 0.025, na.rm = 1),
          "95CIHigh" = apply(x, 2, quantile, 0.975, na.rm = 1)
        ),
        4
      ))
    }
    #mean, sd and confidence intervals for all test and train datasets
    summarydf <- as.data.frame(cbind(
      "test" = stats_summary(df_modelstats_test),
      "train" = stats_summary(df_modelstats_train)
    ))

    output <- list()
    output$test <- df_modelstats_test
    output$train <- df_modelstats_train
    output$testaverage <-
      sapply(df_modelstats_test, mean, na.rm = 1)
    output$trainaverage <-
      sapply(df_modelstats_train, mean, na.rm = 1)
    output$tuned_cv_models <- models_for_each_cv
    output$randomseed <- randomseed
    output$call <- Call
    output$cv <- c(outer_cv, inner_cv, repeat_cv)
    output$summarydf <- summarydf
    class(output) <- "survensemble_cv"
    time_1 <- Sys.time()
    output$time <- time_1 - time_0
    print(time_1 - time_0)
    return(output)
  }
