############# SRF ################

#' Internal function to compute survival probability
#'  by time from a fitted survival random forest
#'
#' @param rfmodel pre-trained survrf_train model
#' @param df_to_predict test data
#' @param fixed_time  at which event probabilities are computed
#' @param oob TRUE/FALSE use out-of-bag prediction
#' @examples
#' df <- simsurv_nonlinear()
#' #params<- c("age", "hyp", "bmi")
#' #s <- survrf_train(df, params)
#' #p <- survrf_predict(s, df, 5)
#' @return output list: output$train, test, testaverage, traintaverage, time
srf_survival_prob_for_time <-
  function(rfmodel,
           df_to_predict,
           fixed_time,
           oob = FALSE) {
    # finds survival prediction from srf model
    if (oob) {
      predicted_matrix <- predict(
        rfmodel,
        newdata = df_to_predict,
        ensemble = "oob",
        outcome = "test"
      )
    } else {
      predicted_matrix <- predict(rfmodel, newdata = df_to_predict)
    }

    j_for_fixedtime <-
      match(1,
            round(predicted_matrix$time.interest, 1) == fixed_time,
            nomatch = -100)

    if (j_for_fixedtime == -100) {
      # print("no fixed time match was found, using closest")
      j_for_fixedtime <-
        which.min(abs(predicted_matrix$time.interest - fixed_time))
    }

    if (oob) {
      y_predicted <- predicted_matrix$survival.oob[, j_for_fixedtime]
    } else {
      y_predicted <- predicted_matrix$survival[, j_for_fixedtime]
    }
    return(y_predicted)
  }

#' Internal function to tune SRF model, in nested CV loop
#'
#' @param df_tune data frame
#' @param predict.factors predictor names
#' @param cv_number k in k-fold CV
#' @param fixed_time  NaN
#' @param randomseed  random seed
#' @param mtry  c(3,4,5) tuning parameter
#' @param nodesize  c(10,20,50) tuning parameter
#' @param nodedepth  100 tuning paramreter
#' @param verbose  FALSE
#' @param oob  TRUE
#' @param nodesize  at which event probabilities are computed
#' @param oob TRUE/FALSE use out-of-bag prediction
#' @examples
#' df <- simsurv_nonlinear()
#' @return  output=list(modelstats, bestbrier, bestauc, bestcindex)
srf_tune <- function(df_tune,
                     predict.factors,
                     cv_number = 3,
                     fixed_time = NaN,
                     randomseed = NULL,
                     mtry = c(3, 4, 5),
                     nodesize = c(10, 20, 50),
                     nodedepth = c(100),
                     verbose = FALSE,
                     oob = TRUE) {
  # function to tune survival random forest by mtry, nodesize and nodedepth grid
  # if oob = TRUE, there is no CV ! as OOB does the job already

  # take out the factors which are not in df_tune or the ones which take 1 value
  predict.factors <- eligible_params(predict.factors, df_tune)

  # set seed
  if (is.null(randomseed)) {randomseed <- round(stats::runif(1)*1e9,0)}
  set.seed(randomseed)

  nodesize <- nodesize[nodesize <= dim(df_tune)[1] / 6]
  # limit mtry with the number of predictors and nodesize by 1/6 of sample size
  if (sum(nodesize > dim(df_tune)[1] / 6) > 0) {
    if (verbose) {
      print("Warning - some min nodesize is > 1/6 of the sample size (1/2 of CV fold)")
    }
  }

  # if all numbers higher than number of factors, only check this factor
  if (sum(mtry > length(predict.factors)) == length(mtry)) {
    mtry <- c(length(predict.factors))
  }
  mtry <- mtry[mtry <= length(predict.factors)]

  # grid of values to tune
  grid_of_values <- expand.grid(
    "mtry" = mtry,
    "nodesize" = nodesize,
    "nodedepth" = nodedepth
  )
  if (verbose) {
    print(paste("Grid size is", dim(grid_of_values)[1]))
  }

  if (dim(grid_of_values)[1] == 0) {
    output <- list()
    return(output)
  }

  # defining output for fixed_time
  if (sum(is.nan(fixed_time) > 0) |
      length(fixed_time) > 1) {
    # not implemented for multiple time tuning yet
    fixed_time <-
      round(quantile(df_tune[df_tune$event == 1, "time"], 0.9), 1)
  }

  df_tune$time_f <-
    ifelse(df_tune$time <= fixed_time, df_tune$time, fixed_time)
  df_tune$event_f <-
    ifelse(df_tune$event == 1 & df_tune$time <= fixed_time, 1, 0)

  # going through combinations
  modelstats <- list()

  if (oob == FALSE) {
    # we do CV instead of using OOB predictions to tune SRF
    set.seed(randomseed)
    cv_folds <-
      caret::createFolds(df_tune$event, k = cv_number, list = FALSE)
    for (i in 1:dim(grid_of_values)[1]) {
      if (verbose) {
        print(grid_of_values[i, ])
      }
      mtry_i <- grid_of_values[i, "mtry"]
      nodesize_i <- grid_of_values[i, "nodesize"]
      nodedepth_i <- grid_of_values[i, "nodedepth"]

      # train grid combination for each cv_iteration
      modelstats_cv <- list()
      for (cv_iteration in 1:cv_number) {
        print(paste("SRF tuning CV step", cv_iteration, "/out of", cv_number))
        df_train_cv <- df_tune[cv_folds != cv_iteration, ]
        df_test_cv <- df_tune[cv_folds == cv_iteration, ]
        # train SRF
        rf.dt <- randomForestSRC::rfsrc(
          as.formula(paste(
            "Surv(time, event) ~",
            paste(predict.factors, collapse = "+")
          )),
          data = df_train_cv,
          nodesize = nodesize_i,
          # this is AVERAGE size, so we want this to be quite high
          ntree = 300,
          mtry = mtry_i,
          nodedepth = nodedepth_i,
          nsplit = 50,
          splitrule = "logrank",
          statistics = FALSE,
          membership = TRUE,
          importance = "none",
          # to speed up by switching off VIMP calculations
          seed = randomseed
        )
        # compute predicted event probability for all people
        y_predicted <-
          1 - srf_survival_prob_for_time(rf.dt, df_test_cv, fixed_time, oob = FALSE)
        validation_stats <-
          survval(y_predicted,
                              fixed_time,
                              df_train_cv,
                              df_test_cv,
                              weighted = TRUE
          )
        modelstats_cv[[cv_iteration]] <- c(
          "mtry" = mtry_i,
          "nodesize" = nodesize_i,
          "nodedepth" = nodedepth_i,
          "time" = validation_stats$T,
          "AUCROC" = validation_stats$AUCROC,
          "BS" = validation_stats$BS,
          "BS_scaled" = validation_stats$BS_scaled,
          "C_score" = validation_stats$C_score,
          "Calib_alpha" = validation_stats$Calib_alpha,
          "Calib_slope" = validation_stats$Calib_slope
        )
      } # end k-fold CV for one grid combination

      # averaging over cv-steps, firs transform to data.frame to use mean()
      modelstats_cv_df <- data.frame(t(modelstats_cv[[1]]))
      for (j in 2:cv_number) {
        modelstats_cv_df <- rbind(modelstats_cv_df, t(modelstats_cv[[j]]))
      }
      modelstats[[i]] <-
        c(
          modelstats_cv[[1]]["mtry"],
          modelstats_cv[[1]]["nodesize"],
          modelstats_cv[[1]]["nodedepth"],
          "AUCROC" = mean(modelstats_cv_df$AUCROC, na.rm = 1),
          "BS" = mean(modelstats_cv_df$BS, na.rm = 1),
          "BS_scaled" = mean(modelstats_cv_df$BS_scaled, na.rm = 1),
          "C_score" = mean(modelstats_cv_df$C_score, na.rm = 1),
          "Calib_alpha" = mean(modelstats_cv_df$Calib_alpha, na.rm = 1),
          "Calib_slope" = mean(modelstats_cv_df$Calib_slope, na.rm = 1),
          "time" = fixed_time
        )
    } # end for grid
  } else {
    # end if(oob==false)
    if (verbose) {
      print(
        "No internal CV for training SRF,
           instead out-of-bag predictions used to assess performance"
      )
    }
    for (i in 1:dim(grid_of_values)[1]) {
      if (verbose) {
        print(grid_of_values[i, ])
      }
      mtry_i <- grid_of_values[i, "mtry"]
      nodesize_i <- grid_of_values[i, "nodesize"]
      nodedepth_i <- grid_of_values[i, "nodedepth"]

      rf.dt <- randomForestSRC::rfsrc(
        as.formula(paste(
          "Surv(time, event) ~",
          paste(predict.factors, collapse = "+")
        )),
        data = df_tune,
        nodesize = nodesize_i,
        # this is AVERAGE size, so we want this to be quite high
        ntree = 300,
        mtry = mtry_i,
        nodedepth = nodedepth_i,
        nsplit = 50,
        splitrule = "logrank",
        statistics = FALSE,
        membership = TRUE,
        importance = "none",
        # to speed up by switching off VIMP calculations
        seed = randomseed
      )

      # compute predicted event probability for all people
      y_predicted <-
        1 - srf_survival_prob_for_time(rf.dt, df_tune, fixed_time, oob = TRUE)
      validation_stats <-
        survval(y_predicted, fixed_time, df_tune, df_tune, weighted = TRUE)
      modelstats[[i]] <- c(
        "mtry" = mtry_i,
        "nodesize" = nodesize_i,
        "nodedepth" = nodedepth_i,
        "time" = validation_stats$T,
        "AUCROC" = validation_stats$AUCROC,
        "BS" = validation_stats$BS,
        "BS_scaled" = validation_stats$BS_scaled,
        "C_score" = validation_stats$C_score,
        "Calib_alpha" = validation_stats$Calib_alpha,
        "Calib_slope" = validation_stats$Calib_slope
      )
    } # end for (i in grid)
  } # end else

  # reshaping into data frame
  df_modelstats <- data.frame("V1" = modelstats[[1]])
  # check if there was more than 1 grid search
  if (dim(grid_of_values)[1] > 1) {
    for (i in 2:dim(grid_of_values)[1]) {
      df_modelstats[i] <- modelstats[[i]]
    }
  }
  df_modelstats <- data.frame(t(df_modelstats))

  if (verbose == TRUE) {

    print(paste("AUC varies from", round(min(
      df_modelstats$AUCROC
    ), 4), "to", round(max(
      df_modelstats$AUCROC
    ), 4)))
    print(paste(
      "Brier score varies from",
      round(min(df_modelstats$BS), 4),
      "to",
      round(max(df_modelstats$BS), 4)
    ))
    print("Combination with highest AUC")
    print(df_modelstats[
      which.max(df_modelstats$AUCROC),
      c("mtry", "nodesize", "nodedepth")
    ])
    print("Combination with lowest Brier Score")
    print(df_modelstats[
      which.min(df_modelstats$BS),
      c("mtry", "nodesize", "nodedepth")
    ])
  }
  output <- list()
  output$modelstats <- df_modelstats
  output$bestbrier <- df_modelstats[which.min(df_modelstats$BS), ]
  output$bestauc <- df_modelstats[which.max(df_modelstats$AUCROC), ]
  output$bestcindex <-
    df_modelstats[which.max(df_modelstats$C_score), ]
  return(output)
}


#' Fits randomForestSRC, with tuning by mtry,nodedepth and nodesize
#' underlying model is by Ishwaran et al(2008)
#' https://www.randomforestsrc.org/articles/survival.html
#' Ishwaran H, Kogalur UB, Blackstone EH, Lauer MS. Random survival forests.
#' The Annals of Applied Statistics. 2008;2:841–60.
#'
#' @param df_train  data, "time" and "event" should describe survival outcome
#' @param predict.factors list of the column names to be used as predictors
#' @param fixed_time time at which performance is maximized
#' @param inner_cv k in k-fold CV for model tuning
#' @param randomseed random seed
#' @param srf_tuning list of mtry, nodedepth and nodesize, default is NULL
#' @param fast_version  TRUE/FALSE, TRUE by default
#' @param oob FALSE/TRUE, TRUE by default
#' @param verbose TRUE/FALSE, FALSE by default
#' @examples
#' df <- simsurv_nonlinear()
#' @return output = list(beststats, allstats, model)
survrf_train <- function(df_train,
                         predict.factors,
                         fixed_time = NaN,
                         inner_cv = 3,
                         randomseed = NULL,
                         srf_tuning = NULL,
                         fast_version = TRUE,
                         oob = TRUE,
                         verbose = FALSE) {
  # take out predictors not in df_train or constant
  predict.factors <- eligible_params(predict.factors, df_train)

  if (is.null(randomseed)) {
    randomseed <- round(stats::runif(1) * 1e9, 0)
  }

  stopifnot(expr = {
    is.data.frame(df_train)
    predict.factors %in% colnames(df_train)
  })

  # setting fixed_time if not given
  if (sum(is.nan(fixed_time)) > 0 | (length(fixed_time) > 1)) {
    fixed_time <-
      round(quantile(df_train[df_train$event == 1, "time"], 0.9), 1)
  }

  # for now only for best AUC but can be amended for brier score or cindex
  # defining the tuning grid for SRF
  p <- length(predict.factors) # number of predictors
  n <- dim(df_train)[1]
  mtry_default <- round(sqrt(p), 0)

  #set default mtry, nodesize and nodedepth
  if (p <= 10) {
    mtry <- c(2, 3, 4, 5)
  } else {
    if (p <= 25) {
      mtry <- c(3, 5, 7, 10, 15)
    } else {
      mtry <- round(c(p / 10, p / 5, p / 3, p / 2, mtry_default), 0)
    }
  }
  mtry <- mtry[mtry<=length(predict.factors)]

  nodesize <-
    seq(min(15, round(n / 6 - 1, 0)), max(min(n / 10, 50), 30), 5)

  nodedepth <- c(5, 25)

  # update them if given in srf_tuning
  if (is.list(srf_tuning)) {
    if (!is.null(srf_tuning$nodesize)) {
      nodesize = srf_tuning$nodesize
    }
    if (!is.null(srf_tuning$mtry)) {
      mtry = c(srf_tuning$mtry, mtry_default)
    }
    if (!is.null(srf_tuning$nodedepth)) {
      nodedepth = srf_tuning$nodedepth
    }
  }

  if (verbose) {
    print(expand.grid(
      "mtry" = mtry,
      "nodesize" = nodesize,
      "nodedepth" = nodedepth
    ))
  }

  if (fast_version == TRUE) {
    # fast version:
    # first, takes recommended mtry and optimizes by depth and node size,
    # second, optimizes by mtry

    tune1 <- srf_tune(
      df_tune = df_train,
      cv_number = inner_cv,
      predict.factors,
      fixed_time = fixed_time,
      randomseed = randomseed,
      mtry = mtry_default,
      nodesize = nodesize,
      nodedepth = nodedepth,
      oob = oob
    )

    nodesize_best <- as.integer(tune1$bestauc["nodesize"])
    nodedepth_best <- as.integer(tune1$bestauc["nodedepth"])
    # using the depth and size check the best mtry

    tune2 <- srf_tune(
      df_tune = df_train,
      cv_number = inner_cv,
      predict.factors,
      fixed_time = fixed_time,
      randomseed = randomseed,
      mtry = mtry,
      nodesize = nodesize_best,
      nodedepth = nodedepth_best,
      oob = oob
    )
    mtry_best <- tune2$bestauc["mtry"]
    best_combo_stat <- tune2$bestauc
    modelstatsall <- rbind(tune1$modelstats, tune2$modelstats)
  } else {
    # slower version:
    # goes through all the grid combinations mtry, nodesize, nodedepth

    print("tuneall")

    tuneall <- srf_tune(
      df_tune= df_train,
      cv_number = inner_cv,
      eligible_params(predict.factors, df_train),
      fixed_time = fixed_time,
      randomseed = randomseed,
      mtry = mtry,
      nodesize = nodesize,
      nodedepth = nodedepth,
      oob = oob
    )
    nodesize_best <- tuneall$bestauc["nodesize"]
    nodedepth_best <- tuneall$bestauc["nodedepth"]
    mtry_best <- tuneall$bestauc["mtry"]
    best_combo_stat <- tuneall$bestauc
    modelstatsall <- tuneall$modelstats
  }

  final.rfs <- randomForestSRC::rfsrc(
    as.formula(paste(
      "Surv(time, event) ~",
      paste(predict.factors, collapse = "+")
    )),
    data = df_train,
    nodesize = nodesize_best,
    # this is AVERAGE size, so we want this to be quite high
    ntree = 500,
    mtry = as.integer(mtry_best),
    nodedepth = as.integer(nodedepth_best),
    nsplit = 50,
    splitrule = "logrank",
    statistics = FALSE,
    membership = TRUE,
    importance = "none",
    # to speed up by switching off VIMP calculations
    seed = randomseed
  )
  output <- list()
  output$beststats <- best_combo_stat
  output$allstats <- modelstatsall
  output$model <- final.rfs
  # calibrate SRF with the best parameters
  return(output)
}

#' Predicts event probability for a fitted SRF model
#'
#' @description
#' Predicts event probability for a fitted SRF model randomForestSRC::rfsrc.
#' Essentially a wrapper of \link{srf_survival_prob_for_time}.
#' @param model_srf trained model
#' @param df_test test data
#' @param fixed_time  time for which probabilities are computed
#' @param oob TRUE/FALSE , default is FALSE
#' @examples
#' df <- simsurv_nonlinear()
#' @return returns matrix(nrow = length(newdata), ncol = length(times))
survrf_predict <- function(model_srf,
                               df_test,
                               fixed_time,
                               oob = FALSE) {
  if (class(model_srf) == "list") {
    srf <- model_srf$model
  } else {
    srf <- model_srf
  }

  if (length(fixed_time) == 1) {
    return(1 - srf_survival_prob_for_time(srf, df_test, fixed_time, oob = oob))
  }

  predicted_event_prob <-
    matrix(nrow = dim(df_test)[1], ncol = length(fixed_time))
  for (t in 1:length(fixed_time)) {
    predicted_event_prob[, t] <-
      1 - srf_survival_prob_for_time(srf, df_test, fixed_time[t], oob = oob)
  }
  colnames(predicted_event_prob) <- round(fixed_time, 3)
  return(predicted_event_prob)
}

#' Cross-validates SRF model
#'
#' @param df data frame with the data, "time" and "event" for survival outcome
#' @param predict.factors list of predictor names
#' @param fixed_time  at which performance metrics are computed
#' @param cv_number k in k-fold CV, default 3
#' @param inner_cv kk in the inner look of kk-fold CV, default 3
#' @param randomseed random seed
#' @param srf_tuning list of mtry, nodedepth and nodesize, default is NULL
#' @param parallel TRUE/FALSE, default = FALSE
#' @param return_models TRUE/FALSE, if TRUE returns all CV objects
#' @param repeat_cv if NULL, runs once, otherwise repeats CV
#' @examples
#' df <- simsurv_nonlinear()
#' @return output list: output$train, test, testaverage, traintaverage, time
survrf_cv <- function(df,
                          predict.factors,
                          fixed_time = NaN,
                          cv_number = 3,
                          inner_cv = 3,
                          randomseed = NULL,
                          srf_tuning = NULL,
                          parallel = FALSE,
                          return_models = FALSE,
                          repeat_cv = NULL) {
  time_0 <- Sys.time()

  if (is.null(randomseed)) {randomseed <- round(stats::runif(1)*1e9,0)}

  stopifnot(expr = {
    is.data.frame(df)
    predict.factors %in% colnames(df)
  })

  if (sum(is.nan(fixed_time)) > 0) {
    fixed_time <- round(quantile(df[df$event == 1, "time"], 0.9), 1)
  }

  predict.factors <- eligible_params(predict.factors, df)
  if (length(predict.factors) == 0) {
    print("No eliible params")
    return(NULL)
  }

  #defining number of repeated cv
  if (is.null(repeat_cv)) {repeat_cv = 1}
  if (is.numeric(repeat_cv) & repeat_cv > 1) {
    repeat_cv = round(repeat_cv, 0)
  }else{
    repeat_cv = 1
  }

  print(paste("Cross-validating SRF ( ",repeat_cv, " repeat(s), ",
              cv_number," outer, ",inner_cv," inner loops)",sep=""))

  modelstats_train <- list()
  modelstats_test <- list()
  models_for_each_cv <- list()
  #progress bar
  pb <- utils::txtProgressBar(0, cv_number*repeat_cv, style = 3)
  utils::setTxtProgressBar(pb, cv_number*repeat_cv / 50)

  for (rep_cv in 1:repeat_cv) {
    set.seed(randomseed + rep_cv)
    if (rep_cv != 1) {df <- df[sample(1:nrow(df)),]}
    cv_folds <-
      caret::createFolds(df$event, k = cv_number, list = FALSE)

    for (cv_iteration in 1:cv_number) {
      # print(paste("CV ",cv_iteration,"/",cv_number, sep =""))
      df_train_cv <- df[cv_folds != cv_iteration, ]
      df_test_cv <- df[cv_folds == cv_iteration, ]
      srf.model.tuned <-
        survrf_train(
          df_train_cv,
          predict.factors,
          fixed_time = fixed_time,
          inner_cv = inner_cv,
          randomseed = randomseed,
          srf_tuning = srf_tuning,
          fast_version = TRUE,
          oob = TRUE
        )
      y_predict_test <- survrf_predict(srf.model.tuned, df_test_cv,
                                       fixed_time, oob = FALSE)
      y_predict_train <- survrf_predict(srf.model.tuned, df_train_cv,
                                        fixed_time, oob = FALSE)
      modelstats_test[[cv_iteration + (rep_cv - 1) * cv_number]] <-
        survval(y_predict_test,fixed_time,df_train_cv,
                df_test_cv,weighted = 1)
      modelstats_train[[cv_iteration + (rep_cv - 1) * cv_number]] <-
        survval(y_predict_train,fixed_time,df_train_cv,
                df_train_cv,weighted = 1)
      if (return_models) {
        models_for_each_cv[[cv_iteration + (rep_cv - 1) * cv_number]] <-
          srf.model.tuned$model
      }
      utils::setTxtProgressBar(pb, cv_iteration + (rep_cv - 1) * cv_number)
    }
  }

  df_modelstats_test <- data.frame(modelstats_test[[1]])
  df_modelstats_train <- data.frame(modelstats_train[[1]])

  for (i in 2:(repeat_cv*cv_number)) {
    df_modelstats_test[i, ] <- modelstats_test[[i]]
    df_modelstats_train[i, ] <- modelstats_train[[i]]
  }
  row.names(df_modelstats_train)<- 1:(cv_number*repeat_cv)
  row.names(df_modelstats_test)<- 1:(cv_number*repeat_cv)

  utils::setTxtProgressBar(pb, cv_number*repeat_cv)
  close(pb)

  df_modelstats_test$test <- 1
  df_modelstats_train$test <- 0
  output <- list()
  output$test <- df_modelstats_test
  output$train <- df_modelstats_train
  output$testaverage <- sapply(df_modelstats_test, mean, na.rm = 1)
  output$trainaverage <-  sapply(df_modelstats_train, mean, na.rm = 1)
  output$pretrained_srf_models <- models_for_each_cv
  time_1 <- Sys.time()
  print(time_1 - time_0)
  output$time <- time_1 - time_0
  return(output)
}

