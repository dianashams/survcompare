############# SRF ################

#' Internal function to compute survival probability
#'  by time from a fitted survival random forest
#'
#' @param rfmodel pre-trained survsrf_train model
#' @param df_to_predict test data
#' @param fixed_time  at which event probabilities are computed
#' @param oob TRUE/FALSE use out-of-bag prediction
#' @examples
#' df <- simulate_nonlinear()
#' #params<- c("age", "hyp", "bmi")
#' #s <- survsrf_train(df, params)
#' #p <- survsrf_predict(s, df, 5)
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
#' @param inner_cv k in k-fold CV, applied if oob=FALSE
#' @param fixed_time  NaN
#' @param randomseed  random seed
#' @param mtry  tuning parameter
#' @param nodesize  tuning parameter
#' @param nodedepth  tuning parameter
#' @param verbose  FALSE
#' @param oob  TRUE
#' @param nodesize  at which event probabilities are computed
#' @param oob TRUE/FALSE use out-of-bag predictions while tuning instead of cross-validation, default is TRUE and is faster
#' @return  output=list(modelstats, bestbrier, bestauc, bestcindex)
survsrf_tune <- function(df_tune,
                         predict.factors,
                         inner_cv = 3,
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
  if (is.null(randomseed)) {
    randomseed <- round(stats::runif(1) * 1e9, 0)
  }
  set.seed(randomseed)

  # limit nodesize by 1/5 of sample size
  if (sum(nodesize > dim(df_tune)[1] / 5) > 0) {
    if (verbose) {
      print("Warning - some nodesize parameters are > 1/5 of the sample size and may prevent effective tree building.")
    }
  }

  # limit mtry by the number of predictors
  # if all numbers higher than number of factors, only check this factor
  if (sum(mtry > length(predict.factors)) == length(mtry)) {
    mtry <- c(length(predict.factors))
    if (verbose) {
      print("Warning - some mtry are > number of params and will be ignored.")
    }
  }
  mtry <- mtry[mtry <= length(predict.factors)]

  # grid of values to tune
  grid_of_values <- expand.grid("mtry" = mtry,
                                "nodesize" = nodesize,
                                "nodedepth" = nodedepth)
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
      caret::createFolds(df_tune$event, k = inner_cv, list = FALSE)
    for (i in 1:dim(grid_of_values)[1]) {
      if (verbose) {
        print(grid_of_values[i,])
      }
      mtry_i <- grid_of_values[i, "mtry"]
      nodesize_i <- grid_of_values[i, "nodesize"]
      nodedepth_i <- grid_of_values[i, "nodedepth"]

      # train grid combination for each cv_iteration
      modelstats_cv <- list()
      for (cv_iteration in 1:inner_cv) {
        print(paste("SRF tuning CV step", cv_iteration, "/out of", inner_cv))
        df_train_cv <- df_tune[cv_folds != cv_iteration,]
        df_test_cv <- df_tune[cv_folds == cv_iteration,]
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
          surv_validate(y_predicted,
                        fixed_time,
                        df_train_cv,
                        df_test_cv,
                        weighted = TRUE)
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
      for (j in 2:inner_cv) {
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
        print(grid_of_values[i,])
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
        surv_validate(y_predicted, fixed_time, df_tune, df_tune, weighted = TRUE)
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
    print(paste(
      "Concordance index (C-score) varies from",
      round(min(df_modelstats$C_score), 4),
      "to",
      round(max(df_modelstats$C_score), 4)
    ))
    print(paste(
      "Brier score varies from",
      round(min(df_modelstats$BS), 4),
      "to",
      round(max(df_modelstats$BS), 4)
    ))
    print("Combination with highest C_score")
    print(df_modelstats[which.max(df_modelstats$C_score),
                        c("mtry", "nodesize", "nodedepth")])
    print("Combination with lowest Brier Score")
    print(df_modelstats[which.min(df_modelstats$BS),
                        c("mtry", "nodesize", "nodedepth")])
  }
  output <- list()
  output$modelstats <- df_modelstats
  output$bestbrier <- df_modelstats[which.min(df_modelstats$BS),]
  output$bestauc <- df_modelstats[which.max(df_modelstats$AUCROC),]
  output$bestcindex <-
    df_modelstats[which.max(df_modelstats$C_score),]
  return(output)
}


#' Fits randomForestSRC, with tuning by mtry, nodedepth, and nodesize.
#' Underlying model is by Ishwaran et al(2008)
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
#' @param oob TRUE/FALSE use out-of-bag predictions while tuning SRF instead of cross-validation, default is TRUE and is faster
#' @param verbose TRUE/FALSE, FALSE by default
#' @return output = list(beststats, allstats, model)
#' @export
survsrf_train <- function(df_train,
                          predict.factors,
                          fixed_time = NaN,
                          inner_cv = 3,
                          randomseed = NULL,
                          srf_tuning = list(),
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
  mtry <-
    ifelse(p <= 10, c(2, 3, 4, 5),
           ifelse(p <= 25, c(3, 5, 7, 10, 15),
                  round(c(
                    p / 10, p / 5, p / 3, p / 2, mtry_default
                  ), 0)))
  mtry <- mtry[mtry <= length(predict.factors)]
  nodesize <-
    seq(min(15, round(n / 6 - 1, 0)), max(min(n / 10, 50), 30), 5)
  nodedepth <- c(5, 25)

  # update them if given in srf_tuning
  if (length(srf_tuning)>0) {
    if (!is.null(srf_tuning$nodesize)) {
      nodesize = srf_tuning$nodesize
    }
    if (!is.null(srf_tuning$mtry)) {
      mtry = c(srf_tuning$mtry)
    }
    if (!is.null(srf_tuning$nodedepth)) {
      nodedepth = srf_tuning$nodedepth
    }
  }

  if (verbose) {
    print("Tuning the following mtry, nodesize, and nodedepth:")
    print(mtry)
    print(nodesize)
    print(nodedepth)
  }

  if (fast_version == TRUE) {
    # fast version:
    # first, takes recommended mtry and optimizes by depth and node size,
    # second, optimizes by mtry
    tune1 <- survsrf_tune(
      df_tune = df_train,
      inner_cv = inner_cv,
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
    tune2 <- survsrf_tune(
      df_tune = df_train,
      inner_cv = inner_cv,
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
    tuneall <- survsrf_tune(
      df_tune = df_train,
      inner_cv = inner_cv,
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
  output$tuning <-
    list("mtry" = mtry,
         "nodesize" = nodesize,
         "nodedepth" = nodedepth)
  # calibrate SRF with the best parameters
  return(output)
}

#' Predicts event probability for a fitted SRF model
#'
#' @description
#' Predicts event probability for a fitted SRF model randomForestSRC::rfsrc.
#' Essentially a wrapper of \link{srf_survival_prob_for_time}.
#' @param trained_model trained model
#' @param newdata test data
#' @param fixed_time  time for which probabilities are computed
#' @param oob TRUE/FALSE use out-of-bag predictions while tuning instead of cross-validation, default is TRUE and is faster
#' @return returns vector of predictions (or matrix if fixed_time is a vector of times)
#' @export
survsrf_predict <- function(trained_model,
                            newdata,
                            fixed_time,
                            oob = FALSE) {
  if (inherits(trained_model, "list")) {
    srf <- trained_model$model
  } else {
    srf <- trained_model
  }

  if (length(fixed_time) == 1) {
    return(1 - srf_survival_prob_for_time(srf, newdata, fixed_time, oob = oob))
  }

  predicted_event_prob <-
    matrix(nrow = dim(newdata)[1], ncol = length(fixed_time))
  for (t in 1:length(fixed_time)) {
    predicted_event_prob[, t] <-
      1 - srf_survival_prob_for_time(srf, newdata, fixed_time[t], oob = oob)
  }
  colnames(predicted_event_prob) <- round(fixed_time, 3)
  return(predicted_event_prob)
}

#' Cross-validates SRF model
#'
#' @param df data frame with the data, "time" and "event" for survival outcome
#' @param predict.factors list of predictor names
#' @param fixed_time  at which performance metrics are computed
#' @param outer_cv k in k-fold CV, default 3
#' @param repeat_cv if NULL, runs once, otherwise repeats CV
#' @param randomseed random seed
#' @param return_models TRUE/FALSE, if TRUE returns all CV objects
#' @param inner_cv k in the inner loop of k-fold CV for SRF hyperparameters tuning, default is 3
#' @param srf_tuning list of tuning parameters for random forest: 1) NULL for using a default tuning grid, or 2) a list("mtry"=c(...), "nodedepth" = c(...), "nodesize" = c(...))
#' @param oob TRUE/FALSE use out-of-bag prediction accuracy while tuning instead of cross-validation, TRUE by default
#' @examples \donttest{
#' \dontshow{rfcores_old<- options()$rf.cores; options(rf.cores = 1)}
#' df <- simulate_nonlinear()
#' srf_cv <- survsrf_cv(df, names(df)[1:4])
#' summary(srf_cv)
#' \dontshow{options(rf.cores=rfcores_old)}
#' }
#' @return list of outputs
#' @export
survsrf_cv <- function(df,
                       predict.factors,
                       fixed_time = NaN,
                       outer_cv = 3,
                       repeat_cv = 2,
                       randomseed = NULL,
                       return_models = FALSE,
                       inner_cv = 3,
                       srf_tuning = list(),
                       oob = TRUE) {
  Call <- match.call()
  inputs <- list(df,
                 predict.factors,
                 fixed_time,
                 outer_cv,
                 repeat_cv,
                 randomseed,
                 return_models,
                 inner_cv,
                 srf_tuning,
                 oob)
  inputclass<- list(df = "data.frame",
                    predict.factors = "character",
                    fixed_time = "numeric",
                    outer_cv = "numeric",
                    repeat_cv = "numeric",
                    randomseed = "numeric",
                    return_models = "logical",
                    inner_cv = "numeric",
                    srf_tuning = "list",
                    oob = "logical")
  cp<- check_call(inputs, inputclass, Call)
  if (cp$anyerror) stop (paste(cp$msg[cp$msg!=""], sep=""))

  if (sum(is.na(df[c("time", "event", predict.factors)])) > 0) {
    stop("Missing data can not be handled. Please impute first.")
  }

  output <- surv_CV(
    df = df,
    predict.factors = predict.factors,
    fixed_time = fixed_time,
    outer_cv = outer_cv,
    inner_cv = inner_cv,
    repeat_cv = repeat_cv,
    randomseed = randomseed,
    return_models = return_models,
    train_function = survsrf_train,
    predict_function = survsrf_predict,
    model_args = list("srf_tuning" = srf_tuning, "oob" = oob),
    model_name = "Survival Random Forest"
  )
  output$call <- Call
  return(output)
}
