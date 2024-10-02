
#' Predicts event probability by a trained Survival Random Forest
#'
#' @param trained_model a trained SRF model, output of survsrf_train(), or randomForestSRC::rfsrc()
#' @param newdata new data for which predictions are made
#' @param fixed_time time of interest for which event probabilities are computed
#' @param extrapsurvival if probabilities are extrapolated beyond trained times (using probability of
#' the lastest available time). Can be helpful for cross-validation of small data, where random
#' split may cause the time of interest being outside of the training set.
#' @return vector of predicted event probabilities
#' @export
survsrf_predict <-
  function(trained_model,
           newdata,
           fixed_time,
           extrapsurvival = TRUE
  )
  {
    if (inherits(trained_model, "list")) {trained_model <- trained_model$model}
    predictions <- predict(trained_model, newdata = newdata)
    train_times <- predictions$time.interest
    s1 <- predictions$survival
    if(fixed_time > max(train_times)) {
      if(!extrapsurvival) {return(rep(NaN, dim(newdata)[1]))
      }else{
        fixed_time = max(train_times)
      }
    }
    f <- function(i) {
      approxfun(train_times, s1[i, ], method = "constant")(fixed_time)}
    predict_eventprob <- 1 - unlist(lapply(1:dim(s1)[1], f))
    return(predict_eventprob)
  }

#' Fits randomForestSRC, with tuning by mtry, nodedepth, and nodesize.
#' Underlying model is by Ishwaran et al(2008)
#' https://www.randomforestsrc.org/articles/survival.html
#' Ishwaran H, Kogalur UB, Blackstone EH, Lauer MS. Random survival forests.
#' The Annals of Applied Statistics. 2008;2:841–60.
#'
#' @param df_train  data, "time" and "event" should describe survival outcome
#' @param predict.factors list of predictor names
#' @param fixed_time time at which performance is maximized
#' @param tuningparams if given, list of hyperparameters, list(mtry=c(), nodedepth=c(),nodesize=c()), otherwise a wide default grid is used
#' @param max_grid_size number of random grid searches for model tuning
#' @param inner_cv number of cross-validation folds for hyperparameters' tuning
#' @param randomseed random seed to control tuning including data splits
#' @param verbose TRUE/FALSE, FALSE by default
#' @return output = list(bestparams, allstats, model)
#' @examples
#' d <-simulate_nonlinear(100)
#' p<- names(d)[1:4]
#' tuningparams = list(
#'  "mtry" = c(5,10,15),
#'  "nodedepth" = c(5,10,15,20),
#'  "nodesize" =    c(20,30,50)
#' )
#' m_srf<- survsrf_train(d,p,tuningparams=tuningparams)
#' @export
survsrf_train <-
  function(df_train,
           predict.factors,
           fixed_time = NaN,
           tuningparams = list(),
           max_grid_size = 10,
           inner_cv = 3,
           randomseed = NaN,
           verbose = TRUE) {

    if (is.nan(randomseed)) {
      randomseed <- round(stats::runif(1) * 1e9, 0)
    }
    grid_of_hyperparams <-
      ml_hyperparams_srf(
        mlparams = tuningparams,
        dftune_size = dim(df_train)[1],
        max_grid_size = max_grid_size,
        randomseed = randomseed,
        p = length(predict.factors)
      )
    if (dim(grid_of_hyperparams)[1] == 1) {
      #if only one option is given, we use it to train
      bestparams = grid_of_hyperparams[1,]
    } else{
      #if there are options, we tune the grid of hyperparams
      tuning_m <- survsrf_tune(
        df_train,
        predict.factors,
        repeat_tune = 1,
        fixed_time = fixed_time,
        tuningparams = tuningparams,
        max_grid_size = max_grid_size,
        inner_cv = inner_cv,
        randomseed = randomseed
      )
      bestparams = tuning_m$bestparams
    }
    if (verbose) {print(bestparams)}
    srfm <- randomForestSRC::rfsrc(
      as.formula(paste("Surv(time, event) ~",paste(predict.factors, collapse = "+"))),
      data = df_train,
      nodesize = bestparams$nodesize, #AVERAGE size, so needs to be large
      ntree = 300,
      mtry = bestparams$mtry,
      nodedepth = bestparams$nodedepth,
      nsplit = 50,
      splitrule = "logrank",
      statistics = FALSE,
      membership = TRUE,
      importance = "none", # to speed up
      seed = randomseed
    )
    output = list()
    output$model_name = "Survival Random Forest"
    output$model = srfm
    output$bestparams = bestparams
    output$grid = grid_of_hyperparams
    return(output)
  }

#'A repeated 3-fold CV over a hyperparameters grid
#'
#' @param repeat_tune number of repeats
#' @param df_tune data
#' @param predict.factors list of predictor names
#' @param fixed_time  not used here, but for some models the time for which performance is optimized
#' @param tuningparams if given, list of hyperparameters, list(mtry=c(), nodedepth=c(),nodesize=c()), otherwise a wide default grid is used
#' @param inner_cv number of cross-validation folds for hyperparameter tuning
#' @param max_grid_size number of random grid searches for model tuning
#' @param randomseed to choose random subgroup of hyperparams
#' @return  output=list(cindex_ordered, bestparams)
survsrf_tune <-
  function(df_tune,
           predict.factors,
           repeat_tune=1,
           fixed_time = NaN,
           tuningparams = list(),
           max_grid_size = 10,
           inner_cv = 3,
           randomseed = NaN) {
    srfgrid <-
      ml_hyperparams_srf(
        mlparams = tuningparams,
        dftune_size = dim(df_tune)[1],
        max_grid_size = max_grid_size,
        randomseed = randomseed,
        p = length(predict.factors)
      )
    means = c()
    for (i in (1:repeat_tune)) {
      srf_tune_fus <-
        survsrf_tune_single(
          df_tune = df_tune,
          predict.factors =  predict.factors,
          fixed_time = fixed_time,
          grid_hyperparams = srfgrid,
          inner_cv = inner_cv,
          randomseed=randomseed
        )
      means <- cbind(means, srf_tune_fus$cindex_mean)
      remove(srf_tune_fus)
    }
    allmeans <- apply(means, 1, mean)
    # grid and average cindex ordered by cindex (highest first)
    output = list()
    output$cindex_ordered <-
      cbind(srfgrid, allmeans)[order(allmeans, decreasing = TRUE),]
    output$bestparams <- output$cindex_ordered[1,]
    return(output)
  }


#' Internal function for survsrf_tune(), performs 1 CV
#' @param df_tune data
#' @param predict.factors list of predictor names
#' @param fixed_time predictions for which time are computed for c-index
#' @param grid_hyperparams  hyperparameters grid (or a default will be used )
#' @param inner_cv number of folds for each CV
#' @param randomseed randomseed
#' @param progressbar FALSE(default)/TRUE
#' @return  output=list(grid, cindex, cindex_mean)
survsrf_tune_single <-
  function(df_tune,
           predict.factors,
           fixed_time = NaN,
           grid_hyperparams=c(),
           inner_cv = 3,
           randomseed = NaN,
           progressbar = FALSE) {
    #fixed_time
    if (is.nan(fixed_time) | length(fixed_time) > 1) {
      # not implemented for multiple time
      fixed_time <-
        quantile(df_tune[df_tune$event == 1, "time"], 0.9)
    }
    grid_size <- dim(grid_hyperparams)[1]

    if (length(grid_hyperparams) == 0) {
      grid_hyperparams <-
        ml_hyperparams_srf(mlparams = list(),
                           dftune_size = dim(df_tune)[1],
                           p = length(predict.factors),
                           randomseed = randomseed)
    }
    #placeholder for c-index
    cind = matrix(NA, nrow = grid_size, ncol = inner_cv)
      #progress bar
    if(progressbar) {
      pb <- utils::txtProgressBar(0, inner_cv * grid_size, style = 3)
      utils::setTxtProgressBar(pb, inner_cv * grid_size / 50)
    }
    # tuning cross-validation loop
    for (cv_iteration in 1:inner_cv) {
      if (!is.nan(randomseed)) {
        set.seed(randomseed + 123 + cv_iteration)
      }
      cv_folds <-
        caret::createFolds(df_tune$event, k = inner_cv, list = FALSE)
      df_train_cv <- df_tune[cv_folds != cv_iteration,]
      df_test_cv <- df_tune[cv_folds == cv_iteration,]
      #Grid search
      for (i in 1:grid_size) {
        srfm <-
          randomForestSRC::rfsrc(
            as.formula(paste(
              "Surv(time, event) ~",
              paste(predict.factors, collapse = "+")
            )),
            data = df_train_cv,
            nodesize =  grid_hyperparams[i, "nodesize"],
            # this is AVERAGE size, so we want this to be quite high
            ntree = 300,
            mtry = grid_hyperparams[i, "mtry"],
            nodedepth = grid_hyperparams[i, "nodedepth"],
            nsplit = 50,
            splitrule = "logrank",
            statistics = FALSE,
            membership = TRUE,
            importance = "none"
          )
        #check test performance
        pp <-
          survsrf_predict(
            trained_model = srfm,
            newdata = df_test_cv,
            fixed_time = fixed_time
          )
        cind[i, cv_iteration] =
          surv_validate(pp, fixed_time, df_train_cv, df_test_cv)[1, "C_score"]

        if (progressbar) {
          utils::setTxtProgressBar(pb, grid_size + (i - 1) * cv_iteration)
        }
      } # here we already have cindex for allgrid for cv_iteration
    }
    if(progressbar) {
      utils::setTxtProgressBar(pb, inner_cv * grid_size)
      close(pb)
    }
    remove(srfm)
    cindex_mean <- apply(cind, 1, mean)
    output = list()
    output$cindex = cind
    output$cindex_mean = cindex_mean
    output$cindex_ordered <-
      cbind(grid_hyperparams, cindex_mean)[order(cindex_mean, decreasing = TRUE), ]
    output$bestparams <- output$cindex_ordered[1, ]
    return(output)
  }


#' Internal function for getting grid of hyperparameters
#' for random or grid search of size = max_grid_size
#' @param mlparams list of params
#' @param dftune_size size of the tuning data to define nodesize options
#' @param p number of predictors to detine mtry options
#' @param max_grid_size grid size for tuning
#' @param randomseed randomseed to select the tuning grid
#' @export
ml_hyperparams_srf <- function(mlparams = list(),
                              p = 10,
                              max_grid_size = 10,
                              dftune_size = 1000,
                              randomseed = NaN) {
  if (p<=10) {mtry = c(2,3,4,5)
  }else if (p<25) {mtry = c(3,5,7,10,15)
  }else{mtry = round(c(p/10, p/5, p/3, p/2, max(1, sqrt(p))),0)  }

  nodesize <- seq(15, max(min(dftune_size/10, 50), 30), 5)
  nodedepth <- c(2,3,5,7,10,15,50)
  default_grid <- list(mtry = mtry, nodesize = nodesize, nodedepth = nodedepth)

  if (length(mlparams) == 0) {
    grid_of_hyperparams <- expand.grid(default_grid)
  } else{
    if (is.null(mlparams$mtry)) {mlparams$mtry = default_grid$mtry}
    if (is.null(mlparams$nodesize)) {mlparams$nodesize = default_grid$nodesize}
    if (is.null(mlparams$nodedepth)) {mlparams$nodedepth = default_grid$nodedepth}

    mlparams$mtry <- mlparams$mtry[mlparams$mtry <= p]
    if (length(mlparams$mtry)==0 ) {mlparams$mtry <- c(max(1,round(sqrt(p),0)))}

    grid_of_hyperparams <- expand.grid(
      "mtry" = mlparams$mtry,
      "nodesize" = mlparams$nodesize,
      "nodedepth" = mlparams$nodedepth
      )
  }
  grid_size <- dim(grid_of_hyperparams)[1]
  final_grid <- grid_of_hyperparams
  if (grid_size > max_grid_size) {
    if (!is.nan(randomseed)) {set.seed(randomseed)}
    final_grid <-
      grid_of_hyperparams[sample(1:grid_size, max_grid_size, replace = FALSE),]
    grid_size <- max_grid_size
  }
  remove(grid_of_hyperparams); remove(default_grid)
  return(final_grid)
}

#' Cross-validates Survival Random Forest
#'
#' @param df  data, "time" and "event" should describe survival outcome
#' @param predict.factors list of predictor names
#' @param fixed_time time at which performance is maximized
#' @param outer_cv number of cross-validation folds for model validation
#' @param inner_cv number of cross-validation folds for hyperparameters' tuning
#' @param repeat_cv number of CV repeats, if NaN, runs once
#' @param randomseed random seed to control tuning including data splits
#' @param return_models if all models are stored and returned
#' @param tuningparams if given, list of hyperparameters, list(mtry=c(), nodedepth=c(),nodesize=c()), otherwise a wide default grid is used
#' @param max_grid_size number of random grid searches for model tuning
#' @param verbose FALSE(default)/TRUE
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
                       inner_cv = 3,
                       repeat_cv = 2,
                       randomseed = NaN,
                       return_models = FALSE,
                       tuningparams = list(),
                       max_grid_size = 10,
                       verbose = FALSE,
                       parallel = FALSE) {

  Call <- match.call()
  inputs <- list(df,
                 predict.factors,
                 fixed_time,
                 outer_cv,
                 inner_cv,
                 repeat_cv,
                 randomseed,
                 return_models,
                 tuningparams,
                 max_grid_size,
                 verbose
                 )
  inputclass<- list(df = "data.frame",
                    predict.factors = "character",
                    fixed_time = "numeric",
                    outer_cv = "numeric",
                    inner_cv = "numeric",
                    repeat_cv = "numeric",
                    randomseed = "numeric",
                    return_models = "logical",
                    tuningparams = "list",
                    max_grid_size = "numeric",
                    verbose = "logical")

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
    model_args = list(tuningparams = tuningparams,
                      max_grid_size= max_grid_size,
                      fixed_time = fixed_time,
                      randomseed = randomseed,
                      verbose = verbose),
    model_name = "Survival Random Forest",
    parallel = parallel
  )
  output$call <- Call
  return(output)
}
