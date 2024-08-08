
################  stack_deephit_train ####################
# Create out-of-bag Cox predictions, then train deephit
#' @export
stack_deephit_train <-
  function(df_train,
           predict.factors,
           fixed_time = NaN,
           inner_cv = 3,
           randomseed = NaN,
           useCoxLasso = FALSE,
           tuningparams = list(),
           max_grid_size =10,
           verbose = FALSE) {

    # setting fixed_time if not given
    if (is.nan(fixed_time)| (length(fixed_time) > 1)) {
      fixed_time <-
        round(quantile(df_train[df_train$event == 1, "time"], 0.9), 1)
    }
    if(verbose) cat("\nTraining baseline learners ...\n")
    #base DeepHit
    ml_base_model <-
      deephit_train(df_train = df_train,
                    predict.factors = predict.factors,
                    fixed_time = fixed_time,
                    tuningparams = tuningparams,
                    max_grid_size = max_grid_size,
                    inner_cv = inner_cv,
                    randomseed = randomseed )

    #base cox model
    cox_base_model <-
      survcox_train(df_train, predict.factors, useCoxLasso = useCoxLasso)

    #------ Computationally demanding fit of the stack model:---
    #first, using k-fold CV, underlying models are trained (using already tuned parameters though)
    #and out-of-sample predictions for Cox and ML are computed
    #then, a regression (outcome ~ cox_predictions + ml_predictions)
    #constitutes the meta-learner
    #meta-learner: compute out-of-sample Cox and ML (DeepHit) predictions
    if(verbose) cat("\nTraining meta-learner... computing out-of-sample predictions...\t")

    if (!is.nan(randomseed)) {set.seed(randomseed)}
    tuningparams_tuned = list(learning_rate = ml_base_model$bestparams$learning_rate,
                              dropout = ml_base_model$bestparams$dropout,
                              num_nodes = ml_base_model$bestparams$num_nodes,
                              batch_size = ml_base_model$bestparams$batch_size,
                              epochs = ml_base_model$bestparams$epochs,
                              mod_alpha = ml_base_model$bestparams$mod_alpha,
                              sigma = ml_base_model$bestparams$sigma,
                              batch_norm = ml_base_model$bestparams$batch_norm,
                              cuts = ml_base_model$bestparams$cuts,
                              weight_decay = ml_base_model$bestparams$weight_decay,
                              early_stopping = ml_base_model$bestparams$early_stopping
    )
    bestparams_base <- ml_base_model$bestparams

    #avoid variance in predictions for small data, use higher number of folds(10)
    k_for_oob = ifelse(dim(df_train)[1]<=350, 10, 5)

    cv_folds <-
      caret::createFolds(df_train$event, k = k_for_oob, list = FALSE)
    cindex_train <- vector(length = k_for_oob)
    cindex_test <- vector(length = k_for_oob)

    for (cv_iteration in 1:k_for_oob) {
      if(verbose) cat("\t", cv_iteration, "/", k_for_oob)
      data_train <- df_train[cv_folds != cv_iteration, ]
      data_oob <- df_train[cv_folds == cv_iteration, ]
      # train cox model on data_train
      cox_m_cv <-
        survcox_train(data_train,
                      predict.factors,
                      useCoxLasso = useCoxLasso)
      # predict for data_oob
      cox_predict_oob <-
        survcox_predict(cox_m_cv, data_oob, fixed_time)[,1]
      # adding Cox prediction to the df_train in the column "cox_predict"
      df_train[cv_folds == cv_iteration, "cox_predict"] <- cox_predict_oob

      #train ML model on data_train
      deephit.cv <-
        deephit_train(df_train = data_train,
                      predict.factors = predict.factors,
                      fixed_time = fixed_time,
                      tuningparams = tuningparams_tuned,
                      max_grid_size = 1,
                      inner_cv = inner_cv,
                      randomseed = randomseed + cv_iteration,
                      verbose = FALSE)
      # predict for data_oob
      ml_predict_oob <-
        deephit_predict(trained_model = deephit.cv,
                        newdata = data_oob,
                        predict.factors = predict.factors,
                        fixed_time= fixed_time)
      # adding ML prediction to the df_train in the column "ml_predict"
      df_train[cv_folds == cv_iteration, "ml_predict"] <- ml_predict_oob

    }

    if(verbose) cat("\t computing meta-learner's weights ...")

    # 1/0 by fixed_time:
    df_train$event_t <-
      ifelse(df_train$time <= fixed_time & df_train$event == 1, 1, 0)

    # Excluding censored observations before fixed_time, leave those with known state
    df_train_in_scope <-
      df_train[(df_train$time >= fixed_time) |
                 (df_train$time < fixed_time & df_train$event == 1),]

    # p = 0+ p1 + lambda*(p2-p1), const==0  OR p=(1-lambda)p1 + lambda*p2
    stack_model <-
      lm(event_t ~ 0 + offset(cox_predict) + I(ml_predict - cox_predict),
         data = df_train_in_scope)

    # same, different fitting for lambda


    #alt1) p = b0 + b1*p1 + b2*p2
    stack_model_2 <-
      lm(event_t ~ cox_predict + ml_predict,
         data = df_train_in_scope)


    #alt2) ln(p/1-p) = b1* ln(p1/1-p1)+b2*ln(p2/1-p2)+c
    # logit = function(x) {
    #   log(pmax(pmin(x, 0.9999), 0.0001) / (1 - pmax(pmin(x, 0.9999), 0.0001)))
    # }
    # df_train$cox_predict_logit <- logit(df_train$cox_predict)
    # df_train$ml_predict_logit <- logit(df_train$ml_predict)
    # stack_model_2 <-
    #   glm(event_t ~ cox_predict_logit + ml_predict_logit,
    #       data = df_train_in_scope,
    #       family = "binomial")

    bestparams_meta <-
      c(stack_model$coefficients[1],
        stack_model$coefficients[2],
        stack_model$coefficients[3],
        stack_model_2$coefficients[1],
        stack_model_2$coefficients[2],
        stack_model_2$coefficients[3]
      )

    #output
    output = list()
    output$model_name <- "Stacked_DeepHit_CoxPH"
    output$oob_predictions <- df_train
    output$model <- stack_model
    output$model_2 <- stack_model_2
    output$model_base_cox <- cox_base_model
    output$model_base_ml <- ml_base_model
    output$randomseed <- randomseed
    output$bestparams <- c(bestparams_base, bestparams_meta)
    output$call <-  match.call()
    class(output) <- "survensemble"
    return(output)
  }

################  stack_deephit_predict ####################

#same as deephit_predict
#' @export
stack_deephit_predict <-
  function(trained_object,
           newdata,
           fixed_time,
           predict.factors,
           use_alternative_model = FALSE
  ) {
    trained_model <- trained_object$model
    if(use_alternative_model) {trained_model <- trained_object$model_2}
    predictdata <- newdata[predict.factors]

    # use model_base with the base Cox model to find cox_predict
    predictdata$cox_predict <-
      survcox_predict(trained_model = trained_object$model_base_cox,
                      newdata = newdata, fixed_time = fixed_time)[,1]
    predictdata$ml_predict <-
      deephit_predict(trained_model = trained_object$model_base_ml,
                      newdata = newdata, predict.factors = predict.factors,
                      fixed_time = fixed_time)
    logit = function(x) {
      log(pmax(pmin(x, 0.9999), 0.0001) / (1 - pmax(pmin(x, 0.9999), 0.0001)))
    }
    predictdata$cox_predict_logit <- logit(predictdata$cox_predict)
    predictdata$ml_predict_logit <- logit(predictdata$ml_predict)

    p <- predict(trained_model, newdata = predictdata)
    return(p)
  }


############### stack_deephit_cv #############
#' @export
stack_deephit_cv <- function(df,
                           predict.factors,
                           fixed_time = NaN,
                           outer_cv = 3,
                           inner_cv = 3,
                           repeat_cv = 2,
                           randomseed = NaN,
                           return_models = FALSE,
                           useCoxLasso = FALSE,
                           tuningparams = list(),
                           max_grid_size =10
) {
  Call <- match.call()

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
    train_function = stack_deephit_train,
    predict_function = stack_deephit_predict,
    model_args = list("tuningparams" = tuningparams,
                      "useCoxLasso" = useCoxLasso,
                      "max_grid_size" = max_grid_size,
                      "randomseed" = randomseed),
    predict_args = list("predict.factors" = predict.factors),
    model_name = "Stacked_DeepHit_CoxPH"
  )
  output$call <- Call
  return(output)
}
