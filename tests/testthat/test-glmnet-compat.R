test_that("Cox-Lasso fits using explicit folds and supported glmnet controls", {
  d <- cox_test_data()
  original_cv <- glmnet::cv.glmnet
  observed <- NULL
  fit_completed <- FALSE
  local_mocked_bindings(
    cv.glmnet = function(x, y, ...) {
      observed <<- list(...)
      result <- original_cv(x = x, y = y, ...)
      fit_completed <<- TRUE
      result
    },
    .package = "glmnet"
  )
  for (requested in c(3, 5)) {
    set.seed(31)
    fit_completed <- FALSE
    messages <- character()
    fit <- withCallingHandlers(
      survcox_train(d, c("x1", "x2"), useCoxLasso = TRUE,
                    inner_cv = requested),
      warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_true(fit_completed)
    expect_s3_class(fit, "coxph")
    expect_true(all(is.finite(stats::predict(fit, d, type = "risk"))))
    expect_identical(observed$nfolds, requested)
    expect_false("nfold" %in% names(observed))
    glmnet_formals <- names(formals(glmnet::glmnet))
    if ("control" %in% glmnet_formals) {
      expect_equal(observed$control$maxit, 5000)
      expect_false("maxit" %in% names(observed))
    } else {
      expect_equal(observed$maxit, 5000)
      expect_false("control" %in% names(observed))
    }
    if ("cox.ties" %in% glmnet_formals) {
      expect_identical(observed$cox.ties, "breslow")
    } else {
      expect_false("cox.ties" %in% names(observed))
    }
    expect_false(any(grepl("deprecated|partial argument match|default.*ties|ties.*default",
                          messages, ignore.case = TRUE)), info = paste(messages, collapse = "\n"))
  }
})
