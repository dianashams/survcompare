test_that("survcompare forwards tuning folds to its Cox comparison", {
  for (requested in c(3, 5)) {
    local({
      observed <- NULL
      local_mocked_bindings(
        survcox_cv = function(df, predict.factors, inner_cv = -1,
                              useCoxLasso = FALSE, ...) {
          observed <<- inner_cv
          expect_true(useCoxLasso)
          stop("Cox wrapper reached")
        },
        .package = "survcompare"
      )
      expect_error(survcompare(cox_test_data(), c("x1", "x2"),
                               fixed_time = 1, inner_cv = requested,
                               useCoxLasso = TRUE, suppresswarn = FALSE),
                   "Cox wrapper reached")
      expect_identical(observed, requested)
    })
  }
})

test_that("each CV wrapper delivers requested folds to the trainer", {
  d <- cox_test_data()
  pairs <- list(c("survcox_cv", "survcox_train"),
                c("survsrf_cv", "survsrf_train"),
                c("survsrfens_cv", "survsrfens_train"),
                c("survsrfstack_cv", "survsrfstack_train"))
  for (pair in pairs) {
    for (requested in c(3, 5)) {
      local({
        observed <- NULL
        trainer <- function(df_train, predict.factors, inner_cv = -1, ...) {
          observed <<- inner_cv
          stop("trainer reached")
        }
        do.call(local_mocked_bindings,
                c(stats::setNames(list(trainer), pair[2]),
                  list(.package = "survcompare", .env = environment())))
        args <- list(df = d, predict.factors = c("x1", "x2"),
                     fixed_time = 1, inner_cv = requested,
                     outer_cv = 2, repeat_cv = 1, randomseed = 5,
                     suppresswarn = FALSE)
        if (pair[1] != "survsrf_cv") args$useCoxLasso <- TRUE
        expect_error(do.call(get(pair[1], asNamespace("survcompare")), args),
                     "trainer reached")
        expect_identical(observed, requested, info = paste(pair, collapse = " -> "))
      })
    }
  }
})

test_that("ensemble Cox fits and forest tuning use the same requested folds", {
  d <- cox_test_data()
  for (method in c("survsrfens_train", "survsrfstack_train")) {
    for (requested in c(3, 5)) {
      local({
        cox_folds <- forest_folds <- numeric()
        local_mocked_bindings(
          survcox_train = function(df_train, predict.factors,
                                   useCoxLasso, inner_cv = -1, ...) {
            expect_true(useCoxLasso)
            cox_folds <<- c(cox_folds, inner_cv)
            survival::coxph(survival::Surv(time, event) ~ x1 + x2,
                            data = df_train, x = TRUE)
          },
          survsrf_train = function(df_train, predict.factors,
                                   inner_cv = -1, ...) {
            forest_folds <<- c(forest_folds, inner_cv)
            list(model = NULL, bestparams = list(mtry = 1))
          },
          survsrf_predict = function(trained_model, newdata, fixed_time, ...) {
            stats::plogis(newdata$x1)
          },
          .package = "survcompare"
        )
        result <- do.call(get(method, asNamespace("survcompare")),
                          list(df_train = d, predict.factors = c("x1", "x2"),
                               fixed_time = 1, inner_cv = requested,
                               randomseed = 21, useCoxLasso = TRUE))
        expect_s3_class(result, "survensemble")
        expected_cox <- if (method == "survsrfens_train") 11L else 6L
        expected_forest <- if (method == "survsrfens_train") 1L else 6L
        expect_equal(cox_folds, rep(requested, expected_cox))
        expect_equal(forest_folds, rep(requested, expected_forest))
      })
    }
  }
})
