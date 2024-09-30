##################################################################

#' Prints trained survensemble object
#'
#'@param x survensemble object
#'@param ... additional arguments to be passed
#'@return x
#'@export
print.survensemble_cv <- function(x, ...) {
  if (!inherits(x, "survensemble")) {stop("Not a \"survensemble\" object")}
  summary.survensemble_cv(x)
}

#' Prints summary of a trained survensemble_cv object
#'
#'@param object survensemble_cv object
#'@param ... additional arguments to be passed
#'@return object
#'@export
summary.survensemble_cv <- function(object, ...) {
  if (!inherits(object, "survensemble")) {
    stop("Not a \"survensemble\" object")
  }
  cat("Survival ensemble", object$model_name ,"\n")
  if (!is.null(cl <- object$call)) {
    cat("Call:\n")
    dput(cl)
  }
  if(!is.null(object$lambda)){
    cat("\n=> Lambda (ML contribution share) :", object$lambda, "\n")
  }
  cat("\n=> Items available as object$item are: ")
  cat(names(object), sep = ", ")
}

##################################################################
#' Prints survensemble_cv object
#'
#'@param x survensemble_cv object
#'@param ... additional arguments to be passed
#'@return x
#'@export
print.survensemble_cv <- function(x, ...) {
  if (!inherits(x, "survensemble_cv")){stop("\nNot a \"survensemble_cv\" object")}
  summary.survensemble_cv(x)
}

#' Prints a summary of survensemble_cv object
#'
#'@param object survensemble_cv object
#'@param ... additional arguments to be passed
#'@return object
#'@export
summary.survensemble_cv <- function(object, ...) {
  if (!inherits(object, "survensemble_cv")) {
    stop("Not a \"survensemble_cv\" object")
  }
  cat("Cross-validation results. Computation time:",
      round(object$time,2), "sec. \n")
  if (!is.null(cl <- object$call)) {
    cat("Call:\n")
    dput(cl)
  }
  print(round(object$testaverage,4))
  cat("\nThe stats are computed from the ", dim(object$test)[1]," data splits.\n")
  print(
    object$summarydf[!rownames(object$summarydf) %in% c("repeat_cv", "outer_cv"), ])
}
