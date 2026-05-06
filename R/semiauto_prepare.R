#' Prepare Data for semiparametric regression
#'
#' This function prepares the cleaned response vector, predictor matrix,
#' selected model data, and time ordering before model fitting.
#'
#' @param data A data frame.
#' @param y A character string specifying the response variable name.
#' @param x A character vector specifying predictor variable names. If NULL,
#' all columns except y, time, and id will be used as predictors.
#' @param time Optional character string specifying the time variable.
#' @param id Optional character string specifying the subject or group ID variable.
#' @param na_action Missing value handling method. Options are "omit" or "fail".
#' @param verbose Logical. If TRUE, prints the preparation summary.
#'
#' @return An object of class \code{semiauto_prepare}.
#'
#' @examples
#' prepared <- semiauto_prepare(
#'   data = mtcars,
#'   y = "mpg",
#'   x = c("wt", "hp", "qsec")
#' )
#'
#' @export
#'
semiauto_prepare <- function(data,
                             y,
                             x = NULL,
                             time = NULL,
                             id = NULL,
                             na_action = c("omit", "fail"),
                             verbose = TRUE) {

  na_action <- match.arg(na_action)


# 1. Run basic input checking using semiauto_check()
  check_result <- semiauto_check(
    data = data,
    y = y,
    x = x,
    time = time,
    id = id,
    allow_missing = (na_action == "omit"),
    verbose = FALSE
  )

  if (!isTRUE(check_result$ok)) {
    stop(
      "Data preparation failed because input checking found problems:\n",
      paste(check_result$problems, collapse = "\n")
    )
  }

# Use predictors from semiauto_check()
  x <- check_result$predictors

# 2. Run missing value summary using semiauto_missing()
  missing_info <- semiauto_missing(
    data = data,
    y = y,
    x = x,
    verbose = FALSE
  )

# 3. Select model variables
  selected_vars <- unique(c(y, x, time, id))
  model_data <- data[selected_vars]

  n_original <- nrow(model_data)

# 4. Handle missing values
  incomplete_rows <- !stats::complete.cases(model_data)
  n_removed <- sum(incomplete_rows)

  if (n_removed > 0 && na_action == "fail") {
    stop(
      "Missing values were found. Use na_action = 'omit' to remove incomplete rows."
    )
  }

  if (n_removed > 0 && na_action == "omit") {
    model_data <- model_data[stats::complete.cases(model_data), , drop = FALSE]
  }

  if (nrow(model_data) == 0) {
    stop("No complete rows remain after missing value handling.")
  }

# 5. Order data by id and/or time
  if (!is.null(id) && !is.null(time)) {
    ord <- order(model_data[[id]], model_data[[time]])
    model_data <- model_data[ord, , drop = FALSE]
  } else if (!is.null(time)) {
    ord <- order(model_data[[time]])
    model_data <- model_data[ord, , drop = FALSE]
  }

# 6. Create response vector and predictor matrix
  y_vector <- model_data[[y]]

  X <- as.matrix(model_data[x])
  X <- cbind(Intercept = 1, X)

# 7. Create result object
  result <- list(
    data = data,
    model_data = model_data,
    check_result = check_result,
    missing_info = missing_info,
    response = y,
    predictors = x,
    time = time,
    id = id,
    n_original = n_original,
    n_used = nrow(model_data),
    n_removed = n_removed,
    na_action = na_action,
    y = y_vector,
    X = X,
    call = match.call()
  )

  class(result) <- "semiauto_prepare"

  if (isTRUE(verbose)) {
    print(result)
  }

  invisible(result)
}
#' Print Method for semiauto_prepare
#'
#' @param x An object of class \code{semiauto_prepare}.
#' @param ... Additional arguments.
#'
#' @export
print.semiauto_prepare <- function(x, ...) {

  cat("\nsemiautoReg Data Preparation\n")
  cat("----------------------------\n")

  cat("Response variable :", x$response, "\n")
  cat("Predictors        :", paste(x$predictors, collapse = ", "), "\n")

  if (!is.null(x$time)) {
    cat("Time variable     :", x$time, "\n")
  }

  if (!is.null(x$id)) {
    cat("ID variable       :", x$id, "\n")
  }

  cat("Original rows     :", x$n_original, "\n")
  cat("Rows used         :", x$n_used, "\n")
  cat("Rows removed      :", x$n_removed, "\n")
  cat("NA action         :", x$na_action, "\n")
  cat("X matrix size     :", nrow(x$X), "x", ncol(x$X), "\n")

  invisible(x)
}
