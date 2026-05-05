#' Check Input Data for semiautoReg
#'
#' This function checks whether the input data is suitable for semiparametric
#' regression modeling in the semiautoReg package.
#'
#' @param data A data frame containing the response and predictor variables.
#' @param y A character string specifying the response variable name.
#' @param x A character vector specifying predictor variable names. If NULL,
#' all columns except y, time, and id will be used as predictors.
#' @param time Optional character string specifying the time variable.
#' @param id Optional character string specifying the subject or group ID variable.
#' @param allow_missing Logical. If FALSE, missing values are reported as a problem.
#' @param min_n Minimum required number of observations.
#' @param min_unique Minimum number of unique values required for the response variable.
#' @param verbose Logical. If TRUE, prints the checking result.
#'
#' @return An object of class \code{semiauto_check}.
#'
#' @examples
#' semiauto_check(mtcars, y = "mpg", x = c("wt", "hp", "qsec"))
#'
#' @export
#'
#'
semiauto_check <- function(data,
                           y,
                           x = NULL,
                           time = NULL,
                           id = NULL,
                           allow_missing = FALSE,
                           min_n = 10,
                           min_unique = 3,
                           verbose = TRUE) {
  problems <- character()
  notes <- character()

# 1. Check data type
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }

  if (nrow(data) == 0) {
    stop("`data` has no rows.")
  }

  if (ncol(data) == 0) {
    stop("`data` has no columns.")
  }

  if (anyDuplicated(names(data)) > 0) {
    stop("`data` contains duplicated column names.")
  }

# 2. Check response variable
  if (missing(y) || !is.character(y) || length(y) != 1) {
    stop("`y` must be a single character string, for example y = 'mpg'.")
  }

  if (!(y %in% names(data))) {
    stop("Response variable `y` was not found in `data`.")
  }

# 3. Check optional time and id variables
  if (!is.null(time)) {
    if (!is.character(time) || length(time) != 1) {
      stop("`time` must be a single character string.")
    }

    if (!(time %in% names(data))) {
      stop("`time` variable was not found in `data`.")
    }
  }

# 4. Define predictors
  if (is.null(x)) {
    x <- setdiff(names(data), c(y, time, id))
  }

  if (!is.character(x)) {
    stop("`x` must be a character vector of predictor variable names.")
  }

  if (length(x) == 0) {
    stop("No predictor variables were selected.")
  }

  missing_x <- setdiff(x, names(data))

  if (length(missing_x) > 0) {
    stop(
      "The following predictor variables were not found in `data`: ",
      paste(missing_x, collapse = ", ")
    )
  }

  model_cols <- c(y, x)

# 5. Check sample size
  if (nrow(data) < min_n) {
    problems <- c(
      problems,
      paste0("The number of observations is too small. Minimum required: ", min_n, ".")
    )
  }

# 6. Check numeric variables
  numeric_check <- vapply(data[model_cols], is.numeric, logical(1))
  non_numeric_cols <- names(numeric_check)[!numeric_check]

  if (!is.numeric(data[[y]])) {
    problems <- c(
      problems,
      paste0("Response variable `", y, "` must be numeric.")
    )
  }

  non_numeric_x <- intersect(non_numeric_cols, x)

  if (length(non_numeric_x) > 0) {
    problems <- c(
      problems,
      paste0(
        "The following predictor variables are not numeric: ",
        paste(non_numeric_x, collapse = ", ")
      )
    )
  }

# 7. Check missing values
  missing_summary <- data.frame(
    variable = model_cols,
    n_missing = vapply(data[model_cols], function(z) sum(is.na(z)), integer(1)),
    pct_missing = round(
      vapply(data[model_cols], function(z) mean(is.na(z)) * 100, numeric(1)),
      2
    ),
    row.names = NULL
  )

  vars_with_missing <- missing_summary$variable[missing_summary$n_missing > 0]

  if (length(vars_with_missing) > 0) {
    msg <- paste0(
      "Missing values were found in: ",
      paste(vars_with_missing, collapse = ", ")
    )

    if (isFALSE(allow_missing)) {
      problems <- c(problems, msg)
    } else {
      notes <- c(notes, msg)
    }
  }

# 8. Check infinite values
  numeric_cols <- model_cols[numeric_check]

  non_finite_summary <- data.frame(
    variable = numeric_cols,
    n_non_finite = vapply(
      data[numeric_cols],
      function(z) sum(!is.finite(z) & !is.na(z)),
      integer(1)
    ),
    row.names = NULL
  )

  vars_with_non_finite <- non_finite_summary$variable[
    non_finite_summary$n_non_finite > 0
  ]

  if (length(vars_with_non_finite) > 0) {
    problems <- c(
      problems,
      paste0(
        "Infinite or non-finite values were found in: ",
        paste(vars_with_non_finite, collapse = ", ")
      )
    )
  }

# 9. Check constant variables
  unique_counts <- vapply(
    data[numeric_cols],
    function(z) length(unique(z[!is.na(z) & is.finite(z)])),
    integer(1)
  )

  constant_cols <- names(unique_counts)[unique_counts <= 1]

  if (y %in% constant_cols) {
    problems <- c(
      problems,
      paste0("Response variable `", y, "` is constant.")
    )
  }

  constant_x <- intersect(constant_cols, x)

  if (length(constant_x) > 0) {
    problems <- c(
      problems,
      paste0(
        "The following predictor variables are constant: ",
        paste(constant_x, collapse = ", ")
      )
    )
  }

# 10. Check response unique values
  if (is.numeric(data[[y]])) {
    y_unique <- length(unique(data[[y]][!is.na(data[[y]]) & is.finite(data[[y]])]))

    if (y_unique < min_unique) {
      problems <- c(
        problems,
        paste0(
          "Response variable `", y,
          "` has too few unique values. Minimum required: ",
          min_unique, "."
        )
      )
    }
  }

# 11. Check duplicated rows in model variables
  duplicated_rows <- sum(duplicated(data[model_cols]))

  if (duplicated_rows > 0) {
    notes <- c(
      notes,
      paste0(duplicated_rows, " duplicated rows were found in the selected variables.")
    )
  }

# 12. Create result object
  result <- list(
    ok = length(problems) == 0,
    data = data,
    response = y,
    predictors = x,
    time = time,
    id = id,
    n_obs = nrow(data),
    n_predictors = length(x),
    missing_summary = missing_summary,
    non_finite_summary = non_finite_summary,
    non_numeric_columns = non_numeric_cols,
    constant_columns = constant_cols,
    duplicated_rows = duplicated_rows,
    problems = problems,
    notes = notes,
    call = match.call()
  )

  class(result) <- "semiauto_check"

  if (isTRUE(verbose)) {
    print(result)
  }

  invisible(result)
}


#' Print Method for semiauto_check
#'
#' @param x An object of class \code{semiauto_check}.
#' @param ... Additional arguments.
#'
#' @export
print.semiauto_check <- function(x, ...) {

  cat("\nsemiautoReg Data Check\n")
  cat("----------------------\n")

  cat("Status      : ")
  if (isTRUE(x$ok)) {
    cat("PASSED\n")
  } else {
    cat("FAILED\n")
  }

  cat("Observations:", x$n_obs, "\n")
  cat("Response    :", x$response, "\n")
  cat("Predictors  :", paste(x$predictors, collapse = ", "), "\n")
  cat("N Predictors:", x$n_predictors, "\n")

  if (!is.null(x$time)) {
    cat("Time variable:", x$time, "\n")
  }

  if (!is.null(x$id)) {
    cat("ID variable  :", x$id, "\n")
  }

  if (length(x$problems) > 0) {
    cat("\nProblems:\n")
    for (p in x$problems) {
      cat("- ", p, "\n", sep = "")
    }
  }

  if (length(x$notes) > 0) {
    cat("\nNotes:\n")
    for (n in x$notes) {
      cat("- ", n, "\n", sep = "")
    }
  }

  cat("\nMissing Summary:\n")
  print(x$missing_summary, row.names = FALSE)

  invisible(x)
}
