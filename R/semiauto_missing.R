#' Missing Value Summary for semiautoReg
#'
#' This function summarizes missing values in the selected variables.
#' It reports the number and percentage of missing values for each variable,
#' as well as the number of complete and incomplete rows.
#'
#' @param data A data frame.
#' @param y Optional character string specifying the response variable name.
#' @param x Optional character vector specifying predictor variable names.
#' If both y and x are NULL, all columns in the data will be checked.
#' @param verbose Logical. If TRUE, prints the missing value summary.
#'
#' @return An object of class \code{semiauto_missing}.
#'
#' @examples
#' semiauto_missing(mtcars)
#' semiauto_missing(mtcars, y = "mpg", x = c("wt", "hp", "qsec"))
#'
#' @export
semiauto_missing <- function(data,
                             y = NULL,
                             x = NULL,
                             verbose = TRUE) {

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

# 2. Define selected variables
  if (is.null(y) && is.null(x)) {
    selected_vars <- names(data)
  } else {
    selected_vars <- c(y, x)
  }

  selected_vars <- unique(selected_vars)

  if (!is.character(selected_vars)) {
    stop("`y` and `x` must be character variable names.")
  }

  missing_vars <- setdiff(selected_vars, names(data))

  if (length(missing_vars) > 0) {
    stop(
      "The following variables were not found in `data`: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  selected_data <- data[selected_vars]

# 3. Missing value summary by variable
  missing_summary <- data.frame(
    variable = selected_vars,
    n_missing = vapply(
      selected_data,
      function(z) sum(is.na(z)),
      integer(1)
    ),
    pct_missing = round(
      vapply(
        selected_data,
        function(z) mean(is.na(z)) * 100,
        numeric(1)
      ),
      2
    ),
    row.names = NULL
  )

# 4. Missing value summary by row
  row_missing_count <- rowSums(is.na(selected_data))

  n_complete_rows <- sum(row_missing_count == 0)
  n_incomplete_rows <- sum(row_missing_count > 0)

  pct_complete_rows <- round(n_complete_rows / nrow(selected_data) * 100, 2)
  pct_incomplete_rows <- round(n_incomplete_rows / nrow(selected_data) * 100, 2)

# 5. Variables with missing values
  variables_with_missing <- missing_summary$variable[
    missing_summary$n_missing > 0
  ]

# 6. Create result object
  result <- list(
    data = data,
    selected_variables = selected_vars,
    missing_summary = missing_summary,
    variables_with_missing = variables_with_missing,
    n_complete_rows = n_complete_rows,
    n_incomplete_rows = n_incomplete_rows,
    pct_complete_rows = pct_complete_rows,
    pct_incomplete_rows = pct_incomplete_rows,
    row_missing_count = row_missing_count,
    has_missing = length(variables_with_missing) > 0,
    call = match.call()
  )
  class(result) <- "semiauto_missing"

  if (isTRUE(verbose)) {
    print(result)
  }

  invisible(result)
}


#' Print Method for semiauto_missing
#'
#' @param x An object of class \code{semiauto_missing}.
#' @param ... Additional arguments.
#'
#' @export
print.semiauto_missing <- function(x, ...) {

  cat("\nsemiautoReg Missing Value Summary\n")
  cat("----------------------------------\n")

  cat("Selected variables :", paste(x$selected_variables, collapse = ", "), "\n")
  cat("Complete rows      :", x$n_complete_rows, "(", x$pct_complete_rows, "% )\n")
  cat("Incomplete rows    :", x$n_incomplete_rows, "(", x$pct_incomplete_rows, "% )\n")

  cat("Status             : ")

  if (isTRUE(x$has_missing)) {
    cat("Missing values found\n")
  } else {
    cat("No missing values found\n")
  }

  if (length(x$variables_with_missing) > 0) {
    cat("Variables with NA  :", paste(x$variables_with_missing, collapse = ", "), "\n")
  }

  cat("\nMissing Summary:\n")
  print(x$missing_summary, row.names = FALSE)

  invisible(x)
}
