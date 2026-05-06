test_that("semiauto_missing works with complete data", {
  result <- semiauto_missing(
    data = mtcars,
    y = "mpg",
    x = c("wt", "hp", "qsec"),
    verbose = FALSE
  )

  expect_s3_class(result, "semiauto_missing")
  expect_false(result$has_missing)
  expect_equal(result$n_complete_rows, nrow(mtcars))
  expect_equal(result$n_incomplete_rows, 0)
})

test_that("semiauto_missing detects missing values", {
  test_data <- mtcars
  test_data$mpg[1:3] <- NA
  test_data$hp[5] <- NA

  result <- semiauto_missing(
    data = test_data,
    y = "mpg",
    x = c("wt", "hp", "qsec"),
    verbose = FALSE
  )

  expect_true(result$has_missing)
  expect_true("mpg" %in% result$variables_with_missing)
  expect_true("hp" %in% result$variables_with_missing)
  expect_equal(result$n_incomplete_rows, 4)
})

test_that("semiauto_missing gives error for missing variable names", {
  expect_error(
    semiauto_missing(
      data = mtcars,
      y = "mpg",
      x = c("wt", "not_exist"),
      verbose = FALSE
    )
  )
})

test_that("semiauto_missing checks all variables when y and x are NULL", {
  result <- semiauto_missing(
    data = mtcars,
    verbose = FALSE
  )

  expect_s3_class(result, "semiauto_missing")
  expect_equal(length(result$selected_variables), ncol(mtcars))
})
