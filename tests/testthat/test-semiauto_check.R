test_that("semiauto_check works with valid data", {
  result <- semiauto_check(
    data = mtcars,
    y = "mpg",
    x = c("wt", "hp", "qsec"),
    verbose = FALSE
  )

  expect_s3_class(result, "semiauto_check")
  expect_true(result$ok)
  expect_equal(result$response, "mpg")
  expect_equal(result$n_predictors, 3)
})

test_that("semiauto_check detects missing response variable", {
  expect_error(
    semiauto_check(
      data = mtcars,
      y = "not_exist",
      x = c("wt", "hp"),
      verbose = FALSE
    )
  )
})

test_that("semiauto_check detects non numeric predictor", {
  test_data <- mtcars
  test_data$group <- as.character(test_data$cyl)

  result <- semiauto_check(
    data = test_data,
    y = "mpg",
    x = c("wt", "group"),
    verbose = FALSE
  )

  expect_false(result$ok)
  expect_true("group" %in% result$non_numeric_columns)
})

test_that("semiauto_check detects missing values", {
  test_data <- mtcars
  test_data$mpg[1] <- NA

  result <- semiauto_check(
    data = test_data,
    y = "mpg",
    x = c("wt", "hp"),
    verbose = FALSE
  )

  expect_false(result$ok)
  expect_true(length(result$problems) > 0)
})
