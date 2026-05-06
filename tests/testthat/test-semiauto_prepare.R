test_that("semiauto_prepare works with complete data", {
  result <- semiauto_prepare(
    data = mtcars,
    y = "mpg",
    x = c("wt", "hp", "qsec"),
    verbose = FALSE
  )

  expect_s3_class(result, "semiauto_prepare")
  expect_equal(result$response, "mpg")
  expect_equal(result$predictors, c("wt", "hp", "qsec"))
  expect_equal(result$n_used, nrow(mtcars))
  expect_equal(ncol(result$X), 4)
})

test_that("semiauto_prepare omits missing values", {
  test_data <- mtcars
  test_data$mpg[1:3] <- NA
  test_data$hp[5] <- NA

  result <- semiauto_prepare(
    data = test_data,
    y = "mpg",
    x = c("wt", "hp", "qsec"),
    na_action = "omit",
    verbose = FALSE
  )

  expect_s3_class(result, "semiauto_prepare")
  expect_equal(result$n_original, 32)
  expect_equal(result$n_used, 28)
  expect_equal(result$n_removed, 4)
  expect_equal(nrow(result$X), 28)
})

test_that("semiauto_prepare fails when missing values exist and na_action is fail", {
  test_data <- mtcars
  test_data$mpg[1] <- NA

  expect_error(
    semiauto_prepare(
      data = test_data,
      y = "mpg",
      x = c("wt", "hp"),
      na_action = "fail",
      verbose = FALSE
    )
  )
})

test_that("semiauto_prepare uses all predictors when x is NULL", {
  result <- semiauto_prepare(
    data = mtcars,
    y = "mpg",
    verbose = FALSE
  )

  expect_s3_class(result, "semiauto_prepare")
  expect_equal(result$response, "mpg")
  expect_equal(length(result$predictors), ncol(mtcars) - 1)
  expect_equal(nrow(result$X), nrow(mtcars))
})
