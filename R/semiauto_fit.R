#internal helper : SCAD derivative

scad_derivative=function(t,lambda,a=3.7) {
  t=abs(t)
  out=numeric(length(t))
  #region 1
  out[t<=lambda]=lambda
  #region 2
  idx=t>lambda&t<=a*lambda
  out[idx]=(a*lambda-t[idx])/(a-1)
  #region 3
  out[t>a*lambda]=0
  return(out)
}

#internal helper : LQA weights
scad_lqa_weight=function(t,lambda,a=3.7,eps=1e-6){
  t_abs=pmax(abs(t),eps)
  weights=scad_derivative(t_abs,lambda=lambda,a=a)/t_abs
  return(weights)
}

# =========================================
# Main function
# =========================================

#' Fit Semiparametric Regression with AR Errors
#'
#' @description
#' Fits a semiparametric regression model using
#' B-spline basis expansion, SCAD penalty,
#' and autoregressive error correction.
#'
#' @param prepared An object returned by \code{semiauto_prepare()}. This argument is required.
#' @param spline_var Character. Name of the predictor to be modeled nonparametrically using spline basis. This argument is required.
#' @param df Numeric. Degrees of freedom for the spline basis. Default is 5.
#' @param degree Numeric. Degree of the spline basis. Default is 3.
#' @param lambda Numeric. SCAD penalty parameter. Default is 0.1.
#' @param max_iter Numeric. Maximum number of iterations. Default is 50.
#' @param tol Numeric. Convergence tolerance. Default is 1e-6.
#'
#' @return A list containing fitted model results.
#'
#' @importFrom splines bs
#' @importFrom MASS ginv
#'
#' @examples
#' prep <- semiauto_prepare(
#'   data = mtcars,
#'   y = "mpg",
#'   x = c("wt", "hp", "qsec")
#' )
#'
#' fit <- semiauto_fit(
#'   prepared = prep,
#'   spline_var = "wt"
#' )
#'
#' fit
#'
#' @export

semiauto_fit <- function(prepared,
                         spline_var,
                         df = 5,
                         degree = 3,
                         lambda = 0.1,
                         max_iter = 50,
                         tol = 1e-6) {

  # =========================
  # 1. Basic checks
  # =========================
  if (!inherits(prepared, "semiauto_prepare")) {
    stop("prepared must be an object returned by semiauto_prepare().")
  }

  model_data <- prepared$model_data
  y_name <- prepared$response
  x <- prepared$predictors

  if (!spline_var %in% x) {
    stop("spline_var must be one of the predictors in prepared$predictors.")
  }

  if (!is.numeric(model_data[[spline_var]])) {
    stop("spline_var must be numeric.")
  }

  # =========================
  # 2. Prepare variables
  # =========================
  Y <- as.numeric(model_data[[y_name]])

  linear_vars <- setdiff(x, spline_var)

  if (length(linear_vars) > 0) {
    X <- stats::model.matrix(
      stats::as.formula(
        paste("~", paste(linear_vars, collapse = " + "))
      ),
      data = model_data
    )
  } else {
    X <- matrix(1, nrow = nrow(model_data), ncol = 1)
    colnames(X) <- "(Intercept)"
  }

  Z <- splines::bs(
    model_data[[spline_var]],
    df = df,
    degree = degree,
    intercept = TRUE
  )

  W <- cbind(X, Z)

  p_x <- ncol(X)
  p_z <- ncol(Z)

  # =========================
  # 3. Initial fit
  # =========================
  beta_theta <- stats::lm.fit(W, Y)$coefficients
  beta_theta[is.na(beta_theta)] <- 0

  beta_old <- beta_theta[1:p_x]
  theta_old <- beta_theta[(p_x + 1):(p_x + p_z)]

  gamma_old <- 0

  # =========================
  # 4. Iteration: SCAD + AR(1)
  # =========================
  for (iter in seq_len(max_iter)) {

    fitted_values <- as.vector(W %*% c(beta_old, theta_old))
    residuals <- Y - fitted_values

    if (length(residuals) > 2) {
      gamma_new <- stats::acf(
        residuals,
        plot = FALSE,
        lag.max = 1
      )$acf[2]
    } else {
      gamma_new <- 0
    }

    if (is.na(gamma_new)) {
      gamma_new <- 0
    }

    Y_star <- Y[-1] - gamma_new * Y[-length(Y)]

    W_star <- W[-1, , drop = FALSE] -
      gamma_new * W[-nrow(W), , drop = FALSE]

    weights_theta <- scad_lqa_weight(
      theta_old,
      lambda = lambda
    )

    penalty <- diag(c(rep(0, p_x), weights_theta))

    A <- crossprod(W_star) + penalty
    b <- crossprod(W_star, Y_star)

    beta_theta_new <- tryCatch(
      solve(A, b),
      error = function(e) MASS::ginv(A) %*% b
    )

    beta_new <- as.numeric(beta_theta_new[1:p_x])
    theta_new <- as.numeric(beta_theta_new[(p_x + 1):(p_x + p_z)])

    diff <- max(abs(c(
      beta_new - beta_old,
      theta_new - theta_old,
      gamma_new - gamma_old
    )))

    beta_old <- beta_new
    theta_old <- theta_new
    gamma_old <- gamma_new

    if (diff < tol) {
      break
    }
  }

  # =========================
  # 5. Final fitted values
  # =========================
  coefficients <- c(beta_old, theta_old)
  fitted_values <- as.vector(W %*% coefficients)
  residuals <- Y - fitted_values

  selected_spline <- abs(theta_old) > 1e-6

  names(beta_old) <- colnames(X)
  names(theta_old) <- colnames(Z)
  names(coefficients) <- c(colnames(X), colnames(Z))

  result <- list(
    call = match.call(),
    y = y_name,
    x = x,
    spline_var = spline_var,
    linear_vars = linear_vars,
    coefficients = coefficients,
    beta = beta_old,
    theta = theta_old,
    gamma = gamma_old,
    fitted_values = fitted_values,
    residuals = residuals,
    selected_spline = selected_spline,
    df = df,
    degree = degree,
    lambda = lambda,
    iterations = iter,
    converged = iter < max_iter,
    model_data = model_data,
    X = X,
    Z = Z,
    W = W
  )

  class(result) <- "semiauto_fit"

  return(result)
}

#' Print semiauto_fit Object
#'
#' @param x Object from semiauto_fit().
#' @param ... Additional arguments.
#' @export

print.semiauto_fit <- function(x, ...) {

  cat("\n")
  cat("semiautoReg Model Fit\n")
  cat("----------------------\n")

  cat("Response variable : ", x$y, "\n", sep = "")
  cat("Spline variable  : ", x$spline_var, "\n", sep = "")
  cat("Iterations       : ", x$iterations, "\n", sep = "")
  cat("Converged        : ", x$converged, "\n", sep = "")
  cat("AR(1) gamma      : ", round(x$gamma, 4), "\n", sep = "")

  cat("\nLinear coefficients:\n")
  print(round(x$beta, 4))

  cat("\nSpline coefficients:\n")
  print(round(x$theta, 4))

  invisible(x)
}
