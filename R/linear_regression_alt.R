#' Implements simple linear regression by gradient descent
#'
#' @param dat A data frame
#' @param response The name of a response variable in the data frame (unquoted)
#' @param explanatory The name of the explanatory variable in the data frame (unquoted)
#'
#' @return A data frame of coefficients
#'
#' @import dplyr
#'
#' @export
slr_gd <- function(dat, response, explanatory){

  ### Compute coefficients by gradient descent
  ### Return a data frame of the same form as in the `simple_linear_regression`
  x <- dat %>% pull({{explanatory}})
  y <- dat %>% pull({{response}})

  explan_name <- dat %>%
    select({{explanatory}}) %>%
    names()
  m = 0
  c = 0
  L = .00001
  yhat =(m*x+c)
  MSEold = sum((y - yhat)^2) / length(x)
  MSEnew = 10000000000
  while(abs(MSEnew - MSEold) != 0){
    Dm = sum(2 * (y-yhat) * (-1*x))/length(x)
    Dc = sum(x*(y - yhat))  * -2/length(x)
    MSEold = MSEnew
    m = m - L*Dm
    c = c - L*Dc
    yhat =(m*x+c)
    MSEnew = sum((y - yhat)^2) / length(x)
    if (abs(MSEnew - MSEold) <.00001){
      results <- tibble::tibble(
        Intercept = c,
        Slope = m
      )
      names(results)[2] <- explan_name
      return(results)
    }
  }

}


#' Implements linear regression with many predictors by gradient descent
#'
#' This function computes coefficients for multiple regression by gradient descent
#' All columns of the provided data frame are used as predictors, except the
#' one specified as a response.
#'
#' No interaction terms are included.
#'
#'
#' @param dat A data frame
#' @param response The name of a response variable in the data frame (unquoted)
#'
#' @return A data frame of coefficients
#'
#' @import dplyr
#'
#'@export
mlr_gd <- function(dat, response) {
  x = as.matrix(select(dat, -{{response}}))
  x = cbind(1, x)
  m = runif(ncol(x), 0, 1)
  y = as.matrix(select(dat, {{response}}))
  L = .00001
  explan_name <- dat %>%
    select(-{{response}}) %>%
    names()
  explan_name = c("Intercept", explan_name)
  yhat = x %*% m
  MSEold = sum((y - yhat)^2) / length(x)
  MSEnew = 10000000000
  converged = FALSE
  while(converged == FALSE){
    Dm = 1/nrow(x) * t(x) %*% (yhat - y)
    MSEold = MSEnew
    m = m - L*Dm
    yhat = x %*% m
    MSEnew = sum((y - yhat)^2) / length(x)
    if (abs(MSEnew - MSEold) <= .0001){
      results = m
      converged = TRUE
      results = as.data.frame(t(m))
      names(results) = explan_name
      return(results)
    }
  }
  ### Compute coefficients by gradient descent
  ### Return a data frame of the same form as in the `multiple_linear_regression`


}
#' Implements linear regression with many predictors by matrix decomposition
#'
#' This function computes coefficients for multiple regression by QR matrix decomposition
#' All columns of the provided data frame are used as predictors, except the
#' one specified as a response.
#'
#' No interaction terms are included.
#'
#'
#' @param dat A data frame
#' @param response The name of a response variable in the data frame (unquoted)
#'
#' @return A data frame of coefficients
#'
#' @import dplyr
#'
#'@export
mlr_qr <- function(dat, response) {
  x = as.matrix(select(dat, -{{response}}))
  x = cbind(1, x)
  y = as.matrix(select(dat, {{response}}))
  QR = qr(x)
  Q = qr.Q(QR)
  R = qr.R(QR)
  results = backsolve(R, crossprod(Q,y))
  ### Compute coefficients by QR decomposition
  ### Return a data frame of the same form as in the `multiple_linear_regression`
  explan_name <- dat %>%
    select(-{{response}}) %>%
    names()
  explan_name = c("Intercept", explan_name)
  results = as.data.frame(t(results))
  names(results) = explan_name
  return(results)

}
