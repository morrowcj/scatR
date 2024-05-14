#' Logit function
#'
#' @description
#' Convert a probability into a quantile of the logistic distribution
#'
#' @param p numeric, the probability to convert
#'
#' @details
#' The equation to convert the probability p into a quantile q is
#' q = log(p / (1 - p))
#'
#' @return numeric, quantile of the logistic regression corresponding to \code{p}
#'
#' @export
logit <- function(p) {log(p / (1 - p))}

#' Inverse logit function
#'
#' @description
#' Convert a quantile of the logistic distribution into a probability
#'
#' @param q numeric, the quantile to convert
#'
#' @details
#' The equation to convert the logistic quantile q into a probability p is
#' p = 1 / (1 + exp(-q))
#'
#'
#' @return numeric, the probability corresponding to \code{q}
inverse_logit <- function(q) {1 / (1 + exp(-q))}
