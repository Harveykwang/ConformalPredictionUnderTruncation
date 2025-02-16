#' Generate the product limiting estimator for survival functions, given a, y, and label, d.
#'
#' @export

survPL <- function(y, d){
  tFail <- sort(unique(y[d]))
  nDeat <- apply(outer(y, tFail, FUN = "=="), 2, sum)
  nRisk <- apply(outer(y, tFail, FUN = ">="), 2, sum)
  nRisk[length(nRisk)] <- nRisk[length(nRisk)] + 1
  surv <- c(1, cumprod(1 - nDeat / nRisk))
  function(t) sapply(t, function(t) surv[sum(tFail <= t) + 1])
}

#' Generate the product limiting estimator for survival functions, given a, y, and label, d.
#'
#' @export

survTrunPL <- function(a, y, d){
  tFail <- sort(unique(y[d]))
  nDeat <- apply(outer(y, tFail, FUN = "=="), 2, sum)
  nRisk <- apply(outer(y, tFail, FUN = ">=") & outer(a, tFail, FUN = "<="), 2, sum)
  surv <- c(1, cumprod(1 - nDeat / nRisk))
  function(t) sapply(t, function(t) surv[sum(tFail <= t) + 1])
}

#' Generate the product limiting estimator for survival functions, given a, y, and label, d.
#'
#' @export

survTrunPLNoCens <- function(a, y){
  tFail <- sort(unique(y))
  nDeat <- apply(outer(y, tFail, FUN = "=="), 2, sum)
  nRisk <- apply(outer(y, tFail, FUN = ">=") & outer(a, tFail, FUN = "<="), 2, sum)
  surv <- c(1, cumprod(1 - nDeat / nRisk))
  function(t) sapply(t, function(t) surv[sum(tFail <= t) + 1])
}

#' Generate the product limiting estimator for survival functions, given a, y, and label, d.
#'
#' @export

weightFunFact <- function(y, a, d = NULL){
  if(is.null(d)) {
    weight <- 1 / survTrunPLNoCens(a, y)(a)
    function(t) sapply(t, function(t) 1 / sum(weight[a <= t]))
  } else {
    survCen <- survPL(y = y - a, d = !d)
    weight <- 1 / survTrunPL(a, y, d)(a)
    function(t) sapply(t, function(t) 1 / sum(survCen(t - a[a <= t]) * weight[a <= t]))
  }
}
