#' Conformal prediction inference for survival data with truncation
#'
#' The predInterval function produces conformal prediction intervals for survival time given the coverates.
#'
#' @param XTest a matrix of covariates with each row corresponding to an instance of covariate vector.
#' @param ATr truncation time.
#' @param YTr the observed event time.
#' @param XTr the observed covariate matrix
#' @param del the censoring indicator vector.
#' @param cens the censoring existence indicator.
#' @param alp miscoverage probability.
#' @param model  a model name to indicate which model will be used to fit the data.
#' @param tau the restricted maximum time.
#'
#' @export

predInterval <- function(XTest, ATr, YTr, XTr, del = NULL, cens = TRUE, alp = 0.1, model = "cox", tau = NULL) {

  pkgs <- c("survival", "flexsurv", "LTRCforests", "randomForestSRC", "pec", "glmnet")
  addPkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(addPkgs)) install.packages(addPkgs, repos='http://cran.us.r-project.org')
  suppressPackageStartupMessages(res <- lapply(X=pkgs, FUN=require, character.only=TRUE))

  if(!is.matrix(XTest)) XTest <- as.matrix(XTest, nrow = 1)
  p <- ncol(XTest)
  colnames(XTest) <- colnames(XTr) <- xName <- paste0("X", 1:p)
  nTest <- nrow(XTest)
  XTest <- as.data.frame(XTest)
  n <- length(ATr)
  if(!cens) {
    data <- as.data.frame(cbind(A = ATr, Y = YTr, del = rep(1, n), X = XTr))
  } else {
    data <- as.data.frame(cbind(A = ATr, Y = YTr, del = del, X = XTr))
  }

  # split the training data
  groupSize <- ceiling(n / 2)
  groupLable <- sample(n, groupSize)
  dataTrai <- data[groupLable, ]
  dataPred <- data[-groupLable, ]

  if(is.null(tau)) {
    failExac <- dataPred[dataPred$del == 1, ]
    if(cens) {
      weigFun <- weightFunFact(y = dataTrai$Y, a = dataTrai$A, d = dataTrai$del)
    } else {
      weigFun <- weightFunFact(y = dataTrai$Y, a = dataTrai$A)
    }
    tRang <- range(failExac$Y)
    tGrid <- seq(tRang[1], tRang[2], length.out = 1e3) ###  # take grid points of the range of T
  } else {
    dataPred$Y <- pmin(dataPred$Y, tau)
    index <- (dataPred$Y == tau | dataPred$del == 1) & dataPred$A < tau
    failExac <- dataPred[index, ]
    if(cens) {
      weigFun <- weightFunFact(y = dataTrai$Y, a = dataTrai$A, d = dataTrai$del)
    } else {
      weigFun <- weightFunFact(y = dataTrai$Y, a = dataTrai$A)
    }
    tRang <- range(failExac$Y)
    tGrid <- seq(tRang[1], min(tRang[2], tau), length.out = 1e3) ###  # take grid points of the range of T
  }
  weigGrid <- weigFun(tGrid)
  weig <- weigFun(failExac$Y)

  if (model == "cox") {
    fml <- as.formula(paste0("Surv(time = A, time2 = Y, event = del) ~ ", paste(xName, collapse = "+")))
    modelFit <- coxph(fml, data = dataTrai)
    newData <- failExac
    newData$A <- 0
    survPred <- predict(modelFit, newdata = newData, type = "survival")
    score <- abs(survPred - .5)
    survCurv <- survfit(modelFit, newdata = XTest)
    survGrid <- apply(rbind(1, survCurv$surv), 2, function(x) x[rowSums(outer(tGrid, survCurv$time, ">=")) + 1])
    scorGrid <- abs(survGrid - .5)
    weig <- weig[order(score)]
    score <- c(sort(score), 0.5)
    weigMat <- rbind(matrix(weig, nrow = length(weig), ncol = length(weigGrid)), weigGrid)
    weigMatNorm <- t(t(weigMat) / colSums(weigMat))
    quanRank <- apply(apply(weigMatNorm, 2, cumsum), 2, function(x) which.max(x >= 1 - alp))
    indi <- scorGrid <= score[quanRank]
    lowerInd <- apply(indi, 2, which.max)
    indi <- apply(indi, 2, rev)
    upperInd <- nrow(indi) + 1 - apply(indi, 2, which.max)
    lower <- tGrid[pmax(lowerInd - 1, 1)]
    upper <- tGrid[pmin(upperInd + 1, 1000)]
  }

  if (model == "lognormal") {
    fml <- as.formula(paste0("Surv(time = A, time2 = Y, event = del) ~ ", paste(xName, collapse = "+")))
    modelFit <- flexsurvreg(fml, data = dataTrai, dist = "lognormal")
    coe <- modelFit$coefficients
    nr <- nrow(failExac)
    survPred <- pnorm(log(failExac$Y), mean = coe[1] + apply(coe[3:(p + 2)] * t(failExac[ , xName]), 2, sum),
                      sd = exp(coe[2]), lower.tail = FALSE)
    score <- abs(survPred - .5)
    meanModel <- coe[1] + apply(coe[3:(p + 2)] * t(XTest), 2, sum)
    survGrid <- sapply(1:nTest, function(i) pnorm(log(tGrid), mean = meanModel[i], sd = exp(coe[2]), lower.tail = FALSE))
    scorGrid <- abs(survGrid - .5)
    weig <- weig[order(score)]
    score <- c(sort(score), 0.5)
    weigMat <- rbind(matrix(weig, nrow = length(weig), ncol = length(weigGrid)), weigGrid)
    weigMatNorm <- t(t(weigMat) / colSums(weigMat))
    quanRank <- apply(apply(weigMatNorm, 2, cumsum), 2, function(x) which.max(x >= 1 - alp))
    indi <- scorGrid <= score[quanRank]
    lowerInd <- apply(indi, 2, which.max)
    indi <- apply(indi, 2, rev)
    upperInd <- nrow(indi) + 1 - apply(indi, 2, which.max)
    lower <- tGrid[pmax(lowerInd - 1, 1)]
    upper <- tGrid[pmin(upperInd + 1, 1000)]
  }

  if (model == "weibull") {
    fml <- as.formula(paste0("Surv(time = A, time2 = Y, event = del) ~ ", paste(xName, collapse = "+")))
    modelFit <- flexsurvreg(fml, data = dataTrai, dist = "weibull", method = "Nelder-Mead")
    coe <- modelFit$coefficients
    nr <- nrow(failExac)
    survPred <- pweibull(failExac$Y, shape = exp(coe[1]),
                         scale = exp(coe[2]) * exp(apply(coe[3:(p + 2)] * t(failExac[ , xName]), 2, sum)), lower.tail = FALSE)
    score <- abs(survPred - .5)
    scaleModel <- exp(coe[2]) * exp(apply(coe[3:(p + 2)] * t(XTest), 2, sum))
    survGrid <- sapply(1:nTest, function(i) pweibull(tGrid, shape = exp(coe[1]), scale = scaleModel[i], lower.tail = FALSE))
    scorGrid <- abs(survGrid - .5)
    weig <- weig[order(score)]
    score <- c(sort(score), 0.5)
    weigMat <- rbind(matrix(weig, nrow = length(weig), ncol = length(weigGrid)), weigGrid)
    weigMatNorm <- t(t(weigMat) / colSums(weigMat))
    quanRank <- apply(apply(weigMatNorm, 2, cumsum), 2, function(x) which.max(x >= 1 - alp))
    indi <- scorGrid <= score[quanRank]
    lowerInd <- apply(indi, 2, which.max)
    indi <- apply(indi, 2, rev)
    upperInd <- nrow(indi) + 1 - apply(indi, 2, which.max)
    lower <- tGrid[pmax(lowerInd - 1, 1)]
    upper <- tGrid[pmin(upperInd + 1, 1000)]
  }

  if (model == "randomForest") {
    fml <- as.formula(paste0("Surv(time = A, time2 = Y, event = del) ~ ", paste(xName, collapse = "+")))
    modelFit <- ltrccif(fml, data = dataTrai, ntree = 100L, stepFactor = 3)
    survPred <- sapply(1:nrow(failExac),
                       function(i) predictProb(modelFit, newdata = failExac[i, ], time.eval = failExac$Y[i])$survival.probs)
    survPred[is.na(survPred)] <- 0
    score <- abs(survPred - .5)
    newData <- cbind(XTest, A = 0, del = 1, Y = Inf)
    survGrid <- predictProb(modelFit, newdata = newData, time.eval = tGrid)$survival.probs
    scorGrid <- abs(survGrid - .5)
    weig <- weig[order(score)]
    score <- c(sort(score), 0.5)
    weigMat <- rbind(matrix(weig, nrow = length(weig), ncol = length(weigGrid)), weigGrid)
    weigMatNorm <- t(t(weigMat) / colSums(weigMat))
    quanRank <- apply(apply(weigMatNorm, 2, cumsum), 2, function(x) which.max(x >= 1 - alp))
    indi <- scorGrid <= score[quanRank]
    lowerInd <- apply(indi, 2, which.max)
    indi <- apply(indi, 2, rev)
    upperInd <- nrow(indi) + 1 - apply(indi, 2, which.max)
    lower <- tGrid[pmax(lowerInd - 1, 1)]
    upper <- tGrid[pmin(upperInd + 1, 1000)]
  }

  if (model == "lasso") {
    fml <- Surv(time = dataTrai$A, time2 = dataTrai$Y, event = dataTrai$del)
    xMat <- as.matrix(dataTrai[ , xName])
    modelFit <- cv.glmnet(x = xMat, y = fml, family = "cox")
    survCurv <- survival::survfit(modelFit, x = xMat, y = fml, s = "lambda.1se", newx = data.matrix(failExac[, xName]))
    survMat <- rbind(1, survCurv$surv)
    survPred <- sapply(1:nrow(failExac),
                       function(i) survMat[tail(rank(c(survCurv$time, failExac$Y[i])), n = 1) , i])
    score <- abs(survPred - .5)

    survCurv <- survival::survfit(modelFit, x = xMat, y = fml, s = "lambda.1se", newx = data.matrix(XTest))
    survMat <- rbind(1, survCurv$surv); survTime <- survCurv$time
    tInd <- colSums(outer(survTime, tGrid, "<=")) + 1
    survGrid <- sapply(1:nTest, function(i) survMat[tInd, i])
    scorGrid <- abs(survGrid - .5)
    weig <- weig[order(score)]
    score <- c(sort(score), 0.5)
    weigMat <- rbind(matrix(weig, nrow = length(weig), ncol = length(weigGrid)), weigGrid)
    weigMatNorm <- t(t(weigMat) / colSums(weigMat))
    quanRank <- apply(apply(weigMatNorm, 2, cumsum), 2, function(x) which.max(x >= 1 - alp))
    indi <- scorGrid <= score[quanRank]
    lowerInd <- apply(indi, 2, which.max)
    indi <- apply(indi, 2, rev)
    upperInd <- nrow(indi) + 1 - apply(indi, 2, which.max)
    lower <- tGrid[pmax(lowerInd - 1, 1)]
    upper <- tGrid[pmin(upperInd + 1, 1000)]
  }

  return(list(lower = lower, upper = upper))

}
