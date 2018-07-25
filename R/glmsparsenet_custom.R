#' Call glmSparseNet model with Cox regression
#'
#' This is an auxiliary method for the analysis. It uses the run.cache method
#' to cache the results and speed up analysis
#'
#' @param xdata input matrix
#' @param ydata Surival data.dataframe with time and status columns
#' @param network network to use, could be a matrix, a degree vector or a string, see ?glmSparseNet
#' @param target.vars number of variables to target in calculations
#' @param xdata.digest sha256 checksum of xdata
#' @param cache.prefix prefix for the cache files
#' @param force.recalc force cache to be recalculated
#' @param ... additional parameters for glmnet
#'
#' @return a glmnet model
#' @export
glmSparseNet.cox <- function(xdata,
                             ydata,
                             target.vars,
                             alpha            = 1,
                             network          = 'correlation',
                             xdata.digest     = NULL,
                             cache.prefix     = 'glmSparseNet.cache',
                             force.recalc     = FALSE,
                             lambda.min.ratio = 0.001,
                             ...) {

  result <- glmSparseNet.cox_(xdata,
                              ydata,
                              target.vars,
                              network,
                              xdata.digest,
                              cache.prefix,
                              force.recalc,
                              alpha            = alpha,
                              nlambda          = 500,
                              lambda.min.ratio = lambda.min.ratio,
                              ...)

  if (sum(result$coef != 0) < (target.vars - 5)) {
    result <- glmSparseNet.cox_(xdata            = xdata,
                                ydata            = ydata,
                                target.vars      = target.vars,
                                network          = network,
                                alpha            = alpha,
                                nlambda          = 1000,
                                xdata.digest     = xdata.digest,
                                cache.prefix     = cache.prefix,
                                lambda.min.ratio = .00001,
                                ...)
  }
  return(result)
}

#' Call glmSparseNet model with Cox regression
#'
#' This is an auxiliary method for the analysis. It uses the run.cache method
#' to cache the results and speed up analysis
#'
#' @param xdata input matrix
#' @param ydata Surival data.dataframe with time and status columns
#' @param network network to use, could be a matrix, a degree vector or a string, see ?glmSparseNet
#' @param target.vars number of variables to target in calculations
#' @param xdata.digest sha256 checksum of xdata
#' @param cache.prefix prefix for the cache files
#' @param force.recalc force cache to be recalculated
#'
#' @return a glmnet model
glmSparseNet.cox_ <- function(xdata,
                              ydata,
                              target.vars,
                              network          = 'correlation',
                              xdata.digest     = NULL,
                              cache.prefix     = 'glmSparseNet.cache',
                              force.recalc     = FALSE,
                              ...) {

  new.model <- run.cache(glmSparseNet,
                         xdata,
                         Surv(ydata$time, ydata$status),
                         family           ='cox',
                         standardize      = F,
                         network          = network,
                         ...,
                         #
                         force.recalc = force.recalc,
                         cache.prefix = gsub('_', '.', gsub('_old', '', prefix)),
                         cache.digest = list(xdata.digest))

  if (any(new.model$df == target.vars)) {
    var.ix <- which(new.model$df == target.vars)
  } else if (any(new.model$df > target.vars)) {
    new.target <- min(new.model$df[new.model$df > target.vars])
    var.ix <- which(new.model$df == new.target)
  } else {
    new.target <- max(new.model$df[new.model$df < target.vars])
    var.ix <- which(new.model$df == new.target)
  }
  new.target.lambda <- new.model$lambda[var.ix[sort(var.ix, decreasing = T, index.return = T)$ix[1]]]
  new.coef          <- as.vector(coef(new.model, s = new.target.lambda))
  names(new.coef)   <- colnames(xdata[,xdata.ix])
  return(list(model = new.model, lambda = new.target.lambda, coef = new.coef))
}
