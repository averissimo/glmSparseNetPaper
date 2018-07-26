#' Title
#'
#' @param prefix
#' @param penalty.factor
#' @param nfolds
#' @param lambda.min.ratio
#' @param show.messages
#'
#' @return
#' @export
#'
#' @examples
my.cv.glmnet <- function(prefix,
                         alpha,
                         xdata,
                         ydata,
                         penalty.factor,
                         xdata.digest      = NULL,
                         nfolds            = 10,
                         lambda.min.ratio  = .001,
                         n.cores           = 1,
                         seed              = 1985) {
  #
  set.seed(seed)
  foldid <- loose.rock::balanced.cv.folds(ydata$status, nfolds = nfolds)$output
  #
  new.model <- loose.rock::run.cache(cv.glmnet,
                                     xdata, Surv(ydata$time, ydata$status),
                                     family           ='cox',
                                     foldid           = foldid,
                                     alpha            = alpha,
                                     nlambda          = 1000,
                                     lambda.min.ratio = lambda.min.ratio,
                                     standardize      = F,
                                     penalty.factor   = penalty.factor,
                                     #
                                     mc.cores         = n.cores,
                                     force.recalc = F,
                                     cache.prefix = prefix,
                                     cache.digest = list(xdata.digest))
  new.coef <- glmnet::coef.cv.glmnet(new.model, s = 'lambda.min')[,1]
  new.target.lambda <- new.model$lambda.min
  #
  return(list(model = new.model, lambda = new.target.lambda, coef = new.coef))
}
