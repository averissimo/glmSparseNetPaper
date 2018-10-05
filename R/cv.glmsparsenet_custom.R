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
my.cv.glmnet <- function(alpha,
                         xdata,
                         ydata,
                         penalty.factor,
                         nlambda           = 1000,
                         nfolds            = 10,
                         lambda.min.ratio  = .001,
                         seed              = 1985,
                         ...) {
  #
  set.seed(seed)
  foldid <- loose.rock::balanced.cv.folds(ydata$status, nfolds = nfolds)$output
  #
  new.model <- glmnet.mclapply::cv.glmnet(xdata, Surv(ydata$time, ydata$status * 1),
                         family           ='cox',
                         penalty.factor   = penalty.factor,
                         #
                         alpha            = alpha,
                         foldid           = foldid,
                         nlambda          = nlambda,
                         lambda.min.ratio = lambda.min.ratio,
                         standardize      = FALSE,
                         ...)

  if (length(unique(new.model$nzero)) < 5) {
    lambda.first <- new.model$lambda[1]

    lambda.nrow <- 150
    lambda.ncol <- 3

    lambda <- (matrix(rep(1 / 10^seq(0, lambda.ncol - 1, 1), lambda.nrow),
                      nrow = lambda.nrow, byrow = TRUE) *
                 array(lambda.first * seq(1/lambda.nrow, 1, 1/lambda.nrow),
                       dim = c(lambda.nrow, lambda.ncol))) %>%
      as.vector %>% unique %>% sort(decreasing = TRUE)

    new.model <- glmnet.mclapply::cv.glmnet(xdata, Surv(ydata$time, ydata$status * 1),
                                            family           ='cox',
                                            penalty.factor   = penalty.factor,
                                            #
                                            alpha            = alpha,
                                            foldid           = foldid,
                                            lambda           = lambda,
                                            standardize      = FALSE,
                                            ...)
  }
  #
  new.coef          <- glmnet::coef.cv.glmnet(new.model, s = 'lambda.min')[,1]
  new.target.lambda <- new.model$lambda.min
  #
  return(list(model = new.model, lambda = new.target.lambda, coef = new.coef))
}
