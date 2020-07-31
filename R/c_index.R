# #' Title
# #'
# #' @param all.pairs
# #' @param fitted.risk
# #' @param ydata
# #' @param n.cores
# #'
# #' @return
# #' @export
# c.index.fun <- function(all.pairs, fitted.risk, ydata, n.cores = params$n.cores) {
#   unlist(mclapply(seq(ncol(all.pairs)), function(ix) {
#     ix.1 <- all.pairs[1,ix]
#     ix.2 <- all.pairs[2,ix]
#     return(my.c.index.cmp(fitted.risk[ix.1], fitted.risk[ix.2],
#                           ydata[ix.1,1], ydata[ix.2,1],
#                           ydata[ix.1,2], ydata[ix.2,2]))
#   }, mc.cores = n.cores))
# }


#' Calculate c-index
#'
#' @param coef.v
#' @param xdata
#' @param ydata
#' @param model.name
#' @param show.message
#' @param n.cores
#'
#' @return
#' @export
fit.risk <- function(coef.v, xdata, ydata, model.name, show.message = FALSE, n.cores = 1) {
  # fitted.risk  <- as.vector(predict(models$glmnet, newx = xdata.train[,xdata.ix], s = lambdas$glmnet, type = 'response'))
  fitted.risk <- exp(as.vector(xdata %*% coef.v))
  c.res <- loose.rock::run.cache(survcomp::concordance.index,
                                 fitted.risk, ydata[,'time'],
                                 ydata[,'status'] * 1,
                                 method="noether", show.message = F,
                                 cache.prefix = 'c-index')

  if (is.na(c.res$c.index)) {
    warning(sprintf('Could not calculate c-index, probably to few events -- event: %d -- censored %d\n',
                    sum(ydata[,'status']), sum(!ydata[,'status'])))
    return(-1)
  }
  #if (show.message)
  #  flog.info(' * %s: %f (se = %f | pvalue = %f)', model.name, c.res$c.index, c.res$se, c.res$pvalue)
  return(c.res$c.index)
}
