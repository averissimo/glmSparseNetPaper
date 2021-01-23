#' Title
#'
#' @param best.ones
#' @param xdata
#' @param ydata
#' @param ntimes.results
#'
#' @return
#' @export
#'
#' @examples
best.model.print <- function(filter.model, best.ones, xdata, ydata, ntimes.results) {
  best.ix    <- (best.ones %>% filter(model.name == filter.model))$ix
  best.model <- (best.ones %>% filter(model.name == filter.model))$model

  colnames(xdata) %>% {
    .[which(as.logical(ntimes.results[[best.ix]]$coefs[[best.model]] != 0))]
  } %>% geneNames() %>% { .$external_gene_name } %>%
    sort %>%
    paste(collapse = ', ') %>% {
      flog.info('Coefs. list', ., capture = TRUE)
    }

  test.ix          <- ntimes.results[[best.ix]]$ixs$test
  coef.list        <- list(ntimes.results[[best.ix]]$coefs[[best.model]])
  names(coef.list) <- filter.model
  cat('\n\n')
  cat('#### Test set {-}')
  cat('\n\n')
  print(
    separate2GroupsCox(coef.list, xdata[test.ix, ], ydata[test.ix, ],
                       plot.title = paste0(prettify.labels(best.model), ' (test set)'),
                       break.x.by = 12)$plot
  )
  cat('\n\n')
  cat('#### Full set {-}')
  cat('\n\n')
  print(
    separate2GroupsCox(coef.list, xdata, ydata,
                       plot.title = paste0(prettify.labels(best.model), ' (full data set)'),
                       break.x.by = 12)$plot
  )
  cat('\n\n')
}
