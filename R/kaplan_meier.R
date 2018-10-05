#' Title
#'
#' @param title
#' @param alpha
#' @param is.network
#'
#' @return
#' @export
#'
#' @examples
#' km.name('This title', .6)
#' km.name('This title', 1)
#' km.name('This title', 0, 'glmnet')
#' km.name('This title', .7, 'orphan')
km.name <- function(title, alpha, is.network = NULL) {
  l1 <- sprintf('%.1f * L1', alpha)
  l2 <- sprintf('%.1f * L2', 1 - alpha)
  sep <- ' + '
  subtitle <- 'Classic elastic net'
  if (alpha == 1) {
    l1 <- 'L1'
    l2 <- ''
    sep <- ''
  } else if (alpha == 0) {
    l1 <- ''
    l2 <- 'L2'
    sep <- ''
  }
  type <- if (is.null(is.network)) {
    subtitle <- 'Classic Elastic Net model'
    ''
  } else {
    is.network <- prettify.labels(is.network, 'model')
    type <- if (is.network == 'Elastic Net') {
      subtitle <- 'Classic Elastic Net model'
      ''
    } else if (is.network == 'Degree_log') {
      subtitle <- 'DegreeCox: Hubs are promoted -- log'
      'Promotes high degree with '
    } else if (is.network == 'Hub') {
      subtitle <- 'HubCox: Hubs are promoted'
      'Promotes high degree with '
    } else if (is.network == 'Degree_old') {
      subtitle <- 'DegreeCox: Hubs are promoted -- old'
      'Promotes high degree with '
    } else if (is.network == 'Orphan') {
      subtitle <- 'OrphanCox: Low connected nodes are promoted'
      'Promotes low degree with '
    }
  }
  return(ggplot2::ggtitle(sprintf('%s: %s%s%s%s', title, type, l1, sep, l2), subtitle = subtitle))
}


#' Title
#'
#' @param coef.l
#' @param plot.title
#' @param xdata
#' @param ydata
#' @param xdata.ix
#'
#' @return
#' @export
my.draw.kaplan <- function(coef.l, xdata, ydata, plot.title = '', no.plot = FALSE) {
  return(glmSparseNet::separate2GroupsCox(coef.l, xdata, ydata, no.plot = no.plot,
                     plot.title = plot.title, legend.outside = F, ylim = c(0,1)))
}
