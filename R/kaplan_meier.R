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
    #is.network <- sub('\\..*', '', is.network)
    is.network <- gsub('\\..*$', '', is.network)
    type <- if (is.network == 'glmnet') {
      subtitle <- 'Classic Elastic Net model'
      ''
    } else if (is.network == 'degree_log') {
      subtitle <- 'DegreeCox: Hubs are promoted -- log'
      'Promotes high degree with '
    } else if (is.network == 'degree_new') {
      subtitle <- 'HubCox: Hubs are promoted'
      'Promotes high degree with '
    } else if (is.network == 'degree_old') {
      subtitle <- 'DegreeCox: Hubs are promoted -- old'
      'Promotes high degree with '
    } else if (is.network == 'orphan') {
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
my.draw.kaplan <- function(coef.l, plot.title, xdata, ydata) {
  return(glmSparseNet::draw.kaplan(coef.l, xdata, ydata,
                     plot.title = plot.title, legend.outside = F, ylim = c(0,1)))
}
