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
#' km.name('This title', 0, 'glmnet--glmnet--0.20--39')
#' km.name('This title', 0, 'glmnet--glmnet--0.20--39', 30)
#' km.name('This title', .7, 'degree_new--orphan--0.20--42')
km.name <- function(title, alpha, is.network = NULL, nvars = NULL) {
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
    if (is.null(nvars)) {
      append.str <- sprintf(' (%s target variables to select)', prettify.labels(is.network, 'target.vars'))
    }
    else {
      append.str <- sprintf(' (%d selected variables with target of %s)', nvars, prettify.labels(is.network, 'target.vars'))
    }
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
  subtitle <- sprintf('%s%s', subtitle, append.str)
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
  ydata <- ydata %>% mutate(time = time / 365 * 12)
  return(glmSparseNet::separate2GroupsCox(coef.l, xdata, ydata, no.plot = no.plot,
                     plot.title = plot.title, legend.outside = F, ylim = c(0,1),
                     break.x.by = 12,
                     xlab = "Time in months",
                     surv.median.line = "hv"))
}
