#' Title
#'
#' @param params
#'
#' @return
#' @export
#'
#' @examples
show.params <- function(params) {
  max.chars <- max(sapply(names(params), nchar))
  for (ix.names in sort(names(params))) {
    prefix <- paste(array(' ', max.chars - nchar(ix.names)), collapse = '')
    if (is.vector(params[[ix.names]]) && length(params[[ix.names]]) == 1) {
      if (is.list(params[[ix.names]])) {
        flog.info('  %s%s: %s', prefix, ix.names, paste(params[[ix.names]], collapse = ', '))
      } else if (is.character(params[[ix.names]])) {
        flog.info('  %s%s: %s', prefix, ix.names, params[[ix.names]])
      } else if (is.infinite(params[[ix.names]])) {
        flog.info('  %s%s: %f', prefix, ix.names, params[[ix.names]])
      } else if (is.logical(params[[ix.names]])) {
        flog.info('  %s%s: %s', prefix, ix.names, params[[ix.names]])
      } else if (is.integer(params[[ix.names]])) {
        flog.info('  %s%s: % 11d', prefix, ix.names, params[[ix.names]])
      } else {
        flog.info('  %s%s: % 11.3f', prefix, ix.names, params[[ix.names]])
      }
    } else if (is.vector(params[[ix.names]])) {
      flog.info('  %s%s: %s', prefix, ix.names, paste(params[[ix.names]], collapse = ', '))
    } else {
      flog.info('  %s%s: (I do not know how to display this)', prefix, ix.names)
    }
  }
}
