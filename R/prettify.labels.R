#' Title
#'
#' @param label
#'
#' @return
#' @export
#'
#' @examples
#' prettify.labels('orphan--glmnet--0.20--83')
#' prettify.labels('orphan--0.01--103', 'base.model')
#' prettify.labels('orphan--degree_new--0.01--83', 'base.model')
#' prettify.labels('orphan--glmnet--0.01--83', 'base.model')
#' prettify.labels(c('orphan--glmnet--0.01--83', 'degree_new--orphan--0.01--83'))
#' prettify.labels('orphan--glmnet--0.01--83', 'string')
#' prettify.labels('orphan--glmnet--0.01--83', 'string.short')
#' prettify.labels('orphan--glmnet--0.01--83', 'string.base')
prettify.labels <- function(label, ret.value = 'string.short') {

  label.old <- label
  if (is.factor(label)) {
    label.old <- levels(label)
  }

  all.splitted <- strsplit(label.old, '--')

  str.t <- sapply(all.splitted, function(splitted) {
    if (length(splitted) == 4) {
      nvars <- splitted[4]
      alpha <- splitted[3]

      l.models <- gsub('glmnet', 'elastic net', splitted[seq_len(2)], ignore.case = TRUE)
      l.models <- gsub('degree_new', 'hub', l.models, ignore.case = TRUE)
      l.models <- proper(l.models)

      my.str.short <- sprintf('%s (base model: %s, nvars: %s alpha: %s)',
                              l.models[1], l.models[2], nvars, alpha)
      my.str       <- sprintf('%s model calculated with a target nvars: %s alpha: %s (base model: %s)',
                              l.models[1], nvars, alpha, l.models[2])
      my.str.base  <- sprintf('%s with target nvars: %s alpha: %s',
                              l.models[2], nvars, alpha)

      return(list(string = my.str, string.base = my.str.base, string.short = my.str.short,
                  model = l.models[1], base.model = l.models[2], target.vars = nvars,
                  alpha = alpha)[[ret.value]])
    } else {
      nvars <- splitted[2]
      alpha <- splitted[3]
alpha <- splitted[3]
      l.models <- gsub('glmnet', 'elastic net', splitted[1], ignore.case = TRUE)
      l.models <- gsub('degree_new', 'hub', l.models, ignore.case = TRUE)
      l.models <- proper(l.models)

      my.str.short <- sprintf('%s | nvars: %s alpha: %s',
                              l.models[1], l.models[2], nvars, alpha)
      my.str.base <- sprintf('%s with target nvars: %s alpha: %s',
                             l.models[1], nvars, alpha)

      return(list(string = NULL, string.base = my.str.base, string.short = my.str.short,
                  model = NULL, base.model = l.models[1], target.vars = nvars,
                  alpha = alpha)[[ret.value]])
    }
  })

  if (is.factor(label)) {
    levels(label) <- str.t
    return(label)
  }
  return(str.t)
}
