build.multiple.runs <-
  function(params, xdata, ydata,
           penalty.factor.degree.new, penalty.factor.degree.old,
           penalty.factor.orphan,
           xdata.digest.cache,
           ydata.digest.cache) {

  set.seed(params$seed)
  seed.vec <- sample(1000 + (max(c(1:1e7,params$ntimes))))[1:params$ntimes]

  varnames <- colnames(xdata)

  if (!exists('global.n.cores')) {
    global.n.cores <- 1
  }

  ntimes.results <- parallel::mclapply(seed.vec, function(seed) {
    calc.params  <-
      params[c('alpha', 'subset', 'train', 'calc.params.old', 'target.vars',
               'balanced.sets')]
    calc.params$seed <- seed
    #
    flog.info(sprintf('Resampling with seed: %d (%d)',
                      calc.params$seed, which(calc.params$seed == seed.vec)))
    my.output <- capture.output(
      ret.val <- run.cache(call.results,
                           xdata, ydata,
                           penalty.factor.degree.new,
                           penalty.factor.degree.old,
                           penalty.factor.orphan,
                           calc.params,
                           no.models = TRUE,
                           no.plot = TRUE,
                           #
                           cache.digest = list(xdata.digest.cache,
                                               ydata.digest.cache),
                           cache.prefix = 'big-diff',
                           show.message = TRUE)
    )

    for (ix.name in names(ret.val$coefs)) {
      ret.val$coefs[[ix.name]] <- as(ret.val$coefs[[ix.name]], 'sparseVector')
    }

    flog.info(sprintf('Finish with seed: %d (%d)',
                      seed, which(seed == seed.vec)), my.output, capture = TRUE)
    flog.info('------------------------ end of iteration')

    # Reduce memory footprint
    ret.val$result <- NULL
    ret.val$models <- NULL
    ret.val$varnames <- NULL # this is returned in upper
    for (ix.coef.names in names(ret.val$coefs)) {
      names(ret.val$coefs[[ix.coef.names]]) <- NULL
    }
    return(ret.val)
  }, mc.cores = min(global.n.cores), mc.allow.recursive = FALSE, mc.preschedule = FALSE)

  return(list(ntimes.results = ntimes.results, varnames = varnames))
}


handle.multiple <- function(ntimes.results) {
  big.df    <- data.frame()
  rank.df   <- tibble()
  values.df <- tibble()

  for (ix in seq_along(ntimes.results)) {
    if (!is.null(ntimes.results[[ix]]) && is.list(ntimes.results[[ix]])) {
      el <- ntimes.results[[ix]]$metrics
      for (ix.el in names(el)) {
        my.values <- sapply(names(el[[ix.el]]), function(ix.model) {el[[ix.el]][[ix.model]]})
        my.names  <- rep(ix.el, length(my.values))
        my.models <- names(el[[ix.el]])
        new.line  <- data.frame(ix     = ix,
                                metric = my.names,
                                model  = my.models,
                                values = as.numeric(my.values),
                                stringsAsFactors = FALSE)
        big.df    <- rbind(big.df, new.line)

        #
        sorted.ix <- sort(my.values)
        base.models <- prettify.labels(my.models, 'string.base')

        for(base.ix in unique(base.models)) {
          ixs <- which(base.models == base.ix)
          my.rank <- my.models[ixs][sort(my.values[ixs], index.return = TRUE)$ix]
          order.fun <- if (grepl('c.index', my.names[1])) {
            function(vec) { rev(vec) }
          } else {
            function(vec) { vec }
          }
          my.rank <- order.fun(my.rank)
          new.el <- c(list(my.names[1], base.ix), prettify.labels(my.rank, 'model'))
          names(new.el) <- c('type', 'base.model', paste0('X', seq_len(length(new.el) - 2)))
          rank.df <- rbind(rank.df, as.tibble(new.el))

          ix.val <- sort(prettify.labels(my.models[ixs], 'model'), index.return = TRUE)$ix
          new.el.val <- c(list(my.names[1], base.ix), my.values[ixs][ix.val])
          names(new.el.val) <- c('type', 'base.model', prettify.labels(my.models[ixs][ix.val], 'model'))
          values.df <- rbind(values.df, as.tibble(new.el.val))
        }
      }
    }
  }
  levels.str <- c('c.index.train', 'c.index.test', 'km.train', 'km.test')
  labels.str <- c('C-Index (Train set)', 'C-Index (Test set)',
                  'Log-rank (Train set)', 'Log-rank (Test set)')

  big.df$metric  <- factor(big.df$metric, levels = levels.str, labels = labels.str)
  rank.df$type   <- factor(rank.df$type, levels = levels.str, labels = labels.str)
  values.df$type <- factor(values.df$type, levels = levels.str, labels = labels.str)

  return(list(big.df         = big.df,
              rank.df        = rank.df,
              values.df      = values.df))
}
