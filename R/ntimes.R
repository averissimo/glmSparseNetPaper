# train.perc <- params$train
# subset <- params$subset
# seed <- 1985
call.results <- function(seed,
                         xdata,
                         ydata,
                         train.perc,
                         subset,
                         penalty.factor.degree.new,
                         penalty.factor.degree.old,
                         penalty.factor.orphan,
                         params) {
  #flog.info('seed: %d', seed)
  set.seed(seed)

  #
  # Build training data
  ixs <- glmSparseNet::balanced.train.and.test(which(ydata$status), which(!ydata$status), train.perc = train.perc)
  #flog.info('ixs:', ixs, capture = T)
  xdata.test <- xdata[ixs$test,]
  ydata.test <- ydata[ixs$test,]
  #
  xdata.train <- xdata[ixs$train,]
  ydata.train <- ydata[ixs$train,]

  xdata.ix <- seq(ncol(xdata))
  xdata.ix.no.added <- xdata.ix

  if (subset < ncol(xdata.train)) {
    xdata.ix <- sample(xdata.ix, params$subset)
  }

  xdata.train.digest <- glmSparseNet::digest.cache(xdata.train[, xdata.ix])

  #
  # MODELS
  #

  models  <- list()
  lambdas <- list()
  coefs   <- list()
  result  <- list()

  glmnet.params <- list()

  for (target.name in names(params$target.vars)) {
    target <- params$target.vars[[target.name]]
    glmnet.params <- c(glmnet.params, list(list(penalty = rep(1, ncol(xdata.train)),
                                                name = 'glmnet',
                                                target = target,
                                                target.name = target.name)))
    glmnet.params <- c(glmnet.params, list(list(penalty = penalty.factor.degree.new,
                                                name = 'degree_new',
                                                target = target,
                                                target.name = target.name)))
    glmnet.params <- c(glmnet.params, list(list(penalty = penalty.factor.orphan,
                                                name = 'orphan',
                                                target = target,
                                                target.name = target.name)))
    if (params$calc.params.old) {
      glmnet.params <- c(glmnet.params, list(list(penalty = penalty.factor.degree.old,
                                                  name = 'degree_old',
                                                  target = target,
                                                  target.name = target.name)))
      glmnet.params <- c(glmnet.params, list(list(penalty = penalty.factor.degree.log,
                                                  name = 'degree_log',
                                                  target = target,
                                                  target.name = target.name)))
    }
  }

  outer.result <- mclapply(seq_along(glmnet.params), function(ix) {
    el       <- glmnet.params[[ix]]
    ix.name  <- sprintf('%s.%s.%d', el$name, el$target.name, el$target)
    ix.cache <- sprintf('%s_models', el$name)
    #
    suppressWarnings(
      result  <- glmSparseNet.cox(xdata        = xdata.train[,xdata.ix],
                                  ydata        = ydata.train,
                                  target.vars  = el$target,
                                  alpha        = params$alpha,
                                  network      = el$penalty,
                                  xdata.digest = xdata.train.digest,
                                  cache.prefix = ix.cache)
    )

    #
    return(list(result = result, name = ix.name))
    #})

  }, mc.cores = 1, mc.allow.recursive = FALSE)


  for (ix in outer.result) {
    result[[ix$name]]  <- ix$result
    models[[ix$name]]  <- ix$result$model
    lambdas[[ix$name]] <- ix$result$lambda
    coefs[[ix$name]]   <- ix$result$coef
  }

  #
  # Kaplan-Meier (p.value) and C-INDEX
  #

  #
  km.train <- list()
  km.test  <- list()
  #
  c.index.train <- list()
  c.index.test  <- list()

  for (ix.name in names(coefs)) {
    km.train[[ix.name]] <- my.draw.kaplan(list(ix.name = coefs[[ix.name]]), plot.title = 'Train set', xdata.train[, xdata.ix], ydata.train)$pvalue
    km.test[[ix.name]]  <- my.draw.kaplan(list(ix.name = coefs[[ix.name]]), plot.title = 'Test set', xdata.test[, xdata.ix], ydata.test)$pvalue
    #
    c.index.train[[ix.name]] <- fit.risk(coefs[[ix.name]], xdata.train[, xdata.ix], ydata.train, ix.name, n.cores = 1, show.message = FALSE)
    c.index.test[[ix.name]]  <- fit.risk(coefs[[ix.name]], xdata.test[, xdata.ix], ydata.test, ix.name, n.cores = 1, show.message = FALSE)
  }

  return(list(metrics = list(km.train      = km.train,
                             km.test       = km.test,
                             c.index.train = c.index.train,
                             c.index.test  = c.index.test),
              coefs = coefs,
              ixs   = ixs))
}
