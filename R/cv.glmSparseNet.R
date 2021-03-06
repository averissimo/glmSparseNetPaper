#' Calculate cross validating GLM model with network-based regularization
#'
#' network parameter accepts:
#'
#'  * string to calculate network based on data (correlation, covariance)
#'  * matrix representing the network
#'  * vector with already calculated penalty weights (can also be used directly with glmnet)
#'
#' @param xdata input data, can be a matrix or MultiAssayExperiment
#' @param ydata response data compatible with glmnet
#' @param network type of network, see below
#' @param network.options options to calculate network
#' @param ... parameters that cv.glmnet accepts
#'
#' @return an object just as cv.glmnet
#' @export
#'
#' @examples
#' xdata <- matrix(rnorm(100), ncol = 20)
#' cv.glmSparseNet(xdata, rnorm(nrow(xdata)), 'correlation', family = 'gaussian')
#' cv.glmSparseNet(xdata, rnorm(nrow(xdata)), 'covariance', family = 'gaussian')
#'
#' #
#' #
#' # Using MultiAssayExperiment
#'
#' # load data
#' xdata <- MultiAssayExperiment::miniACC
#' # build valid data with days of last follow up or to event
#' event.ix <- which(!is.na(xdata$days_to_death))
#' cens.ix <- which(!is.na(xdata$days_to_last_followup))
#' xdata$surv_event_time <- array(NA, nrow(xdata@colData))
#' xdata$surv_event_time[event.ix] <- xdata$days_to_death[event.ix]
#' xdata$surv_event_time[cens.ix] <- xdata$days_to_last_followup[cens.ix]
#' # Keep only valid individuals
#' valid.ix <- as.vector(!is.na(xdata$surv_event_time) &
#'                       !is.na(xdata$vital_status) &
#'                       xdata$surv_event_time > 0)
#' xdata.valid <- xdata[, rownames(xdata@colData)[valid.ix]]
#' ydata.valid <- xdata.valid@colData[,c('surv_event_time', 'vital_status')]
#' colnames(ydata.valid) <- c('time', 'status')
#' cv.glmSparseNet(xdata.valid,
#'                   ydata.valid,
#'                   family          = 'cox',
#'                   network         = 'correlation',
#'                   experiment.name = 'RNASeq2GeneNorm',
#'                   mc.cores = 10)
setGeneric('cv.glmSparseNet.mclapply', function(xdata, ydata, network, network.options = network.options.default(), ...) {
  stop('wrong arguments, see help for cv.glmSparseNet and cv.glmnet')
})


#' Calculate GLM model with network-based regularization
#'
#' @inheritParams cv.glmSparseNet
#' @inherit cv.glmSparseNet return examples
#' @export
#'
setMethod('cv.glmSparseNet.mclapply', signature(xdata = 'matrix'), function(xdata, ydata, network,
                                                                   network.options = network.options.default(), ...) {
  return(glmSparseNet:::glmSparseNet.private(glmnet.mclapply::cv.glmnet, xdata, ydata, network, network.options = network.options, ...))
})

#' Calculate GLM model with network-based regularization
#'
#' @inheritParams cv.glmSparseNet
#' @param experiment.name Name of experiment in MultiAssayExperiment
#'
#' @inherit cv.glmSparseNet return examples
#' @export
setMethod('cv.glmSparseNet.mclapply', signature(xdata = 'MultiAssayExperiment'), function(xdata, ydata, network,
                                                                                 experiment.name = NULL,
                                                                                 network.options = network.options.default(), ...) {
  return(glmSparseNet:::glmSparseNet.private(glmnet.mclapply::cv.glmnet, xdata, ydata, network, experiment.name = experiment.name , network.options = network.options, ...))
})


#' Calculate GLM model with network-based regularization
#'
#' @inheritParams cv.glmSparseNet
#' @inherit cv.glmSparseNet return examples
#' @export
setMethod('cv.glmSparseNet.mclapply', signature(xdata = 'SummarizedExperiment'), function(xdata, ydata, network,
                                                                                 network.options = network.options.default(), ...) {
  return(glmSparseNet:::glmSparseNet.private(glmnet.mclapply::cv.glmnet, xdata, ydata, network, network.options = network.options, ...))
})



