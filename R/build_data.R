#' Prepare FPKM data from TCGA project (specific tissue)
#'
#' This function will load data and pre-process it to be used in
#'  survival models.
#'
#' It will:
#'  * load data
#'  * handle duplicate samples for same individal (default is to keep only first)
#'  * remove individuals with missing vital_status or both follow-up/death time span
#'  * remove individuals with follow-up/death time span == 0
#'  * remove genes from RNASeqData (xdata) with standard deviation == 0
#'
#' @param project tcga project that has a package avaliable see https://github.com/averissimo/tcga.data
#' @param tissue.type type of tissue, can be 'primary.solid.tumor', 'metastatic', etc... depending on project.
#' @param input.type either 'rna' for RNASeq or 'dna' for mutation data
#' @param center scale xdata by subtracting by mean and dividing by standard deviation
#' @param include.negative.survival shift the survival times to include negative survival
#' @param handle.duplicates strategy to handle multiple samples for same individual, can take 'keep_first' or 'keep_all'
#' @param coding.genes filter the genes to only include coding genes, see loose.rock::coding.genes
#'
#' @return a list with data ready to be used in survival analysis, the 'xdata.raw' and 'ydata.raw' elements
#' have the full dataset for the specific tissue and the 'xdata' and 'ydata' have been cleaned by handling
#' patients with multiple samples, removing individuals with event time <= 0, missing and genes that have
#' standard_deviation == 0. It also returns a sha256 checksum for each of the data
#' @export
#'
#' @examples
#' # # Install a data package to load cancer data
#' # source("https://bioconductor.org/biocLite.R")
#' # biocLite(paste0('https://github.com/averissimo/tcga.data/releases/download/',
#' #                 '2016.12.15-brca/brca.data_1.0.tar.gz'))
#' # # Run
#' # prepare.tcga.survival.data('brca', 'primary.solid.tumor', 'keep_first')
#' # prepare.tcga.survival.data('brca', 'primary.solid.tumor', 'keep_first', input.type = 'dna')
prepare.tcga.survival.data <- function(project                      = 'brca',
                                       tissue.type                  = 'primary.solid.tumor',
                                       input.type                   = 'rna',
                                       center                       = TRUE,
                                       log2.normalize               = TRUE,
                                       include.negative.survival    = TRUE,
                                       handle.duplicates            = 'keep_first',
                                       coding.genes                 = FALSE) {

  package.name <- paste0(project, '.data')

  if (!package.name %in% utils::installed.packages()) {
    stop(sprintf('There is no package called \'%s\' installed, please go to https://github.com/averissimo/tcga.data/releases and install the corresponding release.',
                 package.name))
  }

  # An environment is necessary to adhere to best practices of ?data
  dat.env <- new.env()

  if (input.type == 'rna') {
    utils::data("fpkm.per.tissue", package = package.name, envir = dat.env)
    fpkm.per.tissue <- dat.env$fpkm.per.tissue

    futile.logger::flog.info('Loading data from %s package', package.name)
    futile.logger::flog.info('Types of tissue:\n * %s',
                             paste(sprintf('%s (%d)',
                                           names(fpkm.per.tissue),
                                           sapply(fpkm.per.tissue, ncol)), collapse = '\n * '))

    #
    # use only tissue from parameters
    xdata.raw <- t(fpkm.per.tissue[[tissue.type]])

    #
    # remove genes that don't have any variability
    sd.xdata  <- sapply(seq(ncol(xdata.raw)), function(ix) { stats::sd(xdata.raw[,ix]) })

    #
    futile.logger::flog.info('Non-expressed genes to be removed (from %d total genes) : %d', ncol(xdata.raw), sum(sd.xdata == 0))
    futile.logger::flog.info('  Remaining genes : %d', ncol(xdata.raw) - sum(sd.xdata == 0))
    xdata.raw <- xdata.raw[,sd.xdata != 0]

    #
    # Normalize

    if (log2.normalize == TRUE) {
      xdata.raw <- log2(xdata.raw + 1)
    }

    if (center == TRUE) {
      xdata.raw <- scale(xdata.raw, center = TRUE, scale = TRUE)
    } else {
      xdata.raw <- apply(xdata, 1, function(row) {
        if (max(abs(row)) == 0) {
          return(row)
        }
        return(row / max(abs(row)))
      })
    }
  } else if (input.type == 'dna') {
    utils::data("mutation", package = package.name, envir = dat.env)
    xdata.raw <- t(dat.env$mutation$count)

    # remove genes that don't have any variability
    sd.xdata  <- sapply(seq(ncol(xdata.raw)), function(ix) { stats::sd(xdata.raw[,ix]) })
    xdata.raw <- xdata.raw[,sd.xdata != 0]
  }



  if (coding.genes) {
    futile.logger::flog.info('Using only coding genes:')
    coding <- loose.rock::coding.genes()
    xdata.raw <- xdata.raw[,colnames(xdata.raw) %in% coding$ensembl_gene_id]
    futile.logger::flog.info('  * total coding genes: %d', length(coding$ensembl_gene_id))
    futile.logger::flog.info('  * coding genes in data: %d (new size of xdata)', ncol(xdata.raw))
  }

  if (handle.duplicates == 'keep_first') {
    xdata <- xdata.raw[!duplicated(strtrim(rownames(xdata.raw), 12)),]
  } else if (handle.duplicates == 'keep_all') {
    xdata <- xdata.raw
  }

  #
  # YDATA

  # load data
  utils::data('clinical', package = package.name, envir = dat.env)
  utils::data('gdc',      package = package.name, envir = dat.env)
  clinical  <- dat.env$clinical
  follow.up <- dat.env$gdc$follow.up

  # TODO: remove and use from package
  aa <- update.survival.from.followup(clinical$all, follow.up) %>% as.data.frame
  rownames(aa) <- aa$bcr_patient_barcode

  clinical[[tissue.type]][, 'vital_status']    <- aa[clinical[[tissue.type]]$bcr_patient_barcode, 'vital_status']
  clinical[[tissue.type]][, 'surv_event_time'] <- aa[clinical[[tissue.type]]$bcr_patient_barcode, 'surv_event_time']

  # load only patients with valid bcr_patient_barcode (non NA)
  ix.clinical <- !is.na(clinical[[tissue.type]]$bcr_patient_barcode)

  # build ydata data.frame
  ydata <- data.frame(time   = clinical[[tissue.type]]$surv_event_time,
                      status = clinical[[tissue.type]]$vital_status)[ix.clinical,]

  #
  # Shift time value if option to include all is present

  if (include.negative.survival && min(ydata$time, na.rm = TRUE) < 0) {
    ydata <- ydata %>%
      dplyr::mutate(time.tcga = time, time = time + abs(min(time, na.rm = TRUE)) + 1)
  }

  # name each row with patient code
  rownames(ydata) <- clinical[[tissue.type]]$bcr_patient_barcode[ix.clinical]


  # removing patients with:
  #  * negative follow-up
  #  * missing follow-up time
  futile.logger::flog.info('Number of patients removed with: \n * followup time < 0:   %d\n * followup time is.na: %d',
                           sum(!is.na(ydata$time) & ydata$time <= 0), sum(is.na(ydata$time)))
  ydata        <- ydata[!is.na(ydata$time) & ydata$time > 0,]

  # status description:
  #  * == 1 for dead (event happening)
  #  * == 0 for alive (censored)
  ydata$status <- ydata$status != 'Alive'

  # Multiple samples for same individual
  #  rename by appending to name
  xdata <- xdata[strtrim(rownames(xdata), 12) %in% rownames(ydata),]
  if (length(strtrim(rownames(xdata), 12)) != length(unique(strtrim(rownames(xdata), 12)))) {
    warning(sprintf('There are multiple samples for the same individual.. using strategy: \'%s\'', handle.duplicates))
    #
    new.row.names   <- strtrim(rownames(xdata), 12)
    rownames(xdata) <- sapply(seq_along(new.row.names), function(ix) {
      ix.name <- new.row.names[ix]
      count.b4 <- sum(new.row.names[1:ix] == ix.name)
      return(sprintf('%s.%d', ix.name, count.b4))
    })
    ydata <- ydata[strtrim(rownames(xdata), 12),]
    rownames(ydata) <- rownames(xdata)
  } else {
    rownames(xdata) <- strtrim(rownames(xdata), 12)
  }

  #
  # Pre-calculated sha256 to use in run.cache if necessary

  xdata.digest     <- loose.rock::digest.cache(xdata)
  xdata.raw.digest <- loose.rock::digest.cache(xdata.raw)

  ydata.digest     <- loose.rock::digest.cache(ydata)

  return(list(xdata            = xdata,
              xdata.raw        = xdata.raw,
              ydata            = ydata,
              clinical         = clinical,
              xdata.digest     = xdata.digest,
              xdata.raw.digest = xdata.raw.digest,
              ydata.digest     = ydata.digest))

}
