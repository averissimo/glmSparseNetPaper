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
#' \dontrun{
#'     # Install a data package to load cancer data
#'     source("https://bioconductor.org/biocLite.R")
#'     biocLite(paste0('https://github.com/averissimo/tcga.data/releases/download/',
#'                     '2016.12.15-brca/brca.data_1.0.tar.gz'))
#'     prepare.tcga.survival.data('brca', 'primary.solid.tumor', 'keep_first')
#'     prepare.tcga.survival.data('brca', 'primary.solid.tumor', 'keep_first', input.type = 'dna')
#' }
prepare.tcga.survival.data <- function(project                      = 'brca',
                                       tissue.type                  = 'primary.solid.tumor',
                                       input.type                   = 'rna',
                                       normalization                = 'none',
                                       log2.pre.normalize           = FALSE,
                                       include.negative.survival    = FALSE,
                                       handle.duplicates            = 'keep_first',
                                       only.coding.genes            = FALSE) {

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

    if (log2.pre.normalize == TRUE) {
      xdata.raw <- log2(xdata.raw + 1)
    }

    if (is.null(normalization)) {
      normalization = 'none'
    }

    if (normalization == 'center') {
      # scales the columns of the matrix
      xdata.raw2 <- scale(xdata.raw, center = TRUE, scale = TRUE)
    } else if (normalization == 'max') {
      # scales the columns by maximum value
      xdata.raw2 <- apply(xdata.raw, 2, function(feature) {
        if (max(abs(feature)) == 0) {
          return(feature)
        }
        return(feature / max(abs(feature)))
      })
    } else {
      xdata.raw2 <- xdata.raw
    }
  } else if (input.type == 'dna') {
    utils::data("mutation", package = package.name, envir = dat.env)
    xdata.raw <- t(dat.env$mutation$count)

    # remove genes that don't have any variability
    sd.xdata  <- sapply(seq(ncol(xdata.raw)), function(ix) { stats::sd(xdata.raw[,ix]) })
    xdata.raw <- xdata.raw[,sd.xdata != 0]
    xdata.raw2 <- xdata.raw
  }



  if (only.coding.genes) {
    futile.logger::flog.info('Using only coding genes:')
    coding <- glmSparseNet::run.cache(loose.rock::coding.genes)
    xdata.raw2 <- xdata.raw2[,colnames(xdata.raw2) %in% coding$ensembl_gene_id]
    futile.logger::flog.info('  * total coding genes: %d', length(coding$ensembl_gene_id))
    futile.logger::flog.info('  * coding genes in data: %d (new size of xdata)', ncol(xdata.raw2))
  }

  if (handle.duplicates == 'keep_first') {
    xdata <- xdata.raw2[!duplicated(strtrim(rownames(xdata.raw2), 12)),]
  } else if (handle.duplicates == 'keep_all') {
    xdata <- xdata.raw2
  }

  #
  # YDATA

  curated <- glmSparseNet::run.cache(curatedTCGAData::curatedTCGAData,
                                   project, 'RNASeq2GeneNorm', FALSE)

  ydata <- curated@colData %>% as.data.frame %>%
    dplyr::select(patientID, status = vital_status, days_to_last_followup, days_to_death) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(time = max(-Inf, days_to_last_followup, days_to_death, na.rm = TRUE),
                  status = status == 1) %>%
    dplyr::select(patientID, time, status) %>%
    as.data.frame

  rownames(ydata) <- ydata$patientID
  ydata <- ydata %>% dplyr::select(time, status)

  # removing patients with:
  #  * negative follow-up
  #  * missing follow-up time
  futile.logger::flog.info('Number of patients removed with: \n * followup time < 0:   %d\n * followup time is.na: %d',
                           sum(!is.na(ydata$time) & ydata$time <= 0), sum(is.na(ydata$time)))
  ydata        <- ydata[!is.na(ydata$time) & ydata$time > 0,]

  # status description:
  #  * == 1 for dead (event happening)
  #  * == 0 for alive (censored)
  # ydata$status <- ydata$status != 'Alive'

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

  ydata <- ydata[rownames(ydata) %in% strtrim(rownames(xdata), 12),]
  ydata <- ydata[strtrim(rownames(xdata), 12), ]

  #
  # Pre-calculated sha256 to use in run.cache if necessary

  xdata.digest     <- glmSparseNet::digest.cache(xdata)
  xdata.raw.digest <- glmSparseNet::digest.cache(xdata.raw)

  ydata.digest     <- glmSparseNet::digest.cache(ydata)

  return(list(xdata            = xdata,
              xdata.raw        = xdata.raw,
              ydata            = ydata,
              clinical         = colData(curated[rownames(ydata),]),
              xdata.digest     = xdata.digest,
              xdata.raw.digest = xdata.raw.digest,
              ydata.digest     = ydata.digest))

}
