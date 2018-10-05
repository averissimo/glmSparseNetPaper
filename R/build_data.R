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
prepare.tcga.survival.data <- function(project             = 'brca.data.2018.09.13',
                                       tissue.type         = 'primary.solid.tumor',
                                       assay               = 'RNASeqFPKM',
                                       handle.duplicates   = 'keep_first',
                                       coding.genes        = FALSE,
                                       log2.pre.normalize  = TRUE,
                                       normalization       = 'center',
                                       subtract.surv.column = NULL) {

  package.name <- project

  if (!require(package.name, character.only = TRUE)) {
    stop(sprintf('There is no package called \'%s\' installed, please go to https://github.com/averissimo/tcga.data/releases and install the corresponding release.'))
  }
  data("multiAssay", package = package.name)
  fpkm.dat <- build.matrix(assay, multiAssay)

  fpkm.per.tissue <- fpkm.dat$data

  # required to substract survival if that is the case
  for (ix.name in names(fpkm.per.tissue)) {
    colnames(fpkm.per.tissue[[ix.name]]) <- fpkm.dat$original.codes[[ix.name]]
  }

  futile.logger::flog.info('Loading data from %s package', package.name)
  futile.logger::flog.info('Types of tissue:\n * %s', paste(sprintf('%s (%d)', names(fpkm.per.tissue), sapply(fpkm.per.tissue, ncol)), collapse = '\n * '))

  xdata.raw <- fpkm.per.tissue[[tissue.type]]

  # remove genes that don't have any variability
  sd.xdata  <- sapply(seq(nrow(xdata.raw)), function(ix) { sd(xdata.raw[ix,]) })
  #
  futile.logger::flog.info('Non-expressed genes to be removed (from %d total genes) : %d', nrow(xdata.raw), sum(sd.xdata == 0))
  futile.logger::flog.info('  Remaining genes : %d', nrow(xdata.raw) - sum(sd.xdata == 0))
  xdata.raw <- xdata.raw[sd.xdata != 0,]

  #
  # Normalize
  if (log2.pre.normalize) {
    xdata.raw <- log2(1 + xdata.raw)
  }

  xdata.raw <- apply(xdata.raw, 1, function(row) {
    if (normalization == 'max') {
      if (max(abs(row)) == 0) {
        return(row)
      }
      return(row / max(abs(row)))
    } else if (normalization == 'center') {
      return(scale(row)[,1])
    } else {
      return(row)
    }
  })

  if (coding.genes) {
    futile.logger::flog.info('Using only coding genes:')
    coding    <- loose.rock::coding.genes()
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
  clinical <- fpkm.dat$clinical

  # load only patients with valid bcr_patient_barcode (non NA)
  ix.clinical <- !is.na(clinical[[tissue.type]]$bcr_patient_barcode)

  # build ydata data.frame
  if (is.null(subtract.surv.column)) {
    ydata <- data.frame(time   = clinical[[tissue.type]]$surv_event_time,
                        status = clinical[[tissue.type]]$vital_status)[ix.clinical,]
  } else {
    if (any(colnames(clinical[[tissue.type]]) == subtract.surv.column)) {
      ydata <- data.frame(time   = clinical[[tissue.type]]$surv_event_time - clinical[[tissue.type]][[subtract.surv.column]],
                          status = clinical[[tissue.type]]$vital_status)[ix.clinical,]
    } else {
      data("gdc.original", package = package.name)
      if (any(columns(gdc.original$bio.sample) == subtract.surv.column)) {
        bio.sample           <- gdc.original$bio.sample %>%
          dplyr::filter(bcr_sample_barcode %in% strtrim(rownames(xdata), 16))
        rownames(bio.sample) <- bio.sample$bcr_sample_barcode
        bio.sample <- bio.sample %>% select(bcr_patient_barcode, days_to_collection)
        bio.sample <- bio.sample[strtrim(rownames(xdata),16),]
        surv.days <- clinical[[tissue.type]]$surv_event_time - bio.sample$days_to_collection
        bio.sample[which(surv.days <= 0),]
        clinical[[tissue.type]]$days_to_collection <- bio.sample$days_to_collection
        clinical[[tissue.type]]$surv.days <- surv.days
        clinical[[tissue.type]][which(surv.days <= 0), c('bcr_patient_barcode', 'surv.days', 'days_to_collection', 'days_to_last_followup', 'days_to_death', 'surv_event_time', 'vital_status')] %>% as.data.frame %>% as_tibble()
      } else {
        stop('Could not find column to subtract the survival time.')
      }
    }
  }

  # name each row with patient code
  rownames(ydata) <- clinical[[tissue.type]]$bcr_patient_barcode[ix.clinical]

  # removing patients with:
  #  * negative follow-up
  #  * missing follow-up time
  futile.logger::flog.info('Number of patients removed with: \n * followup time <= 0:   %d\n * followup time is.na: %d',
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
    flog.warn('There are multiple samples for the same individual.. using strategy: \'%s\'', handle.duplicates)
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

  xdata.digest     <- loose.rock::digest.cache(xdata)
  xdata.raw.digest <- loose.rock::digest.cache(xdata.raw)

  ydata.digest <- loose.rock::digest.cache(ydata)

  return(list(xdata            = xdata,
              xdata.digest     = xdata.digest,
              ydata            = ydata,
              ydata.digest     = ydata.digest,
              xdata.raw        = xdata.raw,
              xdata.raw.digest = xdata.raw.digest,
              clinical         = clinical))

}
