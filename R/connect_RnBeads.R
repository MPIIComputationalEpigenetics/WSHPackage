#' create annotation
#'
#' This function creates a GRanges object from an input rnb.set object which then is input
#' for the IHS score calculation
#'
#' @param rnb.path path to an RnBSet object stored in disk
#'
#' @return GRanges object containing the CpG sites of interest
#' @noRd
create.annotation <- function(rnb.path){
  if(inherits(rnb.path,"RnBSet")){
    rnb.set <- rnb.path
  }else if(file.exists(rnb.path)){
    rnb.set <- load.rnb.set(rnb.path)
  }else{
    stop("Invalid value for rnb.set")
  }
  anno <- annotation(rnb.set)
  coverage <- covg(rnb.set)
  mean_coverage <- rowMeans(coverage,na.rm=TRUE)
  rm(coverage)
  high_enough <- mean_coverage >= get.option('coverage.threshold')
  anno <- anno[high_enough,]
  if(all(c("chrX","chrY")%in%anno$Chromosome)){
    anno <- anno[anno$Chromosome!="chrX",]
    anno <- anno[anno$Chromosome!="chrY",]
  }else if(all(c("X","Y")%in%anno$Chromosome)){
    anno <- anno[anno$Chromosome!="X",]
    anno <- anno[anno$Chromosome!="Y",]
  }
  anno <- GRanges(Rle(anno$Chromosome),IRanges(start=anno$Start,end=anno$End),anno$Strand)
  names <- paste(as.character(seqnames(anno)),start(ranges(anno)),end(ranges(anno)),sep="_")
  unique <- unique(names)
  match <- match(unique,names)
  anno <- anno[match]
  return(anno)
}

#' rnb.calculate.fdrp
#'
#' This function calculates the FDRP scores for a given RnBSet and a bam file containing the reads.
#'
#' @param rnb.set RnBSet object containing the CpG sites for which calculation should be conducted
#' @param bam.path path to the bam file containing the reads
#' @param log.path path to the log directory, if not existing it is created
#' @param cores cores available for the analysis
#'
#' @return FDRP scores for the CpG sites in the RnBSet with a higher coverage than COVERAGE.THRESHOLD
#' @export
rnb.calculate.fdrp <- function(rnb.set,bam.path,log.path=getwd(),cores=1){
  logger.start("Computing FDRP from RnBSet object")
  anno <- create.annotation(rnb.set)
  fdrps <- calculate.fdrp(bam.path,anno,log.path=log.path,cores=cores)
  logger.completed()
  return(fdrps)
}

#' rnb.calculate.qfdrp
#'
#' This function calculates the qFDRP scores for a given RnBSet and a bam file containing the reads.
#'
#' @param rnb.set RnBSet object containing the CpG sites for which calculation should be conducted
#' @param bam.path path to the bam file containing the reads
#' @param log.path path to the log directory, if not existing it is created
#' @param cores cores available for the analysis
#'
#' @return qFDRP scores for the CpG sites in the RnBSet with a higher coverage than COVERAGE.THRESHOLD
#' @export
rnb.calculate.qfdrp <- function(rnb.set,bam.path,log.path=getwd(),cores=1){
  logger.start("Computing qFDRP from RnBSet object")
  anno <- create.annotation(rnb.set)
  qfdrps <- calculate.qfdrp(bam.path,anno,log.path=log.path,cores=cores)
  logger.completed()
  return(qfdrps)
}

#' rnb.calculate.pdr
#'
#' This function calculates the PDR scores for a given RnBSet and a bam file containing the reads.
#'
#' @param rnb.set RnBSet object containing the CpG sites for which calculation should be conducted
#' @param bam.path path to the bam file containing the reads
#' @param log.path path to the log directory, if not existing it is created
#' @param cores cores available for the analysis
#'
#' @return PDR scores for the CpG sites in the RnBSet with a higher coverage than COVERAGE.THRESHOLD
#' @export
rnb.calculate.pdr <- function(rnb.set,bam.path,log.path=getwd(),cores=1){
  logger.start("Computing PDR from RnBSet object")
  anno <- create.annotation(rnb.set)
  pdrs <- calculate.pdrs(bam.path,anno,log.path=log.path,cores=cores)
  logger.completed()
  return(pdrs)
}

#' rnb.calculate.mhl
#'
#' This function calculates the MHL scores for a given RnBSet and a bam file containing the reads.
#'
#' @param rnb.set RnBSet object containing the CpG sites for which calculation should be conducted
#' @param bam.path path to the bam file containing the reads
#' @param roi.path path to a directory to which output can be added
#'
#' @return MHL scores for the CpG sites in the RnBSet with a higher coverage than COVERAGE.THRESHOLD
#' @export
rnb.calculate.mhl <- function(rnb.set,bam.path,roi.path=getwd()){
  logger.start("Computing MHL from RnBSet object")
  anno <- create.annotation(rnb.set)
  roi.df <- data.frame(Chromosome=seqnames(anno),Start=start(anno),End=end(anno))
  roi.location <- file.path(roi.path,"roi.bed")
  write.table(roi.df,roi.location,sep="\t",quote=F,row.names=F)
  mhl <- calculate.mhl(roi.location,bam.path)
  logger.completed()
  return(mhl)
}
