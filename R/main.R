#' main.R
#'
#' This scripts contains main functions to calculate the scores.
#'
#'
## G L O B A L S #######################################################################################################

VALID.SCORES <- c("fdrp","qfdrp","pdr","epipolymorphism","entropy","mhl")

IHS.OPTIONS <- new.env()
assign('ALL',c('window.size','mapq.filter','max.reads','min.overlap','fdrp.type','coverage.threshold','methclone.methylation.diff','perl.path','samtools.path'),IHS.OPTIONS)
assign('WINDOW.SIZE',50,IHS.OPTIONS)
assign('MAPQ.FILTER',35,IHS.OPTIONS)
assign('MAX.READS',40,IHS.OPTIONS)
assign('MIN.OVERLAP',35,IHS.OPTIONS)
assign('FDRP.TYPE','FDRP',IHS.OPTIONS)
assign('COVERAGE.THRESHOLD',10,IHS.OPTIONS)
assign('METHCLONE.METHYLATION.DIFF',0,IHS.OPTIONS)
assign('PERL.PATH',"/usr/bin/perl",IHS.OPTIONS)
assign('SAMTOOLS.PATH',"/usr/bin",IHS.OPTIONS)

## F U N C T I O N S ###################################################################################################

#' set.option
#'
#' Change global options for IHS score calculation
#'
#' @param window.size Window around the CpG site of interest to consider in FDRP and qFDRP calculation, the higher the value, the more
#' likely it is to find heterogeneity
#' @param mapq.filter mapq filter used to filter out low quality reads
#' @param max.reads Maximum number of reads to be considered in FDRP and qFDRP calculation. The scores compute all pairs, therefore this
#' is a crucial parameter for the running time of the calculation.
#' @param min.overlap Miniumum overlap between two reads to consider it as a read pair in FDRP/qFDRP calculation in bp.
#' @param fdrp.type FDRP type to be used: either FDRP or qFDRP
#' @param perl.path  Path to an installation of perl on the machine
#' @param coverage.threshold Coverage Threshold emloyed to select the sites in a RnBSet annotation that fullfill havin a coverage
#'                   higher than this threshold
#' @param methclone.methylation.diff Methylation difference parameter employed by the methclone software. Only sites are considered
#'                   that have a methylation difference higher than this value in the methclone package.
#' @param samtools.path path to the directory where samtools is located in your machine
#'
#' @export
#'
set.option <- function(window.size,
                       mapq.filter,
                       max.reads,
                       min.overlap,
                       fdrp.type,
                       coverage.threshold,
                       methclone.methylation.diff,
                       perl.path,
                       samtools.path){
  if(!missing(window.size)) IHS.OPTIONS[['WINDOW.SIZE']] <- window.size
  if(!missing(mapq.filter)) IHS.OPTIONS[['MAPQ.FILTER']] <- mapq.filter
  if(!missing(max.reads)) IHS.OPTIONS[['MAX.READS']] <- max.reads
  if(!missing(min.overlap)) IHS.OPTIONS[['MIN.OVERLAP']] <- min.overlap
  if(!missing(coverage.threshold)) IHS.OPTIONS[['COVERAGE.THRESHOLD']] <- coverage.threshold
  if(!missing(methclone.methylation.diff)) IHS.OPTIONS[['METHCLONE.METHYLATION.DIFF']] <- methclone.methylation.diff
  if(!missing(perl.path)) IHS.OPTIONS[['PERL.PATH']] <- perl.path
  if(!missing(samtools.path)){
    example.bam <- system.file(file.path("extData","small_example.bam"),package="ISH")
    cmd <- paste0(samtools.path,"/samtools view -H ",example.bam)
    tryCatch(output <- system(command = cmd,intern = T),error=function(e){
      stop("Invalid value for samtools.path, please specify a directory with a working version of samtools")
    })
    if(substr(output[1],1,1)!="@"){
      stop("Invalid value for samtools.path, please specify a directory with a working version of samtools")
    }
    IHS.OPTIONS[['SAMTOOLS.PATH']] <- samtools.path
  }
  if(!missing(fdrp.type)){
    if(!(fdrp.type%in%c('FDRP','qFDRP'))) stop('Invalid Value for fdrp.type, either FDRP or qFDRP requred')
    IHS.OPTIONS[['FDRP.TYPE']] <- fdrp.type
  }
}

#' get.option
#' Print the value of the global option
#'
#' @param names string or character vector containing the names of the options to be printed
#'
#' @return the option for the specified option
#' @export
get.option <- function(names){
  if(!all(names %in% IHS.OPTIONS[['ALL']])){
    stop(paste0('No option(s) available named: ',names[!(names%in%IHS.OPTIONS[['ALL']])]))
  }
  ret <- c()
  if('window.size'%in%names){
    ret <- c(ret,window.size=IHS.OPTIONS[['WINDOW.SIZE']])
  }
  if('mapq.filter'%in%names){
    ret <- c(ret,mapq.filter=IHS.OPTIONS[['MAPQ.FILTER']])
  }
  if('max.reads'%in%names){
    ret <- c(ret,max.reads=IHS.OPTIONS[['MAX.READS']])
  }
  if('min.overlap'%in%names){
    ret <- c(ret,min.overlap=IHS.OPTIONS[['MIN.OVERLAP']])
  }
  if('coverage.threshold'%in%names){
    ret <- c(ret,coverage.threshold=IHS.OPTIONS[['COVERAGE.THRESHOLD']])
  }
  if('methclone.methylation.diff'%in%names){
    ret <- c(ret,methclone.methylation.diff=IHS.OPTIONS[['METHCLONE.METHYLATION.DIFF']])
  }
  if('fdrp.type'%in%names){
    ret <- c(ret,fdrp.type=IHS.OPTIONS[['FDRP.TYPE']])
  }
  if('perl.path'%in%names){
    ret <- c(ret,perl.path=IHS.OPTIONS[['PERL.PATH']])
  }
  if('samtools.path'%in%names){
    ret <- c(ret,samtools.path=IHS.OPTIONS[['SAMTOOLS.PATH']])
  }
  return(ret[names])
}

#' ish.run.example
#'
#' Returns the scores for a small example data set of size 50kb
#'
#' @param score The ISH score which should be computed, needs to be one of \code{fdrp},\code{qfdrp},\code{pdr},\code{epipolymorphism},
#'               \code{entropy} or \code{mhl}
#'
#' @return A data frame containing the annotation and the corresponding score
#'
#' @author Michael Scherer
#' @export
ish.run.example <- function(score="qfdrp"){
  if(!(score%in%VALID.SCORES)){
    stop(paste("Invalid value for score, must be one of",VALID.SCORES))
  }
  logger.start("ISH score example")
  example.rnb.set <- system.file(file.path("extData","small_rnbSet.zip"),package = "ISH")
  example.bam <- system.file(file.path("extData","small_example.bam"),package = "ISH")
  score <- compute.score.rnb(example.bam,example.rnb.set,score)
  logger.completed()
  return(score)
}

#' compute.score.rnb
#'
#' Main function to compute ISH scores based on a input RnBSet and a bam file containing the reads.
#'
#' @param bam.file path to bam file containing the reads
#' @param rnb.set path to RnBSet contaning methylation, coverage and sample meta information
#' @param score The ISH score which should be computed, needs to be one of \code{fdrp},\code{qfdrp},\code{pdr},\code{epipolymorphism},
#'               \code{entropy} or \code{mhl}
#'
#' @return data frame containing the annotation and the computed ISH scores
#' @author Michael Scherer
#' @export
compute.score.rnb <- function(bam.file,rnb.set,score){
  if(!(score%in%VALID.SCORES)){
    stop(paste("Invalid value for score, must be one of",VALID.SCORES))
  }
  if(!(file.exists(bam.file))){
    stop(paste("File",bam.file,"does not exist"))
  }
  if(!(file.exists(rnb.set))){
    stop(paste("File",rnb.set,"does not exist"))
  }
  if(score=="qfdrp"){
    ret <- rnb.calculate.qfdrp(rnb.set,bam.file)
  }
  if(score=="fdrp"){
    ret <- rnb.calculate.fdrp(rnb.set,bam.file)
  }
  if(score=="pdr"){
    ret <- rnb.calculate.pdr(rnb.set,bam.file)
  }
  if(score=="mhl"){
    ret <- rnb.calculate.mhl(rnb.set,bam.file)
  }
  if(score%in%c("epipolymorphism","entropy")){
    logger.warning("Epipolymorphism and Entropy do not take a RnBSet object as input. The annotation is computed from the bam file directly.")
    if(score=="epipolymorphism"){
      ret <- calculate.epipolymorphism(bam.file)
    }else{
      ret <- calculate.entropy(bam.file)
    }
  }
  return(ret)
}

#' compute.score.GRanges
#'
#' Main function to compute ISH scores based on a input GRanges object and a bam file containing the reads.
#'
#' @param bam.file path to bam file containing the reads
#' @param range GRanges object containing the annotation
#' @param score The ISH score which should be computed, needs to be one of \code{fdrp},\code{qfdrp},\code{pdr},\code{epipolymorphism},
#'               \code{entropy} or \code{mhl}
#'
#' @return data frame containing the annotation and the computed ISH scores
#' @author Michael Scherer
#' @export
compute.score.GRanges <- function(bam.file,range,score){
  if(!(score%in%VALID.SCORES)){
    stop(paste("Invalid value for score, must be one of",VALID.SCORES))
  }
  if(!(file.exists(bam.file))){
    stop(paste("File",bam.file,"does not exist"))
  }
  if(!(file.exists(rnb.set))){
    stop(paste("File",rnb.set,"does not exist"))
  }
  if(!(inherits(range,"GRanges"))){
    stop(paste("Invalid value for range, needs to be a GRanges object"))
  }
  if(score=="qfdrp"){
    ret <- calculate.qfdrp(rnb.set,anno)
  }
  if(score=="fdrp"){
    ret <- calculate.fdrp(rnb.set,anno)
  }
  if(score=="pdr"){
    ret <- rnb.calculate.pdr(rnb.set,bam.file)
  }
  if(score=="mhl"){
    # do transformation of GRanges to bed file
    ret <- calculate.mhl(rnb.set,bam.file)
  }
  if(score%in%c("epipolymorphism","entropy")){
    logger.warning("Epipolymorphism and Entropy do not take a RnBSet object as input. The annotation is computed from the bam file directly.")
    if(score=="epipolymorphism"){
      ret <- calculate.epipolymorphism(bam.file)
    }else{
      ret <- calculate.entropy(bam.file)
    }
  }
  return(ret)
}