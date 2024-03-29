#' main.R
#'
#' This scripts contains main functions to calculate the scores.
#'
#'
## G L O B A L S #######################################################################################################

VALID.SCORES <- c("fdrp","qfdrp","pdr","epipolymorphism","entropy","mhl")
WINDOWS.CAPABLE <- c("fdrp","qfdrp","pdr")

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
#' Change global options for WSH score calculation
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
#' @author Michael Scherer
#' @examples
#' \donttest{
#' get.option("coverage.threshold")
#' set.option(coverage.threshold=42)
#' get.option("coverage.threshold")
#' }
set.option <- function(window.size=50,
                       mapq.filter=35,
                       max.reads=40,
                       min.overlap=35,
                       fdrp.type='FDRP',
                       coverage.threshold=10,
                       methclone.methylation.diff=0,
                       perl.path="/usr/bin/perl",
                       samtools.path="/usr/bin"){
  if(length(window.size)!=1){
    stop("Please specify the options one by one, not as a vector or list.")
  }
  if(!missing(window.size)) IHS.OPTIONS[['WINDOW.SIZE']] <- window.size
  if(!missing(mapq.filter)) IHS.OPTIONS[['MAPQ.FILTER']] <- mapq.filter
  if(!missing(max.reads)) IHS.OPTIONS[['MAX.READS']] <- max.reads
  if(!missing(min.overlap)) IHS.OPTIONS[['MIN.OVERLAP']] <- min.overlap
  if(!missing(coverage.threshold)) IHS.OPTIONS[['COVERAGE.THRESHOLD']] <- coverage.threshold
  if(!missing(methclone.methylation.diff)) IHS.OPTIONS[['METHCLONE.METHYLATION.DIFF']] <- methclone.methylation.diff
  if(!missing(perl.path)){
    cmd <- paste(perl.path,"--help")
    tryCatch(output <- system(command = cmd,intern = T),error=function(e){
      stop("Invalid value for perl.path, please specify the location of a working perl version")
    })
    if(is.null(output)){
      stop("Invalid value for perl.path, please specify the location of a working perl version")
    }
    IHS.OPTIONS[['PERL.PATH']] <- perl.path
  }
  if(!missing(samtools.path)){
    example.bam <- system.file(file.path("extData","small_example.bam"),package="WSH")
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
#' @author Michael Scherer
#' @examples
#' \donttest{
#' get.option("coverage.threshold")
#' set.option(coverage.threshold=42)
#' get.option("coverage.threshold")
#' }
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

#' wsh.run.example
#'
#' Returns the scores for a small example data set of size 50kb
#'
#' @param score The WSH score which should be computed, needs to be one of \code{fdrp},\code{qfdrp},\code{pdr},\code{epipolymorphism},
#'               \code{entropy} or \code{mhl}
#'
#' @return A data frame containing the annotation and the corresponding score
#'
#' @author Michael Scherer
#' @examples
#' \donttest{
#' wsh.run.example("pdr")
#' }
#' @export
wsh.run.example <- function(score="qfdrp"){
  logger.start("WSH score example")
  example.GRanges <- system.file(file.path("extData","example_GRanges.RData"),package = "WSH")
  load(example.GRanges)
  example.bam <- system.file(file.path("extData","small_example.bam"),package = "WSH")
  score <- compute.score.GRanges(example.bam,example.GRanges,score)
  logger.completed()
  return(score)
}

#' compute.score.rnb
#'
#' Main function to compute WSH scores based on a input RnBSet and a bam file containing the reads.
#'
#' @param bam.file path to bam file containing the reads
#' @param rnb.set path to RnBSet contaning methylation, coverage and sample meta information
#' @param score The WSH score which should be computed, needs to be one of \code{fdrp},\code{qfdrp},\code{pdr},\code{epipolymorphism},
#'               \code{entropy} or \code{mhl}
#' @param use.sex.chromosomes Flag indicating if scores are also to be computed for the sex chromosomes
#' @param ... Further parameters passed to \code{calculate.*}
#'
#' @return data frame containing the annotation and the computed WSH scores
#' @author Michael Scherer
#' @examples
#' \dontrun{
#' example.rnb.set <- system.file(file.path("extData","small_rnbSet.zip"),package="WSH")
#' example.bam <- system.file(file.path("extData","small_example.bam"),package="WSH")
#' fdrp <- compute.score.rnb(bam.file=example.bam,rnb.set=example.rnb.set)
#' }
#' @export
compute.score.rnb <- function(bam.file,
                              rnb.set,
                              score="qfdrp",
                              use.sex.chromosomes=use.sex.chromosomes,
                              ...){
  wsh.check.validity(score)
  if(!(file.exists(bam.file))){
    stop(paste("File",bam.file,"does not exist"))
  }
  if(!(inherits(rnb.set,"RnBSet"))){
    if(!(file.exists(rnb.set))){
      stop(paste("File",rnb.set,"does not exist"))
    }
  }
  if(score=="qfdrp"){
    ret <- rnb.calculate.qfdrp(rnb.set,bam.file,use.sex.chromosomes=use.sex.chromosomes,...)
  }
  if(score=="fdrp"){
    ret <- rnb.calculate.fdrp(rnb.set,bam.file,use.sex.chromosomes=use.sex.chromosomes,...)
  }
  if(score=="pdr"){
    ret <- rnb.calculate.pdr(rnb.set,bam.file,use.sex.chromosomes=use.sex.chromosomes,...)
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
#' Main function to compute WSH scores based on a input GRanges object and a bam file containing the reads.
#'
#' @param bam.file path to bam file containing the reads
#' @param range GRanges object containing the annotation
#' @param score The WSH score which should be computed, needs to be one of \code{fdrp},\code{qfdrp},\code{pdr},\code{epipolymorphism},
#'               \code{entropy} or \code{mhl}
#' @param use.sex.chromosomes Flag indicating if scores are also to be computed for the sex chromosomes
#' @param ... Further parameters passed to \code{calculate.*}
#'
#' @return data frame containing the annotation and the computed WSH scores
#' @author Michael Scherer
#' @examples
#' \donttest{
#' load(system.file(file.path("extData","example_GRanges.RData"),package="WSH"))
#' example.bam <- system.file(file.path("extData","small_example.bam"),package="WSH")
#' fdrp <- compute.score.GRanges(bam.file=example.bam,range=example.GRanges)
#' }
#' @export
compute.score.GRanges <- function(bam.file,
                                  range,
                                  score="qfdrp",
                                  use.sex.chromosomes=FALSE,
                                  ...){
  wsh.check.validity(score)
  if(!(file.exists(bam.file))){
    stop(paste("File",bam.file,"does not exist"))
  }
  if(!(inherits(range,"GRanges"))){
    stop(paste("Invalid value for range, needs to be a GRanges object"))
  }
  if(score=="qfdrp"){
    ret <- calculate.qfdrp(bam.file,range,use.sex.chromosomes=use.sex.chromosomes,...)
  }
  if(score=="fdrp"){
    ret <- calculate.fdrp(bam.file,range,use.sex.chromosomes=use.sex.chromosomes,...)
  }
  if(score=="pdr"){
    ret <- calculate.pdr(bam.file,range,use.sex.chromosomes=use.sex.chromosomes,...)
  }
  if(score=="mhl"){
    # do transformation of GRanges to bed file
    range.table <- data.frame(Chromosome=seqnames(range),Start=start(range),End=end(range))
    mhl.file <- paste0(tempfile(pattern="WSH_"),".bed")
    write.table(range.table,mhl.file,sep="\t",row.names=FALSE,quote=FALSE)
    ret <- calculate.mhl(mhl.file,bam.file)
  }
  if(score%in%c("epipolymorphism","entropy")){
    logger.warning("Epipolymorphism and Entropy do not take a GRanges object as input. The annotation is computed from the bam file directly.")
    if(score=="epipolymorphism"){
      ret <- calculate.epipolymorphism(bam.file)
    }else{
      ret <- calculate.entropy(bam.file)
    }
  }
  return(ret)
}

#' compute.score
#'
#' Generic function to call, passes its arugments either to \code{\link{compute.score.rnb}} or
#' \code{\link{compute.score.GRanges}}.
#'
#' @param bam.file path to bam file containing the reads
#' @param ... additional arugment. Either RnBSet, GRanges or empty (only for Epipolymorphism and Entropy)
#' @param score The WSH score which should be computed, needs to be one of \code{fdrp},\code{qfdrp},\code{pdr},\code{epipolymorphism},
#'               \code{entropy} or \code{mhl}
#' @param use.sex.chromosomes Flag indicating if scores are also to be computed for the sex chromosomes
#' @param ignore.strand  The \code{ignore.strand} parameter from the function \code{\link{findOverlaps}}
#'
#' @return data frame containing the annotation and the computed WSH scores
#' @author Michael Scherer
#' @examples
#' \donttest{
#' load(system.file(file.path("extData","example_GRanges.RData"),package="WSH"))
#' example.rnb.set <- system.file(file.path("extData","small_rnbSet.zip"),package="WSH")
#' example.bam <- system.file(file.path("extData","small_example.bam"),package="WSH")
#' fdrp <- compute.score(bam.file=example.bam,example.GRanges,score="fdrp")
#' qfdrp <- compute.score(bam.file=example.bam,example.rnb.set,score="qfdrp")
#' }
#' @export
compute.score <- function(bam.file,
                          ...,
                          score="qfdrp",
                          use.sex.chromosomes=FALSE,
                          ignore.strand=TRUE){
  wsh.check.validity(score)
  optlist <- list(...)
  if(length(optlist)==0&!(score%in%c("epipolymorphism","entropy"))){
    stop("Annotation inference only applicable for epipolymorphism and entropy. Please specify annotation.")
  }
  if(length(optlist)>1){
    stop("Only single argument accepted. Either GRanges or RnBSet.")
  }
  if(length(optlist)==0){
    if(score=="epipolymorphism"){
      res <- calculate.epipolymorphism(bam.file)
    }else if(score=="entropy"){
      res <- calculate.entropy(bam.file)
    }else{
      stop("Invalid configuration. Please specify annotation.")
    }
  }else{
    if(inherits(optlist[[1]],"RnBSet")||is.character(optlist[[1]])){
      res <- compute.score.rnb(bam.file=bam.file,rnb.set=optlist[[1]],score=score,use.sex.chromosomes=use.sex.chromosomes,ignore.strand=ignore.strand)
    }else if(inherits(optlist[[1]],"GRanges")){
      res <- compute.score.GRanges(bam.file = bam.file,range=optlist[[1]],score=score,use.sex.chromosomes=use.sex.chromosomes,ignore.strand=ignore.strand)
    }else{
      stop("Invalid value for additional argument, needs to be GRanges or RnBSet")
    }
  }
  return(res)
}

#' wsh.check.validity
#'
#' Checks if the score is compatible with the current setting, or if a non-valid score was specified.
#'
#' @param score Score to be checked
#' @author Michael Scherer
#' @noRd
wsh.check.validity <- function(score){
  if(!(score%in%VALID.SCORES)){
    stop(paste("Invalid value for score, must be one of",VALID.SCORES))
  }
  sys.name <- Sys.info()['sysname']
  if(sys.name=="Windows"&!(score%in%WINDOWS.CAPABLE)){
    stop(paste("Cannot compute score",score,"on windows, only",paste(WINDOWS.CAPABLE,collapse = ", "),"are possible"))
  }
}

#' remove.sex.chromosomes
#'
#' Removes sex chromosomes from the annotation, since we do not compute scores for those chromosomes, only for the autosomes.
#'
#' @param annotation Annotation (RnBSet or GRanges), for which sex chromosomes are to be removed.
#'
#' @return The annotation without the sex chromosomes
#' @noRd
remove.sex.chromosomes <- function(annotation){
  logger.start("Removing Sex chromosomes")
  if(inherits(annotation,"RnBSet")){
    annotation <- rnb.execute.sex.removal(rnb.set=annotation)
  }else if(inherits(annotation,"GRanges")){
    keep <- !(as.character(seqnames(annotation)) %in% c("chrX","X","chrY","Y"))
    annotation <- annotation[keep]
  }else{
    stop("Invalid value for annotation.")
  }
  logger.completed()
  return(annotation)
}

#' create.genomebrowser.track
#'
#' Converts the output of compute.score und compute.score.GRanges into genome browser tracks, that
#' can be directly imported into e.g. UCSC Genome Browser.
#'
#' @param score.output The output as generated through compute.score or compute.score.GRanges: a
#'        data.frame with chromosome positions and the score values for those positions.
#' @param sample.name The name of the sample to be added to the file.
#' @param store.path The path to which the file is to be stored.
#' @param bin.width Width of bins in bp to aggregate the scores. If NULL, the scores are stored
#'        as single CpGs.
#' @author Michael Scherer
#' @export
create.genomebrowser.track <- function(score.output, sample.name="Sample", store.path=getwd(),bin.width=5000){
	if(!is.data.frame(score.output)){
		stop("Invalid value for score.output, needs to be a data.frame")
	}
	if(!all((c("chromosome","start","end") %in% colnames(score.output)))){
		stop("Neccessary column names not available")
	}
  score <- colnames(score.output)[4]
	if(!is.null(bin.width)){
		agg <- lapply(unique(score.output$chromosome),function(chr){
			start.chr <- min(score.output$start[score.output$chromosome %in% chr])
			end.chr <- max(score.output$end[score.output$chromosome %in% chr])
			cbind(as.character(chr),seq(start.chr,end.chr-bin.width,by=bin.width),seq(start.chr+bin.width,end.chr,by=bin.width))
		})
		agg.frame <- c()
		for(chr in agg){
		  agg.frame <- rbind(agg.frame,chr)
		}
		agg.ranges <- GRanges(Rle(agg.frame[,1]),IRanges(start=as.numeric(agg.frame[,2]),as.numeric(agg.frame[,3])))
		score.anno <- GRanges(Rle(score.output$chromosome),IRanges(start=score.output$start,end=score.output$end))
		op <- findOverlaps(score.anno,agg.ranges)
		agg.score <- aggregate(score.output[queryHits(op),4],by=list(subjectHits(op)),mean,na.rm=T)
		score.output <- data.frame(agg.frame[unique(subjectHits(op)),],agg.score$x)
	}
	final.output <- data.frame(score.output[,1:3],paste0("'",round(as.numeric(score.output[,4])*100),"%'"),round(as.numeric(score.output[,4])*100,2))
	final.output <- apply(final.output,c(1,2),as.character)
	to.write <- c()
	for(i in 1:nrow(final.output)){
	  to.write <- c(to.write,paste(paste(c(final.output[i,]),collapse = " ")))
	}
	to.write <- c(paste0("browser position ",final.output[1,1],":",final.output[1,2],"-",final.output[1,3]),
					paste0("track type=bed name=\"",sample.name,"\" description=\"",score, " scores in ", bin.width," tiles\" useScore=1"),
					to.write)
	writeLines(to.write,file.path(store.path,paste0(sample.name,"_",score,".bed")))
}
