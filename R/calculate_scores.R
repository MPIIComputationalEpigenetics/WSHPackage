########################################################################################################################
## calculate_scores.R
## created: 2017-03-06
## creator: Michael Scherer
## ---------------------------------------------------------------------------------------------------------------------
## Main script of the package. Includes all relevant functions for computing the scores.
########################################################################################################################

#' convert
#' This function converts a sequence representation as a string into a binary (TRUE/FALSE)
#' representation according to the methylation state of the CpG dinucleotide at the given position.
#'
#' @param y position of the CpG nucleotide in the string
#' @param string sequence of the read as a string
#' @return TRUE, if the CpG dinucleotide was methylated at the given position and FALSE otherwise
#'
#' @author Michael Scherer
#' @noRd
convert <- function(y,string){
  if(nchar(string) >= y+2){
    s <- substr(string,y+1,y+2)
    ret <- ifelse(s == 'CG',TRUE,FALSE)
    ret
  }
}

#' toCpGs
#' This function converts the sequence representation as the result of a sequencing technology to a
#' binary representation with TRUE values for the methylated CpGs and FALSE values for the remaining
#' CpGs
#'
#' @param index index of the read that should be converted to the custom representation
#' @param match_read_cpg list containing the mapping from the reads to the CpG sites covered by
#'							the read
#' @param starts_cpgs vector of positions for the CpG sites
#' @param starts_reads starting positions of the reads
#' @param seqs_reads raw sequence of the read
#' @return representation of the read as a logical vector for each CpG site in the read
#'						where TRUE represents a methylated CpG and FALSE a non-methylated (or a sequence
#'						variation)
#'@details This function calls the function convert for each CpG site covered by the read.
#'						Sequence variations and unmethylated CpG sites are treated in the same way.
#'
#' @author Michael Scherer
#' @noRd
toCpGs <- function(index,match_read_cpg,starts_cpgs,starts_reads,seqs_reads){
  covered_cpgs <- match_read_cpg[[index]]
  start_cpgs <- starts_cpgs[covered_cpgs]
  names(start_cpgs) <- start_cpgs
  start_of_read <- starts_reads[index]
  start_cpgs <- start_cpgs-start_of_read
  sequence <- seqs_reads[index]
  representation <- lapply(start_cpgs,convert,sequence)
  representation <- unlist(representation)
  representation
}

#' classify.read
#' This functions classifies a read into either discordant (TRUE) or concordant (FALSE).
#'
#' @param index index of the read that should be converted to the custom representation
#' @param match_read_cpg list containing the mapping from the reads to the CpG sites covered by
#'							the read
#' @param starts_cpgs vector of positions for the CpG sites
#' @param starts_reads starting positions of the reads
#' @param seqs_reads raw sequence of the read
#' @return representation of the read as a logical value, where TRUE represents a discordant
#'                 and FALSE a concordant read
#'@details A read is classified as discordant, if the methylation levels found at all CpGs in the
#'                 read are heterogeneous. That means if there is at least one methylated and one
#'                 umethylated CpG, the read is discordant, otherwise concordant.
#'
#' @noRd
#' @author Michael Scherer
classify.read <- function(index,match_read_cpg,starts_cpgs,starts_reads,seqs_reads){
  covered_cpgs <- match_read_cpg[[index]]
  # This is one of the criteria implemented in the PDR, each read has to contain at least 4 CpGs to
  # be used for the classification into discordant/concordant
  if(length(covered_cpgs)<4){
    return(NA)
  }
  start_cpgs <- starts_cpgs[covered_cpgs]
  names(start_cpgs) <- start_cpgs
  start_of_read <- starts_reads[index]
  start_cpgs <- start_cpgs-start_of_read
  sequence <- seqs_reads[index]
  representation <- lapply(start_cpgs,convert,sequence)
  representation <- unlist(representation)
  #' A read is only concordant if all of the CpGs show the same methylation status
  concordant <- (all(representation) || all(!representation))
  return(!concordant)
}


#' restrict
#' This function restricts a read to those CpG sites that are at most WINDOW.SIZE away from one
#' another.
#'
#' @param positions positions of the CpG sites in the read of interest
#' @param cpg position of the CpG site of interest
#' @return positions of the CpG sites that are at most WINDOW.SIZE away from one
#'						another
#'
#' @noRd
#' @author Michael Scherer
restrict <- function(positions,cpg){
  # We only restrict something, if the read is longer than 50 bp
  if(!any(positions == cpg)) return(NA)
  distances <- abs(as.numeric(positions)-as.numeric(cpg))
  positions <- positions[distances<=get.option('window.size')]
  if((as.numeric(positions)[length(positions)]-as.numeric(positions)[1])>get.option('window.size')){
    distances <- distances[distances<=get.option('window.size')]
    end <- length(distances)
    remove <- rep(FALSE,end)
    pos <- match(cpg,positions)
    if(is.na(pos)){
      return(NA)
    }
    i.left <- pos-1
    i.right <- pos+1
    finished.left <- FALSE
    finished.right <- FALSE
    while((!finished.left) || (!finished.right)){
      left <- distances[i.left]
      right <- distances[i.right]
      if(finished.left){
        i.right <- i.right + 1
        if(i.right > end){
          finished.right <- TRUE
          i.right <- end
          next
        }
        right <- distances[i.right]
        if(get.option('window.size') < (left + right)){
          remove[i.right:end] <- TRUE
          finished.right <- TRUE
          i.right <- max(which(!remove))
        }
      }else if (finished.right){
        i.left <- i.left - 1
        if(i.left < 1){
          finished.left <- TRUE
          i.left <- 1
          next
        }
        left <- distances[i.left]
      }else{
        if((left < right) || all(remove[i.right:end])){
          i.left <- i.left - 1
          if(i.left < 1){
            finished.left <- TRUE
            i.left <- 1
            next
          }
          left <- distances[i.left]
        }else{
          i.right <- i.right + 1
          if(i.right > end){
            finished.right <- TRUE
            i.right <- end
            next
          }
          right <- distances[i.right]
        }
      }
    }
    # select which side to choose
    if((left+right)>get.option('window.size')){
        num.left <- pos - i.left
        num.right <- i.right - pos
        if(num.left<num.right){
            remove[i.left:(pos-1)] <- TRUE
        }else if(num.right<num.left){
            remove[(pos+1):i.right] <- TRUE
        }else{
            select.side <- sample(1:2,1)
            if(select.side==1){
                remove[i.left:(pos-1)] <- TRUE
            }else{
                remove[(pos+1):i.right] <- TRUE
            }
        }
    }
    positions <- positions[!remove]
  }
  positions
}

#' compute.discordant
#' This function decides for a given read pair if it is discordant or not.
#'
#' @param index index of the read pair
#' @param read1 vector of indeces for the query reads for the read pairs
#' @param read2 vector of indeces for the subject reads for the read pairs
#' @param values GRanges object containing the reads from which the read pairs
#'					were calculated
#' @param site position of the CpG site of interest
#' @return classification of the read pair as discordant or concordant
#'
#' @details The given read pair is only concordant, if all CpG sites that overlap between the
#' 				reads reflect the same methylation status. Otherwise the read pair is classified
#' 				as discordant.
#' @noRd
#' @author Michael Scherer
compute.discordant <- function(index,read1,read2,values,site){
  # this is the first read
  q <- read1[index]
  # this is the second
  s <- read2[index]
  # the values are stored in the corresponding part of the GRanges object
  v1 <- values[q]
  # we check if we have any information about the mehtylation state available
  if(length(v1)>0){
    v1 <- v1[[1]]
  }else{
    return(NA)
  }
  v2 <- values[s]
  if(length(v2)>0){
    v2 <- v2[[1]]
  }else{
    return(NA)
  }
  # now we match the corresponding positions to each other
  names1 <- names(v1)
  names2 <- names(v2)
  both <- intersect(names1,names2)
  both <- restrict(both,site)
  if(any(is.na(both))){
    return(NA)
  }
  intersection1 <- v2[both]
  intersection1 <- intersection1[!is.na(intersection1)]
  intersection2 <- v1[both]
  intersection2 <- intersection2[!is.na(intersection2)]
  # two reads are only concordant to each other if they reflect the same methylation
  # state at each position covered by the reads
  # we do a case distinction for either FDRP or qFDRP: while FDRP only classifies each read pair to either
  # discordant or concordant, the qFDRP computes the discordance fraction between the two reads
  if(get.option('fdrp.type')=='FDRP'){
    discordant <- any(intersection1 != intersection2)
  }else{
    discordant <- (sum(intersection1 != intersection2))/length(intersection1)
  }
  discordant
}

#' calculate.fdrp.site
#' This function compute the FDRP score for a given cpg site.
#'
#' @param pos index of the cpg site
#' @param cpg list containing the mapping from each cpg site to the reads that
#'					contain this site
#' @param reads GRanges object containing the reads needed for the calculation
#'					of the FDRP, already converted into the custom representation
#' @param site position of the CpG site of interest as an integer number
#' @return FDRP score for the given CpG site
#'
#' @details This function is called by calculate.fdrps.chromosome each CpG present on
#' 				the chormosome and calls compute.discordant for each pair of reads.
#' @export
#' @author Michael Scherer
calculate.fdrp.site <- function(pos,cpg,reads,site){
  cpg <- cpg[[pos]]
  site <- site[[pos]]
  # we only consider the calculation if there are more than 2 reads for this CpG site
  if(length(cpg)<3){
    return(NA)
  }
  # the maximum number of reads to be calculated is 40, which already are
  # (1/2)*39*40 = 780 read pairs for a single site. Otherwise 40 reads are
  # are sampled from the all reads.
  if(length(cpg)>get.option('max.reads')){
    cpg <- sample(cpg,get.option('max.reads'))
  }
  selected <- reads[cpg]
  rm(reads)

  # here we calculate all pairs of the reads that the corresponding CpG site covers
  # we require the overlap to be at least 35 bp long
  overlap <- findOverlaps(selected,selected,minoverlap=get.option('min.overlap'),ignore.strand=TRUE)
  query <- queryHits(overlap)
  if(length(query)==0){
    return(NA)
  }
  subject <- subjectHits(overlap)
  # we only consider each of the pairs once and the findOverlaps function calculates
  # all pairs in both directions
  smaller <- query<subject
  query <- query[smaller]
  subject <- subject[smaller]

  # now we start with the calculation of the FDRP for the sample
  values <- values(selected)[,"CpG"]
  rm(selected)
  # query is as long as all read pairs we want to consider at this specific position
  ret <- as.list(1:length(query))
  ret <- lapply(ret,compute.discordant,query,subject,values,site)
  ret <- unlist(ret)
  #' we actually calculate the FDRP as frac{#discordant read pairs}{#all read pairs}
  fdrp <- sum(ret,na.rm=TRUE)/length(ret)
  fdrp
}

#' calculate.pdr.site
#' This function computes the PDR score for a given cpg site.
#'
#' @param cpg index of the cpg site
#' @param reads GRanges object containing the reads needed for the calculation
#'					of the PDR, already converted into the custom representation
#' @return PDR score for the given CpG site
#'
#' @export
#' @author Michael Scherer
calculate.pdr.site <- function(cpg,reads){
  # we only consider the calculation if the cpg site is covered by more than 10 reads
  # Another requirement stated in the PDR paper
  if(length(cpg)<10){
    return(NA)
  }

  selected <- reads[cpg]
  rm(reads)

  # We select the reads that contain the given CpG site
  values <- values(selected)[,"isDiscordant"]
  rm(selected)
  # we calculate the PDR as frac{#discordant reads}{#all reads} that contain the given CpG
  values <- unlist(values)
  pdr <- mean(values,na.rm=TRUE)
  pdr
}


#' calculate.fdrp.by.chromosome
#' This function calculates the FDRP scores for the reads given in the bam files
#' for the CpGs present in the annotation.
#'
#' @param bam bam file with the reads from a bisulfite sequencing technique
#'				already aligned to a reference genome
#' @param anno annotation as a GRanges object with the CpG sites to be analyzed
#' @param ignore.strand The \code{ignore.strand} parameter from the function \code{\link{findOverlaps}}
#' @return vector of fdrp scores for the given CpG sites in anno
#'
#' @details This function is called by calculate.fdrps for each chromosome
#'				separately and calls toCpG as well as calculate.fdrp.site.
#'
#' @author Michael Scherer
#'
#' @export
#' @import Rsamtools
#' @import GenomicAlignments
#' @import rtracklayer
calculate.fdrp.by.chromosome <- function(bam,
                                         anno,
                                         ignore.strand=TRUE){
  print(get.option('max.reads'))
  chromosome <- as.character(seqnames(anno))[1]
  is.sex.chromosome <- grepl("X|Y|23|24",chromosome)
  logger.start(paste('Calculation of',chromosome))
  start <- start(ranges(anno)[1])
  end <- end(ranges(anno)[length(anno)])
  if(!(chromosome %in% names(scanBamHeader(bam)[[1]]))){
    if(!is.sex.chromosome){
        chromosome <- unique(na.omit(as.numeric(unlist(strsplit(chromosome,"[^0-9]+")))))
    }else{
        chromosome <- unique(gsub("chr","",chromosome))
    }
  }
  which <- paste0(chromosome,":",start,"-",end)
  which <- GRanges(which)
  param <- ScanBamParam(which=which,what="seq",mapqFilter=get.option('mapq.filter'),flag=scanBamFlag(isNotPassingQualityControls=FALSE,isDuplicate=FALSE))
  reads <- readGAlignments(bam,param=param)

  range_reads <- GRanges(reads)
  rm(reads)
  newStyle <- mapSeqlevels(seqlevels(range_reads),'UCSC')
  newStyle <- newStyle[!is.na(newStyle)]
  range_reads <- renameSeqlevels(range_reads,newStyle)

  # we only analyze those CpGs that are covered (on average) by enough reads in the complete dataset
  overlap <- findOverlaps(range_reads,anno,ignore.strand=ignore.strand)
  query <- queryHits(overlap)
  query <- unique(query)
  range_reads <- range_reads[query]

  ####### REPRESENTATION ############
  # This part converts the raw sequencing reads from the alignment into
  # an representation, where only CpG positions are considered and from whom we
  # can infer discordance or concordance of reads
  logger.start(paste('Representation',chromosome))
  range_cpgs <- ranges(anno)
  starts_cpgs <- start(range_cpgs)
  rm(range_cpgs)
  seqs_reads <- as.character(values(range_reads)$seq)
  starts_reads <- start(ranges(range_reads))
  overlap <- findOverlaps(range_reads,anno,ignore.strand=ignore.strand)
  match_read_cpg <- as(overlap,"list")
  overlap <- findOverlaps(anno,range_reads,ignore.strand=ignore.strand)
  fdrps <- as.list(rep(NA,length(anno)))
  rm(anno)
  # for each read we convert the covered CpG sites into a custom representation
  read_representation <- as.list(1:length(range_reads))
  read_representation <- lapply(read_representation,toCpGs,match_read_cpg,starts_cpgs,starts_reads,seqs_reads)
  rm(match_read_cpg)
  rm(starts_reads)
  rm(seqs_reads)
  values(range_reads) <- DataFrame(cbind(CpG=read_representation))
  rm(read_representation)
  logger.completed()

  ######### FDRP CALCULATION ###########
  # we only calculate the FDRP for the CpGs that are acutally covered by any read in the
  # corresponding sample
  logger.start(paste('FDRP',chromosome))
  match_cpg_reads <- as(overlap,"list")
  rm(overlap)
  null <- lapply(match_cpg_reads,function(x){length(x)>0})
  null <- unlist(null)
  match_cpg_reads <- match_cpg_reads[null]
  starts_cpgs <- starts_cpgs[null]
  toApply <- 1:length(starts_cpgs)
  fdrps_actual <- lapply(toApply,calculate.fdrp.site,match_cpg_reads,range_reads,starts_cpgs)
  # when we do not have a read that covers this site, we set the FDRP for this site to NA
  fdrps[null] <- fdrps_actual
  fdrps <- unlist(fdrps)
  # here we set the name for the corresponding position in the genome
  logger.completed()
  logger.completed()
  fdrps
}

#' calculate.pdr.by.chromosome
#'
#' This function calculates the PDR scores for the reads given in the bam files
#' for the CpGs present in the annotation.
#'
#' @param bam bam file with the reads from a bisulfite sequencing technique
#'				already aligned to a reference genome
#' @param anno annotation as a GRanges object with the CpG sites to be analyzed
#' @param ignore.strand The \code{ignore.strand} parameter from the function \code{\link{findOverlaps}}
#' @return vector of pdr scores for the given CpG sites in anno
#'
#' @author Michael Scherer
#'
#' @export
#' @import Rsamtools
#' @import GenomicAlignments
#' @import rtracklayer
calculate.pdr.by.chromosome <- function(bam, 
                                        anno,
                                        ignore.strand=TRUE){
  chromosome <- as.character(seqnames(anno))[1]
  is.sex.chromosome <- grepl("X|Y|23|24",chromosome)
  logger.start(paste('Calculation of',chromosome))
  start <- start(ranges(anno)[1])
  end <- end(ranges(anno)[length(anno)])
  if(!(chromosome %in% names(scanBamHeader(bam)[[1]]))){
    if(!is.sex.chromosome){
        chromosome <- unique(na.omit(as.numeric(unlist(strsplit(chromosome,"[^0-9]+")))))
    }else{
        chromosome <- unique(gsub("chr","",chromosome))
    }
  }
  which <- paste0(chromosome,":",start,"-",end)
  which <- GRanges(which)
  param <- ScanBamParam(which=which,what="seq",mapqFilter=get.option('mapq.filter'),flag=scanBamFlag(isNotPassingQualityControls=FALSE,isDuplicate=FALSE))
  reads <- readGAlignments(bam,param=param)

  range_reads <- GRanges(reads)
  rm(reads)
  newStyle <- mapSeqlevels(seqlevels(range_reads),'UCSC')
  newStyle <- newStyle[!is.na(newStyle)]
  range_reads <- renameSeqlevels(range_reads,newStyle)

  # we only analyze those CpGs that are covered (on average) by enough reads in the complete dataset
  overlap <- findOverlaps(range_reads,anno,ignore.strand=ignore.strand)
  query <- queryHits(overlap)
  query <- unique(query)
  range_reads <- range_reads[query]

  ####### REPRESENTATION ############
  #This part clasifies all reads into either discordant or concordant
  logger.start(paste('Representation',chromosome))
  range_cpgs <- ranges(anno)
  starts_cpgs <- start(range_cpgs)
  rm(range_cpgs)
  seqs_reads <- as.character(values(range_reads)$seq)
  starts_reads <- start(ranges(range_reads))
  overlap <- findOverlaps(range_reads,anno,ignore.strand=ignore.strand)
  match_read_cpg <- as(overlap,"list")
  overlap <- findOverlaps(anno,range_reads,ignore.strand=ignore.strand)
  pdrs <- as.list(rep(NA,length(anno)))
  rm(anno)
  # we classify each read into either discordant or concordant
  classified_reads <- as.list(1:length(range_reads))
  classified_reads <- lapply(classified_reads,classify.read,match_read_cpg,starts_cpgs,starts_reads,seqs_reads)
  rm(match_read_cpg)
  rm(starts_reads)
  rm(seqs_reads)
  rm(starts_cpgs)
  values(range_reads) <- DataFrame(cbind('isDiscordant'=classified_reads))
  rm(classified_reads)
  logger.completed()

  ######### PDR CALCULATION ###########
  # Starting from the classification of the reads, we now calculate the actual PDR
  # values for the CpG sites of interest
  logger.start(paste('PDR',chromosome))
  match_cpg_reads <- as(overlap,"list")
  rm(overlap)
  null <- lapply(match_cpg_reads,function(x){length(x)>0})
  null <- unlist(null)
  match_cpg_reads <- match_cpg_reads[null]
  pdrs_actual <- lapply(match_cpg_reads,calculate.pdr.site,range_reads)
  # when we do not have a read that covers this site, we set the FDRP for this site to NA
  pdrs[null] <- pdrs_actual
  pdrs <- unlist(pdrs)
  logger.completed()
  logger.completed()
  pdrs
}

#' calculate.fdrp.score
#'
#' This function calculates the FDRPs/qFDRPs (depending on how the options are set) for all CpG sites in
#' the annotation from the reads provided by the bam file.
#'
#' @param bam.file bath to the bam file to be analyzed
#'				already aligned to a reference genome
#' @param anno CpG sites to be analyzed as a GRanges object
#' @param log.path location of the log file
#' @param cores number of cores available for the analysis
#' @param window.size window size used to restrict the concordance/discordance classification of each read pair
#' 						DEFAULT: 50 as the maximum distance
#' @param use.sex.chromosomes Flag indicating if scores are also to be computed for the sex chromosomes
#' @param ignore.strand The \code{ignore.strand} parameter from the function \code{\link{findOverlaps}}
#'
#' @return FDRP scores for the given annotation.
#'
#' @author Michael Scherer
#'
#' @import RnBeads
#' @import doParallel
#' @import parallel
#' @export
calculate.fdrp.score <- function(bam.file,
                                 anno,
                                 log.path=getwd(),
                                 cores=1,
                                 window.size=unname(get.option('window.size')),
                                 max.reads=unname(get.option('max.reads')),
                                 mapq.filter=unname(get.option('mapq.filter')),
                                 coverage.threshold=unname(get.option('coverage.threshold')),
                                 use.sex.chromosomes=FALSE,
                                 ignore.strand=TRUE){
  output.frame <- data.frame(chromosome=seqnames(anno),start=start(anno),end=end(anno))
  bam <- BamFile(bam.file)
  if(!file.exists(file.path(log.path,'log'))){
    dir.create(file.path(log.path,'log'))
  }
  cl <- makeCluster(cores,outfile=file.path(log.path,'log','log_FDRP.log'))
  registerDoParallel(cl)
  anno <- prepare.annotation(anno,use.sex.chromosomes=use.sex.chromosomes)
  if(get.option('fdrp.type')=='FDRP'){
    logger.start("FDRP calculation")
    fdrps <- foreach(chromosome=anno,.combine='c',.packages=c('RnBeads','GenomicAlignments','Rsamtools','rtracklayer'),.export=c('calculate.fdrp.by.chromosome',
        'max.reads',
        'calculate.fdrp.site',
        'compute.discordant',
        'mapq.filter',
        'window.size',
        'coverage.threshold',
        'toCpGs',
        'convert',
        'restrict',
        'set.option',
        'get.option',
        'IHS.OPTIONS')) %dopar%{
      set.option(mapq.filter=mapq.filter)
      set.option(coverage.threshold=coverage.threshold)
      set.option(window.size=window.size)
      set.option(max.reads=max.reads)
      set.option('fdrp.type'='FDRP')
      calculate.fdrp.by.chromosome(bam,chromosome,ignore.strand = ignore.strand)
    }
    logger.completed()
  }else{
    logger.start("qFDRP calculation")
    fdrps <- foreach(chromosome=anno,.combine='c',.packages=c('RnBeads','GenomicAlignments','Rsamtools','rtracklayer'),.export=c('calculate.fdrp.by.chromosome',
        'max.reads',
        'calculate.fdrp.site',
        'compute.discordant',
        'mapq.filter',
        'coverage.threshold',
        'window.size',
        'toCpGs',
        'convert',
        'restrict',
        'set.option',
        'get.option',
        'IHS.OPTIONS')) %dopar%{
      set.option(mapq.filter=mapq.filter)
      set.option(coverage.threshold=coverage.threshold)
      set.option(window.size=window.size)
      set.option(max.reads=max.reads)
      set.option('fdrp.type'='qFDRP')
      calculate.fdrp.by.chromosome(bam,chromosome,ignore.strand = ignore.strand)
    }
    logger.completed()
  }
  stopCluster(cl)
  fdrps <- unlist(fdrps)
  output.frame$FDRP <- fdrps
  return(output.frame)
}

#' calculate.fdrp
#'
#' This function calculates the FDRP scores for the reads given in the bam files
#' for the CpGs present in the annotation.
#'
#' @param bam.file bath to the bam file to be analyzed
#'				already aligned to a reference genome
#' @param anno annotation as a GRanges object with the CpG sites to be analyzed
#' @param log.path location of the log file
#' @param cores number of cores available for the analysis
#' @param window.size window size used to restrict the concordance/discordance classification of each read pair
#' 						DEFAULT: 50 as the maximum distance
#' @param use.sex.chromosomes Flag indicating if scores are also to be computed for the sex chromosomes
#' @param ignore.strand The \code{ignore.strand} parameter from the function \code{\link{findOverlaps}}
#'
#' @return FDRP scores for the given annotation.
#'
#' @author Michael Scherer
#' @export
calculate.fdrp <- function(bam.file,
                           anno,
                           log.path=getwd(),
                           cores=1,
                           window.size=get.option('window.size'),
                           use.sex.chromosomes=FALSE,
                           ignore.strand=TRUE){
  set.option(fdrp.type='FDRP')
  if(!use.sex.chromosomes){
    anno <- remove.sex.chromosomes(anno)
  }
  qfdrp <- calculate.fdrp.score(bam.file,anno,log.path,cores,use.sex.chromosomes = use.sex.chromosomes,ignore.strand = ignore.strand)
  return(qfdrp)
}

#' calculate.qfdrp
#'
#' This function calculates the qFDRP scores for the reads given in the bam files
#' for the CpGs present in the annotation.
#'
#' @param bam.file bath to the bam file to be analyzed
#'				already aligned to a reference genome
#' @param anno annotation as a GRanges object with the CpG sites to be analyzed
#' @param log.path location of the log file
#' @param cores number of cores available for the analysis
#' @param window.size window size used to restrict the concordance/discordance classification of each read pair
#' 						DEFAULT: 50 as the maximum distance
#' @param use.sex.chromosomes Flag indicating if scores are also to be computed for the sex chromosomes
#' @param ignore.strand The \code{ignore.strand} parameter from the function \code{\link{findOverlaps}}
#'
#' @return qFDRP scores for the given annotation.
#'
#' @author Michael Scherer
#' @export
calculate.qfdrp <- function(bam.file,
                            anno,
                            log.path=getwd(),
                            cores=1,
                            window.size=get.option('window.size'),
                            use.sex.chromosomes=FALSE,
                            ignore.strand=TRUE){
  set.option(fdrp.type='qFDRP')
  if(!use.sex.chromosomes){
    anno <- remove.sex.chromosomes(anno)
  }
  qfdrp <- calculate.fdrp.score(bam.file,anno,log.path,cores,use.sex.chromosomes=use.sex.chromosomes,ignore.strand = ignore.strand)
  colnames(qfdrp)[ncol(qfdrp)] <- "qFDRP"
  return(qfdrp)
}

#' calculate.pdr
#'
#'  This function calculates the PDRs for all analyzed CpG sites in
#' the bam file of the corresponding sample
#'
#' @param bam.file bath to the bam file to be analyzed
#'				already aligned to a reference genome
#' @param anno	annotation as a GRanges object with the CpG sites to be analyzed
#' @param log.path location of the log file
#' @param cores number of cores available for the analysis
#' @param use.sex.chromosomes Flag indicating if scores are also to be computed for the sex chromosomes
#' @param ignore.strand The \code{ignore.strand} parameter from the function \code{\link{findOverlaps}}
#'
#' @return PDR scores for the given annotation.
#'
#' @author Michael Scherer
#'
#' @import RnBeads
#' @import doParallel
#' @import parallel
#' @export
calculate.pdr <- function(bam.file,
                          anno,
                          log.path=getwd(),
                          cores=1,
                          window.size=unname(get.option('window.size')),
                          max.reads=unname(get.option('max.reads')),
                          mapq.filter=unname(get.option('mapq.filter')),
                          coverage.threshold=unname(get.option('coverage.threshold')),
                          use.sex.chromosomes=FALSE,
                          ignore.strand=TRUE){
  logger.start("PDR calculation")
  if(!use.sex.chromosomes){
    anno <- remove.sex.chromosomes(anno)
  }  
  output.frame <- data.frame(chromosome=seqnames(anno),start=start(anno),end=end(anno))
  bam <- BamFile(bam.file)
  if(!file.exists(file.path(log.path,'log'))){
    dir.create(file.path(log.path,'log'))
  }
  cl <- makeCluster(cores,outfile=file.path(log.path,'log','log_PDR.log'))
  registerDoParallel(cl)
  anno <- prepare.annotation(anno,use.sex.chromosomes=use.sex.chromosomes)
  pdrs <- foreach(chromosome=anno,.combine='c',.packages=c('RnBeads','GenomicAlignments','Rsamtools','rtracklayer'),.export=c('calculate.pdr.by.chromosome','calculate.pdr.site','classify.read',
        'max.reads',
        'mapq.filter',
        'coverage.threshold',
        'window.size',
        'toCpGs',
        'convert',
        'restrict',
        'set.option',
        'get.option',
        'IHS.OPTIONS')) %dopar%{
      set.option(mapq.filter=mapq.filter)
      set.option(coverage.threshold=coverage.threshold)
      set.option(window.size=window.size)
      set.option(max.reads=max.reads)
    calculate.pdr.by.chromosome(bam,chromosome,ignore.strand = ignore.strand)
  }
  stopCluster(cl)
  pdrs <- unlist(pdrs)
  output.frame$PDR <- pdrs
  logger.completed()
  return(output.frame)
}

#' calculate.epipoly.line
#'
#' This functions computes epipolymorphim for a single line of methclone's output.
#'
#' @param line single line in methclone's output
#' @return epipolymorphism for this region
#'
#' @author Michael Scherer
#' @noRd
calculate.epipoly.line <- function(line){
  percentages <- as.numeric(line[15:30])
  epipoly <- 1-(sum((percentages/100)^2))
  epipoly <- round(epipoly,4)
  epipoly
}


#' calculate.epipolymorphism
#'
#' This function calculates Epipolymorphism by calling the methclone software to compute epiallele frequencies and
#' then uses the definition of epipolymorphism for the calculation.
#'
#' @param bam.file path to bam file containing the reads of the data set
#' @param out.folder folder to store the temporary file produced by methclone
#' @param out.name name of temporary file
#'
#' @return a data frame containing the position in the reference genome and the Epipolymorphism scores
#'
#' @import RnBeads
#' @export
calculate.epipolymorphism <- function(bam.file,out.folder=getwd(),out.name="methclone"){
  logger.start("Epipolymorphism calculation")
  run.methclone(bam.file=bam.file,out.folder=out.folder,out.name=out.name)
  methclone.file <- file.path(out.folder,paste0(out.name,"_tmp.txt.gz"))
  methclone.data <- read.csv(methclone.file,sep='\t')
  system(paste('rm',methclone.file))
  output.frame <- data.frame(chromosome=methclone.data$chr,start=methclone.data$start,end=methclone.data$end,strand=methclone.data$strand)
  values <- apply(methclone.data,1,calculate.epipoly.line)
  output.frame <- data.frame(output.frame,Epipolymorphism=values)
  logger.completed()
  return(output.frame)
}

#' calculate.entropy.line
#'
#' This functions computes methylation entropy for a single line of methclone's output.
#'
#' @param line single line in methclone's output
#' @return entropy for this region
#'
#' @author Michael Scherer
#' @noRd
calculate.entropy.line <- function(line){
  percentages <- as.numeric(line[15:30])
  percentages <- percentages[percentages!=0]
  entropy <- -0.25*(sum((percentages/100)*log2(percentages/100)))
  entropy <- round(entropy,4)
  entropy
}


#' calculate.entropy
#'
#' This function calculates Methylation entropy by calling the methclone software to compute epiallele frequencies and
#' then uses the definition of entropy for the calculation.
#'
#' @param bam.file path to bam file containing the reads of the data set
#' @param out.folder folder to store the temporary file produced by methclone
#' @param out.name name of temporary file
#'
#' @return a data frame containing the position in the reference genome and the Entropy scores
#'
#' @import RnBeads
#' @export
calculate.entropy <- function(bam.file,out.folder=getwd(),out.name="methclone"){
  logger.start("Entropy calculation")
  run.methclone(bam.file=bam.file,out.folder=out.folder,out.name=out.name)
  methclone.file <- file.path(out.folder,paste0(out.name,"_tmp.txt.gz"))
  methclone.data <- read.csv(methclone.file,sep='\t')
  system(paste('rm',methclone.file))
  output.frame <- data.frame(chromosome=methclone.data$chr,start=methclone.data$start,end=methclone.data$end,strand=methclone.data$strand)
  values <- apply(methclone.data,1,calculate.entropy.line)
  output.frame <- data.frame(output.frame,Entropy=values)
  logger.completed()
  return(output.frame)
}

#' calculate.mhl
#'
#' This function computed the Methylation Haplotype Load (MHL) for all sites in a given annotation by employing the scripts
#' downloaded from the paper's website
#'
#' @param roi Region of interest for which haplotype information should be computed
#' @param bam.file path to bam file containing the reads
#' @param out.folder folder to write the output
#' @param out.name name of the output file
#' @param bam.type alignment tools used to create the bam file
#'
#' @return a data frame containing the position in the reference genome and the MHL scores
#' @export
calculate.mhl <- function(roi,bam.file,out.folder=getwd(),out.name="mhl.txt",bam.type="bismark"){
  logger.start("MHL calculation")
  hapinfo.file <- run.haplotype.calculation(roi,bam.file,out.folder,bam.type)
  mhl.file <- run.mhl.calculation(hapinfo.file,out.folder,out.name)
  mhl.data <- read.table(mhl.file,sep='\t',skip=1)
  anno.mhl <- as.character(mhl.data[,1])
  anno.mhl <- strsplit(anno.mhl,'[[:punct:]]')
  anno.mhl <- as.data.frame(anno.mhl)
  anno.mhl <- t(anno.mhl)
  mhl.scores <- as.numeric(mhl.data[,2])
  out.frame <- data.frame(chromosome=anno.mhl[,1],start=anno.mhl[,2],end=anno.mhl[,3],MHL=mhl.scores)
  row.names(out.frame) <- 1:nrow(out.frame)
  logger.completed()
  return(out.frame)
}
