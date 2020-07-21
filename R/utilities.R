########################################################################################################################
## utilities.R
## created: 2020-07-17
## creator: Michael Scherer
## ---------------------------------------------------------------------------------------------------------------------
## Utility functions for the WSH package
########################################################################################################################

#' anno.split
#' This function splits an annotation as a GRanges object into a GRangesList
#' with two entries.
#'
#' @param anno annotation as a GRanges object to be split into two annotations
#' @return GRangesList object with two more or less equally sized parts
#'
#' @details the first with lower starting points as the half of the
#' maximum starting number and the second with the remaining. In the best case,
#' the two annotations are equally sized in the end.
#'
#' @author Michael Scherer
#' @noRd
anno.split <- function(anno){
  start <- start(ranges(anno)[1])
  end <- end(ranges(anno)[length(anno)])
  end <- round(end/2,0)
  part1 <- anno[end(ranges(anno))<=end]
  part2 <- anno[start(ranges(anno))>end]
  return(GRangesList(part1,part2))
}

#' prepare.annotation
#' 
#' This function prepares the annotation for the use for WSH score computation.
#' 
#' @param anno The annotation that is to be prepared
#' @param use.sex.chromosomes Flag indiciating if the sex chromosomes are to be included in the anysis
#' @return The prepared annotation as a \code{GRangesList}
#' 
#' @author Michael Scherer
#' @noRd
prepare.annotation <- function(anno,
                               use.sex.chromosomes){
  sel.chromosomes <- ifelse(use.sex.chromosomes,"[1-9]|X|Y","[1-9]")
  sel.chromosomes <- grepl(sel.chromosomes,seqnames(anno))
  anno <- anno[sel.chromosomes]
  anno <- split(anno,seqnames(anno))
  anno <- anno[lengths(anno)>0]
  if(length(anno)>2){
    first <- anno.split(anno[[1]])
    anno <- c(first,anno[2:length(anno)])
    second <-  anno.split(anno[[3]])
    anno <- c(anno[1:2],second,anno[4:length(anno)])
  }  
  return(anno)
}  