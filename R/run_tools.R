#' run.methclone
#'
#' Runs the methclone software
#'
#' @param bam.file path to bam file containing the reads for defining epialleles
#' @param out.folder path to folder where methclones output should be stored
#' @param out.name name of the output file
#'
#' @author Michael Scherer
#' @export
run.methclone <- function(bam.file,out.folder=getwd(),out.name="methclone"){
  logger.start("Computing epialleles with methclone software")
  location <- system.file(file.path("bin","methclone"),package="ISH")
  out.file <- file.path(out.folder,paste0(out.name,"_tmp.txt.gz"))
  cmd <- paste(location,bam.file,bam.file,out.file,out.name,get.option('methclone.methylation.diff'),get.option('window.size'),get.option('coverage.threshold'))
  logger.info(paste("Executing:",cmd))
  system(command = cmd)
  logger.completed()
}


#' run.haplotype.calculation
#'
#' Runs the haplotype calculation associated with MHL's publication
#'
#' @param roi Region of interest for which haplotype information should be computed
#' @param bam.file path to bam file containing the reads
#' @param out.folder folder to write the output
#' @param out.name name of the output file
#' @param bam.type alignment tools used to create the bam file
#'
#' @author Michael Scherer
#' @export
run.haplotype.calculation <- function(roi,bam.file,out.folder=getwd(),out.name="hapinfo.txt",bam.type="bismark"){
  logger.start("Computing haplotypes with perl scripts (might take several hours/days)")
  Sys.setenv(PATH=paste0(get.option("samtools.path"),":",Sys.getenv("PATH")))
  perl.location <- get.option('perl.path')
  script.location <- system.file(file.path("scripts","mergedBam2hapInfo_RRBS_v1.0.pl"),package = "ISH")
  out.file <- file.path(out.folder,out.name)
  cmd <- paste(perl.location,script.location,roi,bam.file,bam.type,roi,">",out.file)
  logger.info(paste("Exceuting:",cmd))
  system(command = cmd)
  logger.completed()
  return(out.file)
}

#' run.mhl.calculation
#'
#' Runs the MHL calculation given the haplotype information file associated with MHL's publication
#'
#' @param hapinfo.file path to haplotype information file as produced by \code{\link{run.haplotype.calculation}}
#' @param out.folder folder to write the output
#' @param out.name name of the output file
#'
#' @author Michael Scherer
#' @export
run.mhl.calculation <- function(hapinfo.file,out.folder=getwd(),out.name="mhl.txt"){
  logger.start("Computing MHL score from haplotype information")
  perl.location <- get.option('perl.path')
  script.location <- system.file(file.path("scripts","hapinfo2mhl.pl"),package = "ISH")
  out.file <- file.path(out.folder,out.name)
  temp.info <- file.path(out.folder,"info_list.txt")
  writeLines(hapinfo.file,temp.info)
  cmd <- paste(perl.location,script.location,temp.info,">",out.file)
  logger.info(paste("Exceuting:",cmd))
  system(command = cmd)
  logger.completed()
  return(out.file)
}
