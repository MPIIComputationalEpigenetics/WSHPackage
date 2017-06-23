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
  location <- system.file(file.path("bin","methclone"),package="IHS")
  cmd <- paste(location,bam.file,bam.file,out.folder,get.option('methclone.methylation.diff'),get.option('window.size'),get.option('coverage.threshold'))
  print(cmd)
  system(command = cmd)
}
