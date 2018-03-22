#' unit testing for ISH R package

#' checks if the example is correctly working
test.example <- function(){
  qfdrp <- ish.run.example()$qFDRP
  fdrp <- ish.run.example("fdrp")$FDRP
  pdr <- ish.run.example("pdr")$PDR
  #mhl <- ish.run.example("mhl")
  epipoly <- ish.run.example("epipolymorphism")$Epipolymorphism
  entropy <- ish.run.example("entropy")$Entropy
  passes <- is.numeric(qfdrp)&is.numeric(fdrp)&is.numeric(pdr)&is.numeric(epipoly)&is.numeric(entropy)
  checkTrue(passes)
}

#' tests function to compute ISH scores from GRanges objects
test.GRanges <- function(){
  example.bam <- system.file(file.path("extData","small_example.bam"),package="ISH")
  example.GRanges <- GRanges(Rle(rep("chr2",10)),IRanges(start = c(2298361,2298554,2298732,2298743,2298787,2298792,2298827,2298884,
                                                                   2298915,2298921),end=c(2298361,2298554,2298732,2298743,2298787,
                                                                                          2298792,2298827,2298884,2298915,2298921)+1))
  qfdrp <- compute.score(example.bam,example.GRanges)
  passes <- is.numeric(qfdrp$qFDRP)
  checkTrue(passes)
}

#' tests function to compute ISH scores from RnBSet objects
# test.rnbSet <- function(){
#   example.bam <- system.file(file.path("extData","small_example.bam"),package="ISH")
#   example.rnb.set <- system.file(file.path("extData","small_rnbSet.zip"),package="ISH")
#   pdr <- compute.score(example.bam,example.rnb.set,score="pdr")
#   passes <- is.numeric(pdr$PDR)
#   checkTrue(passes)
# }


#' main testing function
execute.unit.test <- function(){
  require("RUnit")
  logger.start("Unit testing")
    logger.start("Test example")
      test.example()
    logger.completed()
    logger.start("Test GRanges function")
      test.GRanges()
    logger.completed()
    #' Only test locally
    # logger.start("Test RnBSet function")
    #   test.rnbSet()
    # logger.completed()
  logger.completed()
}

execute.unit.test()
