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
#  example.bam <- system.file(file.path("extData","small_example.bam"),package="ISH")
#  example.rnb.set <- system.file(file.path("extData","small_rnbSet.zip"),package="ISH")
#  pdr <- compute.score(example.bam,example.rnb.set,score="pdr")
#  passes <- is.numeric(pdr$PDR)
#  checkTrue(passes)
# }

test.options <- function(){
  names.new.options <- c("window.size","mapq.filter","max.reads","min.overlap","fdrp.type","coverage.threshold",
                         "methclone.methylation.diff")
  new.options <- c(window.size=42,mapq.filter=5,max.reads=10,min.overlap=59,fdrp.type="qFDRP",coverage.threshold=5,
                   methclone.methylation.diff=0)
  set.option(window.size=42,mapq.filter=5,max.reads=10,min.overlap=59,fdrp.type="qFDRP",coverage.threshold=5,
             methclone.methylation.diff=0)
  package.options <- get.option(names.new.options)
  passes <- all(package.options == new.options)
  tryCatch(set.option(perl.path="foo"),error=function(e){
    passes <- FALSE
  })
  tryCatch(set.option(samtools.path="foo"),error=function(e){
    passes <- FALSE
  })
  checkTrue(passes)
}

# test.option.influence <- function(){
#   example.bam <- system.file(file.path("extData","small_example.bam"),package="ISH")
#   example.rnb.set <- system.file(file.path("extData","small_rnbSet.zip"),package="ISH")
#   fdrp.default <- rnb.calculate.fdrp(example.rnb.set,example.bam)
#   set.option(coverage.threshold = 50)
#   fdrp.new <- rnb.calculate.fdrp(example.rnb.set,example.bam)
#   passes <- nrow(fdrp.default) != nrow(fdrp.new)
#   set.option(coverage.threshold = 10)
#   fdrp.new <- rnb.calculate.fdrp(example.rnb.set,example.bam)
#   passes <- passes & (nrow(fdrp.default)==nrow(fdrp.new))
#   set.option(mapq.filter = 0)
#   fdrp.new <- rnb.calculate.fdrp(example.rnb.set,example.bam)
#   passes <- passes & (any(fdrp.default$FDRP!=fdrp.new$FDRP))
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
    logger.start("Test package options")
      test.options()
    logger.completed()
    #' Only test locally
    # logger.start("Test package option influence")
    #   test.option.influence()
    # logger.completed()
  logger.completed()
}

execute.unit.test()
