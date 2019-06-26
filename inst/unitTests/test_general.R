#' unit testing for WSH R package

#' checks if the example is correctly working
test.example <- function(){
  qfdrp <- wsh.run.example()$qFDRP
  fdrp <- wsh.run.example("fdrp")$FDRP
  pdr <- wsh.run.example("pdr")$PDR
  #mhl <- wsh.run.example("mhl")
  epipoly <- tryCatch(wsh.run.example("epipolymorphism")$Epipolymorphism,error=function(e){
    if(Sys.info()['sysname']=="Windows"){
     return(0)
    }
  })
  entropy <- tryCatch(wsh.run.example("entropy")$Entropy,error=function(e){
    if(Sys.info()['sysname']=="Windows"){
      return(0)
    }
  })
  passes <- is.numeric(qfdrp)&is.numeric(fdrp)&is.numeric(pdr)&is.numeric(epipoly)&is.numeric(entropy)
  checkTrue(passes)
}

#' tests function to compute WSH scores from GRanges objects
test.GRanges <- function(){
  example.bam <- system.file(file.path("extData","small_example.bam"),package="WSH")
  example.GRanges <- GRanges(Rle(rep("chr2",10)),IRanges(start = c(2298361,2298554,2298732,2298743,2298787,2298792,2298827,2298884,
                                                                   2298915,2298921),end=c(2298361,2298554,2298732,2298743,2298787,
                                                                                          2298792,2298827,2298884,2298915,2298921)+1))
  qfdrp <- compute.score(example.bam,example.GRanges)
  passes <- is.numeric(qfdrp$qFDRP)
  checkTrue(passes)
}

#' tests function to compute WSH scores from RnBSet objects
 test.rnbSet <- function(){
  example.bam <- system.file(file.path("extData","small_example.bam"),package="WSH")
  example.rnb.set <- system.file(file.path("extData","small_rnbSet.zip"),package="WSH")
  pdr <- compute.score(example.bam,example.rnb.set,score="pdr")
  passes <- is.numeric(pdr$PDR)
  checkTrue(passes)
 }

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

 test.option.influence <- function(){
   example.bam <- system.file(file.path("extData","small_example.bam"),package="WSH")
   example.rnb.set <- system.file(file.path("extData","small_rnbSet.zip"),package="WSH")
   fdrp.default <- rnb.calculate.fdrp(example.rnb.set,example.bam)
   set.option(coverage.threshold = 50)
   fdrp.new <- rnb.calculate.fdrp(example.rnb.set,example.bam)
   passes <- nrow(fdrp.default) != nrow(fdrp.new)
   set.option(coverage.threshold = 10)
   fdrp.new <- rnb.calculate.fdrp(example.rnb.set,example.bam)
   passes <- passes & (nrow(fdrp.default)==nrow(fdrp.new))
   set.option(mapq.filter = 0)
   fdrp.new <- rnb.calculate.fdrp(example.rnb.set,example.bam)
   passes <- passes & (any(fdrp.default$FDRP!=fdrp.new$FDRP))
   checkTrue(passes)
 }

test.genomebrowser <- function(){
	qfdrp <- wsh.run.example()
	create.genomebrowser.track(qfdrp)
	create.genomebrowser.track(qfdrp,bin.width=NULL)
	res <- readLines("Sample_qFDRP.bed")
	unlink("Sample_qFDRP.bed")
	passes <- length(res) > 0
	checkTrue(passes)
}

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
    # Only test locally
    # logger.start("Test RnBSet function")
    #    test.rnbSet()
    # logger.completed()
    logger.start("Test package options")
      test.options()
    logger.completed()
    # Only test locally
    # logger.start("Test package option influence")
    #   test.option.influence()
    # logger.completed()
	logger.start("Test Genome Browser conversion")
		test.genomebrowser()
	logger.completed()
  logger.completed()
}

execute.unit.test()
