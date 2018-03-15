#' unit testing for ISH R package

#' checks if the example is correctly working
test.example <- function(){
  qfdrp <- ish.run.example()
  fdrp <- ish.run.example("fdrp")
  pdr <- ish.run.example("pdr")
  mhl <- ish.run.example("mhl")
  epipoly <- ish.run.example("epipolymorphism")
  entropy <- ish.run.example("entropy")
  passes <- is.numeric(qfdrp)&is.numeric(fdrp)&is.numeric(pdr)&is.numeric(mhl)&is.numeric(epipoly)&is.numeric(entropy)
  checkTrue(passes)
}

#' main testing function
execute.unit.test <- function(){
  require("RUnit")
  logger.start("Unit testing")
    logger.start("Test example")
      test.example()
    logger.completed()
  logger.completed()
}

execute.unit.test()
