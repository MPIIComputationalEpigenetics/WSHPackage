# WSH
[![Build Status](https://travis-ci.org/MPIIComputationalEpigenetics/WSHPackage.svg?branch=master)](https://travis-ci.org/schmic05/WSH_package)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI: 10.1093/nar/gkaa120](https://zenodo.org/badge/doi/10.1093/nar/gkaa120.svg)](https://doi.org/10.1093/nar/gkaa120)

## Short description
R package for the calculation of the following Intra-Sample Heterogeneity Scores in Bisulfite Sequencing Data: FDRP, qFDRP, PDR, Epipolymorphism, Methylation Entropy and MHL. The package is distributed under the GNU GPL-3 license, except for the external tools located in *inst/bin* (methclone, GNU Lesser GPL) and *inst/scripts* (MHL scripts, software distribution granted in the source file). All remaining code may be used, copied, changed and further distributed.

## Installation
```r
if(!requireNamespace("devtools")) install.packages("devtools",repos = "https://cloud.r-project.org/")
devtools::install_github("MPIIComputationalEpigenetics/WSHPackage")
```

## Documentation
```r
example.bam <- system.file(file.path("extData","small_example.bam"),
                           package="WSH")
example.rnb.set <- system.file(file.path("extData","small_rnbSet.zip"),
                               package="WSH")
example.rnb.set <- load.rnb.set(example.rnb.set)
set.option(coverage.threshold = 10)
qfdrp <- rnb.calculate.qfdrp(example.rnb.set,example.bam)
```

A more detailed documentation of the WSH R-package is available [here](vignettes/WSH.md).

## Citation
If you are using this package, please cite:

- Scherer M., et al., Quantitative comparison of within-sample heterogeneity scores for DNA methylation data , Nucleic Acids Research, 2020, [10.1093/nar/gkaa120](https://doi.org/10.1093/nar/gkaa120)

## Contact
You can contact [Michael Scherer](mailto:mscherer@mpi-inf.mpg.de) for reporting bugs, feature requests or questions.
