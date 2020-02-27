LPM
===

LPM (Latent Probit Model), is an efficient statistical approach to characterize relationship among complex traits using summary statistics from multiple GWASs and functional annotations. 'LPM' package provides model parameter estimation as well as statistical inference.

To make LPM scalable to millions of SNPs and hundreds of traits, instead of working with a bruteforce algorithm to handle all the data simultaneously, we develope an efficient parameter-expanded EM (PX-EM) algorithm for pair-wise analysis and implement a dynamic threading strategy to enhance its parallel property. This pairwise strategy is guaranteed to give consistent results by our theoretical analysis from the perspective of the composite likelihood approach.

Installation
===========

To install the development version of LPM, it's easiest to use the 'devtools' package. Note that LPM depends on the 'Rcpp' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

```
#install.packages("devtools")
library(devtools)
install_github("mingjingsi/LPM")
```

Windows users need to download [local323](http://www.stats.ox.ac.uk/pub/Rtools/goodies/multilib/local323.zip), unzip it in C:/local323 and set the Windows parameter LIB-GSL to be C:\local323. Mac users need to install GNU Scientific Library (GSL).

Usage
===========

[The 'LPM' vignette](https://github.com/mingjingsi/LPM/blob/master/inst/doc/LPM_package.pdf?raw=true) will provide a good start point for the genetic analysis using LPM package. The following help page will also provide quick references for LPM package and the example command lines:

```
library(LPM)
package?LPM
```

References
==========

Jingsi Ming, Tao Wang, Can Yang, LPM: a latent probit model to characterize the relationship among complex traits using summary statistics from multiple GWASs and functional annotations, Bioinformatics, , btz947, https://doi.org/10.1093/bioinformatics/btz947

[Supplementary Document](https://github.com/mingjingsi/LPM/blob/master/suppl_doc/LPM%20-%20Supplement.pdf?raw=true)


Reproducibility
==========

All the simulation results can be reproduced by using the code at [sim-LPM](https://github.com/mingjingsi/sim-LPM). Real data sets used in the paper have been made publicly available, including functional annotations ([link](https://drive.google.com/file/d/1Jn_MEDJVZR16UB3lXFjNjAuQDCKECfUw/view)) and summary statistics from GWAS ([link](https://drive.google.com/file/d/1vXZX1l5IWXEx9De0psprR5GCniI3VG0l/view)).


Development
==========

This R package is developed by Jingsi Ming and Can Yang (macyang@ust.hk)
