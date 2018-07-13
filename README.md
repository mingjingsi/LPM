LPM
===

LPM (Latent Probit Model), is an efficient statistical approach to characterize relationship among complex traits using summary statistics from multiple GWASs and functional annotations. 'LPM' package provides model parameter estimation as well as statistical inference.

Installation
===========

To install the development version of LPM, it's easiest to use the 'devtools' package. Note that LPM depends on the 'Rcpp' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

```
#install.packages("devtools")
library(devtools)
install_github("mingjingsi/LPM")
```

Windows users need to download [local323](http://www.stats.ox.ac.uk/pub/Rtools/goodies/multilib/local323.zip), unzip it in C:/local323 and set the Windows parameter LIB-GSL to be C:\local323.

Usage
===========

[The 'LPM' vignette](https://github.com/mingjingsi/LPM/blob/master/inst/doc/LPM_package.pdf?raw=true) will provide a good start point for the genetic analysis using LPM package. The following help page will also provide quick references for LPM package and the example command lines:

```
library(LPM)
package?LPM
```

References
==========

Jingsi Ming, Tao Wang and Can Yang; LPM: a latent probit model to characterize relationship among complex traits using summary statistics from multiple GWASs and functional annotations.


<!Reproducibility
==========

All the simulation results can be reproduced by using the code at [sim-LPM](https://github.com/mingjingsi/sim-LPM).>


Development
==========

This R package is developed by Jingsi Ming and Can Yang (macyang@ust.hk)
