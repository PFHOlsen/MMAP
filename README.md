# MMAP
Implementation of robust version of EM algorithm for Markovian Arrival Processes (MAP) from [Breuer (2002)](https://link.springer.com/article/10.1023/A:1020981005544)
and new EM algorithm for Marked Markovian Arrival Processes (MMAP). A numerical example with MAP-fitting and MMAP-fitting to forest fire data is provided.

## Files
**1.** `EM_algorithms.cpp` -  C++ implementation of robust EM algorithms for MAP and MMAP.\
**2.** `MAP_utils.R` - R code with helper functions to generate MAP parameters and compute conditional path-dependent distribution of data generating Markov process.\
**3.** `data.R` - R code with preprocessing of forest fire data from Oregon during 2018-2022. Returns arrival sequences for all regions and the regions EOA, NOA, SOA respectively.  \
**4.** `fitting.R` - R code with statistical fitting of MAP and MMAP to forest fire data using EM algorithms.\
**5.** `diagnostics.R` - R code with construction of uniform pseudo-samples and QQ-plots.\
\
Note: The files should be run in the order assigned above.

## Dependencies
This project requires the following R packages:
- matrixdist
- Rcpp
- RcppArmadillo
- tidyverse
  
Install them with:
```r
install.packages(c("matrixdist","Rcpp", "RcppArmadillo","tidyverse"))
``` 
