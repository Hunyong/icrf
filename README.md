# icrf
Simulations for the interval censored recursive forest (icrf) method by Cho, Jewell, and Kosorok (2020+).  

This set of code is for implementing the simulations and the avalanche victims data analysis.

## I. simulation:  
  * Outline of the main code  
    * I-1. `1test.R`:          compares different methods (icrf, STIC, SFIC, Cox) under different scenarios.  
    * I-2. `1B.largeSample.R`: illustrates the performance under different sample sizes.  
    * I-3. `1C.honesty.R`:     compares quasi-honesty and exploitative approaches.  

  * All the code above is dependent on the following helper functions and code, and they are sourced in the code:  
    * `0functions.R`:         Elementary helper functions and versions of data-generating functions.  
    * `1setting.R`:           Given parameters it generates some global variables   
                        (a data-generating function, file names, etc.). 
  
  * These can be done with parallel computing using bash scripts 
      * `bash-install.sh`           for initial installation of relevant packages  
      * `bash-icrf-sim-run.sh`      for `1test.R`  
      * `bash-icrf-sizeSim-run.sh`  for `1B.largeSample.R`   
      * `bash-icrf-honesty-run.sh`  for `1C.honesty.R`   

  * The parameters for each setting (`scenario`, `sim`, `n.monitor`) are specified by `args`.  
  
  * Once simulations are done, outputs will be generated in the output folders.  
  These will be used in the following analysis code.  
  

## II. analyzing data:  
  * `2analysis.R`               analyzes the outputs of `1test.R` and `1C.honesty.R`.  
  * `2analysis_largeSample.R`   analyzes the outputs of `1B.largeSample.R`.  
  
  * These (I and II) are dependent on the following functions, and those files are sourced in the corresponding code (I and II):  
  `0functions.R` and `2analysis-functions.R`  
  

## III. Avalanche victims data analysis:  
  * `5avalanche.R`            analyzes the avalanche victims data.   
  * The data are not publicly available.  
  