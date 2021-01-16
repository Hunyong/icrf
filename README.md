# icrf
Simulations for the interval censored recursive forest (icrf) method by Cho, Jewell, and Kosorok (2020+).  

This set of code is for implementing the simulations and the avalanche victims data analysis.

## I. simulation:  
  * Outline of the main code  
    * I-1. `scripts/1test.R`:          compares different methods (icrf, STIC, SFIC, Cox) under different scenarios.  
    * I-2. `scripts/1B.largeSample.R`: illustrates the performance under different sample sizes.  
    * I-3. `scripts/1C.honesty.R`:     compares quasi-honesty and exploitative approaches.  
    * (Note that the working directory is the root directory. So the scripts should be run without changing the working directory to `scripts/`.)  

  * All the code above is dependent on the following helper functions and code, and they are sourced in the code:  
    * `scripts/0functions.R`:         Elementary helper functions and versions of data-generating functions.  
    * `scripts/1setting.R`:           Given parameters it generates some global variables   
                        (a data-generating function, file names, etc.). 
  
  * These can be done with parallel computing using bash scripts 
      * `scripts/bash-install.sh`           for initial installation of relevant packages  
      * `scripts/bash-icrf-sim-run.sh`      for `1test.R`  
      * `scripts/bash-icrf-sizeSim-run.sh`  for `1B.largeSample.R`   
      * `scripts/bash-icrf-honesty-run.sh`  for `1C.honesty.R`   

  * The parameters for each setting (`scenario`, `sim`, `n.monitor`) are specified by `args`.  
  
  * Once simulations are done, outputs will be generated in the output folders.  
  These will be used in the following analysis code.  
  

## II. analyzing data:  
  * `scripts/2analysis.R`               analyzes the outputs of `scripts/1test.R` and `scripts/1C.honesty.R`.  
  * `scripts/2analysis_largeSample.R`   analyzes the outputs of `scripts/1B.largeSample.R`.  
  
  * These (I and II) are dependent on the following functions, and those files are sourced in the corresponding code (I and II):  
  `scripts/0functions.R` and `scripts/2analysis-functions.R`  
  

## III. Avalanche victims data analysis:  
  * `scripts/5avalanche.R`            analyzes the avalanche victims data.   
  * The data are not publicly available.  
  