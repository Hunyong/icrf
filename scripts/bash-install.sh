module load r/3.5.2
Rscript -e "install.packages('icrf')"
Rscript -e "install.packages('icenReg',  repos = 'https://cloud.r-project.org')"
#Rscript -e "install.packages('icrf_2.0.0.tar.gz', repos = NULL, INSTALL_opts = c('--no-lock'))"
#Rscript -e "install.packages('LTRCtrees',  repos = 'https://cloud.r-project.org')"
#Rscript -e "install.packages('BiocManager')"
#Rscript -e "BiocManager::install('Icens')"