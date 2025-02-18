#######################################################
# INSTALLING BIOGEOBEARS
#######################################################
# Run the install.packages commands ONCE

# Please run them one at a time.

# Installing devtools
install.packages("devtools", type="binary", repos="https://cloud.r-project.org")

# 
# IF YOU GET A MESSAGE LIKE THIS, TYPE "n" FOR NO:
#  There are binary versions available but the source versions are later:
#           binary source needs_compilation
# jsonlite   1.8.3  1.8.4              TRUE
# htmltools  0.5.3  0.5.4              TRUE
# 
# Do you want to install from sources the packages which need compilation? (Yes/no/cancel)
# 

install.packages("ape", type="binary", repos="https://cloud.r-project.org")
install.packages("Rcpp", type="binary", repos="https://cloud.r-project.org")
install.packages("ape", type="binary", repos="https://cloud.r-project.org")
install.packages("FD", type="binary", repos="https://cloud.r-project.org")
install.packages("snow", type="binary", repos="https://cloud.r-project.org")
install.packages("phytools", type="binary", repos="https://cloud.r-project.org")
install.packages("phangorn", type="binary", repos="https://cloud.r-project.org")
install.packages("phylobase", type="binary", repos="https://cloud.r-project.org")
install.packages("optimx", type="binary", repos="https://cloud.r-project.org")
install.packages("GenSA", type="binary", repos="https://cloud.r-project.org")

# R packages by Nick Matzke -- dependencies of BioGeoBEARS
install.packages("rexpokit", type="binary", repos="https://cloud.r-project.org")
install.packages("cladoRcpp", type="binary", repos="https://cloud.r-project.org")


# Install BioGeoBEARS from GitHub
# (BioGeoBEARS is pure R, so installation is easy *if* the above 
#  packages have been installed)
library(devtools)
devtools::install_github(repo="nmatzke/BioGeoBEARS", INSTALL_opts="--byte-compile", upgrade="never")

##############################################################################################
##############################################################################################
### SE OS PROCEDIMENTOS ACIMA NÃO FUNCIONAREM CORRETAMENTE... ###################################
install.packages("devtools")
library(devtools)

install.packages("Rcpp")
install.packages("ape")
install.packages("phytools")
install.packages("FD")
install.packages("GenSA")
install.packages("snow")
install.packages("rexpokit")
install.packages("cladoRcpp")

library(devtools)

# BioGeoBEARS now lives on GitHub instead of CRAN
# To install
# BioGeoBEARS version 1.1 from GitHub, install with:

devtools::install_github(repo="nmatzke/BioGeoBEARS")


# Check that your installed packages will load:
library(ape)
library(phytools)
library(GenSA)
library(FD)
library(snow)
library(parallel)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

#######################################################

