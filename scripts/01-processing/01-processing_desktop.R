setwd('~/data/shanghai_clinical_scRNAseq/')


#set lib.loc
#.libPaths('/gs/gsfs0/users/aferrena/R/x86_64-pc-linux-gnu-library/4.0')


#increase globals size, for integration step... uses future under the hood...
# i believe this increases the maximum RAM used for each parallel "worker" thread.
# the error reporting for this is quite verbose thankfully
#options(future.globals.maxSize = 1000 * 1024^2) --> 1gb
#options(future.globals.maxSize = 5000 * 1024^2) #--> 5gb
options(future.globals.maxSize = 10000 * 1024^2) #--> 10gb


#close any existing plotting devices
while (!is.null(dev.list()))  dev.off()


#try to check ram
message('check total RAM')
print(benchmarkme::get_ram())

#check cores
message('check cores')
parallel::detectCores()

library(tidyverse)
library(RISC)
library(Seurat)
#library(FerrenaSCRNAseq)

#library(biomaRt)


set.seed(2020)
