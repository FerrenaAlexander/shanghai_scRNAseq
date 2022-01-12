#time things:
sys_start_time <- Sys.time()

library(tidyverse)
#library(RISC)
library(Seurat)
library(FerrenaSCRNAseq)
library(DoubletFinder)


set.seed(2020)

setwd('~/data/shanghai_clinical_scRNAseq/')


#close any existing plotting devices
while (!is.null(dev.list()))  dev.off()


#try to check ram
message('check total RAM: ')
print(benchmarkme::get_ram())

#check cores
message('check cores: ')
parallel::detectCores()







#indir and outdir
indir <- 'data/rawdata_GEO/unzipped/'
outdir <- 'data/processed/'


#get sample names and sort them
samps <- list.files(indir)
samps <- str_sort(samps, numeric=T)

#append path
samppaths <- paste0(indir, samps)


for(sampdex in c(1:length(samps)) ){
  
  #get samples
  samp <- samps[sampdex]
  samppath <- samppaths[sampdex]
  
  #set output dirs
  outputdir <- paste0(outdir, samp)
  dir.create(outputdir)
  dir.create( paste0(outputdir, '/qc') )
  
  message('\nReading in ', samp, '\n')
  
  
  sobj <- CreateSeuratObject(   Read10X(samppath), min.cells= 3, project = samp)
  
  #mito content, add to metadata
  mito.features <- grep(pattern = "^mt-", x = rownames(x = sobj), value = TRUE, ignore.case = T)
  sobj[["percent.mito"]] <- Seurat::PercentageFeatureSet(sobj, features = mito.features)
  
  
  #normalize and cluster
  suppressWarnings(sobj <- Seurat::SCTransform(sobj, verbose = T))
  
  sobj <- Seurat::RunPCA(object = sobj, verbose = F)
  
  sobj <- Seurat::FindNeighbors(object = sobj, dims = 1:20, verbose = F)
  sobj <- Seurat::FindClusters(object = sobj, resolution = 0.1, verbose = F, algorithm = 4)
  
  sobj <- RunUMAP(sobj, dims = 1:20)
  
  
  ### run auto filter ###
  reportlist <- FerrenaSCRNAseq::automatedfiltering(sobj, clusters = 'SCT_snn_res.0.1',
                                                    iterativefilter.mito = F)
  
  #add autofilter results to metadata
  autofilterres <- reportlist[[1]]
  sobj$filteredout <- autofilterres$filteredout
  sobj$filterreason <- autofilterres$filterreason
  
  
  
  #plot autofilter results
  pdf(  paste0(outputdir, '/qc/autofilter.pdf'), 7,7)
  print( DimPlot(sobj, label = T) )
  
  
  print( FeaturePlot(sobj, c('nCount_RNA', 'nFeature_RNA', 'percent.mito'), order = T) + 
           DimPlot(sobj, group.by = 'filteredout')
  )
  
  print(reportlist)
  dev.off()
  
  #save autofilter output
  saveRDS( reportlist, paste0(outputdir, '/qc/reportlist-autofilter.rds') )
  
  #filter
  goodcells <- autofilterres[autofilterres$filteredout == 'No', 'barcodes']
  sobj <- sobj[,goodcells]
  
  #reprocess
  #normalize and cluster
  suppressWarnings(sobj <- Seurat::SCTransform(sobj, verbose = T))
  
  sobj <- Seurat::RunPCA(object = sobj, verbose = F)
  
  sobj <- Seurat::FindNeighbors(object = sobj, dims = 1:20, verbose = F)
  sobj <- Seurat::FindClusters(object = sobj, resolution = 0.1, verbose = F, algorithm = 4)
  
  sobj <- RunUMAP(sobj, dims = 1:20)
  
  # DimPlot(sobj, label = T)
  
  
  #use doublet filtering
  dfdf <- FerrenaSCRNAseq::doubletfinderwrapper(sobj, clusters = 'SCT_snn_res.0.1')
  
  #add doublet info to seurat
  sobj$DoubletFinderClassification <- dfdf$DoubletFinderClassification
  
  pdf(  paste0(outputdir, '/qc/doublets.pdf'), 7,7)
  
  print( Seurat::DimPlot(sobj, group.by = 'DoubletFinderClassification') )
  
  
  dev.off()
  
  
  #remove all doublets
  md <- sobj@meta.data
  md <- md[md$DoubletFinderClassification == 'Singlet' ,]
  goodcells <- rownames(md)
  
  sobj <- sobj[,goodcells]
  
  
  #reprocess without doublets
  suppressWarnings(sobj <- Seurat::SCTransform(sobj, verbose = T))
  
  sobj <- Seurat::RunPCA(object = sobj, verbose = F)
  
  sobj <- Seurat::FindNeighbors(object = sobj, dims = 1:20, verbose = F)
  sobj <- Seurat::FindClusters(object = sobj, resolution = 0.1, verbose = F, algorithm = 4)
  
  sobj <- RunUMAP(sobj, dims = 1:20)
  
  # make sure to add in cell cycle, etc
  
  #save it
  saveRDS(sobj,  paste0(outputdir, '/sobj.rds') )
  
}










message('\n\n\nAll done!\n\n\n')


sessionInfo()

date()
sys_end_time <- Sys.time()
message('Runtime:')
sys_end_time - sys_start_time


rm(list=ls())

q('no')

