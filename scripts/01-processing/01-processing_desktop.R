#time things:
sys_start_time <- Sys.time()



# this is processing for samples on threadripper after cellranger




library(tidyverse)
#library(RISC)
library(Seurat)
library(FerrenaSCRNAseq)
library(DoubletFinder)
library(enrichR)
library(ggVennDiagram)




setwd('/home/Alex/data/shanghai_clinical_scRNAseq/results')

#close any existing plotting devices
while (!is.null(dev.list()))  dev.off()







set.seed(2020)



#list input and output paths. 

indir <- '/home/Alex/data/shanghai_clinical_scRNAseq/data/SRA-reprocessed/'

procdir <- 'processed/'

dir.create(procdir)


samps <- list.files(indir, pattern='BC')
samps <- str_sort(samps, numeric = T)
samppaths <- paste0(indir, '/', samps)

sampdex=1

for(sampdex in c(1:length(samps)) ){
  
  #get samples
  samp <- samps[sampdex]
  samppath <- samppaths[sampdex]
  
  #set output dirs
  outputdir <- paste0(procdir, samp)
  dir.create(outputdir)
  dir.create( paste0(outputdir, '/qc') )
  dir.create( paste0(outputdir, '/celltypeid') )
  
  message('\nReading in ', samp, '\n')
  
  
  
  #read in all data
  # cellranger h5 file
  # souporcell output
  h5file <- paste0(samppath, '/', 'filtered_feature_bc_matrix.h5')
  soupfile <- paste0(samppath, '/', 'clusters.tsv')
  
  #tryCatch({
  
  
  sobj <- CreateSeuratObject(   Seurat::Read10X_h5(h5file), min.cells= 3, project = samp)
  
  #mito content, add to metadata
  mito.features <- grep(pattern = "^mt-", x = rownames(x = sobj), value = TRUE, ignore.case = T)
  sobj[["percent.mito"]] <- Seurat::PercentageFeatureSet(sobj, features = mito.features)
  
  
  #normalize and cluster
  suppressWarnings(sobj <- Seurat::SCTransform(sobj, verbose = T))
  
  sobj <- Seurat::RunPCA(object = sobj, verbose = F)
  
  sobj <- Seurat::FindNeighbors(object = sobj, dims = 1:20, verbose = F)
  sobj <- Seurat::FindClusters(object = sobj, resolution = 0.1, verbose = F, algorithm = 1)
  
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
  sobj <- Seurat::FindClusters(object = sobj, resolution = 0.1, verbose = F, algorithm = 1)
  
  sobj <- RunUMAP(sobj, dims = 1:20)
  
  # DimPlot(sobj, label = T)
  
  
  #use doublet filtering
  dfdf <- FerrenaSCRNAseq::doubletfinderwrapper(sobj, clusters = 'SCT_snn_res.0.1')
  
  #read in souporcell and select only cells in seurat + match order
  soup <- read.table(soupfile, header = T)
  soup <- soup[match(colnames(sobj), soup$barcode),]
  
  #add doublet info to seurat
  sobj$DoubletFinderClassification <- dfdf$DoubletFinderClassification
  sobj$SouporcellDoubletClassification <- stringr::str_to_title( soup$status )
  
  #add souporcell clusters
  sobj$SouporcellClusters <- soup$assignment
  
  pdf(  paste0(outputdir, '/qc/doublets.pdf'), 7,7)
  
  md <- sobj@meta.data
  dfdoub <- rownames( md[md$DoubletFinderClassification == 'Doublet', ] )
  soupdoub <- rownames( md[md$SouporcellDoubletClassification == 'Doublet',] )
  soupunc <- rownames( md[md$SouporcellDoubletClassification == 'Unassigned',] )
  
  v <- list("DoubletFinder\nDoublet" = dfdoub,
            "Souporcell\nDoublet" = soupdoub,
            "Souporcell\nUnclassified" = soupunc)
  
  print( ggVennDiagram::ggVennDiagram(v) )
  
  print( DimPlot(sobj, group.by = 'DoubletFinderClassification') )
  print( DimPlot(sobj, group.by = 'SouporcellDoubletClassification') )
  
  
  dev.off()
  
  
  #remove all doublets
  md <- sobj@meta.data
  md <- md[md$DoubletFinderClassification == 'Singlet' ,]
  md <- md[md$SouporcellDoubletClassification == 'Singlet' ,]
  goodcells <- rownames(md)
  
  sobj <- sobj[,goodcells]
  
  
  #reprocess without doublets
  suppressWarnings(sobj <- Seurat::SCTransform(sobj, verbose = T))
  
  sobj <- Seurat::RunPCA(object = sobj, verbose = F)
  
  sobj <- Seurat::FindNeighbors(object = sobj, dims = 1:20, verbose = F)
  sobj <- Seurat::FindClusters(object = sobj, resolution = 0.1, verbose = F, algorithm = 1)
  
  sobj <- RunUMAP(sobj, dims = 1:20)
  
  #for all other purposes, we should use RNA assay...
  DefaultAssay(sobj) <- 'RNA'
  sobj <- NormalizeData(sobj)
  sobj <- ScaleData(sobj)
  
  
  
  ### try to call cell type ###
  
  
  #get DEGs
  m <- FindAllMarkers(sobj, test.use = 'LR')
  
  n <- 5
  top <- m %>% group_by(cluster) %>% top_n(n = n, wt = avg_log2FC) #this code is changed, in new seurat versions, need "log2FC"...
  
  
  message('\n\n\nRunning EnrichR...')
  
  e_hu_list <- list()
  e_mo_list <- list()
  for(clust in unique(m$cluster)){
    
    message('\nenrichR for cluster ', clust)
    dbs <- c("Human_Gene_Atlas", "Mouse_Gene_Atlas")
    
    enriched <- enrichr(m[m$cluster==clust, 'gene'], dbs)
    e_hu <- enriched[[1]][1:5, c(1,2,3,4,7,8)]
    e_mo <- enriched[[2]][1:5, c(1,2,3,4,7,8)]
    
    e_hu$clust <- clust
    e_mo$clust <- clust
    
    e_hu_list[[clust]] <- e_hu
    e_mo_list[[clust]] <- e_mo
    
  }
  
  
  
  e_hu <- dplyr::bind_rows(e_hu_list)
  e_mo <- dplyr::bind_rows(e_mo_list)
  
  e_hu$Term <- factor(e_hu$Term, levels = rev(unique(e_hu$Term)))
  e_hu$clust <- factor(e_hu$clust, levels = str_sort(unique(e_hu$clust), numeric = T))
  
  
  
  
  pdf(  paste0(outputdir, '/celltypeid/celltype.pdf'), 7,7)
  
  #make a ncie table
  
  avgmd <- data.frame(cluster = c(names(table(sobj$SCT_snn_res.0.1)), 'TOTAL'),
                      cellnum = c(table(sobj$SCT_snn_res.0.1), ncol(sobj)),
                      medianUMI = c(aggregate(nCount_RNA ~ SCT_snn_res.0.1, sobj@meta.data, median)[,2], median(sobj$nCount_RNA)),
                      medianNumGene = c(aggregate(nFeature_RNA ~ SCT_snn_res.0.1, sobj@meta.data, median)[,2], median(sobj$nFeature_RNA)),
                      medianPercentMito = paste0( formatC( c(aggregate(percent.mito  ~ SCT_snn_res.0.1, sobj@meta.data, median)[,2], median(sobj$percent.mito )), 2), '%')
                      
  )
  
  print(plot(gridExtra::tableGrob(avgmd)))
  
  
  
  print( DimPlot(sobj, label = T) + ggtitle('Louvain, Res = 0.1') )
  
  #repelplot
  mapping <- data.frame(cluster = names(table(sobj$SCT_snn_res.0.1)),
                        nCount_RNA = log(aggregate(nCount_RNA ~ SCT_snn_res.0.1, sobj@meta.data, median)[,2]),
                        nFeature_RNA = log(aggregate(nFeature_RNA ~ SCT_snn_res.0.1, sobj@meta.data, median)[,2]))
  
  g1=ggplot(sobj@meta.data, aes(log(nCount_RNA), log(nFeature_RNA), col = SCT_snn_res.0.1))+geom_point()+
    ggrepel::geom_label_repel(inherit.aes = F, data = mapping, aes(nCount_RNA, nFeature_RNA, label = cluster, fill = cluster))
  
  
  FeaturePlot(sobj, c('nCount_RNA', 'nFeature_RNA', 'percent.mito'), order = T, label = T) +g1
  
  print( DimPlot(sobj, label = T, group.by = 'SouporcellClusters') )
  
  
  genepanel <- c('Ptprc', 'Tnfrsf11a', 'Cd68',
                 'Cd3e', 'Cd19', 'Fap',
                 'Sp7', 'Sox9', 'Nes')
  
  #if human...
  genepanel <- toupper(genepanel)
  
  print( FeaturePlot(sobj, genepanel, label = T, order = T) )
  
  print(DoHeatmap(sobj, top$gene, raster = F) + NoLegend() )
  
  
  
  
  # Plot EnrichR
  
  print(
    ggplot(e_hu, aes(x=clust, y = Term, size = Odds.Ratio, col=-log10(P.value) ))+
      geom_point()+
      ggtitle("EnrichR: Human_Gene_Atlas")
  )
  
  e_mo$Term <- factor(e_mo$Term, levels = rev(unique(e_mo$Term)))
  e_mo$clust <- factor(e_mo$clust, levels = str_sort(unique(e_mo$clust), numeric = T))
  
  print(
    ggplot(e_mo, aes(x=clust, y = Term, size = Odds.Ratio, col=-log10(P.value) ))+
      geom_point()+
      ggtitle("EnrichR: Mouse_Gene_Atlas")
  )
  
  
  dev.off()
  
  
  #save sobj and degs
  
  saveRDS(sobj,  paste0(outputdir, '/sobj.rds') )
  write.csv(m, paste0(outputdir, '/markers_lr.csv'), quote = F, row.names = F )
  
  
  
  
  
}









message('\n\n\nAll done!\n\n\n')


sessionInfo()

date()
sys_end_time <- Sys.time()
message('Runtime:')
sys_end_time - sys_start_time


rm(list=ls())

q('no')

