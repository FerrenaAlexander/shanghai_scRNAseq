#time things:
sys_start_time <- Sys.time()



# this is to plot the key markers from the paper
# I also cluster with higher resolution in this script.




library(tidyverse)
#library(RISC)
library(Seurat)
library(FerrenaSCRNAseq)



set.seed(2020)



setwd('/home/Alex/data/shanghai_clinical_scRNAseq/results')

#close any existing plotting devices
while (!is.null(dev.list()))  dev.off()



#define cluster markers, sub-cluster markers from paper

# these are defined in paper supplementary table 2

clustermarkers <- list(Osteoblastic = c('RUNX2', 'COL1A1', 'CDH11', 'IBSP'),
                       Chondroblastic = c('SOX9', 'ACAN', 'PTH1R'),
                       Osteoclast = c('ACP5', 'CTSK', 'MMP9'),
                       Myeloid = c('CD74', 'CD14', 'FCGR3A'),
                       T_NK = c('CD3D', 'IL7R', 'CD8A', 'CD4', 'NKG7', 'GNLY'),
                       # NK_NKT = c('NKG7', 'GNLY', 'CD3D'), #NK is CD3-, NKT is cd3+
                       DC = c('CD1C', 'FCER1A', 'CLEC9A', 'CCR7', 'CD14', 'CD163'),
                       Fibroblast = c('DCN', 'COL1A1'),
                       Pericyte = c('RGS5', 'ACTA2'),
                       MSC = c('MME', 'THY1', 'CXCL12', 'SFRP2'),
                       Endothelial = c('PECAM1', 'VWF'),
                       Myoblast = c('MYL1', 'MYLPF'),
                       Bcell = c('MS4A1', 'CD19', 'JCHAIN'),
                       Proliferative = c('MKI67', 'TOP2A', 'PCNA')
)

#sub-cluster markers
# we will need the files...
supdatfolder <- 'paper_supplement_subclustermarkers'


subclustermarkers <- list()

supdatfiles <- list.files(supdatfolder, full.names = T)
for(fp in supdatfiles){
  
  #get celltype from file name
  ct <- tools::file_path_sans_ext(basename(fp))
  ct <- str_split_fixed(ct, '\\_', 2)[1,2]
  
  #read file
  subm <- readxl::read_excel(fp, skip = 1)
  
  #keep only pdj<0.05 genes
  subm <- subm[subm$p_val_adj<0.05,]
  
  #make list...
  subctlist <- list()
  
  subclusts <- unique(subm$subcluster)
  for(subct in subclusts){
    
    subctlist[[subct]] <- subm[subm$subcluster==subct,]$gene
  }
  
  #add to whole list
  subclustermarkers[[ct]] <- subctlist
}
rm(subctlist, subm, subclusts, subct, supdatfiles, supdatfolder, ct, fp)





#list input and output paths. 

procdir <- 'desktop_minimal_clone/results/processed/'
procdir <- 'processed/'

samps <- list.files(procdir)

sampdex=1

for(sampdex in c(1:length(samps)) ){
  
  #get sample
  samp <- samps[sampdex]
  
  #get path
  fp <- paste0(procdir, '/', samp, '/sobj.rds')
  
  
  
  #read in all data
  sobj <- readRDS(fp)
  
  
  
  
  #clustering and deeper clustering
  d_01 <-  DimPlot(sobj, label = T, group.by = 'SCT_snn_res.0.1') + labs(title = sobj@project.name, 
                                           subtitle= 'Louvain, Res = 0.1')
  
  
  DefaultAssay(sobj) <- 'SCT'
  sobj <- FindClusters(sobj, resolution = 0.8)
  
  d_08 <- DimPlot(sobj, label = T, group.by = 'SCT_snn_res.0.8') + labs(title = sobj@project.name, 
                                          subtitle= 'Louvain, Res = 0.8')
  
  
  #save sobj after this just in case.
  saveRDS(sobj, fp)
  
  
  #sourporcell
  d_soup <- DimPlot(sobj, group.by = 'SouporcellClusters')
  
  ### plot markers ###
  DefaultAssay(sobj) <- 'RNA'
  
  
  #plot global markers...
  genepanel <- c('Ptprc', 'Nfatc1', 'Aif1',
                 'Cd3e', 'Cd19', 'Fap',
                 'Sp7', 'Sox9', 'Pecam1')
  
  #if human...
  genepanel <- toupper(genepanel)

  
  f_gp <- FeaturePlot(sobj, genepanel, label = T, order = T, repel = F)
  
  #plot cell type markers...
  celltypes <- names(clustermarkers)
  celltypeplots <- list()
  for(ct in celltypes){
    genes <- clustermarkers[[ct]]
    
    celltypeplots[[ct]] <- FeaturePlot(sobj, genes, label = T, order = T) + patchwork::plot_annotation(title = ct)
    
  }
  
  
  pdfname <- paste0(procdir, '/', samp, '/celltypeid/', 'paper_clustermarkers.pdf')
  pdf(pdfname, height = 9, width = 9)
  
  print(d_01)
  print(d_08)
  print(d_soup)
  
  print(f_gp)
  
  print(celltypeplots)
  
  dev.off()
  
  
}









message('\n\n\nAll done!\n\n\n')


sessionInfo()

date()
sys_end_time <- Sys.time()
message('Runtime:')
sys_end_time - sys_start_time


rm(list=ls())

q('no')

