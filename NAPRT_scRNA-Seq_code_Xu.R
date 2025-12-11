# first -----
rm(list = ls())
options(stringsAsFactors = F)
gc()
options(future.globals.maxSize = 2000 * 1024^2)
options(bitmapType = 'cairo')

##
library(ggplot2)
library(ggthemes)
library(ggsci)
library(Seurat)
library(VISION)
##
setwd(".../path2data/NaprtKO/scRNA-Seq/")
#Read gene expression matrix from output of cellranger ----
#data dir
data_files = ".../path2data/NaprtKO/scRNA-Seq/"

object_name = "NaprtKO"

samples<-list.files(data_files)

#samples<-samples[-1] 

## read data -----
#read multiple samples data----
#samples = c("WT","KO","WT_DSS","KO_DSS")


{
  sceList = lapply(samples, function(pro){
    folder = paste0(data_files,pro,"/out/")
    obj=CreateSeuratObject(counts = Read10X(folder,gene.column = 1),
                           project = pro)
    return(obj)
  })
  
  seuratdata<-merge(sceList[[1]],sceList[-1])
  
  dir.create("output")
  dir.create("Rdata")
  
  saveRDS(seuratdata, file = paste0("./Rdata/1.0_read10X_", object_name, ".Rds"))
}

seuratdata$Sample <- seuratdata@active.ident
Idents(seuratdata)<-"Sample"
Idents(seuratdata)<-factor(Idents(seuratdata),levels = samples)
seuratdata$Sample <- seuratdata@active.ident

saveRDS(seuratdata, file = paste0("./Rdata/1.1_read10X_", object_name, "_sample.Rds"))

seuratdata<-readRDS("./")
names(seuratdata)
Idents(seuratdata) = "Sample"
levels(seuratdata)

#
{
  plotlength = 1
  
  MT_name <- "^mt-"
  Rib_name <- "^Rp[Sl]"
  Hgb_name <- "^Hb[Ab]"
  number_pca = 50
  number_pca_choose = 1:20
  UMAP_min_dist = 0.4
  reduction_number = c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 0.9,1, 1.2, 1.4, 1.6)
  harmony_nclust = 15
  harmony_max_iter = 15
  cluster_max_iter = 20
  number_harmony_choose = 1:15
  sample_names = c("WT","KO","WT_DSS","KO_DSS")
}


### make QC plot
Idents(seuratdata) <-  "Sample"
{
  seuratdata <- PercentageFeatureSet(seuratdata, pattern = MT_name, col.name = "percent.mt")
  seuratdata <- PercentageFeatureSet(seuratdata, pattern = Rib_name, col.name = "percent.rib")
  seuratdata <- PercentageFeatureSet(seuratdata, pattern = Hgb_name, col.name = "percent.hgb")
  p1 <- VlnPlot(seuratdata, features = "nFeature_RNA", pt.size = 0.01) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_nFeature_RNA.png", width = plotlength*length(levels(seuratdata)), height = 6)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_nFeature_RNA.pdf", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "nFeature_RNA", pt.size = 0) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_nFeature_RNA_nopoint.png", width = plotlength*length(levels(seuratdata)), height = 6)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_nFeature_RNA_nopoint.pdf", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "nCount_RNA", pt.size = 0.01) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_nCount_RNA.png", width = plotlength*length(levels(seuratdata)), height = 6)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_nCount_RNA.pdf", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "nCount_RNA", pt.size = 0) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_nCount_RNA_nopoint.png", width = plotlength*length(levels(seuratdata)), height = 6)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_nCount_RNA_nopoint.pdf", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.mt", pt.size = 0.01) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_percent.mt.png", width = plotlength*length(levels(seuratdata)), height = 6)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_percent.mt.pdf", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.mt", pt.size = 0) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_percent.mt_nopoint.png", width = plotlength*length(levels(seuratdata)), height = 6)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_percent.mt_nopoint.pdf", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.mt", y.max = 30, pt.size = 0.01) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_percent.mt2.png", width = plotlength*length(levels(seuratdata)), height = 6)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_percent.mt2.pdf", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.mt", y.max = 30, pt.size = 0) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_percent.mt2_nopoint.png", width = plotlength*length(levels(seuratdata)), height = 6)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_percent.mt2_nopoint.pdf", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.rib", pt.size = 0.01) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_percent.rib.png", width = plotlength*length(levels(seuratdata)), height = 6)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_percent.rib.pdf", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.hgb", pt.size = 0.01) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_percent.hgb.png", width = plotlength*length(levels(seuratdata)), height = 6)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_percent.hgb.pdf", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.rib", pt.size = 0) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_percent.rib_nopoint.png", width = plotlength*length(levels(seuratdata)), height = 6)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_percent.rib_nopoint.pdf", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.hgb", pt.size = 0) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_percent.hgb_nopoint.png", width = plotlength*length(levels(seuratdata)), height = 6)
  ggsave(p1, file = "./output/1.1_raw_data_VlnPlot_percent.hgb_nopoint.pdf", width = plotlength*length(levels(seuratdata)), height = 6)
  
  summ_data <- cell_summ(seuratdata, type = "Sample", is.Spatial = FALSE)
  
  write.table(summ_data, file = "./output/1.1_raw_data_summ.xls", row.names = F, sep = "\t")
  
  saveRDS(seuratdata, file = paste0("./Rdata/1.1_", object_name, "_raw_data.Rds"))
  
  print("reading data finished.")
}
## DecontX ----- get rid of ambient RNA
library(celda)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(MatrixGenerics)
library(matrixStats)

{
  options(mc.cores = 15) #use multiple cores
  seuratdata <- DecontX_seuratdata(seuratdata = seuratdata, convergence = 0.0005,
                                   min.cells = 3, min.features = 3)  # remove ambient RNA
  
  options(mc.cores = 1)
  Idents(seuratdata) <- factor(Idents(seuratdata), levels = sample_names)
  seuratdata$Sample <- seuratdata@active.ident
  
  seuratdata <- PercentageFeatureSet(seuratdata, pattern = MT_name, col.name = "percent.mt")
  seuratdata <- PercentageFeatureSet(seuratdata, pattern = Rib_name, col.name = "percent.rib")
  seuratdata <- PercentageFeatureSet(seuratdata, pattern = Hgb_name, col.name = "percent.hgb")
  
  p1 <- VlnPlot(seuratdata, features = "nFeature_RNA", pt.size = 0.01) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.2_DecontX_VlnPlot_nFeature_RNA.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "nFeature_RNA", pt.size = 0) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.2_DecontX_VlnPlot_nFeature_RNA_nopoint.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "nCount_RNA", pt.size = 0.01) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.2_DecontX_VlnPlot_nCount_RNA.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "nCount_RNA", pt.size = 0) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.2_DecontX_VlnPlot_nCount_RNA_nopoint.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.mt", pt.size = 0.01) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.2_DecontX_VlnPlot_percent.mt.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.mt", pt.size = 0) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.2_DecontX_VlnPlot_percent.mt_nopoint.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.mt", y.max = 40, pt.size = 0.01) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.2_DecontX_VlnPlot_percent.mt2.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.mt", y.max = 40, pt.size = 0) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.2_DecontX_VlnPlot_percent.mt2_nopoint.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.rib", pt.size = 0.01) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.2_DecontX_VlnPlot_percent.rib.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.rib", pt.size = 0) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.2_DecontX_VlnPlot_percent.rib_nopoint.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.rib", y.max = 30, pt.size = 0.01) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.2_DecontX_VlnPlot_percent.rib2.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.rib", y.max = 30, pt.size = 0) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.2_DecontX_VlnPlot_percent.rib2_nopoint.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.hgb", pt.size = 0.01) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.2_DecontX_VlnPlot_percent.hgb.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.hgb", pt.size = 0) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.2_DecontX_VlnPlot_percent.hgb_nopoint.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  summ_data <- cell_summ(seuratdata, type = "Sample", is.Spatial = FALSE)
  
  write.table(summ_data, file = "./output/1.2_DecontX_data_summ.xls", row.names = F, sep = "\t")
  
  saveRDS(seuratdata, file = paste0("./Rdata/1.2_", object_name, "_DecontX_data.Rds"))
  
  message("DecontX finished.")
}
## filter data -----
library(dplyr)

{
  Idents(seuratdata) <- "Sample"
  seuratdata <- subset(seuratdata, subset = percent.rib < 60) #remove cells with ribosome genes expression % >40 
  seuratdata <- subset(seuratdata, subset = percent.hgb < 5)  #remove cells with hemoglobin genes expression % >5
  
  Idents(seuratdata) <- "Sample"
  seuratdata <- filter_seuratdata(seuratdata = seuratdata, filter_percent.mt = 10,
                                  sd_number = 2)  #remove cells with mitochondrial genes expression % >40
{

  seuratdata <- PercentageFeatureSet(seuratdata, pattern = MT_name, col.name = "percent.mt")
  seuratdata <- PercentageFeatureSet(seuratdata, pattern = Rib_name, col.name = "percent.rib")
  seuratdata <- PercentageFeatureSet(seuratdata, pattern = Hgb_name, col.name = "percent.hgb")
  # seuratdata <- subset(seuratdata, subset = (nCount_RNA > 200))
  # seuratdata <- subset(seuratdata, subset = (nFeature_RNA > 200))
  
  Idents(seuratdata) <- "Sample"
  
  seuratdata$Sample <- seuratdata@active.ident
  
  p1 <- VlnPlot(seuratdata, features = "nFeature_RNA", pt.size = 0.01) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.4_filter_data_VlnPlot_nFeature_RNA.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "nFeature_RNA", pt.size = 0) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.4_filter_data_VlnPlot_nFeature_RNA_nopoint.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "nCount_RNA", pt.size = 0.01) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.4_filter_data_VlnPlot_nCount_RNA.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "nCount_RNA", pt.size = 0) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.4_filter_data_VlnPlot_nCount_RNA_nopoint.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.mt", pt.size = 0.01) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.4_filter_data_VlnPlot_percent.mt.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.mt", pt.size = 0) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.4_filter_data_VlnPlot_percent.mt_nopoint.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.rib", pt.size = 0.01) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.4_filter_data_VlnPlot_percent.rib.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.rib", pt.size = 0) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.4_filter_data_VlnPlot_percent.rib_nopoint.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.hgb", pt.size = 0.01) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.4_filter_data_VlnPlot_percent.hgb.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  p1 <- VlnPlot(seuratdata, features = "percent.hgb", pt.size = 0) & NoLegend() & labs(x = NULL)
  ggsave(p1, file = "./output/1.4_filter_data_VlnPlot_percent.hgb_nopoint.png", width = plotlength*length(levels(seuratdata)), height = 6)
  
  summ_data <- cell_summ(seuratdata, type = "Sample", is.Spatial = FALSE)
  
  write.table(summ_data, file="./output/1.4_filter_data_summ.xls", row.names = F, sep = "\t")
  
  message("filtering data finished.")
}



saveRDS(seuratdata, file = paste0("./Rdata/1.3_", object_name, "_clean_data.Rds"))


## CellCycle -----
{
  seuratdata <- NormalizeData(seuratdata, verbose = F)
  
  seuratdata <- cellcycle_seuratdata(seuratdata = seuratdata)
  
  message("cell cycle running finished.")
}

## SCTransform -----
{
  options(mc.cores = 10)
  seuratdata <- SCTransform(seuratdata, do.scale = F, variable.features.n = 3000, return.only.var.genes = TRUE,
                            vars.to.regress = c("nCount_RNA", "CC.Difference"), conserve.memory = F)
  
  options(mc.cores = 1)
  
  var.features <- VariableFeatures(seuratdata)
  var.features <- var.features[stringr::str_detect(var.features, pattern = "^mt-", negate = T)]
  var.features <- var.features[stringr::str_detect(var.features, pattern = "^Rp[sl]", negate = T)]
  var.features <- var.features[stringr::str_detect(var.features, pattern = "^Hb[ab]", negate = T)]
  var.features <- var.features[stringr::str_detect(var.features, pattern = "^Ig[hkl]", negate = T)]
  var.features <- var.features[stringr::str_detect(var.features, pattern = "Ac233755.1", negate = T)]
  var.features <- var.features[stringr::str_detect(var.features, pattern = "^Tr[abdg]", negate = T)]
  var.features <- var.features[stringr::str_detect(var.features, pattern = "Stmn1", negate = T)]
  var.features <- var.features[stringr::str_detect(var.features, pattern = "Mki67", negate = T)]
  var.features <- var.features[stringr::str_detect(var.features, pattern = "^Ifi", negate = T)]
  var.features <- var.features[stringr::str_detect(var.features, pattern = "^Irf", negate = T)]
  VariableFeatures(seuratdata) <- var.features
  
  message("SCTransform running finished.")
}

## PCA -----
{
  seuratdata <- RunPCA(seuratdata, npcs = number_pca, verbose = F)
  
  message("PCA running finished.")
}

## no correction -----
{
  seuratdata <- FindNeighbors(seuratdata, reduction = "pca", dims = number_pca_choose)
  
  seuratdata <- RunUMAP(seuratdata, n.components = 2, reduction = "pca", n.neighbors = 40,
                        min.dist = UMAP_min_dist, dims = number_pca_choose, seed.use = 3,
                        reduction.name = "umapNC", reduction.key = "umapNC_")
  
  seuratdata <- RunTSNE(seuratdata, dims = number_pca_choose, reduction = "pca",
                        reduction.name = "tsneNC", reduction.key = "tsneNC_")
  
  p1 <- DimPlot(seuratdata, reduction = "umapNC", label = TRUE, group.by = "Sample", pt.size = 0.25, raster = FALSE)
  ggsave(p1, file = "./output/2.1_umap_Sample_NC.png", width = 6.5, height = 6)
  ggsave(p1, file = "./output/2.1_umap_Sample_NC.pdf", width = 6.5, height = 6)
  
  p1 <- DimPlot(seuratdata, reduction = "umapNC", label = FALSE, group.by = "Sample", pt.size = 0.25, raster = FALSE)
  ggsave(p1, file = "./output/2.1_umap_Sample_NC_nolabel.png", width = 6.5, height = 6)
  ggsave(p1, file = "./output/2.1_umap_Sample_NC_nolabel.pdf", width = 6.5, height = 6)
  
  p1 <- DimPlot(seuratdata, reduction = "tsneNC", label = TRUE, group.by = "Sample", pt.size = 0.25, raster = FALSE)
  ggsave(p1, file = "./output/2.2_tsne_Sample_NC.png", width = 6.5, height = 6)
  ggsave(p1, file = "./output/2.2_tsne_Sample_NC.pdf", width = 6.5, height = 6)
  
  p1 <- DimPlot(seuratdata, reduction = "tsneNC", label = FALSE, group.by = "Sample", pt.size = 0.25, raster = FALSE)
  ggsave(p1, file = "./output/2.2_tsne_Sample_NC_nolabel.png", width = 6.5, height = 6)
  ggsave(p1, file = "./output/2.2_tsne_Sample_NC_nolabel.pdf", width = 6.5, height = 6)
  
  for(reduction_n in reduction_number) {
    seuratdata <- FindClusters(seuratdata, resolution = reduction_n)
    
    seuratdata@meta.data[[paste0("SCT_snn_res.", reduction_n, "NC")]] <- seuratdata@active.ident
    
    p1 <- DimPlot(seuratdata, reduction = "umapNC", label = TRUE,
                  group.by = paste0("SCT_snn_res.", reduction_n, "NC"), pt.size = 0.25, raster = FALSE)
    ggsave(p1, file = paste0("./output/2.1_umapNC_", reduction_n, "_NC.png"), width = 7, height = 6)
    ggsave(p1, file = paste0("./output/2.1_umapNC_", reduction_n, "_NC.pdf"), width = 7, height = 6)
    
    p1 <- DimPlot(seuratdata, reduction = "umapNC", label = FALSE,
                  group.by = paste0("SCT_snn_res.", reduction_n, "NC"), pt.size = 0.25, raster = FALSE)
    ggsave(p1, file = paste0("./output/2.1_umapNC_", reduction_n, "_NC_nolabel.png"), width = 7, height = 6)
    ggsave(p1, file = paste0("./output/2.1_umapNC_", reduction_n, "_NC_nolabel.pdf"), width = 7, height = 6)
    
    p1 <- DimPlot(seuratdata, reduction = "tsneNC", label = TRUE,
                  group.by = paste0("SCT_snn_res.", reduction_n, "NC"), pt.size = 0.25, raster = FALSE)
    ggsave(p1, file = paste0("./output/2.2_tsneNC_", reduction_n, "_NC.png"), width = 7, height = 6)
    ggsave(p1, file = paste0("./output/2.2_tsneNC_", reduction_n, "_NC.pdf"), width = 7, height = 6)
    
    p1 <- DimPlot(seuratdata, reduction = "tsneNC", label = FALSE,
                  group.by = paste0("SCT_snn_res.", reduction_n, "NC"), pt.size = 0.25, raster = FALSE)
    ggsave(p1, file = paste0("./output/2.2_tsneNC_", reduction_n, "_NC_nolabel.png"), width = 7, height = 6)
    ggsave(p1, file = paste0("./output/2.2_tsneNC_", reduction_n, "_NC_nolabel.pdf"), width = 7, height = 6)
  }
  
  saveRDS(seuratdata, file = paste0("./Rdata/1.4_", object_name, "_NC.Rds"))
  
  message("UMAP with no correction running finished.")
}
## nCount nFeature NC-----
{
  plotdata <- seuratdata@meta.data
  plotdata <- plotdata[, c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rib", "Phase")]
  plotdata$cell <- rownames(plotdata)
  
  temp_dimplot <- as.data.frame(Embeddings(seuratdata, reduction = "umapNC"))
  temp_dimplot$cell <- rownames(temp_dimplot)
  
  plotdata <- merge(x = plotdata, y = temp_dimplot, by = "cell")
  rownames(plotdata) <- plotdata$cell
  plotdata$cell <- NULL
  
  for (i in c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rib")) {
    tmp <- plotdata[, c(i, "umapNC_1", "umapNC_2")]
    tmp <- tmp[order(tmp[,1]),]
    
    p1 <- ggplot(tmp, aes(umapNC_1, umapNC_2, color = tmp[,1])) + geom_point(size = 0.15) +
      scale_colour_gradientn(colors = c("#3d0e4a", "#32668d", "#3fbca6", "#FFFF00")) +
      theme_classic() + theme(legend.title = element_blank()) + ggtitle(i) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(p1, file = paste0("./output/3.4_umapNC_", i, ".png"), width = 6.5, height = 6)
    ggsave(p1, file = paste0("./output/3.4_umapNC_", i, ".pdf"), width = 6.5, height = 6)
  }
  
  p1 <- DimPlot(seuratdata, reduction = "umapNC", group.by = "Phase", pt.size = 0.25, raster = FALSE)
  ggsave(p1, file = "./output/3.4_umapNC_Phase.png", width = 6.5, height = 6)
  ggsave(p1, file = "./output/3.4_umapNC_Phase.pdf", width = 6.5, height = 6)
}



ALLMARKER<-c()

Idents(seuratdata) = "SCT_snn_res.0.3NC"
for (j in 1:length(markerdata)) {
  
  cellname<-names(markerdata[j])
  marker_gene <- unique(markerdata[[j]])
  
  ALLMARKER<-c(ALLMARKER,marker_gene)
  
  FeaturePlot(seuratdata, features = marker_gene, reduction = "umapNC", raster = FALSE, pt.size = 0.25, order = T,
              cols = c("lightgrey", "red"),max.cutoff = 3)
  
  for(i in marker_gene) {
    p1 <- FeaturePlot(seuratdata, features = i, reduction = "umapNC", raster = FALSE, pt.size = 0.25, order = F,
                      cols = c("lightgrey", "red"),max.cutoff = 3)
    p1
    ggsave(p1, file = paste0("./output/4.1_markerNC_featureplot_", cellname,"_",i, ".png"), width = 6, height = 6)
  }
  
  p1<-DotPlot(seuratdata, features = marker_gene, group.by = "SCT_snn_res.0.3NC",
              cluster.idents = T, assay = "RNA", col.min = -1,col.max = 2) + 
    scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33','#FFCC33'))+
    RotatedAxis() +ggtitle(paste0(cellname)) +theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p1
  ggsave(p1, file = paste0("./output/4.1_markerNC_featureplot_",cellname,"_dotplot.png"), width = 6, height = 6)
  
  print(j)
}

ALLMARKER<-unique(ALLMARKER)
ALLMARKER<-factor(ALLMARKER,levels = ALLMARKER)




