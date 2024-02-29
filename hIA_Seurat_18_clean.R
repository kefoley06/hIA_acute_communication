#Kate Foley, foley@iu.edu
#Indiana University 
#2024 

## 0. Upload Libraries -------------------------------------------------------

library(dplyr)
library(Seurat)
library(SeuratObject)
library(SeuratData)
library(patchwork)
library(RColorBrewer)
library(ggplot2)
library(scCustomize)
library(ComplexHeatmap)
library(ggpattern)
library(beepr)


## 1. Upload data -------------------------------------------------------------
dataset_loc <- "/Users/foley/Documents/Projects/hIA/nfeature_5000/hIA-cellranger/filtered/"
ids<- c("hIA11_F_antiAB_3day",
        "hIA15_M_antiAB_3day",
        "hIA16_F_control_3day",
        "hIA25_F_antiAB_24hr",
        "hIA27_F_control_24hr",
        "hIA28_M_control_3day",
        "hIA37_M_control_24hr",
        "hIA39_M_control_3day",
        "hIA41_M_antiAB_24hr",
        "hIA45_M_antiAB_3day",
        "hIA48_M_antiAB_24hr",
        "hIA49_F_control_24hr")

d10x.data <- sapply(ids, function(i){
  d10x <- Read10X(file.path(dataset_loc,i))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})
beep()

experiment.data <- do.call("cbind", d10x.data)

experiment.aggregate <- CreateSeuratObject(
  experiment.data,
  project = "hIA_2",
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "\\-")

hIA_data<-experiment.aggregate

## 2. QC and Selecting Cells for further analysis -----------------------------
hIA_data[["percent.mt"]] <- PercentageFeatureSet(hIA_data, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(hIA_data, features = c("nFeature_RNA"), ncol = 1, pt.size = 0)

plot1 <- FeatureScatter(hIA_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hIA_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

hIA <- subset(hIA_data, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)

#filtered vs unfiltered
hIA #filtered
hIA_data #original
beep()

## 3. Normalize the data ------------------------------------------------------
hIA <- NormalizeData(hIA, normalization.method = "LogNormalize", scale.factor = 10000)

## 4. Assign variables --------------------------------------------------------
hIA$whole.ident<-hIA$orig.ident
head(hIA$whole.ident)

#assign hIA Study ID
hIA$hIA_ID<-substr(hIA$whole.ident,1,5)
head(hIA$hIA_ID)
#assign Sex
hIA$Sex<-substr(hIA$whole.ident,7,7)
head(hIA$Sex)
#assign Injection
hIA$Injection<-substr(hIA$whole.ident,9,12)
head(hIA$Injection)
tail(hIA$Injection)

#assign Timepoint
hIA$Timepoint<-substr(hIA$whole.ident,16,20)
head(hIA$Timepoint)
tail(hIA$Timepoint)

#assign Timepoint2
hIA$Timepoint2<-substr(hIA$whole.ident,16,20)
head(hIA$Timepoint)
tail(hIA$Timepoint)
hIA$Timepoint2[hIA$Timepoint2 == "3day"]<- 'antiAB3day'
hIA$Timepoint2[hIA$Timepoint2 == "24hr"]<- 'antiAB24hr'
hIA$Timepoint2[hIA$Timepoint2 == "_3day"]<- 'control3day'
hIA$Timepoint2[hIA$Timepoint2 == "_24hr"]<- 'control24hr'

## 5. Identification of highly variable features ------------------------------
hIA <- FindVariableFeatures(hIA, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(hIA), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(hIA)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
beep()

## 6.Scale the data -----------------------------------------------------------
all.genes <- rownames(hIA)
hIA <- ScaleData(hIA, features = all.genes)

## 7. Perform linear dimensional reduction ------------------------------------
hIA <- RunPCA(hIA, features = VariableFeatures(object = hIA))
beep()

# Examine and visualize PCA results a few different ways
print(hIA[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(hIA, dims = 1:2, reduction = "pca")
DimPlot(hIA, reduction = "pca")
DimHeatmap(hIA, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(hIA, dims = 1:11, cells = 500, balanced = TRUE)

#### Determine the 'dimensionality' of the dataset ####
ElbowPlot(hIA)
#we will choose 11 here as thats where it really dips off

## 8. Cluster the Cells -------------------------------------------------------
hIA <- FindNeighbors(hIA, dims = 1:11)

# Determine the clusters for various resolutions                                
#hIA <- FindClusters(object = hIA, resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
hIA <-FindClusters(hIA, resolution = .5)
#how many clusters are formed - Look at cluster IDs of the first 5 cells
head(Idents(hIA), 5)
#resolution of .5 = 22 clusters
#staying with .5 for now 

## 9. Run non-linear dimensional reduction (UMAP/tSNE) ------------------------
hIA <- RunUMAP(hIA, dims = 1:11)

color_pal <-c("#AA0DFEFF", #0
              "#3B00FBFF",
              "#66B0FFFF",
              "#FA0087FF", 
              "#782AB6FF", #4
              "#7ED7D1FF", 
              "#1C8356FF", 
              "#16FF32FF", #7
              "#F7E1A0FF",
              "#B00068FF", 
              "#F6222EFF", #10
              "#DEA0FDFF", 
              "#FE00FAFF", #12
              "#325A9BFF",
              "#D85FF7FF",
              "#F8A19FFF", #15
              "#00868B",
              "#FBE426FF",
              "#FEAF16FF", #18
              "#2ED9FFFF", 
              "#3283FEFF",
              "#FC1CBFFF" , #21
              "#C075A6FF")

set.seed(1)
DimPlot(hIA, reduction = "umap", label = T, seed =1)
DimPlot_scCustom(hIA, reduction = "umap", label.size = 5, label.box = T, colors_use = color_pal, repel = F)
DimPlot_scCustom(hIA, reduction = "umap", label = T, label.box = T, split.by = "Timepoint2", ncol = 2, colors_use = color_pal)
DimPlot_scCustom(hIA, reduction = "umap", label = F, label.box = F, split.by = "hIA_ID", ncol = 2, colors_use = color_pal)

#saveRDS(hIA, file = paste0("/Users/kate.foley/Documents/Projects/hAI/SCseq/Rscripts/Savepoint1",version., ".rds"))
#save.image("~/Documents/Projects/hAI/SCseq/Savepoint1.RData")

## 10. Finding DE features (cluster biomarkers)
#find all markers of cluster 1
#cluster0.markers<-FindMarkers(hIA, ident.1 = 0, ident.2 = NULL, min.pct = .25)

# find all markers distinguishing cluster 5 from clusters 0 and 2
clustermicro_7tomicro.markers <- FindMarkers(hIA, ident.1 = c(7), ident.2 = c(0,1,5,6), min.pct = 0.25)
head(clustermicro_7tomicro.markers, n = 20)
write.csv(clustermicro_7tomicro.markers, "clustermicro_7tomicro.markers.csv")

# find markers for every cluster compared to all remaining cells, both positive and negative
hIA.markers <- FindAllMarkers(hIA, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
beep()

hIA.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

clustertopmarkers<-hIA.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

## 11. Gene Expression per cluster ---------------------------------------
#see how particular genes of interest are expressed per cluster
Stacked_VlnPlot(hIA, features = c("Cebpd", "Ptx3"), pt.size = 0, split.by = "Timepoint2")

# you can plot raw counts as well
VlnPlot(hIA, 
        features = c("Spp1", "Gpnmb", "Fos", "Ch25h", "Cst7"), 
        pt.size = 0, 
        slot = "counts", 
        log = T,
        split.by = "Timepoint2")

FeaturePlot(hIA, features = c("Cst3"), label = T,split.by = "Timepoint2")

#Microglia: 
       #Homeo: Tmem119, Cx3cr1, Hexb, p2ry12, c1qa, c1qc, csf1r,, ctss, cst3
        #DAM: Cst7, Lpl, Apoe, trem2, axl, csf1, timp2, itgax, 
        #IRM: Ifitm3, Ifit3, Irf7
        #proliferative: Top2a, Birc5
        #chemokine expressing: ccl3, ccl4
        #micro-endo binding: Itga6, Itagal
#Astrocytes
        # Aqp4, gfap, aldoc, s100b, slc1a3, slc1a2, aldh1l1
#Endothelial Cells
        #Cldn5, Flt1, pecam1, slc2a1, nos3, ocln, tek, cdh5, bsg
#Pericyte
        # pdgfrb,Vtn, Ptn, ifitm1,kcnj8, rgs5
#Neurons
        # Excitatory: Slc17a7, Vsnl1, Snap25, Dnm1
        # Inhibitory: Gad2, Slc6a1
#OPC
        # Myt1
#Oligodendrocyte
        # Opalin, Plp1, mbp, mog, olig1, mag, cldn11, olig2, sox10
#smooth muscle cell
        #acta2, des, tagln, myl9
#Fibroblast
        #Col1a2
#PBMC
        #CD3e

#marker genes
FeaturePlot(hIA, features = c("Tmem119", "Cx3cr1", "Hexb", "P2ry12", "C1qa", "Csf1r"), label = T)
FeaturePlot(hIA, features = c("Apoe", "Cst7", "Lpl", "Csf1", "Itgax", "Trem2"), label = T)
FeaturePlot(hIA, features = c("Ifitm3", "Ifit3", "Irf7"), label = T)
FeaturePlot(hIA, features = c("Top2a", "Birc5"), label = T)
FeaturePlot(hIA, features = c("Ccl3", "Ccl4"), label = T)
FeaturePlot(hIA, features = c("Itga6", "Itgal"), label = T)
FeaturePlot(hIA, features = c("Gfap", "S100b", "Aldh1l1", "Aqp4", "Aldoc", "Slc1a3"), label = T)
FeaturePlot(hIA, features = c("Pecam1", "Tek", "Flt1", "Cldn5", "Ocln", "Slc2a1"), label = T)
FeaturePlot(hIA, features = c("Pdgfrb", "Vtn", "Ptn", "Ifitm1", "Kcnj8", "Rgs5"), label = T)
FeaturePlot(hIA, features = c("Slc17a7", "Vsnl1", "Snap25", "Dnm1"), label = T)
FeaturePlot(hIA, features = c("Gad2", "Slc6a1"), label = T)
FeaturePlot(hIA, features = c("Myt1"), label = T)
FeaturePlot(hIA, features = c("Mbp", "Mog", "Mag", "Plp1", "Olig1", "Olig2"), label = T)
FeaturePlot(hIA, features = c("Acta2", "Des", "Tagln", "Myl9"), label = T)
FeaturePlot(hIA, features = c("Col1a2"), label = T)
FeaturePlot(hIA, features = c("Cd3e"), label = T)

FeaturePlot_scCustom(hIA, features = c("Tgfbi", "Ifitm1", "Ifitm2", "Ifitm3"), label = T, num_columns = 2)
DimPlot(hIA, reduction = "umap", split.by = "Timepoint2", ncol = 2, label = T)

FeaturePlot(hIA, features = c("Tmem119", "Hexb", "Cx3cr1"), label = T, split.by = "Timepoint", by.col = T)
FeaturePlot(hIA, features = c("Lpl", "Trem2", "Csf1"), label = T, split.by = "Timepoint", by.col = T)
FeaturePlot(hIA, features = c("Aldoc", "Aqp4", "Aldh1l1"), label = T, split.by = "Timepoint", by.col = T)
FeaturePlot(hIA, features = c("Top2a", "Birc5"), label = T, split.by = "Timepoint", by.col = T)
FeaturePlot(hIA, features = c("Pecam1", "Tek", "Flt1"), label = T, split.by = "Timepoint", by.col = T)
FeaturePlot(hIA, features = c("Ccl3", "Ccl4"), label = T, split.by = "Timepoint", by.col = T)
FeaturePlot(hIA, features = c("Pdgfrb", "Vtn"), label = T, split.by = "Timepoint", by.col = T)
FeaturePlot(hIA, features = c("Olig1", "Olig2"), label = T, split.by = "Timepoint", by.col = T)
FeaturePlot(hIA, features = c("Spp1"), label = T, split.by = "Timepoint", by.col = T)

#heatmap for top 10 genes (pos and neg) for each cluster
markers_per_cluster<-hIA.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(hIA, features = top5$gene) + NoLegend()

# Get number of cells per cluster and per sample of origin
table(hIA@meta.data$RNA_snn_res.0.5, hIA@meta.data$Timepoint2)
cellspercluster<-table(hIA@meta.data$RNA_snn_res.0.5, hIA@meta.data$Timepoint2)
write.csv(cellspercluster, "cellspercluster.csv")

percentage_cellspercluster<-prop.table(cellspercluster, margin = 1)*100

table(hIA@meta.data$RNA_snn_res.0.5, hIA@meta.data$hIA_ID)

## 12. Adding Module scores ---------------------------------------------
  #homeo microglia
  microglia_homeo<-list(c('Tmem119', 'Cx3cr1', 'Hexb', 'P2ry12', 'C1qa', 'C1qc', 'Csf1r', 'Ctss', 'Cst3'))
  hIA<-AddModuleScore(hIA,
                      features = microglia_homeo,
                      name = "Micro_Homeo")
  # Plot scores
  FeaturePlot_scCustom(hIA,
                       features = "Micro_Homeo1", 
                       repel = TRUE,
                       order = T,
                       #split.by = c("Timepoint2"), 
                       keep.scale = 'all',
                       #num_columns = 2,
                       label = T,
                       na_cutoff = NA,
                       na_color = "light grey",
                       colors_use = viridis_magma_light_high) 
  #"For instance if plotting module score which contains negative values you will probably want to remove the cutoff value entirely to avoid misconstruing results."
  
  FeaturePlot(hIA,
              features = "Micro_Homeo1", 
              label = TRUE, 
              repel = TRUE)+
              #cols = c("red","white", "blue"))
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

  
  #DAM microglia
  microglia_dam<-list(c('Cst7', 'Lpl', 'Apoe', 'Trem2', 'Axl', 'Csf1', 'Timp2', 'Itgax'))
  #DAM: Cst7, Lpl, Apoe, trem2, axl, csf1, timp2, itgax, 
  hIA<-AddModuleScore(hIA,
                      features = microglia_dam,
                      name = "Micro_Dam")
  # Plot scores
  FeaturePlot_scCustom(hIA,
                       features = "Micro_Dam1", 
                       repel = TRUE,
                       order = T,
                       #split.by = c("Timepoint2"), 
                       #keep.scale = 'all',
                       num_columns = 2,
                       label = T,
                       na_cutoff = 0)
  
  #IRM microglia
  microglia_IRM<-list(c('Ifitm3', 'Ifit3', 'Irf7'))#IRM: Ifitm3, Ifit3, Irf7
  hIA<-AddModuleScore(hIA,
                      features = microglia_IRM,
                      name = "Micro_IRM")
  # Plot scores
  FeaturePlot_scCustom(hIA,
                       features = "Micro_IRM1", 
                       repel = TRUE,
                       order = T,
                       split.by = c("Timepoint2"), 
                       keep.scale = 'all',
                       num_columns = 2,
                       label = T,
                       na_cutoff = NA)
  
  
  #Proliferative microglia
  microglia_prolif<-list(c('Top2a', 'Birc5'))  #proliferative: Top2a, Birc5
  hIA<-AddModuleScore(hIA,
                      features = microglia_prolif,
                      name = "Micro_prolif")
  # Plot scores
  FeaturePlot_scCustom(hIA,
                       features = "Micro_prolif1", 
                       repel = TRUE,
                       order = T,
                       #split.by = c("Timepoint2"), 
                       keep.scale = 'all',
                       #num_columns = 2,
                       label = T,
                       na_cutoff = NA)
  

  #chemokine expressing microglia
  microglia_chemokine<-list(c('Ccl3', 'Ccl4'))          
  #chemokine expressing: ccl3, ccl4
  hIA<-AddModuleScore(hIA,
                      features = microglia_chemokine,
                      name = "Micro_chemokine")
  # Plot scores
  FeaturePlot_scCustom(hIA,
                       features = "Micro_chemokine1", 
                       repel = TRUE,
                       order = T,
                       #split.by = c("Timepoint2"), 
                       keep.scale = 'all',
                       #num_columns = 2,
                       label = T,
                       na_cutoff = NA)
  
  
  #Astrocytes
  astrocytes<-list(c('Aqp4', 'Gfap', 'Aldoc', 'S100b', 'Slc1a3', 'Slc1a2', 'Aldh1l1'))         
  hIA<-AddModuleScore(hIA,
                      features = astrocytes,
                      name = "Astrocytes")
  # Plot scores
  FeaturePlot_scCustom(hIA,
                       features = "Astrocytes1", 
                       repel = TRUE,
                       order = T,
                       #split.by = c("Timepoint2"), 
                       keep.scale = 'all',
                       #num_columns = 2,
                       label = T,
                       na_cutoff = NA)
  
  
  #Endothelial cells
  endos<-list(c('Cldn5', 'Flt1', 'Pecam1', 'Slc2a1', 'Nos3', 'Ocln', 'Tek', 'Cdh5'))         
  hIA<-AddModuleScore(hIA,
                      features = endos,
                      name = "Endos")
  # Plot scores
  FeaturePlot_scCustom(hIA,
              features = "Endos1", 
              repel = TRUE,
              order = T,
              split.by = c("Timepoint2"), 
              keep.scale = 'all',
              num_columns = 2,
              label = T,
              na_cutoff = NA)

  
  #Pericytes
  peris<-list(c('Pdgfrb','Vtn', 'Ptn', 'Ifitm1','Kcnj8', 'Rgs5'))         
  hIA<-AddModuleScore(hIA,
                      features = peris,
                      name = "Pericytes")
  # Plot scores
  FeaturePlot_scCustom(hIA,
                       features = "Pericytes1", 
                       repel = TRUE,
                       order = T,
                       split.by = c("Timepoint2"), 
                       keep.scale = 'all',
                       num_columns = 2,
                       label = T,
                       na_cutoff = NA)
  
  #OPCs/oligos
  Oligos<-list(c('Opalin', 'Plp1', 'Mbp', 'Mog', 'Olig1', 'Mag', 'Cldn11', 'Olig2', 'Sox10'))         
  hIA<-AddModuleScore(hIA,
                      features = Oligos,
                      name = "Oligos")
  # Plot scores
  FeaturePlot_scCustom(hIA,
                       features = "Oligos1", 
                       repel = TRUE,
                       order = T,
                       split.by = c("Timepoint2"), 
                       keep.scale = 'all',
                       num_columns = 2,
                       label = T,
                       na_cutoff = NA)

  #SMCs
  SMCs<-list(c('Acta2', 'Des', 'Tagln', 'Myl9'))         
  hIA<-AddModuleScore(hIA,
                      features = SMCs,
                      name = "SMCs")
  # Plot scores
  FeaturePlot_scCustom(hIA,
                       features = "SMCs1", 
                       repel = TRUE,
                       order = T,
                       split.by = "Timepoint2", 
                       keep.scale = 'all',
                       num_columns = 2,
                       label = T,
                       na_cutoff = NA)
  
  #Keren-Shaul - DAM stage 3
  list_1<-list(c('Cd9', 'Itgax', 'Clec7a', 'Cd63','Apoe', 'B2m', 'Tyrobp', 'Ctsd',
                 'Lpl', 'Cst7', 'Trem2', 'Ctsb','Fth1', 'Spp1', 'Axl', 'Lyz2',
                 'Igf1', 'Gpnmb', 'Lilrb4', 'Fabp5','Lgals3'))         
  hIA<-AddModuleScore(hIA,
                      features = list_1,
                      name = "list_1")
  # Plot scores
  FeaturePlot_scCustom(hIA,
                       features = "list_11", 
                       repel = TRUE,
                       order = T,
                       #split.by = c("Timepoint2"), 
                       keep.scale = 'all',
                       #num_columns = 2,
                       label = T,
                       na_cutoff = NA)
  
  #Mathys - Cluster 3
  mathys_cluster3<-list(c('Top2a', 'Uhrf1', 'Rrm2', 'Rad51','Chaf1b', 'Ccl3', 'Ccl4', 'Cxcl16',
                 'Mif'))         
  hIA<-AddModuleScore(hIA,
                      features = mathys_cluster3,
                      name = "mathys_cluster3")
  # Plot scores
  FeaturePlot_scCustom(hIA,
                       features = "mathys_cluster31", 
                       repel = TRUE,
                       order = T,
                       #split.by = c("Timepoint2"), 
                       keep.scale = 'all',
                       #num_columns = 2,
                       label = T,
                       na_cutoff = NA)
  
  #Hammond Cluster 3 - E14.5 microglia genes
  hammond_cluster3<-list(c('Fabp5','Ldha', 'Lgals1', 'Spp1', 'Mif','Tpi1', 'Pkm',
                           'Aldoa', 'Ftl1', 'Mif'))         
  hIA<-AddModuleScore(hIA,
                      features = hammond_cluster3,
                      name = "hammond_cluster3")
  # Plot scores
  FeaturePlot_scCustom(hIA,
                       features = "hammond_cluster31", 
                       repel = TRUE,
                       order = T,
                       split.by = c("Timepoint2"), 
                       keep.scale = 'all',
                       num_columns = 2,
                       label = T,
                       na_cutoff = NA)
  
  #Hammond Cluster 3 - E14.5 microglia genes
  BAM_inquiry<-list(c("Mrc1", "Maf", "Fli1", "Etv1", "Runx1", "Trps1"))         
  hIA<-AddModuleScore(hIA,
                      features = BAM_inquiry,
                      name = "BAM_inquiry")
  # Plot scores
  FeaturePlot_scCustom(hIA,
                       features = "BAM_inquiry1", 
                       repel = TRUE,
                       order = T,
                       #split.by = c("Timepoint2"), 
                       keep.scale = 'all',
                       num_columns = 2,
                       label = T,
                       na_cutoff = NA)
  

# Violin plot by cluster and gene marker 
  colors_1<-c("grey80", "grey29", "grey60", "grey11")
  
  gene_list_plot <- c("Tnf")
  Stacked_VlnPlot(seurat_object = hIA, features = gene_list_plot, x_lab_rotate = TRUE, split.by = "Timepoint2", plot_legend = T, group.by =)
  
  marsh_gene_list_plot <- c("Slc17a7", "Gad2", "Aqp4", "Myt1", "Col1a2", "Cldn5", "Opalin", "Cx3cr1", "Cd3e")
  Stacked_VlnPlot(seurat_object = hIA, features = marsh_gene_list_plot, x_lab_rotate = TRUE, split.by = "Timepoint2", plot_legend = T)

  micro_gene_list<-c("Tmem119", "Cx3cr1", "Hexb", "P2ry12")
  Stacked_VlnPlot(seurat_object = hIA, 
                  features = micro_gene_list, 
                  x_lab_rotate = TRUE, 
                  #split.by = "Timepoint2",
                  #idents = c(0,1,2,5,6,7,9,12,17),
                  plot_legend = F,
                  colors_use = color_pal)
  
  micro_gene_list2<-c("Cst7", "Lpl", "Apoe", "Trem2", "Axl", "Csf1", "Timp2", "Itgax")
  Stacked_VlnPlot(seurat_object = hIA, 
                  features = micro_gene_list2, 
                  x_lab_rotate = TRUE, 
                  split.by = "Timepoint2", 
                  plot_legend = T)
  
  yang_gene_list<-c("Tmem119", "Cx3cr1", "Cst7", "Clec7a", "Apoe", "Ifitm3", "Hexb", "C3ar1", "Stmn1")
  Stacked_VlnPlot(seurat_object = hIA, 
                  features = yang_gene_list, 
                  x_lab_rotate = TRUE, 
                  split.by = "Timepoint2",
                  plot_legend = T)
  
  micro_gene_list3<-c("Ifitm3", "Ifit3", "Irf7", "Top2a", "Birc5", "Ccl3", "Ccl4", "Itga6", "Itgal")
  Stacked_VlnPlot(seurat_object = hIA, 
                  features = micro_gene_list3, 
                  x_lab_rotate = TRUE, 
                  split.by = "Timepoint2",
                  plot_legend = T)
  
  astro_gene_list1<-c("Gfap", "S100b", "Aqp4", "Aldoc", "Vim")
  Stacked_VlnPlot(seurat_object = hIA, 
                  features = astro_gene_list1, 
                  x_lab_rotate = TRUE, 
                #  split.by = "Timepoint2",
                  plot_legend = F,
                  colors_use = color_pal)
  
  endo_gene_list1<-c("Pecam1", "Tek", "Flt1", "Cldn5", "Ocln", "Slc2a1")
  Stacked_VlnPlot(seurat_object = hIA, 
                  features = endo_gene_list1, 
                  x_lab_rotate = TRUE, 
                  split.by = "Timepoint2",
                  plot_legend = T)
  
  peri_gene_list1<-c("Pdgfrb", "Vtn", "Ptn", "Ifitm1", "Kcnj8", "Rgs5")
  Stacked_VlnPlot(seurat_object = hIA, 
                  features = peri_gene_list1, 
                  x_lab_rotate = TRUE, 
                  split.by = "Timepoint2",
                  plot_legend = T)
  
  oligo_gene_list1<-c("Mbp", "Mog", "Mag", "Plp1", "Olig1", "Olig2", "Opalin")
  Stacked_VlnPlot(seurat_object = hIA, 
                  features = oligo_gene_list1, 
                  x_lab_rotate = TRUE, 
                  split.by = "Timepoint2",
                  plot_legend = T)
  
  test<-c("Cldn11","Enpp2","Mag","Mal", "Ctss", "Hexb")
  Stacked_VlnPlot(seurat_object = hIA, 
                  features = test, 
                  x_lab_rotate = TRUE, 
                  split.by = "Timepoint2",
                  idents = c(0,1,5,6,7,4,8,11,2),
                  colors_use = colors_1,
                  plot_legend = T)
  
  fig_gene_list1<-c("Tmem119","Hexb", "Ctss")
  fig_gene_list2<-c("Aldoc", "Aqp4") 
  fig_gene_list3<-c("Opalin", "Plp1", "Myt1")
  fig_gene_list4<-c("Gad2", "Slc17a7")
  #"Cldn5", "Pdgfrb", "Opalin", "Cd3e")
  #colors_1<-c("Dark Blue", "Light Blue", "Dark Green", "Light Green")
  plot_1<-Stacked_VlnPlot(seurat_object = hIA, 
                  features = fig_gene_list1, 
                  x_lab_rotate = TRUE, 
                  plot_legend = F,
                  colors_use = color_pal)
  plot_2<-Stacked_VlnPlot(seurat_object = hIA, 
                          features = fig_gene_list2, 
                          x_lab_rotate = TRUE, 
                          plot_legend = F,
                          colors_use = color_pal)
  plot_3<-Stacked_VlnPlot(seurat_object = hIA, 
                          features = fig_gene_list3, 
                          x_lab_rotate = TRUE, 
                          plot_legend = F,
                          colors_use = color_pal)
  plot_4<-Stacked_VlnPlot(seurat_object = hIA, 
                          features = fig_gene_list4, 
                          x_lab_rotate = TRUE, 
                          plot_legend = F,
                          colors_use = color_pal)
 (plot_1 | plot_2) / 
    (plot_3 | plot_4)
  
  micro_gene_list<-c("Tmem119", "Cx3cr1", "Hexb", "P2ry12", "C1qa", "C1qc", "Csf1r","Ctss", "Cst3","Ifitm3", "Ifit3", "Irf7", "Top2a", "Birc5", "Ccl3", "Ccl4", "Itga6", "Itgal")
  Stacked_VlnPlot(seurat_object = hIA, 
                  features = micro_gene_list, 
                  x_lab_rotate = TRUE, 
                  split.by = "Timepoint2",
                  idents = c(0,1,2,5,6,7,9,12,17),
                  plot_legend = T,
                  colors = colors_1)
  
  cluster6_probe1<-c("C1qa", "C1qc","Ctss", "Lpl", "Cst7","Apoe", "Ccl3", "Ccl4")
  Stacked_VlnPlot(seurat_object = hIA, 
                  features = cluster6_probe1, 
                  x_lab_rotate = TRUE, 
                  split.by = "Timepoint2",
                  idents = c(6),
                  plot_legend = T,
                  colors = color_pal)
  
    color_pal_micro <-c("#AA0DFEFF","#3B00FBFF","#7ED7D1FF","#1C8356FF","#16FF32FF","#B00068FF","#FE00FAFF" )
    test2<-c("Tmem119",
             "Hexb",
             "Apoe",
             "Tnf",
             "Il1b",
             "Spp1",
             "Lpl",
             "Cst3",
             "Top2a",
             "Clu",
             "Aldoc")
  Stacked_VlnPlot(seurat_object = hIA, 
                  features = test2, 
                  x_lab_rotate = TRUE, 
                  #split.by = "Timepoint2",
                  idents = c(0,1,5,6,7,9,12),
                  colors_use = color_pal_micro,
                  plot_legend = T)
  
  #immune marker genes
color_pal_immune <-c("#D85FF7FF","#00868B","#2ED9FFFF")
"Itgal", "Il1b"
test3<-c("Cx3cr1", "Tmem119", "Pf4", "Mrc1", "Lyz2", "Ifitm1", "Ifitm2","Cxcl2", "Cd74","Cd3e","Itgal" )
test4<-c("Tgfb1", "Tgfb2", "Tgfbr1", "Tgfbr2", "Acvr1")
color_two<-c("#F3766E", "#18BDC2")
Stacked_VlnPlot(seurat_object = hIA_3day, 
                  features = test4, 
                  x_lab_rotate = TRUE, 
                  split.by = "Timepoint2",
                  split.plot = T,
                  idents = c(0,1,5,6,7,9,12,14,16,19),
                  plot_legend = T, 
                colors_use = color_two)

  
percentage_cellspercluster <-as.data.frame(percentage_cellspercluster)
percentage_cellspercluster<-percentage_cellspercluster %>% ## unload and reload dplyr if doens't work
  rename(
    Cluster = Var1,
    Timepoint2 = Var2, 
    Percentage = Freq
  )
  
cluster0_percentages<-percentage_cellspercluster[which(percentage_cellspercluster$Cluster == 0),]
ggplot(cluster0_percentages, aes(x = Timepoint2, y = Percentage))+
  geom_bar(stat = "identity", aes(fill = Cluster))+
  scale_fill_manual(values=color_pal[1])+
  theme_classic()
  

## 13. Perform Default differential expression tests

  
  #write out the markers but only the p<0.05 genes
  hIA.markers_signif<-hIA.markers[which(hIA.markers$p_val_adj<0.05),]
  write.csv(hIA.markers_signif, file = "hIA_DEgenes_allclusters_significant.csv")
  
  
#find differences between clusters nearby that have similar marker genes
setwd("~/Documents/Projects/hIA/nfeature_5000/DEgenelists/")
# between 0 and 1 
C12_C01567_DEgenes<-FindMarkers(hIA, ident.1 = "Unknown2", ident.2 = c("Micro1", "Micro2", "Micro4", "Micro-DAM"))
head(C12_C01567_DEgenes, n = 30)
write.csv(C12_C01567_DEgenes, file = "C12_C01567_DEgenes.csv")

# between X and Y 
DEgenes<-FindMarkers(hIA, ident.1 = "Oligo1", ident.2 = "Oligo2")
head(DEgenes, n = 20)

write.csv(DEgenes, file = "DEgenes_C4_v_C8.csv")

#Seurat7 endpoint Save
