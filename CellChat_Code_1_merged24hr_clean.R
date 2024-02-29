### Cell Chat Code

# https://www.biorxiv.org/content/10.1101/2023.11.05.565674v1.full.pdf 
#https://rdrr.io/github/sqjin/CellChat/f/tutorial/CellChat-vignette.Rmd

#install.packages("remotes")
#remotes::install_github("sqjin/CellChat")

library(CellChat)
library(patchwork)
library(Seurat)
library(umap)
library(patchwork)
library(beepr)
library(ggplot2)
library(ComplexHeatmap)
library(SeuratWrappers)
library(pheatmap)
library(ggplotify)
options(stringsAsFactors = FALSE)

setwd("/Users/foley/Documents/Projects/hIA/nfeature_5000/CellChat")
cellchat_24hrAB<-readRDS("cellchat_24AB_hIA_savepoint1.rds")
cellchat_24hrIgG<-readRDS("hIA_24hr_IgG_hIA_savepoint1.rds")
#cellchat_3dayAB<-readRDS("hIA_3day_AB_hIA_savepoint1.rds")
#cellchat_3dayIgG<-readRDS("hIA_3day_IgG_savepoint1.rds")

cellchat_24hrAB <- netAnalysis_computeCentrality(cellchat_24hrAB)
cellchat_24hrIgG <- netAnalysis_computeCentrality(cellchat_24hrIgG)
#cellchat_3dayAB <- netAnalysis_computeCentrality(cellchat_3dayAB)
#cellchat_3dayIgG <- netAnalysis_computeCentrality(cellchat_3dayIgG)

object.list <- list(cc_24hr_AB = cellchat_24hrAB,
                    cc_24hr_IgG = cellchat_24hrIgG)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#save(object.list, file = "cellchat_object.list_24hr_ABIgG.RData")
#save(cellchat, file = "cellchat_merged_24hr_ABIgG.RData")

setwd("/Users/foley/Documents/Projects/hIA/nfeature_5000/CellChat/24hr_ABIgG/")

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


par(mfrow = c(2,2))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2


group.cellType <- c("Micro",  #0 "#AA0DFEFF", 
                    "Micro",   #1 "#3B00FBFF",
                    "UK",   #2 "#66B0FFFF",
                    "Astro",   #3 "#FA0087FF", 
                    "Oligo/OPC",   #4 "#782AB6FF", #4
                    "Micro",  #5 "#7ED7D1FF", 
                    "Micro",  #6 "#1C8356FF", 
                    "Micro",  #7 "#16FF32FF", #7
                    "Oligo/OPC", #8 "#F7E1A0FF",
                    "Micro", #9 "#B00068FF", 
                    "Endo", #10 "#F6222EFF", #10
                    "Oligo/OPC", #11 "#DEA0FDFF", 
                    "Micro", #12 "#FE00FAFF", #12
                    "Neuron", #13 "#325A9BFF",
                    "Immune", #14 "#D85FF7FF",
                    "Neuron", #15 "#F8A19FFF", #15
                    "Immune", #16 "#7ED7D1FF",
                    "Astro", #17 "#FBE426FF",
                    "UK", #18 "#FEAF16FF", #18
                    "Immune", #19 "#2ED9FFFF", 
                    "UK", #20 "#3283FEFF",
                    "UK", #21  "#FC1CBFFF" , #21
                    "Endo")


group.cellType <- factor(group.cellType, 
                         levels = c("Micro",  #0 "#AA0DFEFF", 
                                    "Astro",   #1 "#3B00FBFF",
                                    "Endo",   #2 "#66B0FFFF",
                                    "Oligo/OPC",   #3 "#FA0087FF", 
                                    "Immune",   #4 "#782AB6FF", #4
                                    "Neuron",  #5 "#7ED7D1FF", 
                                    "UK"  #6 "#1C8356FF", 
                                      ))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))



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


# Show the number of interactions or interaction strength betweenany two cell types in each dataset.
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, 
                   weight.scale = T, 
                   label.edge= T, 
                   edge.weight.max = weight.max[3], 
                   edge.width.max = 12, 
                   shape = "circle",
                   #color.use = color_pal,
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

object.list[[1]]@net$count.merged
object.list[[2]]@net$count.merged


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax, color.use = color_pal)
}
patchwork::wrap_plots(plots = gg)

#isolate the scatterplot data 
AB_scatter<-gg[[1]]$data
IgG_scatter<-gg[[2]]$data

#remove all non-microglia/immune clusters
AB_scatter<-AB_scatter[c(1,2,6,7,8,10,13,15,17,20),]
IgG_scatter<-IgG_scatter[c(1,2,6,7,8,10,13,15,17,20),]

color_pal_microimmune<-c("#AA0DFEFF","#3B00FBFF","#7ED7D1FF","#1C8356FF","#16FF32FF","#B00068FF","#FE00FAFF","#D85FF7FF","#00868B","#2ED9FFFF" )

ggplot(AB_scatter, aes(x = x, y = y, size = Count, label = labels))+
 geom_point(color = color_pal_microimmune)+
  scale_size_area(breaks = c(250,500,750,1000, 1250))+
 geom_text(hjust = 0, nudge_x = .5, size = 5)+
  scale_color_manual(values = color_pal_microimmune)+
  xlim(0,20)+
  ylim(0,20)+
  theme_bw()


# Similarly, CellChat can also show the differential number of interactions or interaction strength between any two cell types using circle plot.
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)



#Identify cell populations with significant changes in sending or receiving signals
num.link <- sapply(object.list, function(x) {rowSums(x@net$count)
  + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],title = names(object.list)[i], weight.MinMax = weight.MinMax, color.use = color_pal)
}
patchwork::wrap_plots(plots = gg)


gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
#https://github.com/sqjin/CellChat/issues/63

#endo:11,23
#micro: 1,2,6,7,8,10,13
#astro: 4,18
#immune: 15, 17, 20
#oligo/OPC: 5,9,12
#neuron: 14,16
#uk:

gg1 <- rankNet(cellchat, slot.name = "net", mode = "comparison", stacked = T, do.stat = TRUE, 
               sources.use = c(4),
               targets.use = c(1:23),
               do.flip = T, 
               thresh = 0.05,
               cutoff.pvalue = 0.05)
gg1

net_infoflow_data<-gg1$data

write.csv(net_infoflow_data, "netinfoflow_24hr_C14_to_endos.csv")


i = 1
# combining all the identified signaling pathways from different datasets
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i], 
                                        width = 10, height= 20,
                                        color.heatmap = "GnBu",
                                        color.use = color_pal,
                                        cluster.rows = F,
                                        cluster.cols = F)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i+1], 
                                        width = 10, height = 20,
                                        color.heatmap = "GnBu",
                                        color.use = color_pal,
                                        cluster.rows = F,
                                        cluster.cols = F)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", 
                                        #signaling = pathway.union, 
                                        title = names(object.list)[i], 
                                        width = 10, height = 20, 
                                        color.heatmap = "GnBu",
                                        color.use = color_pal)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", 
                                        signaling = pathway.union, title = names(object.list)[i+1], 
                                        width = 10, height = 20, 
                                        color.heatmap = "GnBu",
                                        color.use = color_pal)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], 
                                        pattern = "all", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i], 
                                        width = 10, height = 20, 
                                        color.heatmap = "GnBu",
                                        color.use = color_pal)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], 
                                        pattern = "all", 
                                        signaling = pathway.union, 
                                        title = names(object.list)[i+1], 
                                        width = 10, height = 20, 
                                        color.heatmap = "GnBu",
                                        color.use = color_pal)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#play with this one to get ideas of commmuncation 
netVisual_bubble(cellchat, 
                 sources.use = c(15), 
                 targets.use = c(12,6,7,8,10,13),
                 sort.by.target = T,
                 comparison = c(1, 2), angle.x = 45,
                 remove.isolate = T, 
                 thresh = 0.01)
gg3<-netVisual_bubble(cellchat, 
                 signaling = "TGFb",
                 sources.use = c(1,2,6,7,8,10,13), 
                 targets.use = c(15,17,20),  
                 comparison = c(1, 2), angle.x = 90,
                 return.data = F, 
                 thresh = 0.05,
                 remove.isolate = F)
gg3

#micro: c(1,2,6,7,8,10,13,17)
#UK: c(3,19,21,22)
#Astro: c(4,18)
#Oligo/OPC: c(5,9,12)
#endo: c(11,23)
#immune: c(15,20)
#neuron: c(14,16)

gg1 <- netVisual_bubble(cellchat, sources.use = c(15), targets.use = c(11,23),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Control", angle.x = 45, remove.isolate = T, thresh = 0.01)
gg2 <- netVisual_bubble(cellchat, sources.use = c(15), targets.use = c(11,23),  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in 24hr ABinjection", angle.x = 45, remove.isolate = T, thresh = 0.01)
gg2
gg1 + gg2
gg1
gg2


gg2 <- netVisual_bubble(cellchat, sources.use = c(15,17,20), targets.use = c(1:23),
                        signaling = "MIF",
                        comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in 24hr ABinjection", angle.x = 45, remove.isolate = T, thresh = 0.01)
gg2

gg4 <- netVisual_bubble(cellchat, sources.use = c(15,17,20), targets.use = c(1:23),
                        signaling = "TNF",
                        comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in 24hr IgG", angle.x = 45, remove.isolate = T, thresh = 0.01)
gg4


gg2 <- netVisual_bubble(cellchat, 
                        sources.use = c(15, 17, 20), 
                        targets.use = c(1:23),
                        signaling = "MIF",
                        comparison = c(1, 2), 
                        #max.dataset = 1, 
                        #title.name = "Increased signaling in ABinjection",
                        angle.x = 45, remove.isolate = T, thresh = 0.05)
gg2

netVisual_bubble(cellchat, 
                        sources.use = c(1,2,6,7,8,10,13), 
                        targets.use = c(15,17,20),  
                        comparison = c(1, 2), 
                        max.dataset = 1, 
                        #title.name = "Increased signaling in ABinjection", 
                        angle.x = 45, remove.isolate = T,
                 thresh = 0.05)




# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "cc_24hr_IgG"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "cc_24hr_IgG",ligand.logFC = 0.2, receptor.logFC = 0.2)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "cc_24hr_AB",ligand.logFC = -0.2, receptor.logFC = -0.2)


gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)


pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, 
                        sources.use = c(15,17,29), 
                        targets.use = c(1:23), 
                        signaling = "TNF",
                        comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, 
                        sources.use = c(15), 
                        targets.use = c(1:23), 
                        #signaling = "TNF",
                        comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2


par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(1:23), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(1:23), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))


# visualize the enriched ligands in the first condition
computeEnrichmentScore(net.down, species = 'mouse')

# visualize the enriched ligands in the second condition
#computeEnrichmentScore(net.up, species = 'mouse')


pathways.show <- c("TGFb") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "circle", 
                      sources.use = c(1,2,6,7,8,10,13),
                      targets.use = c(15,17,20),
                      edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      color.use = color_pal, 
                      signaling.name = paste(pathways.show, 
                                             names(object.list)[i]))
}

pathways.show <- c("TNF") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], 
                               signaling = pathways.show,
                               slot.name = c("netP"),
                               color.heatmap = "GnBu",
                               color.use = color_pal,
                               title.name = paste(pathways.show, "signaling ",
                                                  names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

AB_pathway_signaling<-ht[[1]]@matrix
IgG_pathway_signaling<-ht[[2]]@matrix

write.csv(AB_pathway_signaling, paste0("AB_", pathways.show, "_Signaling.csv"))
write.csv(IgG_pathway_signaling, paste0("IgG_", pathways.show, "_Signaling.csv"))


# Chord diagram
pathways.show <- c("TNF") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "chord", 
                      #sources.use = c(1,2,6,7,8,10,13),
                      #targets.use = c(1,2,6,7,8,10,13, 15,17,20),
                      color.use = color_pal,
                     weight.scale = T,
                      signaling.name = paste(pathways.show, names(object.list)[i]),
                      show.legend = F, 
                     thresh = 0.05)
}



# Chord diagram
#group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(object.list[[1]]@idents)
pathways.show <- c("TGFb") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]), color.use = color_pal)
}
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show,
                       sources.use = c(11),
                       targets.use = c(7),
                       group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]), color.use = color_pal)
}


par(mfrow = c(1, 2), xpd=TRUE)

# compare all the interactions sending from fibroblast to inflamatory immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], 
                       sources.use = c(1,2,6,7,8,10,13), 
                       targets.use = c(2),  
                       title.name = paste0("Signaling received by Microglia to Endo - ", 
                                           thresh = 0.0001,
                                           names(object.list)[i]), legend.pos.x = 10)}
#micro: c(1,2,6,7,8,10,13,17)
#UK: c(3,19,21,22)
#Astro: c(4,18)
#Oligo/OPC: c(5,9,12)
#endo: c(11,23)
#immune: c(15,20)
#neuron: c(14,16)



saveRDS(cellchat, file = "cellchat_comparisonAnalysis_24hours_ABIgG.rds")




cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("cc_24hr_IgG", "cc_24hr_AB")) # set factor level
plotGeneExpression(cellchat, signaling = "IGF", split.by = "datasets", colors.ggplot = T,
                   enriched.only = T, split.plot = T)



gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = c("C19"))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = c("C10"))
patchwork::wrap_plots(plots = list(gg1,gg2))



# export various datasets 
#weight is the strength
#count is the number of interactions
AB_weightmerged<-as.data.frame(cellchat@net$cc_24hr_AB$weight.merged)
AB_countmerged<-cellchat@net$cc_24hr_AB$count.merged
AB_pval<-as.data.frame(cellchat@net$cc_24hr_AB$pval)
AB_weight<-as.data.frame(cellchat@net$cc_24hr_AB$weight)
AB_count<-as.data.frame(cellchat@net$cc_24hr_AB$count)

IgG_weightmerged<-cellchat@net$cc_24hr_IgG$weight.merged
IgG_countmerged<-cellchat@net$cc_24hr_IgG$count.merged
IgG_pval<-as.data.frame(cellchat@net$cc_24hr_IgG$pval)
IgG_weight<-as.data.frame(cellchat@net$cc_24hr_IgG$weight)
IgG_count<-as.data.frame(cellchat@net$cc_24hr_IgG$count)

write.csv(AB_weightmerged, file = "24hr_AB_weightmerge.csv")
write.csv(AB_countmerged, file = "24hr_AB_countmerge.csv")
write.csv(IgG_weightmerged, file = "24hr_IgG_weightmerge.csv")
write.csv(IgG_countmerged, file = "24hr_IgG_countmerge.csv")

write.csv(AB_weight, file = "24hr_AB_weight.csv")
write.csv(AB_count, file = "24hr_AB_count.csv")
write.csv(IgG_weight, file = "24hr_IgG_weight.csv")
write.csv(IgG_count, file = "24hr_IgG_count.csv")

subtraction<-as.data.frame(AB_weightmerged-IgG_weightmerged)



par(mfrow = c(1,2))


p1<-as.ggplot(pheatmap(AB_weight, 
         display_numbers = T,
         cluster_cols = F,
         cluster_rows = F,
         fontsize_number = 5))

p2<-as.ggplot(pheatmap(IgG_weight, 
         display_numbers = T,
         cluster_cols = F,
         cluster_rows = F,
         fontsize_number = 5))

p1+p2

test2<-ht1@matrix
test3<-ht2@matrix

subtraction<-as.data.frame(test2-test3)
subtraction[is.na(subtraction)]=0

as.ggplot(pheatmap(subtraction, 
                   display_numbers = F,
                   cluster_cols = F,
                   cluster_rows = F,
                   fontsize_number = 0))



p2<-as.ggplot(pheatmap(test2, 
                       display_numbers = T,
                       cluster_cols = F,
                       cluster_rows = F,
                       fontsize_number = 0))

test2[is.na(test2)]=0
test3[]

heatmap(test2)
