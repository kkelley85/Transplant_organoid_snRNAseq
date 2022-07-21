#Code to generate cluster overlap plots from reference datasets
#Requires reference datasets to be downloaded and saved as seurat objects
#Code adopted from Bakken et al. Nature 2021: https://www.nature.com/articles/s41586-021-03465-8
#Original code from Bakken et al. found here: https://github.com/AllenInstitute/BICCN_M1_Evo/blob/master/human_m1_mtg_integration/Figure_5_b_c_e_f___ED_figure_8_all.R

library(Seurat)
library(ggplot2)
library(pheatmap)

#Load query dataset (either t-hCS or hCS)

query=readRDS("GSE190815_t-hCS_processed_SeuratObject.rds")

#downsample for computational speed and to balance dataset sizes
Idents(query)="cluster_label"
query=subset(query, downsample =500)

query[["Condition"]]="t-hCS"

#Load pre-processed reference dataset
#For example: Allen MTG; raw data obtained here: https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-smart-seq

reference=readRDS("SeurObj_Allen_MTG_dims_30_res_0.5.rds")

Idents(reference) <- reference$cluster
reference =subset(reference, downsample =500)

reference[["Condition"]]="MTG"
reference[["Plot.name"]]="MTG"

#Comparison to MTG

ref="MTG"
quer="t-hCS"

#Merge and normalize each sample
data=merge(query,reference)
data <- SplitObject(data, split.by = "Plot.name")

for (i in c(1:length(data))){
	
	data[[i]]=SCTransform(data[[i]], vst.flavor = "v2",verbose = T)
	data[[i]]= RunPCA(data[[i]],npcs=100)

}

#Integrate datasets
features <- SelectIntegrationFeatures(object.list = data, nfeatures = 3000)
data <- PrepSCTIntegration(object.list = data, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = data, normalization.method = "SCT",anchor.features = features)

data.combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims=1:30)
data.combined <- RunPCA(data.combined,npcs=100)
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:30, verbose = FALSE)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:30)

#cluster
clust= FindClusters(data.combined, resolution = 0.5)


#cluster label
clust[["seurat_clusters.new"]]=as.numeric(clust$seurat_clusters)


label=paste0(quer,"_",clust@meta.data$CellNames,"_", clust@meta.data$Full_seuratclusters_res0.5)
label[is.na(clust@meta.data$Full_Original_clusters)]=paste0(ref,"_",clust@meta.data$cluster[is.na(clust@meta.data$Full_Original_clusters)])

clust$label_for_heatmap <- label

# Compare clustering
#Code adopted from Bakken et al. Nature 2021: https://www.nature.com/articles/s41586-021-03465-8
#Original code from Bakken et al. found here: https://github.com/AllenInstitute/BICCN_M1_Evo/blob/master/human_m1_mtg_integration/Figure_5_b_c_e_f___ED_figure_8_all.R

#load Allen functions
source('Allen_Coclustering_fxns.R', chdir = TRUE)

#Setup co-clustering for plotting
heat.colors <- colorRampPalette(c("white", "grey70", "black"))(100)
heat.colors <- colorRampPalette(c("darkblue", "darkgreen", "yellow"))(100)

ref.cl <- clust$label_for_heatmap
cca.cl <- clust$seurat_clusters.new
compare.datasets <- unique(clust$Condition)

cl.conf <- compare_cl(ref.cl, cca.cl)
cocl <- cl.conf$cocl

cocl.subset <- cocl[grepl(compare.datasets[1], row.names(cocl)),
                    grepl(compare.datasets[2], row.names(cocl))]  

row.names(cocl.subset) <- sub(paste0(compare.datasets[1], "_"), "",
                              row.names(cocl.subset))

colnames(cocl.subset) <- sub(paste0(compare.datasets[2], "_"), "",  
                             colnames(cocl.subset))

hCS_order=c("Cyc_Prog_15","Progenitor_3","Astroglia_9","Astroglia_5","Astroglia_4","OPC_1","Oligo_14","UL_GluN_13","UL_GluN_2","UL_GluN_0","UL_GluN_7","UL_GluN_6","DL_GluN_11","DL_GluN_8","DL/SP_GluN_10","DL_ImGluN_12","RELN_16")

cocl.subset2 <- reorder_matrix(cocl.subset[hCS_order,], by.rows = F)

cocl.subset2= cocl.subset2[,c(colnames(cocl.subset2)[c(1:length(colnames(cocl.subset2)))[-grep("Inh", colnames(cocl.subset2))]],colnames(cocl.subset2)[grep("Inh", colnames(cocl.subset2))])]

#create heatmap plot
pheatmap(cocl.subset2, cluster_rows = F, cluster_cols = F,
         color = heat.colors, 
         fontsize = 6)
