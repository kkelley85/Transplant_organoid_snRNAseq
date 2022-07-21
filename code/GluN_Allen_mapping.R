#Excitatory neuron mapping to Allen brain labels: Extended Data Figure 4

library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

#Load processed data

thcs=readRDS("GSE190815_t-hCS_processed_SeuratObject.rds")

#t-hCS mapping example

#subset to GluN neurons

Idents(thcs)="subclass_label"
clust=subset(thcs,idents=c("GluN","GluN_UL","GluN_DL","GluN_DL/SP"))

#re-analyze and integrate subsetted GluN neurons

data=SplitObject(clust,split.by="Plot.name")
for (i in c(1:length(data))){
	print(paste0("Processing ",unique(data[[i]]@meta.data$Plot.name)))
	counts=data[[i]]@assays$RNA@counts
	meta=data[[i]]@meta.data
	meta=meta[,c(1:grep("Plot.name",colnames(meta)), grep("seurat_clusters",colnames(meta)), grep("CellNames",colnames(meta)))]
	colnames(meta)[grep("seurat_clusters",colnames(meta))]="Full_seuratclusters_res0.5"
	
	data[[i]]=CreateSeuratObject(counts = counts, meta.data=meta,project = meta$Plot.name[1], min.cells = 5, min.features = 0)
	data[[i]]=SCTransform(data[[i]], vst.flavor = "v2",verbose = T)
	data[[i]]= RunPCA(data[[i]],npcs=100)
	
}


features <- SelectIntegrationFeatures(object.list = data, nfeatures = 3000)
data <- PrepSCTIntegration(object.list = data, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = data, normalization.method = "SCT",anchor.features = features)

data.combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims=1:30)
data.combined <- RunPCA(data.combined,npcs=100)
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:30, verbose = FALSE)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:30)
data.combined <- PrepSCTFindMarkers(data.combined)

query= data.combined

#Load preprocessed reference datast, e.g. Allen M1; obtained here: https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x

reference=readRDS("SeurObj_Allen_M1_dims_30_res_0.5.rds")

#select only class of interest for integration
Idents(reference) <- reference$class_label
reference <- subset(reference, idents = "Glutamatergic")

#downsample to max of 500 cells per cluster for computational speed
Idents(reference) <- reference$cluster_label
reference.ds <- subset(reference, downsample = 500)

#Normalize subsetted data
reference.ds= SCTransform(reference.ds, vst.flavor = "v2",verbose = T)
reference.ds <- RunPCA(reference.ds,npcs=100)
reference.ds <- RunUMAP(reference.ds, reduction = "pca", dims = 1:30, verbose = FALSE)
reference.ds <- FindNeighbors(reference.ds, reduction = "pca", dims = 1:30)
#reference.ds <- PrepSCTFindMarkers(reference.ds)


	anchors <- FindTransferAnchors(
  		reference = reference.ds,
  		query = query,
  		normalization.method = "SCT",
  		reference.reduction="pca",
  		npcs=30,
  		dims=1:30)

	tmp=reference.ds@meta.data[,is.element(colnames(reference.ds@meta.data),"cluster_label")]
	names(tmp)=colnames(reference.ds)

	predictions <- TransferData(
		anchorset = anchors,
		refdata = tmp,
    	dims = 1:30)

	query <- AddMetaData(query, metadata = predictions[,1],col.name="M1_clusterLabel")

	tmp=reference.ds@meta.data[,is.element(colnames(reference.ds@meta.data),"subclass_label")]
	names(tmp)=colnames(reference.ds)
	
	predictions <- TransferData(
		anchorset = anchors,
		refdata = tmp,
    	dims = 1:30)
	
colnames(predictions)[1]="M1_Subclass"

	query <- AddMetaData(query, metadata = predictions[,1],col.name="M1_Subclass")


#make plots on original UMAP coordinates

Idents(thcs)="subclass_label"
clust=subset(thcs,idents=c("GluN","GluN_UL","GluN_DL","GluN_DL/SP"))

all.equal(colnames(query),colnames(clust))

clust[["M1_clusterLabel"]]=query$M1_clusterLabel
clust[["M1_Subclass"]]=query$M1_Subclass

labels= c("L2/3 IT","L5 IT","L5 ET","L6 IT","L6 IT Car3","L5/6 NP","L6 CT","L6b")
color_labels=c("#9E0142","#F46D43","#FDAE61","#FFFFBF","#E6F598","#66C2A5","#3288BD","#5E4FA2")
names(color_labels)=labels

p1= DimPlot(clust, reduction = "umap", cols= color_labels ,group.by="M1_Subclass",shuffle=T)&NoAxes()& ggtitle(paste0(""))&theme(text = element_text(size=6))

	setwd("output")

	pdf(file=paste0("t-hCS_GluN_M1_Subclass_predictions_UMAP.pdf"),width=4,height=3)
	print(p1)
	dev.off()