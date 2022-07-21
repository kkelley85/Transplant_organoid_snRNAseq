#Code to preprocess 10x single-nuclei RNAseq CellRanger Count v6.1.2 outputs

#load packages
library(Seurat)
library(ggplot2)
library(RColorBrewer)

#download RAW data from GEO:
data.dir="../GSE190815_RAW"

#Load sample information file
info=read.csv("../data/SampleInformation.csv")

#Load pre=processing functions
source('../code/PreProcessing_fxns.R', chdir = TRUE)

#create list for each sample
data=vector(mode="list",length=length(info$GEO_ID))

#Obtain putative human nuclei and process data using seurat workflow for each 10x sample

for (i in c(1:length(data))){
	print(paste0("Processing ",info$GEO_ID[i]))

	#temporarily change file name
	files=paste0(data.dir,"/",info$GEO_ID[i],"_",info$PlotName[i],c("_barcodes.tsv.gz","_matrix.mtx.gz","_features.tsv.gz"))
	filesRename= paste0(data.dir,"/",c("barcodes.tsv.gz","matrix.mtx.gz","features.tsv.gz"))
	file.rename(from = files, to = filesRename)

	#get seurat object
	data[[i]]= getSeurat(bc.dir=data.dir,min.cells=0,min.features=0,proj.name=info$PlotName[i])

	#revert file names back and obtain meta data on human and rat counts per cell
	file.rename(from = filesRename, to = files)

	#Add ribosomal gene meta data
	ribo.genes = rownames(data[[i]]@assays$RNA@counts)[grep("^RP[SL]",rownames(data[[i]]@assays$RNA@counts))]
	percent.ribo <- colSums(data[[i]]@assays$RNA@counts[is.element(rownames(data[[i]]@assays$RNA@counts), ribo.genes),])/Matrix::colSums(data[[i]]@assays$RNA@counts)*100
	data[[i]] <- AddMetaData(data[[i]], percent.ribo, col.name = "Hs.percent.ribo")
	
	#Filter putative rat cell and low quality nuclei
	data[[i]]=QC.filter(seu.object=data[[i]],min.human.percent=95,mito.percent=20,min.genes=1000, max.genes=100000,min.cells=5)	

	#Further filter putative doublets, and low quality nuclei, as described in methods using iterative approach

	#load previously filtered meta data for high quality cells
	metaFinal=read.csv(paste0(data.dir,"/",info$GEO_ID[i],"_",info$PlotName[i],"_cellMetaData.csv.gz"))
	metaFinal$barcode=sapply(strsplit(metaFinal$barcode,"_"),function(x){x[[1]]})
	
	#merge with seurat object
	meta=merge(data.frame(barcode=rownames(data[[i]]@meta.data),data[[i]]@meta.data,Order=c(1:dim(data[[i]]@meta.data)[1])), metaFinal[,c(1,grep("seurat_clusters",colnames(metaFinal)),grep("subclass_label",colnames(metaFinal)),grep("cluster_label",colnames(metaFinal)))],by.x=1,by.y=1,all.x=T)
	meta=meta[order(meta$Order,decreasing=F),]
	meta$Filter="keep"
	meta$Filter[is.na(meta$subclass_label)]="out"
	
	if (all.equal(colnames(data[[i]]),meta$barcode)){
		data[[i]]@meta.data$Filter=meta$Filter
		data[[i]]@meta.data$seurat_clusters =meta$seurat_clusters
		data[[i]]@meta.data$subclass_label =meta$subclass_label
		data[[i]]@meta.data$cluster_label =meta$cluster_label
	}else {stop("meta data file does not match seurat object")}
	
	#filter seurat object to high quality cells
	
	Idents(data[[i]])="Filter"
	data[[i]]=subset(data[[i]],idents="keep")

	data[[i]] =SCTransform(data[[i]], vst.flavor = "v2",verbose = T)
	data[[i]] = RunPCA(data[[i]],npcs=100)

	data[[i]] = RunUMAP(data[[i]], dims = 1:30)
	data[[i]] = FindNeighbors(data[[i]], dims = 1:30)


}

#Integrate t-hCS data

thcs=list(data[[1]],data[[2]],data[[4]])

features <- SelectIntegrationFeatures(object.list = thcs, nfeatures = 3000)
thcs <- PrepSCTIntegration(object.list = thcs, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = thcs, normalization.method = "SCT",anchor.features = features)

data.combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims=1:30)
data.combined <- RunPCA(data.combined,npcs=100)
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:30, verbose = FALSE)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:30)
data.combined <- PrepSCTFindMarkers(data.combined)

#thCS
thcs= FindClusters(data.combined, resolution = 0.5)
DefaultAssay(thcs) <- "SCT"

#Integrate hCS data

hcs=list(data[[3]],data[[5]],data[[6]])

features <- SelectIntegrationFeatures(object.list = hcs, nfeatures = 3000)
hcs <- PrepSCTIntegration(object.list = hcs, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = hcs, normalization.method = "SCT",anchor.features = features)

data.combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims=1:30)
data.combined <- RunPCA(data.combined,npcs=100)
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:30, verbose = FALSE)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:30)
data.combined <- PrepSCTFindMarkers(data.combined)

#hCS
hcs= FindClusters(data.combined, resolution = 0.5)
DefaultAssay(hcs) <- "SCT"

