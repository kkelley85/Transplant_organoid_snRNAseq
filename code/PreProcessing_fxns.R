# Functions for analyzing Transplant 10x single-nuclei RNAseq CellRanger Count v6.1.2 outputs aligned to combine human/rat reference


#Function to read in 10x data and add meta.data
#######################
getSeurat=function(
				bc.dir=NULL,
				min.cells=0,
				min.features=0,
				proj.name=NULL
				){
							
	data <- Read10X(data.dir = bc.dir)
	full <- CreateSeuratObject(counts = data, project = proj.name, min.cells = min.cells, min.features = min.features)
	human.genes=grep("GRCh38-2020A-v5-",rownames(full))
	rat.genes=grep("Rnor-6.0--------",rownames(full))

	human.counts=Matrix::colSums(full[human.genes,])
	rat.counts=Matrix::colSums(full[rat.genes,])
	total.counts=Matrix::colSums(full)
	percent.human=100*(human.counts/(total.counts))

	barcode.info=data.frame(name=colnames(full),human.counts=human.counts,rat.counts=rat.counts,total.counts=total.counts,percent.human=percent.human)


	full[["percent.human.mt"]] <- PercentageFeatureSet(full, pattern = "GRCh38-2020A-v5-MT-")
	full[["percent.rat.mt"]] <- PercentageFeatureSet(full, pattern = "Rnor-6.0--------Mt-")

	if (all.equal(colnames(full),rownames(barcode.info))){
		full[["percent.human"]] <- barcode.info$percent.human
		full[["human.counts"]] <- barcode.info$human.counts
		full[["rat.counts"]] <- barcode.info$rat.counts

	} else {stop("cell IDs do not match")}

#Get unique human geenes
    features <- grep(pattern ="GRCh38-2020A-v5-", x = rownames(x = full))
	human.genes=apply(GetAssayData(object = full, assay = NULL, slot = "counts")[features,],2,function(x){sum(x>0)})
	full[["human.genes"]] <- human.genes

#Remove rat genes

	counts=full@assays$RNA@counts
	counts=counts[grep("GRCh38-2020A-v5-",rownames(counts)),]
	rownames(counts)=gsub("GRCh38-2020A-v5-","",rownames(counts))
	meta= full@meta.data

	full <- CreateSeuratObject(counts = counts, project = proj.name, meta.data=meta)

	return(full)

}
#######################


##Quality Control function

QC.filter=function(
				seu.object=NULL,
				min.human.percent=NULL,
				mito.percent=NULL,
				min.genes=NULL,
				max.genes=NULL,
				min.cells=NULL
){
	
	print("Initial dimensions")
	print(paste0("Number of cells: ", dim(seu.object)[2]))
	print(paste0("Number of features: ", dim(seu.object)[1]))

## Filter out rat cells
	seu.object =subset(seu.object,subset= (percent.human >= min.human.percent))
	print("Remove rat cells")
	print(paste0("Number of cells: ", dim(seu.object)[2]))
	print(paste0("Number of features: ", dim(seu.object)[1]))
	
# Filter low quality
	seu.object =subset(seu.object,subset= (human.genes >= min.genes)&(human.genes<= max.genes)&(percent.human.mt<=mito.percent))
	print("Remove low quality cells")
	print(paste0("Number of cells: ", dim(seu.object)[2]))
	print(paste0("Number of features: ", dim(seu.object)[1]))

#subset to features present in at least X cells
	min.cells.gene=apply(GetAssayData(object = seu.object, assay = NULL, slot = "counts"),1,function(x){sum(x>0)})

	seu.object = seu.object[min.cells.gene>=min.cells,]
	print("Subset features to min cells")
	print(paste0("Number of cells: ", dim(seu.object)[2]))
	print(paste0("Number of features: ", dim(seu.object)[1]))

	return(seu.object)	
}
######################

