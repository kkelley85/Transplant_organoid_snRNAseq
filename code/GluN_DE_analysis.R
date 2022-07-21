#Glutamatergic neuron DE analysis using Libra package: https://github.com/neurorestore/Libra

library(Libra)
library(Seurat)
library(ggplot2)

#Load processed data

thcs=readRDS("GSE190815_t-hCS_processed_SeuratObject.rds")
hcs=readRDS("GSE190815_hCS_processed_SeuratObject.rds")

#Setup meta data appropriately to be handled by Libra

thcs[["cell_type"]]= thcs@meta.data$subclass_label
hcs[["cell_type"]]= hcs@meta.data$subclass_label

thcs[["replicate"]]= thcs@meta.data$Plot.name
hcs[["replicate"]]= hcs@meta.data$Plot.name

thcs[["label"]]="t-hCS"
hcs[["label"]]="hCS"

#Combine data into single object
comb=merge(thcs,hcs)

##subset to glutamatergic neurons

subclass=comb$cell_type
subclass[is.element(subclass,c("GluN_DL","GluN","GluN_DL/SP","GluN_UL"))]="GluN"
comb[["subclass_broad"]]=subclass

Idents(comb)="subclass_broad"
comb.glun=subset(comb,idents="GluN")

#Run Libra DE analysis
DE.glun = run_de(comb.glun, meta = comb.glun@meta.data,cell_type_col="subclass_broad")

#subset up-regulated t-hCS DE genes to those expressed in at least 10% of cells t-hCS GluN cells
Idents(thcs)="subclass_label"
glun=subset(thcs,idents=c("GluN_UL","GluN","GluN_DL","GluN_DL/SP"))
glun= SplitObject(glun, split.by = "Plot.name")

glun_proportion=matrix(NA,nrow=dim(glun[[1]])[1],ncol=3)
rownames(glun_proportion)=rownames(glun[[1]])

for (i in c(1:length(glun))){
	
	tmp=glun[[i]]@assays$RNA@counts
	tmp=apply(tmp,1,function(x){sum(x>0)})
	tmp=tmp/dim(glun[[i]])[2]
	tmp=tmp[is.element(names(tmp),rownames(glun[[i]]))]
	tmp=tmp[rownames(glun_proportion)]
	glun_proportion[,i]=tmp
	
}

min.expr=apply(glun_proportion,1,function(x){min(x)})

#filter genes to 2-fold up-regulated in t-hCS, signifiant (adjusted-pvalue <0.05) & expressed in at least 10 percent of cells in each t-hCS smaple
neuron.de=data.frame(DE.glun)
up_neuron=neuron.de[(neuron.de$avg_logFC<=-1)&(neuron.de$p_val_adj<0.05),]
up_neuron=up_neuron[is.element(up_neuron$gene,names(min.expr[min.expr>=0.1])),]

setwd("../output")
write.table(up_neuron,"Up_thCS_All_GluN_fold_2_min.percent_10.csv",sep=",",col.names=T,row.names=F)
write.table(DE.glun,"thCS_All_GluN_DEanalysis_Libra.csv",sep=",",col.names=T,row.names=F)


