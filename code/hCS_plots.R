#Code to generate hCS plots: Extended Data Fig. 4

library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

#load processed hCS .rds data
#download from here: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE190815&format=file&file=GSE190815%5FhCS%5Fprocessed%5FSeuratObject%2Erds%2Egz
clust=readRDS("GSE190815_hCS_processed_SeuratObject.rds")

#Set colors
colors=c("#D53E4F","#2171B5","#6BAED6","#DD3497","#7A0177","#4292C6","#8DD3C7","#9E0142","#FA9FB5","#35978F","#08519C","#C7E9B4","#08306B","#01665E","#FD8D3C","#EDF8B1","#B3DE69")

cellLabel=unique(as.character(clust@meta.data$cluster_label))

names(colors)=cellLabel

#Generate UMAP plot: Extended Data Fig. 4

p1= DimPlot(clust, reduction = "umap",group.by="cluster_label",shuffle=T,cols=colors, label = TRUE, repel = TRUE)&NoLegend()&NoAxes()& ggtitle(paste0(""))&theme(text = element_text(size=6))


setwd("../output")
pdf(file=paste0("hCS_All_clusters_res0.5_Integrate_UMAP.pdf"),width=3,height=3)
print(p1)
dev.off()

#Generate Integration by sample plot: Extended data Fig. 4

colors= brewer.pal(8, "Set1")[1:3]

clust$Plot.name=factor(as.character(clust$Plot.name),levels=c("hCS_d227_2242-1","hCS_d243_1208-2","hCS_d257_8119-1"))

p1= DimPlot(clust, reduction = "umap",group.by="Plot.name",shuffle=T,cols=colors, label = F, repel = TRUE)&NoLegend()&NoAxes()& ggtitle(paste0(""))&theme(text = element_text(size=6))

pdf(file=paste0("hCS_bySample_Integration_UMAP.pdf"),width=3,height=3)
print(p1)
dev.off()


#CellProportions: Extended Data Figure 4

cellLabel=unique(as.character(clust@meta.data$cluster_label))
colors=c("#D53E4F","#2171B5","#6BAED6","#DD3497","#7A0177","#4292C6","#8DD3C7","#9E0142","#FA9FB5","#35978F","#08519C","#C7E9B4","#08306B","#01665E","#FD8D3C","#EDF8B1","#B3DE69")

names(colors)=cellLabel

cluster_order=c("Cyc_Prog_c9","Progenitor_c2","Astroglia_c3","Astroglia_c8","Astroglia_c1","Astroglia_c0","Astroglia_c6","IPC_c13","GluN_UL_c7","GluN_DL_c5","GluN_DL/SP_c10","RELN_c15","IN_c12","IN_c4","Choroid_c14","Mening._c11","Mening._c16")

clust$cluster_label =factor(as.character(clust$cluster_label),levels=cluster_order)
		
df= clust@meta.data

df2=df %>%
	group_by(cluster_label,Plot.name) %>%
	summarize(n())

colnames(df2)[3]="percentage"	

df2$cluster_label =factor(as.character(df2$cluster_label),levels= cluster_order)

pdf(file=paste0("hCS_All_SeuratClusters_CellProportions.pdf"),width=1.65,height=2)
	ggplot(df2, aes(fill= cluster_label,x=Plot.name, y=percentage)) + 
	geom_bar( stat="identity", position="fill") +
				scale_fill_manual(values=colors[levels(df2$cluster_label)]) +
				ylab("proportion") +
				scale_y_continuous(expand = expansion(mult = c(0, .1)))+
				theme_bw()+
				theme(
					axis.line.x = element_line(colour = "black",),
					axis.line.y = element_line(colour = "black",),
					panel.grid.major = element_blank(),
					panel.grid.minor = element_blank(),
					axis.title.y=element_text(size=5,face="plain"),
					axis.title.x=element_text(size=0,face="plain"),
					axis.text.y=element_text(size=5),
					axis.text.x=element_text(size=5,angle=45,hjust=1),
					plot.title=element_text(size=12),
					legend.position="right",
					legend.text=element_text(size=5),
					legend.title=element_text(size=0),
					legend.key.height=unit(0.15,"cm"),
					legend.key.width=unit(0.15,"cm"),
					panel.border = element_blank(),
					panel.background = element_blank())
dev.off()

#Marker genes: supplementary table 4
Idents(clust)="cluster_label"
markers.out=FindAllMarkers(clust,only.pos=T)

#write.table(markers.out,file="hCS_All_FindMarkers_seuratClusters_Positive.csv",sep=",",col.names=T,row.names=F)


