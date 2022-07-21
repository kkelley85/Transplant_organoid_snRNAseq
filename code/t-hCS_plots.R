#Code to generate t-hCS plots, Fig.1 and Extended Data Fig. 4

library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

#load processed t-hCS .rds data
#download from here: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE190815&format=file&file=GSE190815%5FhCS%5Fprocessed%5FSeuratObject%2Erds%2Egz
clust=readRDS("GSE190815_t-hCS_processed_SeuratObject.rds")

#Set colors
colors=c("#C7E9B4","#FCC5C0","#253494","#225EA8","#41AB5D","#AE017E","#F768A1","#7A0177","#DD3497","#FDE0DD","#FA9FB5","#238443","#E7E1EF","#1D91C0","#FFF7F3","#EDF8B1","#FD8D3C")

cellLabel=unique(as.character(clust@meta.data$cluster_label))

names(colors)=cellLabel
		
#Generate UMAP plot: Figure 1g
p1= DimPlot(clust, reduction = "umap",group.by="cluster_label",shuffle=T,cols=colors, label = TRUE, repel = TRUE)&NoLegend()&NoAxes()& ggtitle(paste0(""))&theme(text = element_text(size=6))

setwd("../output")
pdf(file=paste0("t-hCS_All_clusters_res0.5_Integrate_UMAP.pdf"),width=3,height=3)
print(p1)
dev.off()

#Generate Integration by sample plot: Extended data 4
colors= brewer.pal(8, "Set1")[1:3]
clust$Plot.name=factor(as.character(clust$Plot.name),levels=c("t-hCS_d224_2242","t-hCS_d227_2242","t-hCS_d276_Q3"))

p1= DimPlot(clust, reduction = "umap",group.by="Plot.name",shuffle=T,cols=colors, label = F, repel = TRUE)&NoLegend()&NoAxes()& ggtitle(paste0(""))&theme(text = element_text(size=6))

pdf(file=paste0("t-hCS_bySample_Integration_UMAP.pdf"),width=3,height=3)
print(p1)
dev.off()

#Marker genes: supplementary table 3
Idents(clust)="cluster_label"
markers.out=FindAllMarkers(clust, only.pos=T)

#CellProportions: Extended Data Figure 4

cellLabel=unique(as.character(clust@meta.data$cluster_label))

colors=c("#C7E9B4","#FCC5C0","#253494","#225EA8","#41AB5D","#AE017E","#F768A1","#7A0177","#DD3497","#FDE0DD","#FA9FB5","#238443","#E7E1EF","#1D91C0","#FFF7F3","#EDF8B1","#FD8D3C")

names(colors)=cellLabel

cluster_order=c("Cyc_Prog_c15","Progenitor_c3","Astroglia_c9","Astroglia_c4","Astroglia_c5","OPC_c1","Oligo_c14","GluN_UL_c13","GluN_UL_c2","GluN_UL_c0","GluN_UL_c7","GluN_UL_c6","GluN_DL_c11","GluN_DL_c8","GluN_DL/SP_c10","GluN_c12","RELN_c16")

clust$cluster_label =factor(as.character(clust$cluster_label),levels=cluster_order)
		
df= clust@meta.data

df2=df %>%
	group_by(cluster_label,Plot.name) %>%
	summarize(n())

colnames(df2)[3]="percentage"	

df2$CellNames2=factor(as.character(df2$cluster_label),levels= cluster_order)

pdf(file=paste0("t-hCS_All_SeuratClusters_CellProportions.pdf"),width=1.65,height=2)
ggplot(df2, aes(fill= CellNames2,x=Plot.name, y=percentage)) + 
	geom_bar( stat="identity", position="fill") +
				scale_fill_manual(values=colors[levels(df2$CellNames2)]) +
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



