#code to generate Fig. 2 activity-dependent transcription plots

library(Seurat)
library(ggplot2)
library(Libra)
library(edgeR)

#Load processed data
thcs=readRDS("GSE190815_t-hCS_processed_SeuratObject.rds")
hcs=readRDS("GSE190815_hCS_processed_SeuratObject.rds")

#load GluN DE genes generated from code/GluN_DE_analysis.R 

glun_up=read.csv("../output/Up_thCS_All_GluN_fold_2_min.percent_10.csv")
glun_all=read.csv("../output/thCS_All_GluN_DEanalysis_Libra.csv")

#load activity-dependent gene sets
#mouse response genes from Hrvatin et al. Nature Neuroscience 2018 Table S3 
#human enriched response genes from Ataman et al. Nature 2016 Table S4

lateExc=read.csv("../data/Hrvatin_lateResponse_excitatoryNeuron_genes.csv",header=F)
earlyExc=read.csv("../data/Hrvatin_earlyResponse_excitatoryNeuron_genes.csv",header=F)
humanOnly=read.csv("../data/Ataman_Human_lateResponse_genes.csv",header=F)

#Function to perform hypergeometic overlap test
overlapTest=function(GeneSet,DEgenes,AllGenes){
	shared=length(intersect(AllGenes, GeneSet))
	shared.in.DEgenes=length(intersect(DEgenes, GeneSet))
	shared.out.DEgenes= shared-shared.in.DEgenes
	DEgenes.not.shared=length(DEgenes)-shared.in.DEgenes
	out.DEgenes.not.shared=length(AllGenes)-length(DEgenes)-shared.out.DEgenes
	fisher.test(matrix(c(shared.in.DEgenes, DEgenes.not.shared, shared.out.DEgenes, out.DEgenes.not.shared),ncol=2),alternative="greater")$p.val
	}

#Perform significance testing
fisher.out=data.frame(Name=c("Ms early-activity","Ms late-activity","Hs late-activity"),Pvalue=c(overlapTest(earlyExc[,1], glun_up $gene, glun_all$gene), overlapTest(lateExc[,1],glun_up$gene, glun_all$gene), overlapTest(humanOnly[,1],glun_up$gene, glun_all $gene)))

#Plot p-values: Figure 2

df=fisher.out

df$Log=-log10(df$Pvalue)
FDR=-log10((0.05/3))

df$Name =factor(as.character(df$Name),levels=rev(as.character(df$Name)))

col= brewer.pal(9, "Greys")[c(3,5,7)]


setwd("../output")
pdf(file="UpNeuron_Activity_Enrichment.pdf",width=1.1,height=1)

ggplot(df, aes(Name,fill=Name)) +
geom_bar(aes(Name,Log), width=.6, stat = "identity")+
coord_flip() +
ylab("-log10(p-value)") +
scale_y_continuous(breaks=c(0,5,10))+
geom_hline(yintercept = FDR,color="black",size=0.25) +
theme_bw() +
scale_fill_manual(values =col) +
theme(axis.line.x = element_line(colour = "black",),
axis.line.y = element_line(colour = "black",),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(),
axis.text.y = element_text(size=6,colour="black"),
axis.text.x = element_text(size=6,colour="black"),
axis.title.x = element_text(size =6,vjust=+0.1),
axis.title.y = element_blank(),
axis.ticks.x = element_line(size=0.25),
axis.ticks.y = element_line(size=0.25),
plot.background=element_rect(fill = "transparent",colour = NA),
legend.position="",
legend.title= element_text(size=0,colour="black"),
legend.text= element_text(size=6,colour="black"),
axis.line= element_line(colour = 'black', size = 0.25),
legend.key.height=unit(0.25,"cm"),
legend.key.width=unit(0.25,"cm"))

dev.off()

#Generate heatmaps of significant genes
#use to_pseudobulk() function from Libra Package

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

#use to_pseudobulk() function from Libra Package
counts= to_pseudobulk(comb.glun, meta = comb.glun@meta.data,replicate_col="replicate",cell_type_col="subclass_broad")

#generate counts per million using edgeR cpm() function
counts_edge=DGEList(counts= counts$GluN)
counts_edge =calcNormFactors(counts_edge, method = "TMM")

CPM=cpm(counts_edge,log=F)

colnames(CPM)=gsub(":hCS","",colnames(CPM))
colnames(CPM)=gsub(":t-hCS","",colnames(CPM))

GluN_CPM= CPM

#Get significant up-regulated t-hCS GluN activity genes
ms_ERG= glun_up[is.element(glun_up$gene,earlyExc[,1]),]
ms_LRG= glun_up[is.element(glun_up$gene, lateExc[,1]),]
hs_LRG= glun_up[is.element(glun_up$gene, humanOnly[,1]),]

#scale data and make plots for each category

#early-response genes

library(reshape2)
tmp=GluN_CPM[is.element(rownames(GluN_CPM), ms_ERG$gene),]
tmp =t(scale(t(tmp)))

df=melt(tmp)
df$Var1=factor(as.character(df$Var1),levels=(ms_ERG$gene))
df$Var2=factor(as.character(df$Var2),levels=rev(c("t-hCS_d224_2242","t-hCS_d227_2242","t-hCS_d276_Q3","hCS_d227_2242","hCS_d243_Q2","hCS_d257_Q3")))

myPalette=brewer.pal(11,"RdBu")

p1=ggplot(df,aes(x= Var1,y=Var2,fill=value))+
geom_tile(color="white",size=0.05)+
scale_fill_gradientn(colours=rev(myPalette),limits=c(-1.5,1.75))+
theme_bw()+
xlab("") +
ylab("") +
ggtitle("") +
theme(
	axis.text.x=element_text(size=5,angle=45,hjust=1),
	axis.text.y=element_text(size=5,),
	title= element_text(size=5),
	axis.ticks=element_line(size=0.05),
	axis.ticks.length=unit(0.01,"inch"),
	panel.grid.major= element_blank(),
	panel.grid.minor= element_blank(),
	legend.text= element_text(size=6),
	panel.border = element_rect(colour = "black", fill=NA, size=0.25),
	legend.position="",
	plot.background=element_rect(fill="transparent",colour=NA))


#late-response genes
tmp=GluN_CPM[is.element(rownames(GluN_CPM), ms_LRG$gene),]
tmp =t(scale(t(tmp)))

df=melt(tmp)
df$Var1=factor(as.character(df$Var1),levels=(ms_LRG$gene))
df$Var2=factor(as.character(df$Var2),levels=rev(c("t-hCS_d224_2242","t-hCS_d227_2242","t-hCS_d276_Q3","hCS_d227_2242","hCS_d243_Q2","hCS_d257_Q3")))

myPalette=brewer.pal(11,"RdBu")

p2=ggplot(df,aes(x= Var1,y=Var2,fill=value))+
geom_tile(color="white",size=0.05)+
scale_fill_gradientn(colours=rev(myPalette),limits=c(-1.5,1.75))+
theme_bw()+
xlab("") +
ylab("") +
ggtitle("") +
theme(
	axis.text.x=element_text(size=5,angle=45,hjust=1),
	axis.text.y=element_text(size=5,),
	title= element_text(size=10),
	axis.ticks=element_line(size=0.05),
	axis.ticks.length=unit(0.01,"inch"),
	panel.grid.major= element_blank(),
	panel.grid.minor= element_blank(),
	legend.text= element_text(size=6),
	panel.border = element_rect(colour = "black", fill=NA, size=0.25),
	legend.position="",
	plot.background=element_rect(fill="transparent",colour=NA))


#Human LRG
tmp=(GluN_CPM[is.element(rownames(GluN_CPM), hs_LRG$gene),])

tmp =t(scale(t(tmp)))

df=melt(tmp)
df$Var1=factor(as.character(df$Var1),levels=(hs_LRG$gene))
df$Var2=factor(as.character(df$Var2),levels=rev(c("t-hCS_d224_2242","t-hCS_d227_2242","t-hCS_d276_Q3","hCS_d227_2242","hCS_d243_Q2","hCS_d257_Q3")))

myPalette=brewer.pal(11,"RdBu")

p3=ggplot(df,aes(x= Var1,y=Var2,fill=value))+
geom_tile(color="white",size=0.05)+
scale_fill_gradientn(colours=rev(myPalette),limits=c(-1.5,1.75))+
theme_bw()+
xlab("") +
ylab("") +
ggtitle("") +
theme(
	axis.text.x=element_text(size=5,angle=45,hjust=1),
	axis.text.y=element_text(size=5,),
	title= element_text(size=10),
	axis.ticks=element_line(size=0.05),
	axis.ticks.length=unit(0.01,"inch"),
	panel.grid.major= element_blank(),
	panel.grid.minor= element_blank(),
	legend.text= element_text(size=6),
	panel.border = element_rect(colour = "black", fill=NA, size=0.25),
	legend.position="",
	plot.background=element_rect(fill="transparent",colour=NA))


setwd("../output")

pdf(file=paste0("Ms_ERG_heatmap",".pdf"),width=1.75,height=1.15)
print(p1)
dev.off()

pdf(file=paste0("Ms_LRG_heatmap",".pdf"),width=2.6,height=1.25)
print(p2)
dev.off()

pdf(file=paste0("Hs_LRG_heatmap",".pdf"),width=1.85,height=1.275)
print(p3)
dev.off()

