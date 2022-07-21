#Code used to generate QC Plots: Extended Data Fig. 4

library(Seurat)
library(ggplot2)
library(RColorBrewer)

#download RAW data from GEO:
#https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE190815&format=file
data.dir="../GSE190815_RAW"

#Load sample information file
info=read.csv("../data/SampleInformation.csv")

#Load pre=processing functions
source('../code/PreProcessing_fxns.R', chdir = TRUE)

#create list for each sample
data=vector(mode="list",length=length(info$GEO_ID))

#Loop through samples to obtain seurat objects from 10x Cellranger outputs


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

	#Obtain rat_vs_human_counts plots: Extended Data Figure 4
	
	rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
	r <- rf(32)
	tmp=data[[i]]@meta.data
	
	tmp$human.counts=tmp$human.counts+1
	tmp$rat.counts=tmp$rat.counts+1
	tmp$total.counts=tmp$human.counts+tmp$rat.counts

	setwd("../output")
	pdf(file=paste0(info$PlotName[i],"_cells_QC_rat_vs_human_counts.pdf"),width=1.25,height=1.15)
		print(ggplot(data=tmp,aes(human.counts, rat.counts)) +
		stat_bin2d(bins=100) + 
		scale_fill_gradientn(colours=r,trans="log10") +
		scale_x_continuous("human counts",trans='log10',limits=c(100,max(tmp$total.counts))) +
  		scale_y_continuous("rat counts",trans='log10', limits=c(5, 10000)) +
		theme_bw() +
		labs(fill = "cells")+
		  theme(
    		axis.text.x = element_text(size=6,angle=45,hjust=1,color="black"),
    		axis.text.y = element_text(size=6,color = "black"),
    		axis.title=element_text(size=6),
			axis.ticks.length=unit(0.01,"inch"),
			#panel.grid.major= element_blank(),
			panel.grid.major= element_line(size=0.5),
			panel.grid.minor= element_blank(),
    		legend.text=element_text(size=6),
    		legend.title= element_text(size=6),
    		legend.position="",
    		panel.border = element_rect(colour="black", size=0.5),
    		plot.background=element_rect(fill="transparent",colour=NA)
  		))
	dev.off()
	

	#Obtain putative human nuclei
	data[[i]]=QC.filter(seu.object=data[[i]],min.human.percent=95,mito.percent=100,min.genes=0, max.genes=100000,min.cells=0)
	data[[i]][["Plot.name"]]=info$PlotName[i]

	#construct plotting data frame
	
	if (i==1){
		df=data[[i]]@meta.data
	} else {
		
		df=rbind(df, data[[i]]@meta.data)
	}
	

}

#Make plots: Extended Data Fig 4

df$Condition=sapply(strsplit(df$Plot.name,"_"),function(x){x[[1]]})

df$Condition=factor(as.character(df$Condition),levels=c("t-hCS","hCS"))
df$Plot.name=factor(as.character(df$Plot.name),levels=c("t-hCS_d224_2242","t-hCS_d227_2242","t-hCS_d276_Q3","hCS_d227_2242","hCS_d243_Q2","hCS_d257_Q3"))

	genes=ggplot(df, aes(y= nFeature_RNA, x= Plot.name,fill=Condition)) + 
		geom_violin(color="black",size=0.25) +	
		#geom_boxplot(alpha=0.5,notch=T,outlier.colour=NA) +
				scale_fill_manual(values=c("#007238","grey70"))+
				scale_y_continuous(trans="log10",name="# genes",limits=c(300,10100))+
				scale_x_discrete(name="")+
				theme_bw()+
				geom_hline(yintercept=1000, color = "black", size=0.25)+
				theme(
					axis.title.y=element_text(color="black",size=6,face="plain"),
					axis.title.x=element_text(color="black", size=6,face="plain"),
					axis.text.y=element_text(color="black", size=6),
					axis.text.x=element_text(color="black", size=6, hjust=1,angle=45),
					plot.title=element_text(color="black", size=12),
					legend.position="",
					legend.text=element_text(size=6),
					legend.title=element_text(size=0),
					panel.grid.major = element_blank(),
					panel.grid.minor = element_blank(),
					axis.ticks=element_line(size=0.1),
					axis.line = element_line(colour = "black",size=0.25),
					panel.background = element_blank())

	counts=ggplot(df, aes(y= nCount_RNA, x= Plot.name,fill= Condition)) + 
		geom_violin(color="black",size=0.25) +	
		#geom_boxplot(alpha=0.5,notch=T,outlier.colour=NA) +
				scale_fill_manual(values=c("#007238","grey70"))+
				scale_y_continuous(trans="log10",name="# counts",limits=c(300,54000))+
				scale_x_discrete(name="")+
				theme_bw()+
				#geom_hline(yintercept=1000, color = "black", size=0.5)+
				theme(
					axis.title.y=element_text(color="black",size=6,face="plain"),
					axis.title.x=element_text(color="black", size=6,face="plain"),
					axis.text.y=element_text(color="black", size=6),
					axis.text.x=element_text(color="black", size=6, hjust=1,angle=45),
					plot.title=element_text(color="black", size=12),
					legend.position="",
					legend.text=element_text(size=6),
					legend.title=element_text(size=0),
					panel.grid.major = element_blank(),
					panel.grid.minor = element_blank(),
					axis.ticks=element_line(size=0.1),
					axis.line = element_line(colour = "black",size=0.25),
					panel.background = element_blank())


		df2 <-
  		df %>%
  		group_by(Plot.name) %>%
  		mutate(outlier = percent.human.mt > median(percent.human.mt) + IQR(percent.human.mt) * 1.5) %>%
  		ungroup

	humanmito=ggplot(df2, aes(y= percent.human.mt, x= Plot.name, fill= Condition)) + 
#	geom_violin(color="black",fill="grey80") +	
	geom_boxplot(alpha=0.5,notch=F, outlier.colour=NA ,color="black",size=0.25) +
	geom_point(data = function(x) dplyr::filter_(x, ~ outlier), alpha=0.25,position = 'jitter',size=0.01)+
				scale_fill_manual(values=c("#007238","grey70"))+
				scale_y_continuous(name="Percent MT counts",limits=c(0,50))+
				scale_x_discrete(name="")+
				theme_bw()+
				geom_hline(yintercept=20, color = "black", size=0.25)+
				theme(
					axis.title.y=element_text(color="black",size=6,face="plain"),
					axis.title.x=element_text(color="black", size=6,face="plain"),
					axis.text.y=element_text(color="black", size=6),
					axis.text.x=element_text(color="black", size=6, hjust=1,angle=45),
					plot.title=element_text(color="black", size=12),
					legend.position="",
					legend.text=element_text(size=6),
					legend.title=element_text(size=0),
					panel.grid.major = element_blank(),
					panel.grid.minor = element_blank(),
					axis.ticks=element_line(size=0.1),
					axis.line = element_line(colour = "black",size=0.25),
					panel.background = element_blank())


	pdf(file=paste0("QCplots_unfiltered.pdf"),width=2,height= 4.5)
	print(counts/genes/humanmito)
	dev.off()

