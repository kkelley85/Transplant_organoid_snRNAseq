#Code to generate psychEncode/BrainSpan plots Extended Fig. 5

library(ggplot2)
library(dplyr)
library(patchwork)

#download psychENCODE processed data: 
#expression: http://evolution.psychencode.org/files/processed_data/RNA-seq/nhp_development_RPKM_rmTechRep.txt
#sample info: http://evolution.psychencode.org/files/processed_data/RNA-seq/mRNA-seq_QC.xlsx

#load psychENCODE data
phe <- read.csv("SampleInfo.csv")
rawRPKM <- read.table("nhp_development_RPKM_rmTechRep.txt", as.is=TRUE, header=TRUE)

tmp=merge(phe,data.frame(colnames(rawRPKM),Order=c(1:dim(rawRPKM)[2])),by.x=1,by.y=1)
tmp=tmp[order(tmp$Order,decreasing=F),]

phe=tmp[,c(1:19)]
all.equal(colnames(rawRPKM),phe[,1])

ncxRegOrder <- c("MFC","OFC","DFC","VFC","M1C","S1C","IPC","A1C","STC","ITC","V1C") # "HIP","AMY","STR","MD","CBC"

# subset to neocortical human samples
# subset to Period>=3 as done in original publication
hNCXExpr <- rawRPKM[,(phe$Species == "Human") & (phe$Region %in% ncxRegOrder) & (phe$Period >= 3)]

hNCXPhe <- phe[phe$Species == "Human" & phe$Region %in% ncxRegOrder & phe$Period >= 3,]
all.equal(colnames(hNCXExpr), hNCXPhe$Sample)

#load developmental genes
dev_genes=read.csv("../data/PsychEncode_development_genes.csv")

#subset expression to developmental genes
hNCXSelExpr <- hNCXExpr[is.element(rownames(hNCXExpr),dev_genes[,1]),]
colnames(hNCXSelExpr)= c(paste0("P",hNCXPhe$Period,"_", colnames(hNCXSelExpr)))

#use gene symbols
hNCXSelExpr2= hNCXSelExpr
rownames(hNCXSelExpr2)= sapply(strsplit(rownames(hNCXSelExpr),"|",fixed=T),function(x){unlist(x)[2]})

#load pseudobulk t-hCS/hCS data
hcs=read.csv("../data/All_pseudobulk_RPKM.csv")
rownames(hcs)=hcs[,1]
hcs=hcs[,-c(1)]
hcs=log2(hcs+1)
colnames(hcs)=gsub("t.hCS","t-hCS",colnames(hcs))

#combine data
data= merge(hcs, hNCXSelExpr2,by="row.names")
rownames(data)= data[,1]
data= data[,2:dim(data)[2]]

#principle component analysis
sample_pca=prcomp(t(data),scale=T)
pca= sample_pca$x

pc_scores <- pca %>%as_tibble(rownames = "sample")
pc_scores=data.frame(pc_scores)

pc_scores$Group=sapply(strsplit(pc_scores$sample,"_"),function(x){x[[1]]})

#Setup data for plotting
#psychENCODE data
df= pc_scores[!is.element(pc_scores$Group,c("hCS","t-hCS")),]
df$sample_name=sapply(strsplit(df$sample,"_"),function(x){x[[2]]})
df=merge(df, hNCXPhe[,c(1,8)],by.x="sample_name",by.y=1)
df$Log2Days=log2(df$Days)

#hCS/t-hCS data
df2= pc_scores[is.element(pc_scores$Group,c("hCS","t-hCS")),]
df2$Days=c(227,243,257,224,227,276)
df2$Log2Days=log2(df2$Days)

#combine
df=rbind(df[,2:dim(df)[2]],df2)
df$Group2="Brainspan"
df$Group2[is.element(df$Group,"hCS")]="hCS"
df$Group2[is.element(df$Group,"t-hCS")]="t-hCS"

#plot

#define developmental periods as done by Kang et al. Nature 2011
bound=c(0,c(13*7,16*7,19*7,24*7,38*7,266+6*30,266+365,266+365*6,266+365*12,266+365*20,266+365*40,266+365*60)-1)

#psychENCODE plot of PC1

p1=ggplot() + 
	geom_point(data= df[df$Group2=="Brainspan",], aes(x=Days, y=PC1),color="black",size=0.25,alpha=0.25) +
	geom_smooth(data=df[df$Group2=="Brainspan",],aes(x=Days, y=PC1),se=F,size=0.25,color="black", linetype = 3,alpha=0.25) + 
	scale_y_continuous(name="PC1 of dev. regulated genes (n=5,567)",limits=c(-85,72))+
	scale_x_continuous(name="Post-conception days",trans="log2",breaks=bound) +
	scale_color_manual(values=c("darkgreen","red"))+
	scale_shape_manual(values=c(3,4))+
	theme_bw() +
	theme(
	axis.text.x=element_text(size=6,angle=45, hjust=1),
	axis.text.y=element_text(size=6),
	title= element_text(size=7),
	axis.ticks=element_line(size=0.05),
	axis.ticks.length=unit(0.01,"inch"),
	panel.grid.minor= element_blank(),
	legend.text= element_text(size=6),
	legend.position="",
	plot.background=element_rect(fill="transparent",colour=NA))

#t-hCS/hCS plot of PC1
df2= df[is.element(df$Group,c("hCS","t-hCS")),]
df2$Combine="combine"

p2=ggplot() + 
	geom_point(data= df2, aes(x=Combine, y=PC1,color=Group,fill=Group),shape=21,size=0.25,alpha=1,position=position_jitter(width=0.15,height=NULL)) +
	scale_y_continuous(name="PC1 of dev. regulated genes (n=5,567)",limits=c(-85,72))+
	scale_x_discrete(name="")+
	scale_color_manual(values=c("grey70","#129F49"))+
	scale_fill_manual(values=c("grey70","#129F49"))+
	theme_bw() +
	theme(
	axis.text.x=element_text(size=0,angle=45, hjust=1),
	axis.text.y=element_text(size=6),
	title= element_text(size=7),
	axis.ticks=element_line(size=0.05),
	axis.ticks.length=unit(0.01,"inch"),
	panel.grid.minor= element_blank(),
	panel.grid.major.x= element_blank(),
	legend.text= element_text(size=6),
	legend.position="right",
	plot.background=element_rect(fill="transparent",colour=NA))

setwd("../output")
pdf(file=paste0("PsychENCODE_hCS_t-hCS_PC1",".pdf"),width=5,height=2.5)
print(p1+p2+plot_layout(widths = c(2, 0.25)))
dev.off()




