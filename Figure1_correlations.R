#Code for Figure 1

setwd("<directory>")

#load packages
library(phyloseq)
library("ggplot2")
library(tidyverse)
library("cowplot")#to stack plots

#files
physeq3<-readRDS("physeq3.rds") #ITS data from woodrat_qiime2R.R
FQ<-read.csv(file="C:/Users/Sara/Dropbox/mycobiome/FungiQuant_Results_wMetadata_All_13Dec21.csv")
SG<-read.csv(file="kraken_bigDB_Conf0.05_outputs_30Mar22.csv") #kraken outputs from big DB

#filter
ps4<- transform_sample_counts(physeq3, function(x) replace(x, x<25,0) )  #set <25 to zero
ps5<-prune_taxa(taxa_sums(ps4) > 0, ps4) #remove anything no longer present in database

wild25<-subset_samples(ps5, Sample_type=="wild")
w25<-prune_taxa(taxa_sums(wild25) > 0, wild25)
dfReads<-data.frame(sample_sums (w25))
dfReads$sampleID<-sample_data(w25)$Woodrat.ID  

#fungiquant data
FQ_w<-filter(FQ, Sample_type == "wild")
FQ_w1<-FQ_w[,c(7,32)] #columns of interest

#get fungal proprotions from compiled bigDB kraken2 outputs (confidence = 0.05)
SG1<-SG[,c(4,14)]
SG1$woodrat.ID<-as.factor(SG1$woodrat.ID)

#the shotgun data is not a perfect subset of the ITS/fungiquant data, but 50 animals do overlap
df1<-merge(dfReads, SG1, by.x = "sampleID", by.y = "woodrat.ID", all.x = TRUE)

df2<-merge(df1, FQ_w1, by.x = "sampleID", by.y = "Woodrat_ID", all.x = TRUE)
df2$pCFun<-100*df2$pCFun  #convert proprotion to percent for figure

names(df2)<-c("Woodrat.ID", "ITS_reads", "SG_perc_fungi", "Copies_18S_per_ng")
df2$Woodrat.ID<-factor(df2$Woodrat.ID, levels = unique(df2$Woodrat.ID[order(df2$ITS_reads)]))

#outlier rat=953 (ITS_reads=  441253,  SG_perc_fungi=85.38% (!!) Copies_18S_per_ng=732142.7) 
df2[which(df2$ITS_reads > 200000),]

df3<-df2[which(df2$ITS_reads < 200000),]

#For color coding replace NAs with zeros,for sg, 18s THESE ARE NOT true zeros
df3[is.na(df3)] <- 0 

P1<-ggplot(df3, aes(y=ITS_reads, x=Woodrat.ID, color = ITS_reads > 0)) +
	geom_point()+
	scale_color_manual(name = 'PC1 > 0', values = setNames(c('black','grey50'),c(T, F))) +
	theme_classic()+
	theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.title.x=element_blank(), legend.position = "none")

P2<-ggplot(df3, aes(y=Copies_18S_per_ng, x=Woodrat.ID, color = Copies_18S_per_ng > 0)) +
	geom_point()+
	scale_color_manual(name = 'PC1 > 0', values = setNames(c('black','grey50'),c(T, F))) +
	theme_classic()+
	theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.title.x=element_blank(), legend.position = "none")

P3<-ggplot(df3, aes(y=SG_perc_fungi, x=Woodrat.ID, color=SG_perc_fungi > 0)) +
	geom_point()+
	scale_color_manual(name = 'PC1 > 0', values = setNames(c('black','grey50'),c(T, F))) +
	theme_classic()+
	theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.title.x=element_blank(), legend.position = "none")
pdf("sequencing_ITS_SG_18s_30Mar22.pdf", width=4, height=6)
plot_grid(P3,P2, P1, nrow=3, align = "v", labels = c('A', 'B', 'C'))
dev.off()

#using Kendall’s tau to measure correlations between 
#18s rRNA gene c/ng , number of ITS reads, and percent of shotgun reads assigned to fungi. 
#df2 has outlier, df3 does not

cor.test(df2$ITS_reads,df2$SG_perc_fungi, method="kendall") 