#fungal ecology from ITS amplicon data
#Create figures 2B,2C

#load packages
library(phyloseq)
library("ggplot2")
library(tidyverse)
library(eulerr)

setwd("<directory>")

#files
Fun_eco<-read.csv("fungal_taxa_counts_25eco_wild_24May21.csv") #ASVs with assigned ecology
physeq3<-readRDS("physeq3.rds") #ITS data from woodrat_qiime2R.R

#create figures
#Figure 2B:ecology of ASV diversity, unrarified dataset with 1509 ASVs
eco_tab<-table(Fun_eco[, "ecology.niche"], exclude = NULL)
df1<-data.frame(eco_tab/1509)
df2<- df1[order(df1$Freq),]

df2$Var1 <- factor(df2$Var1, levels=df2$Var1[order(df2$Freq)], ordered=TRUE)

hsize <- 1.5

df3 <- df2 %>% 
  mutate(x = hsize)

df3$perc<- round(df3$Freq*100)

greys<-c("grey10", "grey95",
		"grey30", "grey70",
		"white", "grey50",
		"black", "grey20",
		"grey80")
pdf("Eco_pie_ASVs.pdf", width=6, height=4)
ggplot(df3, aes(x = hsize, y = Freq, fill = Var1)) +
 	 geom_col(color="black") +
  	coord_polar(theta = "y") +
  	xlim(c(0.2, hsize + 0.5))+
  	scale_fill_manual(values=greys)+
	geom_text(aes(label = perc),
             position = position_stack(vjust = 0.5)) +
	#theme_classic()+
	theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
dev.off()


#for ecology based on ASV read counts, use rarified data to reduce individual differences in seq depth
ps4<- transform_sample_counts(physeq3, function(x) replace(x, x<25,0) )  #set <25 to zero
ps5<-prune_taxa(taxa_sums(ps4) > 0, ps4) #remove anything no longer present in database
wild25<-subset_samples(ps5, Sample_type=="wild")
w25<-prune_taxa(taxa_sums(wild25) > 0, wild25)
w25R1<-rarefy_even_depth(w25, sample.size = 1000, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)

#Figure 2C
R2<-data.frame(taxa_sums(w25R1))
R3<-merge(R2, Fun_eco, by.x=0, by.y=1)
R4<-data.frame(niche=R3$ecology.niche, count=R3$taxa_sums.w25R1.)
R5<-data.frame(summarize_all(group_by(R4, niche), sum))
R5$prop<-R5$count/sum(R5$count)
R5$perc<- round(R5$prop*100)

R5$niche <- factor(R5$niche, levels=df2$Var1) #use same order as prev

R7 <- R5 %>% 
  mutate(x = hsize)

pdf("Eco_pie_Reads.pdf", width=6, height=4)
ggplot(R7, aes(x = hsize, y = prop, fill = niche)) +
 	 geom_col(color="black") +
  	coord_polar(theta = "y") +
  	xlim(c(0.2, hsize + 0.5))+
  	scale_fill_manual(values=greys)+
	geom_text(aes(label = perc),
             position = position_stack(vjust = 0.5)) +
	#theme_classic()+
	theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
dev.off()

#make euler diagrams for three trophic modes: inset for 3B, 3C)
troph_tab<-table(Fun_eco[, "Trophic_mode"], exclude = NULL)
troph_tab<-data.frame(troph_tab)
#for 74% of ASVs (1119/1509)  with identifiable trophic mode  
troph_ASV<- c(Path = 139, Sap = 654, Sym = 25, 
		"Path&Sap" =136, "Sym&Path" = 2, "Sym&Sap" =14,
		"Path&Sap&Sym" = 149)

ven<-plot(euler(troph_ASV, input = c("disjoint"),shape =  "ellipse"),
	fills = list(fill = c("grey40","grey5", "grey80" ), alpha = 0.95),
	quantities = FALSE)

pdf("Troph_ASV_venn__17jan22.pdf", width=3, height=3)
ven
dev.off()

#for reads, use merged dataframe from read count eco pie above
R4reads<-data.frame(niche=R3$Trophic_mode, count=R3$taxa_sums.w25R1.)
R5reads<-data.frame(summarize_all(group_by(R4reads, niche), sum))
#vens for 82.7% of reads with assigned trophic mode

troph_reads<- c(Path = 12256, Sap = 59524, Sym = 763, 
		"Path&Sap" =12456, "Sym&Path" = 141, "Sym&Sap" =104,
		"Path&Sap&Sym" = 2382)

ven<-plot(euler(troph_reads, input = c("disjoint"),shape =  "ellipse"),
	fills = list(fill = c("grey40","grey5", "grey80" ), alpha = 0.95),
	quantities = FALSE)

pdf("Troph_reads_venn__17jan22.pdf", width=3, height=3)
ven
dev.off()