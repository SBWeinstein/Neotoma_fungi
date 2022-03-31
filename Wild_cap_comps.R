#Neotoma fungal captive/wild comparisons

#packages
library(phyloseq)
library("ggplot2")
library(tidyverse)
library("cowplot")#to stack plots
library(ggbeeswarm) #prettier version of jittering points
library(eulerr)# for euler diagrams
library(gplots)

setwd("<directory>")

#files
physeq3<-readRDS("physeq3.rds") #ITS amplicon data
FQ<-read.csv(file="FungiQuant_Results_wMetadata_All_13Dec21.csv") #Fungiquant data

########################################
#filter, results in 1672 taxa
ps4<- transform_sample_counts(physeq3, function(x) replace(x, x<25,0) )  #set <25 to zero
ps5<-prune_taxa(taxa_sums(ps4) > 0, ps4) #remove anything no longer present in database

#fix extra space in "N. stephensi " 
sampledf <- data.frame(sample_data(ps5))
sampledf$Rat_sp<-as.character(sampledf$Rat_sp)
sampledf["Rat_sp"][sampledf["Rat_sp"] == "N. stephensi "] <- "N. stephensi"
sample_data(ps5) <- sampledf
unique(sample_data(ps5)$Rat_sp)

cap25<-subset_samples(ps5, Sample_type=="captive")
c25<-prune_taxa(taxa_sums(cap25) > 0, cap25)

wild_25<-subset_samples(ps5, Sample_type=="wild")
w_25<-prune_taxa(taxa_sums(wild_25) > 0, wild_25)

ssc<-data.frame(reads=sample_sums(c25))
#how many of the 107 rats have more than 100 reads? more than 1000?
sum(ssc$reads > 0) #35
sum(ssc$reads > 100) #15
sum(ssc$reads > 1000) #1

c25f <- prune_samples(sample_sums(c25)>0, c25)
ssc0<-data.frame(reads=sample_sums(c25f))
ssc0_Out<- filter(ssc0, reads != 10378) 
mean(ssc0_Out$reads)
sd(ssc0_Out$reads)

#diversity in captive samples with reads: ASVs, genera, families
c_g<-tax_glom(c25f, "Family")
c_g1<-estimate_richness(c_g, measures = c("Observed"))
mean(c_g1$Observed)
sd(c_g1$Observed)

##compare changes in reads, and observed richness (ASV, genera, families) 
#between wild and captive samples, using c25 (107 rats) and w_25 (120 rats)

divC<-data.frame(sample_data(c25))
divC1 <- dplyr::select(divC, Woodrat.ID, Rat_sp,ng.ul)
ASV_C<- estimate_richness(c25, measures="Observed") 
C_gen<-estimate_richness(tax_glom(c25, "Genus"),measures="Observed") 
C_fam<-estimate_richness(tax_glom(c25, "Family"),measures="Observed")
ssc<-data.frame(reads=sample_sums(c25))
divC2<-cbind(ObsC_AVS = ASV_C$Observed, ObsC_gen = C_gen$Observed, 
			ObsC_fam = C_fam$Observed, readsC = ssc$reads)
divC3<-cbind(divC1, divC2) 

#repeat for wilds
divW<-data.frame(sample_data(w_25))
divW1 <- dplyr::select(divW, Woodrat.ID, ng.ul)
ASV_W<- estimate_richness(w_25, measures="Observed") 
W_gen<-estimate_richness(tax_glom(w_25, "Genus"),measures="Observed") 
W_fam<-estimate_richness(tax_glom(w_25, "Family"),measures="Observed")
ssw<-data.frame(reads=sample_sums(w_25))
divW2<-cbind(ObsW_AVS = ASV_W$Observed, ObsW_gen = W_gen$Observed, 
			ObsW_fam = W_fam$Observed, readsW = ssw$reads)
divW3<-cbind(divW1, divW2) 

divCW<-merge(divC3, divW3, by = "Woodrat.ID") #107 rats

#difference in amount of gDNA in wild and cap samples?
#compare: ng.ul.x ng.ul.y
#F test to compare variances, they are unequal
var.test(divCW$ng.ul.x, divCW$ng.ul.y)
#default R t-test is Welch t-test
t.test(divCW$ng.ul.x, divCW$ng.ul.y)

#diff in ASVs between wild and cap. use non-parametric stats
wilcox.test(divCW$ObsC_AVS, divCW$ObsW_AVS, paired = TRUE, alternative = "two.sided")
summary(lm(divCW$ObsW_AVS~divCW$ObsC_AVS))

mean(divCW$ObsC_AVS)
sd(divCW$ObsC_AVS)
mean(divCW$ObsW_AVS)
sd(divCW$ObsW_AVS)

#diff in read counts between wild and cap?
#non-parametric alternative to paired t-test
#note that R returns "V" which is related to the Z statistic 
t.test(divCW$ObsC_AVS, divCW$ObsW_AVS)
wilcox.test(divCW$readsC, divCW$readsW, paired = TRUE, alternative = "two.sided")
mean(divCW$readsC)
sd(divCW$readsC)
mean(divCW$readsW)
sd(divCW$readsW)

########################
#Figure 4: 3 stacked plots with  w/c reads, ASVs, and fungiquant data
CW_reads<- divCW[, c(1,7,12)]
#CW_reads1 <- filter(CW_reads, Woodrat.ID != "953")
gathercols <- c("readsC", "readsW")
CW_long<-gather(CW_reads, type, reads, gathercols )
CW_long$type<-factor(CW_long$type, levels=c("readsW", "readsC"))
P1<-ggplot(CW_long, aes(x=type, y=(reads), color= type, fill =type))+
	geom_boxplot(outlier.shape = NA)+
	geom_quasirandom(method = "pseudorandom", alpha=0.3, size=1.5)+
	 scale_fill_manual(values=c("grey10", "grey80"))+
	scale_color_manual(values=c("grey40", "grey40"))+
	geom_hline(yintercept=62859 , linetype="dashed", color= "grey60")+ #mean chow reads
	theme_classic()+
	labs(y="ITS Reads")+
	ylim(0, 160000)+  #add outlier with y axis break back later 
	theme(legend.position="none", axis.title.x = element_blank(),
			axis.text.x=element_blank())

CW_ASVs<- divCW[, c(1,4,9)]
gathercols <- c("ObsC_AVS", "ObsW_AVS")
ASV_long<-gather(CW_ASVs, type, ASVs, gathercols )
ASV_long$type<-factor(ASV_long$type, levels=c("ObsW_AVS", "ObsC_AVS"))
P2<-ggplot(ASV_long, aes(x=type, y=(ASVs), color= type, fill =type))+
	geom_boxplot(outlier.shape = NA)+
	geom_quasirandom(method = "pseudorandom", alpha=0.3, size=1.5)+
	scale_fill_manual(values=c("grey10", "grey80"))+
	scale_color_manual(values=c("grey40", "grey40"))+
	geom_hline(yintercept=132, linetype="dashed", color= "grey60")+
	theme_classic()+
	labs(y="ITS ASVs")+
	theme(legend.position="none", axis.title.x = element_blank(), 
			axis.text.x=element_blank())

#Fungiquant data
FQ_wc<-filter(FQ, Sample_type == "wild" |
			 Sample_type == "captive")
FQ_wc1 <- filter(FQ_wc, Copies_18S_per_ng < 732142) #excludes outlier 953

P3<-ggplot(FQ_wc1, aes(x=Sample_type, y=Copies_18S_per_ng, fill = Sample_type, color = Sample_type)) +
	geom_boxplot(outlier.shape = NA)+
	geom_quasirandom(method = "pseudorandom", alpha=0.3, size=1.5)+
	scale_fill_manual(values=c("grey10", "grey80"))+
	scale_color_manual(values=c("grey40", "grey40"))+
	geom_point(data=FQ_wc953, aes(x=Sample_type, y=Copies_18S_per_ng))+
	geom_hline(yintercept=247759, linetype="dashed", color= "grey60")+
	theme_classic()+
	scale_y_log10()+
	labs(y="18s c/ng")+
	theme(legend.position="none", axis.title.x = element_blank(), 
		axis.text.x=element_blank())

#one column width in journal= 88mm, convert to inches
pdf("WC_comps_BW_22Jan22.pdf", width=3.46457, height=6)
plot_grid(
 P3,P1, P2,
   ncol = 1, align = "v", labels = "AUTO")
dev.off()

########################
###overlap in ASVs between wild, captive, cage, chow
wTax<-rownames(tax_table(w_25))#1509
cTax<-rownames(tax_table(c25)) #32

controls<-subset_samples(ps5, Sample_type=="environ_control")
cages<-c("17950X232", "17950X233", "17950X234")
chows<-c("17950X235", "17950X236", "17950X237")
cage<-subset_samples(controls, sample_names(controls) %in% cages)
cage1<-prune_taxa(taxa_sums(cage) > 0, cage)#3 ASVs
chow<-subset_samples(controls, sample_names(controls) %in% chows)
chow1<-prune_taxa(taxa_sums(chow) > 0, chow) #183 asvs

chowTax<-rownames(tax_table(chow1))#183
cageTax<-rownames(tax_table(cage1)) #3

#amounts are unique to each catergory. eg "A" is only in A, no where else.
taxList <- list(wTax=wTax, cTax=cTax, chowTax=chowTax, cageTax=cageTax)
ItemsList <- venn(taxList, show.plot = FALSE)

lengths(attributes(ItemsList)$intersections)
overlaps<-c("Captive&Chow"=6, "Wild&Cage"=2, "Wild&Chow"=23, "Wild&Captive&Chow"=7, 
             Captive=8, Wild=1466, Cage=1, Chow=147, "Wild&Captive"=11)
fit3 <- euler(overlaps, shape = "ellipse")

#Figure 5A
pdf("ITS_taxa_venn_BW_30Dec21.pdf", width=3.46457, height=3.46457)
plot(fit3, quantities = TRUE,
		fills=list(fill=c("grey90" ,"grey60","grey10", "white"), alpha=0.9),
		edges=c("black" ,"black","black", "black"))
dev.off()

wild_retained<-subset_taxa(c25, rownames(otu_table(c25)) %in% rownames(tax_table(w_25)))
plot_bar(wild_retained, fill="Genus")

retained_gen<-tax_glom(wild_retained, "Genus")
#how prevalent are these?
# Compute prevalence of each feature, store as dataframe
prevdf = apply(X = otu_table(wild_retained),
                 MARGIN = ifelse(taxa_are_rows(wild_retained), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
percent<-(prevdf/107)*100
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf, Percent = percent,
                      TotalAbundance = taxa_sums(wild_retained),
                      tax_table(wild_retained))
prevdf<-prevdf[order(-prevdf$Percent),]
#look at table
view<-print(prevdf[,-c(4,6,7)], row.names = F)

#there are 8 taxa only seen in captive rats, what are they?
#use overlaps from venn diagram prep
capOnly<-attributes(ItemsList)$intersections$cTax

capOnly1<-subset_taxa(c25, rownames(otu_table(c25)) %in% capOnly )
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(capOnly1),
                 MARGIN = ifelse(taxa_are_rows(capOnly1), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
percent<-(prevdf/107)*100
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf, Percent = percent,
                      TotalAbundance = taxa_sums(capOnly1),
                      tax_table(capOnly1))
prevdf<-prevdf[order(-prevdf$Percent),]
#look at table
view<-print(prevdf[,-c(4,6,7)], row.names = F)

#captive chow overlap:  "cTax:chowTax"
capchow<-subset_taxa(c25, rownames(otu_table(c25)) %in% chowTax )
capsums<-sample_sums(c25)
capchowsums<-sample_sums(capchow)
propchow <-capchowsums/capsums
sums<-data.frame(capsums, capchowsums,propchow )
mean(sums$propchow, na.rm = TRUE)
sd(sums$propchow, na.rm = TRUE)

chowPhy<-tax_glom(chow, "Phylum")
sample_sums(chow)
chp<-data.frame(as(tax_table(chowPhy), "matrix"))
chp$sums<-taxa_sums(chowPhy) 
chp[,3:7]<-NULL
chp$prop<- chp$sum/ sum(chp$sums) 