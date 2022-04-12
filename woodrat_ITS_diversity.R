#Fungal prevalence, diversity, quantity, and composition analyses
#R version 3.6.3 or v4.1.0 

###################
#packages
library(phyloseq)
library("ggplot2")
library(tidyverse)
library("cowplot")#to stack plots
library(MASS) #for regressions
library(vegan) # for community comparisons
library(emmeans) #pairwise comparisons
library(multcomp) #pairwise comparisons, CLD

setwd("<directory>")

#files
w25<-readRDS("phyloseq_fungal_w25_27Dec21.rds") # filtered, wild rat ITS data  from woodrat_qiime2R.R script
Fun_eco<-read.csv("fungal_taxa_counts_25eco_wild_24May21.csv") #amplicon ecology data to get plant-associates
precip<-read.csv("precip.by.site.csv") #precip data from NOAA, with site info added in excel
plant<-readRDS( file = "diet1percent_rar_29May20.rds") 
FQ<-read.csv(file="FungiQuant_Results_wMetadata_All_13Dec21.csv") #fungiquant data

#############
w25_sp<-tax_glom(w25, "Species", NArm=FALSE) # 551 "species"

#calculate prevalence of each species
wb_spP<- transform_sample_counts(w25_sp, function(x) replace(x, x>0,1) ) 
dfP<-data.frame(taxa_sums(wb_spP))
ttP<-data.frame(as(tax_table(wb_spP), "matrix"))

dft<-merge( ttP,dfP, by=0)
colnames(dft)[9]<-"Count"  #rename column
dft$prop<-dft$Count/120
dft_sort <- dft[order(dft$Count),]  #table of "species" prevalence in wild rats

#figure 2A histogram
pdf("species_histo_5jan22.pdf",  width= 7, height=3.5)
ggplot(dft, aes(Count/120)) + 
  geom_histogram(color="black", fill="grey50",binwidth= 0.02) + 
	ylab("Count of fungal species") + 
	xlab("Proportion of animals")+
	theme_classic()
dev.off()

###Predicting fungal diversity in wild rats
#build data.frame, starting with sample_data(w25)
div_dat<-data.frame(sample_data(w25))

set.seed(42) #for rarifying

obs.rich.ITS<- estimate_richness(w25, measures="Observed") #not rarefied
div_dat<-cbind(div_dat, obs.rich.ITS)

#add diversity metrics for rarefied data, use 1000 read cut off- 14 samples removed
w25R<-rarefy_even_depth(w25, sample.size = 1000, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
saveRDS(w25R, file = "phyloseq_fungal_w25Rar_5Jan22.rds")

obs.rich.R <-estimate_richness(w25R, measures="Observed")
obs.shan.R <-estimate_richness(w25R, measures="Shannon")

dd4<-data.frame(sample_data(w25R)$Woodrat.ID,obs.rich.R,obs.shan.R )
colnames(dd4)<-c("Woodrat.ID", "obs.ITS.R", "shan.ITS.R")
dd2<-merge(div_dat, dd4, by = "Woodrat.ID", all.x=TRUE)

div_dat2<-merge(dd2, precip, by="site")  #note that sites FC,WR, LR have extra space after name

####add plant-associated (PA) fungal diversity
eco_plant<-filter(Fun_eco, ecology.niche == "plant associate") 
plant_assoc<-as.vector(eco_plant$X)  #624 taxa, convert from factor to vector for prune function
#subset taxa in w25 to rownames of plant associates
ps_pa<-prune_taxa(plant_assoc, w25)

#calculate richness of PAs
obs.rich.ITS.PA<- estimate_richness(ps_pa, measures="Observed") #not rarefied
dpa<-data.frame(sample_data(w25)$Woodrat.ID,obs.rich.ITS.PA)
names(dpa)<-c("Woodrat.ID", "Obs_plant_assoc")
#merge with div_dat2 dataframe
div_dat3<-merge(div_dat2, dpa, by="Woodrat.ID")

#add diet diversity
psF<-tax_glom(plant, "Family")
df_plant<-data.frame(estimate_richness(psF, measures=c("Observed", "Shannon")), #observed, shannon richness of families per sample
				Woodrat.ID=sample_data(psF)$Woodrat.ID)			
colnames(df_plant)<-c("Obs.PF", "Sh.PF","Woodrat.ID")

div_dat4<-merge(div_dat3, df_plant, by= "Woodrat.ID", all.x=TRUE) #a few missing diets, but mostly there.

div_dat4$site<-as.character(div_dat4$site)
div_dat4["site"][div_dat4["site"] == "Rio Mesa, UT"] <- "Rio Mesa"  #rename to be consistent in plots
#treat Pioneertown and Pioneertown Reserve as one site, as very close to each other
div_dat4["site"][div_dat4["site"] == "Pioneertown Reserve"] <- "Pioneertown"

div_dat4$site = reorder(div_dat4$site, div_dat4$Observed, mean)

#for figure 3A
pdf("ITS_richness_by_site_5Jan22.pdf", width=3.46457, height=4, useDingbats=FALSE)
ggplot(div_dat4, aes(site, Observed)) + 
  geom_boxplot(outlier.shape = NA, fill = "grey90", color="black" )+
	geom_jitter(aes(shape = Rat_sp, , fill = Rat_sp), width=.05, height=0)+
	scale_shape_manual(values = c(21, 22, 24, 25, 21,22,23)) +
  	scale_fill_manual(
    		values = c("grey95", "grey40", "grey40", "black","black", "grey95", "black")
  				)+
	theme_classic()+
	theme(axis.text.x=element_text(angle=90,hjust=1, vjust=0.1), 
		legend.position="bottom", legend.title = element_blank())+
	labs(y = "ASV Richness (ITS)", x= "Site")
dev.off()

#######################################
#what predicts total fungal diversity? and plant-associated fungal diversity?
#Note that residual deviance >> df, using negative binomial error distribution
#5 animals are missing plant data- remove for model comparisons- drops to 115 samples
#as diet not sig predictor, return to full dataset after factor dropped
complete_plants<-div_dat4[complete.cases(div_dat4$Obs.PF), ] 

##for plant associated fungal (observed diversity)- only site and species significant factors
m1 <- glm.nb(Obs_plant_assoc~site+Rat_sp+precip+Obs.PF, data = complete_plants)
m1.p <- glm.nb(Obs_plant_assoc~Rat_sp+precip+Obs.PF, data = complete_plants)#drop precip
m2 <- glm.nb(Obs_plant_assoc~site+Rat_sp+Obs.PF, data = complete_plants)
m3 <- glm.nb(Obs_plant_assoc~site+Rat_sp, data = complete_plants)
m4 <- glm.nb(Obs_plant_assoc~site+Rat_sp, data = div_dat4) #site and species also sig with full dataset
anova(m1,m2)
anova(m2,m3)
anova(m4, test='LRT')
summary(m4)

##model checks, use for all linear models
model= m4
plot(fitted(model), residuals(model, type="pearson"), main= "Residuals vs Fitted", xlab = "Fitted Values", ylab = "Residuals")  
abline(h = 0, lty = 2)  
qqnorm(residuals(model))

###
##for all fungal diversity: rarefied obs.ITS.R- sig site and species
m1 <- glm.nb(obs.ITS.R~site+Rat_sp+precip+Obs.PF, data = complete_plants)
m1.p <- glm.nb(obs.ITS.R~Rat_sp+precip+Obs.PF, data = complete_plants) #colinearities between site and rainfall
m2 <- glm.nb(obs.ITS.R~site+Rat_sp+Obs.PF, data = complete_plants)
m3 <- glm.nb(obs.ITS.R~site+Rat_sp, data = complete_plants)
m4b <- glm.nb(obs.ITS.R~site+Rat_sp, data = div_dat4) #best

anova(m2,m3)#compare between models
anova(m4b, test='LRT') #model comparison for stats table
summary(m4b)

##for all fungal diversity: not rarefied Observed- sig site 
m1 <- glm.nb(Observed~site+Rat_sp+precip+Obs.PF, data = complete_plants)
m1.p <- glm.nb(Observed~Rat_sp+precip+Obs.PF, data = complete_plants) #colinearities between site and rainfall
m2 <- glm.nb(Observed~site+Rat_sp+Obs.PF, data = complete_plants)
m3 <- glm.nb(Observed~site+Rat_sp, data = complete_plants) #species is close (o.o6), but not sig
m3a <- glm.nb(Observed~site+Rat_sp, data = div_dat4)
m4a <- glm.nb(Observed~site, data = div_dat4) #best

anova(m4a) 

#Which sites significantly differ? for Compact Letter Display, Figure 1A
marginal = emmeans(m4a, ~ site)
cld(marginal, Letters=letters)				 

#shan.ITS.R (rarified) not count data, ranges 0-4, surprisingly normal distribution, good fit for lm, only site sig
hist(div_dat4$shan.ITS.R)
m1 <- lm(shan.ITS.R~site+Rat_sp+precip+Obs.PF, data = complete_plants)
m2 <- lm(shan.ITS.R~site+Rat_sp+Obs.PF, data = complete_plants)
m3 <- lm(shan.ITS.R~site+Rat_sp, data = complete_plants)
m4 <- lm(shan.ITS.R~site+Rat_sp, data = div_dat4)
m5<-lm(shan.ITS.R~site, data = div_dat4)

summary(m4)
anova(m5)

#does 18s quantity (funqiquant data) follow the same patterns? 
FQ_w<-filter(FQ, Sample_type == "wild")
#FQ_w1 <- filter(FQ_w, Woodrat_ID != "953") #953 is the outlier, examine with and w/o

FQW<-FQ_w[,c(7,32)]
names(FQW)[1]<-c("Woodrat.ID")
div_dat5<-merge(div_dat4, FQW, by="Woodrat.ID")#add to diversity data frame
div_dat5Out<- filter(div_dat5, Woodrat.ID != "953") #953 is the high outlier, , examine with and w/o

hist(div_dat5Out$Copies_18S_per_ng)#not normaly distrubuted, log10 transform- residuals, etc look pretty good
hist(log10(div_dat5Out$Copies_18S_per_ng))
complete_plants5<-div_dat5Out[complete.cases(div_dat5Out$Obs.PF), ] 
m1 <- lm(log10(Copies_18S_per_ng)~site+Rat_sp+precip+Obs.PF, data = complete_plants5)
m1.2 <- lm(log10(Copies_18S_per_ng)~Rat_sp+precip+Obs.PF, data = complete_plants5)
m2 <- lm(log10(Copies_18S_per_ng)~Rat_sp+precip, data = complete_plants5)
m3 <- lm(log10(Copies_18S_per_ng)~site+Rat_sp, data = complete_plants5)
m4 <- lm(log10(Copies_18S_per_ng)~site+Rat_sp, data = div_dat5Out)
m5 <- lm(log10(Copies_18S_per_ng)~site, data = div_dat5Out)

anova(m5,m4)
anova(m5)

#############################################################
##predictors of composition
#w25R:the 1000 read subsampled dataset w/ 106 samples
#Permanova, adonis function, testing for different centroids
#jaccard distance
psJac <- phyloseq::distance(w25R, method = "jaccard")
# make  dataframe from the sample_data
sampledf <- data.frame(sample_data(w25R))
sampledf["site"][sampledf["site"] == "Pioneertown Reserve"] <- "Pioneertown"

#Adonis test
adonis(psJac ~ site + Rat_sp, data = sampledf)
adonis(psJac ~ Rat_sp+site, data = sampledf)
#sites and species have different centroids.  

# Homogeneity of dispersion test
beta <- betadisper(psJac, sampledf$site)
beta <- betadisper(psJac, sampledf$Rat_sp)
permutest(beta)
# betadisper significant- groups have different dispersions

 
