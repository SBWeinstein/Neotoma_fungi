#Convert qiime files for analyses in R
#R version 4.1.0 
###################
#first time--install qiime2R
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20
###################

#load packages
library("qiime2R") #read in QIIME2 .qza files, convert to phyloseq object
library(phyloseq)
library("ggplot2")
library(tidyverse)

setwd("<directory>")
################################

#files that will be loaded in script below
#table.qza
#taxonomy.qza
#gnomex_metadata_ITS_trim.csv

#qza files to phyloseq object
physeq<-qza_to_phyloseq(
	 features="table.qza",
	 taxonomy="taxonomy.qza")

#add metadata 
#must be a dataframe & rownames must match the sample names in the otu_table
#original metadata contains 6 blanks that are not retained in OTU table, removed from trimmed metadata
metdat<-read.csv("gnomex_metadata_ITS_trim.csv", row.names=1) 
metdf<-data.frame(metdat) #convert to dataframe

physeq1 = merge_phyloseq(physeq, sample_data(metdf))#merge with phyloseq object

#quick look at taxa and taxa prevalence in everything sequenced
table(tax_table(physeq1)[, "Kingdom"], exclude = NULL)# all fungi, good

#count number of taxa in each phyla (4307 taxa distributed in 10 phyla)
# ~200 taxa in "NA" and "unidentified"
table(tax_table(physeq1)[, "Phylum"], exclude = NULL)

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(physeq1),
                 MARGIN = ifelse(taxa_are_rows(physeq1), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(physeq1),
                      tax_table(physeq1))

#look at relationship between prevalence and abundance for taxa in each phylum
ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(physeq1),color=Phylum)) +
 	geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) +  #line is 1%, present in 2-3 samples
	geom_point(size = 2, alpha = 0.7) +
  	scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  	facet_wrap(~Phylum) + theme_classic()+ theme(legend.position="none")
	
#check the positive controls
#nearly all saccharomyces cerevisiae; three ASVs; remove these Sc ASVs from other samples

#where do the positive control yeasts show up?
SC = subset_taxa(physeq1, Species=="Saccharomyces_cerevisiae")
SC_strains<-row.names(tax_table(SC))# the 6 Sc ASVs in the data set
pos_SC<-row.names(tax_table(pos0))# ASVs in the positive controls
pc<-is.element(SC_strains, pos_SC)#is the strain from the PC (T/F)
pc1<- paste0(pc, sep = "", 1:6)

#filter out 3 Sc strains in positive controls 
pcSC<- SC_strains[c(1,2,6)]
keep <- setdiff(taxa_names(physeq1), pcSC)
physeq2 <- prune_taxa(keep, physeq1)

#remove positive control samples
physeq3<-subset_samples(physeq2, Sample_type!="pos_control")

saveRDS(physeq3, file = "physeq3.rds")

#what's in the blanks? Not much...that's good
blanks<-subset_samples(physeq3 , Sample_type=="blank")
b2<-prune_taxa(taxa_sums(blanks) > 0, blanks) #remove any taxa that are no longer present
plot_bar(b2, fill="Genus")

#how many Neocallimastigomycota in unfiltered data?
p3calli = subset_taxa(physeq3, Phylum=="Neocallimastigomycota")
sum(sample_sums(p3calli))#29 reads from 3 ASVs in one sample (none >25, so lost in filtering)

#filter (25 reads), results in 1672 taxa
ps4<- transform_sample_counts(physeq3, function(x) replace(x, x<25,0) )  #set <25 to zero
ps5<-prune_taxa(taxa_sums(ps4) > 0, ps4) #remove anything no longer present in database

#fix extra space in "N. stephensi " in metadata
sampledf <- data.frame(sample_data(ps5))
sampledf$Rat_sp<-as.character(sampledf$Rat_sp)
sampledf["Rat_sp"][sampledf["Rat_sp"] == "N. stephensi "] <- "N. stephensi"
sample_data(ps5) <- sampledf
unique(sample_data(ps5)$Rat_sp)

wild25<-subset_samples(ps5, Sample_type=="wild")
w25<-prune_taxa(taxa_sums(wild25) > 0, wild25)

#the ITS data has field codes for animals from rio mesa, replace with animal IDs for merging w/ diet data set.
fix_codes<-read.csv("rio_mesa_animal_codes.csv")
y<- as.character(sample_data(w25)$Woodrat.ID)
y[15:21]<-fix_codes$Woodrat.ID
sample_data(w25)$Woodrat.ID<-y

saveRDS(w25, file = "phyloseq_fungal_w25_27Dec21.rds")