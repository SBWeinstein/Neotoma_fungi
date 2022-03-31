#distance matrix regression
#generates table 3, and Figure 3B

#setwd("<>") #set working directory

#load packages
library(ggplot2)
library(phyloseq)
library(phangorn)
library(tidyr)
library(geosphere)
library(ecodist)
library(eulerr)#for euler diagrams in figures
library(dplyr) #use to round values in MRM output table

#files
Neo_tree<-read.tree(file="Neotoma_phylosymtree.nwk") #Tree from Matoqc
site<-as.data.frame(read.csv("Site_lat_long.csv")) #lat/long data
plant<-readRDS( file = "diet1percent_rar_29May20.rds") 
Fun_eco<-read.csv("fungal_taxa_counts_25eco_wild_24May21.csv")
w25r<-readRDS("phyloseq_fungal_w25Rar_27Dec21.rds") #rarified wild fungal data from woodrat_ITS_diversity.R

#use phylogeny, diet, geography data to make distance matrices for MRM models
#make distance matrix for phylogeny using rat tree
#replace tip labels to match format of species names in mycobiome data
Neo_tree$tip.label<-c("N. albigula", "N. stephensi", "N. macrotis",
                      "N. cinerea", "N. bryanti", "N. lepida", "N. devia" )
#confirm that this matches names in dataset:
intersect(Neo_tree$tip.label,unique(sample_data(w25r)$Rat_sp))#7 there

#compute the pairwise distances between the pairs of tips using branch lengths
D_evo<-cophenetic(Neo_tree)
#convert into long format, species1, species2, distance
D_evo_long<-gather(as.data.frame(D_evo), key="Species2", value= "distance")
D_evo_long["Species1"]<-c(rep(rownames(D_evo), 7)) #make extra column for compared to species
D_evo_long <- D_evo_long[c(3,1,2)]  #reorder columns

#make distance matrix for geography using lat/long data
#confirm match names in dataset:
intersect(site$site,unique(sample_data(w25r)$site)) #19
D_site<-distm(cbind(site$longitude,site$latitude), fun=distGeo) #longitude first
row.names(D_site) <- site$site
colnames(D_site) <- site$site
D_site_km<-D_site/1000 #convert to km
#convert to long format
D_site_long<-gather(as.data.frame(D_site_km), key="Site2", value= "distance")
D_site_long["Site1"]<-c(rep(rownames(D_site_km), 19))
D_site_long <- D_site_long[c(3,1,2)]

#Make  diet distance matrix
#using individual diet data, at family level
psF<-tax_glom(plant, "Family")
#check name matchs
intersect(sample_data(psF)$Woodrat.ID,unique(sample_data(w25r)$Woodrat.ID)) #102 (106 in w25r, 155 in psF)

D_bray_diet<-phyloseq::distance(psF, "bray")
dietm<-as.matrix(D_bray_diet)
row.names(dietm) <- sample_data(psF)$Woodrat.ID
colnames(dietm) <- sample_data(psF)$Woodrat.ID
#convert to long format
D_diet_long<-gather(as.data.frame(dietm), key="pop2", value= "distance")
D_diet_long["pop1"]<-c(rep(rownames(dietm), dim(dietm)[1]))
D_diet_long <- D_diet_long[c(3,1,2)] #reorder columns


#######################
#create mycobiome distance matrix for rarified wild fungal: w25r 

#restrict to 102 animals with diet data
w25rP<-subset_samples(w25r, Woodrat.ID %in% row.names(dietm))
w25rP<-prune_taxa(taxa_sums(w25rP) > 0, w25rP)# 1221 taxa and 102 samples 
###################################

#analysis variation: rerun with separate plant, and not-plant associate fungal data sets
eco_plant<-filter(Fun_eco, ecology.niche == "plant associate") 
plant_assoc<-as.vector(eco_plant$X)  #624 taxa, convert from factor to vector for prune function
#subset taxa in w25rP to rownames of plant associates
ps_pa<-prune_taxa(plant_assoc, w25rP) #505 taxa

eco_Notplant<-filter(Fun_eco, ecology.niche != "plant associate") 
Notplant_assoc<-as.vector(eco_Notplant$X)  #885 taxa, convert from factor to vector for prune function
#subset taxa in w25rP to rownames of plant associates
ps_NP<-prune_taxa(Notplant_assoc, w25rP) #716 taxa and 102 samples 

#PAs are not in all samples, NPs are in all samples
ps_pa1<-prune_samples(sample_sums(ps_pa) > 0, ps_pa) # 505 taxa and 100 samples
ps_NP1<-prune_samples(sample_sums(ps_NP) > 0, ps_NP) # 716 taxa and 102 samples

########
#MRM models
ps_data<-sample_data(ps_NP1) #full data: w25rP, plant: ps_pa1, Notplant: ps_NP1

###Make big host phylogeny matrix to match microbiome dimensions
#start with blank new matrix with dims and names matching the microbiome matrix
my.matrix <- matrix(0, nrow = length(ps_data$Rat_sp), ncol = length(ps_data$Rat_sp))
row.names(my.matrix) <- ps_data$Rat_sp
colnames(my.matrix) <- ps_data$Rat_sp

#use the table to fill in values in matrix with dimensions matching all the samples
for (i in 1:nrow(D_evo_long)) {
   my.matrix[rownames(my.matrix) == D_evo_long$Species2[i],
	 colnames(my.matrix) == D_evo_long$Species1[i]] <- D_evo_long$distance[i]
}
big_evo_matrix<-my.matrix #rename

#Make big SITE matrix to match microbiome dimensions
#make blank new matrix with dims and names matching the microbiome matrix
big.site.matrix <- matrix(0, nrow = length(ps_data$site), ncol = length(ps_data$site))
row.names(big.site.matrix) <- ps_data$site
colnames(big.site.matrix) <- ps_data$site

#use the table to fill in values in matrix with dimensions matching all the samples
for (i in 1:nrow(D_site_long)) {
   big.site.matrix[rownames(big.site.matrix) == D_site_long$Site2[i],
	 colnames(big.site.matrix) == D_site_long$Site1[i]] <- D_site_long$distance[i]
}

#make big DIET matrix with dims and names matching the microbiome matrix
big.diet.matrix <- matrix(0, nrow = length(ps_data$Woodrat.ID), ncol = length(ps_data$Woodrat.ID))
row.names(big.diet.matrix) <- ps_data$Woodrat.ID
colnames(big.diet.matrix) <- ps_data$Woodrat.ID

#use the table to fill in values in matrix with dimensions matching all the samples
for (i in 1:nrow(D_diet_long)) {
   big.diet.matrix[rownames(big.diet.matrix) == D_diet_long$pop2[i],
	 colnames(big.diet.matrix) == D_diet_long$pop1[i]] <- D_diet_long$distance[i]
}

###################################################
###calculate values for jaccard distance
output<-matrix(0, nrow =1,  ncol =19 )
colnames(output)<-c("dist", "DPSr", "DPr", "DSr", "PSr", "Sr", "Dr", "Pr",
					"Diet","Phy", "Site",
				"PhyDiSi", "DietPhy", "PhySite", "DietSite","DPSp","Sp", "Dp", "Pp")


	#alter myco data input: full: w25rP, plant assoc: ps_pa1, not_plant: ps_NP1
	D_bray_ps<-phyloseq::distance(ps_NP1, method="jaccard")#"jaccard" "bray"
	dm<-as.matrix(D_bray_ps)

	#matrices need to be a "distance matrix"
	DPS<-MRM(as.dist(dm) ~ as.dist(log10(big.site.matrix+1)) +
          as.dist(big_evo_matrix)+ as.dist(big.diet.matrix),  nperm=1000)

	#For variance partioning
	DP<-MRM(as.dist(dm) ~  as.dist(big_evo_matrix)+ as.dist(big.diet.matrix),  nperm=1000)
	DS<-MRM(as.dist(dm) ~  as.dist(log10(big.site.matrix+1))+ as.dist(big.diet.matrix),  nperm=1000)
	PS<-MRM(as.dist(dm) ~  as.dist(log10(big.site.matrix+1))+ as.dist(big_evo_matrix),  nperm=1000)
	S<-MRM(as.dist(dm) ~ as.dist(log10(big.site.matrix+1)),  nperm=1000)
	D<-MRM(as.dist(dm) ~ as.dist(big.diet.matrix),  nperm=1000)
	P<-MRM(as.dist(dm) ~ as.dist(big_evo_matrix), nperm=1000)

	DPSr<-DPS$r.squared[1] #r values are 1
	DPSp<-DPS$r.squared[2] #p values are 2
	DPr<-DP$r.squared[1]
	DSr<-DS$r.squared[1]
	PSr<-PS$r.squared[1]
	Sr<-S$r.squared[1]
	Sp<-S$r.squared[2]
	Dr<-D$r.squared[1]
	Dp<-D$r.squared[2]
	Pr<-P$r.squared[1]
	Pp<-P$r.squared[2]
	
	j<-1
	output[j,1]<- "Distance"
	output[j,2]<-DPS$r.squared[1]
	output[j,3]<-DP$r.squared[1]
	output[j,4]<-DS$r.squared[1]
	output[j,5]<-PS$r.squared[1]
	output[j,6]<-S$r.squared[1]
	output[j,7]<-D$r.squared[1]
	output[j,8]<-P$r.squared[1]

	#get contributions
	output[j,9]<-DPSr-PSr
	output[j,10]<-DPSr-DSr
	output[j,11]<-DPSr-DPr
	#overlaps
	output[j,12]<- DPSr +Dr + Pr+ Sr-DPr-DSr-PSr
	output[j,13]<- PSr +  DSr - Sr - DPSr
	output[j,14]<- DSr +DPr - Dr -DPSr
	output[j,15]<- PSr + DPr - Pr-DPSr
	#pvalues
	output[j,16]<- DPSp
	output[j,17]<- Sp
	output[j,18]<- Dp
	output[j,19]<- Pp
	

op1<-data.frame(output)

op2<-dplyr::mutate_if(op1[,-1], is.character, ~ as.numeric(as.character(.x)))
op3<-round(op2,digits=3)
op4<-data.frame(metric=op1[,1], op3)
Jac<-op4

#for making Euler diagram
D100<-as.vector(Jac$Diet)*100
P100<-as.vector(Jac$Phy)*100
S100<-as.vector(Jac$Site)*100
PS100<-as.vector(Jac$PhySite)*100
PD100<-as.vector(Jac$DietPhy)*100
DS100<-as.vector(Jac$DietSite)*100
PDS100<-as.vector(Jac$PhyDiSi)*100

##plot as venn diagram
#the R circle gives a reference size for 100% of variance
combo_wild <- c(Diet = D100, Phylogeny = P100, Site = S100,
	"Diet&Phylogeny" = PD100, "Diet&Site" = DS100, "Phylogeny&Site" = PS100,
	"Diet&Phylogeny&Site"=PDS100, "R"=100, "R&Phylogeny"=.0001)


ven<-plot(euler(combo_wild, input = c("disjoint"),shape =  "ellipse"),
	fills = list(fill = c("darkgreen","darkmagenta", "blue3" ), alpha = 0.8),
	quantities = FALSE)

pdf("ITS_jac_venn_NotPlantAssoc__21Jan22.pdf", width=4, height=4)
ven
dev.off()

