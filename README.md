<img align="right" src="https://github.com/SBWeinstein/Neotoma_fungi/blob/main/github_graphic-01.svg" width="20%">

# This respository contains code associated with "Wild herbivorous mammals (genus *Neotoma*) host a diverse but transient assemblage of fungi"

### Scripts for manuscript are as follows:

**Zac_scripts_here** sequence processing and taxonomy assignments? Amplicon sequences are available from the [NCBI sequence archive](https://www.ncbi.nlm.nih.gov/sra) under BioProject #####. Inputs are fastq.gz files. 

**woodrat_qiime2R.R:** Convert qiime files from Zac_scripts to phyloseq format, add metadata, check Saccharomyces_cerevisiae positive controls and blanks, remove positive controls. 

**Methods_comps_Fig1.R:** Compare the number of ITS reads, 18s copies per ng DNA, and percent of metagenomic reads assigned to fungi. Uses ITS data from phyloseq outputs from woodrat_qiime2R.R, Fungiquant results from FungiQuant_Results_wMetadata_All_13Dec21.csv, and kraken2 taxonomy summaries combined into phylo66_kraken_outputs.csv. Produces Figure 1.

**NOAA_precip_data.R:** Download precipitation data from NOAA 

**woodrat_ITS_diversity.R:** Calculates prevalence of fungal taxa and examines factors predicting fungal diversity, quantity and composition in wild rats. Script uses phyloseq formatted files produced by woodrat_qiime2R.R, precipitation data from NOAA_precip_data.R, compiled data on the ecology of identified ITS amplicons, and wild woodrat diet data (from Weinstein et al 2021). Generates figures 2A, 3A.

**Ecology_ASVs_Reads.R** Uses compiled data on the ecology of all identified fungal ASVs to examine the distribution of ecological niches and trophic modes across reads and ASVs.  Produces figure 2B, 2C
 
**Wild_cap_comps.R:** Script uses phyloseq formatted files produced by woodrat_qiime2R.R to examine ITS read counts and diversity in captive rats, with comparisons to wild individuals. Generates figures 4 and 5.

**MRM_models.R:** Script uses rarefied phyloseq formatted files produced in woodrat_ITS_diversity.R to examine how host phylogeny, geography and diet contribute to mycobiome composition in wild individuals. Generates Figure 3b and Table 3.
