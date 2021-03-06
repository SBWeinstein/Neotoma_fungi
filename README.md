<img align="right" src="https://github.com/SBWeinstein/Neotoma_fungi/blob/main/github_graphic-01.svg" width="20%">

# This respository contains code associated with "Wild herbivorous mammals (genus *Neotoma*) host a diverse but transient assemblage of fungi"

### Scripts for manuscript are as follows:

**[ITS_seq_preprocess](ITS_seq_preprocess):** Sequence processing to ASV tables and taxonomy calling. Uses QIIME2 framework with DADA2 for denoising and UNITE database for taxonomic calling. Amplicon sequences are available from the [NCBI sequence archive](https://www.ncbi.nlm.nih.gov/sra) under BioProject PRJNA824056. Inputs are paired fastq.gz files. 

**[woodrat_qiime2R.R](woodrat_qiime2R.R):** Convert qiime outputs (table.qza, taxonomy.qza) to phyloseq format, add metadata (gnomex_metadata_ITS_trim.csv), check positive controls and blanks. Removes positive controls to produce physeq3.rds, remove ASVs with <25 read counts and subset to wild samples to produce phyloseq_fungal_w25_27Dec21.rds used in subsequent analyses.

**[Shotgun_bash_commands.sh](Shotgun_bash_commands.sh):** Shotgun metagenomic read processing including quality control, removing host reads, and classification using kraken2.

**[Kraken_outputs_to_R.R](Kraken_outputs_to_R.R):** Creates csv file with proportion of reads assigned to fungi (and other groups) from kraken2 outputs. Produces kraken_bigDB_Conf0.05_outputs_30Mar22.csv.

**[Figure1_correlations.R](Figure1_correlations.R):** Compare the number of ITS reads, 18s copies per ng DNA, and percent of metagenomic reads assigned to fungi. Uses phyloseq formatted ITS amplicon data  from woodrat_qiime2R.R (physeq3.rds), Fungiquant results from FungiQuant_Results_wMetadata_All_13Dec21.csv, and kraken2 taxonomy reports combined into kraken_bigDB_Conf0.05_outputs_30Mar22.csv. Produces Figure 1.

**[NOAA_precip_data.R](NOAA_precip_data.R):** Download precipitation data from NOAA 

**[woodrat_ITS_diversity.R](woodrat_ITS_diversity.R):** Calculates prevalence of fungal taxa and examines factors predicting fungal diversity, quantity and composition in wild rats. Script uses phyloseq formatted files produced by woodrat_qiime2R.R (phyloseq_fungal_w25_27Dec21.rds), precipitation data from NOAA_precip_data.R, compiled data on the ecology of identified ITS amplicons (fungal_taxa_counts_25eco_wild_24May21.csv), Fungiquant results, and wild woodrat diet data (from Weinstein et al 2021, see https://github.com/SBWeinstein/Neotoma2021). Generates Figures 2A, 3A.

**[Ecology_ASVs_Reads.R](Ecology_ASVs_Reads.R)** Uses phyloseq formatted amplicon data (physeq3.rds) and compiled data (fungal_taxa_counts_25eco_wild_24May21.csv) on the ecology of all identified fungal ASVs to examine the distribution of ecological niches and trophic modes across reads and ASVs.  Produces Figures 2B, 2C

**[MRM_models.R](MRM_models.R):** Script uses rarefied phyloseq formatted files produced in woodrat_ITS_diversity.R to examine how host phylogeny, geography and diet contribute to mycobiome composition in wild individuals. To restrict analyses to plant associated fungi requires fungal ecology data (fungal_taxa_counts_25eco_wild_24May21.csv). Analyses also require woodrat phylogeny with branch lengths from Matoqc 2007 (Neotoma_phylosymtree.nwk), coordinates of sampling sites (Site_lat_long.csv), and diet data from from Weinstein et al 2021.  Generates Figure 3b and Table 3.

**[Wild_cap_comps.R](Wild_cap_comps.R):** Script uses phyloseq formatted files (physeq3.rds) produced by woodrat_qiime2R.R and Fungiquant results (FungiQuant_Results_wMetadata_All_13Dec21.csv) to examine how fungal quantities, ITS read counts and ASV diversity differ between captive and wild individuals. Generates Figures 4 and 5.
