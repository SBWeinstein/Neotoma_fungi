#Classify shotgun metagenomic reads from wild woodrats
#bash command list

#folder containing reads in fastq.gz format from 66 wild woodrat samples 	
export ReadDir=/Metagenomics/Neotoma_sequences/Phylo	

#path to the directory for outputs
OutDir=/Metagenomics/	

WorkDir=mgs_phylo66

##################################
### Run FastP  
mkdir -p $WorkDir/reads/filt
mkdir -p $WorkDir/reads/fastp_reports

#run fastp for batch of folders, each with the forward and reverse reads for one sample
for dir in Neotoma_sequences/Phylo/*; do   
	SampleID=$( ls $dir | head -1 | rev | cut -c22- | rev)
	mkdir -p $WorkDir/reads/filt/${SampleID}
	fastp --in1 $dir/*_L001_R1_001.fastq.gz --in2 $dir/*_L001_R2_001.fastq.gz   --out1 $WorkDir/reads/filt/${SampleID}/${SampleID}_1.fastq.gz  --out2 $WorkDir/reads/filt/${SampleID}/${SampleID}_2.fastq.gz  --detect_adapter_for_pe  --cut_front  --cut_tail  --cut_window_size 4  --cut_mean_quality 20  --length_required 60 --json $WorkDir/reads/fastp_reports/$SampleID.json  --html $WorkDir/reads/fastp_reports/$SampleID.html
done

##############################################
### Remove reads that map to the Neotoma lepida genome
## Index host genomes
bowtie2-build -f /Metagenomics/Neotoma_lepida.Arrow.fasta.masked N_lepida  

mkdir -p $WorkDir/reads/nohost

#remove host reads, each sample is a folder with filtered F,R files 
for dir in $WorkDir/reads/filt/*; do 
	SampleID=$( basename $dir )
	mkdir -p $WorkDir/reads/nohost/${SampleID}

bowtie2 --threads 24  -1 $dir/${SampleID}_1.fastq.gz -2 $dir/${SampleID}_2.fastq.gz -x N_lepida --end-to-end --un-conc-gz $WorkDir/reads/nohost/"${SampleID}_"noNlep.fastq.gz > /dev/null

	mv $WorkDir/reads/nohost/"${SampleID}_"noNlep.fastq.1.gz $WorkDir/reads/nohost/${SampleID}/"${SampleID}_"noHost_R1.fastq.gz 
	mv $WorkDir/reads/nohost/"${SampleID}_"noNlep.fastq.2.gz $WorkDir/reads/nohost/${SampleID}/"${SampleID}_"noHost_R2.fastq.gz
done

###################
#Classify reads using Kraken2 using PlusPFP or BigDB, with varying confidence levels
mkdir -p $WorkDir/reads/kraken2

#run kraken2, make kraken2 folder with a .report and .out file for each sample
for dir in $WorkDir/reads/nohost/*; do 
	SampleID=$( basename $dir )
kraken2   --db Kraken2_DB --confidence 0.05 --threads 24 --output $WorkDir/reads/kraken2/${SampleID}_kraken2.out  --report $WorkDir/reads/kraken2/${SampleID}_kraken2.report --paired $dir/"${SampleID}_"noHost_R1.fastq.gz $dir/"${SampleID}_"noHost_R2.fastq.gz 
done

##########################################################
#download Krakentools scripts to extract fungi from kraken reports
wget https://github.com/jenniferlu717/KrakenTools/archive/master.zip 
unzip master.zip -d KrakenTools

#make scripts runnable
chmod +x /Metagenomics/KrakenTools/KrakenTools-master/kreport2krona.py
cd /Metagenomics/KrakenTools/KrakenTools-master/

#can convert to Krona or Metaphlan formats
#krona output needs intermediate ranks to see fungi as a rank and number is reads assigned to that taxonomy (does not include children)
#metaphlan output has fungi as a major rank, parent amounts include children

#use krakentools to extract Neo reads from bigDB classified samples.  blast to see if reasonable (they're not)
#directory of Neo reads (R1, R2) from each sample 
mkdir -p /Metagenomics/mgs_phylo66/Kraken_bigDB/kraken_bg_DB_outputs_Neo

#path to the directory for outputs
OutDir=/Metagenomics/mgs_phylo66/Kraken_bigDB/kraken_bg_DB_outputs_Neo

#path to reads
export Reads=/Metagenomics/mgs_phylo66/reads/nohost

#path to kraken outputs
export Kraken_out=/Metagenomics/mgs_phylo66/Kraken_bigDB/Mycobiome_Slurm_outputs_bigDB

#run extract_kraken_reads script for each sample
for dir in $Reads/*; do 
	SampleID=$( basename $dir )
./extract_kraken_reads.py -k $Kraken_out/${SampleID}_kraken2.out -r $Kraken_out/${SampleID}_kraken2.report -s1 $Reads/${SampleID}/${SampleID}_noHost_R1.fastq.gz -s2 $Reads/${SampleID}/${SampleID}_noHost_R2.fastq.gz -o $OutDir/"${SampleID}_"Neo_R1.fq -o2 $OutDir/"${SampleID}_"Neo_R2.fq -t 29007 --include-children
done

#####get number of fungal families and genera from kraken report files.  
#convert to krona format to pull out unique families and genera
./kreport2krona.py -r /Metagenomics/mgs_phylo66/bigDB_index_kraken2.report -o index_bigDB_krona_IM_index --intermediate-ranks

#get families (start with f__)
grep "x__Fungi"   /Metagenomics/KrakenTools/KrakenTools-master/index_bigDB_krona_IM_index |  grep -v "^0" |grep -o '\bf__\w*' |  uniq | wc -l
#825

#get gen
grep "x__Fungi"   /Metagenomics/KrakenTools/KrakenTools-master/index_bigDB_krona_IM_index |  grep -v "^0" |grep -o '\bg__\w*' |  uniq | wc -l
#6250

################################
#to get numbers of reads classfied as fungi, plants, bacteria, convert to metaphlan format, then read files into R to make single CSV file
#folder containing kraken2 reports from big DB or pfp DB runs (66 phylo samples)		
export ReportDir=/Metagenomics/mgs_phylo66/mycobiome_kraken2_pfpDB_outputs/kraken2_conf0.5_outputs

#path to the directory for outputs
OutDir=/Metagenomics/mgs_phylo66/mycobiome_kraken2_pfpDB_outputs/kraken2_conf0.5_outputs/meta_outputs

mkdir $OutDir

#make new folder of converted mpa format kraken reports
for file in $ReportDir/*.report; do
	 f="$(basename -- $file)"
	 ./kreport2mpa.py -r $file -o $OutDir/${f}  --no-intermediate-ranks
done

