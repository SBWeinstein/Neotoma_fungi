#!/bin/bash

# Requires:
# - qiime2-2020.2
# - ITSxpress plugin in qiime2. Version must be greater than 2019.10  when "reverse primer" was implemented in the plugin.
# - Paired input .fastq.gz files. This script assumes file pair names end with "R1_001.fastq.gz" and "R2_001.fastq.gz".
#      - Note that the paired seq files were still input although ASVs in this script are only called from a single end read. 
# - A trained classifier in QIIME2 artefact format for taxonomic classification.

############ Define variables and directory locations here ###################
# Specify location of the classifier artefact file. A file.
TrainedClassifierArtifact=
# Num of processes for all commands except DADA2. An integer.
Nproc=26
# Define the sequencing run ID. A string.
RUNID=17950R
# Specify location of input paired fastq files. A directory. Only the input fastq files should be present.                              
RAWDIR=
# Specify a scratch directory location. A directory.
SCRATCH=
# Specify the working directory for result files. A directory.
WRKDIR=
# A metatadata map file (optional) for adding metadata to barplots.
MetadataMapFile=
#########################################################################
###########
mkdir -p ${WRKDIR}/metadata
mkdir -p ${WRKDIR}/q2_viz
mkdir -p ${SCRATCH}
###########

cd ${RAWDIR}

# Make manifest file on the fly, place in metadata file in project directory. 
# Notably, we read in both forward and reverse at first even though we only use the forward read in the end. This was because the itsxpress qiime2 plugin command for single end was not working properly at the time, but we can use the 'trim-pair-output-unmerged' command and output the unpaired for single end analyses.
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > ${WRKDIR}/metadata/manifest_${RUNID}_PE.txt
for f in *R1_001.fastq.gz
do
  SAMPLEID=`basename ${f%%_*}`
    echo -e "${SAMPLEID}\t${PWD}/${f}\t${PWD}/${f%R1_001.fastq.gz}R2_001.fastq.gz" >> ${WRKDIR}/metadata/manifest_${RUNID}_PE.txt
done

cd ${SCRATCH}

ManifestFile=${WRKDIR}/metadata/manifest_${RUNID}_PE.txt

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ${ManifestFile} --output-path inseq_PE_demux.qza --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize --i-data inseq_PE_demux.qza --o-visualization inseq_PE_demux.qzv

# Here, we use itsxpress call to use paired end reads, because the single end trimming wass not working, but we output unmerged pairs. 
qiime itsxpress trim-pair-output-unmerged --p-reversed-primers --p-taxa F --i-per-sample-sequences inseq_PE_demux.qza --p-region ITS2 --p-threads ${Nproc} --o-trimmed inseq_PE_demux_ITSx.qza --verbose

qiime demux summarize --i-data inseq_PE_demux_ITSx.qza --o-visualization inseq_PE_demux_ITSx.qzv

# We now output the ITS region and trimmed sequences  (via ITSx) so that we can extract just the forward reads, reimport them, and denoise them alone with dada2 and create a table. 
########
qiime tools export --input-path inseq_PE_demux_ITSx.qza --output-path inseq_PE_demux_ITSx
# From QIIME2: Exported inseq_PE_demux_ITSx.qza as SingleLanePerSamplePairedEndFastqDirFmt to directory inseq_PE_demux_ITSx

# Using the exported MANIFEST file, edit it and create a read1 and read2 only manifest.
# Importantly, because ITSxpress swaps these due to our seq strategy in reverse orientation, the "R1" files exported are actually read2 and vice-versa. Thus, hereafter, I use "for"(forward) for R1 which were originally R2 reads and vice versa. Thes are now forward read in correct orientation relative to standard ITS orientation.
cd inseq_PE_demux_ITSx

echo -e "sample-id\tabsolute-filepath" > ${WRKDIR}/metadata/manifest_${RUNID}_SE_for.txt
for f in *R1_001.fastq.gz
do
  SAMPLEID=`basename ${f%%_*}`
    echo -e "${SAMPLEID}\t${PWD}/${f}" >> ${WRKDIR}/metadata/manifest_${RUNID}_SE_for.txt
done

echo -e "sample-id\tabsolute-filepath" > ${WRKDIR}/metadata/manifest_${RUNID}_SE_rev.txt
for g in *R2_001.fastq.gz
do
  SAMPLEID=`basename ${g%%_*}`
    echo -e "${SAMPLEID}\t${PWD}/${g}" >> ${WRKDIR}/metadata/manifest_${RUNID}_SE_rev.txt
done

cd ../

# Now, re-import and summarize
qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path ${WRKDIR}/metadata/manifest_${RUNID}_SE_for.txt --output-path inseq_SE_for_demux_ITSx.qza --input-format SingleEndFastqManifestPhred33V2
qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path ${WRKDIR}/metadata/manifest_${RUNID}_SE_rev.txt --output-path inseq_SE_rev_demux_ITSx.qza --input-format SingleEndFastqManifestPhred33V2
qiime demux summarize --i-data inseq_SE_for_demux_ITSx.qza --o-visualization inseq_SE_for_demux_ITSx.qzv
qiime demux summarize --i-data inseq_SE_rev_demux_ITSx.qza --o-visualization inseq_SE_rev_demux_ITSx.qzv
########

# R1(original R2)/forward read denoising and table
########
qiime dada2 denoise-single \
--i-demultiplexed-seqs inseq_SE_for_demux_ITSx.qza \
--p-trunc-len 0 \
--o-table table_SE_for.qza \
--o-representative-sequences rep_set_SE_for.qza \
--o-denoising-stats dada2_DenoiseStats_SE_for.qza

qiime feature-table tabulate-seqs \
--i-data rep_set_SE_for.qza \
--o-visualization rep_set_SE_for.qzv

qiime feature-table summarize \
--i-table table_SE_for.qza \
--o-visualization table_SE_for.qzv

qiime feature-classifier classify-sklearn \
--i-classifier ${TrainedClassifierArtifact} \
--i-reads rep_set_SE_for.qza \
--o-classification taxonomy_SE_for.qza \
--p-n-jobs ${Nproc}

qiime metadata tabulate \
--m-input-file taxonomy_SE_for.qza \
--o-visualization taxonomy_SE_for.qzv

qiime taxa barplot \
--i-table table_SE_for.qza \
--i-taxonomy taxonomy_SE_for.qza \
--o-visualization taxbarplots_AllSamplesNoFilter_SE_for.qzv \
--m-metadata-file ${MetadataMapFile}
########

# R2(original R1)/reverse read denoising and table
########
qiime dada2 denoise-single \
--i-demultiplexed-seqs inseq_SE_rev_demux_ITSx.qza \
--p-trunc-len 0 \
--o-table table_SE_rev.qza \
--o-representative-sequences rep_set_SE_rev.qza \
--o-denoising-stats dada2_DenoiseStats_SE_rev.qza

qiime feature-table tabulate-seqs \
--i-data rep_set_SE_rev.qza \
--o-visualization rep_set_SE_rev.qzv

qiime feature-table summarize \
--i-table table_SE_rev.qza \
--o-visualization table_SE_rev.qzv

qiime feature-classifier classify-sklearn \
--i-classifier ${TrainedClassifierArtifact} \
--i-reads rep_set_SE_rev.qza \
--o-classification taxonomy_SE_rev.qza \
--p-n-jobs ${Nproc}

qiime metadata tabulate \
--m-input-file taxonomy_SE_rev.qza \
--o-visualization taxonomy_SE_rev.qzv

qiime taxa barplot \
--i-table table_SE_rev.qza \
--i-taxonomy taxonomy_SE_rev.qza \
--o-visualization taxbarplots_AllSamplesNoFilter_SE_rev.qzv \
--m-metadata-file ${MetadataMapFile}
########

# Copy key outputs back to wrkdir:
cp *.qza ${WRKDIR}; cp *.qzv ${WRKDIR}/q2_viz
