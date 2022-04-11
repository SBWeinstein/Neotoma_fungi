#!/bin/bash

# Requires: 
# - qiime2-2020.2
# - ITSxpress plugin in qiime2. Version must be greater than 2019.10  when "reverse primer" was implemented in the plugin.
# - Paired input .fastq.gz files. This script assumes file pair names end with "R1_001.fastq.gz" and "R2_001.fastq.gz".
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
############
 
cd ${RAWDIR}

# Make manifest file on the fly, place in metadata file in project directory
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > ${WRKDIR}/metadata/manifest_${RUNID}.txt
for f in *R1_001.fastq.gz
do
  SAMPLEID=`basename ${f%%_*}`
    echo -e "${SAMPLEID}\t${PWD}/${f}\t${PWD}/${f%R1_001.fastq.gz}R2_001.fastq.gz" >> ${WRKDIR}/metadata/manifest_${RUNID}.txt
done

cd ${SCRATCH}

ManifestFile=${WRKDIR}/metadata/manifest_${RUNID}.txt

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ${ManifestFile} --output-path inseq_PE_demux.qza --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize --i-data inseq_PE_demux.qza --o-visualization inseq_PE_demux.qzv

qiime itsxpress trim-pair-output-unmerged --p-reversed-primers --p-taxa F --i-per-sample-sequences inseq_PE_demux.qza --p-region ITS2 --p-threads ${Nproc} --o-trimmed inseq_PE_demux_ITSx.qza --verbose

qiime demux summarize --i-data inseq_PE_demux_ITSx.qza --o-visualization inseq_PE_demux_ITSx.qzv

qiime dada2 denoise-paired \
--i-demultiplexed-seqs inseq_PE_demux_ITSx.qza \
--p-trunc-len-f 0 \
--p-trunc-len-r 0 \
--o-table table.qza \
--o-representative-sequences rep_set.qza \
--o-denoising-stats dada2_DenoiseStats.qza \
--verbose

qiime feature-table tabulate-seqs \
--i-data rep_set.qza \
--o-visualization rep_set.qzv

qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv

qiime feature-classifier classify-sklearn \
--i-classifier ${TrainedClassifierArtifact} \
--i-reads rep_set.qza \
--o-classification taxonomy.qza \
--p-n-jobs ${Nproc}

qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy.qza \
--o-visualization taxbarplots_AllSamplesNoFilter.qzv \
--m-metadata-file ${MetadataMapFile}

# Copy key outputs back to wrkdir:
cp *.qza ${WRKDIR}; cp *.qzv ${WRKDIR}/q2_viz
