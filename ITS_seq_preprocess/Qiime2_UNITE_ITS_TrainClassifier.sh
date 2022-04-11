#!/bin/bash

# Requires:
# - qiime2-2020.2
# - URL of UNITE database set. 

# Method: Pulls reference sets from UNITE database, unzip, import to qiime object and train classifier. Likely only works with the specific UNITE database version and naming scheme in that.
# See UNITE donwloads page: https://unite.ut.ee/repository.php
# Here, we used the "dynamic" species hypothesis clustering set.
# I also use the provided set that was trimmed by ITSx. The download includes a "developer" folder with the untrimmed set that has flanking SSU, LSU seqs not trimmed.

############### Define variables ############################
# Specify the directory to keep the reference sequences and artefact files created
REFSEQDIR=
# URL to UNITE seqs. The original URL used is provided.
WebHostedRef=https://files.plutof.ut.ee/public/orig/98/AE/98AE96C6593FC9C52D1C46B96C2D9064291F4DBA625EF189FEC1CCAFCF4A1691.gz
# The version date or number to append to files. The original used is provided.
VersionDate=2020.02.04
##############################################################


cd ${REFSEQDIR}
wget -O DownloadedRefSeqs.tar.gz ${WebHostedRef}
tar -xzf DownloadedRefSeqs.tar.gz
cd sh_qiime_release*
InputRepSeqs=`ls sh_refs_qiime*_dynamic_*.fasta`
InputTaxa=`ls sh_taxonomy*_dynamic_*.txt`

qiime tools import \
--input-path ${InputRepSeqs} \
--output-path ../${InputRepSeqs}.qza \
--type 'FeatureData[Sequence]'

qiime tools import \
--input-path ${InputTaxa} \
--type FeatureData[Taxonomy] \
--input-format HeaderlessTSVTaxonomyFormat \
--output-path ../${InputTaxa}.qza

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ../${InputRepSeqs}.qza \
--i-reference-taxonomy ../${InputTaxa}.qza \
--o-classifier ../naive_bayes_classifier_${InputRepSeqs%.fasta}_${VersionDate}.qza

cd ../
rm -R sh_qiime_release*
