#! /bin/bash

##Pipeline to analyse 16S rDNA from bacterial communities of the octocoral Pacifigorgia Cairnsi
##The steps were performe using QIIME2 v2021.4
##By Sandra Montaño Salazar

#Import data

qiime tools import --type 'SampleDta[PairedEndSequencesWithQuality]' --input-path seccoral/manifest.csv --input-format PairedEndFastqManifestPhred33 --output-path seccoral/demux-paired-end.qza

#Vizualize the data

qiime demux summarize --i-data 16S-paired-end.qza --o-visualization 16S-paired-end.qzv

#Merge the forward and reverse

qiime vsearch join-pairs --i-demultiplexed-seqs demux.qza --o-joined-sequences demux-joined.qza --p-minmergelen 248 --p-maxmergelen 258  --p-maxdiffs 10

#Quality filter based on q-score

qiime quality-filter q-score-joined --i-demux ~/demux-joinedt2.qza --o-filtered-sequences ~/demux-joined-final-filtered-q20.qza --o-filter-stats ~/demux-joined-filter-stats.qza --p-min-quality 20

#Sequences processed using Deblur to determine the ASVs

qiime deblur denoise-16S --i-demultiplexed-seqs ~/demux-joined-final-filtered-q20.qza --p-trim-length 250 --p-sample-stats --o-representative-sequences ~/rep-seqs-deblur-250final.qza --o-table ~/table-deblur250-final.qza --o-stats ~/deblur-stats-250final.qza

#For assign taxonomy, the first step is obtain silva classifier for V4 region

wget -O "silva_132_99_515F_806R_nb_classifier.qza" "https://data.qiime2.org/2018.8/common/silva-132-99-515-806-nb-classifier.qza"

#The second step for taxonomical assign is training the classifier

qiime feature-classifier classify-sklearn --i-classifier silva_132_99_515F_806R_nb_classifier.qza --i-reads rep-seqs-dada2.qza --o-classification taxonomy_silva.qza

#For filter out the mitochondrias, chloroplast and Eukaryota

qiime taxa filter-table --i-table ~/table-deblur250-final.qza --i-taxonomy ~/taxonomy_silva.qza --p-exclude Mitochondria,Chloroplast,Eukaryota --o-filtered-table ~/table-deblur-final-fmce.qza

#Finally, filter the minimum number of samples and minimum frecuency

qiime feature-table filter-features --i-table ~/table-deblur-final-fmce.qza --p-min-frequency 3 --p-min-samples 3 --o-filtered-table ~/table-minsamp3-minfrec3.qza

#To export the tables

qiime tools export ~/table-minsamp3-minfrec3.qza --output-dir ~table-mce-min3

#Convert the table from biom format to txt file

biom convert -i ~/table-mce-min3/feature-table.biom -o ~/table-mce-min3/export-table-mce-min3.txt --to-tsv