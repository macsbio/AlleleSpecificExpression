#!/bin/bash

# This contains all the lines of script needed to generate ASE variant counts using
#  the ASEReadCounter function of GATK 3.8.1
# Newer versions of ASEReadCounter (GATK >4.0) apparently stopped working. 
#  Hence until fixed, will rely on older version (November 2017, so not that old).
#
# Assumes script is run from within BAM file folder, that contains the folders:
#  	- ASE_results (were the ASReadCounter CSVs will end up)
#	- ref (containing all FASTA derived files and the variant reference VCF)
#   - tools (containing picard 2.18.11, GATK 3.8.1, faToTwoBit 16-07-2018)

# 0.  Preliminaries
# ORIGINAL FASTA REQUIRED TO WHICH BAM FILES ARE ALLIGNED
# Needs to be converted to 2bit (for GC correction) and a DICT file needs 
# to be created (for GATK):
#	./faToTwoBit GRCh38.fa GRCh38.2bit
#	java -jar picard.jar CreateSequenceDictionary REFERENCE=GRCh38.fa OUTPUT=GRCh38.dict
#
# Sometimes need to change VCF file to match chromosome notation of BAM and FASTA
# e.g. to change 'chr1' etc notation to '1' etc notation (and here chrM to MT).
#	awk '{gsub(/^chrM/,"chrMT"); print}' 1000G_phase1.snps.high_confidence.hg38.vcf | \
#	awk '{gsub(/^chr/,""); print}' > 1000G_phase1.snps.high_confidence.chrRecode.hg38.vcf
#	bgzip 1000G_phase1.snps.high_confidence.chrRecode.hg38.vcf
#	tabix -p vcf 1000G_phase1.snps.high_confidence.chrRecode.hg38.vcf.gz
# To install tabix (not sure if biobuntu specific): 
#	apt install tabix

# Some variables to set:
FASTA=GRCh38

# MAKE SURE TO CHECK TO MATCH FILENAME PATTERNS IN YOUR DATASET
# Start of loop:
for f in $PWD/*.bam ; do
  FILE=`basename ${f}`
  FILE=${FILE%.*}

# 2.a Reorder BAM file to match exactly with FASTA reference file contigs. This is ok as long as
#     the main chromosomes (i.e. autosomes and M, X, Y) contigs are identical between the 
#     FASTA reference file and the BAM files. In this case I have manually checked.
#     Do not care about the random, alternative and uncharacterized contigs, as they
#     do not matter in the ASE variant counting (since almost all SNPs will be on the
#     main chromosomes, and any that are on the 'junk' I will ignore):
java -jar $PWD/tools/picard.jar ReorderSam INPUT=${FILE}.bam \
 OUTPUT=${FILE}.r.bam CREATE_INDEX=true REFERENCE=$PWD/ref/${FASTA}.fa \
 ALLOW_INCOMPLETE_DICT_CONCORDANCE=true ALLOW_CONTIG_LENGTH_DISCORDANCE=true
 
 # The error it gives is due to more contigs in FASTA file compared to BAM.
 #   See: https://goo.gl/YE1k6D
 
# 2.b A better solution would be to extract all main autosomes and M,X,Y from BAMs,
#     i.e. filter out all the uncharacterized, alternative and random contigs:
#bash -c "samtools view -b ${FILE}.r.bam \
# chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 \
# chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrM chrX chrY \
# > ${FILE}.f.bam"
# But this isn't working at the moment, generates errors later on in step 2

# 3.  Apparently the below is needed to generate additional read group characteristics
#     and generate a new index file. The former is required by GATK:
java -jar $PWD/tools/picard.jar AddOrReplaceReadGroups \
       I=${FILE}.r.bam \
       O=${FILE}.rg.bam \
       RGID=dcm \
       RGLB=dcm \
       RGPL=illumina \
       RGPU=dcm \
       RGSM=dcm \
       CREATE_INDEX=true
rm -f ${FILE}.r.ba*

# 4. ASE variant calling:
# Does not work currently: 
#  ./gatk-4.0.7.0/gatk ASEReadCounter -I ${FILE}.reordered.rg.bam \
#  -O ${FILE}_ASEresults.csv -V 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
#  -R ${FASTA}.fa -min-depth 10 -mmq 20 -mbq 5

# Using GATK 3.8.1 works: (45 minutes)
java -jar $PWD/tools/GenomeAnalysisTK.jar -R $PWD/ref/${FASTA}.fa -T ASEReadCounter \
 -I ${FILE}.rg.bam -o $PWD/ASE_results/${FILE}_ASEresults.csv -sites \
 $PWD/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
 --minDepthOfNonFilteredBase 10 --minMappingQuality 20 --minBaseQuality 5 \
 -U ALLOW_N_CIGAR_READS --outputFormat CSV
rm -f ${FILE}.rg.ba*

done

echo $SECONDS

# END
