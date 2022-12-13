# AlleleSpecificExpression

THIS PIPELINE CAN BE USED TO ANALYSE RNA SEQUENCING DATA TO FIND IMBALANCED ALLELES THAT FUNCTION AS A PROXY FOR CIS-REGULATORY VARIATION.

RNA PREPROCESSING
Script in bash, "ASEReadCounter.sh".
Alignment to reference genome.
Necessary files and tools are noted in this script.
No GC-correction.
Because this removes reads, and we don't know which ones it removes.
This might change the imbalance scores.
Also, we don't care about total expression, only relative expression on the same locus.

ASE TO R
Get all the files from the bash script in R-readable format.
Combine certain data features and subset others.

ASE ANALYSIS
Calculate ASE scores
Setting a threshold for likely homozygous with misread in RNA-seq.
Individual analysis/population level.
Group comparisons.
