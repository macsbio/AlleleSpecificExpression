This folder contains one directory for the initial processing of RNA-sequencing data (ASEReadCounter)
ASEReadCounter is a GATK tool used to count allele specific reads
    The script also utilizes a few other files (reference genome) and tools (picard/jar)
    These are all descried in the .sh script!

The other folder (ASE_analysis) contains the actual analysis of the allele specific read counts
    ASEReadCounter_OutputProcessing takes the output of ASEReadCounter and turns it into several .csv files
    These files contain read counts per locus, total counts, and ratios for further analysis by using 
    the scripts in that folder
