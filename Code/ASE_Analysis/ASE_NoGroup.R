#=============================================================================#
# DCM_MUMC_ASE_analysis.R                                                     #
#                                                                             #
# Version: 4.0                                                                #
# Date:  Fall 2020 - Spring 2021                                              #
# Author: Michiel Adriaens, PhD; MaCSBio, Maastricht University               #
# Additional: Daan van Beek; MaCSBio                                          #
# History:                                                                    #
#  1.0: Creation                                                              #
#  2.0: Updated for new data and objects                                      #
#  2.1: Removed all ideas, thoughts, etc. from start of script and merged     #
#       into notes (Onenote)                                                  #
#  3.0: Likely homozygote filtering using WES data                            #
#  3.1: Updated and streamlined SNP <-> gene functionality (fast now!)        #
#  3.2: Bug fixing in SNP <-> gene lookup functions, added LD lookup function #
#  3.3: Updated with semi-complete data (n = 80) and phen. clustering         #
#  4.0: Update with n = 87 data, absolute deviation, new plotting             #
#                                                                             #
#=============================================================================#

#-----------------------------------------------------------------------------#
# 1. Initialization
#-----------------------------------------------------------------------------#

# Some information about necessary files and column headers
# Clinical file should contain an SID header and Group (for cluster number or case/control) header
# EnsgID header --> ensemble_gene_id
# HGNC header --> hgnc_symbol
# SNP id header --> rsid
# Stat values header --> pvalue, qvalue
# Other input files are created by running ASEReadCounter_OutputProcessing
#   on the output from the ASEReadCounter pre-processing in bash/Linux
# If interested in comparisons with known related genes
#   change the name of DCM_genes.txt
#   It should have the header ensembl_gene_id and hgnc_symbol

#Clear working space
rm(list = ls(all.names = T))

ASEReadCounterFile <- "ASEReadCounterOutput_Processed_DCM.RData"
genes_of_interest <- "DCM_genes.txt"

# If TRUE, provide a filename for the genotype data
geno <- TRUE
genotypedata <- "wesHetDataWES_HG38HG19.RData"

# Load packages
if(!require(rstudioapi)){
  install.packages("rstudioapi")
  library(rstudioapi)}
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)}
if(!require(ggsci)){
  install.packages("ggsci")
  library(ggsci)}
if(!require(dplyr)) {
  install("dplyr")
  library(dplyr) }
if(!require(tidyr)) {
  install("tidyr")
  library(tidyr) }
if(!require(rJava)) {
  install.packages("rJava")
  library(rJava) }
if(!require(reshape2)){
  install.packages("reshape2")
  library(reshape2)}
if(!require(readxl)){
  install.packages("readxl")
  library(readxl)}
if(!require(pROC)){
  install.packages("pROC")
  library(pROC)}
if(!require(BiocManager)) {
  install.packages("BiocManager") }
if(!require(biomaRt)){
  BiocManager::install("biomaRt")
  library(biomaRt)}
if(!require(qvalue)){
  BiocManager::install("qvalue")
  library(qvalue)}
if(!require(limma)) {
  BiocManager::install("limma")
  library(limma) }
if(!require(edgeR)) {
  BiocManager::install("edgeR") }
library(edgeR)
if(!require(RColorBrewer)) {
  BiocManager::install("RColorBrewer") }
library(RColorBrewer)
if(!require(ggrepel)) {
  BiocManager::install("ggrepel") }
library(ggrepel)
if(!require(pcaMethods)) {
  BiocManager::install("pcaMethods") }
library(pcaMethods)
if(!require(org.Hs.eg.db)){
  BiocManager::install("org.Hs.eg.db") }
library(org.Hs.eg.db)
if(!require(topGO)){
  BiocManager::install("topGO") }
library(topGO)
if(!require("RCy3")) {
  BiocManager::install("RCy3") }
library(RCy3)

# Set directories
DATA.DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
if(!dir.exists(paste0(DATA.DIR, "/Output"))) {
  dir.create(paste0(DATA.DIR, "/Output"))
} 
RES.DIR <- paste0(DATA.DIR, "/Output") 

options(stringsAsFactors = F)

# Load processed ASEReadCounter output data
setwd(DATA.DIR)
load(ASEReadCounterFile)

# Check if all orders of columns and rows are equal across the dataframes
all(dosageData[, 1] == altRatioData[, 1]) # TRUE
all(colnames(dosageData) == colnames(altRatioData)) # TRUE
all(dosageData[, 1] == totalCountData[, 1]) # TRUE
all(colnames(dosageData) == colnames(totalCountData)) # TRUE
all(dosageData[, 1] == refCountData[, 1]) # TRUE
all(colnames(dosageData) == colnames(refCountData)) # TRUE
all(dosageData[, 1] == altCountData[, 1]) # TRUE
all(colnames(dosageData) == colnames(altCountData)) # TRUE
all(dosageData[, 1] == refRatioData[, 1]) # TRUE
all(colnames(dosageData) == colnames(altCountData)) # TRUE

# Restructure data
rwNms <- dosageData[, 1]
dosageData   <- dosageData[, -1]
altRatioData <- altRatioData[, -1]
totalCountData <- totalCountData[, -1]
refCountData <- refCountData[, -1]
refRatioData <- refRatioData[, -1]
altCountData <- altCountData[, -1]
rownames(altRatioData) <- rownames(dosageData) <- rownames(totalCountData) <- rwNms
rownames(refCountData) <- rownames(altCountData) <- rownames(refRatioData) <- rwNms

# Calculate the ASE score according to what is common practice:
# Ranging from 0.5 for no imbalance to 1 for full imbalance
# This is achieved by taking the alternative allele count divided by the
# total allele count (a value between 0 and 1), minus 0.5 in absolute value, + 0.5.
absData <- (abs(altRatioData - 0.5) + 0.5)

# This part is only data cleaning for OUR PROJECT
colnames(altRatioData) <- colnames(dosageData) <- colnames(totalCountData) <- colnames(absData) <- gsub("SID.","", colnames(absData))
colnames(altRatioData) <- colnames(dosageData) <- colnames(totalCountData) <- colnames(absData) <- gsub("D","", colnames(absData))

colnames(altCountData) <- colnames(refCountData) <- colnames(altRatioData)
colnames(refRatioData) <- colnames(altCountData)

#-----------------------------------------------------------------------------#
# 2. Filter the data based on likely homozygosity
#-----------------------------------------------------------------------------#

# Derive threshold for calling homozygotes from ASE ratios
# Can be done using genotype data (geno == TRUE)
# Or with other methodologies
# The one we presented is based on including only loci with >= 10 reads
#                                                                 per allele
#-----------------------------------------------------------------------------#
if (geno == TRUE) {
  load('wesHetDataWES_HG38HG19.RData')
  colnames(wesHetData)[-1] <- gsub("\\.hg19", "", colnames(wesHetData)[-1])
  
  # Take only loci with ASE scores
  absData.sel <- absData
  absData.sel <- absData.sel[apply(absData.sel, 1, function(x) any(!is.na(x))),]
  
  # Define SNPs shared between WES and ASE datasets
  sharedSnps <- intersect(wesHetData[, 1], rownames(absData.sel))
  
  # Keep only SNPs where alt and ref are defined the same (less than 0.2% is lost)
  sharedSnps <- sharedSnps[which(
    annData[match(sharedSnps, annData$id), 'refAllele'] ==
      wesAnnData[match(sharedSnps, wesAnnData$id), 'refAllele'] &
      annData[match(sharedSnps, annData$id), 'altAllele'] ==
      wesAnnData[match(sharedSnps, wesAnnData$id), 'altAllele'])]
  
  # Define shared samples
  commonSamples <- intersect(colnames(absData.sel), colnames(wesHetData))
  
  # Align WES and ASE datasets for easy comparison 
  wesHetData.sel <- wesHetData[, -1]
  rownames(wesHetData.sel) <- wesHetData[, 1]
  absData.sel <- absData.sel[sharedSnps, commonSamples]
  totalCountData.sel <- totalCountData[sharedSnps, commonSamples]
  
  # Recode WES data allele dosage to text (for easy plotting)
  wesHetData.sel <- wesHetData.sel[sharedSnps, commonSamples]
  wesHetData.sel[wesHetData.sel == 1] <- "het"
  wesHetData.sel[wesHetData.sel == 2] <- "altHom"
  wesHetData.sel[is.na(wesHetData.sel)]  <- "refHom" # not correct at all*
  
  # Create density of ASE as a function of zygosity
  # Uses all rsIDs and all sIDs overlapping in WES and RNA-seq data
  setwd(RES.DIR)
  absData.sel1 <- absData.sel
  absData.sel1$rsid <- rownames(absData.sel1) 
  wds <- wesHetData.sel
  wds$rsid <- rownames(wds)
  abstest <- pivot_longer(absData.sel1, c(1:(ncol(absData.sel1)-1)), values_to = "ASE", names_to = "SID")
  wdstest <- pivot_longer(wds, c(1:(ncol(wds)-1)), values_to = "ZYG", names_to = "SID")
  
  # Recode homozygotes to 1 and heterozygotes to 0
  wdstest$ZYG[wdstest$ZYG=="altHom"] <- 1
  wdstest$ZYG[wdstest$ZYG=="refHom"] <- 1
  wdstest$ZYG[wdstest$ZYG=="het"] <- 0
  
  # Join ASE scores with zygosity
  wdstest$ZYG <- as.factor(wdstest$ZYG)
  homtest <- full_join(abstest, wdstest, by = c("rsid","SID"))
  homtest <- na.omit(homtest)
  
  # Create ROC and calculate threshold using Youden's J statistic
  threshTest <- list()
  homtesthom <- homtest[homtest$ZYG==1,]
  homtesthet <- homtest[homtest$ZYG==0,]
  
  for (i in 1:100) {
    homtesthomCurr <- homtesthom[sample(nrow(homtesthom), size = 167329),]
    homtestCurr <- rbind(homtesthomCurr,homtesthet)
    roc1Curr <- roc(homtestCurr$ZYG,homtestCurr$ASE,plot=T)
    threshTest[[i]] <- coords(roc1Curr, x="best", ret="threshold", 
                              best.method="youden")
  }
  
  roc1 <- roc(homtest$ZYG,homtest$ASE,plot=T)
  thresh <- as.numeric(coords(roc1, x="best", ret="threshold", 
                              best.method="youden"))
  
  homtest$ZYG <- as.character(homtest$ZYG)
  homtest$ZYG[homtest$ZYG==1] <- "Homozygote"
  homtest$ZYG[homtest$ZYG==0] <- "Heterozygote"
  
# Density plot for ASE values in all genes
density <- ggplot(homtest, aes(x = ASE)) + 
  geom_density(aes(color = ZYG, fill = ZYG, y=after_stat(scaled)),
               linewidth = 0.8, alpha = 0.2) +
  labs(y = "Density", x = "ASE") +
  scale_color_lancet(name = "Zygosity") +
  scale_fill_lancet(name = "Zygosity") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.line = element_line(colour = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    title = element_text(face = "bold"),
    axis.text.x = element_text(size = 8, margin = margin(l=0,r=0,t=10,b=0)),
    axis.text.y = element_text(size = 8, margin = margin(l=0,r=10,t=0,b=0)),
    axis.title.x = element_text(size = 10, margin = margin(l=0,r=0,t=10,b=0)),
    axis.title.y = element_text(size = 10, margin = margin(l=0,r=10,t=0,b=0)),
    axis.ticks = element_blank(), 
  )
tiff("ASE_Density.tif", res = 300, width = 12, height = 9, units = "cm")
print(density)
dev.off()
  # Remove homozygote samples (based on genotyping data)
  absData.sel[wesHetData.sel=="altHom" | wesHetData.sel=="refHom"] <- NA
  
  # Apply zygosity ASE threshold to non-sequenced samples
  absData[dosageData == 0 | dosageData == 2] <- NA
  absData[absData>thresh] <- NA
  absData1 <- absData
  absData[sharedSnps, commonSamples] <- absData.sel
} else {
  altCountData <- (altRatioData * totalCountData)
  refCountData <- (refRatioData * totalCountData)
  absData[altCountData<10] <- NA
  absData[refCountData<10] <- NA
}

# Not necessary: remove data dat is not needed anymore
rm(list=setdiff(ls(), c("absData","absData1","thresh","annData","DATA.DIR","RES.DIR","totalCountData", "genes_of_interest")))

# Filter out SNPs with only NAs in absData (i.e. not measured or only homozygotes)
selection2 <- which(apply(absData, 1, function(x) all(is.na(x))))
totalCountData <- totalCountData[-selection2, ]
absData <- absData[-selection2, ]
absData1 <- absData1[-selection2, ]

nSamples <- ncol(absData)

# Calculate the median ASE value among all measured loci for binomial model
absSnpMedians <- apply(absData, 2, median, na.rm = T)
IS.MEDIAN = median(absSnpMedians)

# Export all SNPs for future reference
setwd(RES.DIR)
write.table(rownames(absData), file = "allAseSnps.txt", sep = "\t",
            row.names = F, col.names = F, quote = F)

#-----------------------------------------------------------------------------#
# 3. From SNPs to genes: functions and annotation
#-----------------------------------------------------------------------------#

# Create giant table, once, linking all ASE SNPs to all genes
#-----------------------------------------------------------------------------#

# Save latest sample IDs for next step
if (!file.exists("previousIDs.txt")) {
  write.table(list(colnames(absData)), "previousIDs.txt", 
              col.names = F, row.names = F)
}

# IMPORTANT: THIS PART NEEDS TO BE RERUN FOR EVERY CHANGE IN INCLUDED SAMPLES
if (!all(colnames(absData) %in% read.table("previousIDs.txt")$V1)) {
  grch38.snp <<- useMart(biomart = "ENSEMBL_MART_SNP", 
                         dataset = "hsapiens_snp")
  grch38.gene <<- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                          dataset = "hsapiens_gene_ensembl")
  bmSnpData <- getBM(attributes = c("refsnp_id", "ensembl_gene_stable_id"), 
                     filters = "snp_filter", 
                     values = rownames(absData),
                     mart = grch38.snp)
  bmGeneData <- getBM(attributes = c("hgnc_symbol", "external_gene_name", 
                                     "ensembl_gene_id","description"), 
                      filters = "ensembl_gene_id", 
                      values = unique(bmSnpData$ensembl_gene_stable_id), 
                      mart = grch38.gene)
  bmGeneData[bmGeneData[, 1] == "", 1] <- bmGeneData[bmGeneData[, 1] == "", 2]
  colnames(bmSnpData)[2] <- "ensembl_gene_id"
  snp2geneData <- merge(bmSnpData, bmGeneData[, c(1, 3)], 
                        by = "ensembl_gene_id")
  colnames(snp2geneData)[2] <- "rsid"
  setwd(RES.DIR)
  save(list = "snp2geneData", file = "snp2geneData.RData")
  write.table(list(colnames(absData)), "previousIDs.txt", 
              col.names = F, row.names = F)
}

load("snp2geneData.RData")

#-----------------------------------------------------------------------------#
# 5. Individual based prioritization metrics
#-----------------------------------------------------------------------------#
# Population median as expected!

# Binomial test for significance of ASE columns are samples, rows are rsIDs
ase.binom.test <- function(countdata, asedata, expected) {
  asePvalueMatrix <- c()
  for (i in colnames(countdata)) {
    nr  <- countdata[, i]
    absind <- asedata[, i]
    p <- apply(cbind(absind, nr), 1, function(r) {
      if (is.na(r[1])) {
        NA
      } else {
        binom.test(round(r[1] * r[2]), r[2], 
                   expected)$p.value
      }
    })
    asePvalueMatrix <- cbind(asePvalueMatrix, p)
    paste0(i)
  }
  rownames(asePvalueMatrix) <- rownames(asedata)
  colnames(asePvalueMatrix) <- colnames(asedata)
  return(asePvalueMatrix)
}

asePvalueMatrix <- ase.binom.test(countdata = totalCountData,
                                  asedata = absData,
                                  expected = IS.MEDIAN)

# Correct p-values
aseQvalueMatrix <- asePvalueMatrix
aseQvalueMatrix[!is.na(aseQvalueMatrix)] <- 
  qvalue(as.vector(as.matrix(aseQvalueMatrix[!is.na(aseQvalueMatrix)])))$qvalues

# Check distributions of ASE and q-values
absScores <- absData
absScores$rsid <- rownames(absScores)
absScores <- pivot_longer(absScores, c(1:(ncol(absScores)-1)), values_to = "ASE", names_to = "SID")
absScores <- na.omit(absScores)
absScores <- full_join(absScores,snp2geneData)

# Density plot for ASE values in all genes
tiff("ASE_Density_Total.tif", res = 300, width = 85, height = 85, units = "mm")
ggplot(absScores, aes(x = ASE)) +
  geom_density(linewidth = 0.8, alpha = 0.2) +
  labs(y = "Density", x = "ASE") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_fixed(ratio = 0.5) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.line = element_line(colour = "black", linewidth = 1),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "white"),
    title = element_text(face = "bold"),
    axis.text.x = element_text(size = 8, margin = margin(l=0,r=0,t=10,b=0)),
    axis.text.y = element_text(size = 8, margin = margin(l=0,r=10,t=0,b=0)),
    axis.title.x = element_text(size = 10, margin = margin(l=0,r=0,t=10,b=0)),
    axis.title.y = element_text(size = 10, margin = margin(l=0,r=10,t=0,b=0)),
    axis.ticks = element_blank(), 
  )
dev.off()

# CREAT QQ PLOT FOR ALL SNPs WITH THRESHOLD APPLIED ALSO TO WES SNPS
# For a fairer comparison
# Calculate p-values to all SNPs
asePvalueMatrix1 <- ase.binom.test(countdata = totalCountData,
                                   asedata = absData1,
                                   expected = IS.MEDIAN)

# Create matrix and merge with gene annotations
apm1 <- as.data.frame(asePvalueMatrix1)
apm1$rsid <- rownames(apm1)
apm1 <- pivot_longer(apm1, c(1:(ncol(apm1)-1)),names_to = "SID")
colnames(apm1)[3] <- "pvalue"
apm1 <- na.omit(apm1)
apm1 <- apm1[!apm1$pvalue==0,]
snps2genesP <- full_join(apm1,snp2geneData)

rm(apm1, asePvalueMatrix1)
snps2genesP <- snps2genesP[!is.na(snps2genesP$pvalue),]

# Get top 10 known DCM-associated genes and classify each p-value
setwd(DATA.DIR)
genesOI <- read.table(genes_of_interest, header = T)

snps2genesP$top <- "All other genes"
snps2genesP$top[snps2genesP$ensembl_gene_id %in% genesOI$ensembl_gene_id] <- "Robust DCM genes"

snps2genesP$logp <- (-log10(snps2genesP$pvalue))

setwd(RES.DIR)

# Make the plot
qqplot <- ggplot(snps2genesP, aes(sample = logp)) +
  geom_qq(aes(color = top),shape = 1, size = 0.5, distribution = qchisq, dparams = list(df = 1)) +
  scale_color_brewer(name = "Gene status", palette = "Dark2") +
  xlab("Expected -log10(pvalue)") +
  ylab("Observed -log10(pvalue)") +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.line = element_line(colour = "black", linewidth = 1),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "white"),
    title = element_text(face = "bold"),
    axis.text.x = element_text(size = 8, margin = margin(l=0,r=0,t=10,b=0)),
    axis.text.y = element_text(size = 8, margin = margin(l=0,r=10,t=0,b=0)),
    axis.title.x = element_text(size = 10, margin = margin(l=0,r=0,t=10,b=0)),
    axis.title.y = element_text(size = 10, margin = margin(l=0,r=10,t=0,b=0)),
    axis.ticks = element_blank(),
  )
tiff("QQ_plot_logp_SNP.tif", res = 300, width = 175, height = 85, units = "mm")
print(qqplot)
dev.off()

# Link all q-values to ENSGID for later analyses
aqm <- as.data.frame(aseQvalueMatrix)
aqm$rsid <- rownames(aqm)
aqm <- pivot_longer(aqm, c(1:(ncol(aqm)-1)),names_to = "SID")
colnames(aqm)[3] <- "qvalue"
aqm <- na.omit(aqm)
aqm$qvalue[aqm$qvalue==0] <- min(aqm$qvalue[!aqm$qvalue==0])
snps2genesQ <- full_join(aqm,snp2geneData)

# For individual q-value ASE Manhattan plot
# Get annotation data necessary for plot
annData <- annData[,c(1,2,3)]
colnames(annData)[3] <- "rsid"
annData[,1] <- gsub("chr", "", annData[,1])
annData[,1] <- gsub("_.*", "", annData[,1])
annData <- annData[annData$rsid %in% rownames(aseQvalueMatrix),]
annData[annData$chromosome=="X", "chromosome"] <- 23
annData$chromosome <- as.integer(annData$chromosome)
annData$position <- as.integer(annData$position)
annData <- full_join(annData,snp2geneData)

# Function to make Manhattan plot
# Expects data frame of rsid SID and pvalue/qvalue
# Expects data frame of chromosome, position, rsid, hgnc_symbol
# Standard statvalue = "pvalue" can be changed to "qvalue"
aseManhattan <- function(input, 
                         statvalue = "pvalue", 
                         annotation, 
                         gws, 
                         filename) {
  if (statvalue == "qvalue") {
    input <- full_join(input, annotation, by = "rsid")
    input <- input %>% arrange(chromosome, position)
    
    # Create cumulative positions for ordering x-axis
    input$BPcum <- NA
    input$BPcum <- as.numeric(input$BPcum)
    s <- 0
    nbp <- c()
    chrnum <- unique(input$chromosome)
    for (i in 1:length(chrnum)) {
      nbp[i] <- max(input[input$chromosome == chrnum[i],]$position)
      input[input$chromosome == chrnum[i],"BPcum"] <- input[input$chromosome == chrnum[i],"position"] + s
      s <- s + nbp[i]
    }
    
    #Create centering for each chromosome x-axis location
    axis.set <- input %>%
      group_by(chromosome) %>%
      summarize(center = (max(BPcum) + min(BPcum)) / 2)
    ylim <- (-log10(min(input$qvalue)) + 1)
    sig = gws
    
    # ggplot
    manhattanplot <- ggplot(input, aes(x = BPcum, y = -log10(qvalue), 
                                       colour = as.factor(chromosome))) +
      geom_point() +
      geom_point(data = ~head(input[order(input$qvalue),],5), aes(x = BPcum, y = -log10(qvalue)), 
                 color = "#D95F02") +
      geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
      geom_text_repel(aes(label=hgnc_symbol), size = 2,
                      data = (slice_min(input, n = 20, qvalue))) +
      scale_x_continuous(expand = c(0.05,0.05), 
                         breaks = axis.set$center, labels = axis.set$chromosome,
                         guide = guide_axis(check.overlap = TRUE)) +
      scale_y_continuous(expand = c(0,0), limits = c(0,ylim)) +
      scale_color_manual(values = rep(c("#276EBF", "#183059"), 23)) + 
      labs(x = "Chromosome", y = "-log10(q)") +
      theme_minimal() + 
      theme( 
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        title = element_text(face = "bold"),
        axis.title.x = element_text(size = 10, margin = margin(l=0,r=0,t=10,b=0)),
        axis.title.y = element_text(size = 10, margin = margin(l=0,r=10,t=0,b=0)),
        axis.ticks = element_blank()
      )
    tiff(filename = paste0(filename,".tif"), res = 300, width = 175, height = 85, units = "mm")
    print(manhattanplot)  
    dev.off()
    
  } else {
    input <- full_join(input, annotation, by = "rsid")
    input <- input %>% arrange(chromosome, position)
    
    # Create cumulative positions for ordering x-axis
    input$BPcum <- NA
    input$BPcum <- as.numeric(input$BPcum)
    s <- 0
    nbp <- c()
    chrnum <- unique(input$chromosome)
    for (i in 1:length(chrnum)) {
      nbp[i] <- max(input[input$chromosome == chrnum[i],]$position)
      input[input$chromosome == chrnum[i],"BPcum"] <- input[input$chromosome == chrnum[i],"position"] + s
      s <- s + nbp[i]
    }
    
    #Create centering for each chromosome x-axis location
    axis.set <- input %>%
      group_by(chromosome) %>%
      summarize(center = (max(BPcum) + min(BPcum)) / 2)
    ylim <- (-log10(min(input$pvalue)) + 1)
    sig = gws
    
    # ggplot
    manhattanplot <- ggplot(input, aes(x = BPcum, y = -log10(pvalue), 
                                       colour = as.factor(chromosome))) +
      geom_point() +
      geom_point(data = ~head(input[order(input$pvalue),],5), aes(x = BPcum, y = -log10(pvalue)), 
                 color = "#D95F02") +
      geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
      geom_text_repel(aes(label=hgnc_symbol), size = 2,
                      data = (slice_min(input, n = 20, pvalue))) +
      scale_x_continuous(expand = c(0.05,0.05), 
                         breaks = axis.set$center, labels = axis.set$chromosome,
                         guide = guide_axis(check.overlap = TRUE)) +
      scale_y_continuous(expand = c(0,0), limits = c(0,ylim)) +
      scale_color_manual(values = rep(c("#276EBF", "#183059"), 23)) + 
      labs(x = "Chromosome", y = "-log10(p)") +
      theme_minimal() +
      theme( 
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        title = element_text(face = "bold"),
        axis.title.x = element_text(size = 10, margin = margin(l=0,r=0,t=10,b=0)),
        axis.title.y = element_text(size = 10, margin = margin(l=0,r=10,t=0,b=0)),
        axis.ticks = element_blank()
      )
    tiff(filename = paste0(filename,".tif"), res = 300, width = 175, height = 85, units = "mm")
    print(manhattanplot)  
    dev.off()
  }
}

# Make directory
if(!dir.exists("Individual_Manhattan_plots")) {
  dir.create("Individual_Manhattan_plots")
} 
setwd("Individual_Manhattan_plots")

# For all samples, make a Manhattan plot
for (i in unique(aqm$SID)) {
  dataManhattan <- aqm[aqm$SID==i,]
  annDataManhattan <- annData[annData$rsid %in% dataManhattan$rsid,]
  annDataManhattan <- unique(annDataManhattan)
  aseManhattan(input = dataManhattan,
               statvalue = "qvalue",
               annotation = annDataManhattan,
               gws = 0.05,
               filename = paste0(i,"_Manhattan"))
}

setwd(RES.DIR)

# Calculate gene counts
# Create a data frame of all IDs with corresponding lowest qvalues for each gene
genesQsample <- snps2genesQ[,c("SID","qvalue","ensembl_gene_id")]

genesQsample <- genesQsample[order(genesQsample$SID, 
                                   genesQsample$ensembl_gene_id, 
                                   genesQsample$qvalue),]

# Take most significant hit per person per gene
genesQsample <- genesQsample %>% distinct(SID, ensembl_gene_id, .keep_all= TRUE)
genesQsample <- na.omit(genesQsample)

# Count amount of times gene pops up as < 0.05
AseGeneCounts <- as.data.frame(table(genesQsample[genesQsample$qvalue<0.05,"ensembl_gene_id"]))
colnames(AseGeneCounts) <- c("ensembl_gene_id","Frequency")

# Barplot gene frequencies across individuals
genecounts <- ggplot(AseGeneCounts, aes(x = Frequency)) +
  geom_bar(color = "#00468BFF", fill = "#00468BFF", alpha = 0.2) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Number of individuals", 
       y = "Number of shared gene imbalances") +
  theme(
    legend.position = "none",
    legend.title = element_text(size = 8),
    axis.title.y = element_text(margin = margin(l=0,r=10,t=0,b=0)),
    axis.line = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 8, margin = margin(l=0,r=0,t=10,b=0)),
    axis.text.y = element_text(size = 8, margin = margin(l=0,r=10,t=0,b=0)),
    axis.title = element_text(face = "bold", size = 12),
    axis.ticks = element_blank(), 
  )
tiff("Gene_counts_histogram.tif", res = 300, width = 85, height = 85, units = "mm")
print(genecounts)
dev.off()

# Gene frequencies top DCM genes
write.table(full_join(AseGeneCounts[AseGeneCounts$ensembl_gene_id %in% 
                                      genesOI$ensembl_gene_id,],
                      genesOI)[,c(1,3,2)], 
            "Gene_counts_topDCM.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")


# Get genes with a significant ASE measurement in at least third of samples
AseGenes33perc <- as.character(AseGeneCounts[AseGeneCounts$Freq>=(ncol(absData)/3),1])

# FUNCTION FOR TOPGO
ase.topGO <- function(inputData, statvalue = "pvalue", filename, genelist = FALSE) {
  if (genelist == TRUE) {
    # Run the topGO analysis
    for (j in c("BP", "CC", "MF")) {			
      message(j)		
      xx <- annFUN.org(j, mapping = "org.Hs.eg.db", ID = "ensembl")
      allGenes <- unique(unlist(xx))		
      allGenesFactor <- factor(as.integer(allGenes %in% inputData))		
      names(allGenesFactor) <- allGenes		
      GOdata <- new("topGOdata",		
                    ontology = j,
                    allGenes = allGenesFactor,
                    nodeSize = 5, # recommended setting
                    annot = annFUN.org,
                    mapping = "org.Hs.eg.db",
                    ID = "ensembl")
      goRes <- runTest(GOdata, algorithm = "parentchild", statistic = "fisher") # recommended setting		
      allRes <- GenTable(GOdata, Pvalue = goRes, 		
                         topNodes = length(goRes@score),
                         orderBy = "Pvalue", ranksOf = "Pvalue",
                         numChar = 1E9)
      if (all(as.numeric(allRes$Pvalue) >= 0.05)) {
        allRes <- allRes[1, ]
        allRes[, 1:6] <- NA
      } else {
        allRes <- allRes[1:max(which(as.numeric(allRes[,6]) < 0.05)), ] # leave only results with p < 0.05
      }
      write.xlsx(allRes,		
                 file = filename, 
                 sheet = ifelse(j == "BP", "BiologicalProcess", 
                                ifelse(j == "CC", "CellularComponent", "MolecularFunction")),
                 row.names = FALSE,
                 append = TRUE)
    }
  } else {
    # Create a data frame of all IDs with corresponding lowest qvalues for each gene
    genesQsample <- inputData[,c(statvalue,"ensembl_gene_id")]
    
    genesQsample <- genesQsample[order(genesQsample$ensembl_gene_id, 
                                       genesQsample[[statvalue]]),]
    
    # Take most significant hit per person per gene
    genesQsample <- genesQsample %>% distinct(ensembl_gene_id, .keep_all= TRUE)
    genesQsample <- na.omit(genesQsample)
    
    # Get list of all the significant genes
    AseGenes <- na.omit(unique(genesQsample$ensembl_gene_id[genesQsample[[statvalue]]<0.05]))
    
    # Run the topGO analysis
    for (j in c("BP", "CC", "MF")) {			
      message(j)		
      xx <- annFUN.org(j, mapping = "org.Hs.eg.db", ID = "ensembl")
      allGenes <- unique(unlist(xx))		
      allGenesFactor <- factor(as.integer(allGenes %in% AseGenes))		
      names(allGenesFactor) <- allGenes		
      GOdata <- new("topGOdata",		
                    ontology = j,
                    allGenes = allGenesFactor,
                    nodeSize = 5, # recommended setting
                    annot = annFUN.org,
                    mapping = "org.Hs.eg.db",
                    ID = "ensembl")
      goRes <- runTest(GOdata, algorithm = "parentchild", statistic = "fisher") # recommended setting		
      allRes <- GenTable(GOdata, Pvalue = goRes, 		
                         topNodes = length(goRes@score),
                         orderBy = "Pvalue", ranksOf = "Pvalue",
                         numChar = 1E9)
      if (all(as.numeric(allRes$Pvalue) >= 0.05)) {
        allRes <- allRes[1, ]
        allRes[, 1:6] <- NA
      } else {
        allRes <- allRes[1:max(which(as.numeric(allRes[,6]) < 0.05)), ] # leave only results with p < 0.05
      }
      write.xlsx(allRes,		
                 file = filename, 
                 sheet = ifelse(j == "BP", "BiologicalProcess", 
                                ifelse(j == "CC", "CellularComponent", "MolecularFunction")),
                 row.names = FALSE,
                 append = TRUE)
    }
  }
}

ase.topGO(inputData = snps2genesQ,
          statvalue = "qvalue",
          filename = "topGO_total.xlsx")

# topGO for genes shared in at least 33% of samples
ase.topGO(inputData = AseGenes33perc,
          genelist = TRUE,
          filename = "topGO_33perc.xlsx")

#-----------------------------------------------------------------------------#
# 7. eQTL and sQTL enrichment analysis (and topGO per cluster)
#-----------------------------------------------------------------------------#

setwd(DATA.DIR)

gtexeqtl <- 
  read.table("GTEx_Analysis_v8_eQTL/Heart_Left_Ventricle.v8.egenes.txt",
             header = T, sep = "\t")
gtexeqtl$gene_id <- gsub("\\..*","",gtexeqtl$gene_id)

gtexsqtl <- 
  read.table("GTEx_Analysis_v8_sQTL/Heart_Left_Ventricle.v8.sgenes.txt",
             header = T, sep = "\t")
gtexsqtl$gene_id <- gsub("\\..*","",gtexsqtl$gene_id)

setwd(RES.DIR)

# Function for sQTL and eQTL enrichment analysis
# Input is a dataframe of ensembl_gene_id and pvalue
# Needs to know name of gtex sQTL and eQTL data frame
# keep qval column name for eQTL and sQTL data
ase.qtl.enrich <- function(inputData, statvalue = "pvalue", sqtl, eqtl, cutoff = 0.05) {
  genesQsample <- inputData[,c(statvalue,"ensembl_gene_id")]
  
  genesQsample <- genesQsample[order(genesQsample$ensembl_gene_id, 
                                     genesQsample[[statvalue]]),]
  
  # Take most significant hit per person per gene
  genesQsample <- genesQsample %>% distinct(ensembl_gene_id, .keep_all= TRUE)
  genesQsample <- na.omit(genesQsample)
  
  # Get list of all the significant genes
  AseGenes <- genesQsample %>% filter(genesQsample[[statvalue]]<cutoff)
  NonAseGenes <- genesQsample %>% filter(genesQsample[[statvalue]]>=cutoff)
  
  # Get eGene and non-eGene counts for ASE and non-ASE genes
  AseEgenes <- nrow(unique(AseGenes[AseGenes$ensembl_gene_id %in% 
                                      eqtl[eqtl$qval<0.05,"gene_id"],]))
  
  AseNonEgenes <- nrow(unique(AseGenes[AseGenes$ensembl_gene_id %in% 
                                         eqtl[eqtl$qval>=0.05,"gene_id"],]))
  
  NonAseEgenes <- nrow(unique(NonAseGenes[NonAseGenes$ensembl_gene_id %in% 
                                            eqtl[eqtl$qval<0.05,"gene_id"],]))
  
  NonAseNonEgenes <- nrow(unique(NonAseGenes[NonAseGenes$ensembl_gene_id %in% 
                                               eqtl[eqtl$qval>=0.05,"gene_id"],]))
  
  # Get sGene and non-sGene counts for ASE and non-ASE genes
  AseSgenes <- nrow(unique(AseGenes[AseGenes$ensembl_gene_id %in% 
                                      sqtl[sqtl$qval<0.05,"gene_id"],]))
  
  AseNonSgenes <- nrow(unique(AseGenes[AseGenes$ensembl_gene_id %in% 
                                         sqtl[sqtl$qval>=0.05,"gene_id"],]))
  
  NonAseSgenes <- nrow(unique(NonAseGenes[NonAseGenes$ensembl_gene_id %in% 
                                            sqtl[sqtl$qval<0.05,"gene_id"],]))
  
  NonAseNonSgenes <- nrow(unique(NonAseGenes[NonAseGenes$ensembl_gene_id %in% 
                                               sqtl[sqtl$qval>=0.05,"gene_id"],]))
  
  # Create 2x2 matrix and perform Fisher test
  contE <- matrix(c(AseEgenes, AseNonEgenes, NonAseEgenes, NonAseNonEgenes), 
                  nrow = 2, byrow = T)
  contS <- matrix(c(AseSgenes, AseNonSgenes, NonAseSgenes, NonAseNonSgenes), 
                  nrow = 2, byrow = T)
  cat(paste0("The p-value for the eQTL Fisher's exact test is ",
             fisher.test(contE)$p.value, "\n",
             "The p-value for the sQTL Fisher's exact test is ", fisher.test(contS)$p.value))
}

# QTL enrichment for the imbalances over all samples
ase.qtl.enrich(inputData = snps2genesQ,
               statvalue = "qvalue",
               sqtl = gtexsqtl,
               eqtl = gtexeqtl)

#-----------------------------------------------------------------------------#
# 8. Save stuff
#-----------------------------------------------------------------------------#

setwd(DATA.DIR)
save.image("ASE_RESULTS.RData")
# END
