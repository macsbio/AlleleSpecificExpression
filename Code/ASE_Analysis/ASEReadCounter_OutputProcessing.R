#=============================================================================#
# ASEsampleDataReadCounter_OutputProcessing.R                                 #
#                                                                             #
# Version: 1.3                                                                #
# Date: November 10, 2020                                                     #
# Author: Michiel Adriaens, PhD; MaCSBio, Maastricht University               #
# History:                                                                    #
#  1.0: Creation                                                              #
#  1.1: Updated a few lines based on Kasper's suggestions                     #
#  1.2: Updated for semi-complete dataset (n = 80)                            # 
#                                                                             #
#=============================================================================#

# Load libraries and set some stuff
#----------------------------------
options(stringsAsFactors = FALSE)
require(dplyr)
rm(list = ls(all.names = T))

# Set data directory, gather filenames
#-------------------------------------
DATA.DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
PROJECT.NAME <- "DCM"
setwd(DATA.DIR)
# Change sample ID prefix in the next line if necessary
fileNames <- dir()[grep("SID.+\\.csv$", dir())]
sampleNames <- make.names(gsub("\\_ASEresults\\.csv", "", fileNames))
sampleNames <- t(as.data.frame(sampleNames))
colnames(sampleNames) <- sampleNames[1,]

sampleNames <- as.character(sampleNames[1,])

# Import data and calculate genotype data, alt an ref ratios, allele dosage 
#  (w.r.t. alternative allele), etc.
#--------------------------------------------------------------------------
arcData <- list()
for (i in seq(along = sampleNames)) {
  # Select sampleData from complete dataset
  sampleData <- read.delim(fileNames[i], sep = ",", as.is = T)
  sampleData$variantID[sampleData$variantID=="."] <- 
    paste0(sampleData$contig[sampleData$variantID=="."],
           ":",
           sampleData$position[sampleData$variantID=="."])
  # Remove all data that does not have a dbSNP ID associated:
  # FIXME: or should I add dummy names a la 1KG? e.g. 'chr1:10000'
  # Create genotype matrix for allele1 and allele2
  al <- cbind(rep(NA, nrow(sampleData)), rep(NA, nrow(sampleData)))
  # Create indexes for homozygous alt, homozygous ref and heterozygous
  #  variants
  altHomIndex <- sampleData[, "refCount"] == 0
  refHomIndex <- sampleData[, "altCount"] == 0
  hetIndex <- sampleData[, "altCount"] != 0 & sampleData[, "refCount"] != 0
  # Fill out the genotype matrix using these indexes
  al[altHomIndex, ] <- as.character(sampleData[altHomIndex, "altAllele"])
  al[refHomIndex, ] <- as.character(sampleData[refHomIndex, "refAllele"])
  al[hetIndex, ] <- cbind(as.character(sampleData[hetIndex, "refAllele"]),
                          as.character(sampleData[hetIndex, "altAllele"]))
  # Create alt allele dosage vector
  dosage <- rep(0, nrow(sampleData))
  dosage[altHomIndex] <- 2
  dosage[hetIndex] <- 1 
  # Append new data to sampleData and replace in complete dataset
  sampleData$al1 <- al[, 1]
  sampleData$al2 <- al[, 2]
  sampleData$genotype <- paste0(al[, 1], al[, 2])
  sampleData$dosage <- dosage
  sampleData$refRatio <- sampleData$refCount / sampleData$totalCount
  sampleData$altRatio <- sampleData$altCount / sampleData$totalCount
  arcData[[i]] <- sampleData
  cat("\r", i)
}; cat("\n")

# Create altRatio, refRatio, totalCount and dosage matrices
#----------------------------------------------------------
date()
for (k in 2:length(sampleNames)) {
  if (k == 2) {
    dosageData <- merge(arcData[[1]][, c("variantID", "dosage")], 
                        arcData[[2]][, c("variantID", "dosage")], by = "variantID", 
                        all = TRUE)
    altRatioData <- merge(arcData[[1]][, c("variantID", "altRatio")], 
                          arcData[[2]][, c("variantID", "altRatio")], by = "variantID", 
                          all = TRUE)
    refRatioData <- merge(arcData[[1]][, c("variantID", "refRatio")], 
                          arcData[[2]][, c("variantID", "refRatio")], by = "variantID", 
                          all = TRUE)
    gtData <- merge(arcData[[1]][, c("variantID", "genotype")], 
                    arcData[[2]][, c("variantID", "genotype")], by = "variantID", 
                    all = TRUE)
    totalCountData <- merge(arcData[[1]][, c("variantID", "totalCount")], 
                            arcData[[2]][, c("variantID", "totalCount")], by = "variantID",
                            all = TRUE)
    refCountData <- merge(arcData[[1]][, c("variantID", "refCount")], 
                          arcData[[2]][, c("variantID", "refCount")], by = "variantID",
                          all = TRUE)
    altCountData <- merge(arcData[[1]][, c("variantID", "altCount")], 
                          arcData[[2]][, c("variantID", "altCount")], by = "variantID",
                          all = TRUE)
    annData.temp <- arcData[[1]][, 1:5]
    annData.temp <- distinct(rbind(annData.temp, arcData[[k]][, 1:5])) 
  } else {
    dosageData <- merge(dosageData, 
                        arcData[[k]][, c("variantID", "dosage")], by = "variantID", 
                        all = TRUE)
    altRatioData <- merge(altRatioData, 
                          arcData[[k]][, c("variantID", "altRatio")], by = "variantID", 
                          all = TRUE)
    refRatioData <- merge(refRatioData, 
                          arcData[[k]][, c("variantID", "refRatio")], by = "variantID", 
                          all = TRUE)
    gtData <- merge(gtData, 
                    arcData[[k]][, c("variantID", "genotype")], by = "variantID", 
                    all = TRUE)
    totalCountData <- merge(totalCountData, 
                            arcData[[k]][, c("variantID", "totalCount")], by = "variantID",
                            all = TRUE)
    refCountData <- merge(refCountData, 
                          arcData[[k]][, c("variantID", "refCount")], by = "variantID",
                          all = TRUE)
    altCountData <- merge(altCountData, 
                          arcData[[k]][, c("variantID", "altCount")], by = "variantID",
                          all = TRUE)
    annData.temp <- distinct(rbind(annData.temp, arcData[[k]][, 1:5]))
  }
  cat("\r", k)
}; cat("\n")
date()
colnames(dosageData)[-1] <- colnames(altRatioData)[-1] <- 
  colnames(refRatioData)[-1] <- colnames(gtData)[-1] <- 
  colnames(totalCountData)[-1] <- colnames(altCountData)[-1] <- 
  colnames(refCountData)[-1] <- sampleNames

# Align annData with other objects
#----------------------------------
annData <- annData.temp[match(altRatioData[, 1], annData.temp[, 3]), ]
colnames(annData) <- c("chromosome", "position", "id", "refAllele",
                       "altAllele")

# Some checks:
all(dosageData[, 1]       ==     altRatioData[, 1])
all(altRatioData[, 1]     ==     refRatioData[, 1])
all(refRatioData[, 1]     ==     gtData[, 1])
all(gtData[, 1]           ==     totalCountData[, 1])
all(totalCountData[, 1]   ==     annData[, 3])

# Add MAF, minor allele
#----------------------
# NB: dosage is w.r.t. alternative allele
n0 <- apply(dosageData == 0, 1, sum, na.rm = T)
n1 <- apply(dosageData == 1, 1, sum, na.rm = T)
n2 <- apply(dosageData == 2, 1, sum, na.rm = T)
n <- n0 + n1 + n2
p <- (2 * n0 + n1) / (2 * n) # ref allele freq.
q <- 1 - p # alt. allele freq
annData$maf <- pmin(p, q)
annData$minorAllele <- NA
annData$minorAllele[p < q] <- annData$refAllele[p < q]
annData$minorAllele[q < p] <- annData$altAllele[q < p]

# Save data to disk
#------------------
setwd(DATA.DIR)
save(list = c('dosageData', 'altRatioData',
              'totalCountData', 'annData','refRatioData', 'refCountData', 'altCountData'), 
     file = paste('ASEReadCounterOutput_Processed_', PROJECT.NAME, 
                  '.RData', sep = ""))
# END