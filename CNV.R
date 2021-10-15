###################################################COPY NUMBER VARIANT ANALYSIS#############################################
#R CNV.R bedfile Normal.bam Tumor.Bam outputt
library(CNVPanelizer)
library(data.table)

#reading the arguments

args = commandArgs(TRUE)
bedfile <-args[1]
Bamfile_Normal <-args[2]
Bamfile_Tumor <-args[3]
output <-args[4]

genomicRangesFromBed <-  BedToGenomicRanges(paste0(bedfile),ampliconColumn = 4,split = "_")
metadataFromGenomicRanges <-  elementMetadata(genomicRangesFromBed)
geneNames = metadataFromGenomicRanges["geneNames"][, 1]
ampliconNames = metadataFromGenomicRanges["ampliconNames"][, 1]
referenceReadCounts <-  ReadCountsFromBam(paste0(Bamfile_Normal), gr = genomicRangesFromBed,ampliconNames = ampliconNames, sampleNames = "Normal",removeDup = FALSE)
sampleReadCounts <-  ReadCountsFromBam(paste0(Bamfile_Tumor),gr = genomicRangesFromBed,ampliconNames = ampliconNames,sampleNames = "Tumor",removeDup = FALSE)
normalizedReadCounts <-  CombinedNormalizedCounts(sampleReadCounts,referenceReadCounts,ampliconNames = ampliconNames)
samplesNormalizedReadCounts = normalizedReadCounts["samples"][[1]]
referenceNormalizedReadCounts = normalizedReadCounts["reference"][[1]]
bootList <- BootList(geneNames,samplesNormalizedReadCounts,referenceNormalizedReadCounts,replicates = 100)
backgroundNoise <- Background(geneNames,samplesNormalizedReadCounts,referenceNormalizedReadCounts,bootList, replicates = 100,significanceLevel = 0.1)
reportTables <-  ReportTables(geneNames,samplesNormalizedReadCounts,referenceNormalizedReadCounts,bootList, backgroundNoise)
cnvs =as.data.frame(reportTables)
cnvs$Tumor.Gene = rownames(cnvs)
write.csv(cnvs, file = paste0(output,"_CNV.csv"),row.names=FALSE)
