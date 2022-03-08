# minfi tutorial
# 
# Analysis of 450k data using minfi
#
# https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html

library(minfi)
library(minfiData)
library(sva)

baseDir <- system.file("extdata", package="minfiData")
targets <- read.metharray.sheet(baseDir)
RGSet <- read.metharray.exp(targets = targets)

phenoData <- pData(RGSet)
phenoData[,1:6]

manifest <- getManifest(RGSet)
manifest

head(getProbeInfo(manifest))

MSet <- preprocessRaw(RGSet) 
MSet

head(getMeth(MSet)[,1:3])

head(getUnmeth(MSet)[,1:3])

RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
RSet

beta <- getBeta(RSet)

GRset <- mapToGenome(RSet)
GRset

beta <- getBeta(GRset)
M <- getM(GRset)
CN <- getCN(GRset)

sampleNames <- sampleNames(GRset)
probeNames <- featureNames(GRset)
pheno <- pData(GRset)

gr <- granges(GRset)
head(gr, n= 3)

annotation <- getAnnotation(GRset)
names(annotation)

# Quality control ####

head(getMeth(MSet)[,1:3])

head(getUnmeth(MSet)[,1:3])

qc <- getQC(MSet)
head(qc)

plotQC(qc)

densityPlot(MSet, sampGroups = phenoData$Sample_Group)

densityBeanPlot(MSet, sampGroups = phenoData$Sample_Group)

controlStripPlot(RGSet, controls="BISULFITE CONVERSION II")

# qcReport(RGSet, pdf= "qcReport.pdf")

## Infer sex ####

predictedSex <- getSex(GRset, cutoff = -2)$predictedSex
head(predictedSex)

# error here :(
# error in evaluating the argument 'table' in selecting a method for 
# function '%in%': error in evaluating the argument 'x' in selecting 
# a method for function 'colnames': unable to find an inherited method 
# for function ‘colData’ for signature ‘"DFrame"’
plotSex(getSex(GRset, cutoff = -2))

# Preprocessing and normalization ####

MSet.illumina <- preprocessIllumina(RGSet, bg.correct = TRUE,
                                    normalize = "controls")

MSet.swan <- preprocessSWAN(RGSet)

GRset.quantile <- preprocessQuantile(RGSet, fixOutliers = TRUE,
                                     removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                     quantileNormalize = TRUE, stratified = TRUE, 
                                     mergeManifest = FALSE, sex = NULL)

# Another error :(

GRset.funnorm <- preprocessFunnorm(RGSet)

# Genetic variants and cell type composition ####
snps <- getSnpInfo(GRset)
head(snps,10)

GRset <- addSnpInfo(GRset)

GRset <- dropLociWithSnps(GRset, snps=c("SBE","CpG"), maf=0)

# library(FlowSorted.Blood.450k)
# cellCounts <- estimateCellCounts(RGset)

# Identifying DMRs and DMPs ####
beta <- getBeta(GRset.funnorm)
age  <- pData(GRset.funnorm)$age
dmp <- dmpFinder(beta, pheno = age  , type = "continuous")
head(dmp)
