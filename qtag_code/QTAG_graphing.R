#! /usr/bin/Rscript
# Gradient Binning Output- graphing
# setwd("~")
# setwd("/Volumes/MCDataDrive/Users/melissachen/Documents/Masters_data/Undergrad_Environmental/Roughwork/ALLDATASETS_7march2017_16s/2_analysis_ALL/SALBIN_16s_FRASER_13march2017/")
########################## OPTPARSE ####################################
library(optparse)
library(MASS)
library(tidyverse)

option_list = list(
  make_option(c("-b", "--boundaries"), type="character",
              help="Boundaries.txt output from binning script", default = "boundaries.txt"),
  make_option(c("-M", "--model_boundaries"),
              help="Model Boundaries output from binning script", type="character", default = 'modelBoundaries_type.txt'),
  make_option(c("-t", "--taxa_abund"), 
              help="Taxa abundances across gradient file output from binning script", type="character", default = 'taxa_abundances_across_gradient.txt'),
  make_option(c("-T", "--types"),
              help="Types across gradient output file from binning script", type="character", default = "types_across_gradient"),
  make_option(c("-N","--Namesgradient"), 
              help= "Gradient text output from binning script", type = "character", default = "gradient.txt"),
  make_option(c("-B", "--biom"), type="character",
              help="biom-- already collapsed, but in biom format.", metavar="character", default = "OTUTableText.txt"),
  make_option(c("-n", "--minthreshold"), type="character",
              help="Minimum relabundance to be plotted", metavar="numeric", default = 0.05),
  make_option(c("-O", "--OTUTable"), type = "character",
              help="Type OTU Table output from Classification script", metavar = "character", default = "OTUTablebyType.txt"),
  make_option(c("-R", "--toleranceRanges"), type = "character",
              help="Tolerance Ranges across gradient", metavar = "character", default = "toleranceRanges.txt"),
  make_option(c("-L","--LOG"), 
              help= "LOG output from binning script", type = "character", default = "LOG.txt")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

boundariesPWD = opt$boundaries
modelboundariesPWD = opt$model_boundaries
taxaAbundPWD = opt$taxa_abund
typesPWD = opt$types
gradientNamesPWD = opt$Namesgradient
biomPWD = opt$biom
minthreshold = opt$minthreshold
LOGPWD = opt$LOG
toleranceRangesPWD = opt$toleranceRanges
typetablePWD = opt$OTUTable


########################## For Testing ####################################
# setwd('/Users/melissachen/Documents/Masters/Project_QTAG/laptop_current_24oct2019/SALBIN_20april2020/18sFraser/')
# 
# setwd('/Users/melissachen/Documents/Masters/Project_QTAG/testing/testing_b16')
# wd <- getwd()
# boundariesPWD = unlist(paste0(wd, "/boundaries.txt"))
# modelboundariesPWD = paste0(wd, "/modelBoundaries_type.txt")
# taxaAbundPWD = paste0(wd, "/taxa_abundances_across_gradient.txt")
# typesPWD = paste0(wd, "/types_across_gradient")
# gradientNamesPWD = paste0(wd, "/gradient.txt")
# biomPWD = paste0(wd, "/OTUTableText.txt")
# minthreshold = 0.05
# LOGPWD = paste0(wd, "/LOG.txt")
# toleranceRangesPWD=paste0(wd,"/toleranceRanges.txt")
# typetablePWD = paste0(wd, "/OTUTablebyType.txt")

# setwd('/Users/melissachen/Documents/Masters/Project_Environmental/FromBotaClust_1feb2018/SALBIN_7June2018_vsearch_silva132/16sFRASER_7june2018')
# wd <- getwd()
# boundariesPWD = unlist(paste0(wd, "/boundaries.txt"))
# modelboundariesPWD = paste0(wd, "/modelBoundaries_type.txt")
# taxaAbundPWD = paste0(wd, "/taxa_abundances_across_gradient.txt")
# typesPWD = paste0(wd, "/types_across_gradient")
# gradientNamesPWD = paste0(wd, "/gradient.txt")
# biomPWD = paste0(wd, "/OTUTableText.txt")
# minthreshold = 0.05
# LOGPWD = paste0(wd, "/LOG.txt")
# typetablePWD = paste0(wd, "/OTUTablebyType.txt")

# 
# setwd('/Users/melissachen/Documents/Masters/Project_Environmental/FromBotaClust_1feb2018/SALBIN_7June2018_vsearch_silva132/18sFraser_7june2018')
# wd <- getwd()
# boundariesPWD = unlist(paste0(wd, "/boundaries.txt"))
# modelboundariesPWD = paste0(wd, "/modelBoundaries_type.txt")
# taxaAbundPWD = paste0(wd, "/taxa_abundances_across_gradient.txt")
# typesPWD = paste0(wd, "/types_across_gradient")
# gradientNamesPWD = paste0(wd, "/gradient.txt")
# biomPWD = paste0(wd, "/OTUTableText.txt")
# minthreshold = 0.05
# LOGPWD = paste0(wd, "/LOG.txt")
# typetablePWD = paste0(wd, "/OTUTablebyType.txt")


########################## LOAD FILES ####################################

# Load files
print("Loading files for graphing...")

boundaries <- read.delim(paste0(boundariesPWD), header = FALSE, strip.white = TRUE)
# histBoundaries <- as.numeric(boundaries$V1)

logvalues <- read.delim(paste0(LOGPWD), header = FALSE, strip.white = TRUE)
minGrad <- as.character(logvalues[grep("minX", logvalues$V1),])
minGrad <- floor(as.numeric(gsub("minX =","", minGrad)))
maxGrad <- as.character(logvalues[grep("maxY", logvalues$V1),])
maxGrad <- ceiling(as.numeric(gsub("maxY =","", maxGrad)))
metadataPWD <- as.character(logvalues[grep("metadata", logvalues$V1),])
metadataPWD <- gsub("metadata = ","../",metadataPWD)

threshold <- as.character(logvalues[grep("^threshold:", logvalues$V1),])
threshold <- as.numeric(gsub("threshold:", "", threshold))

modelBoundaries <- read.delim(paste0(modelboundariesPWD), header = TRUE, row.names = 1, strip.white = TRUE, stringsAsFactors = FALSE, na.strings = c('','na','NA'))

taxaAbundances <- read.delim(paste0(taxaAbundPWD), header = FALSE, row.names = 1, strip.white = TRUE, stringsAsFactors = FALSE, na.strings = c('','na','NA'))
gradient <- as.numeric(taxaAbundances[1,])
taxaAbundances <- taxaAbundances[-1,]

typesAll <- read.delim(paste0(typesPWD,"_all.txt"), header = TRUE, strip.white = TRUE, stringsAsFactors = FALSE, row.names = 1, na.strings = c('','na','NA'))
colnames(typesAll) <- gsub("X", "", colnames(typesAll))
typesCondensed <- read.delim(paste0(typesPWD,"_condensed.txt"), header = TRUE, strip.white = TRUE, stringsAsFactors = FALSE, row.names = 1, na.strings = c('','na','NA'))
colnames(typesCondensed) <- gsub("X", "", colnames(typesCondensed))

gradientNames <- read.delim(paste0(gradientNamesPWD), header = FALSE, strip.white = TRUE, sep = ",", stringsAsFactors = FALSE)

toleranceRanges <- read.delim(paste0(toleranceRangesPWD), header = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
toleranceRanges$OTU <- as.character(toleranceRanges$OTU)
########################## HISTOGRAM (Single) ####################################

# Plot histogram
print("Plotting single histogram")
ggsave("BoundaryDistributions.pdf"
       ,ggplot(data=boundaries) + geom_histogram(aes(x=V1), bins=25) + geom_vline(aes(xintercept=mean(V1)), col="red") +
         xlab("Salinity") + ylab("Frequency") +theme_bw())

########################## TAXA ABUND (Individual) ############

# Taxa abundances for each taxa
print("Plotting individual taxa abund")
system("mkdir INDIV_PLOTS")
taxaTotal=nrow(taxaAbundances)
pb <- txtProgressBar(min=0, max=taxaTotal, style=3)
for (x in (1:nrow(taxaAbundances))) {
  Sys.sleep(0.1)
  setTxtProgressBar(pb,x)
  taxa <- rownames(taxaAbundances)[x]
  Avalue <- modelBoundaries[taxa,'A']
  Bvalue <- modelBoundaries[taxa,'B']
  testB <- modelBoundaries[taxa,'sigAB']
  meanA <- modelBoundaries[taxa,'meanA']
  meanB <- modelBoundaries[taxa,'meanB']
  meanC <- modelBoundaries[taxa,'meanC']
  typeMB <- modelBoundaries[taxa,'type']
  bloom <- modelBoundaries[taxa, 'bloom']
  bounds <- modelBoundaries[taxa,c('boundaries','boundariestwo')]
  pdf(file = paste0("INDIV_PLOTS/",typeMB,'_Indiv',taxa,'.pdf'))
  plot(as.numeric(taxaAbundances[x,]) ~gradient
       , xlab = 'Salinity'
       , ylab = 'Relative Abundance'
       , main = paste(taxa)
       , sub = paste(typeMB,"bloom?",bloom,signif(modelBoundaries[taxa,'sigAB'],3)," ",signif(modelBoundaries[taxa,'sigBC'],3), " ",signif(modelBoundaries[taxa,'sigAC'],1)))
  if (testB != 'None') {
  lines(c(minGrad,Avalue), c(meanA,meanA)
        , lwd = 3
        , col = 'blue')
  lines(c(Avalue,Bvalue), c(meanB,meanB)
        , lwd = 3
        , col = 'purple')
  lines(c(Bvalue,maxGrad), c(meanC,meanC)
        , lwd = 3
        , col = 'red')
  } else if (testB == 'None') {
    lines(c(minGrad,Avalue), c(meanA,meanA)
          , lwd = 3
          , col = 'blue')
    lines(c(Avalue,Bvalue), c(meanA,meanC)
          , lwd = 3
          , col = 'purple')
    lines(c(Bvalue,maxGrad), c(meanC,meanC)
          , lwd = 3
          , col = 'red')
        }
  for (b in as.numeric(bounds)) {
    abline(v=b, col="grey")
  }
  dev.off()
}
close(pb)

########################## BARPLOT: Taxa types across gradient ####################################

# Taxa types across gradient
print("Plotting barplots for taxa types")

colorsGradient <- c("blue" #lorestr
                    ,"darkblue" # lobloom
                    ,"cyan" #lopeak
                    ,"slateblue1" #bracklotoler
                    ,"purple" #brackrestr
                    ,"magenta" # brack bloom
                    ,"pink4" #brack peak
                    ,"violetred" #brackpeakhitoler
                    ,"orange" #hipeak
                    ,"red" #hirestr
                    ,"darkred"#hibloom
                    ,"black" #noclass
                    , "grey" #ubiq
                    ,"white" #invbrack
                    )
colorsGradientCondensed <- c("blue","purple","red","black")


typesAll <- as.matrix(typesAll)
pdf(file = "taxatypesacrossgradient_all.pdf")
par(mar=c(5.1,4.1,4.1,9), xpd = TRUE)
barplot(typesAll, beside = FALSE
      , col = colorsGradient
      , main = 'Taxa types across gradient'
      , xlab = 'Gradient'
      , ylab = 'Relative Abundance')
legend('right', inset = c(-0.4,0), legend = rev(rownames(typesAll)), pch = 21, col = c("black",rev(colorsGradient)[-1]), pt.bg = rev(colorsGradient), bty="n")
dev.off()


typesCondensed <- as.matrix(typesCondensed)
pdf(file = "taxatypesacrossgradient_condensed.pdf")
par(mar=c(5.1,4.1,4.1,9), xpd = TRUE)
barplot(typesCondensed, beside = FALSE
      , col = colorsGradientCondensed
      , main = 'Taxa types across gradient'
      , xlab = 'Gradient'
      , ylab = 'Relative Abundance')
legend('right', inset = c(-0.4,0), legend = rev(rownames(typesCondensed)), pch = 21, pt.bg = rev(colorsGradientCondensed), col = rev(colorsGradientCondensed), bty="n")
dev.off()

  
  ########################## TAXA ABUND (by type) ####################################
system('mkdir TaxaAbund_bytype_all')
system('mkdir TaxaAbund_bytype_condensed')
print("Making taxa abundance plots by type")


for (i in 1:nrow(typesAll)){
pdf(file = paste0('./TaxaAbund_bytype_all/Indiv_', rownames(typesAll)[i], '_all.pdf'))
plot(typesAll[i,] ~ colnames(typesAll)
     , type = 'h'
     , lwd = 10
     , col = 'darkgrey'
     , main = paste0(rownames(typesAll)[i],' abundance across gradient')
     , xlab = 'Gradient'
     , ylab = 'Relative abundance')
dev.off()
}



for (i in 1:nrow(typesCondensed)){
pdf(file = paste0('./TaxaAbund_bytype_condensed/Indiv_', rownames(typesCondensed)[i], '.pdf'))
plot(typesCondensed[i,] ~ colnames(typesCondensed)
     , type = 'h'
     , lwd = 10
     , col = 'darkgrey'
     , main = paste0(rownames(typesCondensed)[i],' abundance across gradient')
     , xlab = 'Gradient'
     , ylab = 'Relative abundance')
dev.off()
}


########################## HISTOGRAM (by type) ####################################
# Histogram but separated by 'type'
print ("Making histograms, by type")

lowOnly <- modelBoundaries[grep(paste0(gradientNames[1]), modelBoundaries$typeSimple),]
lowOnlyTrue <- modelBoundaries[grep(paste0("^",gradientNames[1],"Restricted$"), modelBoundaries$type),]
lowOnlyHalf <- modelBoundaries[grep(paste0(gradientNames[1],"Peak"), modelBoundaries$type),]
lowOnlyBloom <- modelBoundaries[grep(paste0(gradientNames[1],"Bloom"), modelBoundaries$type),]


highOnly <- modelBoundaries[grep(paste0(gradientNames[3]), modelBoundaries$typeSimple),]
highOnlyTrue <- modelBoundaries[grep(paste0("^",gradientNames[3],"Restricted$"), modelBoundaries$type),]
highOnlyHalf <- modelBoundaries[grep(paste0(gradientNames[3], "Peak"), modelBoundaries$type),]
highOnlyBloom <- modelBoundaries[grep(paste0(gradientNames[3], "Bloom"), modelBoundaries$type),]

interOnly <- modelBoundaries[grep(paste0(gradientNames[2]), modelBoundaries$typeSimple),]
interOnlyTrue <- modelBoundaries[grep(paste0("^",gradientNames[2],"Restricted$"), modelBoundaries$type),]
interOnlyHalf <- modelBoundaries[grep(paste0(gradientNames[2], "Peak"), modelBoundaries$type),]
interOnlyBloom <- modelBoundaries[grep(paste0(gradientNames[2], "Bloom"), modelBoundaries$type),]
interOnlyHi <- modelBoundaries[grep(paste0(gradientNames[2], "PeakHiToler"), modelBoundaries$type),]
interOnlylo <- modelBoundaries[grep(paste0(gradientNames[2], "PeakLoToler"), modelBoundaries$type),]

# Set even bin widths 
is.even <- function(x) { x%%2 == 0}
if (is.even(maxGrad-minGrad)){
  binBreaks <- seq(minGrad,maxGrad, by = 2)
} else {
  binBreaks <- seq(minGrad,maxGrad+1, by = 2)
}

system("mkdir Histograms")

########################## Low ####################################

allLow <- c(lowOnly$boundariestwo)
TrueLow <- c(lowOnlyTrue$boundariestwo)
HalfLow <- c(lowOnlyHalf$boundariestwo)
BloomLow <- c(lowOnlyBloom$boundariestwo)
nonBloomLow <- c(HalfLow, TrueLow)



lowBoundaries <- c()
if (length(allLow) > 0){
for (i in 1:length(allLow)) {
  if (!is.na(allLow[i])) {
    lowBoundaries <- c(lowBoundaries, allLow[i])
  } 
}
meanLow <- round(mean(lowBoundaries),2)
pdf(file = paste0("./Histograms/Histogram_",gradientNames[1],"Boundaries_all.pdf"))
hist(as.numeric(lowBoundaries)
     , xlim = range(minGrad,maxGrad)
     , xlab = "Gradient"
     , breaks = binBreaks
     , col = "grey"
     , main = paste0("Boundary distribution for ALL ", gradientNames[1],": u=",as.character(meanLow)))
abline(v=meanLow
       , col = "red"
       , lwd = 2)
dev.off()
}

lowBoundaries <- c()
if (length(nonBloomLow) > 0){
for (i in 1:length(nonBloomLow)) {
  if (!is.na(nonBloomLow[i])) {
    lowBoundaries <- c(lowBoundaries, nonBloomLow[i])
  } 
}
meanLow <- round(mean(lowBoundaries),2)
pdf(file = paste0("./Histograms/Histogram_",gradientNames[1],"Boundaries_nonBloom.pdf"))
hist(lowBoundaries
     , xlim = range(minGrad,maxGrad)
     , xlab = "Gradient"
     , col = "grey"
     , breaks = binBreaks
     , main = paste0("Boundary distribution for nonBLOOM ", gradientNames[1],": u=",as.character(meanLow)))
abline(v=meanLow
       , col = "red"
       , lwd = 2)
dev.off()
}


# True
lowBoundaries <- c()
if (length(TrueLow) > 0){
for (i in 1:length(TrueLow)) {
  if (!is.na(TrueLow[i])) {
    lowBoundaries <- c(lowBoundaries, TrueLow[i])
  } 
}
meanLow <- round(mean(lowBoundaries),2)
pdf(file = paste0("./Histograms/Histogram_",gradientNames[1],"Boundaries_True.pdf"))
hist(lowBoundaries
     , xlim = range(minGrad,maxGrad)
     , xlab = "Gradient"
     , breaks = binBreaks
     , col = "grey"
     , main = paste0("Boundary distribution for TRUE ", gradientNames[1],": u=",as.character(meanLow)))
abline(v=meanLow
       , col = "red"
       , lwd = 2)
dev.off()
}

# Half
if (length(HalfLow) > 0){
for (i in 1:length(HalfLow)) {
  if (!is.na(HalfLow[i])) {
    lowBoundaries <- c(lowBoundaries, HalfLow[i])
  } 
}
meanLow <- round(mean(lowBoundaries),2)
pdf(file = paste0("./Histograms/Histogram_",gradientNames[1],"Boundaries_Half.pdf"))
hist(as.numeric(lowBoundaries)
     , xlim = range(minGrad,maxGrad)
     , xlab = "Gradient"
     , col = "grey"
     , breaks = binBreaks
     , main = paste0("Boundary distribution for HALF ", gradientNames[1],": u=",as.character(meanLow)))
abline(v=meanLow
       , col = "red"
       , lwd = 2)
dev.off()
}

# Bloom
lowBoundaries <- c()
if (length(BloomLow) > 0){
for (i in 1:length(BloomLow)) {
  if (!is.na(BloomLow[i])) {
    lowBoundaries <- c(lowBoundaries, BloomLow[i])
  } 
}
meanLow <- round(mean(lowBoundaries),2)
pdf(file = paste0("./Histograms/Histogram_",gradientNames[1],"Boundaries_Bloom.pdf"))
hist(lowBoundaries
     , xlim = range(minGrad,maxGrad)
     , xlab = "Gradient"
     , breaks = binBreaks
     , col = "grey"
     , main = paste0("Boundary distribution for BLOOM ", gradientNames[1],": u=",as.character(meanLow)))
abline(v=meanLow
       , col = "red"
       , lwd = 2)
dev.off()
}


########################## High ####################################

# Now pull out all boundaries for fresh
allHigh <- c(highOnly$boundaries)
TrueHigh <- c(highOnlyTrue$boundaries)
HalfHigh <- c(highOnlyHalf$boundaries)
BloomHigh <- c(highOnlyBloom$boundaries)
nonBloomHigh <- c(HalfHigh,TrueHigh)

HighBoundaries <- c()
if (length(allHigh) > 0){
for (i in 1:length(allHigh)) {
  if (!is.na(allHigh[i])) {
    HighBoundaries <- c(HighBoundaries, allHigh[i])
  }
}
meanHigh <- round(mean(HighBoundaries),2)
pdf(file = paste0("./Histograms/Histogram_",gradientNames[3],"Boundaries_all.pdf"))
hist(HighBoundaries
     , xlim = range(minGrad,maxGrad)
     , xlab = "Gradient"
     , breaks = binBreaks
     , col = "grey"
     , main = paste0("Boundary distribution for ALL ", gradientNames[3],": u=",as.character(meanHigh)))
abline(v=meanHigh
       , col = "red"
       , lwd = 2)
dev.off()
}
HighBoundaries <- c()
if (length(nonBloomHigh) > 0){
for (i in 1:length(nonBloomHigh)) {
  if (!is.na(nonBloomHigh[i])) {
    HighBoundaries <- c(HighBoundaries, nonBloomHigh[i])
  }
}
meanHigh <- round(mean(HighBoundaries),2)
pdf(file = paste0("./Histograms/Histogram_",gradientNames[3],"Boundaries_nonBloom.pdf"))
hist(HighBoundaries
     , xlim = range(minGrad,maxGrad)
     , xlab = "Gradient"
     , breaks = binBreaks
     , col = "grey"
     , main = paste0("Boundary distribution for nonBlOOM ", gradientNames[3],": u=",as.character(meanHigh)))
abline(v=meanHigh
       , col = "red"
       , lwd = 2)
dev.off()
}


# True
HighBoundaries <- c()
if (length(TrueHigh) > 0){
for (i in 1:length(TrueHigh)) {
  if (!is.na(TrueHigh[i])) {
    HighBoundaries <- c(HighBoundaries, TrueHigh[i])
  }
}
meanHigh <- round(mean(HighBoundaries),2)
pdf(file = paste0("./Histograms/Histogram_",gradientNames[3],"Boundaries_True.pdf"))
hist(HighBoundaries
     , xlim = range(minGrad,maxGrad)
     , xlab = "Gradient"
     , breaks = binBreaks
     , col = "grey"
     , main = paste0("Boundary distribution for TRUE ", gradientNames[3],": u=",as.character(meanHigh)))
abline(v=meanHigh
       , col = "red"
       , lwd = 2)
dev.off()
}

# Half
HighBoundaries <- c()
if (length(HalfHigh) > 0) {
for (i in 1:length(HalfHigh)) {
  if (!is.na(HalfHigh[i])) {
    HighBoundaries <- c(HighBoundaries, HalfHigh[i])
  }
}
meanHigh <- round(mean(HighBoundaries),2)
pdf(file = paste0("./Histograms/Histogram_",gradientNames[3],"Boundaries_Half.pdf"))
hist(HighBoundaries
     , xlim = range(minGrad,maxGrad)
     , xlab = "Gradient"
     , breaks = binBreaks
     , col = "grey"
     , main = paste0("Boundary distribution for HALF ", gradientNames[3],": u=",as.character(meanHigh)))
abline(v=meanHigh
       , col = "red"
       , lwd = 2)
dev.off()
}

# Bloom
HighBoundaries <- c()
if (length(BloomHigh) > 0) {
for (i in 1:length(BloomHigh)) {
  if (!is.na(BloomHigh[i])) {
    HighBoundaries <- c(HighBoundaries, BloomHigh[i])
  }
}
meanHigh <- round(mean(HighBoundaries),2)
pdf(file = paste0("./Histograms/Histogram_",gradientNames[3],"Boundaries_Bloom.pdf"))
hist(HighBoundaries
     , xlim = range(minGrad,maxGrad)
     , xlab = "Gradient"
     , breaks = binBreaks
     , col = "grey"
     , main = paste0("Boundary distribution for BLOOM ", gradientNames[3],": u=",as.character(meanHigh)))
abline(v=meanHigh
       , col = "red"
       , lwd = 2)
dev.off()
}

# }
########################## Inter ####################################

# Now pull out all boundaries for fresh
allInter <- cbind(interOnly$boundaries, interOnly$boundariestwo)
TrueInter <- cbind(interOnlyTrue$boundaries,interOnlyTrue$boundariestwo)
HalfInter <- cbind(interOnlyHalf$boundaries,interOnlyHalf$boundariestwo)
BloomInter <- cbind(interOnlyBloom$boundaries, interOnlyBloom$boundariestwo)
HiInter <- cbind(interOnlyHi$boundaries, interOnlyHi$boundariestwo)
LoInter <- cbind(interOnlylo$boundaries, interOnlylo$boundariestwo)
nonBloomInter <- rbind(HalfInter,TrueInter)

InterBoundaries <- c()
if (length(allInter) > 0){
InterBoundaries <- allInter
meanInterOne <- round(mean(InterBoundaries[,1]),2)
meanInterTwo <- round(mean(InterBoundaries[,2]),2)
maximumY <- max(max(table(InterBoundaries[,1])),max(table(InterBoundaries[,2])))*1.5
pdf(file = paste0("./Histograms/Histogram_",gradientNames[2],"Boundaries_all.pdf"))
hist(InterBoundaries[,1]
     , xlim = range(minGrad,maxGrad)
     , xlab = "Gradient"
     , breaks = binBreaks
     , ylim = c(0,maximumY)
     , col = rgb(0,0,1,0.5)
     , main = paste0("Boundary distribution for ALL ", gradientNames[2],": u=",as.character(meanInterOne),",",as.character(meanInterTwo)))
abline(v=meanInterOne
       , col = "blue"
       , lwd = 2)
hist(InterBoundaries[,2]
     , xlim = range(minGrad,maxGrad)
     , breaks = binBreaks
     , ylim = c(0,maximumY)
     , col = rgb(1,0,0,0.5)
     , add = TRUE)
abline(v=meanInterTwo
       , col = "red"
       , lwd = 2)
dev.off()
}
InterBoundaries <- c()
if (length(nonBloomInter) > 0){
for (i in 1:length(nonBloomInter)) {
InterBoundaries <- nonBloomInter
meanInterOne <- round(mean(InterBoundaries[,1]),2)
meanInterTwo <- round(mean(InterBoundaries[,2]),2)
maximumY <- max(max(table(InterBoundaries[,1])),max(table(InterBoundaries[,2])))*1.5
pdf(file = paste0("./Histograms/Histogram_",gradientNames[2],"Boundaries_nonBloom.pdf"))
hist(InterBoundaries[,1]
     , xlim = range(minGrad,maxGrad)
     , xlab = "Gradient"
     , ylim = c(0,maximumY)
     , col = rgb(0,0,1,0.5)
     , breaks = binBreaks
     , main = paste0("Boundary distribution for nonBLOOM ", gradientNames[2],": u=",as.character(meanInterOne),",",as.character(meanInterTwo)))
abline(v=meanInterOne
       , col = "blue"
       , lwd = 2)
hist(InterBoundaries[,2]
     , xlim = range(minGrad,maxGrad)
     , ylim = c(0,maximumY)
     , breaks = binBreaks
     , col = rgb(1,0,0,0.5)
     , add = TRUE)
abline(v=meanInterTwo
       , col = "red"
       , lwd = 2)
dev.off()
}
}


# True
InterBoundaries <- c()
if (length(TrueInter) > 0){
InterBoundaries <- TrueInter
meanInterOne <- round(mean(InterBoundaries[,1]),2)
meanInterTwo <- round(mean(InterBoundaries[,2]),2)
maximumY <- max(max(table(InterBoundaries[,1])),max(table(InterBoundaries[,2])))*1.5
pdf(file = paste0("./Histograms/Histogram_",gradientNames[2],"Boundaries_True.pdf"))
hist(InterBoundaries[,1]
     , xlim = range(minGrad,maxGrad)
     , xlab = "Gradient"
     , ylim = c(0,maximumY)
     , breaks = binBreaks
     , col = rgb(0,0,1,0.5)
     , main = paste0("Boundary distribution for TRUE ", gradientNames[2],": u=",as.character(meanInterOne),",",as.character(meanInterTwo)))
abline(v=meanInterOne
       , col = "blue"
       , lwd = 2)
hist(InterBoundaries[,2]
     , xlim = range(minGrad,maxGrad)
     , ylim = c(0,maximumY)
     , breaks = binBreaks
     , col = rgb(1,0,0,0.5)
     , add = TRUE)
abline(v=meanInterTwo
       , col = "red"
       , lwd = 2)
dev.off()
}

# Half
InterBoundaries <- c()
if (length(HalfInter) > 0){
InterBoundaries <- HalfInter
meanInterOne <- round(mean(InterBoundaries[,1]),2)
meanInterTwo <- round(mean(InterBoundaries[,2]),2)
maximumY <- max(max(table(InterBoundaries[,1])),max(table(InterBoundaries[,2])))*1.5
pdf(file = paste0("./Histograms/Histogram_",gradientNames[2],"Boundaries_Half.pdf"))
hist(InterBoundaries[,1]
     , xlim = range(minGrad,maxGrad)
     , xlab = "Gradient"
     , ylim = c(0,maximumY)
     , breaks = binBreaks
     , col = rgb(0,0,1,0.5)
     , main = paste0("Boundary distribution for Half ", gradientNames[2],": u=",as.character(meanInterOne),",",as.character(meanInterTwo)))
abline(v=meanInterOne
       , col = "blue"
       , lwd = 2)
hist(InterBoundaries[,2]
     , xlim = range(minGrad,maxGrad)
     , ylim = c(0,maximumY)
     , breaks = binBreaks
     , col = rgb(1,0,0,0.5)
     , add = TRUE)
abline(v=meanInterTwo
       , col = "red"
       , lwd = 2)
dev.off()
}

# Bloom
InterBoundaries <- c()
if (length(BloomInter) > 0){
InterBoundaries <- BloomInter
meanInterOne <- round(mean(InterBoundaries[,1]),2)
meanInterTwo <- round(mean(InterBoundaries[,2]),2)
maximumY <- max(max(table(InterBoundaries[,1])),max(table(InterBoundaries[,2])))*1.5
pdf(file = paste0("./Histograms/Histogram_",gradientNames[2],"Boundaries_Bloom.pdf"))
hist(InterBoundaries[,1]
     , xlim = range(minGrad,maxGrad)
     , xlab = "Gradient"
     , ylim = c(0,maximumY)
     , breaks = binBreaks
     , col = rgb(0,0,1,1)
     , main = paste0("Boundary distribution for BLOOM ", gradientNames[2],": u=",as.character(meanInterOne),",",as.character(meanInterTwo)))
abline(v=meanInterOne
       , col = "blue"
       , lwd = 2)
hist(InterBoundaries[,2]
     , xlim = range(minGrad,maxGrad)
     , ylim = c(0,maximumY)
     , breaks = binBreaks
     , col = rgb(1,0,0,0.5)
     , add = TRUE)
abline(v=meanInterTwo
       , col = "red"
       , lwd = 2)
dev.off()
}

# Hi
InterBoundaries <- c()
if (length(HiInter) > 0){
InterBoundaries <- HiInter
meanInterOne <- round(mean(InterBoundaries[,1]),2)
meanInterTwo <- round(mean(InterBoundaries[,2]),2)
maximumY <- max(max(table(InterBoundaries[,1])),max(table(InterBoundaries[,2])))*1.5
pdf(file = paste0("./Histograms/Histogram_",gradientNames[2],"Boundaries_Hi.pdf"))
hist(InterBoundaries[,1]
     , xlim = range(minGrad,maxGrad)
     , xlab = "Gradient"
     , ylim = c(0,maximumY)
     , breaks = binBreaks
     , col = rgb(0,0,1,0.5)
     , main = paste0("Boundary distribution for HI ", gradientNames[2],": u=",as.character(meanInterOne),",",as.character(meanInterTwo)))
abline(v=meanInterOne
       , col = "blue"
       , lwd = 2)
hist(InterBoundaries[,2]
     , xlim = range(minGrad,maxGrad)
     , ylim = c(0,maximumY)
     , breaks = binBreaks
     , col = rgb(1,0,0,0.5)
     , add = TRUE)
abline(v=meanInterTwo
       , col = "red"
       , lwd = 2)
dev.off()
}

# Lo
InterBoundaries <- c()
if (length(LoInter) > 0){
InterBoundaries <- LoInter
meanInterOne <- round(mean(InterBoundaries[,1]),2)
meanInterTwo <- round(mean(InterBoundaries[,2]),2)
maximumY <- max(max(table(InterBoundaries[,1])),max(table(InterBoundaries[,2])))*1.5
pdf(file = paste0("./Histograms/Histogram_",gradientNames[2],"Boundaries_Lo.pdf"))
hist(InterBoundaries[,1]
     , xlim = range(minGrad,maxGrad)
     , xlab = "Gradient"
     , ylim = c(0,maximumY)
     , breaks = binBreaks
     , col = rgb(0,0,1,0.5)
     , main = paste0("Boundary distribution for LO ", gradientNames[2],": u=",as.character(meanInterOne),",",as.character(meanInterTwo)))
abline(v=meanInterOne
       , col = "blue"
       , lwd = 2)
hist(InterBoundaries[,2]
     , xlim = range(minGrad,maxGrad)
     , ylim = c(0,maximumY)
     , breaks = binBreaks
     , col = rgb(1,0,0,0.5)
     , add = TRUE)
abline(v=meanInterTwo
       , col = "red"
       , lwd = 2)
dev.off()
}


########################## TAXASUMMARIES SECTION ####################################
print("In taxa summaries section")

# library(outliers)
# library(stats)
output <- "TaxaSummaries"

system(paste0("mkdir ",output))
# system(paste0("biom convert -i ",biom," -o ./",output,"/OTUTable_text.txt --to-tsv --header-key taxonomy"))
############## Pre-amble; loading files, formatting, relative abund ###################

taxa <- read.delim(paste0(biomPWD)
                   , strip.white = TRUE
                   , stringsAsFactors = FALSE
                   , header = TRUE
                   , row.names = 1
                   # , skip = 1
                   )
colnames(taxa) <- gsub(".","-", colnames(taxa), fixed = TRUE)

# # Order by abundance
# orderAbund <- order(rowSums(taxa[,-ncol(taxa)]), decreasing = TRUE)
# taxa <- taxa[orderAbund,]
# OR order by alphabetical
orderAlpha <- order(taxa[,ncol(taxa)])
taxa <- taxa[orderAlpha,]
# Make taxonomy ref
taxonomyRef <- data.frame(taxa[,ncol(taxa)])
rownames(taxonomyRef) <- rownames(taxa)

# Make taxa summary separate
taxa <- taxa[,-(ncol(taxa))]

# Make taxa summary relative abundance
relAbund <- function(x) {
  total <- sum(x)
  newx <- x/total
  return(newx)
}

# Make get taxa names function
getTaxaNames <- function(x, delim) {
    tempname <- strsplit(as.character(x),split = paste0(delim))[[1]][c(3,6,7)]
    if (length(grep("Proteobacteria",tempname[1])) >0) {
        tempname[1] <- strsplit(as.character(x),split=paste0(delim))[[1]][4]
    }
    tempname <- gsub("^.*__","",tempname)
    return(paste0(tempname[1],": ",tempname[2],"_", tempname[3]))
}

# Apply function to taxa
taxaAbund <- data.frame(apply(taxa, 2, FUN = relAbund))
names(taxaAbund) <- names(taxa)

#--

# Load metadata
metadata <- read.delim(paste0(metadataPWD), stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""), row.names = 1)
rownames(metadata) <- gsub(".","-",rownames(metadata), fixed = TRUE)

# Make sure only samples in metadata are ones in taxa summaries
metadata <- metadata[which(rownames(metadata) %in% names(taxaAbund)),]

# Use re-order with gradientNames[4] to make in order of gradient
metadata <- metadata[order(as.numeric(metadata[[paste0(gradientNames[4])]])),]
# then make taxa the same
taxaAbund <- taxaAbund[,rownames(metadata)]

# Load modelBoundaries
modelBoundaries <- read.delim(paste0(modelboundariesPWD)
                              , header = TRUE
                              , row.names = 1
                              , strip.white = TRUE
                              , stringsAsFactors = FALSE
                              , na.strings = c('','na','NA'))

# Then, extract gradient values for these
ordered.gradient <- metadata[[paste0(gradientNames[4])]]

# Now, I want to replace the headers of the euktaxa.ordered with these values.
taxaAbund.grad <- taxaAbund
names(taxaAbund.grad) <- ordered.gradient


# Now, I want to take these 'trend' lines and plot them all on single plots so I can see general trends.

# Changing the dataframe into a matrix so that I can plug it directly into barplot. 
taxaAbund.matrix <- as.matrix(taxaAbund.grad, stringsAsFactors = FALSE, header = TRUE, rownames.force = TRUE)

# Want to look at ONLY things that are grthan the 0.05 threshold
# Create an empty vector
grthanTF <- vector()
# Now, for each row in taxa.matrix, I want to test whether ANY value in that line is greater than 0.05.
# I take this vector (T/F,T/F,T/F) and ask whether any of them are true; if one of them are, then that line is saved to the vector grthan.
for (i in 1:nrow(taxaAbund.matrix)) {
  line <- as.factor(taxaAbund.matrix[i,] >= as.numeric(minthreshold))
  if (any(line == "TRUE")) {
    grthanTF[i] <- TRUE
  } else{
    grthanTF[i] <- FALSE
  }
}


# Get only major taxa
taxaAbund.grad.filt <- taxaAbund.grad[grthanTF,]

######## Making barplot; taxa summaries ##########
print("Making barplots")

# For the data matrix, I want to colour the bars so that all taxanomic groups that represent more than 5% of the total OTUs found in any sample are coloured. Note that this is slightly different than in the paper; they claim to use a threshold of 1%, but when I looked at the data it seemed they were wrong. There were samples with more than 1% (but less than 5%) representation that they did not label. Thus, I will be using the 5% threshold.

# Then, I use the vector 'grthan' to do some things:

# Try to do a calculation-based way to space out taxa so they are spread out by gradient.
# First, list all unique salinities
unique.num <- unique(names(taxaAbund.grad.filt))

# # Get random colours
set.seed(1989)
colors.lines <- sample(colors()[grep('gr(a|e)y|white|snow', colors(), invert = TRUE)], nrow(taxaAbund.grad.filt), replace = FALSE)

# Then, count how many of each unique gradient there is. 
# This loop counts how many reps there are of each 'unique' gradient in the taxa file
spacing <- data.frame()
for (i in 1:length(unique.num)) {
  spacing[i,1] <- as.numeric(unique.num[i])
  spacing[i,2] <- sum(unique.num[i] == names(taxaAbund.grad.filt))
}
names(spacing) <- c("Sal","count")

# Let us imagine that each "gradient" is worth 1 unit. each separation of gradient, then, is also worth 1 unit. 
# Therefore, the barwidth should be calculated as 1/#ofsites@thesamegradient, then multiplied by a constant to make the average bar width 1
# We need the average bar width to be 1 because the 'space' in the barplot function is calculated as a multiple of the average bar width. Thus, we want a '1' 'to actually equal '1'
# The distance between clusters should be calculated as 1*difference, where the differences are scaled to the average bar width.


# Calculating proportional band widths.
spacing.width <- vector()
for (i in 1:length(spacing[,2])) {
  
  spacing.width <- c(spacing.width, rep((1/spacing[i,2]), spacing[i,2]))
  
}
# And now multiplying by 1/barwidth so the mean bar width is one

spacing.wid <- spacing.width*(1/mean(spacing.width))

# Calculating distances
spacing.num <- vector()
for (i in 1:length(spacing[,2])) {
  if (i == 1) {
    spacing.num <- c(0, rep(0, spacing[i,2]-1))
  } else {
    spacing.num <- c(spacing.num, (spacing[i,1]-spacing[i-1,1])*1 , rep(0, spacing[i,2]-1))
  }
}

# COLORS
grthan <- gsub(FALSE,NA,grthanTF)
grthan <- as.logical(grthan)

# I make a factor that tells me which row names should have coloured bars. Bars with less than 5% representation across all samples are listed as "NA".
# coloredbars <- as.factor(rownames((taxaAbund.matrix[grthan,])))
coloredbarsvec <- as.vector(rownames((taxaAbund.matrix[grthan,])))
coloredbars <- coloredbarsvec[-which(is.na(coloredbarsvec))]

# taxonomyRef[coloredbars,]
# coloredbars <- as.factor(as.vector(rownames((taxaAbund.matrix[grthan,]))))

#get taxa names function
getTaxaNames <- function(x, delim) {
    tempname <- strsplit(as.character(x),split = paste0(delim))[[1]][c(3,6,7)]
    # if (length(grep("Proteobacteria",tempname[1])) >0) {
    #     tempname[1] <- strsplit(as.character(x),split=paste0(delim))[[1]][4]
    # }
    tempname <- gsub("^.*__","",tempname)
    return(paste0(tempname[1],": ",tempname[2],"_", tempname[3]))
}

# get names of things with colored bars
# namesOTU <- levels(coloredbars)
# namesTaxa <- taxonomyRef[match(namesOTU,rownames(taxonomyRef)),]
# abbrTaxa <- sapply(namesTaxa, function(x) {getTaxaNames(x,delim="; ")})
namesTaxa <- taxonomyRef[match(coloredbars,rownames(taxonomyRef)),]
abbrTaxa <- sapply(namesTaxa, function(x) {getTaxaNames(x,delim="; ")})

# get alphabetical order
# alphaOrder <- order(abbrTaxa)

# I make a factor that tells me how many different coloured bars there are. Since the factor above has NA values as well, I count the ones that are NOT "NA".
ncoloredbars <-sum(!is.na(coloredbars))

# assign colors
colorsassigned <- data.frame("color"=colors.lines, "bars"=coloredbars, stringsAsFactors =FALSE)
# assign colors for barplot
colorsPlot <- rbind(colorsassigned,c("white","less"))[match(rownames(taxaAbund.matrix),colorsassigned[,2], nomatch = (nrow(colorsassigned)+1)),1]

# version with adjust bar widths
pdf(paste0("./",output,"/TaxaSummariesPartitioned_byOTU.pdf"), width = 10, height = 5)
# quartz(,10,5)
par(mar=c(5.1,4.1,3,1), xpd=TRUE)
barplot(taxaAbund.matrix 
        , col = colorsPlot
        # , col = c(colors.lines)[coloredbars] # Using the generated colours from above, I colour each bar that has >5% representation a different colour.
        , space = spacing.num # This lines the bar plot up with the boxplot in spacing
        , cex.axis = 0.6 # scaling axis
        , cex.lab = 1 # scaling label
        , ylab = "% Representation" # y axis label
        , xlab = "Salinity"
        , xaxt = 'n'
        , width = spacing.wid
        , xlim = c(0,sum(spacing.wid,spacing.num))
        , border = NA
        , main = 'Taxa Summaries; partitioned (by OTU)'
        , sub = paste0('Minimum max abund = ', minthreshold)
        
)
# Add manual x axis
mindif35 <- as.numeric(colnames(taxaAbund.matrix)) - 35
which35 <- which.min(abs(mindif35))
scalepos <- mindif35[which35]
axis(side = 1
     , at = seq(0,sum(spacing.wid[1:which35],spacing.num[1:which35])-scalepos
                # , by = sum(spacing.wid,spacing.num)/(30/5)
                , length.out=35/5 +1)
     , labels = seq(0,35,by = 5))
dev.off()


# not adjusted; even spaces
pdf(paste0("./",output,"/TaxaSummarieseven_byOTU.pdf"), width = 10, height = 5)
# quartz(,10,5)
par(mar=c(5.1,4.1,3,1), xpd=TRUE)
barplot(taxaAbund.matrix 
        , col = colorsPlot
        , cex.axis = 0.6 # scaling axis
        , cex.lab = 1 # scaling label
        , ylab = "% Representation" # y axis label
        , xlab = "Salinity"
        , xaxt = 'n'
        , border = NA
        # , space = spacing.num
        , main = 'Taxa Summaries; partitioned (by OTU)'
        , sub = paste0('Minimum max abund = ', minthreshold)
        
)
# Add manual x axis: space out properly
axis(side = 1
     , at = seq(0.6
                # ,ncol(taxaAbund.matrix)-0.5
                , length.out=ncol(taxaAbund.matrix)
                , by=1.2)
     , labels = colnames(taxaAbund.matrix)
     , cex.axis = 0.5
     , line=-1
     , tick=FALSE
     )
dev.off()

# not adjusted; clustered
pdf(paste0("./",output,"/TaxaSummariesevenclustered_byOTU.pdf"), width = 10, height = 5)
# quartz(,10,5)
par(mar=c(5.1,4.1,3,1), xpd=TRUE)
barplot(taxaAbund.matrix 
        , col = colorsPlot
        , cex.axis = 0.6 # scaling axis
        , cex.lab = 1 # scaling label
        , ylab = "% Representation" # y axis label
        , xlab = "Salinity"
        , xaxt = 'n'
        , border = NA
        , space = spacing.num
        , main = 'Taxa Summaries; partitioned (by OTU)'
        , sub = paste0('Minimum max abund = ', minthreshold)
        
)
# Add manual x axis: space out properly
totalDist <- 0:length(spacing.num)
for ( i in 1:length(spacing.num) ) {
    totalDist[i+1] <- totalDist[i] + spacing.num[i] + 1
}
totalDist <- totalDist-0.5
axis(side = 1
     , at = totalDist[-1]
     , labels = colnames(taxaAbund.matrix)
     , cex.axis = 0.5
     , line=-1
     , tick=FALSE
)
dev.off()

#### MAKE LEGEND HERE

# Type out legend labels because I want them a certain way

# This makes an 'empty' plot with no labels or axis.
pdf(paste0("./",output,"/Legend_byOTU_partitionedTaxaSummaries.pdf"), 15,15)
plot(0,0,xlab="",ylab="",xaxt="n",yaxt="n",bty="n",pch="")
legend("center"
       , legend = c(abbrTaxa, "Other (<5%)")
       , pch = 22 # squares
       , col= "black" # black outline
       , pt.bg = c(as.character(colorsassigned[,1]), "white") #fill squares with colours
       , bty = "n" # Get rid of outer box
       , title = "Taxa"
       , cex = 1
       , xpd = TRUE
)
dev.off()

########### Taxa summaries partitioned; by class ###########

combineGroupAbund <- function(g,allOTUnames,groups,taxaAbund.matrix) {
    tomerge <- allOTUnames[which(groups==g)]
    if (length(tomerge) > 1) {
        return(colSums(taxaAbund.matrix[tomerge,]))
    } else {
        return(taxaAbund.matrix[tomerge,])
    }
}

levelNames <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
for ( level in 3:5) {
    allOTUnames <- rownames(taxaAbund.matrix)
    allTaxanames <- taxonomyRef[allOTUnames,]
    groups <- c()
    for ( i in allTaxanames) {
        grouptemp <- strsplit(as.character(i),split="; ")[[1]][level]
        # if ( ( level == 3 ) & (length(grep("Proteobacteria",grouptemp))) ){
        #     grouptemp <-  strsplit(as.character(i),split="; ")[[1]][4]
        # }
        leveltemp <- level-1
        while (is.na(grouptemp)) {
            grouptemp <- strsplit(as.character(i),split="; ")[[1]][leveltemp]
            leveltemp <- leveltemp-1
        }
        grouptemp <- gsub("^.*__","",grouptemp)
        groups <- c(groups,grouptemp)
    }
    uniqueGroups <- unique(groups)
    newMatrix <- matrix(ncol=ncol(taxaAbund.matrix), nrow=length(unique(groups)), dimnames = list(uniqueGroups,colnames(taxaAbund.matrix)))
    for ( g in uniqueGroups ) {
        newMatrix[g,] <- combineGroupAbund(g,allOTUnames,groups,taxaAbund.matrix)
    }
    
    #Now, filter by abundance
    groupBelowThresh <- c()
    for ( r in 1:nrow(newMatrix) ) {
        if (max(newMatrix[r,]) < minthreshold) {
            groupBelowThresh <- c(groupBelowThresh,r)
        }
    }
    #get colors
    set.seed(1993)
    groupsWithColors <- uniqueGroups[-groupBelowThresh]
    getColors <-sample(colors()[-grep("gr[a|e]y|white", colors())],length(groupsWithColors), replace=FALSE)
    colorsGroup <- rep("white",nrow(newMatrix))
    #make things that that need color colored
    rowsNeedColors <- match(groupsWithColors,rownames(newMatrix))
    for ( r in 1:length(rowsNeedColors) ) {
        row <- rowsNeedColors[r]
        colorsGroup[row] <- getColors[r]
    }
    
    pdf(paste0("./",output,"/TaxaSummariesPartitioned_by",levelNames[level],".pdf"), width = 10, height = 5)
    # quartz(,10,5)
    par(mar=c(5.1,4.1,3,1), xpd=TRUE)
    barplot(newMatrix 
            , col = colorsGroup # Using the generated colours from above, I colour each bar that has >5% representation a different colour.
            , space = spacing.num # This lines the bar plot up with the boxplot in spacing
            , cex.axis = 0.6 # scaling axis
            , cex.lab = 1 # scaling label
            , ylab = "% Representation" # y axis label
            , xlab = "Salinity"
            , xaxt = 'n'
            , width = spacing.wid
            , xlim = c(0,sum(spacing.wid,spacing.num))
            , border = NA
            , main = 'Taxa Summaries; partitioned (by OTU)'
            , sub = paste0('Minimum max abund = ', minthreshold)
            
    )
    # Add manual x axis
    mindif35 <- as.numeric(colnames(newMatrix)) - 35
    which35 <- which.min(abs(mindif35))
    scalepos <- mindif35[which35]
    axis(side = 1
         , at = seq(0,sum(spacing.wid[1:which35],spacing.num[1:which35])-scalepos
                    # , by = sum(spacing.wid,spacing.num)/(30/5)
                    , length.out=35/5 +1)
         , labels = seq(0,35,by = 5))
    # 
    # pos30 <- max(which(colnames(taxaAbund.matrix) == "30"))
    # axis(side = 1
    #      , at = seq(0,sum(spacing.wid[1:pos30],spacing.num[1:pos30])
    #                 # , by = sum(spacing.wid,spacing.num)/(30/5)
    #                 , length.out=30/5+1)
    #      , labels = seq(0,30,by = 5))
    dev.off()
    
    # Make a legend for this
    
    #get groups that were > 5
    namesForLegend <- uniqueGroups[-groupBelowThresh]
    colorsForLegend <- colorsGroup[-groupBelowThresh]
    # This makes an 'empty' plot with no labels or axis.
    pdf(paste0("./",output,"/Legend_by",levelNames[level],"_partitionedTaxaSummaries.pdf"), 15,15)
    plot(0,0,xlab="",ylab="",xaxt="n",yaxt="n",bty="n",pch="")
    legend("center"
           , legend = c(namesForLegend, "Other (<5%)")
           , pch = 22 # squares
           , col= "black" # black outline
           , pt.bg = c(as.character(colorsForLegend), "white") #fill squares with colours
           , bty = "n" # Get rid of outer box
           , title = "Taxa"
           , cex = 1
           , xpd = TRUE
    )
    dev.off()
    
    
    # not adjusted; even spaces
    pdf(paste0("./",output,"/TaxaSummarieseven_by",levelNames[level],".pdf"), width = 10, height = 5)
    # quartz(,10,5)
    par(mar=c(5.1,4.1,3,1), xpd=TRUE)
    barplot(newMatrix 
            , col = colorsGroup # Using the generated colours from above, I colour each bar that has >5% representation a different colour.
            # , space = spacing.num # This lines the bar plot up with the boxplot in spacing
            , cex.axis = 0.6 # scaling axis
            , cex.lab = 1 # scaling label
            , ylab = "% Representation" # y axis label
            , xlab = "Salinity"
            , xaxt = 'n'
            , border = NA
            , main = 'Taxa Summaries; partitioned (by OTU)'
            , sub = paste0('Minimum max abund = ', minthreshold)
            
    )
    # Add manual x axis: space out properly
    axis(side = 1
         , at = seq(0.6
                    # ,ncol(taxaAbund.matrix)-0.5
                    , length.out=ncol(newMatrix)
                    , by=1.2)
         , labels = colnames(newMatrix)
         , cex.axis = 0.5
         , line=-1
         , tick=FALSE
    )
    dev.off()
    
    # not adjusted; clustered
    pdf(paste0("./",output,"/TaxaSummariesevenclustered_by",levelNames[level],".pdf"), width = 10, height = 5)
    # quartz(,10,5)
    par(mar=c(5.1,4.1,3,1), xpd=TRUE)
    barplot(newMatrix 
            , col = colorsGroup
            , cex.axis = 0.6 # scaling axis
            , cex.lab = 1 # scaling label
            , ylab = "% Representation" # y axis label
            , xlab = "Salinity"
            , xaxt = 'n'
            , border = NA
            , space = spacing.num
            , main = 'Taxa Summaries; partitioned (by OTU)'
            , sub = paste0('Minimum max abund = ', minthreshold)
            
    )
    # Add manual x axis: space out properly
    totalDist <- 0:length(spacing.num)
    for ( i in 1:length(spacing.num) ) {
        totalDist[i+1] <- totalDist[i] + spacing.num[i] + 1
    }
    totalDist <- totalDist-0.5
    axis(side = 1
         , at = totalDist[-1]
         , labels = colnames(newMatrix)
         , cex.axis = 0.5
         , line=-1
         , tick=FALSE
    )
    dev.off()
    
    
}

####### OTU TOLERANCE RANGES ##########
# 

otuBoundaries <- modelBoundaries %>%
  mutate(OTU=as.character(rownames(modelBoundaries))) %>%
  filter(typeSimple!="noclass") %>%
  mutate(boundaries=ifelse(is.na(boundaries), minGrad,boundaries)
         , boundariestwo=ifelse(is.na(boundariestwo), maxGrad,boundariestwo)
         , midBoundary=(boundaries+boundariestwo)/2
         , breadth=abs(c(boundaries-boundariestwo))) %>%
  select(OTU, typeSimple,type,boundaries, boundariestwo, midBoundary, breadth)
# nOTUs <- nrow(otuBoundaries)
# lwd_adj <- 100/nOTUs
nOTUs <- toleranceRanges %>%
  left_join(otuBoundaries) %>%
  filter(!is.na(typeSimple), abundances_nozero !=0, position=="withinBounds") %>% select(OTU)%>%pull() %>%unique() %>% length()
outliers_points <- toleranceRanges %>%
  left_join(otuBoundaries) %>%
  filter(!is.na(typeSimple), abundances_nozero !=0, position=="outlier") %>%
  mutate(typeSimple=gsub("Restricted","",typeSimple)) %>%
  mutate(typeSimple=factor(typeSimple, levels=c("fresh","brackish","marine"))) %>%
  arrange(typeSimple, midBoundary) %>%
  mutate(OTU=factor(OTU, levels=unique(OTU)))
ggsave("tolerance_ranges.pdf", height=10, width=6.5
       ,toleranceRanges %>%
         left_join(otuBoundaries) %>%
         filter(!is.na(typeSimple), abundances_nozero !=0, position=="withinBounds") %>%
         mutate(typeSimple=gsub("Restricted","",typeSimple)) %>%
         mutate(typeSimple=factor(typeSimple, levels=c("fresh","brackish","marine"))) %>%
         arrange(typeSimple, midBoundary) %>%
         mutate(OTU=factor(OTU, levels=unique(OTU))) %>%
         ggplot() + geom_segment(aes(x=startSal, xend=endSal, y=OTU, yend=OTU, col=typeSimple, alpha=abundances_nozero_maxstand), lwd=100/(nOTUs*1.2)) +
         geom_point(data=outliers_points,aes(x=startSal, y=OTU,col=typeSimple, alpha=abundances_nozero_maxstand), cex=1/(nOTUs*1.2) )+
         scale_color_manual(values=c("blue","purple","red")) +
         labs(alpha = paste0("Relative abundance \n(standardized by max OTU)"), col="Specialist Type") +
         theme_classic() +
         theme(axis.text.y = element_blank(), axis.ticks.y=element_blank()) + xlab("Salinity")
       )
# use taxaAbund to get 'extreme' points as well.

##### CDFs #####

### Boundary CDFs ###
allBoundaries <- c(modelBoundaries$boundaries, modelBoundaries$boundariestwo)
freshBoundaries <- c(modelBoundaries$boundaries[modelBoundaries$typeSimple=="freshRestricted"]
                   ,modelBoundaries$boundariestwo[modelBoundaries$typeSimple=="freshRestricted"])
brackBoundaries1  <- c(modelBoundaries$boundaries[modelBoundaries$typeSimple=="brackishRestricted"])
brackBoundaries2  <- c(modelBoundaries$boundariestwo[modelBoundaries$typeSimple=="brackishRestricted"])
marineBoundaries <-  c(modelBoundaries$boundaries[modelBoundaries$typeSimple=="marineRestricted"]
                     ,modelBoundaries$boundariestwo[modelBoundaries$typeSimple=="marineRestricted"])

# Overall CDF
boundaries.cdf <- ecdf(allBoundaries)
Salinity <- seq(0,40,length.out = length(allBoundaries))
`F(Salinity)` <- boundaries.cdf(v=Salinity)

# Fresh CDF
boundaries.fresh.cdf <- ecdf(freshBoundaries)
fresh.b.cdf <- boundaries.fresh.cdf(v=Salinity)

# Brack CDF
boundaries.brack1.cdf <- ecdf(brackBoundaries1)
brack1.b.cdf <- boundaries.brack1.cdf(v=Salinity)
boundaries.brack2.cdf <- ecdf(brackBoundaries2)
brack2.b.cdf <- boundaries.brack2.cdf(v=Salinity)
# Marine CDF
boundaries.marine.cdf <- ecdf(marineBoundaries)
marine.b.cdf <- boundaries.marine.cdf(v=Salinity)


pdf("CDF_boundaries.pdf")
plot(`F(Salinity)` ~ Salinity, type="l", lwd=2, ylim=c(0,1))
lines(fresh.b.cdf ~ Salinity, col="blue", lty=2, lwd=2)
lines(brack1.b.cdf ~ Salinity, col="purple4", lty=2, lwd=2)
lines(brack2.b.cdf ~ Salinity, col="magenta", lty=2, lwd=2)
lines(marine.b.cdf ~ Salinity, col="red", lty=2, lwd=2)
legend("bottomright", legend = c("All","Fresh (end)","Brackish (start)", "Brackish (end)", "Marine (start)"), lty=c(1,3,3,3,3), col=c("black","blue","purple4","magenta","blue"), bty="n")
dev.off()

