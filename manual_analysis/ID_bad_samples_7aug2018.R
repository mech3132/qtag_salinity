#!/bin/bash
library(MASS) # for NMDS plotting (isoMDS)
library(vegan) # bray curtis
##### This script takes otutables an stuff from salinity project and plots them to see if there are differences in extr/seq methods #######

#### Set FP #####
# setwd("/Users/melissachen/Documents/Masters/Project_Environmental/FromBotaClust_1feb2018/Removing_contamined_samples/")
# dm16sFP <- "/beta_fr16/bray_curtis_fraser16S_otu_table_raw.txt"
# dm18sFP <- "beta_fr18/bray_curtis_OTU_Table_wtaxa.txt"
otu16FP <- "raw_otu_and_mf/Fraser_16s/OUTPUT_FILES/OTUTable_text.txt"
otu18FP <- "raw_otu_and_mf/Fraser_18S/OUTPUT_FILES/otu_table_temp.txt"
mf16FP <- "raw_otu_and_mf/Fraser_16S/MANUAL_INPUT_FILES/Fraser16s_mappingfile_merged.txt"
mf18FP <- "raw_otu_and_mf/Fraser_18S/MANUAL_INPUT_FILES/MF_Fraser18s_control_readded.txt"
# badsamp <- c("E-TB-2016-2","E-BA-2016-3","E-BL-2016-1")
badsamp <- c("E.TB.2016.2","E.BA.2016.3","E.BL.2016.1")
badsampreps <- c("E-TB-2016","E-BA-2016","E-BL-2016")

badsamp2 <- c("S-st-1-6-19-r3","S-st-1-6-19-r2","S-st-1-6-19-r1"
              # ,"S-st-3-5-27-r3","S-st-3-5-27-r2","S-st-3-5-27-r1"
              # ,"E-TB-2016-2","E-BA-2016-3","E-BL-2016-1"
              )
badsamp2 <- gsub("-",".", badsamp2, fixed=TRUE)
badsampColors <- c("red","green","purple")
badsampSize <- 2

#############

# dm16 <- read.delim(dm16sFP, header=TRUE, stringsAsFactors = FALSE, row.names=1)
# dm18 <- read.delim(dm18sFP, header=TRUE, stringsAsFactors = FALSE, row.names=1)

otu16 <- read.delim(otu16FP, header=TRUE, skip=1, stringsAsFactors = FALSE, row.names = 1)
otu16_notax <- otu16[,colnames(otu16)!="taxonomy"]
otu18<- read.delim(otu18FP, header=TRUE, skip=1, stringsAsFactors = FALSE, row.names = 1)
otu18_notax <- otu18[,colnames(otu18)!="taxonomy"]

mf16 <- read.delim(mf16FP, header=TRUE, stringsAsFactors = FALSE, row.names=1)
rownames(mf16) <- gsub("-",".", rownames(mf16))
mf18 <- read.delim(mf18FP, header=TRUE, stringsAsFactors = FALSE, row.names=1)
rownames(mf18) <- gsub("-",".", rownames(mf18))

mf16 <- mf16[match(colnames(otu16), rownames(mf16)),]
mf18 <- mf18[match(colnames(otu18), rownames(mf18)),]


dm16 <- vegdist(t(otu16_notax), method = "bray")
dm18 <- vegdist(t(otu18_notax), method = "bray")

set.seed(1234012934)
nmds16 <- isoMDS(dist(dm16), k=2)
nmds18 <- isoMDS(dist(dm18), k=2)

# Plot via colors

salGrad <- colorRampPalette(colors = c("white","blue","darkblue"))
colorSal <- salGrad(34)
salRound <- round(as.numeric(mf16$SalinityEnviron))
salRound18s <- round(as.numeric(mf18$SalinityEnviron))

#### NMDS 16s #####

# get things from replicate of bad samples
outline <- rep("black",length=nrow(nmds16$points))
size <- rep(1, length=nrow(nmds16$points))
# highlight bad samples' replicates
badrepsTB <- rownames(mf16)[ mf16$ColRep %in% badsampreps[1]]
badrepsBA <- rownames(mf16)[ mf16$ColRep %in% badsampreps[2]]
badrepsBL <- rownames(mf16)[ mf16$ColRep %in% badsampreps[3]]

badTB <- match(badrepsTB, rownames(nmds16$points))
badBA <- match(badrepsBA, rownames(nmds16$points))
badBL <- match(badrepsBL, rownames(nmds16$points))

outline[badTB] <- badsampColors[1]
outline[badBA] <- badsampColors[2]
outline[badBL] <- badsampColors[3]

bad <- match(badsamp, rownames(nmds16$points))
size[bad] <- badsampSize


pdf("FR16s_removed_samples.pdf",height=7,width=5)
par(fig=c(0,1,0,0.8), mar=c(4.1,4.1,2.1,2.1))
plot(nmds16$points, xlab="NMDS1", ylab="NMDS2"
     , pch=21
     , bg=colorSal[factor(salRound)]
     , col=outline
     , cex=size
)
par(fig=c(0,1,0.8,1), mar=c(0,0,0,0), new=TRUE)
plot(NULL,xlim=c(0,1), ylim=c(0,1), xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
legend("center", legend=c("E-TB-2016","E-BA-2016","E-BL-2016","Outliers removed"), pch=21, col=c(badsampColors,"grey"), pt.bg = c("white","white","white","grey"), pt.cex = c(1,1,1,2))
dev.off()

#### NMDS 18s ####
# highlight bad samples
bad <- match(badsamp, rownames(nmds18$points))
outline <- rep("black",length=nrow(nmds18$points))
outline[bad] <- badsampColors


pdf("FR18s_removed_samples.pdf",height=7,width=5)
par(fig=c(0,1,0,0.8), mar=c(4.1,4.1,2.1,2.1))
plot(nmds18$points, xlab="NMDS1", ylab="NMDS2"
     , pch=21
     , bg=colorSal[factor(salRound18s)]
     , col=outline
)
par(fig=c(0,1,0.8,1), mar=c(0,0,0,0), new=TRUE)
plot(NULL,xlim=c(0,1), ylim=c(0,1), xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
legend("center", legend=c("E-TB-2016","E-BA-2016","E-BL-2016"), pch=21, col=badsampColors, pt.bg = c("white","white","white"), pt.cex = c(1,1,1))
dev.off()


########## Questionable samples ppt 28?
# highlight bad samples
bad <- match(badsamp2, rownames(nmds18$points))
outline <- rep("black",length=nrow(nmds18$points))
outline[bad] <- badsampColors


pdf("FR18s_candidate_bad.pdf",height=7,width=5)
par(fig=c(0,1,0,0.8), mar=c(4.1,4.1,2.1,2.1))
plot(nmds18$points, xlab="NMDS1", ylab="NMDS2"
     , pch=21
     , bg=colorSal[factor(salRound18s)]
     , col=outline
)
dev.off()


bad <- match(badsamp2, rownames(nmds16$points))
outline <- rep("black",length=nrow(nmds16$points))
outline[bad] <- badsampColors

pdf("FR16s_candidate_bad.pdf",height=7,width=5)
par(fig=c(0,1,0,0.8), mar=c(4.1,4.1,2.1,2.1))
plot(nmds16$points, xlab="NMDS1", ylab="NMDS2"
     , pch=21
     , bg=colorSal[factor(salRound)]
     , col=outline
)
dev.off()


###### CHECKING DELETED SAMPLES #######



