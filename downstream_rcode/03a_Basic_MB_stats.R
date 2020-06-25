#!/bin/bash Rscript
library(tidyverse)
##### Get some basic stats from the modelBoundaries files #####
dirNames = unlist(read.delim("./OUTPUT/dirNames.txt", header=FALSE))

mbB16FP=paste0(dirNames[1],"/modelBoundaries_type.txt")
mbF16FP=paste0(dirNames[2],"/modelBoundaries_type.txt")
mbF18FP=paste0(dirNames[3],"/modelBoundaries_type.txt")

taxaRefB16FP = paste0(dirNames[1],"/taxaIDLegend.txt")
taxaRefF16FP = paste0(dirNames[2],"/taxaIDLegend.txt")
taxaRefF18FP = paste0(dirNames[3],"/taxaIDLegend.txt")

mbB16 = read.delim(mbB16FP, header=TRUE, row.names=1, stringsAsFactors = FALSE)
mbF16 = read.delim(mbF16FP, header=TRUE, row.names=1, stringsAsFactors = FALSE)
mbF18 = read.delim(mbF18FP, header=TRUE, row.names=1, stringsAsFactors = FALSE)

taxaRefB16 = read.delim(taxaRefB16FP,header=FALSE, row.names=1)
taxaRefF16 = read.delim(taxaRefF16FP,header=FALSE, row.names=1)
taxaRefF18 = read.delim(taxaRefF18FP,header=FALSE, row.names=1)

newDirName <- "./OUTPUT/03a_SummaryStats"
dir.create(newDirName)

## Updated after adjusting boundary calculations in original script 2may2020
# # get rid of noclass
mbB16.filt = mbB16[mbB16$typeSimple!='noclass',]
mbF16.filt = mbF16[mbF16$typeSimple!='noclass',]
mbF18.filt = mbF18[mbF18$typeSimple!='noclass',]

B16.ranges <- mbB16$boundariestwo-mbB16$boundaries
F16.ranges <- mbF16$boundariestwo-mbF16$boundaries
F18.ranges <- mbF18$boundariestwo-mbF18$boundaries

maxrange = max(c(B16.ranges, F16.ranges, F18.ranges))

sink(file = paste0(newDirName,"/basic_mb_stats_summary.txt"))
#baltic range limits
print("========16sBaltic========")
print("range limits:")
print(c(min(B16.ranges),max(B16.ranges)))
print("median: ")
print(median(B16.ranges))
#fraser range limits
print("========16sFraser========")
print("range limits:")
print(c(min(F16.ranges),max(F16.ranges)))
print("median: ")
print(median(F16.ranges))
print("========18sFraser========")
print("range limits:")
#fraser18 range limits
print(c(min(F18.ranges), max(F18.ranges)))
print("median: ")
print(median(F18.ranges))
sink()


pdf(paste0(newDirName,"/range_sizes.pdf"),5,10)
# quartz(,5,10)
par(mfrow=c(3,1))
hist(B16.ranges, xlim=c(0,maxrange), main="(A) Baltic Sea (16S) salinity tolerance range", xlab="Salinity Tolerance Range")
abline(v=median(B16.ranges), col="red", lty=3, sub=paste0("Median: ",median(B16.ranges)))
hist(F16.ranges, xlim=c(0,maxrange), main="(B) Fraser River Estuary (16S) salinity tolerance range", xlab="Salinity Tolerance Range")
abline(v=median(F16.ranges), col="red", lty=3, sub=paste0("Median: ",median(F16.ranges)))
hist(F18.ranges, xlim=c(0,maxrange), main="(C) Fraser River Estuary (18S) salinity tolerance range", xlab="Salinity Tolerance Range")
abline(v=median(F18.ranges), col="red", lty=3, sub=paste0("Median: ",median(F18.ranges)))
dev.off()

sink(paste0(newDirName,"/basic_mb_stats_summary.txt"), append = TRUE)
wilcox.test(B16.ranges, F16.ranges)
wilcox.test(F18.ranges, F16.ranges)
wilcox.test(B16.ranges, F18.ranges)

#count number in each
print("========16sBaltic========")
print("brackish proportion")
sum(mbB16$typeSimple=="brackishRestricted")/nrow(mbB16) # brackish Baltic
print("fresh proportion")
sum(mbB16$typeSimple=="freshRestricted")/nrow(mbB16) # fresh Baltic
print("marine proportion")
sum(mbB16$typeSimple=="marineRestricted")/nrow(mbB16) # marine Baltic
print("noclass proportion")
sum(mbB16$typeSimple=="noclass")/nrow(mbB16) # marine Baltic
print("total_reads")
nrow(mbB16)

print("========16sFraser========")
print("brackish proportion")
sum(mbF16$typeSimple=="brackishRestricted")/nrow(mbF16) # brackish Fraser
print("fresh proportion")
sum(mbF16$typeSimple=="freshRestricted")/nrow(mbF16) # fresh Fraser
print("marine proportion")
sum(mbF16$typeSimple=="marineRestricted")/nrow(mbF16) # marine Fraser
print("noclass proportion")
sum(mbF16$typeSimple=="noclass")/nrow(mbF16) # noclass Fraser
print("total_reads")
nrow(mbF16)

print("========18sFraser========")
print("brackish proportion")
sum(mbF18$typeSimple=="brackishRestricted")/nrow(mbF18) # brackish Fraser
print("fresh proportion")
sum(mbF18$typeSimple=="freshRestricted")/nrow(mbF18) # fresh Fraser
print("marine proportion")
sum(mbF18$typeSimple=="marineRestricted")/nrow(mbF18) # marine Fraser
print("noclass proportion")
sum(mbF18$typeSimple=="noclass")/nrow(mbF18) # noclass Fraser
print("total_reads")
nrow(mbF18)
sink()

## Finding mid-point in tolerance range
B16.mid = apply(mbB16[,c("boundaries","boundariestwo")], MARGIN = 1, FUN =mean)
F16.mid = apply(mbF16[,c("boundaries","boundariestwo")], MARGIN = 1, FUN =mean)
F18.mid = apply(mbF18[,c("boundaries","boundariestwo")], MARGIN = 1, FUN =mean)

maxmid = max(c(F16.mid,B16.mid,F18.mid))+2

pdf(paste0(newDirName,"/midpoints.pdf"),5,10)
# quartz(,5,10)
par(mfrow=c(3,1))
hist(B16.mid, xlim=c(0,maxmid), main="Baltic 16s range midpoint", xlab="Midpoint of tolerance range", breaks = 20)
hist(F16.mid, xlim=c(0,maxmid), main="Fraser 16s range midpoint", xlab="Midpoint of tolerance range", breaks = 20)
hist(F18.mid, xlim=c(0,maxmid), main="Fraser 18s range midpoint", xlab="Midpoint of tolerance range", breaks = 20)
dev.off()


##### Comparing shared OTUs ######

####### Shared otus by classification (complex) #####
commonOTUs <- rownames(mbB16[na.omit(match(rownames(mbF16), rownames(mbB16))),])

combinedmb <- matrix(ncol=6, nrow=length(commonOTUs),dimnames=list(commonOTUs,c("btype","ftype","bb1","bb2","fb1","fb2")))
for ( otu in commonOTUs ) {
    ftemp <- as.vector(mbF16[otu,c('type','boundaries','boundariestwo')])
    btemp <- as.vector(mbB16[otu,c('type','boundaries','boundariestwo')])
    combinedmb[otu,] <- c(as.character(btemp[1]),as.character(ftemp[1]),as.numeric(btemp[c(2,3)]), as.numeric(ftemp[c(2,3)]))
}

sum(combinedmb[,c("btype")] == combinedmb[,c("ftype")])/nrow(combinedmb)

categories <- c("freshBloom" # 1
                ,"freshRestricted" #1
                ,"freshPeak" #2
                ,"brackishPeakLoToler" #4
                ,"brackishBloom" #5
                ,"brackishRestricted" #5
                ,"brackishPeakHiToler" #6
                ,"marinePeak" #8
                ,"marineRestricted" #9
                ,"marineBloom") #9
posCat <- c(1,1,2,4,5,5,6,8,9,9)

# lcat <- length(categories)
lcat <- max(posCat)

fpos <- sapply(combinedmb[,c("ftype")],function(x) {
    posCat[match(x,categories)]
})
bpos <- sapply(combinedmb[,c("btype")],function(x) {
    posCat[match(x,categories)]
})

combpos <- as.data.frame(cbind(bpos,fpos, comb=paste0(bpos,"-",fpos)))
# combpos[which(combpos$comb == '9-4'),]
counts_combpos <- table(combpos[,'comb'])
maxCount <- max(counts_combpos)
minCount <- min(counts_combpos)
dividend <- (maxCount)/40
thickness_combpos <- counts_combpos/dividend

pdf(paste0(newDirName,"/shared_otus_lineplot_complex.pdf"))
# quartz()
par(mar=c(4.2,10.2,2.2,10.2))
plot(NULL, xlim=c(0,1), ylim=c(0,lcat+2), xaxt="n", yaxt="n", xlab="",ylab="", bty="n")
for ( i in unique(posCat) ) {
    for (j in unique(posCat) ) {
        combpos_temp <- paste0(i,"-",j)
        thick_temp <- thickness_combpos[combpos_temp]
        if ( !is.na(thick_temp)) {
            # print(paste0(i,"-",j))
            lines(x=c(0,1), y=c(i,j), col=rgb(0.5,0.5,0.5,0.3), lwd=thick_temp)
        }
    }
}
axis(1, at = c(0,1), labels = c("Baltic","Fraser"), las=1, line=-1, cex.axis=2, tick = FALSE)
axis(side = 2, at=unique(posCat), labels=c("Fresh","Fresh/ marinetolerant","Brack/ freshtolerant","Brackish","Brack / marinetolerant","Marine/ freshtolerant","Marine"), las=2)
axis(side = 4, at=unique(posCat), labels=c("Fresh","Fresh/ marinetolerant","Brack/ freshtolerant","Brackish","Brack / marinetolerant","Marine/ freshtolerant","Marine"), las=2)
legend("top", legend=c(paste0("       ",maxCount," OTUs"),"",paste0("       ",minCount," OTU")), lwd=c(max(thickness_combpos),NA,min(thickness_combpos)), bty="n", col=rgb(0.5,0.5,0.5,0.3))
dev.off()

####### Shared otus by classification (simple) #####

combinedmb.simple <- matrix(ncol=6, nrow=length(commonOTUs),dimnames=list(commonOTUs,c("btype","ftype","bb1","bb2","fb1","fb2")))
for ( otu in commonOTUs ) {
    ftemp <- as.vector(mbF16[otu,c('typeSimple','boundaries','boundariestwo')])
    btemp <- as.vector(mbB16[otu,c('typeSimple','boundaries','boundariestwo')])
    combinedmb.simple[otu,] <- c(as.character(btemp[1]),as.character(ftemp[1]),as.numeric(btemp[c(2,3)]), as.numeric(ftemp[c(2,3)]))
}

categories.simple <- c("freshRestricted","brackishRestricted","marineRestricted")
lcat <- length(categories.simple)

fpos <- sapply(combinedmb.simple[,c("ftype")],function(x) {
    match(x,categories.simple)
})
bpos <- sapply(combinedmb.simple[,c("btype")],function(x) {
    match(x,categories.simple)
})

combpos <- as.data.frame(cbind(bpos,fpos,comb=paste0(bpos, "-",fpos)))
counts_combpos <- table(combpos[,'comb'])
maxCount <- max(counts_combpos)
minCount <- min(counts_combpos)
dividend <- (maxCount)/40
thickness_combpos <- counts_combpos/dividend

pdf(paste0(newDirName,"/shared_otus_lineplot_simple.pdf"))
par(mar=c(4.2,10.2,2.2,10.2))
plot(NULL, xlim=c(0,1), ylim=c(0,lcat+1), xaxt="n", yaxt="n", xlab="",ylab="")
for ( i in 1:lcat ) {
    for (j in 1:lcat) {
        combpos_temp <- paste0(i,"-",j)
        thick_temp <- thickness_combpos[combpos_temp]
        if (!is.na(thick_temp)) {
            lines(x=c(0,1), y=c(i,j), col=rgb(0.5,0.5,0.5,0.3), lwd=thick_temp)
        }
    }
}
axis(1, at = c(0,1), labels = c("Baltic","Fraser"), las=2)
axis(side = 2, at=seq(1,lcat), labels=categories.simple, las=2)
axis(side = 4, at=seq(1,lcat), labels=categories.simple, las=2)
dev.off()

########## shared otus by peaking location #########
blim <- c(0,36)
flim <- c(0,34)

categories <- c("freshbloom","freshRestricted","freshPeak","brackishPeakLoToler","brackishBloom","BrackishPeakAllToler","brackishRestricted","brackishPeakHiToler","marinePeak","marineRestricted","marineBloom")
lcat <- length(categories)

peaksTable <- matrix(ncol=2,nrow=length(commonOTUs), dimnames=list(commonOTUs,c("b","f")))
for ( i in commonOTUs)  {
    #BALTIC
    if ( combinedmb.simple[i,'btype'] == 'freshRestricted' ) {
        peakpos <- (as.numeric(combinedmb.simple[i,'bb1']) - flim[1])/2
    } else if ( combinedmb.simple[i,'btype'] == 'brackishRestricted' ) {
        peakpos <- (as.numeric(combinedmb.simple[i,'bb2']) + as.numeric(combinedmb.simple[i,'bb1']))/2
    } else if ( combinedmb.simple[i,'btype'] == 'marineRestricted' ) {
        peakpos <- (flim[2] - as.numeric(combinedmb.simple[i,'bb1']))/2
    }
    peaksTable[i,"b"] <- peakpos
    
    #FRASER
    if ( combinedmb.simple[i,'ftype'] == 'freshRestricted' ) {
        peakpos <- (as.numeric(combinedmb.simple[i,'fb1']) - flim[1])/2
    } else if ( combinedmb.simple[i,'ftype'] == 'brackishRestricted' ) {
        peakpos <- (as.numeric(combinedmb.simple[i,'fb2']) + as.numeric(combinedmb.simple[i,'fb1']))/2
    } else if ( combinedmb.simple[i,'ftype'] == 'marineRestricted' ) {
        peakpos <- (flim[2] - as.numeric(combinedmb.simple[i,'fb1']))/2
    }
    peaksTable[i,"f"] <- peakpos
}
# 
# 
# combpos <- as.data.frame(cbind(fpos,bpos, comb=paste0(fpos,"-",bpos)))
# counts_combpos <- table(combpos[,'comb'])
# maxCount <- max(counts_combpos)
# minCount <- min(counts_combpos)
# dividend <- (maxCount)/40
# thickness_combpos <- counts_combpos/dividend

pdf(paste0(newDirName,"/shared_otus_lineplot_bysal.pdf"))
par(mar=c(4.2,4.2,2.2,4.2))
plot(NULL, xlim=c(0,1), ylim=c(0,36), xaxt="n", yaxt="n", xlab="",ylab="")
for ( r in 1:nrow(peaksTable)) {
    lines(x=c(0,1),y=peaksTable[r,], col=rgb(0.5,0.5,0.5,0.2))
}
axis(1, at = c(0,1), labels = c("Baltic","Fraser"), las=2)
axis(side = 2, at=seq(0,35,by=5), labels=seq(0,35,by=5), las=1)
axis(side = 4, at=seq(0,35,by=5), labels=seq(0,35,by=5), las=1)
mtext(side=2, text="Salinity", line=3)
mtext(side=4, text="Salinity", line=3)
dev.off()

####### COMPOSITION OF SPECIALIST TYPES #########

### CLASS LEVEL FOR BALTIC
B16.fresh.taxa = as.character(taxaRefB16[rownames(mbB16.filt[mbB16.filt$typeSimple == "freshRestricted",]),])
B16.marine.taxa = as.character(taxaRefB16[rownames(mbB16.filt[mbB16.filt$typeSimple == "marineRestricted",]),])
B16.brackish.taxa = as.character(taxaRefB16[rownames(mbB16.filt[mbB16.filt$typeSimple == "brackishRestricted",]),])

sink(paste0(newDirName,"/basic_mb_stats_summary.txt"), append=TRUE)
print("---------- BREAKDOWN BY CLASS-----------")
# Breakdown by class
print("========16sBaltic========")
print("Fresh")
# FRESH
B16.fresh.taxa.split = sapply(1:length(B16.fresh.taxa), function(x) {strsplit(B16.fresh.taxa[x], split="; ", fixed=TRUE)})
round(table(sapply(B16.fresh.taxa.split, FUN="[[", INDEX=3))/length(B16.fresh.taxa.split),2)
print("----")
round(table(sapply(B16.fresh.taxa.split, FUN="[[", INDEX=4))/length(B16.fresh.taxa.split),2)
print("----")
round(table(sapply(B16.fresh.taxa.split, FUN="[[", INDEX=5))/length(B16.fresh.taxa.split),2)
print("Marine")
# MARINE
B16.marine.taxa.split = sapply(1:length(B16.marine.taxa), function(x) {strsplit(B16.marine.taxa[x], split="; ", fixed=TRUE)})
round(table(sapply(B16.marine.taxa.split, FUN="[[", INDEX=3))/length(B16.marine.taxa.split),2)
print("----")
round(table(sapply(B16.marine.taxa.split, FUN="[[", INDEX=4))/length(B16.marine.taxa.split),2)
print("----")
round(table(sapply(B16.marine.taxa.split, FUN="[[", INDEX=5))/length(B16.marine.taxa.split),2)
print("Brackish")
# BRACKISH
B16.brackish.taxa.split = sapply(1:length(B16.brackish.taxa), function(x) {strsplit(B16.brackish.taxa[x], split="; ", fixed=TRUE)})
round(table(sapply(B16.brackish.taxa.split, FUN="[[", INDEX=3))/length(B16.brackish.taxa.split),2)
print("----")
round(table(sapply(B16.brackish.taxa.split, FUN="[[", INDEX=4))/length(B16.brackish.taxa.split),2)
print("----")
round(table(sapply(B16.brackish.taxa.split, FUN="[[", INDEX=5))/length(B16.brackish.taxa.split),2)



### CLASS LEVEL FOR FRASER
F16.fresh.taxa = as.character(taxaRefF16[rownames(mbF16.filt[mbF16.filt$typeSimple == "freshRestricted",]),])
F16.marine.taxa = as.character(taxaRefF16[rownames(mbF16.filt[mbF16.filt$typeSimple == "marineRestricted",]),])
F16.brackish.taxa = as.character(taxaRefF16[rownames(mbF16.filt[mbF16.filt$typeSimple == "brackishRestricted",]),])

print("========16sFraser========")
# Breakdown by class
print("fresh")
# FRESH
F16.fresh.taxa.split = sapply(1:length(F16.fresh.taxa), function(x) {strsplit(F16.fresh.taxa[x], split="; ", fixed=TRUE)})
round(table(sapply(F16.fresh.taxa.split, FUN="[[", INDEX=3))/length(F16.fresh.taxa.split),2)
print("----")
round(table(sapply(F16.fresh.taxa.split, FUN="[[", INDEX=4))/length(F16.fresh.taxa.split),2)
print("----")
round(table(sapply(F16.fresh.taxa.split, FUN="[[", INDEX=5))/length(F16.fresh.taxa.split),2)
print("marine")
# MARINE
F16.marine.taxa.split = sapply(1:length(F16.marine.taxa), function(x) {strsplit(F16.marine.taxa[x], split="; ", fixed=TRUE)})
round(table(sapply(F16.marine.taxa.split, FUN="[[", INDEX=3))/length(F16.marine.taxa.split),2)
print("----")
round(table(sapply(F16.marine.taxa.split, FUN="[[", INDEX=4))/length(F16.marine.taxa.split),2)
print("----")
round(table(sapply(F16.marine.taxa.split, FUN="[[", INDEX=5))/length(F16.marine.taxa.split),2)
print("brackish")
# BRACKISH
F16.brackish.taxa.split = sapply(1:length(F16.brackish.taxa), function(x) {strsplit(F16.brackish.taxa[x], split="; ", fixed=TRUE)})
round(table(sapply(F16.brackish.taxa.split, FUN="[[", INDEX=3))/length(F16.brackish.taxa.split),2)
print("----")
round(table(sapply(F16.brackish.taxa.split, FUN="[[", INDEX=4))/length(F16.brackish.taxa.split),2)
print("----")
round(table(sapply(F16.brackish.taxa.split, FUN="[[", INDEX=5))/length(F16.brackish.taxa.split),2)

### CLASS LEVEL FOR FRASER
F18.fresh.taxa = as.character(taxaRefF18[rownames(mbF18.filt[mbF18.filt$typeSimple == "freshRestricted",]),])
F18.marine.taxa = as.character(taxaRefF18[rownames(mbF18.filt[mbF18.filt$typeSimple == "marineRestricted",]),])
F18.brackish.taxa = as.character(taxaRefF18[rownames(mbF18.filt[mbF18.filt$typeSimple == "brackishRestricted",]),])

print("========18sFraser========")
# Breakdown by class
print("fresh")
# FRESH
F18.fresh.taxa.split = sapply(1:length(F18.fresh.taxa), function(x) {strsplit(F18.fresh.taxa[x], split="; ", fixed=TRUE)})
# add command because 18s has some "Unassigned"
F18.fresh.taxa.split <- lapply(F18.fresh.taxa.split, `length<-`, max(lengths(F18.fresh.taxa.split)))
round(table(sapply(F18.fresh.taxa.split, FUN="[[", INDEX=3))/length(F18.fresh.taxa.split),2)
print("----")
round(table(sapply(F18.fresh.taxa.split, FUN="[[", INDEX=4))/length(F18.fresh.taxa.split),2)
print("----")
round(table(sapply(F18.fresh.taxa.split, FUN="[[", INDEX=5))/length(F18.fresh.taxa.split),2)

print("marine")
# MARINE
F18.marine.taxa.split = sapply(1:length(F18.marine.taxa), function(x) {strsplit(F18.marine.taxa[x], split="; ", fixed=TRUE)})
F18.marine.taxa.split <- lapply(F18.marine.taxa.split, `length<-`, max(lengths(F18.marine.taxa.split)))
round(table(sapply(F18.marine.taxa.split, FUN="[[", INDEX=3))/length(F18.marine.taxa.split),2)
print("----")
round(table(sapply(F18.marine.taxa.split, FUN="[[", INDEX=4))/length(F18.marine.taxa.split),2)
print("----")
round(table(sapply(F18.marine.taxa.split, FUN="[[", INDEX=5))/length(F18.marine.taxa.split),2)

print("brackish")# BRACKISH
F18.brackish.taxa.split = sapply(1:length(F18.brackish.taxa), function(x) {strsplit(F18.brackish.taxa[x], split="; ", fixed=TRUE)})
F18.brackish.taxa.split <- lapply(F18.brackish.taxa.split, `length<-`, max(lengths(F18.brackish.taxa.split)))
round(table(sapply(F18.brackish.taxa.split, FUN="[[", INDEX=3))/length(F18.brackish.taxa.split),2)
print("----")
round(table(sapply(F18.brackish.taxa.split, FUN="[[", INDEX=4))/length(F18.brackish.taxa.split),2)
print("----")
round(table(sapply(F18.brackish.taxa.split, FUN="[[", INDEX=5))/length(F18.brackish.taxa.split),2)

print("18 Class merged")
fresh18.class <- data.frame(round(table(sapply(F18.fresh.taxa.split, FUN="[[", INDEX=3))/length(F18.fresh.taxa.split),2))
names(fresh18.class) <- c("Class","Fresh")
marine18.class <- data.frame(round(table(sapply(F18.marine.taxa.split, FUN="[[", INDEX=3))/length(F18.marine.taxa.split),2))
names(marine18.class) <- c("Class","Marine")
brack18.class <- data.frame(round(table(sapply(F18.brackish.taxa.split, FUN="[[", INDEX=3))/length(F18.brackish.taxa.split),2))
names(brack18.class) <- c("Class","Brackish")

fresh18.class %>% full_join(marine18.class) %>% full_join(brack18.class)

sink()



# Get reverse table: object to class; get distribution of fresh/brack/marine
taxaRefB16.split <- as.data.frame(taxaRefB16) %>%
    rename(Taxa=V2) %>% separate(Taxa, into=paste0("D_",seq(1,7)), sep="; ", remove=FALSE)
taxaRefF16.split <- as.data.frame(taxaRefF16) %>%
    rename(Taxa=V2) %>% separate(Taxa, into=paste0("D_",seq(1,7)), sep="; ", remove=FALSE)
taxaRefF18.split <- as.data.frame(taxaRefF18) %>%
    rename(Taxa=V2) %>% separate(Taxa, into=paste0("D_",seq(1,7)), sep="; ", remove=FALSE)
for (ds in c("B16","F16","F18")) {
    # ds = "F18"
    tempDS <- get(paste0("taxaRef",ds,".split"))
    tempMB <- get(paste0("mb",ds,".filt"))
    for ( lvl in c("D_3","D_4","D_5")) {
        # lvl="D_4"
        allClasses <- unique(tempDS[,lvl]) 
        allClasses.byType <- data.frame()
        total <- length(allClasses)
        pb <- txtProgressBar(min=0, max=total, style=3, )
        for ( r in 1:length(allClasses)) {
            # which(allClasses=="D_3__Rhodobacterales")
            # r = 4
            if (is.na(allClasses[r])) {
                next 
            }
            classOTUsTemp <- rownames(tempDS[which(tempDS[,lvl] == allClasses[r]),])
            if ( lvl == "D_4") {
                nameWrite = paste0(tempDS[rownames(tempDS)==classOTUsTemp[1],c("D_3")],tempDS[rownames(tempDS)==classOTUsTemp[1],c("D_4")])
            } else if (lvl=="D_5") {
                nameWrite =  paste0(tempDS[rownames(tempDS)==classOTUsTemp[1],c("D_3")],tempDS[rownames(tempDS)==classOTUsTemp[1],c("D_4")],tempDS[rownames(tempDS)==classOTUsTemp[1],c("D_5")])
            } else {
                nameWrite = paste0(tempDS[rownames(tempDS)==classOTUsTemp[1],c("D_3")])
                
            }
            
            tableTemp <- table(tempMB[rownames(tempMB) %in%  classOTUsTemp, "typeSimple"])/length(tempMB[rownames(tempMB) %in%  classOTUsTemp, "typeSimple"])
            allClasses.byType <- rbind(allClasses.byType, data.frame(Taxa=nameWrite, Fresh=as.numeric(tableTemp["freshRestricted"]), Brackish=as.numeric(tableTemp["brackishRestricted"]), Marine=as.numeric(tableTemp["marineRestricted"]),total=length(classOTUsTemp)) )
            Sys.sleep(0.1)
            setTxtProgressBar(pb, r)
        }
        close(pb)
        write.table(allClasses.byType, file = paste0(newDirName,"/",ds,"_",lvl,".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    }
}



##### FIND OTUS that are classified as opposites #######
sink(paste0(newDirName,"/basic_mb_stats_summary.txt"), append=TRUE)
print("Number of classifications that are EXACT same")
# Number of classifications that are EXACT same
sum(combinedmb[,'btype'] == combinedmb[,'ftype'])
print("total shared")
nrow(combinedmb)

print("Number of classifications that are EXACT samem noclass ommitted")
# Number of classifications that are EXACT same; noclass ommitted
combinedmb.filt <- combinedmb[-grep("noclass", combinedmb),]
sum(combinedmb.filt[,'btype'] == combinedmb.filt[,'ftype'])
print("total shared")
nrow(combinedmb.filt)
31/108

print("Number of classifications that same broad type")
# Number of classifications that are same broad type
sum(combinedmb.simple[,'btype'] == combinedmb.simple[,'ftype'])
print("total shared")
nrow(combinedmb.simple)
61/129

print("Number of classifications that same broad type; noclass ommitted")
# Number of classifications that are same broad type; noclass ommitted
combinedmb.simple.filt <- combinedmb.simple[-grep("noclass", combinedmb.simple),]
sum(combinedmb.simple.filt[,'btype'] == combinedmb.simple.filt[,'ftype'])
print("total shared")
nrow(combinedmb.simple.filt)
58/108

sink()
# find the one that is opposite
oppositeClass <- c()
oppositeClass <- c(oppositeClass, names(which((combinedmb[,'btype'] == 'brackishPeakLoToler') & (combinedmb[,'ftype'] == 'marineRestricted'))))

