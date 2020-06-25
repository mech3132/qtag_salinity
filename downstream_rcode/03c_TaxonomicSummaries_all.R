#!/bin/bash
library(RColorBrewer) # For palette colors
############################## Taxa summaries plots #############################
 
############ 16s COMPARATIVE ##################
output <- "./OUTPUT/03c_TaxaSummaries"
dir.create(output)

dirNames <- unlist(read.delim("./OUTPUT/dirNames.txt",header=FALSE))

biomFPWD <- paste0(dirNames[2],"/OTUTableText.txt")
biomBPWD <- paste0(dirNames[1],"/OTUTableText.txt")
metadataFPWD <- "./INPUT_DATA/Fraser_16S/MF_16sFraser_noConCOL.txt"
metadataBPWD <- "./INPUT_DATA/Baltic_16S/metadata_table.tsv"

biomF18PWD <- paste0(dirNames[3],"/OTUTableText.txt")
metadataF18PWD <- "./INPUT_DATA/Fraser_18S/MF_18sFraser_noConCOL.txt"

############## Pre-amble; loading files, formatting, relative abund ###################

taxaF <- read.delim(paste0(biomFPWD)
                   , strip.white = TRUE
                   , stringsAsFactors = FALSE
                   , header = TRUE
                   , row.names = 1
                   # , skip = 1
                   )
taxaB <- read.delim(paste0(biomBPWD)
                    , strip.white = TRUE
                    , stringsAsFactors = FALSE
                    , header = TRUE
                    , row.names = 1
                    # , skip = 1
                    )
colnames(taxaF) <- gsub(".","-", colnames(taxaF), fixed = TRUE)
colnames(taxaB) <- gsub(".","-", colnames(taxaB), fixed = TRUE)

# Load metadata
metadataF <- read.delim(paste0(metadataFPWD), stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""), row.names = 1)
metadataB <- read.delim(paste0(metadataBPWD), stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""), row.names = 1)
rownames(metadataF) <- gsub(".","-",rownames(metadataF), fixed = TRUE)
rownames(metadataB) <- gsub(".","-",rownames(metadataB), fixed = TRUE)

# Order by taxa
orderAlpha <- order(taxaF[,ncol(taxaF)])
taxaF <- taxaF[orderAlpha,]
orderAlpha <- order(taxaB[,ncol(taxaB)])
taxaB <- taxaB[orderAlpha,]


# Make taxonomy ref
taxonomyRefF <- data.frame(taxaF[,ncol(taxaF)])
rownames(taxonomyRefF) <- rownames(taxaF)
taxonomyRefB <- data.frame(taxaB[,ncol(taxaB)])
rownames(taxonomyRefB) <- rownames(taxaB)

#combine ref
alltaxa <- unique(rownames(taxonomyRefB), rownames(taxonomyRefF))
synthTaxaRef <- 1:length(alltaxa)
n <- 1
for ( i in alltaxa ) {
    balttaxa <- NA
    frastaxa <- NA
    frasID <- NA
    baltID <- NA
    
    baltID <- match(i,rownames(taxonomyRefB))
    frasID <- match(i, rownames(taxonomyRefF))
    if (!is.na(baltID)) {
        balttaxa <- as.character(taxonomyRefB[baltID,])
    } 
    if (!is.na(frasID)) {
        frastaxa <- as.character(taxonomyRefF[frasID,])
    } 
    if (is.na(balttaxa) & (!is.na(frastaxa))) {
        synthTaxaRef[n] <- frastaxa
    } else if (is.na(frastaxa) & (!is.na(balttaxa))) {
        synthTaxaRef[n] <- balttaxa
    } else if ( balttaxa == frastaxa ) {
        synthTaxaRef[n] <- balttaxa
    } else {
        print("ERROR: NOT SAME TAXA")
        print(frastaxa)
        print(balttaxa)
        synthTaxaRef[n] <- frastaxa
    }
    n <- n+1
}

# Make taxa summary separate
taxaF <- taxaF[,-(ncol(taxaF))]
taxaB <- taxaB[,-(ncol(taxaB))]

# Make taxa summary relative abundance
relAbund <- function(x) {
    total <- sum(x)
    newx <- x/total
    return(newx)
}

# Make get taxa names function
getTaxaNames <- function(x, delim) {
    tempname <- strsplit(as.character(x),split = paste0(delim))[[1]][c(3,6,7)]
    tempname <- gsub("^.*__","",tempname)
    return(paste0(tempname[1],": ",tempname[2],"_", tempname[3]))
}

# Apply function to taxa
taxaAbundF <- data.frame(apply(taxaF, 2, FUN = relAbund))
names(taxaAbundF) <- names(taxaF)
taxaAbundB <- data.frame(apply(taxaB, 2, FUN = relAbund))
names(taxaAbundB) <- names(taxaB)

# Metadata stuff

# Make sure only samples in metadata are ones in taxa summaries
metadataF <- metadataF[rownames(metadataF) %in% names(taxaAbundF),]
metadataB <- metadataB[rownames(metadataB) %in% names(taxaAbundB),]

# Sort taxa summaries sites and metadata so they are increasing in gradient
# Use re-order with gradientNames[4] to make in order of gradient
metadataF <- metadataF[order(as.numeric(metadataF[["SalinityEnviron"]])),]
metadataB <- metadataB[order(as.numeric(metadataB[["SalinityEnviron"]])),]

# then make taxa the same
taxaAbundF <- taxaAbundF[,rownames(metadataF)]
taxaAbundB <- taxaAbundB[,rownames(metadataB)]
                            , na.strings = c('','na','NA'))

# Then, extract gradient values for these
ordered.gradientF <- metadataF[["SalinityEnviron"]]
ordered.gradientB <- metadataB[["SalinityEnviron"]]

# Now, I want to replace the headers of the euktaxa.ordered with these values.
taxaAbundF.grad <- taxaAbundF
names(taxaAbundF.grad) <- ordered.gradientF
taxaAbundB.grad <- taxaAbundB
names(taxaAbundB.grad) <- ordered.gradientB

# Changing the dataframe into a matrix so that I can plug it directly into barplot. 
taxaAbundF.matrix <- as.matrix(taxaAbundF.grad, stringsAsFactors = FALSE, header = TRUE, rownames.force = TRUE)
taxaAbundB.matrix <- as.matrix(taxaAbundB.grad, stringsAsFactors = FALSE, header = TRUE, rownames.force = TRUE)

### Function for combining levels ####
combineGroupAbund <- function(g,allOTUnames,groups,taxaAbund.matrix) {
    tomerge <- allOTUnames[which(groups==g)]
    if (length(tomerge) > 1) {
        return(colSums(taxaAbund.matrix[tomerge,]))
    } else {
        return(taxaAbund.matrix[tomerge,])
    }
}

#### Get spacing ####
######## Making barplot; taxa summaries ##########
print("Making barplots")

#F first
# First, list all unique salinities
uniqueF.num <- unique(names(taxaAbundF.grad))

# Then, count how many of each unique gradient there is. 
# This loop counts how many reps there are of each 'unique' gradient in the taxa file
spacingF <- data.frame()
for (i in 1:length(uniqueF.num)) {
    spacingF[i,1] <- as.numeric(uniqueF.num[i])
    spacingF[i,2] <- sum(uniqueF.num[i] == names(taxaAbundF.grad))
}
names(spacingF) <- c("Sal","count")

# Calculating distances
spacingF.num <- vector()
for (i in 1:length(spacingF[,2])) {
    if (i == 1) {
        spacingF.num <- c(0, rep(0, spacingF[i,2]-1))
    } else {
        spacingF.num <- c(spacingF.num, (spacingF[i,1]-spacingF[i-1,1])*1 , rep(0, spacingF[i,2]-1))
    }
}




#B second
# First, list all unique salinities
uniqueB.num <- unique(names(taxaAbundB.grad))

# Then, count how many of each unique gradient there is. 
# This loop counts how many reps there are of each 'unique' gradient in the taxa file
spacingB <- data.frame()
for (i in 1:length(uniqueB.num)) {
    spacingB[i,1] <- as.numeric(uniqueB.num[i])
    spacingB[i,2] <- sum(uniqueB.num[i] == names(taxaAbundB.grad))
}
names(spacingB) <- c("Sal","count")

# Calculating distances
spacingB.num <- vector()
for (i in 1:length(spacingB[,2])) {
    if (i == 1) {
        spacingB.num <- c(0, rep(0, spacingB[i,2]-1))
    } else {
        spacingB.num <- c(spacingB.num, (spacingB[i,1]-spacingB[i-1,1])*1 , rep(0, spacingB[i,2]-1))
    }
}

levelNames <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
for ( level in 4:6) {
    # For F
    allOTUnamesF <- rownames(taxaAbundF.matrix)
    allTaxanamesF <- taxonomyRefF[allOTUnamesF,]
    groupsF <- c()
    for ( i in allTaxanamesF) {
        classtemp <- strsplit(as.character(i),split="; ")[[1]][3]
        grouptemp <- strsplit(as.character(i),split="; ")[[1]][level]
        
        if (length(grep("Clade ",grouptemp))>0) {
            grouptempadd <- strsplit(as.character(i),split="; ")[[1]][level-1]
            grouptemp <- gsub("^.*__","",grouptemp)
            grouptemp <- paste0(grouptempadd, " ",grouptemp)
            # print(grouptemp)
        }
        
        classtemp <- gsub("^.*__", "", classtemp)
        grouptemp <- gsub("^.*__","",grouptemp)
        groupsF <- c(groupsF,paste0(classtemp,": ",grouptemp))
    }
    uniqueGroupsF <- unique(groupsF)
    
    # For B
    allOTUnamesB <- rownames(taxaAbundB.matrix)
    allTaxanamesB <- taxonomyRefB[allOTUnamesB,]
    groupsB <- c()
    for ( i in allTaxanamesB) {
        classtemp <- strsplit(as.character(i),split="; ")[[1]][3]
        grouptemp <- strsplit(as.character(i),split="; ")[[1]][level]
        
        if (length(grep("Clade ",grouptemp))>0) {
            grouptempadd <- strsplit(as.character(i),split="; ")[[1]][level-1]
            grouptemp <- gsub("^.*__","",grouptemp)
            grouptemp <- paste0(grouptempadd, " ",grouptemp)
            # print(grouptemp)
        }

        classtemp <- gsub("^.*__", "", classtemp)
        grouptemp <- gsub("^.*__","",grouptemp)
        groupsB <- c(groupsB,paste0(classtemp,": ",grouptemp))
    }
    uniqueGroupsB <- unique(groupsB)
    # F mat
    newMatrixF <- matrix(ncol=ncol(taxaAbundF.matrix), nrow=length(unique(groupsF)), dimnames = list(uniqueGroupsF,colnames(taxaAbundF.matrix)))
    for ( g in uniqueGroupsF ) {
        newMatrixF[g,] <-combineGroupAbund(g,allOTUnamesF,groupsF,taxaAbundF.matrix)
    }
    # B mat
    newMatrixB <- matrix(ncol=ncol(taxaAbundB.matrix), nrow=length(unique(groupsB)), dimnames = list(uniqueGroupsB,colnames(taxaAbundB.matrix)))
    for ( g in uniqueGroupsB ) {
        newMatrixB[g,] <-combineGroupAbund(g,allOTUnamesB,groupsB,taxaAbundB.matrix)
    }
    
    #Now, find top 10 most abundant in each group, discounting shared ones.
    top10F.order <- order(rowSums(newMatrixF), decreasing=TRUE)
    top10B.order <- order(rowSums(newMatrixB), decreasing=TRUE)
    top10.taxa.all <- c()
    n <- 1
    while (length(top10.taxa.all) < 20) {
        tempFtaxa <- rownames(newMatrixF)[top10F.order[n]]
        tempBtaxa <- rownames(newMatrixB)[top10B.order[n]]
        top10.taxa.all <- c(top10.taxa.all, tempFtaxa, tempBtaxa)
        top10.taxa.all <- unique(top10.taxa.all)
        n <- n + 1
    }
    # if we happen to have 21:
    #get colors
    col.pal <- c("maroon",brewer.pal(n=8, name=c("Set2")),brewer.pal(n=12, name=c("Paired")))
    # compine with top10 taxa
    top10.col <- cbind(top10.taxa.all,"col.pal"=col.pal[1:length(top10.taxa.all)])
    # quartz()
    # barplot(rep(1,21), col=col.pal)
    
    # get vector of colors for each matrix
    colF <- rep("white", nrow(newMatrixF))
    posColF <- cbind(top10.col,"pos"=match(top10.taxa.all, rownames(newMatrixF)))
    for ( i in 1:nrow(posColF) ) {
        colF[as.numeric(posColF[i,"pos"])] <- posColF[i,"col.pal"]
    }
    
    colB <- rep("white", nrow(newMatrixB))
    posColB <- cbind(top10.col,"pos"=match(top10.taxa.all, rownames(newMatrixB)))
    for ( i in 1:nrow(posColB) ) {
        colB[as.numeric(posColB[i,"pos"])] <- posColB[i,"col.pal"]
    }

    ###### ADJ SPACING FIRST
    pdf(paste0("./",output,"/TaxaSummariesPartitioned_by",levelNames[level],".pdf"), width = 10, height = 5)
    # quartz(,10,5)
    par(mar=c(2,4.1,1,1), xpd=TRUE)
    
    ### BALTIC ###
    par(fig = c(0,0.7,0.5,1))
    barplot(newMatrixB
            , col = colB # Using the generated colours from above, I colour each bar that has >5% representation a different colour.
            , space = spacingB.num # This lines the bar plot up with the boxplot in spacing
            , cex.axis = 0.6 # scaling axis
            , cex.lab = 1 # scaling label
            , xaxt = 'n'
            , border = NA
            
    )
    # Add manual x axis: space out properly
    totalDist <- 0:length(spacingB.num)
    for ( i in 1:length(spacingB.num) ) {
        totalDist[i+1] <- totalDist[i] + spacingB.num[i] + 1
    }
    totalDist <- totalDist-0.5
    axis(side = 1
         , at = totalDist[-1]
         , labels = colnames(taxaAbundB.matrix)
         , cex.axis = 0.5
         , line=-1
         , tick=FALSE
    )
    title(xlab="Salinity",line=1,cex.lab=0.5)
    title(ylab="Relative Abundance", line=2, cex.lab=0.5)
    ### FRASER ####
    par(fig = c(0,0.7,0,0.5), new=TRUE)
    barplot(newMatrixF
            , col = colF # Using the generated colours from above, I colour each bar that has >5% representation a different colour.
            , space = spacingF.num # This lines the bar plot up with the boxplot in spacing
            , cex.axis = 0.6 # scaling axis
            , cex.lab = 1 # scaling label
            , xaxt = 'n'
            , border = NA
            
    )
    # Add manual x axis: space out properly
    totalDist <- 0:length(spacingF.num)
    for ( i in 1:length(spacingF.num) ) {
        totalDist[i+1] <- totalDist[i] + spacingF.num[i] + 1
    }
    totalDist <- totalDist-0.5
    axis(side = 1
         , at = totalDist[-1]
         , labels = colnames(taxaAbundF.matrix)
         , cex.axis = 0.5
         , line=-1
         , tick=FALSE
    )
    title(xlab="Salinity",line=1,cex.lab=0.5)
    title(ylab="Relative Abundance", line=2, cex.lab=0.5)
    ## legend
    par(fig = c(0.7,1,0,1), mar=c(1,1,1,1),new=TRUE)
    plot(0,0,xlab="",ylab="",xaxt="n",yaxt="n",bty="n",pch="")
    # re-order legend so alphabetical
    alphaOrder <- c(length(top10.col[,1])+1,order(top10.col[,1], decreasing = TRUE))
    legend("left"
           , legend = c(top10.col[,1], "Other (<5%)")[alphaOrder]
           , pch = 22 # squares
           , col= "black" # black outline
           , pt.bg = c(as.character(top10.col[,2]), "white")[alphaOrder] #fill squares with colours
           , bty = "n" # Get rid of outer box
           , title = paste0("Class: ",levelNames[level])
           , cex = 0.7
           , pt.cex = 1.5
           , xpd = TRUE
    )
    dev.off()
    
}

############ ~~~~~~18s ALONE~~~~~~ ##################
############## Pre-amble; loading files, formatting, relative abund ###################

taxaF <- read.delim(paste0(biomF18PWD)
                    , strip.white = TRUE
                    , stringsAsFactors = FALSE
                    , header = TRUE
                    , row.names = 1
                    # , skip = 1
                    )

colnames(taxaF) <- gsub(".","-", colnames(taxaF), fixed = TRUE)

# Load metadata
metadataF <- read.delim(paste0(metadataF18PWD), stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""), row.names = 1)
rownames(metadataF) <- gsub(".","-",rownames(metadataF), fixed = TRUE)

# Order by taxa
orderAlpha <- order(taxaF[,ncol(taxaF)])
taxaF <- taxaF[orderAlpha,]


# Make taxonomy ref
taxonomyRefF <- data.frame(taxaF[,ncol(taxaF)])
rownames(taxonomyRefF) <- rownames(taxaF)

# Make taxa summary separate
taxaF <- taxaF[,-(ncol(taxaF))]

# Apply function to taxa
taxaAbundF <- data.frame(apply(taxaF, 2, FUN = relAbund))
names(taxaAbundF) <- names(taxaF)


# # Order by abundance
# orderAbundF <- order(rowSums(taxaAbundF), decreasing = TRUE)
# taxaAbundF <- taxaAbundF[orderAbundF,]
# orderAbundB <- order(rowSums(taxaAbundB), decreasing = TRUE)
# taxaAbundB <- taxaAbundB[orderAbundB,]

#--
# Metadata stuff

# Make sure only samples in metadata are ones in taxa summaries
metadataF <- metadataF[rownames(metadataF) %in% names(taxaAbundF),]

# Sort taxa summaries sites and metadata so they are increasing in gradient
# Use re-order with gradientNames[4] to make in order of gradient
metadataF <- metadataF[order(as.numeric(metadataF[["SalinityEnviron"]])),]

# then make taxa the same
taxaAbundF <- taxaAbundF[,rownames(metadataF)]

# Then, extract gradient values for these
ordered.gradientF <- metadataF[["SalinityEnviron"]]

# Now, I want to replace the headers of the euktaxa.ordered with these values.
taxaAbundF.grad <- taxaAbundF
names(taxaAbundF.grad) <- ordered.gradientF

# Changing the dataframe into a matrix so that I can plug it directly into barplot. 
taxaAbundF.matrix <- as.matrix(taxaAbundF.grad, stringsAsFactors = FALSE, header = TRUE, rownames.force = TRUE)

#### Get spacing ####
######## Making barplot; taxa summaries ##########
print("Making barplots")

#F first
# First, list all unique salinities
uniqueF.num <- unique(names(taxaAbundF.grad))

# Then, count how many of each unique gradient there is. 
# This loop counts how many reps there are of each 'unique' gradient in the taxa file
spacingF <- data.frame()
for (i in 1:length(uniqueF.num)) {
    spacingF[i,1] <- as.numeric(uniqueF.num[i])
    spacingF[i,2] <- sum(uniqueF.num[i] == names(taxaAbundF.grad))
}
names(spacingF) <- c("Sal","count")

# Calculating distances
spacingF.num <- vector()
for (i in 1:length(spacingF[,2])) {
    if (i == 1) {
        spacingF.num <- c(0, rep(0, spacingF[i,2]-1))
    } else {
        spacingF.num <- c(spacingF.num, (spacingF[i,1]-spacingF[i-1,1])*1 , rep(0, spacingF[i,2]-1))
    }
}


levelNames <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
for ( level in 4:6) {
    # For F
    allOTUnamesF <- rownames(taxaAbundF.matrix)
    allTaxanamesF <- taxonomyRefF[allOTUnamesF,]
    groupsF <- c()
    for ( i in allTaxanamesF) {
        classtemp <- strsplit(as.character(i),split="; ")[[1]][3]
        if (is.na(classtemp)) {
            classtemp <- "Unassigned"
        }
        grouptemp <- strsplit(as.character(i),split="; ")[[1]][level]
        if (is.na(grouptemp)) {
            grouptemp <- "Unassigned"
        }
        # 
        # if (length(grep("Clade ",grouptemp))>0) {
        #     grouptempadd <- strsplit(as.character(i),split="; ")[[1]][level-1]
        #     grouptemp <- gsub("^.*__","",grouptemp)
        #     grouptemp <- paste0(grouptempadd, " ",grouptemp)
        #     # print(grouptemp)
        # }
        # 
        classtemp <- gsub("^.*__", "", classtemp)
        grouptemp <- gsub("^.*__","",grouptemp)
        groupsF <- c(groupsF,paste0(classtemp,": ",grouptemp))
    }
    uniqueGroupsF <- unique(groupsF)
    
    # F mat
    newMatrixF <- matrix(ncol=ncol(taxaAbundF.matrix), nrow=length(unique(groupsF)), dimnames = list(uniqueGroupsF,colnames(taxaAbundF.matrix)))
    for ( g in uniqueGroupsF ) {
        newMatrixF[g,] <-combineGroupAbund(g,allOTUnamesF,groupsF,taxaAbundF.matrix)
    }

    
    #Now, find top 10 most abundant in each group, discounting shared ones.
    top10F.order <- order(rowSums(newMatrixF), decreasing=TRUE)
    top10.taxa.all <- c()
    n <- 1
    while (length(top10.taxa.all) < 20) {
        tempFtaxa <- rownames(newMatrixF)[top10F.order[n]]
        top10.taxa.all <- c(top10.taxa.all, tempFtaxa)
        n <- n + 1
    }
    # if we happen to have 21:
    #get colors
    col.pal <- c(brewer.pal(n=8, name=c("Set2")),brewer.pal(n=12, name=c("Paired")))
    # compine with top10 taxa
    top10.col <- cbind(top10.taxa.all,"col.pal"=col.pal[1:length(top10.taxa.all)])
    # quartz()
    # barplot(rep(1,21), col=col.pal)
    
    # get vector of colors for each matrix
    colF <- rep("white", nrow(newMatrixF))
    posColF <- cbind(top10.col,"pos"=match(top10.taxa.all, rownames(newMatrixF)))
    for ( i in 1:nrow(posColF) ) {
        colF[as.numeric(posColF[i,"pos"])] <- posColF[i,"col.pal"]
    }
    
    
    ###### ADJ SPACING FIRST
    pdf(paste0("./",output,"/18sTaxaSummariesPartitioned_by",levelNames[level],".pdf"), width = 10, height = 4)
    # quartz(,10,5)
    par(mar=c(2,4.1,2,1), fig = c(0,0.7,0.1,0.9), xpd=TRUE)
    
    ### FRASER ####
    barplot(newMatrixF
            , col = colF # Using the generated colours from above, I colour each bar that has >5% representation a different colour.
            , space = spacingF.num # This lines the bar plot up with the boxplot in spacing
            , cex.axis = 0.6 # scaling axis
            , cex.lab = 1 # scaling label
            , xaxt = 'n'
            , border = NA
            
    )
    # Add manual x axis: space out properly
    totalDist <- 0:length(spacingF.num)
    for ( i in 1:length(spacingF.num) ) {
        totalDist[i+1] <- totalDist[i] + spacingF.num[i] + 1
    }
    totalDist <- totalDist-0.5
    axis(side = 1
         , at = totalDist[-1]
         , labels = colnames(taxaAbundF.matrix)
         , cex.axis = 0.5
         , line=-1
         , tick=FALSE
    )
    title(xlab="Salinity",line=1,cex.lab=0.5)
    title(ylab="Relative Abundance", line=2, cex.lab=0.5)
    ## legend
    par(fig = c(0.7,1,0,1), mar=c(1,1,1,1),new=TRUE)
    alphaOrder <- c(length(top10.col[,1])+1,order(top10.col[,1], decreasing = TRUE))
    plot(0,0,xlab="",ylab="",xaxt="n",yaxt="n",bty="n",pch="")
    legend("left"
           , legend = c(top10.col[,1], "Other (<5%)")[alphaOrder]
           , pch = 22 # squares
           , col= "black" # black outline
           , pt.bg = c(as.character(top10.col[,2]), "white")[alphaOrder] #fill squares with colours
           , bty = "n" # Get rid of outer box
           , title = paste0("Class: ",levelNames[level])
           , cex = 0.7
           , pt.cex = 1.5
           , xpd = TRUE
    )
    dev.off()
    
    
    ## separate plots by season
    onlySOG <- grep("SOG",rownames(metadataF))
    alldates <- gsub("E-SOG-[0-9]","",rownames(metadataF)[onlySOG])
    uniqueDates <- unique(alldates)
    sepSamplesSOG <- list()
    for ( d in uniqueDates) {
        sepSamplesSOG[[paste0(d)]] <- rownames(metadataF[onlySOG[d == alldates],])
    }
    # Now, plot each one individually
    
    pdf(file=paste0("./",output,"/SimilarComposition_byDate_at",levelNames[level],".pdf"), 10, 5)
    par(mfrow=c(1,3))
    for (d in sort(names(sepSamplesSOG))) {
        colToKeep <- match(sepSamplesSOG[[d]], colnames(taxaAbundF))
        barplot(newMatrixF[,colToKeep]
                , col = colF # Using the generated colours from above, I colour each bar that has >5% representation a different colour.
                , cex.axis = 0.6 # scaling axis
                , cex.lab = 0.5 # scaling label
                , border = NA
                , names.arg = as.character(colnames(newMatrixF[,colToKeep]))
                , ylab="Relative Abundance"
                , main=paste0("Date: ", d)
        )
    }
    dev.off()
    
}

