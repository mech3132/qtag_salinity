#!/bin/bash
library(tidyverse)
library(nlme)
library(phytools)
library(ape)
library(picante)
library(ggtree)
library(grid)
library(splancs)

########### This script takes modelboundaries file and a tree to plot how many are brackish/marine/fresh at each level in each taxonomic group

####### LOOP THROUGH DATASETS #######
dirNames <- unlist(read.delim("./OUTPUT/dirNames.txt", header=FALSE))
newDir = "./OUTPUT/03d_TypesAcrossTaxonomy"
dir.create(newDir)
############### FUNCTIONS #####################
source("./downstream_rcode/00_sourcecode2.R")

dataset <- c("B16","F16","F18")
for ( d in dataset) {
    if (d == "B16") {
        print("STARTING BALTIC DATASET")
        treepwd <- './OUTPUT/tree_B16_filt.tre'
        #Baltic16
        mbpwd <- paste0(dirNames[1],'/modelBoundaries_type.txt')
        otupwd <- paste0(dirNames[1],'/OTUTableText.txt')
        uniqueotupwd <- paste0(dirNames[1],'/OTUs.txt')
        output <- paste0(newDir, "/Baltic16S")
    } else if (d == "F16") {
        print("STARTING FRASER 16 DATASET")
        treepwd <- './OUTPUT/tree_F16_filt.tre'
        #Fraser16
        mbpwd <- paste0(dirNames[2],'/modelBoundaries_type.txt')
        otupwd <- paste0(dirNames[2],'/OTUTableText.txt')
        uniqueotupwd <- paste0(dirNames[2],'/OTUs.txt')
        output <- paste0(newDir, "/Fraser16S")
    } else if (d == "F18") {
        print("STARTING FRASER 18 DATASET")
        treepwd <- './OUTPUT/tree_F18_filt.tre'
        #Fraser18
        mbpwd <- paste0(dirNames[3],'/modelBoundaries_type.txt')
        otupwd <- paste0(dirNames[3],'/OTUTableText.txt')
        uniqueotupwd <- paste0(dirNames[3],'/OTUs.txt')
        output <- paste0(newDir, "/Fraser18S")
    } 
    ######## START ########
    tree16 <- read.tree(file = paste0(treepwd))
    print("DONE READING TREE")
    mb <- read.delim(file=paste0(mbpwd))
    otu <- read.delim(file=paste0(otupwd), skip=1, row.names=1)
    uniqueotu <- read.delim(file=uniqueotupwd, header=FALSE)
    dir.create(output)
    ####################################
   # Get taxonomy
    taxa <- cbind(rownames(otu), as.character(otu[,ncol(otu)]))
    # filter taxa out of OTU
    otu <- otu[,-ncol(otu)]

    # Filter OTU table to get rid of all OTUs not in uniqueotu; should automatically filter out low abund OTU
    otu.filt <- otu[match(uniqueotu$V1, rownames(otu)),]
    # Get rid of NAs
    otu.filt <- otu.filt[-c(rownames(otu.filt)=="NA"),]
    
    # Prune tree so there are only things in the table.
    print("pruning...")
    tree.pruned <- prune.sample(t(otu.filt), tree16)
    print("done pruning")
    
    # if dataset is 18S, root tree [ not already rooted ]
    if (d == "F18") {
        taxaleg <- taxa[taxa[,1] %in% tree.pruned$tip.label,]
        
        taxa.split <- sapply(taxaleg[,2], FUN = function(x) {strsplit(as.character(x), split="; __", fixed = TRUE)})
        
        class_list2 <- vector(length=length(taxa.split))
        for ( i in 1:length(taxa.split) ) {
            if ( is.na(taxa.split[[i]][2]) ) {
                class_list2[i] <- taxa.split[[i]][1]
            } else {
                class_list2[i] <- taxa.split[[i]][2]
            }
        }
        
        names(class_list2) <- as.character(taxaleg[,1])
        # unique(class_list2)
        opistikonts <- names(class_list2)[grep(pattern="Opisthokonta", x=class_list2)]
        # There are two fungi that are outliers that we want to remove.Otherwise, the rooting gets messed up.
        opistikonts <- opistikonts[-match(c("87590", "87593"),opistikonts)]
        
        opist_node <- findMRCA(tree16, tips=opistikonts)
        all_opist <- getDescendants(tree.pruned, node = opist_node)
        all_opist_otus <- tree.pruned$tip.label[all_opist]
        all_opist_otus <- all_opist_otus[!is.na(all_opist_otus)]
        tree.pruned <-unroot(tree.pruned)
        tree.pruned <- root(tree.pruned, outgroup = all_opist_otus, resolve.root = TRUE)
        # tree.pruned <- reroot(tree.pruned, node = opist_node)
        
    }
    
    # get order of taxa that GGtree will plot
    P <- ggtree::ggtree(tree.pruned
                ,branch.length="none"
                ,layout='rectangular')
    # getRoot(tree.pruned)
    toOrderGGtree <- P$data[which(P$data$isTip == TRUE),c("y","label")]
   
    tips <- rev(toOrderGGtree[order(toOrderGGtree$y),"label"])
    
    # get full taxa name
    namesFull <- lapply(1:length(tips$label), function(x) {
        fullname <- taxa[match(tips$label[x],taxa[,1]),2]
        return(fullname)
    })
    names(namesFull) <- tips$label
    
    # For 18S: since names are not evenly 7 levels even though I used the 7 level version?
    if ( d == "F18" ) {
        for ( i in 1:length(namesFull) ) {
            currentName <- unlist(strsplit(namesFull[[i]], "; "))
            newName <- NULL
            while ( length(currentName) < 7 ) {
                newName <- paste(c(currentName, currentName[length(currentName)]), sep = "; " , collapse = "; ") 
                currentName <- unlist(strsplit(newName, "; "))
            }
            if ( !is.null(newName) ) {
                namesFull[[i]] <- newName
            }
        }
    }

    # Get composition of each
    taxaNames <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    if ( d == "F18") {
        taxaNames <- c("D_0","D_1","D_2","D_3","D_4","D_5","D_6")
    }
    taxaLevels <- c(3,4,5,6) # Don't do species bc it doesn't add a lot of information
    alltaxalevelcolors <- matrix(ncol=length(taxaLevels),nrow=length(tips$label),dimnames = list(tips$label,taxaNames[taxaLevels]))
    for (l in taxaLevels) {
        # l = 3
        listClassComp <- propType(level=l,namesFull,mb)
        tableTypes <- tableAllTypes(level=l,namesFull)
        colors <- getColor(listClassComp)
        alltaxalevelcolors[,taxaNames[l]] <- as.vector(makebarcolor(colors,tableTypes,taxaNames[l]))
        if (l == 3) {
            listClassCompClass <- listClassComp
            tableTypesClass <- tableTypes
        } else if (l == 4) {
            listClassCompOrder <- listClassComp
            tableTypesOrder <- tableTypes
        } else if (l == 5 ) {
            listClassCompFam <- listClassComp
            tableTypesFam <- tableTypes
        } else if (l == 6) {
            listClassCompGen <- listClassComp
            tableTypesGen <- tableTypes
        }
    }
    # Get class sizes so I can make big-ass blocks
    keepname <- (lapply(listClassCompClass, FUN=sum) != 0)
    ClassNames <- gsub("_","",names(listClassCompClass)[keepname])
    ttc.sorted <- tableTypesClass[match(ClassNames, gsub("_","",names(tableTypesClass)))]
    blocksizeNames <- cbind(ClassNames,ttc.sorted)

    #Also, get composition of OTU
    alltypes <- mb[match(tips$label, mb[,"taxa"]),"typeSimple"]
    tempTab <- lapply(alltypes,function(x) {table(x,exclude="noclass")})
    OTU <- getColor(tempTab)
    
    #get class level to label
    listClassComp <- propType(level=3,namesFull,mb)
    tableTypes <- tableAllTypes(level=3,namesFull)
    colors <- getColor(listClassComp)
    
    # make empty matrix of 1's
    emptyBarMat <- matrix(ncol=length(taxaLevels)+1,nrow=length(tips$label),dimnames = list(tips$label,c(taxaNames[taxaLevels],"OTU")), data = 1)
    # add OTU colors
    alltaxalevelcolors.final <- cbind(alltaxalevelcolors,"OTU" = OTU)
    
    #####################
    # legend for class
    # Get class of each; then plot each color and all things in that block
    ##### AT CLASS LEVEL LEGEND #####
    alltaxalevelpercol <- list()
    count <- 0
    prevCol <- ""
    for ( i in 1:nrow(alltaxalevelcolors)) {
        currentCol <- alltaxalevelcolors[i,1]
        if (currentCol != prevCol) {
            count  <- count + 1
            alltaxalevelpercol[[count]] <- vector()
        }
        currentotu <- rownames(alltaxalevelcolors)[i]
        currentName <- namesFull[[currentotu]]
        currentClass <- gsub("^.*__","",strsplit(currentName,"; ")[[1]][3])
        alltaxalevelpercol[[count]] <- c(alltaxalevelpercol[[count]],currentClass)
        prevCol <- currentCol
    }
    x <- alltaxalevelcolors[,1]
    allCol <-x[x!=c(x[-1], FALSE)]
    
    #get names of groups, esp if multiple
    groupNames <- vector(length=length(allCol))
    groupAbund <- vector(length=length(allCol))
    for ( i in 1:length(alltaxalevelpercol)) {
        namesTemp <- unique(alltaxalevelpercol[[i]])
        temp <- ""
        for ( j in namesTemp) {
            temp <- paste(temp,j,sep=",")
        }
        groupNames[i] <- temp
        groupAbund[i] <- length(alltaxalevelpercol[[i]])
    }
    
    nameAndAbund <- paste0(groupNames,"(",groupAbund,")")
    
    # ## OLD; for some reason colours are backward
    # pdf(file="LegendForClass.pdf", height=length(nameAndAbund)/5)
    # par(mar=c(4.1,20,2,1))
    # barplot(rep(1,length(allCol)), col=rev(allCol), names.arg = rev(nameAndAbund), horiz = TRUE, las=2, cex.names=0.5, width=rev(groupAbund), border = rev(allCol))
    # dev.off()
    pdf(file=paste0(output,"/LegendForClass.pdf"), height=length(nameAndAbund)/5)
    par(mar=c(4.1,20,2,1))
    barplot(rep(1,length(allCol)), col=(allCol), names.arg = (nameAndAbund), horiz = TRUE, las=2, cex.names=0.5, width=(groupAbund), border = (allCol))
    dev.off()
    
    # Now, make a grey barplot with data above
    blankBar <- rev(sapply(alltaxalevelpercol, function(x) length(x)))
    
    ##### AT ORDER LEVEL LEGEND #####
    alltaxalevelpercol.order <- list()
    count <- 0
    prevCol <- ""
    for ( i in 1:nrow(alltaxalevelcolors)) {
        currentCol <- alltaxalevelcolors[i,2]
        if (currentCol != prevCol) {
            count  <- count + 1
            alltaxalevelpercol.order[[count]] <- vector()
        }
        currentotu <- rownames(alltaxalevelcolors)[i]
        currentName <- namesFull[[currentotu]]
        currentClass <- gsub("^.*__","",strsplit(currentName,"; ")[[1]][4])
        alltaxalevelpercol.order[[count]] <- c(alltaxalevelpercol.order[[count]],currentClass)
        prevCol <- currentCol
    }
    x <- alltaxalevelcolors[,2]
    allCol.order <-x[x!=c(x[-1], FALSE)]
    
    #get names of groups, esp if multiple
    groupNames.order <- vector(length=length(allCol.order))
    groupAbund <- vector(length=length(allCol.order))
    for ( i in 1:length(alltaxalevelpercol.order)) {
        namesTemp <- unique(alltaxalevelpercol.order[[i]])
        temp <- ""
        for ( j in namesTemp) {
            temp <- paste(temp,j,sep=",")
        }
        groupNames.order[i] <- temp
        groupAbund[i] <- length(alltaxalevelpercol.order[[i]])
    }
    
    nameAndAbund <- paste0(groupNames.order,"(",groupAbund,")")
    
    # pdf(file="LegendForOrder.pdf", height=length(nameAndAbund)/10)
    # par(mar=c(4.1,20,2,1))
    # barplot(rep(1,length(allCol.order)), col=rev(allCol.order), names.arg = rev(nameAndAbund), horiz = TRUE, las=2, cex.names=0.5, width=rev(groupAbund), border = rev(allCol.order))
    # dev.off()
    pdf(file=paste0(output,"/LegendForOrder.pdf"), height=length(nameAndAbund)/10)
    par(mar=c(4.1,20,2,1))
    barplot(rep(1,length(allCol.order)), col=(allCol.order), names.arg = (nameAndAbund), horiz = TRUE, las=2, cex.names=0.5, width=(groupAbund), border = (allCol.order))
    dev.off()

    
    ##### AT FAM LEVEL LEGEND #####
    alltaxalevelpercol.fam <- list()
    count <- 0
    prevCol <- ""
    for ( i in 1:nrow(alltaxalevelcolors)) {
        currentCol <- alltaxalevelcolors[i,3]
        if (currentCol != prevCol) {
            count  <- count + 1
            alltaxalevelpercol.fam[[count]] <- vector()
        }
        currentotu <- rownames(alltaxalevelcolors)[i]
        currentName <- namesFull[[currentotu]]
        currentClass <- gsub("^.*__","",strsplit(currentName,"; ")[[1]][5])
        alltaxalevelpercol.fam[[count]] <- c(alltaxalevelpercol.fam[[count]],currentClass)
        prevCol <- currentCol
    }
    x <- alltaxalevelcolors[,3]
    allCol.fam <-x[x!=c(x[-1], FALSE)]
    
    #get names of groups, esp if multiple
    groupNames.fam <- vector(length=length(allCol.fam))
    groupAbund <- vector(length=length(allCol.fam))
    for ( i in 1:length(alltaxalevelpercol.fam)) {
        namesTemp <- unique(alltaxalevelpercol.fam[[i]])
        temp <- ""
        for ( j in namesTemp) {
            temp <- paste(temp,j,sep=",")
        }
        groupNames.fam[i] <- temp
        groupAbund[i] <- length(alltaxalevelpercol.fam[[i]])
    }
    
    nameAndAbund <- paste0(groupNames.fam,"(",groupAbund,")")
    
    # pdf(file="LegendForFam.pdf", height=length(nameAndAbund)/10)
    # par(mar=c(4.1,20,2,1))
    # barplot(rep(1,length(allCol.fam)), col=rev(allCol.fam), names.arg = rev(nameAndAbund), horiz = TRUE, las=2, cex.names=0.5, width=rev(groupAbund), border = rev(allCol.fam))
    # dev.off()
    pdf(file=paste0(output,"/LegendForFam.pdf"), height=length(nameAndAbund)/10)
    par(mar=c(4.1,20,2,1))
    barplot(rep(1,length(allCol.fam)), col=(allCol.fam), names.arg = (nameAndAbund), horiz = TRUE, las=2, cex.names=0.5, width=(groupAbund), border = (allCol.fam))
    dev.off()
    
    ##### AT GEN LEVEL LEGEND #####
    alltaxalevelpercol.gen <- list()
    count <- 0
    prevCol <- ""
    for ( i in 1:nrow(alltaxalevelcolors)) {
        currentCol <- alltaxalevelcolors[i,4]
        if (currentCol != prevCol) {
            count  <- count + 1
            alltaxalevelpercol.gen[[count]] <- vector()
        }
        currentotu <- rownames(alltaxalevelcolors)[i]
        currentName <- namesFull[[currentotu]]
        currentClass <- gsub("^.*__","",strsplit(currentName,"; ")[[1]][6])
        alltaxalevelpercol.gen[[count]] <- c(alltaxalevelpercol.gen[[count]],currentClass)
        prevCol <- currentCol
    }
    x <- alltaxalevelcolors[,4]
    allCol.gen <-x[x!=c(x[-1], FALSE)]
    
    #get names of groups, esp if multiple
    groupNames.gen <- vector(length=length(allCol.gen))
    groupAbund <- vector(length=length(allCol.gen))
    for ( i in 1:length(alltaxalevelpercol.gen)) {
        namesTemp <- unique(alltaxalevelpercol.gen[[i]])
        temp <- ""
        for ( j in namesTemp) {
            temp <- paste(temp,j,sep=",")
        }
        groupNames.gen[i] <- temp
        groupAbund[i] <- length(alltaxalevelpercol.gen[[i]])
    }
    
    nameAndAbund <- paste0(groupNames.gen,"(",groupAbund,")")
    
    # pdf(file="LegendForGen.pdf", height = length(nameAndAbund)/10)
    # par(mar=c(4.1,20,2,1))
    # barplot(rep(1,length(allCol.gen)), col=rev(allCol.gen), names.arg = rev(nameAndAbund), horiz = TRUE, las=2, cex.names=0.5, width=rev(groupAbund), border = rev(allCol.gen))
    # dev.off()
    pdf(file=paste0(output,"/LegendForGen.pdf"), height = length(nameAndAbund)/10)
    par(mar=c(4.1,20,2,1))
    barplot(rep(1,length(allCol.gen)), col=(allCol.gen), names.arg = (nameAndAbund), horiz = TRUE, las=2, cex.names=0.5, width=(groupAbund), border = (allCol.gen))
    dev.off()
    
    ###### Triangle legend ##########
    # Coordinates of the triangle
    invtri <- rbind(sin(0:2*2/3*pi), cos(0:2*2/3*pi))
    tri <- rbind(sin(0:2*2/3*pi), c(-1, 0.5,0.5))

    # Function for calculating the color of a set of points `pt`
    # in relation to the triangle
    tricol <- function(pt, sharpness=2){
        # require(splancs)
        RGB <- sapply(1:3, function(i){
            a <- sweep(pt, 2, invtri[,i])
            # b <- apply(invtri[,-i], 1, mean) - invtri[,i]
            b <- apply(invtri[,-i], 1, mean) - invtri[,i]
            output <- sharpness*((a %*% b) / sum(b^2))-sharpness+1
            return(output)
        })
        
        
        # RGB <- rbind(c(0.4,0.53,1.43),c(0.4,-1,0.55),c(-0.4,-1.4,0.6))
        ###### MY ADDITIONS
        RGB.alt <- pmin(pmax(RGB, 0), 1)
        intensity <- apply(RGB.alt, 1, function(x) {
            max(x)
        })
        RGB.alt <- cbind(RGB.alt, intensity)
        RGB.alt[-inpip(pt,t(tri)),] <- 1    # Color points outside the triangle white
        RGB.final <- unname(as.data.frame(RGB.alt))
        do.call(rgb, RGB.final)
        
        #####
        
        # RGB[-inpip(pt,t(tri)),] <- 1    # Color points outside the triangle white
        # do.call(rgb, unname(as.data.frame(pmin(pmax(RGB, 0), 1))))
    }
    
    # Plot
    res <- 1000                         # Resolution
    xi <- seq(-1, 1, length=res)        # Axis points
    yi <- seq(-1.2, 0.8, length=res)
    x <- xi[1] + cumsum(diff(xi))       # Midpoints between axis points
    y <- yi[1] + cumsum(diff(yi))
    xy <- matrix(1:(length(x)*length(y)), length(x))
    
    
    
    ############# PLOT ###########
    print("about to plot...")
    # phylogenetic plot
    # P <- ggtree(tree.pruned
    #             # ,branch.length="none"
    #             ,layout='rectangular')
    # NOTE; took the above and using it to sort taxa.
    # create an apporpriate viewport.  Modify the dimensions and coordinates as needed
    
    labelsLevels <- c("Phylum","Class","Order","Family","Genus","Species")
    if ( d == "F18") {
        labelsLevels <- c("D_01","D_02","D_03","D_04","D_05","D_06")
    }
    
    pdf(file = paste0(output,"/Composition_of_taxa_types.pdf"),5,7)
    # quartz(,5,7)
    par(fig = c(0,1,0,0.4), mar = c(0,0,0,0), oma = c(2,4,2,4))
    image(xi, yi, xy, col=tricol(as.matrix(expand.grid(x,y))), useRaster=TRUE
          , axes = FALSE
          , xlab = ""
          , ylab = "")
    lines(tri[1,c(1:3,1)], tri[2,c(1:3,1)], type="l")
    
    par(fig = c(0,1,0,0.4), mar = c(0,0,0,0), oma = c(0,0,0,0), new = TRUE)
    plot(0,0
         , xaxt = "n"
         , yaxt = "n"
         , bty = "n"
         , pch = ""
         , xlab = ""
         , ylab = "")
    text(x = c(-0.8,0.8,0)
         , y = c(0.9,0.9,-0.75)
         , labels = c("BRACK","FRESH","MARINE")
    )
    par(fig=c(0,0.4,0.4,1), mar=c(6,1,1,1), new=TRUE)
    plot(0,0,bty="n",xaxt="n",xlab="",ylab="",pch="",yaxt="n")
    vp.TopLeft <- viewport(height=unit(0.42, "npc"), width=unit(0.425, "npc"), 
                           just=c("left","top"), 
                           y=0.98, x=0)
    print(P, vp=vp.TopLeft)
    
    
    par(fig=c(0.2+0.2,0.3+0.2,0.4,1), mar=c(6,.1,1,.1), new=TRUE)
    barplot(as.matrix(rev(blankBar)), col="lightgrey", border = "white",axes = FALSE)
    axis(1,at = 0.7, labels="Taxa", las=2, tick=FALSE)
    par(fig=c(0.3+0.2,0.4+0.2,0.4,1), mar=c(6,.1,1,.1),new=TRUE)
    barplot(as.matrix(rev(emptyBarMat[,1])), col=(alltaxalevelcolors.final[,1]),border=NA, axes = FALSE)
    axis(1,at = 0.7, labels=labelsLevels[2], las=2, tick=FALSE)
    par(fig=c(0.4+0.2,0.5+0.2,0.4,1), mar=c(6,.1,1,.1),new=TRUE)
    barplot(as.matrix(rev(emptyBarMat[,2])), col=(alltaxalevelcolors.final[,2]),border=NA, axes = FALSE)
    axis(1,at = 0.7, labels=labelsLevels[3], las=2, tick=FALSE)
    par(fig=c(0.5+0.2,0.6+0.2,0.4,1), mar=c(6,.1,1,.1),new=TRUE)
    barplot(as.matrix(rev(emptyBarMat[,3])), col=(alltaxalevelcolors.final[,3]),border=NA, axes = FALSE)
    axis(1,at = 0.7, labels=labelsLevels[4], las=2, tick=FALSE)
    par(fig=c(0.6+0.2,0.7+0.2,0.4,1), mar=c(6,.1,1,.1),new=TRUE)
    barplot(as.matrix(rev(emptyBarMat[,4])), col=(alltaxalevelcolors.final[,4]),border=NA, axes = FALSE)
    axis(1,at = 0.7, labels=labelsLevels[5], las=2, tick=FALSE)
    par(fig=c(0.7+0.2,0.8+0.2,0.4,1), mar=c(6,.1,1,.1),new=TRUE)
    barplot(as.matrix(rev(emptyBarMat[,5])), col=(alltaxalevelcolors.final[,5]),border=NA, axes = FALSE)
    axis(1,at = 0.7, labels="OTU", las=2, tick=FALSE)
    
    dev.off()

    ### NO LEGEND PLOT
     pdf(file = paste0(output,"/Composition_of_taxa_types_nolegend.pdf"),7,5)
     # quartz(7,5)
     par(fig=c(0,0.4,0,1), mar=c(6,.1,1,.1))
     plot(0,0,bty="n",xaxt="n",xlab="",ylab="",pch="",yaxt="n")
     vp.TopLeft <- viewport(height=unit(0.75, "npc"), width=unit(0.425, "npc"), 
                            just=c("left","top"), 
                            y=0.975, x=0)
     print(P, vp=vp.TopLeft)
     
     par(fig=c(0.2+0.2,0.3+0.2,0,1), mar=c(6,.1,1,.1), new=TRUE)
     barplot(as.matrix(rev(blankBar)), col="lightgrey", border = "white",axes = FALSE)
     axis(1,at = 0.7, labels="Taxa", las=2, tick=FALSE)
     par(fig=c(0.3+0.2,0.4+0.2,0,1), mar=c(6,.1,1,.1),new=TRUE)
     barplot(as.matrix(rev(emptyBarMat[,1])), col=(alltaxalevelcolors.final[,1]),border=NA, axes = FALSE)
     axis(1,at = 0.7, labels=labelsLevels[2], las=2, tick=FALSE)
     par(fig=c(0.4+0.2,0.5+0.2,0,1), mar=c(6,.1,1,.1),new=TRUE)
     barplot(as.matrix(rev(emptyBarMat[,2])), col=(alltaxalevelcolors.final[,2]),border=NA, axes = FALSE)
     axis(1,at = 0.7, labels=labelsLevels[3], las=2, tick=FALSE)
     par(fig=c(0.5+0.2,0.6+0.2,0,1), mar=c(6,.1,1,.1),new=TRUE)
     barplot(as.matrix(rev(emptyBarMat[,3])), col=(alltaxalevelcolors.final[,3]),border=NA, axes = FALSE)
     axis(1,at = 0.7, labels=labelsLevels[4], las=2, tick=FALSE)
     par(fig=c(0.6+0.2,0.7+0.2,0,1), mar=c(6,.1,1,.1),new=TRUE)
     barplot(as.matrix(rev(emptyBarMat[,4])), col=(alltaxalevelcolors.final[,4]),border=NA, axes = FALSE)
     axis(1,at = 0.7, labels=labelsLevels[5], las=2, tick=FALSE)
     par(fig=c(0.7+0.2,0.8+0.2,0,1), mar=c(6,.1,1,.1),new=TRUE)
     barplot(as.matrix(rev(emptyBarMat[,5])), col=(alltaxalevelcolors.final[,5]),border=NA, axes = FALSE)
     axis(1,at = 0.7, labels="OTU", las=2, tick=FALSE)
     # barplot(as.matrix(emptyBarMat[,6]), col=rev(alltaxalevelcolors.final[,6]),border=NA, axes = FALSE)
     # axis(1,at = 0.7, labels="OTU", las=2, tick=FALSE)
     
     dev.off()
     # ## OLD; for some reason the colors are backward
     # pdf(file = "Composition_of_taxa_types_nolegend.pdf",7,5)
     # # quartz(7,5)
     # par(fig=c(0,0.4,0,1), mar=c(6,.1,1,.1))
     # plot(0,0,bty="n",xaxt="n",xlab="",ylab="",pch="",yaxt="n")
     # vp.TopLeft <- viewport(height=unit(0.82, "npc"), width=unit(0.425, "npc"), 
     #                        just=c("left","top"), 
     #                        y=1.01, x=0)
     # print(P, vp=vp.TopLeft)
     # 
     # par(fig=c(0.2+0.2,0.3+0.2,0,1), mar=c(6,.1,1,.1), new=TRUE)
     # barplot(as.matrix(blankBar), col="lightgrey", border = "white",axes = FALSE)
     # axis(1,at = 0.7, labels="Taxa", las=2, tick=FALSE)
     # par(fig=c(0.3+0.2,0.4+0.2,0,1), mar=c(6,.1,1,.1),new=TRUE)
     # barplot(as.matrix(emptyBarMat[,1]), col=rev(alltaxalevelcolors.final[,1]),border=NA, axes = FALSE)
     # axis(1,at = 0.7, labels="Class", las=2, tick=FALSE)
     # par(fig=c(0.4+0.2,0.5+0.2,0,1), mar=c(6,.1,1,.1),new=TRUE)
     # barplot(as.matrix(emptyBarMat[,2]), col=rev(alltaxalevelcolors.final[,2]),border=NA, axes = FALSE)
     # axis(1,at = 0.7, labels="Order", las=2, tick=FALSE)
     # par(fig=c(0.5+0.2,0.6+0.2,0,1), mar=c(6,.1,1,.1),new=TRUE)
     # barplot(as.matrix(emptyBarMat[,3]), col=rev(alltaxalevelcolors.final[,3]),border=NA, axes = FALSE)
     # axis(1,at = 0.7, labels="Family", las=2, tick=FALSE)
     # par(fig=c(0.6+0.2,0.7+0.2,0,1), mar=c(6,.1,1,.1),new=TRUE)
     # barplot(as.matrix(emptyBarMat[,4]), col=rev(alltaxalevelcolors.final[,4]),border=NA, axes = FALSE)
     # axis(1,at = 0.7, labels="Genus", las=2, tick=FALSE)
     # par(fig=c(0.7+0.2,0.8+0.2,0,1), mar=c(6,.1,1,.1),new=TRUE)
     # barplot(as.matrix(emptyBarMat[,5]), col=rev(alltaxalevelcolors.final[,5]),border=NA, axes = FALSE)
     # axis(1,at = 0.7, labels="OTU", las=2, tick=FALSE)
     # # barplot(as.matrix(emptyBarMat[,6]), col=rev(alltaxalevelcolors.final[,6]),border=NA, axes = FALSE)
     # # axis(1,at = 0.7, labels="OTU", las=2, tick=FALSE)
     # 
     # dev.off()
     print("finished plotting")
    ####### NEAREST NEIGHBOUR #########
    
    # get distance matrix
    PD.dist <- cophenetic(tree.pruned)
    diag(PD.dist) <- NA
    
    allBrack <- as.character(mb[mb[,"typeSimple"] == "brackishRestricted",1])
    
    nearestNeigh <- rep(NA,length(allBrack))
    n <- 1
    for ( otu in (allBrack) ) {
        if ( otu %in% rownames(PD.dist)) {
            nearestType <- "brackishRestricted"
            # order closest to farthest
            orderNeigh <- order(PD.dist[otu,])
            count <- 1
            while ( nearestType == "brackishRestricted" | nearestType == "noclass" ) {
                nearestOTU <- colnames(PD.dist)[orderNeigh[count]]
                nearestType <- as.character(mb[mb[,1]==nearestOTU,"typeSimple"])
                count <- count + 1
            }
            
            nearestNeigh[n] <- nearestType
        }
        n <- n+1
    }
    ratioFreshMar <- table(nearestNeigh)
    
    write(capture.output(binom.test(x=ratioFreshMar[1], n=sum(ratioFreshMar), p=0.5)), file=paste0(output,"/binotest.txt"))
    
    write.table(blocksizeNames, file = paste0(output,"/LargeBlockNames.txt"), sep = "\t")
   
   
    
    #### FINISH UP ####
    
    print("ready to start next dataset")
    

}


########## COMBINE BALTIC AND FRASER TREE ############

treepwd <- 'tree_BF16_filt.tre'
#Readtree
tree16 <- read.tree(file = paste0(treepwd))
#filter tree

#Baltic16
mbpwdB <- paste0(dirNames[1],'/modelBoundaries_type.txt')
otupwdB <- paste0(dirNames[1],'/OTUTableText.txt')
uniqueotupwdB <- paste0(dirNames[1],'/OTUs.txt')
#Fraser16
mbpwdF <- paste0(dirNames[2],'/modelBoundaries_type.txt')
otupwdF <- paste0(dirNames[2],'/OTUTableText.txt')
uniqueotupwdF <- paste0(dirNames[2],'/OTUs.txt')
output <- paste0(newDir,"/COMBO16")

# load mb
mbB <- read.delim(mbpwdB, header=TRUE, row.names=1)
mbF <- read.delim(mbpwdF, header=TRUE, row.names=1)
# load OTU tab
otuB <- read.delim(otupwdB, header=TRUE, row.names=1, skip=1)
otuF <- read.delim(otupwdF, header=TRUE, row.names=1, skip=1)
# get all unique OTUs
allOTUs <- unique(c(rownames(mbB), rownames(mbF)))

# get all shared OTUs; alternative
allOTUs <- c(rownames(mbB), rownames(mbF))[duplicated(c(rownames(mbB), rownames(mbF)))]

# get consenses types
conType <- matrix(ncol=1, nrow=length(allOTUs), dimnames = list(allOTUs, c("typeSimple")))
countNoMatch <- 0
countMatch <- 0
NoMatchTrack <- matrix(ncol=2)
for ( otu in allOTUs ) {
    typeB <- NA
    typeF <- NA
    if ( otu %in% rownames(mbB)) {
        typeB <- as.character(mbB[otu,"typeSimple"])
    }
    if ( otu %in% rownames(mbF)) {
        typeF <- as.character(mbF[otu,"typeSimple"])
    }
    # Now get consensus type
    if (is.na(typeB)) {
        conType[otu,] <- typeF
    } else if (is.na(typeF)) {
        conType[otu,] <- typeB
    } else if (typeF==typeB) {
        conType[otu,] <- typeB
        countMatch <- countMatch + 1
    } else if (typeF != typeB) {
        conType[otu,] <- "noclass"
        countNoMatch <- countNoMatch + 1 # keep track of how many don't match
        NoMatchTrack <- rbind(NoMatchTrack,c(typeB,typeF))
    }
}
# conType
# countMatch
# countNoMatch
# match(rownames(conType),tree16$tip.label)
# rownames(conType)[24]

conType <- conType[-which(conType[,1]=="noclass"),]
conType <- as.matrix(conType)

# prune tree
tree.pruned <- prune.sample(t(conType), tree16)

# get OTU taxonomies
taxaRefB <- cbind(as.character(otuB[,ncol(otuB)]))
rownames(taxaRefB) <- rownames(otuB)
taxaRefF <- cbind(as.character(otuF[,ncol(otuF)]))
rownames(taxaRefF) <- rownames(otuF)

# filter OTUtaxonomies
otusCombo <- rownames(conType)
taxaRefComb <- matrix(ncol=1,nrow=length(otusCombo),dimnames=list(otusCombo, "taxonomy"))
for ( otu in otusCombo ) {
    if ( otu %in% rownames(taxaRefB) ) {
        taxaRefComb[otu,1] <- taxaRefB[otu,]
    } else if ( otu %in% rownames(taxaRefF) ) {
        taxaRefComb[otu,1] <- taxaRefF[otu,]
    }
}

#

# get order of taxa that GGtree will plot
P <- ggtree(tree.pruned
            ,branch.length="none"
            ,layout='rectangular')
toOrderGGtree <- P$data[which(P$data$isTip == TRUE),c("y","label")]
tips <- rev(toOrderGGtree[order(toOrderGGtree$y),"label"])

# tips <- tree.pruned$tip.label

###### BEGIN PLOTTING #######

# get full taxa name
namesFull <- lapply(1:length(tips$label), function(x) {
    fullname <- taxaRefComb[match(tips$label[x],rownames(taxaRefComb)),1]
    return(fullname)
})
names(namesFull) <- tips$label


# Get composition of each
mb <- data.frame("taxa"=rownames(conType),"typeSimple"=conType[,1])
taxaNames <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
taxaLevels <- c(3,4,5,6) # Don't do species bc it doesn't add a lot of information
alltaxalevelcolors <- matrix(ncol=length(taxaLevels),nrow=length(tips$label),dimnames = list(tips$label,taxaNames[taxaLevels]))
for (l in taxaLevels) {
    listClassComp <- propType(level=l,namesFull,mb)
    tableTypes <- tableAllTypes(level=l,namesFull)
    colors <- getColor(listClassComp)
    alltaxalevelcolors[,taxaNames[l]] <- as.vector(makebarcolor(colors,tableTypes,taxaNames[l]))
    if (l == 3) {
        listClassCompClass <- listClassComp
        tableTypesClass <- tableTypes
    } else if (l == 4) {
        listClassCompOrder <- listClassComp
        tableTypesOrder <- tableTypes
    } else if (l == 5 ) {
        listClassCompFam <- listClassComp
        tableTypesFam <- tableTypes
    }
}
# Get class sizes so I can make big-ass blocks
keepname <- (lapply(listClassCompClass, FUN=sum) != 0)
ClassNames <- gsub("_","",names(listClassCompClass)[keepname])
ttc.sorted <- tableTypesClass[match(ClassNames, gsub("_","",names(tableTypesClass)))]
blocksizeNames <- cbind(ClassNames,ttc.sorted)

#Also, get composition of OTU
alltypes <- mb[match(tips, mb[,"taxa"]),"typeSimple"]
tempTab <- lapply(alltypes,function(x) {table(x,exclude="noclass")})
OTU <- getColor(tempTab)

#get class level to label
listClassComp <- propType(level=3,namesFull,mb)
tableTypes <- tableAllTypes(level=3,namesFull)
colors <- getColor(listClassComp)

# make empty matrix of 1's
emptyBarMat <- matrix(ncol=length(taxaLevels)+1,nrow=length(tips$label),dimnames = list(tips$label,c(taxaNames[taxaLevels],"OTU")), data = 1)
# add OTU colors
alltaxalevelcolors.final <- cbind(alltaxalevelcolors,"OTU" = OTU)

#####################
# legend for class


# Get class of each; then plot each color and all things in that block
##### AT CLASS LEVEL LEGEND #####
alltaxalevelpercol <- list()
count <- 0
prevCol <- ""
for ( i in 1:nrow(alltaxalevelcolors)) {
    currentCol <- alltaxalevelcolors[i,1]
    if (currentCol != prevCol) {
        count  <- count + 1
        alltaxalevelpercol[[count]] <- vector()
    }
    currentotu <- rownames(alltaxalevelcolors)[i]
    currentName <- namesFull[[currentotu]]
    currentClass <- gsub("^.*__","",strsplit(currentName,"; ")[[1]][3])
    alltaxalevelpercol[[count]] <- c(alltaxalevelpercol[[count]],currentClass)
    prevCol <- currentCol
}
x <- alltaxalevelcolors[,1]
allCol <-x[x!=c(x[-1], FALSE)]

#get names of groups, esp if multiple
groupNames <- vector(length=length(allCol))
groupAbund <- vector(length=length(allCol))
for ( i in 1:length(alltaxalevelpercol)) {
    namesTemp <- unique(alltaxalevelpercol[[i]])
    temp <- ""
    for ( j in namesTemp) {
        temp <- paste(temp,j,sep=",")
    }
    groupNames[i] <- temp
    groupAbund[i] <- length(alltaxalevelpercol[[i]])
}

nameAndAbund <- paste0(groupNames,"(",groupAbund,")")

pdf(file=paste0(output,"/LegendForClass.pdf"), height=length(nameAndAbund)/2)
par(mar=c(4.1,20,2,1))
barplot(rep(1,length(allCol)), col=(allCol), names.arg = (nameAndAbund), horiz = TRUE, las=2, cex.names=0.5, width=(groupAbund), border = (allCol))
dev.off()

# Now, make a grey barplot with data above
blankBar <- rev(sapply(alltaxalevelpercol, function(x) length(x)))

##### AT ORDER LEVEL LEGEND #####
alltaxalevelpercol.order <- list()
count <- 0
prevCol <- ""
for ( i in 1:nrow(alltaxalevelcolors)) {
    currentCol <- alltaxalevelcolors[i,2]
    if (currentCol != prevCol) {
        count  <- count + 1
        alltaxalevelpercol.order[[count]] <- vector()
    }
    currentotu <- rownames(alltaxalevelcolors)[i]
    currentName <- namesFull[[currentotu]]
    currentClass <- gsub("^.*__","",strsplit(currentName,"; ")[[1]][4])
    alltaxalevelpercol.order[[count]] <- c(alltaxalevelpercol.order[[count]],currentClass)
    prevCol <- currentCol
}
x <- alltaxalevelcolors[,2]
allCol.order <-x[x!=c(x[-1], FALSE)]

#get names of groups, esp if multiple
groupNames.order <- vector(length=length(allCol.order))
groupAbund <- vector(length=length(allCol.order))
for ( i in 1:length(alltaxalevelpercol.order)) {
    namesTemp <- unique(alltaxalevelpercol.order[[i]])
    temp <- ""
    for ( j in namesTemp) {
        temp <- paste(temp,j,sep=",")
    }
    groupNames.order[i] <- temp
    groupAbund[i] <- length(alltaxalevelpercol.order[[i]])
}

nameAndAbund <- paste0(groupNames.order,"(",groupAbund,")")

pdf(file=paste0(output,"/LegendForOrder.pdf"), height=length(nameAndAbund)/2)
par(mar=c(4.1,20,2,1))
barplot(rep(1,length(allCol.order)), col=(allCol.order), names.arg = (nameAndAbund), horiz = TRUE, las=2, cex.names=0.5, width=(groupAbund), border = (allCol.order))
dev.off()



##### AT FAM LEVEL LEGEND #####
alltaxalevelpercol.fam <- list()
count <- 0
prevCol <- ""
for ( i in 1:nrow(alltaxalevelcolors)) {
    currentCol <- alltaxalevelcolors[i,3]
    if (currentCol != prevCol) {
        count  <- count + 1
        alltaxalevelpercol.fam[[count]] <- vector()
    }
    currentotu <- rownames(alltaxalevelcolors)[i]
    currentName <- namesFull[[currentotu]]
    currentClass <- gsub("^.*__","",strsplit(currentName,"; ")[[1]][5])
    alltaxalevelpercol.fam[[count]] <- c(alltaxalevelpercol.fam[[count]],currentClass)
    prevCol <- currentCol
}
x <- alltaxalevelcolors[,3]
allCol.fam <-x[x!=c(x[-1], FALSE)]

#get names of groups, esp if multiple
groupNames.fam <- vector(length=length(allCol.fam))
groupAbund <- vector(length=length(allCol.fam))
for ( i in 1:length(alltaxalevelpercol.fam)) {
    namesTemp <- unique(alltaxalevelpercol.fam[[i]])
    temp <- ""
    for ( j in namesTemp) {
        temp <- paste(temp,j,sep=",")
    }
    groupNames.fam[i] <- temp
    groupAbund[i] <- length(alltaxalevelpercol.fam[[i]])
}

nameAndAbund <- paste0(groupNames.fam,"(",groupAbund,")")

pdf(file=paste0(output,"/LegendForFam.pdf"), height=length(nameAndAbund)/2)
par(mar=c(4.1,20,2,1))
barplot(rep(1,length(allCol.fam)), col=(allCol.fam), names.arg = (nameAndAbund), horiz = TRUE, las=2, cex.names=0.5, width=(groupAbund), border = (allCol.fam))
dev.off()

##### AT GEN LEVEL LEGEND #####
alltaxalevelpercol.gen <- list()
count <- 0
prevCol <- ""
for ( i in 1:nrow(alltaxalevelcolors)) {
    currentCol <- alltaxalevelcolors[i,4]
    if (currentCol != prevCol) {
        count  <- count + 1
        alltaxalevelpercol.gen[[count]] <- vector()
    }
    currentotu <- rownames(alltaxalevelcolors)[i]
    currentName <- namesFull[[currentotu]]
    currentClass <- gsub("^.*__","",strsplit(currentName,"; ")[[1]][6])
    alltaxalevelpercol.gen[[count]] <- c(alltaxalevelpercol.gen[[count]],currentClass)
    prevCol <- currentCol
}
x <- alltaxalevelcolors[,4]
allCol.gen <-x[x!=c(x[-1], FALSE)]

#get names of groups, esp if multiple
groupNames.gen <- vector(length=length(allCol.gen))
groupAbund <- vector(length=length(allCol.gen))
for ( i in 1:length(alltaxalevelpercol.gen)) {
    namesTemp <- unique(alltaxalevelpercol.gen[[i]])
    temp <- ""
    for ( j in namesTemp) {
        temp <- paste(temp,j,sep=",")
    }
    groupNames.gen[i] <- temp
    groupAbund[i] <- length(alltaxalevelpercol.gen[[i]])
}

nameAndAbund <- paste0(groupNames.gen,"(",groupAbund,")")

pdf(file=paste0(output,"/LegendForGen.pdf"), height = length(nameAndAbund)/2)
par(mar=c(4.1,20,2,1))
barplot(rep(1,length(allCol.gen)), col=(allCol.gen), names.arg = (nameAndAbund), horiz = TRUE, las=2, cex.names=0.5, width=(groupAbund), border = (allCol.gen))
dev.off()


###### Triangle legend ##########
# Coordinates of the triangle
invtri <- rbind(sin(0:2*2/3*pi), cos(0:2*2/3*pi))
tri <- rbind(sin(0:2*2/3*pi), c(-1, 0.5,0.5))

# Function for calculating the color of a set of points `pt`
# in relation to the triangle
tricol <- function(pt, sharpness=2){
    # require(splancs)
    RGB <- sapply(1:3, function(i){
        a <- sweep(pt, 2, invtri[,i])
        # b <- apply(invtri[,-i], 1, mean) - invtri[,i]
        b <- apply(invtri[,-i], 1, mean) - invtri[,i]
        output <- sharpness*((a %*% b) / sum(b^2))-sharpness+1
        return(output)
    })
    
    
    # RGB <- rbind(c(0.4,0.53,1.43),c(0.4,-1,0.55),c(-0.4,-1.4,0.6))
    ###### MY ADDITIONS
    RGB.alt <- pmin(pmax(RGB, 0), 1)
    intensity <- apply(RGB.alt, 1, function(x) {
        max(x)
    })
    RGB.alt <- cbind(RGB.alt, intensity)
    RGB.alt[-inpip(pt,t(tri)),] <- 1    # Color points outside the triangle white
    RGB.final <- unname(as.data.frame(RGB.alt))
    do.call(rgb, RGB.final)
    
    #####
    
    # RGB[-inpip(pt,t(tri)),] <- 1    # Color points outside the triangle white
    # do.call(rgb, unname(as.data.frame(pmin(pmax(RGB, 0), 1))))
}

# Plot
res <- 1000                         # Resolution
xi <- seq(-1, 1, length=res)        # Axis points
yi <- seq(-1.2, 0.8, length=res)
x <- xi[1] + cumsum(diff(xi))       # Midpoints between axis points
y <- yi[1] + cumsum(diff(yi))
xy <- matrix(1:(length(x)*length(y)), length(x))


############# PLOT ###########

# phylogenetic plot
# P <- ggtree(tree.pruned,branch.length="none",layout='rectangular')
# create an apporpriate viewport.  Modify the dimensions and coordinates as needed

pdf(file = paste0(output,"/Composition_of_taxa_types.pdf"),5,7)
# quartz(,5,7)
par(fig = c(0,1,0,0.4), mar = c(0,0,0,0), oma = c(2,4,2,4))
image(xi, yi, xy, col=tricol(as.matrix(expand.grid(x,y))), useRaster=TRUE
      , axes = FALSE
      , xlab = ""
      , ylab = "")
lines(tri[1,c(1:3,1)], tri[2,c(1:3,1)], type="l")

par(fig = c(0,1,0,0.4), mar = c(0,0,0,0), oma = c(0,0,0,0), new = TRUE)
plot(0,0
     , xaxt = "n"
     , yaxt = "n"
     , bty = "n"
     , pch = ""
     , xlab = ""
     , ylab = "")
text(x = c(-0.8,0.8,0)
     , y = c(0.9,0.9,-0.75)
     , labels = c("BRACK","FRESH","MARINE")
)
par(fig=c(0,0.4,0.4,1), mar=c(6,1,1,1), new=TRUE)
plot(0,0,bty="n",xaxt="n",xlab="",ylab="",pch="",yaxt="n")
vp.TopLeft <- viewport(height=unit(0.425, "npc"), width=unit(0.425, "npc"), 
                       just=c("left","top"), 
                       y=0.9835, x=0)

print(P, vp=vp.TopLeft)

par(fig=c(0.2+0.2,0.3+0.2,0.4,1), mar=c(6,.1,1,.1), new=TRUE)
barplot(as.matrix(rev(blankBar)), col="lightgrey", border = "white",axes = FALSE)
axis(1,at = 0.7, labels="Taxa", las=2, tick=FALSE)
par(fig=c(0.3+0.2,0.4+0.2,0.4,1), mar=c(6,.1,1,.1),new=TRUE)
barplot(as.matrix(rev(emptyBarMat[,1])), col=(alltaxalevelcolors.final[,1]),border=NA, axes = FALSE)
axis(1,at = 0.7, labels="Class", las=2, tick=FALSE)
par(fig=c(0.4+0.2,0.5+0.2,0.4,1), mar=c(6,.1,1,.1),new=TRUE)
barplot(as.matrix(rev(emptyBarMat[,2])), col=(alltaxalevelcolors.final[,2]),border=NA, axes = FALSE)
axis(1,at = 0.7, labels="Order", las=2, tick=FALSE)
par(fig=c(0.5+0.2,0.6+0.2,0.4,1), mar=c(6,.1,1,.1),new=TRUE)
barplot(as.matrix(rev(emptyBarMat[,3])), col=(alltaxalevelcolors.final[,3]),border=NA, axes = FALSE)
axis(1,at = 0.7, labels="Family", las=2, tick=FALSE)
par(fig=c(0.6+0.2,0.7+0.2,0.4,1), mar=c(6,.1,1,.1),new=TRUE)
barplot(as.matrix(rev(emptyBarMat[,4])), col=(alltaxalevelcolors.final[,4]),border=NA, axes = FALSE)
axis(1,at = 0.7, labels="Genus", las=2, tick=FALSE)
par(fig=c(0.7+0.2,0.8+0.2,0.4,1), mar=c(6,.1,1,.1),new=TRUE)
barplot(as.matrix(rev(emptyBarMat[,5])), col=(alltaxalevelcolors.final[,5]),border=NA, axes = FALSE)
axis(1,at = 0.7, labels="OTU", las=2, tick=FALSE)
# barplot(as.matrix(emptyBarMat[,6]), col=rev(alltaxalevelcolors.final[,6]),border=NA, axes = FALSE)
# axis(1,at = 0.7, labels="OTU", las=2, tick=FALSE)

dev.off()

### NO LEGEND PLOT
pdf(file = paste0(output,"/Composition_of_taxa_types_nolegend.pdf"),7,5)
# quartz(7,5)
par(fig=c(0,0.4,0,1), mar=c(6,.1,1,.1))
plot(0,0,bty="n",xaxt="n",xlab="",ylab="",pch="",yaxt="n")
vp.TopLeft <- viewport(height=unit(0.76, "npc"), width=unit(0.425, "npc"), 
                       just=c("left","top"), 
                       y=0.98, x=0)
print(P, vp=vp.TopLeft)

par(fig=c(0.2+0.2,0.3+0.2,0,1), mar=c(6,.1,1,.1), new=TRUE)
barplot(as.matrix(rev(blankBar)), col="lightgrey", border = "white",axes = FALSE)
axis(1,at = 0.7, labels="Taxa", las=2, tick=FALSE)
par(fig=c(0.3+0.2,0.4+0.2,0,1), mar=c(6,.1,1,.1),new=TRUE)
barplot(as.matrix(rev(emptyBarMat[,1])), col=(alltaxalevelcolors.final[,1]),border=NA, axes = FALSE)
axis(1,at = 0.7, labels="Class", las=2, tick=FALSE)
par(fig=c(0.4+0.2,0.5+0.2,0,1), mar=c(6,.1,1,.1),new=TRUE)
barplot(as.matrix(rev(emptyBarMat[,2])), col=(alltaxalevelcolors.final[,2]),border=NA, axes = FALSE)
axis(1,at = 0.7, labels="Order", las=2, tick=FALSE)
par(fig=c(0.5+0.2,0.6+0.2,0,1), mar=c(6,.1,1,.1),new=TRUE)
barplot(as.matrix(rev(emptyBarMat[,3])), col=(alltaxalevelcolors.final[,3]),border=NA, axes = FALSE)
axis(1,at = 0.7, labels="Family", las=2, tick=FALSE)
par(fig=c(0.6+0.2,0.7+0.2,0,1), mar=c(6,.1,1,.1),new=TRUE)
barplot(as.matrix(rev(emptyBarMat[,4])), col=(alltaxalevelcolors.final[,4]),border=NA, axes = FALSE)
axis(1,at = 0.7, labels="Genus", las=2, tick=FALSE)
par(fig=c(0.7+0.2,0.8+0.2,0,1), mar=c(6,.1,1,.1),new=TRUE)
barplot(as.matrix(rev(emptyBarMat[,5])), col=(alltaxalevelcolors.final[,5]),border=NA, axes = FALSE)
axis(1,at = 0.7, labels="OTU", las=2, tick=FALSE)
# barplot(as.matrix(emptyBarMat[,6]), col=rev(alltaxalevelcolors.final[,6]),border=NA, axes = FALSE)
# axis(1,at = 0.7, labels="OTU", las=2, tick=FALSE)

dev.off()
write.table(blocksizeNames, file = paste0(output,"/LargeBlockNames.txt"), sep = "\t")

pdf(file = paste0(output,"/LEGEND_for_tree.pdf"),7,5)
# quartz(,5,7)
par(fig = c(0,1,0,1), mar = c(0,0,0,0), oma = c(2,4,2,4))
image(xi, yi, xy, col=tricol(as.matrix(expand.grid(x,y))), useRaster=TRUE
      , axes = FALSE
      , xlab = ""
      , ylab = "")
lines(tri[1,c(1:3,1)], tri[2,c(1:3,1)], type="l")

par(fig = c(0,1,0,1), mar = c(0,0,0,0), oma = c(0,0,0,0), new = TRUE)
plot(0,0
     , xaxt = "n"
     , yaxt = "n"
     , bty = "n"
     , pch = ""
     , xlab = ""
     , ylab = "")
text(x = c(-0.8,0.8,0)
     , y = c(0.75,0.75,-0.8)
     , labels = c("BRACK","FRESH","MARINE")
)
dev.off()


