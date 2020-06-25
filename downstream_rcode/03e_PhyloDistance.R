################
# This script takes SalBin output and calculates the pairwise phylogenetic distances between OTU types

# Getting PD distances
library(picante) # for phylogenetic tree
library(MASS) # for stepAIC
library(tidyverse)


###### Iterate through datasets #########

dirNames <- unlist(read.delim("./OUTPUT/dirNames.txt", header=FALSE))
newDirName <- "./OUTPUT/03e_PhyloDistance"
dir.create(newDirName)
dataset <- c("B16","F16","F18")
for ( d in dataset) {
    if (d == "B16") {
        print("STARTING BALTIC DATASET")
        #Baltic
        mbpwd <- paste0(dirNames[1],'/modelBoundaries_type.txt')
        otupwd <- paste0(dirNames[1],'/OTUTableText.txt')
        uniqueotupwd <- paste0(dirNames[1],'/OTUs.txt')
        output <- paste0(newDirName,"/BALTIC16_PD")
        treepwd <- './OUTPUT/tree_B16_filt.tre'
        tolRpwd <- paste0(dirNames[1],'/toleranceRanges.txt')
        
        midBrackThresh <- 15
        
    } else if (d == "F16") {
        print("STARTING FRASER 16 DATASET")
        #Fraser
        mbpwd <- paste0(dirNames[2],'/modelBoundaries_type.txt')
        otupwd <- paste0(dirNames[2],'/OTUTableText.txt')
        uniqueotupwd <- paste0(dirNames[2],'/OTUs.txt')
        output <- paste0(newDirName,"/FRASER16_PD")
        treepwd <- './OUTPUT/tree_F16_filt.tre'
        tolRpwd <- paste0(dirNames[2],'/toleranceRanges.txt')
        
        midBrackThresh <- 15
        
    } else if (d == "F18") {
        print("STARTING FRASER 18 DATASET")
        #Fraser18
        mbpwd <- paste0(dirNames[3],'/modelBoundaries_type.txt')
        otupwd <- paste0(dirNames[3],'/OTUTableText.txt')
        uniqueotupwd <- paste0(dirNames[3],'/OTUs.txt')
        output <- paste0(newDirName,"/FRASER18_PD")
        treepwd <- './OUTPUT/tree_F18_filt.tre'
        tolRpwd <- paste0(dirNames[3],'/toleranceRanges.txt')
        midBrackThresh <- 15
        
    } 
    
    
    ###### Load tree and data and functions #########
    
    tree16 <- read.tree(file = paste0(treepwd))
    
    mb <- read.delim(file=paste0(mbpwd))
    otu <- read.delim(file=paste0(otupwd), row.names=1)
    uniqueotu <- read.delim(file=uniqueotupwd, header=FALSE)
    tolR <- read.delim(file=tolRpwd)
    
    ####################################
    # make new dir
    dir.create(paste0(output))

    taxa <- cbind(rownames(otu), as.character(otu[,ncol(otu)]))
    otu <- otu[,-ncol(otu)]
    
    # Filter OTU table to get rid of all OTUs not in uniqueotu; should automatically filter out low abund OTU
    # I really shouldn't NEED this; but here it is in case.
    # otu.filt <- otu[match(uniqueotu$V1, rownames(otu)),]
    # for 18, need to filter out things that aren't in the tree-- some alignments failed.
    
    otu.filt <- otu[tree16$tip.label,]
    
    # get rid of any NA rows
    
    mb <- mb[match(rownames(otu.filt),mb$taxa),]
    mb <- mb[which(!is.na(mb$taxa)),]
    # Prune tree so there are only things in the table.
    tree.pruned <- prune.sample(t(otu.filt), tree16)
    
    # Find pairwise distances between all
    taxa.dist <- cophenetic(tree.pruned)
    
    # # First, adjust taxa.dist so that diagonal is NA and only includes things in mb
    # taxa.dist.2 <- taxa.dist
    # diag(taxa.dist.2) <- NA
    # # Reorder
    # taxa.dist.2 <- taxa.dist.2[as.character(mb[,1]), as.character(mb[,1])]
    # 
    # # Basics
    fresh.taxa <- mb[mb[,3] == "freshRestricted","taxa"]
    brackish.taxa <- mb[mb[,3] == "brackishRestricted","taxa"]
    marine.taxa <- mb[mb[,3] == "marineRestricted","taxa"]

    ## Split brackish; upper and lower
    # upper.brackish.taxa <-  mb[mb[,2] == "brackishPeakHiToler","taxa"]
    # lower.brackish.taxa <-  mb[mb[,2] == "brackishPeakLoToler","taxa"]
    
    # Split brackish; upper and lower
    brackishSplit <- tolR %>% group_by(OTU) %>% summarize(peak=mean(peakSal_withzero), minSal=min(startSal), maxSal=max(endSal)) %>%
        group_by(OTU) %>%
        mutate(midSal=mean(c(minSal, maxSal))
               , lower = midSal<midBrackThresh
               , uppr = midSal>=midBrackThresh) 
    brackishLower <- brackishSplit%>%
        filter(lower) %>% select(OTU) %>% pull()
    brackishUpper<- brackishSplit%>%
        filter(uppr) %>% select(OTU) %>% pull()
    # lower.brackish.taxa <- mb[mb$taxa %in% brackishLower,]
    # upper.brackish.taxa <- mb[mb$taxa %in% brackishUpper,]
    
    # Get midpoint of range for all brackish
    # also, separate
    temp.boundaries <- mb[match(brackish.taxa,mb[,1]),c("taxa","boundaries","boundariestwo")]
    temp.midpoints <- c(-100,apply(temp.boundaries[,c(2,3)], MARGIN = 1, FUN = mean)) # -100 is to double check this is removed
    names(temp.midpoints) <- c("CON",temp.boundaries[,1])
    high.brackish.taxa <- names(which(temp.midpoints >= 15))
    low.brackish.taxa <- names(which(temp.midpoints < 15))
    
    ######### Get PD for nearest fresh, marine, brack for each OTU #########
    typesSimple <- c("freshRestricted","marineRestricted","brackishRestricted")
    
    # NOTE: I always start with a 2 so that it always has something to remove; makes it easier to go between
    # 16S and 18S
    freshOTUPD <- list()
    for (t in typesSimple) {
        # t <- typesSimple[1]
        freshOTUPD[[t]] <- 2
        for (n in 1:length(fresh.taxa)) {
            # n <- 1
            otuTestList <- as.character(mb[mb[,3]==t,1])
            otu <- as.character(fresh.taxa[n])
            # if it's fresh and fresh need to get rid of OTU that is the same as it
            if ( otu %in% otuTestList ) {
                otuTestList <- otuTestList[-match(otu, otuTestList)]
            }
            freshOTUPD[[t]] <- c(freshOTUPD[[t]],min(taxa.dist[otu,otuTestList]))
            # match(otuTestList, colnames(taxa.dist))
        }
    }
    
    marineOTUPD <- list()
    for (t in typesSimple) {
        marineOTUPD[[t]] <- 2
        for (n in 1:length(marine.taxa)) {
            otuTestList <- as.character(mb[mb[,3]==t,1])
            otu <- as.character(marine.taxa[n])
            # if it's fresh and fresh need to get rid of OTU that is the same as it
            if ( otu %in% otuTestList ) {
                otuTestList <- otuTestList[-match(otu, otuTestList)]
            }
            marineOTUPD[[t]] <- c(marineOTUPD[[t]],min(taxa.dist[otu,otuTestList]))
        }
    }
    
    brackishOTUPD <- list()
    for (t in typesSimple) {
        brackishOTUPD[[t]] <- 2
        for (n in 1:length(brackish.taxa)) {
            otuTestList <- as.character(mb[mb[,3]==t,1])
            otu <- as.character(brackish.taxa[n])
            # if it's fresh and fresh need to get rid of OTU that is the same as it
            if ( otu %in% otuTestList ) {
                otuTestList <- otuTestList[-match(otu, otuTestList)]
            }
            brackishOTUPD[[t]] <- c(brackishOTUPD[[t]],min(taxa.dist[otu,otuTestList]))
        }
    }
    
    # For each analysis, REMOVE pairs where one of them has a distance greater than 1
    
    # BRACKISH statistical comparisons
    remove_temp_brack <- which(apply(cbind(brackishOTUPD$freshRestricted, brackishOTUPD$marineRestricted), MARGIN=1, FUN = function(x) any(x>1)))
    capture.output(t.test(brackishOTUPD$freshRestricted[-remove_temp_brack], brackishOTUPD$marineRestricted[-remove_temp_brack], paired=TRUE),file = paste0(output,"/brackish_compare_t.test_paired.txt"))
    capture.output(shapiro.test(brackishOTUPD$freshRestricted[-remove_temp_brack]- brackishOTUPD$marineRestricted[-remove_temp_brack]),file=paste0(output,"/brackish_compare_shapiro.txt"))
    capture.output(wilcox.test(brackishOTUPD$freshRestricted[-remove_temp_brack], brackishOTUPD$marineRestricted[-remove_temp_brack]), file=paste0(output,"/brackish_compare_wilcox.txt"))
    
    # FRESH statistical comparisons
    remove_temp_fresh <- which(apply(cbind(freshOTUPD$freshRestricted, freshOTUPD$marineRestricted), MARGIN=1, FUN = function(x) any(x>1)))
    capture.output(t.test(freshOTUPD$freshRestricted[-remove_temp_fresh], freshOTUPD$marineRestricted[-remove_temp_fresh], paired=TRUE),file = paste0(output,"/fresh_compare_t.test_paired.txt"))
    capture.output(shapiro.test(freshOTUPD$freshRestricted[-remove_temp_fresh]-freshOTUPD$marineRestricted[-remove_temp_fresh]),file=paste0(output,"/fresh_compare_shapiro.txt"))
    capture.output(wilcox.test(freshOTUPD$freshRestricted[-remove_temp_fresh], freshOTUPD$marineRestricted[-remove_temp_fresh]), file=paste0(output,"/fresh_compare_wilcox.txt"))
    
    # MARINE statistical comparisons
    remove_temp_marine <- which(apply(cbind(marineOTUPD$freshRestricted, marineOTUPD$marineRestricted), MARGIN=1, FUN = function(x) any(x>1)))
    capture.output(t.test(marineOTUPD$freshRestricted[-remove_temp_marine], marineOTUPD$marineRestricted[-remove_temp_marine], paired=TRUE),file = paste0(output,"/marine_compare_t.test_paired.txt"))
    capture.output(shapiro.test(marineOTUPD$freshRestricted[-remove_temp_marine]-marineOTUPD$marineRestricted[-remove_temp_marine]),file=paste0(output,"/marine_compare_shapiro.txt"))
    capture.output(wilcox.test(marineOTUPD$freshRestricted[-remove_temp_marine], marineOTUPD$marineRestricted[-remove_temp_marine]), file=paste0(output,"/marine_compare_wilcox.txt"))
    
    ### First, plot fresh and marine OTUs in relation to each other
    # find limits of fresh and  marine OTUs
    rangeAxis <- range(c(freshOTUPD$freshRestricted[-remove_temp_fresh],marineOTUPD$freshRestricted[-remove_temp_marine],freshOTUPD$marineRestricted[-remove_temp_fresh],marineOTUPD$marineRestricted[-remove_temp_marine]))
    pdf(paste0(output,"/marinevsfresh_PD_to_nearestMarineFresh.pdf"))
    plot(NULL, xlim=c(0,rangeAxis[2]), ylim=c(0,rangeAxis[2]), xlab="Distance to nearest fresh OTU", ylab="Distance to nearest marine OTU")
    points(freshOTUPD$freshRestricted[-remove_temp_fresh], freshOTUPD$marineRestricted[-remove_temp_fresh], pch=19, col="blue")
    points(marineOTUPD$freshRestricted[-remove_temp_marine], marineOTUPD$marineRestricted[-remove_temp_marine], pch=19, col="red")
    abline(a=0,b=1)
    legend("top", legend=c("Fresh OTUs","Marine OTUs"), pch=c(19,19), col=c("blue","red"))
    dev.off()

    maxAxis <- max(c(brackishOTUPD$freshRestricted[-remove_temp_brack],brackishOTUPD$marineRestricted[-remove_temp_brack]))
    pdf(paste0(output,"/brackish_PD_to_nearestMarineFresh.pdf"))
    plot(brackishOTUPD$freshRestricted[-remove_temp_brack], brackishOTUPD$marineRestricted[-remove_temp_brack], pch=19, col="purple"
         , xlab="Distance to nearest fresh OTU"
         , ylab="Distance to nearest marine OTU"
         , xlim=c(0,maxAxis)
         , ylim=c(0,maxAxis))
    abline(a=0,b=1)
    dev.off()
    
    pdf(paste0(output,"/brackish_PD_correlation_to_midpoint.pdf"), 10,5)
    par(mfrow=c(1,2))
    plot(brackishOTUPD$freshRestricted[-remove_temp_brack]~as.numeric(temp.midpoints)[-remove_temp_brack], xlab="Midpoint of tolerance range", ylab="Phyl. Dist. to nearest freshwater specialist", col="blue")
    # cor.to.fresh.summary <- summary(lm(brackishOTUPD$freshRestricted~as.numeric(temp.midpoints)))
    cor.to.fresh.summary = cor.test(as.numeric(temp.midpoints)[-remove_temp_brack],brackishOTUPD$freshRestricted[-remove_temp_brack], method=c("spearman"))
    legend("topright", legend=c(paste0("P-value = ",signif(cor.to.fresh.summary$p.value,2)), paste0("Spearmen's Rho = ",signif(as.numeric(cor.to.fresh.summary$estimate),2))), pch=c("",""), bty="n")
    abline(lm(brackishOTUPD$freshRestricted[-remove_temp_brack]~as.numeric(temp.midpoints)[-remove_temp_brack]))
    
    plot(brackishOTUPD$marineRestricted[-remove_temp_brack]~as.numeric(temp.midpoints)[-remove_temp_brack], xlab="Midpoint of tolerance range", ylab="Phyl. Dist. to nearest marine specialist", col = "red")
    # cor.to.marine.summary <- summary(lm(brackishOTUPD$marineRestricted~as.numeric(temp.midpoints)))
    cor.to.marine.summary = cor.test(as.numeric(temp.midpoints)[-remove_temp_brack],brackishOTUPD$marineRestricted[-remove_temp_brack], method=c("spearman"))
    legend("topright", legend=c(paste0("P-value = ",signif(cor.to.marine.summary$p.value,2)), paste0("Spearman's Rho = ",signif(as.numeric(cor.to.marine.summary$estimate),2))), pch=c("",""), bty="n")
    abline(lm(brackishOTUPD$marineRestricted[-remove_temp_brack]~as.numeric(temp.midpoints)[-remove_temp_brack]))
    dev.off()

    pdf(paste0(output,"/brackish_PD_correlation_to_midpoint_combined.pdf"))
    plot(brackishOTUPD$freshRestricted[-remove_temp_brack]-brackishOTUPD$marineRestricted[-remove_temp_brack]~as.numeric(temp.midpoints)[-remove_temp_brack], xlab="Midpoint of tolerance range", ylab="Phyl. Dist. to nearest freshwater specialist")
    # cor.to.fresh.summary <- summary(lm(brackishOTUPD$freshRestricted~as.numeric(temp.midpoints)))
    cor.to.freshmarine.summary = cor.test(as.numeric(temp.midpoints)[-remove_temp_brack],brackishOTUPD$freshRestricted[-remove_temp_brack]-brackishOTUPD$marineRestricted[-remove_temp_brack], method=c("spearman"))
    legend("topright", legend=c(paste0("P-value = ",signif(cor.to.fresh.summary$p.value,2)), paste0("Spearmen's Rho = ",signif(as.numeric(cor.to.freshmarine.summary$estimate),2))), pch=c("",""), bty="n")
    abline(lm(brackishOTUPD$freshRestricted[-remove_temp_brack]-brackishOTUPD$marineRestricted[-remove_temp_brack]~as.numeric(temp.midpoints)[-remove_temp_brack]))
    dev.off()
    
    brack.fresh <- lm(as.numeric(temp.midpoints)[-remove_temp_brack]~brackishOTUPD$freshRestricted[-remove_temp_brack])
    capture.output(stepAIC(brack.fresh),file=paste0(output,"/AIC_brackfresh.txt"))
    brack.marine <- lm(as.numeric(temp.midpoints)[-remove_temp_brack]~brackishOTUPD$marineRestricted[-remove_temp_brack])
    capture.output(stepAIC(brack.marine),file=paste0(output,"/AIC_brackmarine.txt"))
    
    ## Compare lower brack and high brack
    
    brackishTypes <- c("brackishUpper","brackishLower","fresh.taxa","marine.taxa")

    brackishLowerPD <- list()
    for (t in brackishTypes) {
        brackishLowerPD[[t]] <- 2
        for (n in 1:length(brackishLower)) {
            otuTestList <- as.character(get(t))
            
            otu <- as.character(brackishLower[n])
            # if it's fresh and fresh need to get rid of OTU that is the same as it
            if ( otu %in% otuTestList ) {
                otuTestList <- otuTestList[-match(otu, otuTestList)]
            }
            brackishLowerPD[[t]] <- c(brackishLowerPD[[t]],min(taxa.dist[otu,otuTestList]))
        }
    }
   
    brackishUpperPD <- list()
    for (t in brackishTypes) {
        brackishUpperPD[[t]] <- 2
        for (n in 1:length(brackishUpper)) {
            otuTestList <- as.character(get(t))
            
            otu <- as.character(brackishUpper[n])
            # if it's fresh and fresh need to get rid of OTU that is the same as it
            if ( otu %in% otuTestList ) {
                otuTestList <- otuTestList[-match(otu, otuTestList)]
            }
            brackishUpperPD[[t]] <- c(brackishUpperPD[[t]],min(taxa.dist[otu,otuTestList]))
        }
    }
    
    # Get rid of points with PD more than 1; probably not resolved well.
    anyMore1Lower <- unique(unlist(lapply(brackishLowerPD, function(x) {which(x>1)})))
    brackishLowerPD <- lapply(brackishLowerPD, function(x) {x[-anyMore1Lower]})
    anyMore1Upper <- unique(unlist(lapply(brackishUpperPD, function(x) {which(x>1)})))
    brackishUpperPD <- lapply(brackishUpperPD, function(x) {x[-anyMore1Upper]})
    
    rangeLimBrack <- range(c(unlist(brackishLowerPD), unlist(brackishUpperPD)))
    pdf(paste0(output,"/brackish_upperlowerneighbours.pdf"), height=5, width=6)
    plot(brackishLowerPD$brackishUpper, brackishLowerPD$brackishLower
         , xlim=rangeLimBrack, ylim=rangeLimBrack
         , xlab="PD from Nearest Upper Brackish Neighbour"
         , ylab="PD from Nearest Lower Brackish Neighbour", col=alpha("darkslateblue",0.5) , pch=19)
    points(brackishUpperPD$brackishUpper, brackishUpperPD$brackishLower, col=alpha("magenta",0.5), pch=19)
    abline(a=0, b=1, col="black", lwd=2)
    legend("topleft", legend=c("Lower brackish OTUs","Upper brackish OTUs"), pch=19, col=c(alpha("darkslateblue",0.5),alpha("magenta",0.5)), bty = "n")
    dev.off()
    
    pdf(paste0(output,"/brackish_lowerfreshneighbours.pdf"), height=5, width=6)
    plot(brackishLowerPD$fresh.taxa, brackishLowerPD$brackishLower
         , xlim=rangeLimBrack, ylim=rangeLimBrack
         , xlab="PD from Nearest Fresh Neighbour"
         , ylab="PD from Nearest Lower Brackish Neighbour", col="blue" , pch=19)
    abline(a=0, b=1, col="black", lwd=2)
    dev.off()
    pdf(paste0(output,"/brackish_lowermarineneighbours.pdf"), height=5, width=6)
    plot(brackishLowerPD$marine.taxa, brackishLowerPD$brackishLower
         , xlim=rangeLimBrack, ylim=rangeLimBrack
         , xlab="PD from Nearest Marine Neighbour"
         , ylab="PD from Nearest Lower Brackish Neighbour", col="red" , pch=19)
    abline(a=0, b=1, col="black", lwd=2)
    dev.off()
    
    
    pdf(paste0(output,"/brackish_uppermarineneighbours.pdf"), height=5, width=6)
    plot(brackishLowerPD$brackishUpper, brackishLowerPD$marine.taxa
         , xlim=rangeLimBrack, ylim=rangeLimBrack
         , xlab="PD from Nearest Upper Brackish Neighbour"
         , ylab="PD from Nearest Marine Neighbour", col="red" , pch=19)
    abline(a=0, b=1, col="black", lwd=2)
    dev.off()
    
    pdf(paste0(output,"/brackish_upperfreshneighbours.pdf"), height=5, width=6)
    plot(brackishLowerPD$brackishUpper, brackishLowerPD$fresh.taxa
         , xlim=rangeLimBrack, ylim=rangeLimBrack
         , xlab="PD from Nearest Upper Brackish Neighbour"
         , ylab="PD from Nearest Fresh Neighbour", col="blue", pch=19)
    abline(a=0, b=1, col="black", lwd=2)
    dev.off()
    
    
    pdf(paste0(output,"/brackish_lower_freshorupperneighbours.pdf"), height=5, width=6)
    plot(brackishLowerPD$fresh.taxa, brackishLowerPD$brackishUpper
         , xlim=rangeLimBrack, ylim=rangeLimBrack
         , xlab="PD from Nearest Fresh Neighbour"
         , ylab="PD from Nearest Upper Brackish Neighbour", col="grey" , pch=19)
    abline(a=0, b=1, col="black", lwd=2)
    dev.off()
    pdf(paste0(output,"/brackish_upper_marineorlowerneighbours.pdf"), height=5, width=6)
    plot(brackishUpperPD$marine.taxa, brackishUpperPD$brackishLower
         , xlim=rangeLimBrack, ylim=rangeLimBrack
         , xlab="PD from Nearest Marine Neighbour"
         , ylab="PD from Nearest Lower Brackish Neighbour", col="grey" , pch=19)
    abline(a=0, b=1, col="black", lwd=2)
    dev.off()
    
    ### Try boxplots
    # reorganize
    brackishPD_dat <- rbind(data.frame(PD=brackishUpperPD$brackishUpper, compare="UpperBrack", OTU_type="UpperBrack")
                                 ,data.frame(PD=brackishUpperPD$brackishLower, compare="LowerBrack", OTU_type="UpperBrack")
                                 ,data.frame(PD=brackishUpperPD$fresh.taxa, compare="Fresh", OTU_type="UpperBrack")
                                 ,data.frame(PD=brackishUpperPD$marine.taxa, compare="Marine", OTU_type="UpperBrack")
                                 ,data.frame(PD=brackishLowerPD$brackishUpper, compare="UpperBrack", OTU_type="LowerBrack")
                                 ,data.frame(PD=brackishLowerPD$brackishLower, compare="LowerBrack", OTU_type="LowerBrack")
                                 ,data.frame(PD=brackishLowerPD$fresh.taxa, compare="Fresh", OTU_type="LowerBrack")
                                 ,data.frame(PD=brackishLowerPD$marine.taxa, compare="Marine", OTU_type="LowerBrack")
    ) %>%
        mutate(compare=factor(compare, levels=c("Fresh","LowerBrack","UpperBrack","Marine")))
    
    ggsave(paste0(output,"/brackish_all_boxplot.pdf")
           , brackishPD_dat %>% ggplot(aes(x=compare, y=sqrt(PD))) + geom_boxplot() +
               facet_grid(OTU_type~.)+
               xlab("Compared to")+ylab("Phylogenetic Distance (sqrt)"))
   ## Statistics
    # Are brackish upper/lower closer to themselves, or to others?
    capture.output(wilcox.test(brackishUpperPD$brackishUpper, brackishUpperPD$brackishLower, paired = TRUE)
                   ,file=paste0(output,"/brackish_wilcoxpaired_upper_uppervslower.txt"))
    capture.output(wilcox.test(brackishLowerPD$brackishUpper, brackishLowerPD$brackishLower, paired = TRUE)
                   ,file=paste0(output,"/brackish_wilcoxpaired_lower_uppervslower.txt"))
    temp.onlybrack_upper <- brackishPD_dat %>% filter(compare%in% c("LowerBrack","UpperBrack"), OTU_type=="UpperBrack")
    temp.onlybrack_lower <- brackishPD_dat %>% filter(compare%in% c("LowerBrack","UpperBrack"), OTU_type=="LowerBrack")
    capture.output(stepAIC(lm(PD ~ compare, data=temp.onlybrack_lower)),file=paste0(output,"/brackish_AIC_lower_uppervslower.txt"))
    capture.output(stepAIC(lm(PD ~ compare, data=temp.onlybrack_upper)),file=paste0(output,"/brackish_AIC_upper_uppervslower.txt"))
    
    # Are brackish upper/lower closer to other brackish, or to fresh/marine?
    
    capture.output(wilcox.test(brackishUpperPD$brackishLower, brackishUpperPD$marine.taxa, paired = TRUE)
                   , file=paste0(output,"/brackish_wilcoxpaired_uppervslowerormarine.txt"))
    capture.output(wilcox.test(brackishLowerPD$brackishUpper, brackishLowerPD$fresh.taxa, paired = TRUE)
                   , file=paste0(output,"/brackish_wilcoxpaired_lowervsupperorfresh.txt"))
    
    
}

