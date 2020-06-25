#!bin/bash/Rscript

library("phyloseq")# For unifrac
library("MASS") # for beta div plotting
library("picante") # For unifrac
library("ape") # for tree?
library("vegan") # for vegdist
library("mgcv") # for fitting splines
library("tidyverse")
### Testing for seeing if lower/upper brackish are phylogenetically distinct

### FUNCTIONS ####
source("./downstream_rcode/00_sourcecode.R")

### Loading dirNames ###
dirNames <- unlist(read.delim("./OUTPUT/dirNames.txt", header=FALSE))

#### START ####

newDirName <- "./OUTPUT/03b_ToleranceRangePlots"
dir.create(newDirName)
for ( d in c("16sBaltic","16sFraser","18sFraser")) {
    dshort <- ifelse(d=="16sBaltic","B16",ifelse(d=="16sFraser","F16","F18"))
    temppwd <- paste0(dirNames[match(d, c("16sBaltic","16sFraser","18sFraser"))],"/")
    TR <- read.delim(paste0(temppwd,"toleranceRanges.txt"))
    mb <- read.delim(paste0(temppwd,"modelBoundaries_type.txt"))
    treePWD <- ifelse(d=="16sBaltic","B16",ifelse(d=="16sFraser","F16","F18"))
    tr <- read.tree(paste0("./OUTPUT/tree_",treePWD,"_filt.tre"))
    tax <- read.delim(paste0(temppwd,"taxaIDLegend.txt"), row.names = 1, header = FALSE)
    otu <- read.delim(paste0(temppwd, "OTUTableText.txt"))
    mfPWD <- ifelse(d=="16sBaltic","./INPUT_DATA/Baltic_16S/metadata_table.tsv"
                    ,ifelse(d=="16sFraser","./INPUT_DATA/Fraser_16S/MF_16sFraser_noConCOL.txt"
                            ,"./INPUT_DATA/Fraser_18S/MF_18sFraser_noConCOL.txt"))
    mf <- read.delim(paste0(mfPWD))
    if ( d %in% c("16sFraser", "18sFraser")) {
        mf$X.SampleID <- gsub("-",".",mf$X.SampleID, fixed = TRUE)
    }
    TR$OTU <- as.character(TR$OTU)
    
    # Filter OTU tables
    otu.adj <- otu %>% select(-c(X.OTUID,taxonomy))
    rownames(otu.adj) <- otu.adj$X.OTUID
    # Filter mf and otu
    mf <- mf[mf$X.SampleID %in% colnames(otu),]
    otu.adj <- otu.adj[,colnames(otu.adj) %in% mf$X.SampleID]
    # Make into relative abundance
    otu.adj.relAbund <- apply(otu.adj, 2, FUN=function(c){c/sum(c)})
    
    ## Add extra column of descriptions
    mb.adj <- mb %>% 
        mutate(OTU=as.character(taxa)) %>%
        filter(typeSimple!="noclass") %>%
        mutate(boundaries=ifelse(is.na(boundaries), minGrad,boundaries)
               , boundariestwo=ifelse(is.na(boundariestwo), maxGrad,boundariestwo)
               , midBoundary=(boundaries+boundariestwo)/2
               , breadth=abs(c(boundaries-boundariestwo))) %>%
        select(OTU, typeSimple,type,boundaries, boundariestwo, midBoundary, breadth)
    
    tipDist <- cophenetic.phylo(tr)
    try <- length(mb.adj$taxa)
    
    # Nearest Neighbour ratio- marine/fresh only
    if ( !file.exists(paste0(newDirName,"/",d,"_nnR.RData")) ) {
        allnnR <- c()
        allnnR_brack <- c()
        total <- length(mb.adj$OTU)
        pb <- txtProgressBar(min=0, max=total, style = 3)
        s <- 0
        for ( t in mb.adj$OTU ) {
            tempDist <- tipDist[t,]
            tempDist <- tempDist[-which(t==names(tempDist))]
            nnR <- data.frame(distance=tempDist,ty=mb[match(names(tempDist), mb$taxa),"typeSimple"]) %>%
                as_tibble() %>% 
                group_by(ty) %>%
                arrange(ty,distance) %>%
                mutate(rank=1:n()) %>%
                ungroup() %>%
                filter(ty %in% c("freshRestricted","marineRestricted"), rank %in% 1) %>% 
                spread(key=ty, value=distance) %>%
                mutate(nnR=marineRestricted/mean(c(marineRestricted,freshRestricted))-1) %>%
                pull(nnR)
            nnR_brack <- data.frame(distance=tempDist,ty=mb[match(names(tempDist), mb$taxa),"typeSimple"]) %>%
                as_tibble() %>% 
                group_by(ty) %>%
                arrange(ty,distance) %>%
                mutate(rank=1:n()) %>%
                ungroup() %>%
                filter(ty %in% c("freshRestricted","marineRestricted", "brackishRestricted"), rank %in% 1) %>%
                spread(key=ty, value=distance) %>%
                mutate(maxFM = max(c(freshRestricted, marineRestricted), na.rm = FALSE)
                    , nnR_brack=brackishRestricted/(sum(c(brackishRestricted,maxFM)))) %>%
                pull(nnR_brack)
            allnnR <- c(allnnR, nnR)
            allnnR_brack <- c(allnnR_brack, nnR_brack)
            s <- s+1
            Sys.sleep(0.1)
            setTxtProgressBar(pb, s)
        }
        close(pb)
        save(allnnR, file=paste0(newDirName,"/", d,"_nnR.RData"))
        save(allnnR_brack, file=paste0(newDirName,"/",d,"_nnR_brack.RData"))
        
    } else {
        load(paste0(newDirName,"/",d,"_nnR.RData"))
        load(paste0(newDirName,"/",d,"_nnR_brack.RData"))
    }
    mb.adj$nnR <- allnnR
    mb.adj$nnR_brack <- allnnR_brack
    
    # Nearest Neighbour ratio- marine/fresh/brack
    if ( !file.exists(paste0(newDirName,"/",d,"_nnR_wbrack.RData")) ) {
        allnnR_rgb <- c()
        total <- length(mb.adj$OTU)
        pb <- txtProgressBar(min=0, max=total, style = 3)
        s <- 0
        for ( t in mb.adj$OTU ) {
            # t <- mb.adj$OTU[1]
            tempDist <- tipDist[t,]
            tempDist <- tempDist[-which(t==names(tempDist))]
            nnR_rgb <- data.frame(distance=tempDist,ty=mb[match(names(tempDist), mb$taxa),"typeSimple"]) %>%
                as_tibble() %>% 
                group_by(ty) %>%
                arrange(ty,distance) %>%
                mutate(rank=1:n()) %>%
                ungroup() %>%
                filter(ty %in% c("freshRestricted","marineRestricted", "brackishRestricted"), rank %in% 1) %>% 
                spread(key=ty, value=distance) %>%
                mutate(maxDist=max(c(brackishRestricted, freshRestricted, marineRestricted), na.rm = TRUE)
                       , minDist=min(c(brackishRestricted, freshRestricted, marineRestricted), na.rm = TRUE)
                       , diffDist = maxDist-minDist
                       , green = (diffDist-(brackishRestricted-minDist))/diffDist
                       , red = (diffDist-(marineRestricted-minDist))/diffDist
                       , blue = (diffDist-(freshRestricted-minDist))/diffDist) %>%  
                select(red, green, blue) %>% unlist()
            allnnR_rgb <- rbind(allnnR_rgb, nnR_rgb)
            s <- s+1
            Sys.sleep(0.1)
            setTxtProgressBar(pb, s)
        }
        close(pb)
        save(allnnR_rgb, file=paste0(newDirName,"/",d,"_nnR_wbrack.RData"))
    } else {
        load(paste0(newDirName,"/",d,"_nnR_wbrack.RData"))
    }
    
    mb.adj <- cbind(mb.adj, allnnR_rgb) %>%
        filter(!is.na(red)) %>%
        mutate(rgb_col = rgb(red,green,blue))
    #### Get beta div
    
    if ( !file.exists(paste0(newDirName,"/dm_bray.",dshort,".RData")) ) {
        temp.bray <- as.matrix(vegdist(t(otu.adj.relAbund), method = "bray"))
        
        save(temp.bray, file = paste0(newDirName,"/","dm_bray.",dshort,".RData"))
        
    } else {
        load(paste0(newDirName,"/","dm_bray.",dshort,".RData"))
    }
    
    #### Getting jumps ####
    mf <- mf %>%
        filter(!is.na(as.numeric(SalinityEnviron))) %>%
        mutate(SalRound = as.numeric(SalinityEnviron))
    SalRound <- sort(unique(mf$SalRound))
    
    if ( !file.exists(paste0(newDirName,"/jumps_",dshort,",RData"))) {
        jumps <- data.frame()
        total <- (length(SalRound)-1)
        # create progress bar
        pb <- txtProgressBar(min = 0, max = total, style = 3)
        for ( s in 1:(length(SalRound)-1) ) {
            sstart <- SalRound[s]
            send <- SalRound[s+1]
            starting <- mf[mf$SalRound==sstart,"X.SampleID"]
            ending <- mf[mf$SalRound== send,"X.SampleID"]
            
            # distUWUnifrac <- as.vector(temp.uwunifrac[starting,ending])
            distBray <- as.vector(temp.bray[starting,ending])
            # distJaccard <- as.vector(temp.jaccard[starting,ending])
            core <- coreKept(otu.adj.relAbund, starting, ending)
            turnover <- OTUturnover(otu.adj.relAbund, starting, ending)
            
            jumps <- rbind(jumps, data.frame(sstart,send
                                             # ,distUWUnifrac
                                             ,distBray
                                             # ,distJaccard
                                             ,core, gained=turnover[[1]], lost=turnover[[2]], totalRich=turnover[[3]]))
            
            Sys.sleep(0.1)
            setTxtProgressBar(pb, s)
        }
        close(pb)
        save(jumps, file = paste0(newDirName,"/","jumps_",dshort,".RData")) 
    } else {
        load(paste0(newDirName,"/","jumps_",dshort,".RData"))
    }
    
    ### Fitting GAMs ####
    aveSal <- rowMeans(cbind(jumps$send, jumps$sstart))
    
    gam_bray <- gam(distBray ~ s(aveSal), data=jumps, family = gaussian)
    gam_core <- gam(core ~ s(aveSal), data=jumps, family = gaussian)
    gam_gained <- gam(gained ~ s(aveSal), data=jumps, family = gaussian)
    gam_lost <- gam(lost ~ s(aveSal), data=jumps, family = gaussian)
    
    ggsave(filename = paste0(newDirName,"/",dshort,"turnover.pdf")
           ,data.frame(Bray_curtis_distance=gam_bray$fitted.values, Core_unchanged=gam_core$fitted.values, Salinity=aveSal) %>%
               select(Salinity, Bray_curtis_distance, Core_unchanged) %>%
               gather(-Salinity, key=Metric, value=Change) %>%
               ggplot() +         
               geom_line(aes(x=Salinity, y=Change, col=Metric)) +
               scale_color_manual(values=c("blue","orange"))+
               # geom_line(aes(x=Salinity, y=core), col="orange") +
               ylab("Bray Curtis Distance") +
               scale_y_continuous(sec.axis=sec_axis(
                   trans = ~.*1
                   , name="Percent community unchanged"
               ))  +
               theme(
                   axis.title.y = element_text(color = "blue")
                   ,axis.title.y.right = element_text(color = "orange")) +
               geom_point(data=data.frame(raw=jumps$distBray, Salinity=aveSal), aes(x=Salinity, y=raw),col="blue")+
               geom_point(data=data.frame(raw=jumps$core, Salinity=aveSal), aes(x=Salinity, y=raw),col="orange")
    )
    
    ggsave(filename = paste0(newDirName,"/",dshort,"gained_lost.pdf")
           ,data.frame(lost=gam_lost$fitted.values, gained=gam_gained$fitted.values, Salinity=aveSal) %>%
               gather(-Salinity, key=Metric, value=Change) %>%
               ggplot() +         
               geom_line(aes(x=Salinity, y=Change, col=Metric)) +
               ylab("Number of OTUs lost/gained") 
    )
    
    ggsave(filename = paste0(newDirName,"/",dshort,"turnover_all.pdf")
           ,data.frame(lost=gam_lost$fitted.values/max(gam_lost$fitted.values), gained=gam_gained$fitted.values/max(gam_gained$fitted.values), bray=gam_bray$fitted.values, core=gam_core$fitted.values, Salinity=aveSal) %>%
               gather(-Salinity, key=Metric, value=Change) %>%
               ggplot() +         
               geom_line(aes(x=Salinity, y=Change, col=Metric)) 
    )
    
    ###### Minima/maxima for change ####
    if (d=="16sBaltic") {
        sink(file=paste0(newDirName,"/maxima.txt"))
    } else {
        sink(file=paste0(newDirName,"/maxima.txt", append=TRUE))
    }
    print(d)
    
    predictBray <- predict(gam_bray, newdata = data.frame(aveSal=unique(aveSal)))
    localMax_bray <- unique(aveSal)[localMax(predictBray)]
    print("local max, bray")
    print(localMax_bray)
    
    predictGained <- predict(gam_gained, newdata = data.frame(aveSal=unique(aveSal)))
    localMax_gained <- unique(aveSal)[localMax(predictGained)]
    print("local max, gained")
    print(localMax_gained)
    
    predictLost <- predict(gam_lost, newdata = data.frame(aveSal=unique(aveSal)))
    localMax_lost <- unique(aveSal)[localMax(predictLost)]
    print("local max, lost")
    print(localMax_lost)
    sink()
    #######
    combinedTR_mb <- TR %>%
        left_join(mb.adj)%>%
        filter(!is.na(typeSimple), abundances_nozero !=0, position=="withinBounds") %>%
        mutate(typeSimple=gsub("Restricted","",typeSimple)
               ,typeSimple=gsub("^","QTAG-assigned ",typeSimple)) %>%
        mutate(typeSimple=factor(typeSimple, levels=c("QTAG-assigned fresh","QTAG-assigned brackish","QTAG-assigned marine"))) %>%
        arrange(typeSimple, midBoundary) %>%
        mutate(OTU=factor(OTU, levels=unique(OTU)), typeSimple=factor(typeSimple, levels=c("QTAG-assigned marine","QTAG-assigned brackish","QTAG-assigned fresh")))
    
    nOTUs <- combinedTR_mb %>% select(OTU)%>%pull() %>%unique() %>% length()
    outliers_points <- TR %>%
        left_join(mb.adj) %>%
        filter(!is.na(typeSimple), abundances_nozero !=0, position=="outlier") %>%
        mutate(typeSimple=gsub("Restricted","",typeSimple)
               ,typeSimple=gsub("^","QTAG-assigned ",typeSimple)) %>%
        mutate(typeSimple=factor(typeSimple, levels=c("QTAG-assigned fresh","QTAG-assigned brackish","QTAG-assigned marine"))) %>%
        arrange(typeSimple, midBoundary) %>%
        mutate(OTU=factor(OTU, levels=unique(OTU)))
    lwd_temp <- ifelse(d=="16sBaltic",1,ifelse(d=="16sFraser",0.0002, 0.0001))
    ggsave(paste0(newDirName,"/",d,"_nnR_tolerancerange.pdf"), height=10, width=6.5
           ,combinedTR_mb %>%
               ggplot() + geom_segment(aes(x=startSal, xend=endSal, y=OTU, yend=OTU, col=nnR, alpha=abundances_nozero_maxstand), lwd= lwd_temp #100/nOTUs*1.2)
                                       ) +
               geom_point(data=outliers_points,aes(x=startSal, y=OTU,col=nnR, alpha=abundances_nozero_maxstand), cex=lwd_temp )+
               scale_color_gradient2(low="red",mid="darkgrey",high="blue")+
               labs(alpha = paste0("Relative abundance \n(standardized by max reads)"), col="Nearest Neighbour Ratio") +
               theme_classic() +
               theme(axis.text.y = element_blank(), axis.ticks.y=element_blank()) + xlab("Salinity") +
               geom_vline(data=data.frame(TurnoverPoints=localMax_bray), mapping=aes(xintercept=TurnoverPoints), lty=1, lwd=2, col="black", alpha=0.5) +
               facet_grid(typeSimple~. ,scales="free", space = "free_y")
           
    )
    ggsave(paste0(newDirName,"/",d,"_nnR_tolerancerange_darker.pdf"), height=10, width=6.5
           ,combinedTR_mb %>%
               ggplot() + geom_segment(aes(x=startSal, xend=endSal, y=OTU, yend=OTU, col=nnR, alpha=abundances_nozero_maxstand), lwd= lwd_temp #100/nOTUs*1.2)
               ) +
               geom_point(data=outliers_points,aes(x=startSal, y=OTU,col=nnR, alpha=abundances_nozero_maxstand), cex=lwd_temp )+
               scale_color_gradient2(low="red",mid="black",high="blue")+
               labs(alpha = paste0("Relative abundance \n(standardized by max reads)"), col="Nearest Neighbour Ratio") +
               theme_classic() +
               theme(axis.text.y = element_blank(), axis.ticks.y=element_blank()) + xlab("Salinity") +
               geom_vline(data=data.frame(TurnoverPoints=localMax_bray), mapping=aes(xintercept=TurnoverPoints), lty=1, lwd=2, col="black", alpha=0.5) +
               facet_grid(typeSimple~. ,scales="free", space = "free_y")
           
    )
    ggsave(paste0(newDirName,"/",d,"_nnR_tolerancerange_contrastx5.pdf"), height=10, width=6.5
           ,combinedTR_mb %>%
               ggplot() + geom_segment(aes(x=startSal, xend=endSal, y=OTU, yend=OTU, col=nnR, alpha=abundances_nozero_maxstand^5), lwd= lwd_temp #100/nOTUs*1.2)
               ) +
               geom_point(data=outliers_points,aes(x=startSal, y=OTU,col=nnR, alpha=abundances_nozero_maxstand^5), cex=lwd_temp )+
               scale_color_gradient2(low="red",mid="darkgrey",high="blue")+
               labs(alpha = paste0("Relative abundance \n(standardized by max reads)"), col="Nearest Neighbour Ratio") +
               theme_classic() +
               theme(axis.text.y = element_blank(), axis.ticks.y=element_blank()) + xlab("Salinity") +
               geom_vline(data=data.frame(TurnoverPoints=localMax_bray), mapping=aes(xintercept=TurnoverPoints), lty=1, lwd=2, col="black", alpha=0.5) +
               facet_grid(typeSimple~. ,scales="free", space = "free_y")
           
    )
    ggsave(paste0(newDirName,"/",d,"_nnR_tolerancerange_noabundance.pdf"), height=10, width=6.5
           ,combinedTR_mb %>%
               ggplot() + geom_segment(aes(x=startSal, xend=endSal, y=OTU, yend=OTU, col=nnR), lwd= lwd_temp #100/nOTUs*1.2)
               ) +
               geom_point(data=outliers_points,aes(x=startSal, y=OTU,col=nnR), cex=lwd_temp )+
               scale_color_gradient2(low="red",mid="darkgrey",high="blue")+
               labs(alpha = paste0("Relative abundance \n(standardized by max reads)"), col="Nearest Neighbour Ratio") +
               theme_classic() +
               theme(axis.text.y = element_blank(), axis.ticks.y=element_blank()) + xlab("Salinity") +
               geom_vline(data=data.frame(TurnoverPoints=localMax_bray), mapping=aes(xintercept=TurnoverPoints), lty=1, lwd=2, col="black", alpha=0.5) +
               facet_grid(typeSimple~. ,scales="free", space = "free_y")
           
    )
    # 
    # combinedTR_mb %>%
    #     ggplot() + geom_segment(aes(x=startSal, xend=endSal, y=OTU, yend=OTU, col=nnR), lwd= lwd_temp #100/nOTUs*1.2)
    #     ) +
    #     geom_point(data=outliers_points,aes(x=startSal, y=OTU,col=nnR), cex=100/nOTUs*1.2 )+
    #     scale_color_gradient2(low="red",mid="grey",high="blue")+
    #     labs(alpha = paste0("Relative abundance \n(standardized by max reads)"), col="Nearest Neighbour Ratio") +
    #     theme_classic() +
    #     theme(axis.text.y = element_blank(), axis.ticks.y=element_blank()) + xlab("Salinity") +
    #     geom_vline(data=data.frame(TurnoverPoints=localMax_bray), mapping=aes(xintercept=TurnoverPoints), lty=1, lwd=2, col="black", alpha=0.5)+
    #     facet_grid(typeSimple~. ,scales="free_y")
    
    ggsave(paste0(newDirName,"/",d,"_nnR_tolerancerange2.pdf"), height=10, width=6.5
           , combinedTR_mb  %>%
               ggplot() + geom_segment(aes(x=startSal, xend=endSal, y=OTU, yend=OTU),show.legend = FALSE, col=combinedTR_mb$rgb_col, lwd= 100/(nOTUs*1.2)
               ) +
               geom_point(data=outliers_points,aes(x=startSal, y=OTU),show.legend = FALSE, cex=1/(nOTUs*1.2), col=outliers_points$rgb_col)+
               labs(alpha = paste0("Relative abundance \n(standardized by max reads)"), col="Nearest Neighbour Ratio") +
               theme_classic() +
               theme(axis.text.y = element_blank(), axis.ticks.y=element_blank()) + xlab("Salinity") +
               geom_vline(data=data.frame(TurnoverPoints=localMax_bray), mapping=aes(xintercept=TurnoverPoints), lty=1, lwd=2, col="black", alpha=0.5) +
               facet_grid(typeSimple~. ,scales="free", space = "free_y")
    )
    
    ggsave(paste0(newDirName,"/",d,"_nnR_tolerancerange3.pdf"), height=10, width=6.5
           ,combinedTR_mb %>%
               ggplot() + geom_segment(aes(x=startSal, xend=endSal, y=OTU, yend=OTU, col=nnR, alpha=nnR_brack), lwd= 100/(nOTUs*1.2)
               ) +
               geom_point(data=outliers_points,aes(x=startSal, y=OTU,col=nnR, alpha=nnR_brack), cex=1/(nOTUs*1.2) )+
               scale_color_gradient2(low="red",mid="grey",high="blue")+
               labs(alpha = paste0("Relative Phylogenetic \ndistance \n(Brack vs Fresh/Marine)"), col="Nearest Neighbour Ratio") +
               theme_classic() +
               theme(axis.text.y = element_blank(), axis.ticks.y=element_blank()) + xlab("Salinity") +
               geom_vline(data=data.frame(TurnoverPoints=localMax_bray), mapping=aes(xintercept=TurnoverPoints), lty=1, lwd=2, col="black", alpha=0.5) +
               facet_grid(typeSimple~. ,scales="free", space = "free_y")
           
    )
    
}





