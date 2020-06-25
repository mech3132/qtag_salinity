#!bin/bash/Rscript

library("phyloseq")# For unifrac
library(MASS) # for beta div plotting
library("picante") # For unifrac
library("ape") # for tree?
library("tidyverse")
library("vegan") # for vegdist
library("mgcv") # for fitting splines

#### Load ####
# setwd("../SALBIN_20april2020/")
t.B16PWD <- "./tree_B16_filt.tre"
t.F16PWD <- "./tree_F16_filt.tre"
t.F18PWD <- "./tree_F18_filt.tre"

otu.B16PWD <- "./16sBaltic/OTUTableText.txt"
otu.F16PWD <- "./16sFraser/OTUTableText.txt"
otu.F18PWD <- "./18sFraser/OTUTableText.txt"

mf.B16PWD <- "../raw_otu_and_mf/Baltic_16S/MANUAL_INPUT_FILES/metadata_table.tsv"
mf.F16PWD <- "../raw_otu_and_mf/Fraser_16S/OUTPUT_FILES/MF_16sFraser_noConCOL.txt"
mf.F18PWD <- "../raw_otu_and_mf/Fraser_18S/OUTPUT_FILES/MF_18sFraser_noConCOL.txt"

mb.B16PWD <- "./16sBaltic/modelBoundaries_type.txt"
mb.F16PWD <- "./16sFraser/modelBoundaries_type.txt"
mb.F18PWD <- "./18sFraser/modelBoundaries_type.txt"

t.B16 <- read.tree(t.B16PWD)
t.F16 <- read.tree(t.F16PWD)
t.F18 <- read.tree(t.F18PWD)

otu.B16 <- read.delim(otu.B16PWD)
otu.F16 <- read.delim(otu.F16PWD)
otu.F18 <- read.delim(otu.F18PWD)

mf.B16 <- read.delim(mf.B16PWD)
mf.F16 <- read.delim(mf.F16PWD)
mf.F16$X.SampleID  <- gsub("-",".",mf.F16$X.SampleID, fixed = TRUE)
mf.F18 <- read.delim(mf.F18PWD)
mf.F18$X.SampleID  <- gsub("-",".",mf.F18$X.SampleID, fixed = TRUE)

mb.B16 <- read.delim(mb.B16PWD)
mb.F16 <- read.delim(mb.F16PWD)
mb.F18 <- read.delim(mb.F18PWD)

## Filter
tax.B16 <- otu.B16 %>% select(X.OTUID,taxonomy) %>%
    tidyr::separate(taxonomy, sep="; ", into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) %>%
    data.frame(row.names=1) 
tax.F16 <- otu.F16 %>% select(X.OTUID,taxonomy) %>%
    tidyr::separate(taxonomy, sep="; ", into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) %>%
    data.frame(row.names=1)
tax.F18 <- otu.F18 %>% select(X.OTUID,taxonomy) %>%
    tidyr::separate(taxonomy, sep="; __", into=c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"), fill="left") %>%
    data.frame(row.names=1)

otu.B16.adj <- otu.B16 %>% select(-c(X.OTUID,taxonomy))
rownames(otu.B16.adj) <- otu.B16$X.OTUID
otu.F16.adj <- otu.F16 %>% select(-c(X.OTUID,taxonomy))
rownames(otu.F16.adj) <- otu.F16$X.OTUID
otu.F18.adj <- otu.F18 %>% select(-c(X.OTUID,taxonomy))
rownames(otu.F18.adj) <- otu.F18$X.OTUID    

# Filter mf and otu
mf.B16 <- mf.B16[mf.B16$X.SampleID %in% colnames(otu.B16.adj),]
otu.B16.adj <- otu.B16.adj[,colnames(otu.B16.adj) %in% mf.B16$X.SampleID]

mf.F16 <- mf.F16[mf.F16$X.SampleID %in% colnames(otu.F16.adj),]
otu.F16.adj <- otu.F16.adj[,colnames(otu.F16.adj) %in% mf.F16$X.SampleID]

mf.F18 <- mf.F18[mf.F18$X.SampleID %in% colnames(otu.F18.adj),]
otu.F18.adj <- otu.F18.adj[,colnames(otu.F18.adj) %in% mf.F18$X.SampleID]


# make into relative abundance
otu.B16.adj.relAbund <- apply(otu.B16.adj, 2, FUN=function(c){c/sum(c)})
otu.F16.adj.relAbund <- apply(otu.F16.adj, 2, FUN=function(c){c/sum(c)})
otu.F18.adj.relAbund <- apply(otu.F18.adj, 2, FUN=function(c){c/sum(c)})

## Function to count percent lost/gained
coreKept <- function(otu.temp, starting, ending) {
    allChange = c()
    otuNames <- rownames(otu.temp)
    for ( s in starting) {
        for ( e in ending) {
            tempDat <- data.frame(start=otu.temp[,s], end=otu.temp[,e], row.names = otuNames) %>%
                filter(start+end!=0)%>%
                mutate(change=end-start
                       , lost = ifelse(change>0,0,change)
                       , gain = ifelse(change<0,0,change)
                       , coreKept = ifelse(change<0, end, start)) %>%
                colSums() 
            allChange = c(allChange, as.numeric(tempDat['coreKept']))
            
            
        }
    }
    return(allChange)
}


OTUturnover <- function(otu.temp, starting, ending) {
    otu.PA <- otu.temp
    otu.PA[otu.PA>0] <- 1
    allNew = c()
    allLost = c()
    totCount = c()
    otuNames <- rownames(otu.PA)
    for ( s in starting) {
        for ( e in ending) {
            tempDat_gained <- data.frame(start=otu.PA[,s], end=otu.PA[,e], row.names = otuNames) %>%
                filter(start==0 & end == 1) %>%
                select(end) %>% colSums() %>% as.vector()
            tempDat_lost <- data.frame(start=otu.PA[,s], end=otu.PA[,e], row.names = otuNames) %>%
                filter(start==1 & end == 0) %>%
                select(start) %>% colSums() %>% as.vector()
            totDat <-  sum(otu.PA[,e])
            allNew <- c(allNew, tempDat_gained)
            allLost <- c(allLost, tempDat_lost)
            totCount <- c(totCount, totDat)
            
        }
    }
    return(list(allNew,allLost, totCount))
}

# Create beta div
for ( ty in c("B16","F16","F18")) {
    otu.temp <- get(paste0("otu.",ty,".adj.relAbund"))
    mf.temp <- get(paste0("mf.",ty))
    tree.temp <- get(paste0("t.",ty))
    
    temp.uwunifrac <- as.matrix(unifrac(t(otu.temp), tree.temp))
    temp.bray <- as.matrix(vegdist(t(otu.temp), method = "bray"))
    # temp.jaccard <- as.matrix(vegdist(t(otu.temp), method = "jaccard"))
    
    assign(paste0("dm_uwunifrac.",ty),temp.uwunifrac)
    assign(paste0("dm_bray.",ty),temp.bray)
    # assign(paste0("dm_jaccard.",ty),temp.jaccard)
}

#### Start ####
dir.create("CommunityTurnover")
for ( ty in c("B16","F16","F18")) {
    otu.temp <- get(paste0("otu.",ty,".adj.relAbund"))
    mf.temp <- get(paste0("mf.",ty))
    tree.temp <- get(paste0("t.",ty))
    
    temp.uwunifrac <- get(paste0("dm_uwunifrac.",ty))
    temp.bray <- get(paste0("dm_bray.",ty))
    # temp.jaccard <- get(paste0("dm_jaccard.",ty))
    
    mf.temp <- mf.temp %>%
        filter(!is.na(as.numeric(SalinityEnviron))) %>%
        mutate(SalRound = as.numeric(SalinityEnviron))
        # mutate(SalRound = round(as.numeric(SalinityEnviron)))
    # SalRound <- sort(unique(mf.temp$SalinityEnviron))
    SalRound <- sort(unique(mf.temp$SalRound))
    
    jumps <- data.frame()
    total <- (length(SalRound)-1)
    # create progress bar
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    for ( s in 1:(length(SalRound)-1) ) {
        sstart <- SalRound[s]
        send <- SalRound[s+1]
        starting <- mf.temp[mf.temp$SalRound==sstart,"X.SampleID"]
        ending <- mf.temp[mf.temp$SalRound== send,"X.SampleID"]

        distUWUnifrac <- as.vector(temp.uwunifrac[starting,ending])
        distBray <- as.vector(temp.bray[starting,ending])
        # distJaccard <- as.vector(temp.jaccard[starting,ending])
        core <- coreKept(otu.temp, starting, ending)
        turnover <- OTUturnover(otu.temp, starting, ending)
        
        jumps <- rbind(jumps, data.frame(sstart,send
                                             ,distUWUnifrac
                                         ,distBray
                                         # ,distJaccard
                                         ,core, gained=turnover[[1]], lost=turnover[[2]], totalRich=turnover[[3]]))
        
        Sys.sleep(0.1)
        setTxtProgressBar(pb, s)
    }
    close(pb)
    assign(paste0(ty,".jumps"),jumps)
    
}
save(B16.jumps, file = "CommunityTurnover/B16.jumps.RData")
save(F16.jumps, file = "CommunityTurnover/F16.jumps.RData")
save(F18.jumps, file = "CommunityTurnover/F18.jumps.RData")

B16.jumps %>%
    as_tibble() %>%
    gather(-c(sstart,send), key=metrics, value=dist) %>% 
    ggplot() + geom_point(aes(x=(send), y=dist, col=metrics)) +
    geom_smooth(aes(x=(send), y=dist, col=metrics)) +
    facet_grid(metrics~., scales="free_y")
F16.jumps %>%
    as_tibble() %>%
    gather(-c(sstart,send), key=metrics, value=dist) %>% 
    ggplot() + geom_point(aes(x=(send), y=dist, col=metrics)) +
    geom_smooth(aes(x=(send), y=dist, col=metrics)) +
    facet_grid(metrics~., scales="free_y")
F18.jumps %>%
    as_tibble() %>%
    gather(-c(sstart,send), key=metrics, value=dist) %>% 
    ggplot() + geom_point(aes(x=(send), y=dist, col=metrics)) +
    geom_smooth(aes(x=(send), y=dist, col=metrics)) +
    facet_grid(metrics~., scales="free_y")

### Replace as need; check out all plots

gam_bray.B16 <- gam(distBray ~ s(send), data=B16.jumps, family = gaussian)
# summary(gam_bray.B16)
# gam.check(gam_bray.B16)
plot(gam_bray.B16$fitted.values ~ B16.jumps$send, type="l")

gam_bray.F16 <- gam(distBray ~ s(send), data=F16.jumps, family=gaussian)
# summary(gam_bray.F16)
# gam.check(gam_bray.F16)
plot(gam_bray.F16$fitted.values ~ F16.jumps$send, type="l")

gam_bray.F18 <- gam(distBray ~ s(send), data=F18.jumps, family=gaussian)
# summary(gam_bray.F18)
# gam.check(gam_bray.F18)
plot(gam_bray.F18$fitted.values ~ F18.jumps$send, type="l")



gam_core.B16 <- gam(core ~ s(send), data=B16.jumps, family = gaussian)
# summary(gam_core.B16)
# gam.check(gam_core.B16)
plot(gam_core.B16$fitted.values ~ B16.jumps$send, type="l", col="darkgreen")
lines(gam_bray.B16$fitted.values ~ B16.jumps$send,col="blue")

gam_core.F16 <- gam(core ~ s(send), data=F16.jumps, family = gaussian)
# summary(gam_core.F16)
# gam.check(gam_core.F16)
plot(gam_core.F16$fitted.values ~ F16.jumps$send, type="l")

gam_core.F18 <- gam(core ~ s(send), data=F18.jumps, family = gaussian)
# summary(gam_core.F18)
# gam.check(gam_core.F18)
plot(gam_core.F18$fitted.values ~ F18.jumps$send, type="l")




gam_unifrac.B16 <- gam(distUWUnifrac ~ s(send), data=B16.jumps, family = gaussian)
# summary(gam_unifrac.B16)
# gam.check(gam_unifrac.B16)
plot(gam_unifrac.B16$fitted.values ~ B16.jumps$send, type="l", col="darkgreen")
lines(gam_bray.B16$fitted.values ~ B16.jumps$send,col="blue")

gam_unifrac.F16 <- gam(distUWUnifrac ~ s(send), data=F16.jumps, family = gaussian)
# summary(gam_unifrac.F16)
# gam.check(gam_unifrac.F16)
plot(gam_unifrac.F16$fitted.values ~ F16.jumps$send, type="l")
lines(gam_bray.F16$fitted.values ~ F16.jumps$send,col="blue")

gam_unifrac.F18 <- gam(distUWUnifrac ~ s(send), data=F18.jumps, family = gaussian)
# summary(gam_unifrac.F18)
# gam.check(gam_unifrac.F18)
plot(gam_unifrac.F18$fitted.values ~ F18.jumps$send, type="l")


# All of them together
plot(NULL, xlim = c(0,35), ylim=c(0,0.8))
lines(gam_bray.B16$fitted.values ~ B16.jumps$send, col="darkgreen")
lines(gam_bray.F16$fitted.values ~ F16.jumps$send, col="orange")
# lines(gam_core.F18$fitted.values ~ F18.jumps$send, col="purple")

# All of them together, core
plot(NULL, xlim = c(0,35), ylim=c(0,0.8))
lines(gam_core.B16$fitted.values ~ B16.jumps$send, col="darkgreen")
lines(gam_core.F16$fitted.values ~ F16.jumps$send, col="orange")

## Get local minima for "core", or local maxima for "distance"
localMax <- function(x) {
  return(which(diff(sign(diff(x)))==-2)+1)
}
localMin <- function(x) {
  return(which(diff(sign(diff(x)))==+2)+1)
}

B16.jumps$send[localMax(gam_core.B16$fitted.values)]

# 
# ### If we filter communities to JUST brackish species, maybe we can see signal of turnover
# # Without fresh and marine making weird lumps at either end?
# brackish.B16 <- mb.B16 %>% filter(typeSimple=="brackishRestricted") %>% select(taxa) %>% pull()
# otu.B16.adj.relAbund.brackOnly <- otu.B16.adj.relAbund[brackish.B16,]
# 
# brackish.F16 <- mb.F16 %>% filter(typeSimple=="brackishRestricted") %>% select(taxa) %>% pull()
# otu.F16.adj.relAbund.brackOnly <- otu.F16.adj.relAbund[brackish.F16,]
# 
# brackish.F18 <- as.character(mb.F18 %>% filter(typeSimple=="brackishRestricted") %>% select(taxa) %>% pull())
# otu.F18.adj.relAbund.brackOnly <- otu.F18.adj.relAbund[brackish.F18,]
# 
# # Create beta div for only brackish
# for ( ty in c("B16","F16","F18")) {
#     otu.temp <- get(paste0("otu.",ty,".adj.relAbund.brackOnly"))
#     mf.temp <- get(paste0("mf.",ty))
#     tree.temp <- get(paste0("t.",ty))
#     
#     temp.uwunifrac <- as.matrix(unifrac(t(otu.temp), tree.temp))
#     temp.bray <- as.matrix(vegdist(t(otu.temp), method = "bray"))
#     temp.jaccard <- as.matrix(vegdist(t(otu.temp), method = "jaccard"))
#     
#     assign(paste0("dm_uwunifrac_brackonly.",ty),temp.uwunifrac)
#     assign(paste0("dm_bray_brackonly.",ty),temp.bray)
#     assign(paste0("dm_jaccard_brackonly.",ty),temp.jaccard)
# }
# 
# #### Filter ONLY brackish ####
# for ( ty in c("B16","F16","F18")) {
#     otu.temp <- get(paste0("otu.",ty,".adj.relAbund.brackOnly"))
#     mf.temp <- get(paste0("mf.",ty))
#     tree.temp <- get(paste0("t.",ty))
#     
#     temp.uwunifrac <- get(paste0("dm_uwunifrac_brackonly.",ty))
#     temp.bray <- get(paste0("dm_bray_brackonly.",ty))
#     temp.jaccard <- get(paste0("dm_jaccard_brackonly.",ty))
#     
#     mf.temp <- mf.temp %>%
#         filter(!is.na(as.numeric(SalinityEnviron))) %>%
#         mutate(SalRound = as.numeric(SalinityEnviron))
#     # mutate(SalRound = round(as.numeric(SalinityEnviron)))
#     # SalRound <- sort(unique(mf.temp$SalinityEnviron))
#     SalRound <- sort(unique(mf.temp$SalRound))
#     
#     jumps <- data.frame()
#     total <- (length(SalRound)-1)
#     # create progress bar
#     pb <- txtProgressBar(min = 0, max = total, style = 3)
#     for ( s in 1:(length(SalRound)-1) ) {
#         sstart <- SalRound[s]
#         send <- SalRound[s+1]
#         starting <- mf.temp[mf.temp$SalRound==sstart,"X.SampleID"]
#         ending <- mf.temp[mf.temp$SalRound== send,"X.SampleID"]
#         
#         distUWUnifrac <- as.vector(temp.uwunifrac[starting,ending])
#         distBray <- as.vector(temp.bray[starting,ending])
#         distJaccard <- as.vector(temp.jaccard[starting,ending])
#         core <- coreKept(otu.temp, starting, ending)
#         turnover <- OTUturnover(otu.temp, starting, ending)
#         
#         jumps <- rbind(jumps, data.frame(sstart,send
#                                          ,distUWUnifrac,distBray,distJaccard,core, gained=turnover[[1]], lost=turnover[[2]], totalRich=turnover[[3]]))
#         
#         Sys.sleep(0.1)
#         setTxtProgressBar(pb, s)
#     }
#     close(pb)
#     assign(paste0(ty,".jumpsBrack"),jumps)
#     
# }
# 
# 
# B16.jumpsBrack %>%
#     as_tibble() %>%
#     gather(-c(sstart,send), key=metrics, value=dist) %>% 
#     ggplot() + geom_point(aes(x=(send), y=dist, col=metrics)) +
#     geom_smooth(aes(x=(send), y=dist, col=metrics)) +
#     facet_grid(metrics~., scales="free_y")
# F16.jumpsBrack %>%
#     as_tibble() %>%
#     gather(-c(sstart,send), key=metrics, value=dist) %>% 
#     ggplot() + geom_point(aes(x=(send), y=dist, col=metrics)) +
#     geom_smooth(aes(x=(send), y=dist, col=metrics)) +
#     facet_grid(metrics~., scales="free_y")
# F18.jumpsBrack %>%
#     as_tibble() %>%
#     gather(-c(sstart,send), key=metrics, value=dist) %>% 
#     ggplot() + geom_point(aes(x=(send), y=dist, col=metrics)) +
#     geom_smooth(aes(x=(send), y=dist, col=metrics)) +
#     facet_grid(metrics~., scales="free_y")
# 
# ### Replace as need; check out all plots
# test.gam <- gam(core ~ s(send), data=F16.jumpsBrack)
# plot(test.gam, residuals=TRUE, cex=2)
# 
# summary(test.gam)
# gam.check(test.gam)
# 
# 
# ### PCOA plots
# dm.temp <- dm_bray.B16
# mf.temp <- mf.B16
# otu.temp <- otu.B16.adj.relAbund.brackOnly
# 
# temp.nmds <- isoMDS(dist(dm.temp))
# cbind(mf.temp[match(rownames(temp.nmds$points),mf.temp$X.SampleID),], temp.nmds$points) %>%
#     ggplot() + geom_point(aes(x=`1`, y=`2`, col=as.numeric(SalinityEnviron)))
# 
# 
# kmeans.2 <- kmeans(t(otu.temp), centers = 2)
# kmeans.3 <- kmeans(t(otu.temp), centers = 3)
# kmeans.4 <- kmeans(t(otu.temp), centers = 4)
# kmeans.5 <- kmeans(t(otu.temp), centers = 5)
# kmeans.6 <- kmeans(t(otu.temp), centers = 6)
# 
# #### SS method
# # function to compute total within-cluster sum of square 
# wss <- function(k) {
#     kmeans(t(otu.temp), k, nstart = 10 )$tot.withinss
# }
# 
# # Compute and plot wss for k = 1 to k = 15
# k.values <- 1:15
# 
# # extract wss for 2-15 clusters
# wss_values <- map_dbl(k.values, wss)
# 
# plot(k.values, wss_values,
#      type="b", pch = 19, frame = FALSE, 
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares")
# 
# 
# cbind(mf.temp[match(rownames(temp.nmds$points),mf.temp$X.SampleID),], temp.nmds$points
#       , kmeans2=kmeans.2$cluster[match(rownames(temp.nmds$points),names(kmeans.2$cluster))]
#       , kmeans3=kmeans.3$cluster[match(rownames(temp.nmds$points),names(kmeans.3$cluster))]
#       , kmeans4=kmeans.4$cluster[match(rownames(temp.nmds$points),names(kmeans.4$cluster))]
#       , kmeans5=kmeans.5$cluster[match(rownames(temp.nmds$points),names(kmeans.5$cluster))]
#       , kmeans6=kmeans.6$cluster[match(rownames(temp.nmds$points),names(kmeans.6$cluster))]) %>%
#     ggplot() + geom_point(aes(x=`1`, y=`2`, col=as.numeric(SalinityEnviron), pch=factor(kmeans2)))
# 
# 
# otu.B16.adj.relAbund[,match(mf.B16[order(mf.B16$SalinityEnviron),"X.SampleID"], colnames(otu.B16.adj.relAbund))]
# heatmap(otu.B16.adj.relAbund[,match(mf.B16[order(mf.B16$SalinityEnviron),"X.SampleID"], colnames(otu.B16.adj.relAbund))]
#         , Colv = NA, scale="row")
# 
