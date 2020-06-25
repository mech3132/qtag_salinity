#!bin/bash

###### Community cross-sections ########
library(tidyverse)

#### Load ####
dirNames <- unlist(read.delim("OUTPUT/dirNames.txt", header=FALSE))
output <- "./OUTPUT/03g_BrackishOTUs"

otu.B16PWD <- paste0(dirNames[1],"/OTUTableText.txt")
otu.F16PWD <- paste0(dirNames[2],"/OTUTableText.txt")
otu.F18PWD <- paste0(dirNames[3],"/OTUTableText.txt")

mf.B16PWD <- "./INPUT_DATA/Baltic_16S/metadata_table.tsv"
mf.F16PWD <- "./INPUT_DATA/Fraser_16S/MF_16sFraser_noConCOL.txt"
mf.F18PWD <- "./INPUT_DATA/Fraser_18S/MF_18sFraser_noConCOL.txt"

mb.B16PWD <- paste0(dirNames[1],"/modelBoundaries_type.txt")
mb.F16PWD <- paste0(dirNames[2],"/modelBoundaries_type.txt")
mb.F18PWD <- paste0(dirNames[3],"/modelBoundaries_type.txt")

tolRange.B16PWD <- paste0(dirNames[1],"/toleranceRanges.txt")
tolRange.F16PWD <- paste0(dirNames[2],"/toleranceRanges.txt")
tolRange.F18PWD <- paste0(dirNames[3],"/toleranceRanges.txt")

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

tolR.B16 <- read.delim(tolRange.B16PWD)
tolR.F16 <- read.delim(tolRange.F16PWD)
tolR.F18 <- read.delim(tolRange.F18PWD)

## Filter
tax.B16 <- otu.B16 %>% select(X.OTUID,taxonomy) %>%
    tidyr::separate(taxonomy, sep="; ", into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) %>%
    data.frame(row.names=1) 
tax.F16 <- otu.F16 %>% select(X.OTUID,taxonomy) %>%
    tidyr::separate(taxonomy, sep="; ", into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) %>%
    data.frame(row.names=1)
tax.F18 <- otu.F18 %>% select(X.OTUID,taxonomy) %>%
    tidyr::separate(taxonomy, sep="; __", into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species","LAST")) %>%
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
# 

### quick test to see what brackish they shared

sharedOTUs <- mb.B16[mb.B16$taxa %in% mb.F16$taxa,"taxa"]
sharedTax <- tax.B16[match(sharedOTUs, rownames(tax.B16)),] %>% 
    rownames_to_column("OTU")%>%
    mutate(Class=gsub("D.*__","",Class)
           ,Family=gsub("D.*__","",Family) ) %>% 
    unite(Class, Family, col=TaxLabel, sep=": ") %>%
    select(OTU,TaxLabel) 

# what about NOT shared
notsharedOTUs.B16 <- mb.B16[!(mb.B16$taxa %in% mb.F16$taxa),"taxa"]
notsharedOTUs.F16 <- mb.F16[!(mb.F16$taxa %in% mb.B16$taxa),"taxa"]
NOTsharedTax.B16 <- tax.B16[match(notsharedOTUs.B16, rownames(tax.B16)),] %>% 
    rownames_to_column("OTU")%>%
    mutate(Class=gsub("D.*__","",Class)
           ,Family=gsub("D.*__","",Family) ) %>% 
    unite(Class, Family, col=TaxLabel, sep=": ") %>%
    select(OTU,TaxLabel) 
NOTsharedTax.F16 <- tax.F16[match(notsharedOTUs.F16, rownames(tax.F16)),] %>% 
    rownames_to_column("OTU")%>%
    mutate(Class=gsub("D.*__","",Class)
           ,Family=gsub("D.*__","",Family) ) %>% 
    unite(Class, Family, col=TaxLabel, sep=": ") %>%
    select(OTU,TaxLabel) 

# Change typeSimple to better names
mb.B16 <- mb.B16 %>% mutate(typeSimple=ifelse(typeSimple=="freshRestricted","Fresh",ifelse(typeSimple=="marineRestricted","Marine",ifelse(typeSimple=="brackishRestricted","Brackish","Unclassified"))))
mb.F16 <- mb.F16 %>% mutate(typeSimple=ifelse(typeSimple=="freshRestricted","Fresh",ifelse(typeSimple=="marineRestricted","Marine",ifelse(typeSimple=="brackishRestricted","Brackish","Unclassified"))))
mb.F18 <- mb.F18 %>% mutate(typeSimple=ifelse(typeSimple=="freshRestricted","Fresh",ifelse(typeSimple=="marineRestricted","Marine",ifelse(typeSimple=="brackishRestricted","Brackish","Unclassified"))))

dir.create(paste0(output))
dir.create(paste0(output,"/IndividualPlots"))
total <- length(sharedOTUs)
pb <- txtProgressBar(min = 0, max = total, style = 3)
i <- 0
for ( otu in sharedOTUs) {
    B.type <- mb.B16[mb.B16$taxa==otu,"typeSimple"]
    F.type <- mb.F16[mb.F16$taxa==otu,"typeSimple"]
    
    B16.temp <- otu.B16.adj.relAbund[rownames(otu.B16.adj.relAbund)==otu,] %>% as.data.frame() %>%
        rownames_to_column(var="SampleID") %>% rename(relAbund=".") %>%
        mutate(Salinity=mf.B16[match(SampleID, mf.B16$X.SampleID),"SalinityEnviron"]
               , Dataset=paste0("Baltic Sea (",B.type,")"))
    F16.temp <- otu.F16.adj.relAbund[rownames(otu.F16.adj.relAbund)==otu,] %>% as.data.frame() %>%
        rownames_to_column(var="SampleID") %>% rename(relAbund=".") %>%
        mutate(Salinity=mf.F16[match(SampleID, mf.F16$X.SampleID),"SalinityEnviron"]
               , Dataset=paste0("Fraser River (",F.type,")"))
    
    ratioBF <- max(B16.temp$relAbund)/max(F16.temp$relAbund)*1
    ggsave(filename = paste0(output,"/IndividualPlots/",otu,"_B-",B.type,"_F-",F.type,".pdf"), width=5, height=3
        ,rbind(B16.temp, F16.temp) %>%
        mutate(Salinity=as.numeric(Salinity)
               , adjRatio = ifelse(Dataset==paste0("Baltic Sea (",B.type,")"), 1, ratioBF)
               , adjRelAbund = relAbund*adjRatio) %>%
        ggplot(aes(x=Salinity, y=adjRelAbund, col=Dataset)) +
        geom_point(pch=19) + scale_color_manual(values=c("purple","darkgreen")) +ylab("Relative Abundance (Baltic)") +
            scale_y_continuous(sec.axis=sec_axis(trans=~./ratioBF, name = "Relative Abundance (Fraser)")) +
            theme(axis.title.y.right = element_text(color = "darkgreen")
                  ,axis.title.y.left = element_text(color = "purple")
                  )+
            ggtitle(paste0(sharedTax[sharedTax$OTU==otu,"TaxLabel"]))
    )
    i <- i+1
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, i)
}

#### What proportion of the community is brackish and shared ####
sharedOTUs_types <- data.frame(OTU=sharedOTUs,B.type=mb.B16[match(sharedOTUs,mb.B16$taxa),"typeSimple"]
           , F.type=mb.F16[match(sharedOTUs,mb.F16$taxa),"typeSimple"])
brack_sharedOTUs <- sharedOTUs_types %>% filter(B.type=="Brackish",F.type=="Brackish")

brack_shared_tax <- tax.B16[match(brack_sharedOTUs$OTU, rownames(tax.B16)),] %>% 
    rownames_to_column("OTU")%>%
    mutate(Class=gsub("D.*__","",Class)
           ,Family=gsub("D.*__","",Family) ) %>% 
    unite(Class, Family, col=TaxLabel, sep=": ") %>%
    select(OTU,TaxLabel) 

### BALTIC
otu_bytype.B16 <- otu.B16.adj %>%
    mutate(Type=mb.B16[match(rownames(otu.B16.adj), mb.B16$taxa),"typeSimple"]
           , taxa = brack_shared_tax[match(rownames(otu.B16.adj), brack_shared_tax$OTU),"TaxLabel"]
           , Group=ifelse(is.na(taxa),Type,"SHARED")) %>%
    gather(-c(Type,taxa,Group), key=SampleID, value=Count) %>%
    group_by(Group, SampleID) %>%
    summarize(reads=sum(Count)) %>%
    mutate(Salinity=mf.B16[match(SampleID, mf.B16$X.SampleID),"SalinityEnviron"]) %>%
    mutate(Salinity_round=round(Salinity)) %>%
    ungroup() %>% group_by(Salinity_round) %>%
    mutate(relAbund=reads/sum(reads)) %>% ungroup() %>%
    arrange(Salinity) %>% mutate(SampleID = factor(SampleID, levels=unique(SampleID))) %>%
    mutate(Group=ifelse(Group=="freshRestricted","Fresh"
                       ,ifelse(Group=="marineRestricted","Marine"
                               ,ifelse(Group=="brackishRestricted","Brackish"
                                       , ifelse(Group=="noclass","Unclassified",Group))))) %>%
    mutate(Group=factor(Group, levels=c("Marine","Brackish","SHARED","Fresh","Unclassified")))

ggsave(filename = paste0(output,"/Baltic_proportionShared.pdf")
       ,otu_bytype.B16 %>%
           ggplot() + geom_col(aes(x=Salinity_round, y=relAbund, fill=Group)) + ylab("Relative Abundance")+xlab("Salinity")+
           scale_fill_manual(values=c("red","purple","magenta","blue","black")))





## FRASER
otu_bytype.F16 <- otu.F16.adj %>%
    mutate(Type=mb.F16[match(rownames(otu.F16.adj), mb.F16$taxa),"typeSimple"]
           , taxa = brack_shared_tax[match(rownames(otu.F16.adj), brack_shared_tax$OTU),"TaxLabel"]
           , Group=ifelse(is.na(taxa),Type,"SHARED")) %>%
    gather(-c(Type,taxa,Group), key=SampleID, value=Count) %>%
    group_by(Group, SampleID) %>%
    summarize(reads=sum(Count)) %>%
    mutate(Salinity=mf.F16[match(SampleID, mf.F16$X.SampleID),"SalinityEnviron"]) %>%
    mutate(Salinity_round=round(as.numeric(Salinity))) %>%
    ungroup() %>% group_by(Salinity_round) %>%
    mutate(relAbund=reads/sum(reads)) %>% ungroup() %>%
    arrange(Salinity) %>% mutate(SampleID = factor(SampleID, levels=unique(SampleID))) %>%
    mutate(Group=ifelse(Group=="freshRestricted","Fresh"
                        ,ifelse(Group=="marineRestricted","Marine"
                                ,ifelse(Group=="brackishRestricted","Brackish"
                                        , ifelse(Group=="noclass","Unclassified",Group))))) %>%
    mutate(Group=factor(Group, levels=c("Marine","Brackish","SHARED","Fresh","Unclassified")))

ggsave(filename = paste0(output,"/Fraser_proportionShared.pdf")
       ,otu_bytype.F16 %>%
           ggplot() + geom_col(aes(x=Salinity_round, y=relAbund, fill=Group)) +ylab("Relative Abundance")+xlab("Salinity")+
           scale_fill_manual(values=c("red","purple","magenta","blue","black"))
)


#### Taxa summaries, taxonomy of shared ####

### Baltic
ncolors <- nrow(brack_sharedOTUs)
filtColors <- colors()[-grep("white|gr(a|e)y", colors())]
set.seed(2345)
taxacolors <- sample(filtColors, size=ncolors, replace=FALSE)
names(taxacolors) <- brack_shared_tax$TaxLabel

otu_bytype.B16_sharedOnly <- otu.B16.adj %>%
    mutate(Type=mb.B16[match(rownames(otu.B16.adj), mb.B16$taxa),"typeSimple"]
           , taxa = brack_shared_tax[match(rownames(otu.B16.adj), brack_shared_tax$OTU),"TaxLabel"]
           , Group=ifelse(is.na(taxa),Type,taxa)) %>%
    gather(-c(Type,taxa,Group), key=SampleID, value=Count) %>%
    group_by(Group, SampleID) %>%
    summarize(reads=sum(Count)) %>%
    mutate(Salinity=mf.B16[match(SampleID, mf.B16$X.SampleID),"SalinityEnviron"]) %>%
    mutate(Salinity_round=round(as.numeric(Salinity))) %>%
    ungroup() %>% group_by(Salinity_round) %>%
    mutate(relAbund=reads/sum(reads)) %>% ungroup() %>%
    filter(!(Group %in% c("Fresh","Marine","Brackish","Unclassified"))) %>%
    arrange(Salinity) %>% mutate(SampleID = factor(SampleID, levels=unique(SampleID))) 

ggsave(filename = paste0(output,"/Baltic_proportionShared_brackonlytaxonomy.pdf"), width=10, height=6
       , otu_bytype.B16_sharedOnly %>%
           ggplot() + geom_col(aes(x=Salinity_round, y=relAbund, fill=Group)) +ylab("Relative Abundance")+xlab("Salinity")+
           scale_fill_manual(values=taxacolors)
)


### Fraser
otu_bytype.F16_sharedOnly <- otu.F16.adj %>%
    mutate(Type=mb.F16[match(rownames(otu.F16.adj), mb.F16$taxa),"typeSimple"]
           , taxa = brack_shared_tax[match(rownames(otu.F16.adj), brack_shared_tax$OTU),"TaxLabel"]
           , Group=ifelse(is.na(taxa),Type,taxa)) %>%
    gather(-c(Type,taxa,Group), key=SampleID, value=Count) %>%
    group_by(Group, SampleID) %>%
    summarize(reads=sum(Count)) %>%
    mutate(Salinity=mf.F16[match(SampleID, mf.F16$X.SampleID),"SalinityEnviron"]) %>%
    mutate(Salinity_round=round(as.numeric(Salinity))) %>%
    ungroup() %>% group_by(Salinity_round) %>%
    mutate(relAbund=reads/sum(reads)) %>% ungroup() %>%
    filter(!(Group %in% c("Fresh","Marine","Brackish","Unclassified"))) %>%
    arrange(Salinity) %>% mutate(SampleID = factor(SampleID, levels=unique(SampleID))) 

ggsave(filename = paste0(output,"/Fraser_proportionShared_brackonlytaxonomy.pdf"),width=10, height=6
, otu_bytype.F16_sharedOnly %>%
           ggplot() + geom_col(aes(x=Salinity_round, y=relAbund, fill=Group)) +ylab("Relative Abundance")+xlab("Salinity")+
           scale_fill_manual(values=taxacolors)
)


#### Tolerance range comparison ####
### Plotting quick and dirty tolerance range plot

tolR.B16_shared <- tolR.B16 %>%
    filter(OTU %in% sharedOTUs, position=="withinBounds") %>%
    group_by(OTU) %>%
    summarize(minSal=min(startSal), maxSal=max(endSal)) %>% ungroup() %>%
    mutate(Dataset="Baltic"
           , Type=mb.B16[match(OTU, mb.B16$taxa),"typeSimple"]) 

tolR.F16_shared <- tolR.F16 %>%
    filter(OTU %in% sharedOTUs, position=="withinBounds") %>%
    group_by(OTU) %>%
    summarize(minSal=min(startSal), maxSal=max(endSal)) %>% ungroup() %>%
    mutate(Dataset="Fraser"
           , Type=mb.F16[match(OTU, mb.F16$taxa),"typeSimple"]) 

ggsave(filename = paste0(output,"/ToleranceRangeCompare.pdf"), width=7, height=10
       ,rbind(tolR.B16_shared, tolR.F16_shared) %>%
           group_by(OTU) %>%
           mutate(midSal=mean(c(minSal, maxSal))) %>% ungroup() %>%
           arrange(Dataset,-midSal) %>% mutate(OTU=factor(OTU, levels=unique(OTU))) %>%
           # filter(OTU%in% sharedOTUs[1:50]) %>%
           ggplot() + 
           geom_segment(aes(x=minSal, xend=maxSal, y=OTU, yend=OTU, col=Type))+
           scale_color_manual(values=c("purple","blue","red","black")) +
           facet_grid(.~Dataset, scales = "fixed")+xlab("Salinity")
       # theme(strip.text.y = element_text(angle = 0)
       # , axis.text.y = element_blank())
)

#### Looking at all brackish OTUs, even those not shared ####


### Baltic
balticBrack_OTUs <- mb.B16 %>% filter(typeSimple=="Brackish") %>% select(taxa) 
baltic_brack_tax <- data.frame(balticBrack_OTUs, tax.B16[match(balticBrack_OTUs$taxa, rownames(tax.B16)),])%>% 
    rownames_to_column("OTU")%>%
    mutate(Class=gsub("D.*__","",Class)
           ,Family=gsub("D.*__","",Family) ) %>% 
    unite(Class, Family, col=TaxLabel, sep=": ") 

otu_bytype.B16_all <- otu.B16.adj %>%
    mutate(Type=mb.B16[match(rownames(otu.B16.adj), mb.B16$taxa),"typeSimple"]
           , taxa = baltic_brack_tax[match(rownames(otu.B16.adj), baltic_brack_tax$OTU),"TaxLabel"]
           , Order = baltic_brack_tax[match(rownames(otu.B16.adj), baltic_brack_tax$OTU),"Order"]
           , Group=ifelse(is.na(taxa),Type,taxa)) %>%
    gather(-c(Type,taxa,Group,Order), key=SampleID, value=Count) %>%
    group_by(Group, Order, SampleID) %>%
    summarize(reads=sum(Count, na.rm = TRUE)) %>%
    mutate(Salinity=mf.B16[match(SampleID, mf.B16$X.SampleID),"SalinityEnviron"]) %>%
    mutate(Salinity_round=round(as.numeric(Salinity))) %>%
    ungroup() %>% group_by(Salinity_round) %>%
    mutate(relAbund=reads/sum(reads)) %>% ungroup() %>%
    filter(!(Group %in% c("Fresh","Marine","Brackish","Unclassified"))) %>%
    arrange(Salinity) %>% mutate(SampleID = factor(SampleID, levels=unique(SampleID))) %>%
    separate(Group, into=c("Class","Family"), sep=": ", remove=FALSE) %>%
    mutate(Order=gsub("D.*__","",Order)) %>%
    unite(Class, Order, col=Group_order, sep=": ", remove=FALSE)%>%
    mutate(Dataset="Baltic")
### Fraser
fraserBrack_OTUs <- mb.F16 %>% filter(typeSimple=="Brackish") %>% select(taxa) 
fraser_brack_tax <- data.frame(fraserBrack_OTUs, tax.F16[match(fraserBrack_OTUs$taxa, rownames(tax.F16)),])%>% 
    rownames_to_column("OTU")%>%
    mutate(Class=gsub("D.*__","",Class)
           ,Family=gsub("D.*__","",Family) ) %>% 
    unite(Class, Family, col=TaxLabel, sep=": ", remove=FALSE) 

otu_bytype.F16_all <- otu.F16.adj %>%
    mutate(Type=mb.F16[match(rownames(otu.F16.adj), mb.F16$taxa),"typeSimple"]
           , taxa = fraser_brack_tax[match(rownames(otu.F16.adj), fraser_brack_tax$OTU),"TaxLabel"]
           , Order = fraser_brack_tax[match(rownames(otu.F16.adj), fraser_brack_tax$OTU),"Order"]
           , Group=ifelse(is.na(taxa),Type,taxa)) %>%
    gather(-c(Type,taxa,Group,Order), key=SampleID, value=Count) %>%
    group_by(Group, Order, SampleID) %>%
    summarize(reads=sum(Count, na.rm = TRUE)) %>%
    mutate(Salinity=mf.F16[match(SampleID, mf.F16$X.SampleID),"SalinityEnviron"]) %>%
    mutate(Salinity_round=round(as.numeric(Salinity))) %>%
    ungroup() %>% group_by(Salinity_round) %>%
    mutate(relAbund=reads/sum(reads)) %>% ungroup() %>%
    filter(!(Group %in% c("Fresh","Marine","Brackish","Unclassified"))) %>%
    arrange(Salinity) %>% mutate(SampleID = factor(SampleID, levels=unique(SampleID))) %>%
    separate(Group, into=c("Class","Family"), sep=": ", remove=FALSE) %>%
    mutate(Order=gsub("D.*__","",Order)) %>%
    unite(Class, Order, col=Group_order, sep=": ", remove=FALSE) %>%
    mutate(Dataset="Fraser")

## Combine datasets
otu_bytype_all <- rbind(otu_bytype.B16_all, otu_bytype.F16_all)

filtColors <- colors()[-grep("white|gr(a|e)y", colors())]
listClasses <- unique(c(otu_bytype.B16_all$Class, otu_bytype.F16_all$Class))
Class_colors <- sample(filtColors, size=length(listClasses))
names(Class_colors) <- listClasses

ggsave(filename = paste0(output,"/Brackonlytaxonomy_class.pdf"), width=10, height=6
       , otu_bytype_all %>%
           ggplot() + geom_col(aes(x=Salinity_round, y=relAbund, fill=Class)) +ylab("Relative Abundance")+xlab("Salinity")+
           scale_fill_manual(values=Class_colors) + 
           facet_grid(Dataset~.)
)
# ggsave(filename = "BrackishOTUs_expanded/Baltic_brackonlytaxonomy_class.pdf", width=10, height=6
#        , otu_bytype.B16_all %>%
#            ggplot() + geom_col(aes(x=Salinity_round, y=relAbund, fill=Class)) +ylab("Relative Abundance")+xlab("Salinity")+
#            scale_fill_manual(values=Class_colors)
# )
# ggsave(filename = "BrackishOTUs_expanded/Fraser_brackonlytaxonomy_class.pdf", width=10, height=6
#        , otu_bytype.F16_all %>%
#            ggplot() + geom_col(aes(x=Salinity_round, y=relAbund, fill=Class)) +ylab("Relative Abundance")+xlab("Salinity")+
#            scale_fill_manual(values=Class_colors)
# )

listOrder <- unique(c(otu_bytype.B16_all$Group_order, otu_bytype.F16_all$Group_order))
Order_colors <- sample(filtColors, size=length(listOrder))
names(Order_colors) <- listOrder

ggsave(filename = paste0(output,"/Brackonlytaxonomy_order.pdf"), width=20, height=6
       , otu_bytype_all %>%
           ggplot() + geom_col(aes(x=Salinity_round, y=relAbund, fill=Group_order)) +ylab("Relative Abundance")+xlab("Salinity")+
           scale_fill_manual(values=Order_colors) + 
           facet_grid(Dataset~.)
)


## Fam
listFamily <- unique(c(otu_bytype.B16_all$Group, otu_bytype.F16_all$Group))
Family_colors <- sample(filtColors, size=length(listFamily))
names(Family_colors) <- listFamily

ggsave(filename = paste0(output,"/Brackonlytaxonomy_Family.pdf"), width=30, height=6
       , otu_bytype_all %>%
           ggplot() + geom_col(aes(x=Salinity_round, y=relAbund, fill=Group)) +ylab("Relative Abundance")+xlab("Salinity")+
           scale_fill_manual(values=Family_colors) + 
           facet_grid(Dataset~.)
)

# Fam; collapsed a bit
otu_bytype_all_collapsed <- otu_bytype_all %>%
    unite(Class, Order, Family, sep=": ", col="Group", remove=FALSE) %>%
    group_by(Group, Order, Salinity_round, Dataset) %>%
    summarize(reads=sum(reads,na.rm=TRUE)
              , relAbund=sum(relAbund,na.rm=TRUE)) 

abundantGroups <- otu_bytype_all_collapsed %>%
    group_by(Group) %>%
    summarize(maxRelAbund=max(relAbund, na.rm=TRUE)) %>%
    filter(maxRelAbund>0.03) %>% select(Group) %>% pull()
# restGroups <- otu_bytype_all_collapsed %>%
#     group_by(Group) %>%
#     summarize(maxRelAbund=max(relAbund, na.rm=TRUE)) %>%
#     filter(maxRelAbund<=0.03) %>% select(Group) %>% pull()
## Fam
set.seed(4234)
abundFamily_colors <- sample(filtColors, size=length(abundantGroups))
names(abundFamily_colors) <- abundantGroups
# lowFamily_colors <- rep(NA, length(restGroups))
# names(lowFamily_colors) <- restGroups
allFamily_colors <- c(abundFamily_colors,`_Other`="grey")

# Filter out low-abund 
otu_bytype_all_collapsed_filt <- otu_bytype_all_collapsed %>%
    mutate(highAbund=ifelse(Group %in% abundantGroups, TRUE, FALSE)
           , Taxonomy=ifelse(highAbund, Group, "_Other")) %>%
    group_by(Taxonomy, Salinity_round, Dataset) %>%
    summarize(relAbund=sum(relAbund, na.rm=TRUE))
    

ggsave(filename = paste0(output,"/Brackonlytaxonomy_Family.pdf"), width=15, height=6
       , otu_bytype_all_collapsed_filt %>%
           ggplot() + geom_col(aes(x=Salinity_round, y=relAbund, fill=Taxonomy)) +ylab("Relative Abundance")+xlab("Salinity")+
           scale_fill_manual(values=allFamily_colors) + 
           facet_grid(Dataset~.)
)



#### 18s Brackish ####
### fraser18
fraser18Brack_OTUs <- mb.F18 %>% filter(typeSimple=="Brackish") %>% select(taxa) 
fraser18_brack_tax <- data.frame(fraser18Brack_OTUs, tax.F18[match(fraser18Brack_OTUs$taxa, rownames(tax.F18)),])%>% 
    rownames_to_column("OTU")%>%
    mutate(Class=gsub("D.*__","",Class)
           ,Family=gsub("D.*__","",Family) ) %>% 
    unite(Class, Family, col=TaxLabel, sep=": ", remove=FALSE) 

otu_bytype.F18_all <- otu.F18.adj %>%
    mutate(Type=mb.F18[match(rownames(otu.F18.adj), mb.F18$taxa),"typeSimple"]
           , taxa = fraser18_brack_tax[match(rownames(otu.F18.adj), fraser18_brack_tax$OTU),"TaxLabel"]
           , Order = fraser18_brack_tax[match(rownames(otu.F18.adj), fraser18_brack_tax$OTU),"Order"]
           , Group=ifelse(is.na(taxa),Type,taxa)) %>%
    gather(-c(Type,taxa,Group,Order), key=SampleID, value=Count) %>%
    group_by(Group, Order, SampleID) %>%
    summarize(reads=sum(Count, na.rm = TRUE)) %>%
    mutate(Salinity=mf.F18[match(SampleID, mf.F18$X.SampleID),"SalinityEnviron"]) %>%
    mutate(Salinity_round=round(as.numeric(Salinity))) %>%
    ungroup() %>% group_by(Salinity_round) %>%
    mutate(relAbund=reads/sum(reads)) %>% ungroup() %>%
    filter(!(Group %in% c("Fresh","Marine","Brackish","Unclassified"))) %>%
    arrange(Salinity) %>% mutate(SampleID = factor(SampleID, levels=unique(SampleID))) %>%
    separate(Group, into=c("Class","Family"), sep=": ", remove=FALSE) %>%
    mutate(Order=gsub("D.*__","",Order)) %>%
    unite(Class, Order, col=Group_order, sep=": ", remove=FALSE) %>%
    mutate(Dataset="fraser18")

listClasses <- unique(c(otu_bytype.F18_all$Class))
Class_colors <- sample(filtColors, size=length(listClasses))
names(Class_colors) <- listClasses

ggsave(filename = paste0(output,"/Fraser18_Brackonlytaxonomy_class.pdf"), width=10, height=6
       , otu_bytype.F18_all %>%
           ggplot() + geom_col(aes(x=Salinity_round, y=relAbund, fill=Class)) +ylab("Relative Abundance")+xlab("Salinity")+
           scale_fill_manual(values=Class_colors) 
)

#### Not shared B16 and F16 #####
otu_bytype_all%>%
    filter(Dataset=="Baltic", Group %in% NOTsharedTax.B16$TaxLabel)
## Fam
set.seed(4234)
Family_colors <- sample(filtColors, size=length(unique(otu_bytype_all$Class)))
names(Family_colors) <- unique(otu_bytype_all$Class)
# # lowFamily_colors <- rep(NA, length(restGroups))
# # names(lowFamily_colors) <- restGroups
# allFamily_colors <- c(abundFamily_colors,`_Other`="grey")
# 
# # Filter out low-abund 
# otu_bytype_all_collapsed_filt <- otu_bytype_all_collapsed %>%
#     mutate(highAbund=ifelse(Group %in% abundantGroups, TRUE, FALSE)
#            , Taxonomy=ifelse(highAbund, Group, "_Other")) %>%
#     group_by(Taxonomy, Salinity_round, Dataset) %>%
#     summarize(relAbund=sum(relAbund, na.rm=TRUE))


ggsave(filename = paste0(output,"/Baltic_NOTSHARED_Brackonlytaxonomy_Class.pdf"), width=15, height=6
       , otu_bytype_all %>%
           ggplot() + geom_col(aes(x=Salinity_round, y=relAbund, fill=Class)) +ylab("Relative Abundance")+xlab("Salinity")+
           scale_fill_manual(values=Family_colors) 
)

