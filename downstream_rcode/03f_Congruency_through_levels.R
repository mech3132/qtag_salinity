#!/bin/bash

library(tidyverse)
library(vegan) # for diversity func
library(lme4) # glmer
library(car) # Anova


#### CONGRUENCY THROUGH LEVELS ####
# This calculates what proportion of each taxonomy level are "pure" (aka how robustly are clades either fresh or marine)
invLogit <- function(x) { return(exp(x)/(1+exp(x)))}

dirNames <- unlist(read.delim("./OUTPUT/dirNames.txt", header=FALSE))
output <- "./OUTPUT/03f_CongruencyThroughLevels"
dir.create(output)
allPEI <- data.frame()
dataset <- c("B16","F16","F18")
for ( d in dataset) {
    if (d == "B16") {
        print("STARTING BALTIC DATASET")
        #Baltic16
        mbpwd <- paste0(dirNames[1],'/modelBoundaries_type.txt')
        otupwd <- paste0(dirNames[1],'/OTUTableText.txt')
        uniqueotupwd <- paste0(dirNames[1],'/OTUs.txt')
        # output <- paste0(newDirName,"/BALTIC16")
    } else if (d == "F16") {
        print("STARTING FRASER 16 DATASET")
        # treepwd <- 'tree_F16_filt.tre'
        #Fraser16
        mbpwd <- paste0(dirNames[2],'/modelBoundaries_type.txt')
        otupwd <- paste0(dirNames[2],'/OTUTableText.txt')
        uniqueotupwd <- paste0(dirNames[2],'/OTUs.txt')
        # output <- paste0(newDirName,"/FRASER16")
    } else if (d == "F18") {
        print("STARTING FRASER 18 DATASET")
        #Fraser18
        mbpwd <- paste0(dirNames[3],'/modelBoundaries_type.txt')
        otupwd <- paste0(dirNames[3],'/OTUTableText.txt')
        uniqueotupwd <- paste0(dirNames[3],'/OTUs.txt')
        # output <- paste0(newDirName,"/FRASER18")
    } 
    ######## START ########
    # dir.create(output)
    mb <- read.delim(file=paste0(mbpwd))
    otu <- read.delim(file=paste0(otupwd), skip=1, row.names=1)
    uniqueotu <- read.delim(file=uniqueotupwd, header=FALSE)
    
    ####################################
    # set new dir

    taxa <- cbind(rownames(otu), as.character(otu[,ncol(otu)]))
    otu <- otu[,-ncol(otu)]
    
    # Filter OTU table to get rid of all OTUs not in uniqueotu; should automatically filter out low abund OTU
    otu.filt <- otu[match(uniqueotu$V1, rownames(otu)),]
    # Get rid of NAs
    otu.filt <- otu.filt[-c(rownames(otu.filt)=="NA"),]
    
    # For 18s, make sure taxa is character not int
    mb <- mb %>% mutate(taxa=as.character(taxa))
    sep_temp <- ifelse(d=="F18","; __", "; ")
    allOTU_taxa_mb <- taxa %>%
        as_tibble() %>%
        mutate(V2=gsub("D_[0-9]__","",V2)) %>%
        separate(V2, into=c("D1","D2","D3","D4","D5","D6","D7"), sep=sep_temp) %>%
        rename(taxa=V1) %>%
        left_join(as_tibble(mb[,c("taxa","typeSimple")])) %>%
        filter(typeSimple%in%c("freshRestricted","marineRestricted")) 
    
    
    for ( lvl in 6:3) {
        group1 <- paste0("D",seq(1:lvl))
        allOTU_temp <- allOTU_taxa_mb %>%
            group_by_at(group1) %>%
            summarize(nFresh=sum(typeSimple=="freshRestricted"), nMarine=sum(typeSimple=="marineRestricted")) %>%
            ungroup() %>% 
            unite(paste0(group1), col="taxa_string", remove=FALSE, sep="; ") %>%
            filter(!(is.na(get(paste0("D",lvl))) | grepl("uncultured",get(paste0("D",lvl))) |  grepl("metagenome",get(paste0("D",lvl))))) %>%
            group_by(taxa_string) %>% 
            filter(sum(c(nFresh,nMarine))>0) %>%
            mutate(SH = diversity(c(nFresh, nMarine), index="shannon")
                   , maxSH = log(sum(nFresh, nMarine))
                   , ntotal = nFresh+nMarine
                   , PEI = ifelse(is.infinite(maxSH), NA, ifelse(SH==0,0, SH/maxSH))) %>%
            select(taxa_string, nFresh, nMarine, SH ,maxSH, ntotal, PEI)
        
        allPEI <- rbind(allPEI,data.frame(allOTU_temp, level=lvl, dataset=d))
        
    }
    
    print("ready to start next dataset")
    
}

allPEI <- allPEI %>%
    mutate(PEI_bin = ifelse(PEI>0,1,0)
           , minority = ifelse(nFresh<nMarine, nFresh, nMarine)
           , majority = ifelse(nFresh<nMarine, nMarine, nFresh)
           , pureClade = ifelse(minority==0,1,0)
           , propMin = minority/ntotal) %>%
    mutate(Eukaryote = ifelse(dataset!="F18",1,0)
           ,factoredLevel = paste0("D",level)
           ,factoredLevel = factor(factoredLevel, levels=c("D6","D5","D4","D3")))

# 
# allPEI %>% 
#     group_by(dataset, level) %>%
#     summarize(wm=weighted.mean(PEI, w=maxSH)
#               , wv=var(rep.int(PEI, times=ntotal))
#               , wse = sd(rep.int(PEI, times=ntotal)/sqrt(sum(ntotal)))) %>%
#     ggplot() + geom_point(aes(x=level, y=wm, col=dataset), cex=3) +
#     geom_line(aes(x=level, y=wm, group=dataset, col=dataset)) +
#     geom_ribbon(aes(x=level, ymin=wm-wse, ymax=wm+wse, group=dataset, fill=dataset), alpha=0.2)
# 
# 
# # Number of pure clades
# allPEI %>%
#     group_by(dataset,factoredLevel) %>%
#     summarize(propPure=sum(pureClade)/n()
#               , CI2.5 = qbinom(c(0.025), size=n(), prob=propPure )/n()
#               , CI97.5 = qbinom(c(0.975), size=n(), prob=propPure )/n()) %>%
#     ggplot(aes(x=factoredLevel, y=propPure)) + geom_point() +
#     geom_segment(aes(x=factoredLevel, xend=factoredLevel, y=CI2.5, yend=CI97.5)) +
#     facet_wrap(.~dataset )
# 
# # Size of clades
# allPEI %>%
#     ggplot(aes(x=factoredLevel, y=log(ntotal))) + geom_boxplot() +
#     facet_wrap(.~dataset )

### Just binomial
allPEI_fonly <- allPEI %>%
    filter(dataset !="B16")

# glm_bin <- glm(pureClade ~ factoredLevel + factoredLevel:Eukaryote + dataset
#                    , data=allPEI
#                    , family = binomial)
# summary(glm_bin)

glm_bin_16 <- glm(pureClade ~ factoredLevel*Eukaryote 
                   , data=allPEI_fonly
                   , family = binomial)

# summary(glm_bin_16)
capture.output(summary(glm_bin_16), file=paste0(output,"/glm_bin_16.txt"))
# Get variance
coefvar_bin <- diag(vcov(glm_bin_16))
# Get coefficient means 
coefmean_bin <- glm_bin_16$coefficients

var_bin_adj <- coefvar_bin %>% t() %>%as_tibble() %>%
    mutate(Bact_D6 =`(Intercept)`
           , Bact_D5 = `(Intercept)` + factoredLevelD5
           , Bact_D4 = `(Intercept)` + factoredLevelD4
           , Bact_D3 = `(Intercept)` + factoredLevelD3
           , Euk_D6 = `(Intercept)` + Eukaryote
           , Euk_D5 = `(Intercept)` + factoredLevelD5 + Eukaryote + `factoredLevelD5:Eukaryote`
           , Euk_D4 = `(Intercept)` + factoredLevelD4 + Eukaryote + `factoredLevelD4:Eukaryote`
           , Euk_D3 = `(Intercept)` + factoredLevelD3 + Eukaryote + `factoredLevelD3:Eukaryote`) %>%
    select(Bact_D6, Bact_D5, Bact_D4, Bact_D3, Euk_D6, Euk_D5, Euk_D4, Euk_D3)
glm_bin_16_summary <- coefmean_bin %>% t() %>%as_tibble() %>%
        mutate(Bact_D6 =`(Intercept)`
               , Bact_D5 = `(Intercept)` + factoredLevelD5
               , Bact_D4 = `(Intercept)` + factoredLevelD4
               , Bact_D3 = `(Intercept)` + factoredLevelD3
               , Euk_D6 = `(Intercept)` + Eukaryote
               , Euk_D5 = `(Intercept)` + factoredLevelD5 + Eukaryote + `factoredLevelD5:Eukaryote`
               , Euk_D4 = `(Intercept)` + factoredLevelD4 + Eukaryote + `factoredLevelD4:Eukaryote`
               , Euk_D3 = `(Intercept)` + factoredLevelD3 + Eukaryote + `factoredLevelD3:Eukaryote`) %>%
        select(Bact_D6, Bact_D5, Bact_D4, Bact_D3, Euk_D6, Euk_D5, Euk_D4, Euk_D3) %>%
    gather(key=Parameter, value=E) %>%
    mutate(Variance=t(var_bin_adj)[,1]
           , sd = sqrt(Variance)) %>%
    mutate(E_adj = invLogit(E)
           , lwrsd = invLogit(E-sd)
           , upprsd = invLogit(E+sd)
           , CI2.5_adj = invLogit(E-sd*2)
           , CI97.5_adj = invLogit(E+sd*2))


### Superimposing summary onto data plots
# Number of pure clades
ggsave(paste0(output,"/proportion_pure_clades.pdf"),
glm_bin_16_summary %>%
    separate(Parameter, into=c("Type","Level"), remove=FALSE) %>%
    mutate(
      # TaxaLevel= c("Genus","Family","Order","Class")[match(Level, c("D6","D5","D4","D3"))]
           # , TaxaLevel=factor(TaxaLevel, levels=c("Genus","Family","Order","Class"))
            Type = ifelse(Type=="Bact","Bacteria","Eukaryote")) %>%
    ggplot(aes(x=Level, y=E_adj)) + geom_point() + 
    geom_segment(aes(xend=Level, y=lwrsd, yend=upprsd)) +
    facet_wrap(.~Type) + ylab("Proportion of fresh/marine exclusive clades") + xlab("Taxonomic Level")
)

# Simulate dataset-- yields same result as just using sd from glm fit
sim_dat <- data.frame(Intercept = rnorm(n=1000, mean=coefmean_bin[1], sd = sqrt(coefvar_bin[1]))
      , D5 = rnorm(n=1000, mean=coefmean_bin[2], sd = sqrt(coefvar_bin[2]))
      , D4 = rnorm(n=1000, mean=coefmean_bin[3], sd = sqrt(coefvar_bin[3]))
      , D3 = rnorm(n=1000, mean=coefmean_bin[4], sd = sqrt(coefvar_bin[4]))
      , Eukaryote = rnorm(n=1000, mean=coefmean_bin[5], sd = sqrt(coefvar_bin[5]))
      , Eukaryote_D5 = rnorm(n=1000, mean=coefmean_bin[6], sd = sqrt(coefvar_bin[6]))
      , Eukaryote_D4 = rnorm(n=1000, mean=coefmean_bin[7], sd = sqrt(coefvar_bin[7]))
      , Eukaryote_D3 = rnorm(n=1000, mean=coefmean_bin[8], sd = sqrt(coefvar_bin[8]))
)
# Checking significance
glm_bin_16_sig <- sim_dat %>%
    mutate(Bact_D6 =Intercept
           , Bact_D5 = Intercept + D5
           , Bact_D4 = Intercept + D4
           , Bact_D3 = Intercept + D3
           , Euk_D6 = Intercept + Eukaryote
           , Euk_D5 = Intercept + D5 + Eukaryote + Eukaryote_D5
           , Euk_D4 = Intercept + D4 + Eukaryote + Eukaryote_D4
           , Euk_D3 = Intercept + D3 + Eukaryote + Eukaryote_D3) %>%
    select(Bact_D6, Bact_D5, Bact_D4, Bact_D3, Euk_D6, Euk_D5, Euk_D4, Euk_D3) %>%
    gather(key=Parameter, value=Coefficient) %>%
    group_by(Parameter) %>%
    summarize(E=mean(Coefficient), sd = sd(Coefficient), CI2.5 = quantile(Coefficient,prob=0.025), CI97.5 = quantile(Coefficient,prob=0.975)) %>% ungroup() 
capture.output(Anova(glm_bin_16, type = 3), file = paste0(output,"/anova_glm_bin_16.txt"))

#### Try doing a model to check even-ness

allPEI_fonly_nozeros <- allPEI_fonly %>%
  filter(propMin>0)

glm_bin2_16 <- glm(cbind(minority,ntotal) ~ factoredLevel*Eukaryote
                   # , offset = log(ntotal)
                  , data=allPEI_fonly_nozeros
                  , family = binomial)

# summary(glm_bin2_16)
capture.output(summary(glm_bin2_16), file=paste0(output,"/glm_bin2_16.txt"))

# Get variance
coefvar_bin2 <- diag(vcov(glm_bin2_16))
# Get coefficient means 
coefmean_bin2 <- glm_bin2_16$coefficients

var_bin2_adj <- coefvar_bin2 %>% t() %>%as_tibble() %>%
  mutate(Bact_D6 =`(Intercept)`
         , Bact_D5 = `(Intercept)` + factoredLevelD5
         , Bact_D4 = `(Intercept)` + factoredLevelD4
         , Bact_D3 = `(Intercept)` + factoredLevelD3
         , Euk_D6 = `(Intercept)` + Eukaryote
         , Euk_D5 = `(Intercept)` + factoredLevelD5 + Eukaryote + `factoredLevelD5:Eukaryote`
         , Euk_D4 = `(Intercept)` + factoredLevelD4 + Eukaryote + `factoredLevelD4:Eukaryote`
         , Euk_D3 = `(Intercept)` + factoredLevelD3 + Eukaryote + `factoredLevelD3:Eukaryote`) %>%
  select(Bact_D6, Bact_D5, Bact_D4, Bact_D3, Euk_D6, Euk_D5, Euk_D4, Euk_D3)
glm_bin2_16_summary <- coefmean_bin2 %>% t() %>%as_tibble() %>%
  mutate(Bact_D6 =`(Intercept)`
         , Bact_D5 = `(Intercept)` + factoredLevelD5
         , Bact_D4 = `(Intercept)` + factoredLevelD4
         , Bact_D3 = `(Intercept)` + factoredLevelD3
         , Euk_D6 = `(Intercept)` + Eukaryote
         , Euk_D5 = `(Intercept)` + factoredLevelD5 + Eukaryote + `factoredLevelD5:Eukaryote`
         , Euk_D4 = `(Intercept)` + factoredLevelD4 + Eukaryote + `factoredLevelD4:Eukaryote`
         , Euk_D3 = `(Intercept)` + factoredLevelD3 + Eukaryote + `factoredLevelD3:Eukaryote`) %>%
  select(Bact_D6, Bact_D5, Bact_D4, Bact_D3, Euk_D6, Euk_D5, Euk_D4, Euk_D3) %>%
  gather(key=Parameter, value=E) %>%
  mutate(Variance=t(var_bin2_adj)[,1]
         , sd = sqrt(Variance)) %>%
  mutate(E_adj = invLogit(E)
         , lwrsd = invLogit(E-sd)
         , upprsd = invLogit(E+sd)
         , CI2.5_adj = invLogit(E-sd*2)
         , CI97.5_adj = invLogit(E+sd*2))


### Superimposing summary onto data plots
# Number of pure clades
ggsave(paste0(output,"/proportion_minority_clades.pdf"),
       glm_bin2_16_summary %>%
         separate(Parameter, into=c("Type","Level"), remove=FALSE) %>%
         mutate(
           # TaxaLevel= c("Genus","Family","Order","Class")[match(Level, c("D6","D5","D4","D3"))]
                # , TaxaLevel=factor(TaxaLevel, levels=c("Genus","Family","Order","Class"))
                 Type = ifelse(Type=="Bact","Bacteria","Eukaryote")) %>%
         ggplot(aes(x=Level, y=E_adj)) + geom_point() + 
         geom_segment(aes(xend=Level, y=lwrsd, yend=upprsd)) +
         facet_wrap(.~Type) + ylab("Expected proportion of minority within mixed clades") + xlab("Taxonomic Level")
)

