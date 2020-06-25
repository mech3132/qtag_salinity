#!bin/bash


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
localMax <- function(x) {
  return(which(diff(sign(diff(x)))==-2)+1)
}
localMin <- function(x) {
  return(which(diff(sign(diff(x)))==+2)+1)
}
