#!/bin/bash

#### Get unique list of OTUs
dirNames <- unlist(read.delim("./OUTPUT/dirNames.txt", header = FALSE))

baltic <- paste0(dirNames[1],"/OTUs.txt")
fraser <- paste0(dirNames[2],"/OTUs.txt")

botus <- read.delim(baltic, header=FALSE, stringsAsFactors = FALSE)
fotus <- read.delim(fraser, header=FALSE, stringsAsFactors = FALSE)

allotus <- unique(c(botus$V1,fotus$V1))

write.table(allotus, file="./OUTPUT/allOTUs_combo_BF16.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
