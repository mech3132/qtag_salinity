#!/bin/bash
### Scan for contaminants ###
###################### DESCRPTION OF SCRIPT #######################
# This script scans the whole OTU table using the control sequences and figures out if there are contaminants that we should remove.
# Typical characteristics of a contaminant that was introduced from the environment are:
  # It is present in the blank or control sequences (duh!)
  # It is relatively abundant in low-read samples and relatively spares in high-read samples

# Note that this does NOT identify contaminants that may have been introduced between wells. It only looks for contaminants that have been introduced to "everything".

# This script takes the raw OTU table and plots suspected contaminants but:
# it does NOT remove OTUs from the table; it simply helps you identify potential contaminants so you can further investigate them.

library(optparse)

option_list <- list(
  make_option(c("-i", "--otu_table"), default = "EMPTY", help="OTU table (raw) in text format. MUST HAVE TAXONOMY AS LAST COLUMN"),
  make_option(c("-m", "--mapping_file"), default = "EMPTY", help="Mapping file or metadata. Columns are categories; rows are sample IDs."),
  make_option(c("-b", "--blanks"), help="The column and treatment for blanks. Eg.'Category:treatment'"),
  make_option(c("-s", "--split_table"), help="If you want to split table into multiple tables, indicate column and table to include for each table from mapping file. Usually occurs when OTU table is composed of multiple runs with different controls. Eg. 'Column:table1,table2'. MUST BE INCLUDED. If all samples are from the same run, create a metadata category and make them all belong to one treatment."),
  make_option(c("--filter_low_abund_per_sample"), default=0, help="Prior to scanning, change all OTUs in each sample below this threshold to zero. Default: 0"),
  make_option(c("--filter_low_abund_overall"), default=0, help="Prior to scanning, remove all OTUs whose maximum read count in any sample is less than this threshold. (Note that this is done BEFORE filtering by sample) Default: 0"),
  make_option(c("-o", "--output"), default="contam_scan", help="Folder name for output files"),
  make_option(c("--dashes"), default="FALSE", help="If your sample names have dashes, mark this as TRUE. This is important because R sometimes converts dashes to periods in header names and the mapping file names MUST match OTU table names."),
  make_option(c("--numberstart", default="FALSE", help="If your sampe names start with numbers, mark this as TRUE. This is important because R will put an 'X' in front of header names that start with a number.")),
  make_option(c("--delim", help="Delimiter to use for taxonomy. This will help split taxonomy names."))
  )
opt = parse_args(OptionParser(option_list=option_list))
otu_table = opt$otu_table
mf_fp = opt$mapping_file
blanks = opt$blanks
splittable = opt$split_table
thresh_sample = opt$filter_low_abund_per_sample
thresh_all = opt$filter_low_abund_overall
output = opt$output
dashes = opt$dashes
numberstart = opt$numberstart
delim = opt$delim

# Exit if something is missing
if (otu_table == "EMPTY") {
  print("Please provide OTU Table")
  quit()
}
if (mf_fp == "EMPTY") {
  print("Please provide MF")
  quit()
}

############# FOR TESTING###################
# 
# setwd("/Users/melissachen/Documents/Masters/Project_Environmental/Project_Salinity/Removing_contamined_samples")
# otu_table = "otu_table_wtaxa_16_text.txt"
# mf_fp = "../raw_otu_and_mf/Fraser_16S/MANUAL_INPUT_FILES/Fraser16s_mappingfile_merged.txt"
# blanks = "control:y"
# splittable = "Year:2014,2015"
# thresh_sample = 5
# thresh_all = 10
# output = "TEST"
# dashes = TRUE
# numberstart = FALSE
# delim = ";"

################# READ IN DATA ###############
print("Reading in data...")

# Make output folder
dir.create(output)
setwd(output)

### Check intput pathways for all files
otu_table = paste0("../", otu_table)
mf_fp = paste0("../",mf_fp)

if ( !is.null(grep(".txt", otu_table)) ) {
    file.copy(from=otu_table, to=paste0(getwd(),"/OTU_Table_text_forcontamscan.txt"))
} else {
    print("ERROR: OTU table not recognized as text file.")
    quit()
}

## Find out the line that the OTU table starts on

# system("tr -d '\r' < OTU_Table_text_forcontamscan.txt > OTU_Table_text_forcontamscan_convert.txt")
# system("echo '\n' >> OTU_Table_text_forcontamscan_convert.txt")
# file.remove("OTU_Table_text_forcontamscan.txt")
fileTemp <- file("OTU_table_text_forcontamscan.txt", "r")
lineStart <- 0
notfound <-  TRUE
while ( notfound ) {
    print(lineStart)
    lineTemp <- grep("#OTU ID", readLines(fileTemp, n=1))
    if ( (length(lineTemp)>0)  ) {
        # print("REACH HERE")
        notfound <- FALSE
    } else {
        lineStart <- lineStart + 1
    }
}
close(fileTemp)

print("Reading in OTU table...")
otu <- read.delim("OTU_Table_text_forcontamscan.txt", header=TRUE, row.names=1, skip=lineStart,stringsAsFactors = FALSE)
mf <- read.delim(paste0(mf_fp), header=TRUE, row.names=1, stringsAsFactors = FALSE, na.strings = "")

# Get blank names
blanks_split <- unlist(strsplit(blanks, split = c(":")))
blanks_split <- c(blanks_split[1], strsplit(blanks_split[2], ","))

# Check that blanks split is in correct format
if ( (length(blanks_split[[1]]) == 1) & is.character(blanks_split[[1]]) ) {
    if ( !(blanks_split[[1]] %in% colnames(mf)) ) {
        print("ERROR: BLANKS flag is not in correct format. Must be category:treatment. Category not found in metadata file.")
        quit()
    } 
} else {
    print("ERROR: BLANKS flag is not in correct format. Must be category:treatment. Category not provided.")
}
if ( !is.na(blanks_split[[2]]) & (length(blanks_split[[1]]) == 1) ) {
    if ( !(blanks_split[[2]] %in% mf[,blanks_split[[1]]]) ) {
        print("ERROR: BLANKS flag is not in correct format. Must be category:treatment. Treatment not found in metadata file.")
        quit()
    } 
} else {
    print("ERROR: BLANKS flag is not in correct format. Must be category:treatment. Treatment not provided.")
}

# Get items to split table by
table_split <- unlist(strsplit(splittable, split = c(":")))
table_split <- c(table_split[1], strsplit(table_split[2], ","))

# Check that split table is in correct format.
if ( (length(table_split[[1]]) == 1) & is.character(table_split[[1]]) ) {
    if ( !(table_split[[1]] %in% colnames(mf)) ) {
        print("ERROR: SPLITTABLE flag is not in correct format. Must be category:treatment. Category not found in metadata file.")
        quit()
    } 
} else {
    print("ERROR: SPLITTABLE flag is not in correct format. Must be category:treatment. Category not provided.")
    quit()
}

if ( !any(is.na(table_split[[2]])) ) {
    if ( !any(table_split[[2]] %in% mf[,table_split[[1]]]) ) {
        print("ERROR: SPLITTABLE flag is not in correct format. Must be category:treatment. Treatment not found in metadata file.")
        quit()
    } 
} else {
    print("ERROR: SPLITTABLE flag is not in correct format. Must be category:treatment. Treatment not provided.")
    quite()
}

# Now lastly, check that all split tables have a control. If not, send out warning.
for ( i in table_split[[2]] ) {
    if ( ! any((mf[,table_split[[1]]] == i) & (mf[,blanks_split[[1]]] == blanks_split[[2]])) ) {
        print(paste0("WARNING: The split table ",i," does not have a control. You may get no output because there are no contaminants to scan."))
    }
}


###################### SETTING UP TABLES #########################
print("Setting up tables...")
# Get rid of taxonomy
## Checking that taxonomy is last column
if ( colnames(otu)[ncol(otu)] == "taxonomy" ) {
    taxa <- data.frame(otu[,ncol(otu)])
    otu <- otu[,-ncol(otu)]
    
} else {
    taxa <- rep("TAXA_NA", nrow(otu))
}
taxa <- data.frame(rownames(otu), taxa)

# Fix header names if applicable
if (numberstart == "TRUE") {
  colnames(otu) <- gsub("^X","",colnames(otu))
}
if (dashes == "TRUE") {
  rownames(mf) <- gsub("-",".", rownames(mf), fixed=TRUE)
}

# Check to make sure mf and otu are the same; first, make sure they are shared

tofilt <- which(rownames(mf) %in% colnames(otu)) 
mf.filt <- mf[tofilt,]

# Now, split up into different studies
otu.bystudy <- list()
mf.bystudy <- list()
for ( study in table_split[[2]] ) {
  mf.bystudy[[study]] <- mf.filt[which(mf.filt[[table_split[[1]]]] == study),]
  otu.bystudy[[study]] <- otu[,match(rownames(mf.bystudy[[study]]), colnames(otu))]
}

###### (1) CHANGE ALL COUNTS WITH LESS THAN 5 PER SAMP TO 0 ######
# print("REACH HERE")
print("DELETING COUNTS WITH LESS THAN 5 PER SAMP")
# print("REACH HERE2")
for ( study in names(otu.bystudy) ) {
    print(paste0("Filtering study ",study))
    
  otu.temp <- otu.bystudy[[study]]
  otu.temp[otu.temp < thresh_sample]<- 0
  
  otu.bystudy[[study]] <- otu.temp
}

###### (2) REMOVE ALL OTUS WITH LESS THAN 10 IN WHOLE OTU TABLE ####
# remove all zeros from otu table
print("DELETEING COUTNS WITH LESS THAN 10 per SAMP")
for ( study in names(otu.bystudy) ) {
  otu.temp <- otu.bystudy[[study]]
  todel <- as.vector(which(rowSums(otu.temp) < thresh_all))
  
  if (length(todel) > 0) {
      otu.bystudy[[study]] <- otu.temp[-todel,]
  }
  
}


###### (3) FILTER MF AND OTU TO MAKE SURE ALL SAMPLES INCLUDED ######
# Check to make sure everything is in everything
for ( study in names(otu.bystudy) ) {
    print("The following logicals should be FALSE to indicate that there are no samples in the data that are not found in metadata")
  # all otu in mf?
  print(any(!(colnames(otu.bystudy[[study]]) %in% rownames(mf.bystudy[[study]]))))
  # all mf in otu?
  print(any(!(rownames(mf.bystudy[[study]]) %in% colnames(otu.bystudy[[study]]))))
}


## Now, begin iterating through studies
for ( study in names(otu.bystudy) ) {
    # study = names(otu.bystudy)[6]
  mf.temp <- mf.bystudy[[study]]
  otu.temp <- otu.bystudy[[study]]
  
  # Get sum of all sample reads, and also filter based on which are actually present in the otu table
  # colnames(otu.temp)
  if ( any(mf.temp[blanks_split[[1]]]==blanks_split[[2]]) ) {
      contam.temp <- otu.temp[,which(mf.temp[blanks_split[[1]]]==blanks_split[[2]])]
      contam.temp <- as.data.frame(contam.temp, row.names=rownames(otu.temp))
  } else {
      print(paste0("NO CONTROLS FOUND IN STUDY ",study,"-- DISCONTINUING LOOP"))
      next
  }

  names.contam.temp <- rownames(contam.temp)[ which(rowSums(contam.temp) >0)]
  otu.temp.filt.nocon <- otu.temp[,!c(mf.temp[blanks_split[[1]]]==blanks_split[[2]])]
  RPS.temp <- colSums(otu.temp.filt.nocon)
  otu.temp.filt <- otu.temp.filt.nocon[names.contam.temp,]
  if (nrow(otu.temp.filt)>0 ) {
      # make list of lm of reads in sample vs size of sample
      all.lm.temp <- list()
      for ( r in 1:nrow(otu.temp.filt)) {
          temp <- lm(as.numeric(otu.temp.filt[r,])/RPS.temp ~ as.vector(RPS.temp))
          m <- coef(temp)[2]
          p <- summary(temp)$coefficients[2,4]
          all.lm.temp[[paste0(rownames(otu.temp.filt)[r])]] <- c(m,p)
      }
      # if the slope is negative AND it's significant, then same as tent contam
      tent.contam.temp <- c()
      not.contam.temp <- c()
      for ( o in names(all.lm.temp) ) {
          # o="JN975971.1.1410"
          temp <- all.lm.temp[[paste0(o)]]
          if ( (temp[2] < 0.05) & (temp[1] < 0) ) {
              tent.contam.temp <- c(tent.contam.temp, o)
          } else {
              not.contam.temp <- c(not.contam.temp, o)
          }
      }

      # now plot contam and get names
      names.contam.temp <- taxa[match(tent.contam.temp,taxa[,1]),2]
      # print(names.contam.temp)
      split.names.temp <- strsplit(as.character(names.contam.temp), split = delim)

      if (length(names.contam.temp) > 0) {
          # get abundance in controls
          contam.tab <- contam.temp[match(tent.contam.temp, rownames(contam.temp)),]
          contam.tab <- as.data.frame(contam.tab, rownames=tent.contam.temp)
          contam.counts <- c()
          for ( r in 1:nrow(contam.tab) ) {
              toprint <- ""
              for ( c in 1:ncol(contam.tab) ) {
                  toprint <- paste0(toprint,", ",names(contam.tab)[c],":",contam.tab[r,c])
              }
              contam.counts[r] <- toprint
          }
          mat.dim <- ceiling(sqrt(length(tent.contam.temp))) # make matrix
          pdf(file = paste0("Possible_contam_",study,".pdf"), width = 10*(mat.dim/3), height = 10*(mat.dim/3))
          par(mfrow=c(mat.dim, mat.dim))
          count <- 1
          toPrintContam <- c()
          for (n in tent.contam.temp) {
              lengthName <- length(split.names.temp[[count]])
              nameToPrint <- paste0(split.names.temp[[count]][lengthName],split.names.temp[[count]][lengthName-1])
              # print(split.names.temp[[count]])
              titleToPrint <- paste0(strtrim(nameToPrint,30), "(",n,")")
              truncTitle <- paste0(strtrim(titleToPrint,37),"...")
              
              plot(as.numeric(otu.temp.filt.nocon[n,])/RPS.temp ~ RPS.temp, ylab="Rel Abund", main=truncTitle, xlab="Reads/sample", sub=paste0(contam.counts[count]), cex.sub=0.5, cex.main=0.5)
              count <- count + 1
              toPrintContam <- c(toPrintContam, paste0(nameToPrint, "(",n,")"))
          }
          dev.off()
          write.table(toPrintContam, file=paste0("Possible_contam",study,".txt"))
      }
      
      # plot non-contam and get names
      names.noncontam.temp <- taxa[match(not.contam.temp,taxa[,1]),2]
      split.names.non.temp <- strsplit(as.character(names.noncontam.temp), split = delim)
      # get abundance in controls
      ncontam.tab <- as.data.frame(contam.temp[not.contam.temp,], row.names =not.contam.temp)
      ncontam.counts <- c()

      for ( r in 1:nrow(ncontam.tab) ) {
          toprint <- ""
          for ( c in 1:ncol(ncontam.tab) ) {
              toprint <- paste0(toprint,", ",names(ncontam.tab)[c],":",ncontam.tab[r,c])
          }
          ncontam.counts[r] <- toprint
      }
      
      mat.dim <- ceiling(sqrt(length(not.contam.temp)))# make matrix
      pdf(file = paste0("Other_contam_",study,".pdf"), width = 10*(mat.dim/3), height = 10*(mat.dim/3))
      par(mfrow=c(mat.dim, mat.dim))
      count <- 1
      toPrintContam <- c()
      for (n in not.contam.temp) {

          lengthName <- length(split.names.non.temp[[count]])
          if ( lengthName > 1 ) {
              nameToPrint <- paste0(split.names.non.temp[[count]][lengthName],split.names.non.temp[[count]][lengthName-1])
              
          } else {
              nameToPrint <- paste0(split.names.non.temp[[count]][lengthName])
          }
          titleToPrint <- paste0(strtrim(nameToPrint,30), "(",n,")")
          truncTitle <- paste0(strtrim(titleToPrint,37),"...")
          plot(as.numeric(otu.temp.filt.nocon[n,])/RPS.temp ~ RPS.temp, ylab="Rel Abund", main=truncTitle, xlab="Reads/sample",sub=paste0(ncontam.counts[count]), cex.sub=0.5, cex.main=0.5)
          count <- count + 1
          toPrintContam <-c(toPrintContam,paste0(nameToPrint, "(",n,")"))
      }
      dev.off()
      write.table(toPrintContam, file=paste0("Other_contam",study,".txt"))
      
  }
  
  
}

system("rm OTU_Table_text_forcontamscan.txt")


