#!bin/bash
# Given an level, names.split list, and modelboundaries file, find out proportion of fresh/marine/brackish.
propType <- function(level,namesFull,mb) {
  # level <- 3
  namesSplittemp <- strsplit(as.character(namesFull),"; ")
  namesSplit <- lapply(namesSplittemp, function(x) gsub("^.*__", "", x))
  names(namesSplit) <- names(namesFull)
  allLevelNames <- list()
  for ( n in 1:length(namesSplit) ) {
    # n = 1922
    # Need to find neighbours that have real assignments
    prevfullname <- NULL
    prevcounter <- 1
    nextfullname <- NULL
    nextcounter <- 1
    # While neighbours are unassigned, try to find the flanking taxonomy names
    while ( ( (length(grep("Unassigned", prevfullname)) > 0) | is.null(prevfullname) ) & ( (n-prevcounter) > 0 ) ) {
      prevt <- names(namesSplit)[n-prevcounter]
      prevfullname <- namesSplit[[prevt]]
      prevcounter <- prevcounter + 1
    }  
    # if we run out of spots and it's unassigned, just make it null
    if ( length(grep("Unassigned", prevfullname))>0 ) {
      prevfullname <- NULL
    }
    # while the enighours are unassigned, try to find the flanking taxonomy names
    while ( ( (length(grep("Unassigned", nextfullname)) > 0) | is.null(nextfullname) ) & ( (n+nextcounter) <= length(namesSplit)) ) {
      nextt <- names(namesSplit)[n+nextcounter]
      nextfullname <- namesSplit[[nextt]]
      nextcounter <- nextcounter + 1
    }
    # if we run out of spots and it's unassigned, just make it null
    if ( length(grep("Unassigned", nextfullname))>0 ) {
      nextfullname <- NULL
    }
    
    # t <- names(namesSplit)[25]
    
    prevtopaste <- ""
    nexttopaste <- ""
    
    # now get the CURRENT taxonomy name
    t <- names(namesSplit)[n]
    fullname <- namesSplit[[t]]
    topaste <- ""
    for ( l in 3:level) {
      topaste <- paste0(topaste, fullname[l],"_")
      prevtopaste <- paste0(prevtopaste, prevfullname[l],"_")
      nexttopaste <- paste0(nexttopaste, nextfullname[l],"_")
      # change the flanking names if they are NULL to each other; usually means you reached the end of the set
      if ( prevtopaste == paste0(rep("_", l-2)) ) {
        prevtopaste <- nexttopaste
      }
      if (nexttopaste == paste0(rep("_",l-2)) ) {
        nexttopaste <- prevtopaste
      }
      # Now, if current name is unassigned, then change its name to the flanking names, assuming they are the same.
      if ( length(grep("Unassigned",topaste))>0 & (prevtopaste == nexttopaste) ) {
        topaste <- prevtopaste
      }
    }
    allLevelNames[[t]] <- topaste
    print( paste0( "Finished ",n," out of ",length(namesSplit)))
  }
  # uniqueLevelNames <- unique(allLevelNames)
  uniqueLevelNames <- uniqueconsec(allLevelNames)
  uniqueLevelNames.subscript <- make.names(uniqueLevelNames, unique = TRUE)
  # listClassComp <- sapply(uniqueLevelNames, function(x) NULL)
  listClassComp <- list()
  remainingallLevelNames <- allLevelNames
  for ( n in 1:length(uniqueLevelNames) ) {
    i <- uniqueLevelNames[n]
    i.subscript <- uniqueLevelNames.subscript[n]
    # go through and get all consequtive ones
    last <- FALSE
    firstdone <- FALSE
    otusForLevel <- c()
    while ( !last ) {
      if ( remainingallLevelNames[1] == i ) {
        firstdone <- TRUE
        otusForLevel <- c(otusForLevel, names(remainingallLevelNames[1]))
        remainingallLevelNames <- remainingallLevelNames[-1]
      } else if ( remainingallLevelNames[1] != i & firstdone ) {
        last <- TRUE
        # print(paste0("finished ", i.subscript))
      }
    }
    # otusForLevel <- names(which(allLevelNames == i))
    classes <- mb[match(otusForLevel, mb$taxa),"typeSimple"]
    tabletemp <- table(classes, exclude = "noclass")/sum(table(classes, exclude = "noclass"))
    if (is.na(sum(tabletemp))) {
      tabletemp[] <- c(0,0,0)
    }
    listClassComp[[i.subscript]] <- tabletemp
  }
  return(listClassComp)
}

tableAllTypes <- function(level,namesFull) {
  namesSplittemp <- strsplit(as.character(namesFull),"; ")
  namesSplit <- lapply(namesSplittemp, function(x) gsub("^.*__", "", x))
  names(namesSplit) <- names(namesFull)
  allLevelNames <- c()
  ######
  
  for ( n in 1:length(namesSplit) ) {
    # n = 1922
    # Need to find neighbours that have real assignments
    prevfullname <- NULL
    prevcounter <- 1
    nextfullname <- NULL
    nextcounter <- 1
    # While neighbours are unassigned, try to find the flanking taxonomy names
    while ( ( (length(grep("Unassigned", prevfullname)) > 0) | is.null(prevfullname) ) & ( (n-prevcounter) > 0 ) ) {
      prevt <- names(namesSplit)[n-prevcounter]
      prevfullname <- namesSplit[[prevt]]
      prevcounter <- prevcounter + 1
    }  
    # if we run out of spots and it's unassigned, just make it null
    if ( length(grep("Unassigned", prevfullname))>0 ) {
      prevfullname <- NULL
    }
    # while the enighours are unassigned, try to find the flanking taxonomy names
    while ( ( (length(grep("Unassigned", nextfullname)) > 0) | is.null(nextfullname) ) & ( (n+nextcounter) <= length(namesSplit)) ) {
      nextt <- names(namesSplit)[n+nextcounter]
      nextfullname <- namesSplit[[nextt]]
      nextcounter <- nextcounter + 1
    }
    # if we run out of spots and it's unassigned, just make it null
    if ( length(grep("Unassigned", nextfullname))>0 ) {
      nextfullname <- NULL
    }
    
    # t <- names(namesSplit)[25]
    
    prevtopaste <- ""
    nexttopaste <- ""
    
    # now get the CURRENT taxonomy name
    t <- names(namesSplit)[n]
    fullname <- namesSplit[[t]]
    topaste <- ""
    for ( l in 3:level) {
      topaste <- paste0(topaste, fullname[l],"_")
      prevtopaste <- paste0(prevtopaste, prevfullname[l],"_")
      nexttopaste <- paste0(nexttopaste, nextfullname[l],"_")
      # change the flanking names if they are NULL to each other; usually means you reached the end of the set
      if ( prevtopaste == paste0(rep("_", l-2)) ) {
        prevtopaste <- nexttopaste
      }
      if (nexttopaste == paste0(rep("_",l-2)) ) {
        nexttopaste <- prevtopaste
      }
      # Now, if current name is unassigned, then change its name to the flanking names, assuming they are the same.
      if ( length(grep("Unassigned",topaste))>0 & (prevtopaste == nexttopaste) ) {
        topaste <- prevtopaste
      }
    }
    allLevelNames[t] <- topaste
    print( paste0( "Finished ",n," out of ",length(namesSplit)))
  }
  
  #####
  # for ( t in names(namesSplit) ) {
  #     fullname <- namesSplit[[t]]
  #     topaste <- ""
  #     for ( l in 3:level) {
  #         topaste <- paste0(topaste, fullname[l],"_")
  #     }
  #     allLevelNames[t] <- topaste
  # }
  # make table with subscripted names
  # uniqueLevelNames <- unique(allLevelNames)
  uniqueLevelNames <- uniqueconsec(allLevelNames)
  uniqueLevelNames.subscript <- make.names(uniqueLevelNames, unique = TRUE)
  # listClassComp <- sapply(uniqueLevelNames, function(x) NULL)
  allLevelNames.subscript <- list()
  remainingallLevelNames <- allLevelNames
  for ( n in 1:length(uniqueLevelNames) ) {
    i <- uniqueLevelNames[n]
    i.subscript <- uniqueLevelNames.subscript[n]
    # go through and get all consequtive ones
    last <- FALSE
    firstdone <- FALSE
    otusForLevel <- c()
    while ( !last & length(remainingallLevelNames)>0 ) {
      if ( remainingallLevelNames[1] == i ) {
        firstdone <- TRUE
        otusForLevel <- c(otusForLevel, names(remainingallLevelNames[1]))
        remainingallLevelNames <- remainingallLevelNames[-1]
      } else if ( remainingallLevelNames[1] != i & firstdone ) {
        last <- TRUE
        # print(paste0("Finished ", i.subscript))
      }
    }
    for ( otu in otusForLevel ) {
      for ( pos in names(allLevelNames) ) {
        if ( otu == pos ) {
          allLevelNames[[pos]] <- i.subscript
        }
      }
    }
    # match(otusForLevel,allLevelNames)
    # classes <- mb[match(otusForLevel, mb$taxa),"typeSimple"]
    # tabletemp <- table(classes, exclude = "noclass")/sum(table(classes, exclude = "noclass"))
    # if (is.na(sum(tabletemp))) {
    #     tabletemp[] <- c(0,0,0)
    # }
    # listClassComp[[i.subscript]] <- tabletemp
  }
  
  return(table(allLevelNames))
}

getColor <- function(listClassComp) {
  color <- sapply(listClassComp, function(x) {
    r <- ifelse(is.na(x['marineRestricted']), 0,x['marineRestricted'])
    b <- ifelse(is.na(x['freshRestricted']), 0,x['freshRestricted'])
    g <- ifelse(is.na(x['brackishRestricted']), 0,x['brackishRestricted'])
    intensity <- max(c(r,b,g), na.rm = TRUE)
    return(rgb(r,g,b,intensity))
  })
  return(color)
}

makebarcolor <- function(colors,tableTypes,levelname) {
  barmatrix <- matrix(ncol=1,nrow=sum(tableTypes),dimnames = list(NULL,levelname))
  allColors <- unlist(sapply(names(colors), function(x) {
    rep(colors[x],tableTypes[x])
  }))
  return(allColors)
}

uniqueconsec <- function(vec) {
  finalvec <- c("begin")
  n <- 1
  for ( i in vec ) {
    if ( i != finalvec[n] ) {
      n <- n+1
      finalvec <- c(finalvec,i)
    } 
  }
  return(finalvec[-1])
}

# getTaxaNameAtLevel <- function(fullname, )
