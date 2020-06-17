########################################################################################################################
# Convenience functions for pulling, plotting, and comparing the FIBSEM and FAFB data
########################################################################################################################


### Load the FAFB data in a dataframe
FAFBLoad <- function(region){
  #' Get data frame of FAFB data in a given region
  #' @return Returns a data frame of the FAFB data in the specified region
  #' @param region: The chosen ROI
  
  # Load the csv file
  FAFBDat = read.csv(paste0(
    'connectivity_matrices\\matrix_',region,'.txt')
    , header=TRUE, sep ="")
  
  # Get out the names of the neurons within the dataset and rename where appropriate
  nrons = rownames(FAFBDat)
  nrons <- gsub("D7-","Delta7-",nrons)
  nrons[which(nrons == "Ra")] <- "R4d-a"
  nrons[which(nrons == "Rb")] <- "R4d-b"
  nrons[which(nrons == "Rc")] <- "R4d-c"
  nrons[which(nrons == "Rd")] <- "R4d-d"
  nrons[which(nrons == "Re")] <- "R4d-e"
  nrons[which(nrons == "Rf")] <- "R4d-f"
  nrons[which(nrons == "Rg")] <- "R4d-g"
  
  # Create nameids and partnerids based off of the constituent neurons
  nameids = c()
  partnerids = c()
  for (rpt in (1:length(nrons))){
    nameids = append(nameids,nrons)
    for (r in (1:length(nrons))){
      partnerids = append(partnerids,nrons[rpt])
    }
  }
  
  # Populate a dataframe from the weight matrix
  colNames = names(FAFBDat)
  colNames <- gsub("\\.","-",colNames)
  colNames <- gsub("D7-","Delta7-",colNames)
  colNames[which(colNames == "Ra")] <- "R4d-a"
  colNames[which(colNames == "Rb")] <- "R4d-b"
  colNames[which(colNames == "Rc")] <- "R4d-c"
  colNames[which(colNames == "Rd")] <- "R4d-d"
  colNames[which(colNames == "Re")] <- "R4d-e"
  colNames[which(colNames == "Rf")] <- "R4d-f"
  colNames[which(colNames == "Rg")] <- "R4d-g"
  FAFB_df = data.frame(nameid.from = nameids,nameid.to = partnerids, ROIweight = integer(length(nrons)))
  for (i in 1:length(colNames)){
    for (j in 1:length(nrons)){
      FAFB_df$ROIweight[which(FAFB_df$nameid.to == colNames[i] & FAFB_df$nameid.from == nrons[j])] <- FAFBDat[j,i]
    }
  }
  
  return(FAFB_df)
}

### Neuprint search
getBodyIdsForList = function (neuronList,prefix="",postfix=".*",...){
  #' Get one dataframe of bodyIDs for all search strings in neuronList
  #' @param neuronList: A list of search strings to be passed.
  #' @param prefix: String to be added before each query (default "")
  #' @param postfix: String to be added after each query (default ".*")
  #' @param ...: Parameters to be passed to neuprint_search. Note that meta=FALSE won't work for now.
  #' @return A data frame of metadata for the results of all the queries
  #' @examples
  #' \dontrun{
  #' # Both will return the same
  #' getBodyIdsForList(c("PFL1","PFL2"))
  #' getBodyIdsForList(c("PFL1","PFL2"),postfix="",field="type")
  #' }
  
  neuronList <-  paste0(prefix,neuronList,postfix)
  bodiesList <- lapply(neuronList,neuprint_search,...)
  return(bind_rows(bodiesList))
}

### Load the FIBSEM data
FIBSEMLoad <- function(neuronList,region){
  #' Get data frame of FIBSEM data in a given region
  #' @return Returns a data frame of the FIBSEM data in the specified region
  #' @param region: The chosen ROI
  #' @param neuronList: A list of neuron names to pull
  
  # Get bodyids associated with the neurons
  bIDs <- getBodyIdsForList(neuronList)$bodyid
  bNames <- neuprint_get_meta(bIDs)$name
  
  # Pull the data for the given set of neurons
  FIBSEMDat = getConnectionTable(bIDs, "POST", region) %>%
    filter(to %in% bIDs) %>% select(name.from,from,name.to,to,ROIweight)
  
  # Fill in zeros where no values exist
  for (f in 1:length(bIDs)){
    for (t in 1:length(bIDs)){
      if (length(which(FIBSEMDat$from == bIDs[f] & FIBSEMDat$to == bIDs[t])) == 0){
        FIBSEMDat <- rbind(FIBSEMDat,
                           data.frame(name.from = bNames[f], from = bIDs[f],
                                      name.to = bNames[t], to = bIDs[t],
                                      ROIweight = 0))
      }
    }
  }
  
  # Establish unique names for each neuron that are consistent with the FAFB names
  FIBSEMDat$nameid.from <- paste(as.character(FIBSEMDat$name.from), as.character(FIBSEMDat$from), sep = "_")
  FIBSEMDat$nameid.to <- paste(as.character(FIBSEMDat$name.to), as.character(FIBSEMDat$to), sep = "_")
  
  
  oldNameIds <- c(FIBSEMDat$nameid.from,FIBSEMDat$nameid.to) %>% unique() 
  newNameIds <- FIBSEMRename(oldNameIds)
  FIBSEMDat$nameid.from <- sapply(FIBSEMDat$nameid.from,
                                  function(n){newNameIds[which(oldNameIds == n)] })
  FIBSEMDat$nameid.to <- sapply(FIBSEMDat$nameid.to,
                                function(n){newNameIds[which(oldNameIds == n)] })
  FIBSEMDat <- FAFBSort(FIBSEMDat,region)
  
  return(FIBSEMDat %>% select(nameid.from, nameid.to, ROIweight))
}

### Rename FIBSEM neurons
FIBSEMRename <- function(nameids){
  #' Rename FIBSEM neurons to match the FIBSEM names
  #' @return Returns a renamed vector of neurons
  #' @param nameids: FIBSEM names as returned from neurprintr
  #' @param newNameids: A vector of renamed neuron names
  
  newNameids <- nameids
  newNameids <- gsub("\\([^\\)]+\\)","",newNameids)
  newNameids <- gsub("_L_","_",newNameids)
  newNameids <- gsub("_R_","_",newNameids)
  newNameids <- gsub("_a_","1_",newNameids)
  newNameids <- gsub("_b_","2_",newNameids)
  newNameids <- gsub("_","-",newNameids)
  
  glomNames <- sapply(newNameids, function(n){strsplit(n,'-')[[1]][2]}) %>% as.character()
  for (g in 1:length(glomNames)){
    if (nchar(glomNames[g]) == 2){
      if (substr(glomNames[g],1,1) == "R"){
        glomNames[g] <- paste0(substr(glomNames[g],2,2),"L")
      } else {
        glomNames[g] <- paste0(substr(glomNames[g],2,2),"R")
      }
    }
    if (nchar(glomNames[g]) == 4){
      if (substr(glomNames[g],1,1) == "R"){
        glomNames[g] <- paste0(substr(glomNames[g],2,2),"L",
                               substr(glomNames[g],4,4),"R")
      } else {
        glomNames[g] <- paste0(substr(glomNames[g],2,2),"R",
                               substr(glomNames[g],4,4),"L")
      }
    }
    if (nchar(glomNames[g]) == 6){
      if (substr(glomNames[g],1,1) == "R"){
        glomNames[g] <- paste0(substr(glomNames[g],2,2),"L",
                               substr(glomNames[g],4,4),"L",
                               substr(glomNames[g],6,6),"R")
      } else {
        glomNames[g] <- paste0(substr(glomNames[g],2,2),"R",
                               substr(glomNames[g],4,4),"R",
                               substr(glomNames[g],6,6),"L")
      }
    }
    newNameids[g] <- paste(strsplit(newNameids[g],'-')[[1]][1],glomNames[g],strsplit(newNameids[g],'-')[[1]][3],sep='-');
  }
  
  names <- sapply(newNameids, function(n) paste(strsplit(n,'-')[[1]][1],strsplit(n,'-')[[1]][2],sep='-')) %>%
    as.character() %>% unique() %>% sort()
  for (n in 1:length(names)){
    nronLocs <- which(grepl(names[n], newNameids))
    if (length(nronLocs) > 1){
      newNameids[nronLocs] <-paste0(names[n],letters[1:length(nronLocs)])
    } else {
      newNameids[nronLocs] <-names[n]
    }
  }
  
  return(newNameids)
}

### Sort the FIBSEM data to match the FAFB data
FAFBSort <- function(FIBSEMDat,region){
  #' Sort the names of the neurons in a FIBSEM connectivity matrix to match the FAFB data
  #' @return Returns a data frame of the FIBSEM data with the names sorted to match the FAFB data
  #' @param region: The chosen ROI
  #' @param FIBSEMDat: A FIBSEM dataframe that holds the connectivity data
  
  if (region == "PB"){
    nameOrder <- c("Delta7-1R9R8L","Delta7-2R7L","Delta7-5R4L","Delta7-6R3L","Delta7-7R2L",
                   "EPG-5L","EPG-5R","EPG-6L","EPG-6R",
                   "PEG-5L","PEG-5R","PEG-6L","PEG-6R",
                   "PEN1-5L","PEN1-5R","PEN1-6L","PEN1-6R","PEN1-7L","PEN1-7R",
                   "PEN2-5L","PEN2-5R","PEN2-6L","PEN2-6R","PEN2-7L", "PEN2-7R")
  }
  if (region == "PB-subset"){
    nameOrder <- c("EPG-5R","EPG-5L","EPG-6R","EPG-6L",
                   "Delta7-5R4L","Delta7-6R3L","Delta7-7R2L","Delta7-2R7L","Delta7-1R9R8L",
                   "PEN1-5R","PEN1-5L","PEN1-6R","PEN1-6L","PEN1-7R","PEN1-7L",
                   "PEN2-5R","PEN2-5L","PEN2-6R","PEN2-6L","PEN2-7R", "PEN2-7L",
                   "PEG-5R","PEG-5L","PEG-6R","PEG-6L")
  }
  if (region == "EB"){
    nameOrder <- c("EPG-5R","EPG-5L","EPG-6R","EPG-6L",
                   "PEN1-5R","PEN1-5L","PEN1-6R","PEN1-6L","PEN1-7R","PEN1-7L",
                   "PEN2-5R","PEN2-5L","PEN2-6R","PEN2-6L","PEN2-7R", "PEN2-7L",
                   "PEG-5R","PEG-5L","PEG-6R","PEG-6L",
                   "R4d")
  }
  if (region == "NO"){
    nameOrder <- c("PEN1-5R","PEN1-6R","PEN1-7R","PEN1-5L","PEN1-6L","PEN1-7L",
                   "PEN2-5R","PEN2-6R","PEN2-7R","PEN2-5L","PEN2-6L","PEN2-7L")
  }
  
  names <- c(FIBSEMDat$nameid.from %>% as.character(),FIBSEMDat$nameid.to %>% as.character()) %>% unique()
  nameLevels <- c()
  for (n in 1:length(nameOrder)){
    nameLevels <- c(nameLevels,
                    sort(names[which(grepl(nameOrder[n],names))]))
  }
  
  FIBSEMDat$nameid.from <- factor(FIBSEMDat$nameid.from, levels = nameLevels)
  FIBSEMDat$nameid.to <- factor(FIBSEMDat$nameid.to, levels = nameLevels)
  
  return(FIBSEMDat)
}

### Downselect neurons of interest
FIBSEMDownSel <- function(FIBSEMDat, names, nums){
  #' Remove randomly selected neurons of a given name from a dataframe
  #' @return Returns a downselected data frame of the FIBSEM data
  #' @param FIBSEMDat: A FIBSEM dataframe that holds the connectivity data
  #' @param names: A list of names of neurons to downselect
  #' @param nums: The number of neurons to keep from each named neuron group
  
  nronNames <- c(FIBSEMDat$nameid.to %>% as.character(),FIBSEMDat$nameid.from%>% as.character()) %>% unique()
  
  for (n in 1:length(names)){
    nronsNow <- nronNames[which(grepl(names[n],nronNames))]
    nronsToExclude <- sample(nronsNow,length(nronsNow)-nums[n])
    FIBSEMDat <- FIBSEMDat[which(!(FIBSEMDat$nameid.from %in% nronsToExclude) & !(FIBSEMDat$nameid.to %in% nronsToExclude)),]
  }
  return(FIBSEMDat)
}



### Function to compare the FAFB and FIBSEM data
datComp <- function(FAFBDat,FIBSEMDat){
  #' For a FAFB and a FIBSEM dataframe (that feature comprable neurons), find 
  #' 1) all instances where one dataset has connections that are are absent in the other
  #' 2) Where both datasets show a connection between neurons, compare the synapse counts by calculating: 
  #' (# of FAFB synapses - # of FIBSEM synapses) / (# of FAFB synapses + # of FIBSEM synapses).
  #' @return Returns a dataframe holding the metrics for comparison
  #' @param FIBSEMDat: A FIBSEM dataframe that holds the connectivity data
  #' @param names: A list of names of neurons to downselect
  #' @param nums: The number of neurons to keep from each named neuron group
  
  # Convert the dataframes to matrices
  FAFBMat <- FAFBDat %>% dcast(nameid.from~nameid.to)
  rownames(FAFBMat) <- FAFBMat$nameid.from
  FAFBMat <- FAFBMat[,2:ncol(FAFBMat)] %>% data.matrix
  
  FIBSEMMat <- FIBSEMDat %>% dcast(nameid.from~nameid.to)
  rownames(FIBSEMMat) <- FIBSEMMat$nameid.from
  FIBSEMMat <- FIBSEMMat[,2:ncol(FIBSEMMat)] %>% data.matrix
  
  # Find the sum and difference of the two matrices
  matDiff <- FAFBMat - FIBSEMMat
  matSum <- FAFBMat + FIBSEMMat
  
  # Calculate the two metrics
  perMat <- matrix(0L,nrow=nrow(matDiff),ncol=ncol(matDiff))
  diffMat <- matrix(0L,nrow=nrow(matDiff),ncol=ncol(matDiff))
  equalMat <- matrix(0L,nrow=nrow(matDiff),ncol=ncol(matDiff))
  for (r in 1:nrow(perMat)){
    for (c in 1:ncol(perMat)){
      if (FAFBMat[r,c] == 0){
        diffMat[r,c] <- -FIBSEMMat[r,c]
        next
      }
      if (FIBSEMMat[r,c] == 0){
        diffMat[r,c] <- FAFBMat[r,c]
        next
      }
      if (matSum[r,c] > 0)
        perMat[r,c] <- matDiff[r,c]/matSum[r,c]
    }
    if (FAFBMat[r,c] == FIBSEMMat[r,c]){
      equalMat[r,c] <- FAFBMat[r,c]
    }
  }
  
  # Name the rows and columns appropriately
  rownames(perMat) <- rownames(FAFBMat)
  rownames(diffMat) <- rownames(FAFBMat)
  rownames(equalMat) <- rownames(FAFBMat)
  colnames(perMat) <- colnames(FAFBMat)
  colnames(diffMat) <- colnames(FAFBMat)
  colnames(equalMat) <- colnames(FAFBMat)
  
  # Convert everything to a dataframe
  perMat <- melt(perMat,value.name="per",varnames=c("nameid.from","nameid.to"))
  diffMat <- melt(diffMat,value.name="diff",varnames=c("nameid.from","nameid.to"))
  equalMat <- melt(equalMat,value.name="equal",varnames=c("nameid.from","nameid.to"))
  metrics <- perMat
  metrics$diff <- diffMat$diff
  metrics$equal <- equalMat$equal
  
  return(metrics)
}

### Plot connectivity matrices
conMatPlot <- function(conMat,maxCts){
  g <- ggplot(conMat) + 
    theme_classic() + theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_gradient2(low="white", mid="#ff694a", high="#6b000d", midpoint =0.5*maxCts, limits=c(0,maxCts)) +
    geom_tile(aes(nameid.to,nameid.from,fill=ROIweight))+coord_fixed(ratio=1)
  return(g)
}




### Function to find neuron names and the # of neurons with that name
namesAndNums <- function(df){
  #' Find the anatomical names of different individual neurons
  #' then find the number of anatomically similar neurons of that name
  #' @return Returns a dataframe with the names and the numbers of neurons with that name
  #' @param df: A FIBSEM dataframe that holds the connectivity data
  
  names <- c(df$nameid.from %>% as.character(),df$nameid.to %>% as.character()) %>% unique()
  uniqueNames <- c(names[which(substr(names,nchar(names),nchar(names))=="L")],
                   names[which(substr(names,nchar(names),nchar(names))=="R")])
  otherNames <- names[which(substr(names,nchar(names),nchar(names))!="L" &
                              substr(names,nchar(names),nchar(names)) !="R")]
  otherNames <- substr(otherNames,1,nchar(otherNames)-1) %>% unique()
  allNames <- c(uniqueNames,otherNames) %>% sort()
  nameNums <- sapply(allNames, function(n){length(which(grepl(n,names)))}) %>% as.numeric()
  nameMatch <- data.frame(name=allNames, num=nameNums)
  
  return(nameMatch)
}

### Function to permute the possible FIBSEM to FAFB matchups and calculate comparison metrics
permuteComp <- function(FAFBDat, FIBSEMDat,region,maxIter){
  #' Compare a FAFB connectivity matrix to a corresponding FIBSEM matrix,
  #' randomly permuting the FIBSEM matrix each time if
  #' the total number of permutations is greater than maxIter.
  #' Otherwise, all parameterizations are considered
  #' @return Returns a dataframe with the metrics for each iteration
  #' @param FAFBDat: A FAFM dataframe that holds the connectivity data
  #' @param FIBSEMDat: A FIBSEM dataframe that holds the connectivity data
  #' @param region: The ROI to consider
  #' @param maxIter: The maximum number of permutations to run

  # Find the neuron names and the number of each name from the FAFB data
  FAFBNames <- c(FAFBDat$nameid.from,FAFBDat$nameid.to) %>% unique()
  uniqueNames <- c(FAFBNames[which(substr(FAFBNames,nchar(FAFBNames),nchar(FAFBNames))=="L")],
                   FAFBNames[which(substr(FAFBNames,nchar(FAFBNames),nchar(FAFBNames))=="R")])
  otherNames <- FAFBNames[which(substr(FAFBNames,nchar(FAFBNames),nchar(FAFBNames))!="L" &
                                  substr(FAFBNames,nchar(FAFBNames),nchar(FAFBNames)) !="R")]
  otherNames <- substr(otherNames,1,nchar(otherNames)-1) %>% unique()
  allNames <- c(uniqueNames,otherNames)
  nameNums <- sapply(allNames, function(n){length(which(grepl(n,FAFBNames)))}) %>% as.numeric()
  nameMatch <- data.frame(name=allNames, num=nameNums)
  
  # Find the neuron names and the number of each name for each dataset
  nameMatch <- namesAndNums(FAFBDat)
  colnames(nameMatch) <- c("name","num_FAFB")
  nameMatch$num_FIBSEM <- namesAndNums(FIBSEMDat)$num
  nameMatch <- nameMatch %>% mutate(diff = num_FIBSEM - num_FAFB)
  nameMatch <- nameMatch[which(nameMatch$num_FIBSEM != 1),]
  nameMatch <- nameMatch %>% mutate(possibleCombs = factorial(num_FIBSEM)/factorial(num_FIBSEM-num_FAFB))
  totalCombs <- prod(nameMatch$possibleCombs)
  
  allMetrics <- data.frame(nameid.from <- c(), nameid.to <- c(),
                           per <- c(), diff <- c(), equal <- c(),
                           iteration <- c())
  
  FAFBDat <- FAFBSort(FAFBDat,region)
  if (totalCombs > maxIter){
    for (it in 1:maxIter){
      # Downselect the FIBSEM data
      FIBSEMDat_Cut <- FIBSEMDownSel(FIBSEMDat, nameMatch$name, nameMatch$num_FAFB)
      
      # Convert data to a matrix
      FIBSEMMat <- FIBSEMDat_Cut %>% dcast(nameid.from~nameid.to)
      rownames(FIBSEMMat) <- FIBSEMMat$nameid.from
      FIBSEMMat <- FIBSEMMat[,2:ncol(FIBSEMMat)] %>% data.matrix
      
      # Shuffle the rows of the matrix for the appropraite rows/columns
      FIBSEMMatShuf <- FIBSEMMat
      for (n in 1:nrow(nameMatch)){
        rows <- which(grepl(nameMatch$name[n],rownames(FIBSEMMat)))
        if (length(rows) > 1){
          FIBSEMMatShuf[sample(rows),]=FIBSEMMatShuf[rows,] %>% as.numeric()
        }
        cols <- which(grepl(nameMatch$name[n],colnames(FIBSEMMat)))
        if (length(cols) > 1){
          FIBSEMMatShuf[,sample(cols)]=FIBSEMMatShuf[,cols] %>% as.numeric()
        }
      }
      
      # Convert the data back to a dataframe
      FIBSEMDat_Cut_Shuf <- melt(FIBSEMMatShuf,value.name="ROIweight",varnames=c("nameid.from","nameid.to"))
      FIBSEMDat_Cut_Shuf <- FAFBSort(FIBSEMDat_Cut_Shuf,region)
      
      # Calculate the metrics
      metrics <- datComp(FAFBDat,FIBSEMDat_Cut_Shuf)
      metrics$iteration <- it
      
      allMetrics <- rbind(allMetrics,metrics)
    }
  } else {
    # Get the names of all of the FIBSEM neurons
    FIBSEMNames <- c(FIBSEMDat$nameid.from %>% as.character(),FIBSEMDat$nameid.to %>% as.character()) %>% unique()
    
    # Find the possible permuations for each name
    members <- list()
    for (n in 1:length(nameMatch$name)){
      posCombs <- combn(FIBSEMNames[which(grepl(nameMatch$name[n],FIBSEMNames))],nameMatch$num_FAFB)
      if (is.null(ncol(posCombs))){
        members[[n]] <- permn(posCombs)
      } else {
        permList <- list()
        for (c in 1:ncol(posCombs)){
          permList <- c(permList,permn(posCombs[,c]))
        }
        members[[n]] <- permList
      }
    }
    nameMatch$members <- members
    
    # Convert data to a matrix
    FIBSEMMat <- FIBSEMDat %>% dcast(nameid.from~nameid.to)
    rownames(FIBSEMMat) <- FIBSEMMat$nameid.from
    FIBSEMMat <- FIBSEMMat[,2:ncol(FIBSEMMat)] %>% data.matrix
    rowNamesOrig <- rownames(FIBSEMMat)
    colNamesOrig <- colnames(FIBSEMMat)
    
    # Step through the permutations
    for (it in 1:totalCombs){
      
      # Make a new FIBSEM data matrix
      FIBSEMMatPerm <- FIBSEMMat
      
      # Initialize new row and column names
      rowNamesPerm <- rowNamesOrig
      colNamesPerm <- colNamesOrig
      
      # Modify the row and column names according to the current permutation 
      fact <- 1
      order <- c()
      step <- 1
      while(step <= nrow(nameMatch)){
        orderNow <- nameMatch$members[[step]][1+((ceiling(it/fact)-1) %% nameMatch$num_FAFB[step])][[1]]
        rowIds <- which(grepl(nameMatch$name[step],rowNamesPerm))
        colIds <- which(grepl(nameMatch$name[step],colNamesPerm))
        if (length(rowIds)>0) {
          rowNamesPerm[rowIds] <- orderNow
        }
        if (length(colIds)>0) {
          colNamesPerm[colIds] <- orderNow
        }
        fact <- fact * nameMatch$num_FAFB[step]
        step <- step+1
      }
      
      # Assign the new row and column names to the matrix
      rownames(FIBSEMMatPerm) <- rowNamesPerm
      colnames(FIBSEMMatPerm) <- colNamesPerm
      
      # Convert the data back to a dataframe
      FIBSEMDat_Perm <- melt(FIBSEMMatPerm,value.name="ROIweight",varnames=c("nameid.from","nameid.to"))
      FIBSEMDat_Perm <- FAFBSort(FIBSEMDat_Perm,region)
      
      # Calculate the metrics
      metrics <- datComp(FAFBDat,FIBSEMDat_Perm)
      metrics$iteration <- it
      
      allMetrics <- rbind(allMetrics,metrics)
    }
  }
  return(allMetrics)
}