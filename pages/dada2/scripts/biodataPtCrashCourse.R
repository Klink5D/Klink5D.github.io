plotFastQC <- function(fwdPathList, revPathList, sampleNamesTag, pathDir) {
  
  sampNr <- length(fwdPathList)
  dir.create(pathDir)
  
  for (fastq in 1:sampNr) {
    sampFwdName <- strsplit(basename(fwdPathList[fastq]), split=sampleNamesTag)[[1]][1]
    sampRevName <- strsplit(basename(revPathList[fastq]), split=sampleNamesTag)[[1]][1]
    
    if (sampFwdName == sampRevName) {
      
      pdf(file = paste0(pathDir, '/', sampFwdName, "_","plotQualityProfile.pdf"), width = 13.5, height = 10)
      fwdFastqPlot <- plotQualityProfile(fwdPathList[fastq])
      revFastqPlot <- plotQualityProfile(revPathList[fastq])
      gridExtra::grid.arrange(fwdFastqPlot, revFastqPlot, ncol = 2)
      dev.off()
      
    }
  
    else {
      
      print(paste0("Forward -", sampFwdName, "- and reverse -", sampRevName, "- fastq sample names do not match!"))
      print("Did you sorted the 'fwdPathList' and 'revPathList' vector lists?")
      
    }
  }
}

#-------------------------------------------------------------------------------------------------------------------------------------------

compareSampleNames <- function(fwdPathList, revPathList, splitSymbol, pickElement) {
  
  # `compareSamplesNames()` function takes as input two path lists
  # it retrieves the `basename()`, i.e., the last directory folder or file
  # it split it by `splitSymbol` and picks up `splitSymbol` one of the elements
  
  fwdNames <- sapply(X = strsplit(x = basename(path = fwdPathList), split = splitSymbol), FUN = `[`, pickElement)
  revNames <- sapply(X = strsplit(x = basename(path = revPathList), split = splitSymbol), FUN = `[`, pickElement)
  
  if (length(fwdNames) == length(revNames)) {
    for (name in 1:length(revNames)) {
      if (fwdNames[name] != revNames[name]) {
        #print("")
        print(paste0(fwdNames[name], " does not match ", revNames[name], "!"))
        print("`fwdPathList` and `revPathList` are not sorted!")
        #print("")
      }
      else {
        #print("")
        print(paste0(fwdPathList[name], " and ", revPathList[name], " are sorted!"))
        #print("")
      } 
    }
  }
  else {
    #print("")
    print("`fwdPathList` and `revPathList` have different length!")
    #print("")
    }
}

#-------------------------------------------------------------------------------------------------------------------------------------------


tax2biom <- function(taxTable) {
 
  lenNrow = nrow(taxTable)
  lenNcol = ncol(taxTable)
  id = c()
  fullTax = c()
  
  
  for (irow in 1:lenNrow) {
    rowID = as.character(rownames(taxTable)[irow])
    id = append(id, rowID)
    tax = ""
    
    for (icol in 1:lenNcol) {
      colName <- colnames(taxTable)[icol]
      colSign <- paste0(tolower(strtrim(colName,1)), "__")
      
      if (tax == "") {
        tax <- paste0(colSign, taxTable[irow,icol])
        #tax <- taxTable[irow,icol]
      }
      else {
        tax <- paste(tax, paste0(colSign, taxTable[irow,icol]), sep = "; ")
        #tax <- paste(tax, taxTable[irow,icol], sep = "; ")
      }
    }
    fullTax <- append(fullTax, tax)
  }
  
  taxTable2Biom <- data.frame(id, fullTax)
  colnames(taxTable2Biom) <- c("ASV_ID", "taxonomy")
  
  return(taxTable2Biom)
}
  

#-------------------------------------------------------------------------------------------------------------------------------------------

countUniqueAndTotalFromDerepObjcList <- function(deRepObj) {
  
  ### count the total no. of unique seqs in the dereplication dada2 object
  
  listOut <- list()
  
  uniq <- sapply(X = deRepObj, FUN = summary)[1,]
  uniqOut <- as.numeric(uniq)
  names(uniqOut) <- names(deRepObj)
  
  total <- sapply(X = deRepObj, FUN = summary)[3,]
  totalOut <- as.numeric(total)
  names(totalOut) <- names(deRepObj)

  listOut[["unique"]] <- uniqOut
  listOut[["total"]] <- totalOut
  
  return(listOut)
  
}



#-------------------------------------------------------------------------------------------------------------------------------------------

countUniqueFromDadaObjcList <- function(dadaObj) {
  
  ### count the total no. of unique seqs in the dada2 object

  uniq <- sapply(X = dadaObj, FUN = summary)[1,]
  uniqOut <- as.numeric(uniq)
  names(uniqOut) <- names(dadaObj)
  
  return(uniqOut)
  
}


#-------------------------------------------------------------------------------------------------------------------------------------------

countMergedSeqsFromDadaObjcList <- function(mergedObj) {
  
  ### count the total no. of merged seqs in the dada2 object
  
  merged <- sapply((sapply(X = mergedObj, `[`, 2)), sum)
  mergedOut <- as.numeric(merged)
  names(mergedOut) <- names(mergedObj)
  
  return(mergedOut)
  
}

  
  
  

