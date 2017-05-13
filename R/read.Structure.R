#' Read Structure format file
#' 
#' This function allows you to import a STRUCTURE format file into R. The first row should be locus name (either with or withour column names for sample ID and population label); the first column should be sample ID; the second column should be population label; the rest are genotype. Use "-9" for missing alleles. 
#' @param x STRUCTURE file or path to the file. The filename extension (e.g., .txt) should be included.
#' @param haploid A logical variable (TRUE or FALSE) to specify whether your dataset is haploid data. Default is FALSE.
#' @return This function returns a list comprising three elements. 1. YOU_NAME_IT$DataMatrix: A matrix of genetic data with a population name label ($popNameVector) in the last column. 2. YOU_NAME_IT$SampleID: A vector of sample ID. 3. YOU_NAME_IT$LocusName: A vector of locus name.
#' @examples # infile <- read.Structure("Your_Structure_File.txt")
#' @references Pritchard, J.K., Stephens, M. and Donnelly, P., 2000. Inference of population structure using multilocus genotype data. Genetics, 155(2), pp.945-959.
#' @import stringr
#' @importFrom utils setTxtProgressBar txtProgressBar packageVersion
#' @export
#' 
read.Structure <- function(x, haploid = FALSE){
  dataType <- NULL
  df <- readLines(x)
  #Check if there are sample id and pop label in first row
  firstrow <- df[1]
  firstrow <- str_trim(firstrow, side="both")
  firstrow <- gsub("\\s+"," ",firstrow)
  firstrow <- strsplit(firstrow, " ")[[1]]
  secondrow <- df[2]
  secondrow <- str_trim(secondrow, side="both")
  secondrow <- gsub("\\s+"," ",secondrow)
  secondrow <- strsplit(secondrow, " ")[[1]]
  #Check length of first and second row and then extract locus names
  if(length(secondrow) - length(firstrow) == 0){ #there're ID and Pop label 
    locusNames <- firstrow[3:length(firstrow)]
  }else if(length(secondrow) - length(firstrow) == 1){ #there's ID or Pop label
    locusNames <- firstrow[2:length(firstrow)]
  }else if(length(secondrow) - length(firstrow) == 2){ #there's no ID or Pop label
    locusNames <- firstrow
  }else {
    stop("The numbers of locus and locus name may not be the same.")
  }
  noLocus <- length(locusNames)
  #Process for haploid data
  if(haploid){
    dataType <- "haploid"
    #Get number of individuals
    noInds <- length(df) - 1
    #Create sample ID container
    sampleID_vec <- NULL
    #Create pop label container
    popNames_vector <- NULL
    #Create empty list to save genotype
    geno_list <- list()
    #Extrace sample ID, pop label, and genotype
    for(i in 1:noInds){
      eachIndRow <- df[i+1]
      eachIndRow <- str_trim(eachIndRow, side="both")
      eachIndRow <- gsub("\\s+"," ",eachIndRow)
      eachIndRow <- strsplit(eachIndRow," ")[[1]]
      #Get sample id
      sampleID <- eachIndRow[1]
      #Append to sampleID_vec
      sampleID_vec <- c(sampleID_vec, sampleID)
      #Get pop label
      popLabel <- eachIndRow[2]
      #Append to popNames_vector
      popNames_vector <- c(popNames_vector, popLabel)
      #Get genotype
      genotype <- eachIndRow[3:length(eachIndRow)]
      #Append a genotype list into a master list
      genotype <- list(genotype)
      geno_list <- c(geno_list, genotype)
    }
    
    #Get alleles of each locus across individuals
    #create an empty data frame
    genoMatrix <- data.frame(matrix(ncol=0,nrow=noInds))
    #Setup progress bar
    pb <- txtProgressBar(min = 0, max = noLocus, style = 3)
    #Save missing locus index in a vector
    missLocusIndex <- NULL
    #
    for(m in 1:noLocus){
      oneLocus_vector <- NULL
      setTxtProgressBar(pb, m)
      for(n in 1:noInds){
        eachLocus <- geno_list[[n]][m]
        if(eachLocus == "-9"){
          eachLocus = NA
        }
        oneLocus_vector <- c(oneLocus_vector, eachLocus)
      }#for(j in 1:noInds)
      #Check if a locus is missing data across individuals
      if(all(is.na(oneLocus_vector)) | length(unique(oneLocus_vector)) == 1){
        missLocusIndex <- c(missLocusIndex, m)
      }else {
        #Convert one locus dataset to dummy locus variables
        oneLocusDf <- data.frame(oneLocus_vector)
        dummyLocusMatrix <- as.data.frame(model.matrix( ~ oneLocus_vector-1,data = oneLocusDf))
        #Insert rows if there is missing data
        dummyLocusMatrix <- dummyLocusMatrix[match(rownames(oneLocusDf),rownames(dummyLocusMatrix)),]
        rownames(dummyLocusMatrix) <- rownames(oneLocusDf)
        #Edit locus name
        names(dummyLocusMatrix) <- substring(names(dummyLocusMatrix), 16, 1000L)
        names(dummyLocusMatrix) <- sub("\\b", paste0(locusNames[m],"_"), names(dummyLocusMatrix))
        #Append new locus variables to genoMatrix
        genoMatrix <- cbind(genoMatrix, dummyLocusMatrix)
      }
    }#for(i in 1:length(locusNames))
    #Convert NA to 0
    genoMatrix[is.na(genoMatrix)] <- 0
    #close progrss bar
    close(pb)
    #Process for diploid data  
  }else if(!haploid){
    dataType <- "diploid"
    #Get number of individuals
    noInds <- (length(df) - 1)/2
    #Create sample ID container
    sampleID_vec <- NULL
    #Create pop label container
    popNames_vector <- NULL
    #Create empty list to save genotype
    geno_list <- list()
    #Extrace sample ID, pop label, and genotype
    for(i in 1:(noInds*2)){
      eachHapRow <- df[i+1]
      eachHapRow <- str_trim(eachHapRow, side="both")
      eachHapRow <- gsub("\\s+"," ",eachHapRow)
      eachHapRow <- strsplit(eachHapRow," ")[[1]]
      #Get sample ID and pop label for odd row
      if(i %% 2 ==1){
        sampleID <- eachHapRow[1]
        sampleID_vec <- c(sampleID_vec, sampleID)
        popLabel <- eachHapRow[2]
        popNames_vector <- c(popNames_vector, popLabel)
      }
      #Get genotype data for each row (haploid) as a list and append to geno_list
      genotype <- eachHapRow[3: length(eachHapRow)]
      genotype <- list(genotype)
      geno_list <- c(geno_list, genotype)
    }#for(i in 1:noInds)
    
    #create an empty data frame
    genoMatrix <- data.frame(matrix(ncol=0,nrow=noInds))
    #Setup progress bar
    pb <- txtProgressBar(min = 0, max = noLocus, style = 3)
    #Save missing locus index in a vector
    missLocusIndex <- NULL
    #Process genotype data
    for(m in 1:noLocus){
      oneLocus_vector <- NULL
      setTxtProgressBar(pb, m)
      for(n in 1:(noInds*2)){
        eachlocus <- geno_list[[n]][m]
        if(eachlocus == "-9"){
          eachlocus = NA
        }
        oneLocus_vector <- c(oneLocus_vector, eachlocus)
      }#for(n in 1:(noInds*2))
      #Check if a locus is missing data across individuals
      if(all(is.na(oneLocus_vector)) | length(unique(oneLocus_vector)) == 1){
        missLocusIndex <- c(missLocusIndex, m)
      }else{
        #Convert oneLocus_vector to data frame and create dummy variables
        oneLocusDF <- data.frame(oneLocus_vector)
        dummyLocusMatrix <- as.data.frame(model.matrix( ~ oneLocus_vector - 1, data=oneLocusDF))
        #Insert rows for missing data
        dummyLocusMatrix <- dummyLocusMatrix[match(rownames(oneLocusDF),rownames(dummyLocusMatrix)),]
        rownames(dummyLocusMatrix) <- rownames(oneLocusDF)
        #Replace NA with 0
        dummyLocusMatrix[is.na(dummyLocusMatrix)] <- 0
        #Add n and n+1 rows (each row is haploid of diploid individual) to one row in a new data frame; 
        newDummyLocusMatrix <- NULL
        for(i in 1:nrow(dummyLocusMatrix)){
          if( i %% 2 == 1){ #process 1st, 3rd, 5th...row (add to i+1th row)
            eachIndGenotype <- dummyLocusMatrix[i,] + dummyLocusMatrix[i+1,]
            newDummyLocusMatrix <- rbind(newDummyLocusMatrix, eachIndGenotype)
          }
        }
        #Reorder rownames
        rownames(newDummyLocusMatrix) <- seq(1:noInds)
        newDummyLocusMatrix <- newDummyLocusMatrix / 2
        #Edit column(locus) names
        names(newDummyLocusMatrix) <- substring(names(newDummyLocusMatrix), 16, 1000L)
        names(newDummyLocusMatrix) <- sub("\\b", paste0(locusNames[m],"_"), names(newDummyLocusMatrix))
        #Append to genoMatrix
        genoMatrix <- cbind(genoMatrix, newDummyLocusMatrix)
      }
    }#for(m in 1:noLocus)
    close(pb)
  }#else if(!haploid)
  
  #Get population size
  popDF <- data.frame(table(popNames_vector))
  noPops <- nrow(popDF)
  #count number of columns (alleles) in genetic data matrix
  noLociVar <- ncol(genoMatrix)
  #Remove locus name if it's missing data across individuals
  noMissLocusIndex <- length(missLocusIndex)
  if(!noMissLocusIndex==0){
    locusNames <- locusNames[-missLocusIndex]
  }
  #Create pop name vector and concatenate to the genoMatrix
  genoMatrix <- cbind(genoMatrix, popNames_vector)
  
  #Print some message to console
  cat(paste0("\n  ################ assignPOP v",packageVersion("assignPOP")," ################\n"))
  cat("\n  A STRUCTURE format file was successfully imported!\n")
  cat(paste0("\n  Imported Data Info: ",noInds," obs. by ",noLocus," loci (",dataType,")"))
  cat(paste0("\n  Number of pop: ",noPops))
  for(i in 1:noPops){
    cat(paste0("\n  Number of inds (",popDF$popNames_vector[i],"): ",popDF$Freq[i] ) )
  }
  cat(paste0("\n  DataMatrix: ",nrow(genoMatrix)," rows by ",ncol(genoMatrix), " columns, with ",noLociVar," allele variables"))
  
  cat("\n")
  cat("\n  Data output in a list comprising the following three elements:")
  cat("\n  YOUR_LIST_NAME$DataMatrix")
  cat("\n  YOUR_LIST_NAME$SampleID")
  cat("\n  YOUR_LIST_NAME$LocusName")
  cat("\n\n")
  
  #Remove variables from GlobalEnv.
  #rm(list = ls(pattern="^pop_.*_index$", envir = .GlobalEnv), envir = .GlobalEnv)
  
  finalList <- list(genoMatrix, sampleID_vec, locusNames)
  names(finalList) <- c("DataMatrix", "SampleID" , "LocusName")
  return(finalList)
  
} 