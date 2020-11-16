#' Read Structure format file
#' 
#' This function allows you to import a STRUCTURE format file into R. The first row should be locus name (either with or withour column names for sample ID and population label); the first column should be sample ID; the second column should be population label; the rest are genotype. Use "-9" for missing alleles. 
#' @param x STRUCTURE file or path to the file. The filename extension (e.g., .txt) should be included.
#' @param ploidy An integer of 1, 2, 3, or 4, to indicate haploid, diploid, triploid, or tetraploid data. Default is 2 (diploid).
#' @return This function returns a list comprising three elements. 1. YOU_NAME_IT$DataMatrix: A matrix of genetic data with a population name label ($popNameVector) in the last column. 2. YOU_NAME_IT$SampleID: A vector of sample ID. 3. YOU_NAME_IT$LocusName: A vector of locus name.
#' @examples # infile <- read.Structure("Your_Structure_File.txt")
#' @references Pritchard, J.K., Stephens, M. and Donnelly, P., 2000. Inference of population structure using multilocus genotype data. Genetics, 155(2), pp.945-959.
#' @import stringr
#' @importFrom utils setTxtProgressBar txtProgressBar packageVersion
#' @export
#' 
read.Structure <- function(x, ploidy = 2){
  dataType <- NULL
  #check ploidy entry
  if(!ploidy %in% c(1,2,3,4)){
    stop("Ploidy should be an integer of 1, 2, 3, or 4.")
  }
  df <- readLines(x)
  #check if there're sample ID and pop label in col names
  firstrow <- str_trim(df[1], side="both")
  firstrow <- strsplit(firstrow, split="[ \t]+", perl=T)[[1]]
  
  secondrow <- str_trim(df[2], side="both")
  secondrow <- strsplit(secondrow, split="[ \t]+", perl=T)[[1]]
  
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
  
  #remove first row (col names)
  df <- df[-1]
  #separate all elements
  df_split <- strsplit(df, split="[ \t]+", perl=T)
  #convert to a matrix
  master_mx <- matrix(unlist(df_split), nrow=length(df), byrow=T)
  #get sample ID 
  id_vector <- master_mx[,1]
  #get pop label 
  pop_vector <- master_mx[,2]
  #get genetic data
  genoMatrix <- master_mx[,3:ncol(master_mx)]
  
  #one-hot encoding genetic data
  genoMatrix <- apply(genoMatrix, 2, structure_onehot, ploidy=ploidy)
  
  #check if entire NA data
  if(all(is.na(genoMatrix))){
    stop("Entire NA data due to identical genotype across samples.")
  }
  #check and remove locus that is NA
  LocusNA_idx <- which(is.na(genoMatrix))
  #remove NA locus if exists
  if(length(LocusNA_idx)>0){
    genoMatrix <- genoMatrix[-LocusNA_idx]
    #get locus name
    locusNames <- locusNames[-LocusNA_idx]
  }
  
  #change dataframe's colnames 
  if(length(genoMatrix) == length(locusNames)){
    for(i in 1:length(locusNames)){
      names(genoMatrix[[i]]) <- paste0(locusNames[i],"_",names(genoMatrix[[i]]))
    }
  }else{
    stop("Oops, lengths of genoMatrix and locusNames differ.")
  }
  #concatenate dataframe in onehot_list
  genoMatrix <- do.call(cbind, genoMatrix)
  #count number of columns (alleles) in genetic data matrix
  noLociVar <- ncol(genoMatrix)
  
  #get population label
  popNames_vector <- NULL
  noPops <- length(unique(pop_vector))
  
  #for haploid data
  if(ploidy == 1){
    dataType <- "haploid"
    popNames_vector <- pop_vector
    sampleID_vec <- id_vector
  }else if(ploidy == 2){
    dataType <- "diploid"
    popNames_vector <- pop_vector[c(T,F)]
    sampleID_vec <- id_vector[c(T,F)]
  }else if(ploidy == 3){
    dataType <- "triploid"
    popNames_vector <- pop_vector[c(T,F,F)]
    sampleID_vec <- id_vector[c(T,F,F)]
  }else if(ploidy == 4){
    dataType <- "tetraploid"
    popNames_vector <- pop_vector[c(T,F,F,F)]
    sampleID_vec <- id_vector[c(T,F,F,F)]
  }
  #combine genetic data and pop label
  genoMatrix <- cbind(genoMatrix, popNames_vector)
  #get number of individuals
  noInds <- nrow(genoMatrix)
  #get each pop size
  popDF <- data.frame(table(popNames_vector))
  
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
  
  finalList <- list(genoMatrix, sampleID_vec, locusNames)
  names(finalList) <- c("DataMatrix", "SampleID" , "LocusName")
  return(finalList)
}

########################################
#Structure genetic data one-hot encoding
#######################################
structure_onehot <- function(oneLoc, ploidy=NULL){
  #oneLoc <- genoMatrix[,1]
  
  #if a locus across samples has one allele
  if(length(unique(oneLoc))==1){
    #if all NA
    if("-9" %in% oneLoc){
      onehotMX <- NA
    #if only one allele
    }else{
      onehotMX <- as.data.frame(rep(1,length(oneLoc)))
      names(onehotMX) <- oneLoc[1]
    }
    
  }else{
    #below also process for haploid
    oneLocDF <- data.frame(oneLoc, stringsAsFactors=T)
    onehotMX <- as.data.frame(model.matrix(~0+oneLoc, data=oneLocDF))
    #for diploid
    if(ploidy==2){
      onehotMX <- as.data.frame(onehotMX[c(T,F),] + onehotMX[c(F,T),])/2  
    }else if(ploidy==3){
      onehotMX <- (onehotMX[c(T,F,F),] + onehotMX[c(F,T,F),] + onehotMX[c(F,F,T),])/3
    }else if(ploidy==4){
      onehotMX <- (onehotMX[c(T,F,F,F),]+onehotMX[c(F,T,F,F),]+onehotMX[c(F,F,T,F),]+onehotMX[c(F,F,F,T),])/4
    }
    #reset row names 
    rownames(onehotMX) <- NULL
    #rename column names
    names(onehotMX) <- levels(oneLocDF$oneLoc)
    #remove missing alleles (-9)
    if("-9" %in% names(onehotMX)){
      onehotMX <- onehotMX[, -which(names(onehotMX)=="-9"), drop=FALSE]
    }
  }
  return(onehotMX)
}
