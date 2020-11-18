#' Read GENEPOP format file
#'
#' This function allows you to import a GENEPOP format file into R. Population names can be specified in the argument. See http://genepop.curtin.edu.au/help_input.html for details about GENEPOP format.
#' @param x GENEPOP file or path to the file. The filename extension (e.g., .txt) should be included.
#' @param pop.names A character string vector for population names. The order of the name should be the same with the order (top to down) in your GENEPOP file.
#' @param haploid A logical variable (TRUE or FALSE) to specify whether your dataset is haploid data. Default is FALSE.
#' @param pos A parameter for program development use; users can ignore it.
#' @return This function returns a list comprising three elements. 1. YOU_NAME_IT$DataMatrix: A matrix of genetic data with a population name label ($popNameVector) in the last column. 2. YOU_NAME_IT$SampleID: A vector of sample ID. 3. YOU_NAME_IT$LocusName: A vector of locus name.
#' @examples # infile <- read.Genepop("Your_Genepop_File.txt", pop.names=c("pop_A", "pop_B", "pop_C"))
#' @references Rousset, F. 2008. Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. Mol. Ecol. Resources 8: 103-106
#' @import stringr
#' @importFrom reshape2 melt
#' @importFrom utils packageVersion
#' @export
#'
read.Genepop <- function(x, pop.names=NULL, haploid = FALSE, pos=1){
  dataType <- NULL
  df <- readLines(x)
  df <- df[-1] #remove first line of file description
  popIndex <- grep("pop", df, ignore.case=T)#get "pop" tag index
  noPops <- length(popIndex)
  
  #check out pop name and verify number of pops
  if(length(pop.names)==0){
    pop.names <- paste0("pop.", seq_along(1:noPops))  
  }else if(length(pop.names)>0){
    if(!length(pop.names) == noPops){
      stop("Lengths of 'pop.names' and 'pop' in file not match...")
    }
  }
  
  #Extract locus names and save them in locusNames 
  if(popIndex[1]==2){ #for one-row locus format
    locusNames <- str_trim(df[1],side="both")
    locusNames <- strsplit(locusNames, ",")[[1]]
  }else if(popIndex[1] > 2){ #for multi-row locus format
    index <- popIndex[1] - 1
    locusNames <- df[1:index]
  }
  locusNames <- str_trim(locusNames, side="both")
  noLocus <- length(locusNames)
  
  #Get index for individuals
  for(i in 1:noPops){
    start <- popIndex[i] + 1
    if(i < noPops){
      end <- popIndex[i+1] - 1
      
    }else if(i == noPops){
      end <- length(df)
    }
    assign(paste0("pop_", i, "_index"), start:end, envir=as.environment(pos))
  }
  pop_all <- lapply(paste0("pop_",seq_along(1:noPops),"_index"), FUN=get)
  
  #save individal index in one vector 
  ind_all_index <- unlist(pop_all)
  #count total number of individuals
  noInds <- length(ind_all_index)
  
  #extract individual genetic data
  ind_df <- df[ind_all_index]
  #separate individual ID and genetic data
  ind_df <- strsplit(ind_df, split = ",")
  #get individual ID
  id_vector <- unlist(lapply(ind_df,`[[`,1))
  #clear extra space on both side of id, if exist
  id_vector <- str_trim(id_vector, side="both")
  
  #get genotype data
  cat("\n  Converting data format...\n")
  geno_list <- unlist(lapply(ind_df,`[[`,2))
  geno_list <- str_trim(geno_list, side="both")
  #separate each locus by spaces or tabs 
  geno_list <- strsplit(geno_list, split="[ \t]+", perl=T)
  #convert nested list to matrix; high computing step
  geno_mx <- matrix(unlist(geno_list), nrow=noInds, byrow=T)
  
  #check number of digits in one locus (if fewer or equal to 3, set haploid=T)
  locusCharSize <- nchar(geno_mx[1,1])
  if(locusCharSize <= 3){
    haploid = T
  }
  
  #apply one-hot encoding for genetic data; high computing step
  cat("\n  Encoding genetic data...\n")
  if(haploid){
    dataType <- "haploid"
    onehot_list <- apply(geno_mx,2 ,genepop_onehot, ploidy=1, noChar=locusCharSize)
  }else{
    dataType <- "diploid"
    onehot_list <- apply(geno_mx,2 ,genepop_onehot, ploidy=2, noChar=locusCharSize)
  }
  #check if entire list is NA
  if(all(is.na(onehot_list))){
    stop("Entire NA data due to identical genotype across samples.")
  }
  #check and remove locus that is NA
  LocusNA_idx <- which(is.na(onehot_list))
  #remove NA locus if exists
  if(length(LocusNA_idx)>0){
    onehot_list <- onehot_list[-LocusNA_idx]
    #get locus name
    locusNames <- locusNames[-LocusNA_idx]
  }
  #change dataframe's colnames 
  if(length(onehot_list) == length(locusNames)){
    for(i in 1:length(locusNames)){
      names(onehot_list[[i]]) <- paste0(locusNames[i],"_",names(onehot_list[[i]]))
    }
  }else{
    stop("Oops, lengths of onehot_list and locusNames differ.")
  }
  #concatenate dataframe in onehot_list
  genoMatrix <- do.call(cbind, onehot_list)
  #count number of columns (alleles) in genetic data matrix
  noLociVar <- ncol(genoMatrix)
  
  #Create pop name vector and concatenate to the genoMatrix
  popNames_vector <- NULL
  for (i in 1:noPops){
    popsize <- length(pop_all[[i]])
    popNameVector <- rep(pop.names[i], popsize)
    popNames_vector <- c(popNames_vector, popNameVector)
  }
  genoMatrix <- cbind(genoMatrix, popNames_vector)
  
  #Print some message to console
  cat(paste0("\n  ################ assignPOP v",packageVersion("assignPOP")," ################\n"))
  cat("\n  A GENEPOP format file was successfully imported!\n")
  cat(paste0("\n  Imported Data Info: ",noInds," obs. by ",noLocus," loci (",dataType,")"))
  cat(paste0("\n  Number of pop: ",noPops))
  for(i in 1:noPops){
    popSize <- length(get(paste0("pop_",i,"_index")))
    cat(paste0("\n  Number of inds (",pop.names[i],"): ",popSize ) )
  }
  cat(paste0("\n  DataMatrix: ",nrow(genoMatrix)," rows by ",ncol(genoMatrix), " columns, with ",noLociVar," allele variables"))
  
  cat("\n")
  cat("\n  Data output in a list comprising the following three elements:")
  cat("\n  YOUR_LIST_NAME$DataMatrix")
  cat("\n  YOUR_LIST_NAME$SampleID")
  cat("\n  YOUR_LIST_NAME$LocusName")
  cat("\n\n")
  
  #Remove variables from GlobalEnv.
  rm(list = ls(pattern="^pop_.*_index$", envir = .GlobalEnv), envir = .GlobalEnv)
  
  finalList <- list(genoMatrix, id_vector, locusNames)
  names(finalList) <- c("DataMatrix", "SampleID" , "LocusName")
  return(finalList)
  
}

########################################
# Genepop genetic data one-hot encoding
########################################
genepop_onehot <- function(oneLoc, ploidy=NULL, noChar=NULL){
  #x is character string vector of a locus
  #sample test
  #oneLoc <- geno_mx[,3] #multi-alleles 
  #oneLoc <- geno_mx[,4] #all NA
  #oneLoc <- geno_mx[,5] #single allele
  
  #check if only one allele
  oneAllele <- FALSE
  if(length(unique(oneLoc))==1){
    oneAllele <- TRUE
  }
  
  #if dataset is haploid
  if(ploidy==1){
    #if only one allele
    if(oneAllele){
      #if all NA 
      if(any(c("0","00","000","0000","000000") %in% oneLoc)){
        onehotDF <- NA
      }else{
        onehotDF <- as.data.frame(rep(1, length(oneLoc)))
        names(onehotDF) <- oneLoc[1]
      }
      #multi-alleles    
    }else{
      #convert one locus vector to dataframe
      oneLocDF <- data.frame(oneLoc, stringsAsFactors = T)
      #get one-hot encoding dataframe
      onehotDF <- as.data.frame(model.matrix(~0+oneLocDF[,1]))
      names(onehotDF) <- levels(oneLocDF$oneLoc)
      #remove missing data
      if(any(c("0","00","000","0000","000000") %in% names(onehotDF))){
        onehotDF <- onehotDF[ , -which(names(onehotDF) %in% c("0","00","000","0000","000000")), drop=FALSE]
      }
    }
    #if dataset is diploid
  }else if(ploidy==2){
    #if only one allele
    if(oneAllele){
      #if all NA
      if(any(c("0","00","000","0000","000000") %in% oneLoc)){
        onehotDF <- NA
      }else{
        onehotDF <- as.data.frame(rep(1, length(oneLoc)))
        names(onehotDF) <- substr(oneLoc[1],1,as.integer(noChar/2))
      }
      #if multi-alleles  
    }else{
      #separate alleles
      alleles <- strsplit(oneLoc, split=paste0("(?<=.{",noChar/2,"})"), perl=T)
      alleles <- unlist(alleles)
      #convert one locus to dataframe
      alleles_DF <- data.frame(alleles, stringsAsFactors = T)
      onehotMX <- model.matrix(~0+alleles_DF[,1])
      onehotDF <- as.data.frame((onehotMX[c(T,F),] + onehotMX[c(F,T),])/2)
      rownames(onehotDF) <- NULL
      names(onehotDF) <- levels(alleles_DF$alleles)
      #remove missing data
      if(any(c("0","00","000","0000","000000") %in% names(onehotDF))){
        onehotDF <- onehotDF[ , -which(names(onehotDF) %in% c("0","00","000","0000","000000")), drop=FALSE]
        #drop=FALSE allows to keep it as data frame. Otherwise, only one column will become a vector
      }
    }
  }
  return(onehotDF)
}
