#' Read GENEPOP format input file
#'
#' This function allows you to import a GENEPOP format file into R. Population names can be specified in the argument. See http://genepop.curtin.edu.au/help_input.html for details about GENEPOP format.
#' @param x The GENEPOP file name or path to the file. The filename extension should be included.
#' @param pop.names A character string vector for population names. The order of the name should be the same with the order (top to down) in your GENEPOP file.
#' @param pos defaults to 1 which equals an assingment to global environment
#' @return This function returns a list that comprises three items. [[1]] A matrix of genetic data with a population name label ($popNameVector) in the last column. [[2]] A vector of sample ID. [[3]] A vector of locus name.
#' @examples # infile <- read.genpop("Your_Genepop_File.txt", pop.names=c("pop_A", "pop_B", "pop_C"))
#' @references Rousset, F., 2008. Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. Mol. Ecol. Resources 8: 103-106.
#' @import stringr
#' @importFrom reshape2 melt
#' @importFrom utils setTxtProgressBar txtProgressBar packageVersion
#' @examples 
#' genin <- read.genpop(system.file("extdata/TinyGenepop.txt", package="assignPOP"))
#' #Change file 'TinyGenepop' to 'simGenepop' to get the example used in the tutorial.
#' @export
#'
read.genpop <- function(x, pop.names = NULL, pos=1){
  df <- readLines(x)
  df <- df[-1]#skip the first line file info
  popIndex <- grep("pop",df,ignore.case=TRUE)
  noPops <- length(popIndex)
  if(length(pop.names)==0){ #check if pop.names is specified, if not, assign pop names
    pop.names <- paste0("pop.",seq_along(1:noPops))
  } else if(length(pop.names)>0){ #check if pop.names and number of pop match
    if(!length(pop.names)==noPops){
      cat("\nError: Pop.names and number of pop in data not match...")
      break } }

  #Extract locus name, save locus names in vector "locusNames"
  if (popIndex[1] == 2) {
    locusNames <- str_trim(df[1],side="both")
    locusNames <- strsplit(locusNames,",")[[1]]
  } else if (popIndex[1] > 2) {
    index <- popIndex[1] - 1
    locusNames <- df[1:index]
  }
  locusNames <- str_trim(locusNames,side="both")
  noLocus <- length(locusNames)

  #Get index for individuals and save in a nested list "pop_all"
  for (i in 1:noPops){
    if (i < noPops){
      start <- popIndex[i] + 1; end <- popIndex[i+1] - 1
      assign(paste0("pop_",i,"_index"), start : end, envir=as.environment(pos) )
    } else if (i == noPops){
      start <- popIndex[i] + 1; end <- length(df)
      assign(paste0("pop_",i,"_index"), start : end, envir=as.environment(pos) )
    }
  }
  pop_all <- lapply(paste0("pop_", seq_along(1:noPops),"_index"), FUN = get )

  #save individual index in one vector "ind_all_index"
  ind_all_index <- NULL
  for (i in 1:noPops){
    ind_all_index <- c(ind_all_index, pop_all[[i]])
  }
  noInds <- length(ind_all_index)

  #extract individual data
  ind_df <- df[ind_all_index]
  #split individual ID and genotype data
  id_vector <- NULL
  geno_list <- list()
  for(i in 1:noInds){
    id_n_genotype <- strsplit(ind_df[i],",")[[1]]
    id <- str_trim(id_n_genotype[1], side="both")
    id_vector <- c(id_vector,id)
    geno <- str_trim(id_n_genotype[2], side="both")
    geno <- gsub("\\s+"," ",geno)#clean extra space or change tabs to one single space between loci
    geno <- strsplit(geno," ")[[1]]#make each locus an element
    geno <- list(geno)
    geno_list <- c(geno_list, geno)
  }

  #Get alleles of each locus across individuals
  #create an empty data frame
  genoMatrix <- data.frame(matrix(ncol=0,nrow=noInds))
  #Setup progress bar
  pb <- txtProgressBar(min = 0, max = noLocus, style = 3)
  #Save missing locus index in a vector
  missLocusIndex <- NULL
  for(m in 1:noLocus){
    oneLocus_vector <- NULL
    setTxtProgressBar(pb, m)
    for(n in 1:noInds){
      eachlocus <- geno_list[[n]][m]
      noChar <- nchar(eachlocus)
      diploid <- substring(eachlocus, c(1,(noChar/2)+1), c(noChar/2,noChar))
      oneLocus_vector <- union(oneLocus_vector, diploid)
    }
    oneLocus_vector <- sort(oneLocus_vector)
    #remove "00" or "000" data
    oneLocus_vector <- oneLocus_vector[! oneLocus_vector %in% "00"]
    oneLocus_vector <- oneLocus_vector[! oneLocus_vector %in% "000"]
    #create multi variables for locus
    noVarOfLocus <- length(oneLocus_vector)
    locusDf <- data.frame(matrix(ncol=noVarOfLocus, nrow=0))

    #check if a locus is a missing data across individuals, if so, save the locus index
    if(noVarOfLocus==0){
      missLocusIndex <- c(missLocusIndex,m)
    }
    noMissLocusIndex <- length(missLocusIndex)
    #check individual genotype and convert genotype data to binary-like
    for(p in 1:noInds){
      eachlocus <- geno_list[[p]][m]
      noChar <- nchar(eachlocus)#count number of character (usually either 4 or 6 digits)
      diploid <- substring(eachlocus, c(1,(noChar/2)+1), c(noChar/2,noChar))
      binaryVector <- oneLocus_vector %in% diploid
      #count number of TRUE
      noTrue <- table(binaryVector)["TRUE"]
      noTrue[is.na(noTrue)] <- 0
      if(noTrue==1){ #if there is only one TRUE, meaning a homozygote
        binaryVector <- binaryVector*1
      } else {
        binaryVector <- binaryVector*0.5
      }
      locusDf <- rbind(locusDf,binaryVector)
    }
    if(!noVarOfLocus==0){
      names(locusDf) <- paste0(locusNames[m],"_",seq_along(1:noVarOfLocus))
      genoMatrix <- cbind(genoMatrix,locusDf)
    }
  }
  #count number of columns (alleles) in genetic data matrix
  noLociVar <- ncol(genoMatrix)
  #Remove locus name if it's missing data across individuals
  if(!noMissLocusIndex==0){
    locusNames <- locusNames[-missLocusIndex]
  }
  #close progrss bar
  close(pb)
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
  cat(paste0("\n  DataInfo: ",noInds," obs. by ",noLocus," loci (with ",noLociVar," var.)"))
  cat(paste0("\n  DataMatrix: ",nrow(genoMatrix)," rows by ",ncol(genoMatrix), " columns"))
  cat(paste0("\n  Number of pop: ",noPops))
  for(i in 1:noPops){
    popSize <- length(get(paste0("pop_",i,"_index")))
    cat(paste0("\n  Number of inds (",pop.names[i],"): ",popSize ) )
  }
  cat("\n")
  cat("\n  Data output in a list comprising the following three elements:")
  cat("\n  DataMatrix @ YOUR_LIST_NAME[[1]]")
  cat("\n  IndividualNames @ YOUR_LIST_NAME[[2]]")
  cat("\n  LocusNames @ YOUR_LIST_NAME[[3]]")
  cat("\n\n")

  #Remove variables from GlobalEnv.
  rm(list = ls(pattern="^pop_.*_index$", envir = .GlobalEnv), envir = .GlobalEnv)
  
  return(list(genoMatrix, id_vector, locusNames))
  
}
