#' Remove low variance alleles (dimensionality reduction)
#'
#' This function helps remove alleles that have low variance in the data set such that it can speed up further analyses for a large data set (e.g., > 10K SNPs).
#' @param x A returned object (a list) from the function read.genpop().
#' @param p A threshold of variance for the alleles to be removed. For example, if p = 0.95 (default setting), an allele occupied more than 95 percents across all the samples will be removed.
#' @return This function return the same object as the function read.genpop() except that the number of columns in the matrix [[1]] is reduced and so is the locus name [[3]].
#' @export
#'
reduce.allele <- function(x, p = 0.95){
  genoMatrix <- x[[1]]
  locusNames <- x[[3]]
  noAlleles <- ncol(genoMatrix)-1
  noInds <- nrow(genoMatrix)
  colToRemove <- NULL
  for(i in 1:noAlleles){ #for loop each allele column
    alleleVector <- genoMatrix[,i] #extract allele column to a vector
    if(length(unique(alleleVector))==1){ #if that allele column has only one kind, remove it
      colToRemove <- c(colToRemove, i) #save column #th
    } else if(length(unique(alleleVector)) > 1){ #if allele column has more than one kind,
      mostAllele <- max(table(alleleVector)) #count number of an allele that has the largest number across individuals
      if( (mostAllele/noInds) > p){ #if the most allele has more than p% across individuals
        colToRemove <- c(colToRemove, i) #save column #th
      }
    }
  }
  #check if all loci variance less than p
  if(length(colToRemove)==0){
    stop("All loci have most abudnant alleles less than p threshold. Nothing to be removed.")
  }
  #Remove columns (alleles) that have low variance
  genoMatrix <- genoMatrix[-colToRemove]
  #Count new number of columns
  newNoAlleles <- ncol(genoMatrix) - 1
  #Count how many columns are removed
  noColumnsRemoved <- noAlleles - newNoAlleles
  
  #check if any locus is gone (all alleles gone) and edit the locus name vector
  newLocusNames <- NULL
  newAlleleLeft <- names(genoMatrix) #get remaining allele names
  for(j in 1:length(locusNames)){
    checking <- grep(paste0("^",locusNames[j],"_"), newAlleleLeft)
    if(!length(checking)==0){
      newLocusNames <- c(newLocusNames,locusNames[j])
    }
  }
  newNoLocus <- length(newLocusNames)
  #Print some message
  cat("\n  New data matrix has created! :)")
  cat(paste0("\n  New DataMatrix size: ",noInds," rows by ",ncol(genoMatrix)," columns"))
  cat(paste0("\n  ",noColumnsRemoved," columns (alleles) have been removed"))
  cat(paste0("\n  ",newNoLocus," loci remaining" ))
  cat("\n\n")
  return(list(genoMatrix, x[[2]], newLocusNames))#x[[2]]is individual ID vector, x[[3]] is locus name
  
}