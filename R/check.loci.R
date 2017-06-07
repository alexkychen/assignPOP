#' Check which loci frequently have high Fst across training sets
#'
#' This function reads through training locus file for each assignment test and counts the frequency of those loci and outputs the results in a text file.
#' @param dir A character string to specify the folder with your cross-validation results. A slash should be entered at the end.
#' @param top.loci An integer to specify how many top informative loci to output.
#' @return This function output the results in a text file. It includes the top N informative loci in N rows, and each row has a list of loci sorted by its occurrence.
#' @import stringr
#' @export
#'
check.loci <- function(dir=NULL, top.loci=20){
  #Read all "Out_*" file names in a specified directory
  fileName_vec <- list.files(path=dir, pattern="Loci_*")
  fileName_vec <- sort(fileName_vec)
  noFiles <- length(fileName_vec)#count number of files

  ana_info <- readLines(paste0(dir,"AnalysisInfo.txt"))
  locusSampleMethod <- ana_info[8]
  if(grepl(pattern="fst", locusSampleMethod)){ #check if results from prior fst sampling method
    train_level <- NULL
    for(i in 1:noFiles){
      oneFileName <- unlist(strsplit(fileName_vec[i], split="_")) #split file name to 4 elements (e.g.,"Loci"  "4"  "0.1"  "1.txt");"4"is train.level,"0.1"is fst.level,"1"is iter
      train_level[i] <- oneFileName[2]
      #fst.level[i] <- oneFileName[3]
      #iters[i] <- unlist(strsplit(oneFileName[4],split=".txt"))
    }
    trainLevels <- unique(train_level)#check out unique training individual levels
    #see if only one level of training individuals, if so, skip question
    if(length(trainLevels)==1){
      ans_0 <- "all"
    }else{ # if levels of training individuals are more than one
      cat(paste0("\n  ",length(trainLevels)," levels of training individuals are found."))
      cat("\n  Which levels would you like to check? (separate levels by a whitespace if multiple)")
      cat("\n  Options: ");cat(trainLevels, sep=", ");cat(", or all")
      ans_0 <- readline("\n  enter here: ")
    }
    
    if(grepl(pattern="ALL",toupper(ans_0))){ #if answer is "all"
      cat(paste0("Loci occur in top ",top.loci," high Fst across all training data\n"), file=paste0(dir,"High_Fst_Locus_Freq.txt")) #Output the following result to text file
      lociMatrix <- NULL
      #Read through each Loci_ file.
      for(i in 1:noFiles){
        oneFileLoci <- readLines(paste0(dir,fileName_vec[i]))[1:top.loci] #Read one Loci_ file for [top.loci]
        lociMatrix <- cbind(lociMatrix, oneFileLoci)
      }
      #colnames(lociMatrix) <- fileName_vec #rename column for lociMatrix (minor)
      #lociMatrix_freq <- as.data.frame(matrix(nrow=0, ncol=3), stringsAsFactors=F)
      #Check out top # high Fst loci
      for(n in 1:top.loci){
        locName_sort <- names(sort(table(lociMatrix[n,]), decreasing=T)) #Read lociMatrix each row (n-th row) from the top row (highest Fst), and get locus name sorted by freq
        locFreq_sort <- as.character(sort(table(lociMatrix[n,]), decreasing=T)) #Get freq of locus, sorted
        cat(paste0("top.",n,"(",length(locName_sort),"): "), file=paste0(dir,"High_Fst_Locus_Freq.txt"), append=T)
        for(j in 1:length(locName_sort)){
          cat(locName_sort[j], file=paste0(dir,"High_Fst_Locus_Freq.txt"), append=T)
          cat(paste0("(",locFreq_sort[j],"), ") , file=paste0(dir,"High_Fst_Locus_Freq.txt"), append=T)
        }
        cat("\n",file=paste0(dir,"High_Fst_Locus_Freq.txt"), append=T )
      }
      #print some message to console when done
      cat("\n  Results were saved in a 'High_Fst_Locus_Freq.txt' file in the directory.")

      #Process data if only a portion is selected
    }else if(!grepl(pattern="ALL",toupper(ans_0))){ #else if answer contains no "all"
      #parse user input
      ans_0 <- str_trim(ans_0, side="both")
      ans_0 <- unlist(strsplit(ans_0," "))
      newFileName_vec <- NULL
      #Check out selected train level files
      for(k in ans_0){
        #Look through Loci_ files
        for(i in 1:noFiles){
          oneFileName <- unlist(strsplit(fileName_vec[i], split="_"))
          if(oneFileName[2]==k){ # if train level matches, append the file name to newFileName_vec
            newFileName_vec <- c(newFileName_vec, fileName_vec[i])
          }
        }
      }
      noNewFiles <- length(newFileName_vec)#count number of selected files
      #Output results to text file and print title
      cat(paste0("Loci occur in top ",top.loci," high Fst across selected ( train inds level: "), file=paste0(dir,"High_Fst_Locus_Freq.txt"))
      for(a in ans_0){
        cat(paste0(a," "), file=paste0(dir,"High_Fst_Locus_Freq.txt"), append=T)
      }
      cat(") training data\n", file=paste0(dir,"High_Fst_Locus_Freq.txt"), append=T )
      lociMatrix <- NULL
      #Read through selected Loci_ file.
      for(s in 1:noNewFiles){
        oneFileLoci <- readLines(paste0(dir,newFileName_vec[s]))[1:top.loci] #Read one Loci_ file for [top.loci]
        lociMatrix <- cbind(lociMatrix, oneFileLoci)
      }
      #Check out top # high Fst loci
      for(n in 1:top.loci){
        locName_sort <- names(sort(table(lociMatrix[n,]), decreasing=T)) #Read lociMatrix each row (n-th row) from the top row (highest Fst), and get locus name sorted by freq
        locFreq_sort <- as.character(sort(table(lociMatrix[n,]), decreasing=T)) #Get freq of locus, sorted
        cat(paste0("top.",n,"(",length(locName_sort),"): "), file=paste0(dir,"High_Fst_Locus_Freq.txt"), append=T)
        for(j in 1:length(locName_sort)){
          cat(locName_sort[j], file=paste0(dir,"High_Fst_Locus_Freq.txt"), append=T)
          cat(paste0("(",locFreq_sort[j],"), ") , file=paste0(dir,"High_Fst_Locus_Freq.txt"), append=T)
        }
        cat("\n",file=paste0(dir,"High_Fst_Locus_Freq.txt"), append=T )
      }
      #Print some message in console when done
      cat("\n  Results were saved in a 'High_Fst_Locus_Freq.txt' file in the directory.")
    }

  }else if(grepl(pattern="random", locusSampleMethod)){
    #save for future development
    message('Checking loci of results generated from random sampling is not useful.\nUse loci.sample = "Fst" in your analyses.')
  }
}
