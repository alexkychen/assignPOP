#' Make an assignment maxtrix from cross-validation results
#' 
#' This function generates a pairwise assignment matrix with mean and variation of assignment accuracies estimated across all assignment tests.
#' @param dir A character string to specify the folder that has your cross-validation assignment results.
#' @param train.loci Choose your proportions of training loci used in Monte-Carlo or K-fold cross-validation. Default is "all".
#' @param train.inds Choose your numbers or proportions of training individuals used in Monte-Carlo cross-validation. Default is "all".
#' @param k.fold Choose the k fold values used in K-fold cross-validation. Default is "all".
#' @return The function returns a matrix in R console as well as a file named "assignment_matrix.txt" in the folder.
#' @export
#' 
assign.matrix <- function(dir=NULL, train.loci="all", train.inds="all", k.fold="all"){
  #read files in the folder
  fileName_vec <- list.files(path=dir, pattern="Out_*")
  fileName_vec <- sort(fileName_vec)
  noFiles <- length(fileName_vec)#count number of files
  #Read one of files and get pop names
  result01 <- read.table(paste0(dir,fileName_vec[1]), header=T)
  pops <- names(result01)[4:length(names(result01))] #read column name and get the pop names between 4th to last column
  noPops <- length(pops)#Number of pops
  
  #claim variables
  fileName_select <- NULL
  
  #check whether results are from MC or kfold
  if(grepl("K", fileName_vec[1])){ #check if it is from kfold results
    #Grap selected train.loci & kfold files
    if((train.loci != "all")&(k.fold != "all")){
      for(i in 1:noFiles){
        oneFileName <- unlist(strsplit(fileName_vec[i], split="_"))
        if(oneFileName[2] %in% train.loci){
          if(gsub("K","",oneFileName[3]) %in% k.fold ){
            fileName_select <- c(fileName_select, fileName_vec[i])
          }
        }
      }
    }else if((train.loci == "all")&(k.fold != "all")){
      for(i in 1:noFiles){
        oneFileName <- unlist(strsplit(fileName_vec[i], split="_"))
        if(gsub("K","",oneFileName[3]) %in% k.fold){
          fileName_select <- c(fileName_select, fileName_vec[i])
        }
      }
    }else if((train.loci != "all")&(k.fold == "all")){
      for(i in 1:noFiles){
        oneFileName <- unlist(strsplit(fileName_vec[i], split="_"))
        if(oneFileName[2] %in% train.loci){
          fileName_select <- c(fileName_select, fileName_vec[i])
        }
      }
    }else{
      fileName_select <- fileName_vec
    }
    
  }else{ #else it is from MC results
    #Grap selected train.loci & train.inds files
    if((train.loci != "all") & (train.inds != "all")){
      for(i in 1:noFiles){
        oneFileName <- unlist(strsplit(fileName_vec[i], split="_"))
        if(oneFileName[2] %in% train.inds){
          if(oneFileName[3] %in% train.loci){
            fileName_select <- c(fileName_select, fileName_vec[i])
          }
        }
      }
    }else if((train.loci == "all")&(train.inds != "all")){
      for(i in 1:noFiles){
        oneFileName <- unlist(strsplit(fileName_vec[i], split="_"))
        if(oneFileName[2] %in% train.inds){
          fileName_select <- c(fileName_select, fileName_vec[i])
        }
      }
    }else if((train.loci != "all")&(train.inds == "all")){
      for(i in 1:noFiles){
        oneFileName <- unlist(strsplit(fileName_vec[i], split="_"))
        if(oneFileName[3] %in% train.loci){
          fileName_select <- c(fileName_select, fileName_vec[i])
        }
      }
    }else{
      fileName_select <- fileName_vec
    }
  }
  #check number of selected files
  noFiles_select <- length(fileName_select)
  
  freq_df <- as.data.frame(matrix(nrow=noPops*noPops,ncol=0),stringsAsFactors=F)
  #Read each file and process data
  for(j in 1:noFiles_select){
    df <- read.table(paste0(dir,fileName_select[j]),header=T)
    df$pred.pop <- factor(df$pred.pop, levels=levels(df$origin.pop))#ensure $pred.pop column has the same levels with $origin.pop
    ctable <- table(df$origin.pop,df$pred.pop)
    #calcuate assignment rate;convert number to rate
    for(k in 1:noPops){
      ctable[k,] <- round(ctable[k,]/sum(ctable[k,]), digits=2)
    }
    ftable <- as.data.frame(ctable)#convert table to data frame with frequency column
    freq_df <- cbind(freq_df, ftable$Freq)
  }
  #estimate mean and sd of assignment rate
  assign_mean <- round(apply(freq_df,1,mean),digits=2)
  assign_sd <- round(apply(freq_df,1,sd),digits=2)
  #create dataframe for saving mean and sd
  assign_df <- ftable[c(1,2)];colnames(assign_df) <- c("origin","assignment")
  #create dataframe of assignment mean
  assign_df <- cbind(assign_df,assign_mean, assign_sd)
  
  #print contingency table (assignment mean)
  cat("Mean of assignment rates across ");cat(noFiles_select);cat(" tests.\n")
  xtabs(assign_mean ~ origin + assignment, data=assign_df)
  #print contingency table (assignment sd)
  cat("\nStandard deviation of assignment rates across ");cat(noFiles_select);cat(" tests.")
  xtabs(assign_sd ~ origin + assignment, data=assign_df)
  
  
}#end