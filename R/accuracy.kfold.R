#' Estimate assignment accuracies of K-fold cross-validation results
#'
#' This function allows you to estimate assignment accuracies of K-fold cross-validation results. The output results can be used to make assignment accuracy plots (use function accuracy.plot) and membership probability plot (use function membership.plot)
#' @param dir A character string to specify the folder that has your K-fold cross-validation results. A slash should be included at the end (e.g., dir="YourFolderName/").
#' @return This function outputs the results in a text file (a table). It can return a data frame when a returning object is specified.
#' @export
accuracy.kfold <- function(dir=NULL){
  #Read all "Out_*" file names in a specified directory
  fileName_vec <- list.files(path=dir, pattern="Out_*")
  fileName_vec <- sort(fileName_vec)
  noFiles <- length(fileName_vec)#count number of files
  #Read one of files and get pop names
  result01 <- read.table(paste0(dir,fileName_vec[1]), header=T)
  pops <- names(result01)[4:length(names(result01))] #read column name and get the pop names between 4th to last column
  noPops <- length(pops)#Number of pops
  #create vectors for saving data
  Var1 <- NULL; Var2 <- NULL
  KF <- NULL
  fold <- NULL
  train.loci <- NULL
  assign.rate.all <- NULL
  assign.rate.each <- as.data.frame(matrix(nrow=0,ncol=noPops),stringsAsFactors=F) #this will be an N tests by M pops dataframe
  #Analyze each assignment test result
  for(i in 1:noFiles){
    oneFileName <- unlist(strsplit(fileName_vec[i], split="_")) #split file name to 4 elements (e.g.,"Out"  "0.1"  "K3"  "1.txt");"0.1"is train.loci,"K3"is k fold,"1"is fold
    train.loci[i] <- oneFileName[2]
    KF[i] <- gsub("K","",oneFileName[3])#Remove "k" character in the "K3" string
    fold[i] <- unlist(strsplit(oneFileName[4],split=".txt"))
    df <- read.table(paste0(dir,fileName_vec[i]),header=T)
    #Calculate overall correct assignment rate
    #df$pred.pop <- factor(df$pred.pop, levels=levels(df$origin.pop))#ensure $pred.pop column has the same levels with $origin.pop
    levels(df$pred.pop) <- pops
    levels(df$origin.pop) <- pops
    ctable <- table(df$origin.pop,df$pred.pop)#make contingency table
    ftable <- as.data.frame(ctable)#convert table to data frame with frequency column
    totalSample <- sum(ftable$Freq)
    AllcorrectNo <- sum(subset(ftable,Var1==Var2)$Freq)
    assign.rate.all[i] <- AllcorrectNo/totalSample
    #calculate correct assignment rate each pop
    popCorrectRate_vec <- NULL
    for(p in pops){
      pop_size <- sum(subset(ftable,Var1==p)$Freq)
      popCorrectNo <- subset(subset(ftable,Var1==Var2), Var1==p)$Freq
      popCorrectRate <- popCorrectNo / pop_size
      popCorrectRate_vec <- c(popCorrectRate_vec, popCorrectRate)
    }
    #append correct assign rate of each pop as one row to data frame
    assign.rate.each[i,] <- popCorrectRate_vec
  }#for(i in 1:noFiles)
  #concatenate all data
  assign_rate_df <- cbind(KF, fold, train.loci, assign.rate.all, assign.rate.each)
  names(assign_rate_df)[5:ncol(assign_rate_df)] <- paste0("assign.rate.",pops)
  #output result
  write.table(assign_rate_df, file=paste0(dir,"Rate_of_",nrow(assign_rate_df),"_tests_",noPops,"_pops.txt"), quote=F, row.names=F )
  #Print some message to console
  cat("\n  Correct assignment rates were estimated!!")
  cat(paste0("\n  A total of ",nrow(assign_rate_df)," assignment tests for ",noPops," pops."))
  cat(paste0("\n  Results were also saved in a 'Rate_of_",nrow(assign_rate_df),"_tests_",noPops,"_pops.txt' file in the directory."))
  
  return(assign_rate_df)
  
}
