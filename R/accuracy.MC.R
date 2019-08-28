#' Estimate assignment accuracies of Monte-Carlo cross-validation results
#'
#' This function allows you to estimate assignment accuracies of Monte-Carlo cross-validation results. The output results can be used to make assignment accuracy plots (use function accuracy.plot).
#' @param dir A character string to specify the folder that has your Monte-Carlo cross-validation results. A slash should be included at the end (e.g., dir="YourFolderName/").
#' @return This function outputs the results in a text file (a table). It can return a data frame when a returning object is specified.
#' @export
accuracy.MC <- function(dir=NULL){
  #Read all "Out_*" file names in a specified directory
  fileName_vec <- list.files(path=dir, pattern="Out_*")
  fileName_vec <- sort(fileName_vec)
  noFiles <- length(fileName_vec)#count number of files
  #Read one of files and get pop names
  result01 <- read.table(paste0(dir,fileName_vec[1]), header=T, check.names = F)
  pops <- names(result01)[4:length(names(result01))] #read column name and get the pop names between 4th to last column
  noPops <- length(pops)#Number of pops
  #create vectors for saving data
  Var1 <- NULL; Var2 <- NULL
  train.inds <- NULL
  train.loci <- NULL
  iters <- NULL
  assign.rate.all <- NULL
  assign.rate.each <- as.data.frame(matrix(nrow=0,ncol=noPops),stringsAsFactors=F) #this will be an N tests by M pops dataframe
  #Analyze each assignment test result
  for(i in 1:noFiles){
    oneFileName <- unlist(strsplit(fileName_vec[i], split="_")) #split file name to 4 elements (e.g.,"Out"  "4"  "0.1"  "1.txt");"4"is train.level,"0.1"is fst.level,"1"is iter
    train.inds[i] <- oneFileName[2]
    train.loci[i] <- oneFileName[3]
    iters[i] <- unlist(strsplit(oneFileName[4],split=".txt"))
    df <- read.table(paste0(dir,fileName_vec[i]),header=T)
    #Calculate overall correct assignment rate
    #set levels of df$pred.pop and df$origin.pop to pops; handle cases when test individuals are not from every pop. 
    #df$pred.pop <- factor(df$pred.pop, levels=levels(df$origin.pop)) -- use in ver.1.1.4 and earlier; see issue about accurac.MC error 
    #levels(df$origin.pop) <- pops -- use in ver1.1.5/6 this causes issues that level of pred pop is ordered by pops and mismatch pop when test individuals not assigned to every pop
    #levels(df$pred.pop) <- pops
    df$origin.pop <- factor(df$origin.pop, levels=levels(factor(pops)))
    df$pred.pop <- factor(df$pred.pop, levels=levels(factor(pops)))
    ctable <- table(df$origin.pop,df$pred.pop)#make contingency table
    ftable <- as.data.frame(ctable)#convert table to data frame with frequency column
    totalSample <- sum(ftable$Freq)
    AllcorrectNo <- sum(subset(ftable,Var1==Var2)$Freq)
    assign.rate.all[i] <- AllcorrectNo/totalSample
    #calculate correct assignment rate each pop
    popCorrectRate_vec <- NULL
    for(p in pops){
      pop_size <- sum(subset(ftable,Var1==p)$Freq)
      if(pop_size==0){ 
        #if no individual from a pop was assigned to test set, set correct rate to 0
        popCorrectRate = 0
      }else{
        popCorrectNo <- subset(subset(ftable,Var1==Var2), Var1==p)$Freq
        popCorrectRate <- popCorrectNo / pop_size
      }
      popCorrectRate_vec <- c(popCorrectRate_vec, popCorrectRate)
    }
    #append correct assign rate of each pop as one row to data frame
    assign.rate.each[i,] <- popCorrectRate_vec
  }
  #concatenate all data
  assign_rate_df <- cbind(train.inds, train.loci, iters, assign.rate.all, assign.rate.each)
  names(assign_rate_df)[5:ncol(assign_rate_df)] <- paste0("assign.rate.",pops)
  #output result
  write.table(assign_rate_df, file=paste0(dir,"Rate_of_",nrow(assign_rate_df),"_tests_",noPops,"_pops.txt"), quote=F, row.names=F )
  #Print some message to console
  cat("\n  Correct assignment rates were estimated!!")
  cat(paste0("\n  A total of ",nrow(assign_rate_df)," assignment tests for ",noPops," pops."))
  cat(paste0("\n  Results were also saved in a 'Rate_of_",nrow(assign_rate_df),"_tests_",noPops,"_pops.txt' file in the directory."))
  
  return(assign_rate_df)
  
}