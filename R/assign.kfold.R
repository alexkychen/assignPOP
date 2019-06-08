#' Population assignment test using K-fold cross-validation
#'
#' This function employs K-fold cross-validation for assignment tests. The results help estimate membership probabilities of every individual. It accepts genetic-only [object returned from read.genpop() or reducel.allele()], integrated [object returned from compile.data()], or non-genetic [R data frame with header] data as input, and outputs results to text files. Several built-in options are provided. See below for more details.
#' @param x An input object which should be the object (list) returned from the function read.genpop(), reduce.allele(), or compile.data(). It could also be a data frame (with column name) returned from read.csv() or read.table() if you're analyzing non-genetic data, such as morphormetrics, chemistry data. The non-genetic data frame should have sample ID in the first column and population label in the last column. 
#' @param k.fold The number of groups to be divided for each population. Use a numeric vector to specify multiple sets of k-folds.
#' @param train.loci The proportion (float between 0 and 1) of loci to be used as training data. Use a numeric vector to specify multiple sets of training loci. This argument will be ignored if you're analyzing non-genetic data.
#' @param loci.sample Locus sampling method, "fst" or "random". If loci.sample="fst" (default) and train.loci=0.1, it means that top 10 percent of high Fst loci will be sampled as training loci. On the other hand, if loci.sample="random", then random 10 percent of loci will be sampled as training loci. This argument will be ignored if you're analyzing non-genetic data.
#' @param dir A character string to specify the folder name for saving output files. A slash at the end must be included (e.g., dir="YourFolderName/"). Otherwise, the files will be saved under your working directory.
#' @param scaled A logical variable (TRUE or FALSE) to specify whether to center (make mean of each feature to 0) and scale (make standard deviation of each feature to 1) the entire dataset before performing PCA and cross-validation. Default is FALSE. As genetic data has converted to numeric data between 0 and 1, to scale or not to scale the genetic data should not be critical. However, it is recommended to set scaled=TRUE when integrated data contains various scales of features.
#' @param pca.method Either a character string ("mixed", "independent", or "original") or logical variable (TRUE or FALSE) to specify how to perform PCA on non-genetic data (PCA is always performed on genetic data). The character strings are used when analyzing integrated (genetic plus non-genetic) data. If using "mixed" (default), PCA is perfromed across the genetic and non-genetic data, resulting in each PC summarizing mixed variations of genetic and non-genetic data. If using "independent", PCA is independently performed on non-genetic data. Genetic PCs and non-genetic PCs are then used as new features. If using "original", original non-genetic data and genetic PCs are used as features. The logical variable is used when analyzing non-genetic data alone. If TRUE, it performs PCA on the training data and applys the loadings to the test data. Scores of training and test data will be used as new features.
#' @param pca.PCs A criterion to retain number of PCs. By default, it uses Kaiser-Guttman criterion that any PC has the eigenvalue greater than 1 will be retained as the new variable/feature. Users can set an integer to specify the number of PCs to be retained.
#' @param pca.loadings A logical variable (False or True) to determine whether it prints the loadings of training data to output text files. Default is False, if set True, the overall output files could be large.
#' @param model A character string to specify which classifier to use for creating predictive models. The current options include "lda", "svm", "naiveBayes", "tree", and "randomForest".
#' @param svm.kernel A character string to specify which kernel to be used when using "svm" classifier.
#' @param svm.cost A number to specify the cost for "svm" method.
#' @param ntree A integer to specify how many trees to build when using "randomForest" method.
#' @param multiprocess A logical variable to determine whether using multiprocess. Default is TRUE. If set FALSE, it will only use single core to run the program. 
#' @param processors The number of processors to be used for parallel running. By default, it uses N-1 processors in your computer.
#' @param skipQ A logical variable to determine whether prompting interactive dialogue when analyzing non-genetic data. If set TRUE, default data type and original values of non-genetic data will be used.
#' @param ... Other arguments that could be potentially used for various models
#' @return You don't need to specify a name for the returned object when using this function. It automatically outputs results in text files to your designated folder.
#' @import stringr
#' @import foreach
#' @importFrom reshape2 melt
#' @importFrom caret createFolds
#' @importFrom MASS lda
#' @importFrom e1071 svm naiveBayes
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom tree tree
#' @importFrom randomForest randomForest importance
#' @importFrom utils write.table packageVersion
#' @importFrom stats model.matrix prcomp predict
#' @export
#'
assign.kfold <- function(x, k.fold = c(3,4,5), train.loci=c(0.1,0.25,0.5, 1), loci.sample="fst", dir=NULL, scaled=FALSE,
                         pca.method="mixed", pca.PCs="kaiser-guttman", pca.loadings=F,
                         model="svm", svm.kernel="linear", svm.cost=1, ntree=50, processors=999, multiprocess=FALSE, skipQ=FALSE, ...){
  #check if dir is correctly entered
  if(is.null(dir)){
    stop("Please provide a folder name ending with '/' in argument 'dir' ")
  }else if(substr(dir, start=nchar(dir), stop=nchar(dir))!="/"){
    stop("Please put a forward slash '/' in the end of your folder name (in argument 'dir'). ")
  }
  #check data type
  if(!is.data.frame(x)){#check if input x is a list returned from read.genpop(), reduce.allele(), or compile.data()
    #Analyze genetic or integrated data
    #checking pca.method
    if(!is.character(pca.method)){ #if pca.method is not character string, print message and stop analyzing
      stop("Please specify a correct parameter, 'mixed' or 'independent' or 'original', for argument 'pca.method' ")
    }
    popNames_vector <- NULL; i <- NULL #some NULL variable to handle R CMD check
    genoMatrix <- x[[1]]
    popSizes <- table(genoMatrix$popNames_vector)#get number of individual for each pop in table
    pops <- names(popSizes)#Get what pops in data
    noPops <- length(popSizes)#count number of pops
    locusNames <- x[[3]]#Get locus name
    noLocus <- length(locusNames)
    alleleName <- colnames(genoMatrix) #This includes last column (popNames_vector), and possible non-genetic data
    #Check if pop size and k fold value can work
    if(max(k.fold) > min(popSizes)){
      stop("Max. K value is greater than a small pop. Please adjust your k.fold setting or increase sample size.")
    }
    #Create a folder to save outfiles
    dir.create(file.path(dir))
    #Detect CPU core/thread numbers and use n-1 threads for parallel ananlysis
    if(multiprocess){
      maxCores <- detectCores()-1
      if (processors <= maxCores & processors > 0){
        #cl <- makeCluster(processors)
        #registerDoParallel(cl,cores=processors)
        cat("\n  Parallel computing is on. Analyzing data using",processors,"cores/threads of CPU...\n")
      }else {
        #cl <- makeCluster(maxCores)
        #registerDoParallel(cl,cores=maxCores)
        cat("\n  Parallel computing is on. Analyzing data using",maxCores,"cores/threads of CPU...\n")
      }
    }else {
      cat("\n  Parallele computing is off. Analyzing data using 1 CPU core...\n")
    }
    #Determine if object x includes non-genetic data. If so, length of x (list) will be 5. If not (only genetic data), the length will be 3.
    if(length(x)==3){#Process genetics-only data
      datatype <- "genetics";noVars <- 0
      #loop each k fold
      for(k in k.fold){
        #create fold index based on label (genoMatrix$popNames_vector)
        fold_index <- createFolds(genoMatrix$popNames_vector, k=k) #this is a list containing k folds of index, e.g., fold_index$Fold1, $Fold2...
        #take the K-th fold as test set; remaing folds as training set
        if(multiprocess){
          if (processors <= maxCores & processors > 0){
            cl <- makeCluster(processors)
            registerDoParallel(cl,cores=processors)
          }else {
            cl <- makeCluster(maxCores)
            registerDoParallel(cl,cores=maxCores)
          }
          foreach(i=1:k, .export=c("Fsts","perform.PCA"), .packages=c("e1071","klaR","MASS","tree","randomForest")) %dopar% {
            trainSetMatrix <- genoMatrix[-fold_index[[i]],]
            testIndID <- x[[2]][fold_index[[i]]]
            #check if fst prior for sampling loci
            if(loci.sample=="fst"){
              train_X <- list(trainSetMatrix, row.names(trainSetMatrix), x[[3]])#create a new x list for Fsts function; x[[3]] is locus name
              fstTable <- Fsts(train_X)#Estimate locus Fst for training data
              orderFst <- order(fstTable$Fst, decreasing=T)#Order the loci index by Fst value(from highest)
              for(f in train.loci){
                trainLocusIndex_fstTable <- orderFst[1:round(noLocus*f)]#Get training locus index (from fstTable)
                trainLocusName <- as.character(fstTable[trainLocusIndex_fstTable,]$Locus)#Get training locus name
                trainLocusIndex_genoMatrix <- NULL #create a train locus index for genoMatrix to be extracted
                for(m in 1:length(trainLocusName)){
                  tempAlleleIndex <- grep(pattern=paste0(trainLocusName[m],"_"), alleleName)#alleleName is colnames(genoMatrix)
                  trainLocusIndex_genoMatrix <- c(trainLocusIndex_genoMatrix,tempAlleleIndex )
                }
                trainLocusIndex_genoMatrix <- sort(trainLocusIndex_genoMatrix)
                #genoMatrix_trainLoci <- genoMatrix[,c(trainLocusIndex_genoMatrix, ncol(genoMatrix))]#resahpe the genoMatrix data, comprising only training loci
                genoMatrix_trainLoci <- genoMatrix[,c(trainLocusIndex_genoMatrix)]
                #Scale and center the data if "scaled=T"
                if(scaled){
                  genoMatrix_trainLoci <- as.data.frame(scale(genoMatrix_trainLoci))
                }
                #Add the pop.Names to the last column
                genoMatrix_trainLoci <- cbind(genoMatrix_trainLoci, genoMatrix[,ncol(genoMatrix)]);colnames(genoMatrix_trainLoci)[ncol(genoMatrix_trainLoci)] <- "popNames_vector"
                ##
                trainSetData <- genoMatrix_trainLoci[-fold_index[[i]],]#Get the training set data (training individuals/loci)
                testSetData <- genoMatrix_trainLoci[fold_index[[i]],]#Get the test set data (test individuals/loci)
                ##
                #Peform PCA on training
                PCA_results <- perform.PCA(trainSetData[,1:ncol(trainSetData)-1], method=pca.PCs) #Run PCA without label column
                loadings <- PCA_results[[1]] #loadings (coefficience) of variables and PCs; apply this to test data
                trainSetData_PC <- as.data.frame(PCA_results[[2]])
                trainSetData_PC <- cbind(trainSetData_PC, trainSetData$popNames_vector) ##Will be used for building predicting models
                colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popNames_vector"
                #Convert test data to PC variables based on training's loadings
                testSetData_matrix <- as.matrix(testSetData[,1:ncol(testSetData)-1])
                testSetData_PC <- as.data.frame(testSetData_matrix %*% loadings)
                testSetData_PC <- cbind(testSetData_PC, testSetData$popNames_vector)
                colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popNames_vector"
                #Use training to build models and test on test individuals
                if(model=="svm"){
                  svm.fit <- svm(popNames_vector ~ ., data=trainSetData_PC, kernel=svm.kernel, cost=svm.cost, prob=T, scale=F)
                  svm.pred <- predict(svm.fit, testSetData_PC, type="class",prob=T)
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(svm.pred), attr(svm.pred,"probabilities"))#combine output to data frame
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="knn"){
                  
                }else if (model=="lda"){
                  lda.fit <- lda(popNames_vector ~ ., data=trainSetData_PC)
                  lda.pred <- predict(lda.fit, testSetData_PC)
                  lda.pred.class <- lda.pred$class
                  lda.pred.prob <- lda.pred$posterior
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(lda.pred.class), as.data.frame(lda.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="naiveBayes"){
                  nby.model <- naiveBayes(popNames_vector ~ ., data=trainSetData_PC)
                  nby.pred.class <- predict(nby.model,testSetData_PC,type="class")
                  nby.pred.prob <- predict(nby.model,testSetData_PC,type="raw")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(nby.pred.class), as.data.frame(nby.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="tree"){
                  tree.model <- tree(popNames_vector ~ ., data=trainSetData_PC)
                  tree.pred.class <- predict(tree.model,testSetData_PC,type="class")
                  tree.pred.prob <- predict(tree.model,testSetData_PC,type="vector")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(tree.pred.class), as.data.frame(tree.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  tree_node <- as.character(summary(tree.model)$used)#Get loci used at tree node
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  cat("Loci used at tree node",tree_node,file=paste0(dir,"Loci_treenode_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="randomForest"){
                  rf.model <- randomForest(popNames_vector ~ ., data=trainSetData_PC, ntree=ntree, importance=T)
                  rf.pred.class <- predict(rf.model,testSetData_PC,type="response")
                  rf.pred.prob <- predict(rf.model,testSetData_PC,type="prob")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(rf.pred.class), as.data.frame(rf.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  #Output "importance"(Gini index) of loci to files;see randomForest package's "importance" function
                  write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Loci_importance_",f,"_K",k,"_",i,".txt"), quote=F, row.names=T)
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }
              }#for(f in train.loci)
            }else if(loci.sample=="random"){
              for(f in train.loci){
                tempLocusIndex <- sort(sample(1:noLocus, round(noLocus*f)))
                trainLocusName <- locusNames[tempLocusIndex]
                trainLocusIndex_genoMatrix <- NULL #create a train locus index for genoMatrix to be extracted
                for(m in 1:length(trainLocusName)){
                  tempAlleleIndex <- grep(pattern=paste0(trainLocusName[m],"_"), alleleName)
                  trainLocusIndex_genoMatrix <- c(trainLocusIndex_genoMatrix,tempAlleleIndex )
                }
                trainLocusIndex_genoMatrix <- sort(trainLocusIndex_genoMatrix)
                #resahpe the genoMatrix data, comprising only training loci
                genoMatrix_trainLoci <- genoMatrix[,c(trainLocusIndex_genoMatrix)]
                #Scale and center the data if "scaled=T"
                if(scaled){
                  genoMatrix_trainLoci <- as.data.frame(scale(genoMatrix_trainLoci))
                }
                #Add the pop.Names to the last column
                genoMatrix_trainLoci <- cbind(genoMatrix_trainLoci, genoMatrix[,ncol(genoMatrix)]);colnames(genoMatrix_trainLoci)[ncol(genoMatrix_trainLoci)] <- "popNames_vector"
                ##
                trainSetData <- genoMatrix_trainLoci[-fold_index[[i]],]#Get the training set data (training individuals/loci)
                testSetData <- genoMatrix_trainLoci[fold_index[[i]],]#Get the test set data (test individuals/loci)
                ##
                #Peform PCA on training
                PCA_results <- perform.PCA(trainSetData[,1:ncol(trainSetData)-1], method=pca.PCs) #Run PCA without label column
                loadings <- PCA_results[[1]] #loadings (coefficience) of variables and PCs; apply this to test data
                trainSetData_PC <- as.data.frame(PCA_results[[2]])
                trainSetData_PC <- cbind(trainSetData_PC, trainSetData$popNames_vector) ##Will be used for building predicting models
                colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popNames_vector"
                #Convert test data to PC variables based on training's loadings
                testSetData_matrix <- as.matrix(testSetData[,1:ncol(testSetData)-1])
                testSetData_PC <- as.data.frame(testSetData_matrix %*% loadings)
                testSetData_PC <- cbind(testSetData_PC, testSetData$popNames_vector)
                colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popNames_vector"
                #Use training to build models and test on test individuals
                if(model=="svm"){
                  svm.fit <- svm(popNames_vector ~ ., data=trainSetData_PC, kernel=svm.kernel, cost=svm.cost, prob=T, scale=F)
                  svm.pred <- predict(svm.fit, testSetData_PC, type="class",prob=T)
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(svm.pred), attr(svm.pred,"probabilities"))#combine output to data frame
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="knn"){
                  
                }else if (model=="lda"){
                  lda.fit <- lda(popNames_vector ~ ., data=trainSetData_PC)
                  lda.pred <- predict(lda.fit, testSetData_PC)
                  lda.pred.class <- lda.pred$class
                  lda.pred.prob <- lda.pred$posterior
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(lda.pred.class), as.data.frame(lda.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="naiveBayes"){
                  nby.model <- naiveBayes(popNames_vector ~ ., data=trainSetData_PC)
                  nby.pred.class <- predict(nby.model,testSetData_PC,type="class")
                  nby.pred.prob <- predict(nby.model,testSetData_PC,type="raw")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(nby.pred.class), as.data.frame(nby.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="tree"){
                  tree.model <- tree(popNames_vector ~ ., data=trainSetData_PC)
                  tree.pred.class <- predict(tree.model,testSetData_PC,type="class")
                  tree.pred.prob <- predict(tree.model,testSetData_PC,type="vector")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(tree.pred.class), as.data.frame(tree.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  tree_node <- as.character(summary(tree.model)$used)#Get loci used at tree node
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  cat("Loci used at tree node",tree_node,file=paste0(dir,"Loci_treenode_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="randomForest"){
                  rf.model <- randomForest(popNames_vector ~ ., data=trainSetData_PC, ntree=ntree, importance=T)
                  rf.pred.class <- predict(rf.model,testSetData_PC,type="response")
                  rf.pred.prob <- predict(rf.model,testSetData_PC,type="prob")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(rf.pred.class), as.data.frame(rf.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  #Output "importance"(Gini index) of loci to files;see randomForest package's "importance" function
                  write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Loci_importance_",f,"_K",k,"_",i,".txt"), quote=F, row.names=T)
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }
              }# for(f in train.loci)
            }#else if(loci.sample="random")
          }# foreach(i=1:k,...
          stopCluster(cl)
        }else{
          for(i in 1:k){
            trainSetMatrix <- genoMatrix[-fold_index[[i]],]
            testIndID <- x[[2]][fold_index[[i]]]
            #check if fst prior for sampling loci
            if(loci.sample=="fst"){
              train_X <- list(trainSetMatrix, row.names(trainSetMatrix), x[[3]])#create a new x list for Fsts function; x[[3]] is locus name
              fstTable <- Fsts(train_X)#Estimate locus Fst for training data
              orderFst <- order(fstTable$Fst, decreasing=T)#Order the loci index by Fst value(from highest)
              for(f in train.loci){
                trainLocusIndex_fstTable <- orderFst[1:round(noLocus*f)]#Get training locus index (from fstTable)
                trainLocusName <- as.character(fstTable[trainLocusIndex_fstTable,]$Locus)#Get training locus name
                trainLocusIndex_genoMatrix <- NULL #create a train locus index for genoMatrix to be extracted
                for(m in 1:length(trainLocusName)){
                  tempAlleleIndex <- grep(pattern=paste0(trainLocusName[m],"_"), alleleName)#alleleName is colnames(genoMatrix)
                  trainLocusIndex_genoMatrix <- c(trainLocusIndex_genoMatrix,tempAlleleIndex )
                }
                trainLocusIndex_genoMatrix <- sort(trainLocusIndex_genoMatrix)
                #genoMatrix_trainLoci <- genoMatrix[,c(trainLocusIndex_genoMatrix, ncol(genoMatrix))]#resahpe the genoMatrix data, comprising only training loci
                genoMatrix_trainLoci <- genoMatrix[,c(trainLocusIndex_genoMatrix)]
                #Scale and center the data if "scaled=T"
                if(scaled){
                  genoMatrix_trainLoci <- as.data.frame(scale(genoMatrix_trainLoci))
                }
                #Add the pop.Names to the last column
                genoMatrix_trainLoci <- cbind(genoMatrix_trainLoci, genoMatrix[,ncol(genoMatrix)]);colnames(genoMatrix_trainLoci)[ncol(genoMatrix_trainLoci)] <- "popNames_vector"
                ##
                trainSetData <- genoMatrix_trainLoci[-fold_index[[i]],]#Get the training set data (training individuals/loci)
                testSetData <- genoMatrix_trainLoci[fold_index[[i]],]#Get the test set data (test individuals/loci)
                ##
                #Peform PCA on training
                PCA_results <- perform.PCA(trainSetData[,1:ncol(trainSetData)-1], method=pca.PCs) #Run PCA without label column
                loadings <- PCA_results[[1]] #loadings (coefficience) of variables and PCs; apply this to test data
                trainSetData_PC <- as.data.frame(PCA_results[[2]])
                trainSetData_PC <- cbind(trainSetData_PC, trainSetData$popNames_vector) ##Will be used for building predicting models
                colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popNames_vector"
                #Convert test data to PC variables based on training's loadings
                testSetData_matrix <- as.matrix(testSetData[,1:ncol(testSetData)-1])
                testSetData_PC <- as.data.frame(testSetData_matrix %*% loadings)
                testSetData_PC <- cbind(testSetData_PC, testSetData$popNames_vector)
                colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popNames_vector"
                #Use training to build models and test on test individuals
                if(model=="svm"){
                  svm.fit <- svm(popNames_vector ~ ., data=trainSetData_PC, kernel=svm.kernel, cost=svm.cost, prob=T, scale=F)
                  svm.pred <- predict(svm.fit, testSetData_PC, type="class",prob=T)
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(svm.pred), attr(svm.pred,"probabilities"))#combine output to data frame
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="knn"){
                  
                }else if (model=="lda"){
                  lda.fit <- lda(popNames_vector ~ ., data=trainSetData_PC)
                  lda.pred <- predict(lda.fit, testSetData_PC)
                  lda.pred.class <- lda.pred$class
                  lda.pred.prob <- lda.pred$posterior
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(lda.pred.class), as.data.frame(lda.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="naiveBayes"){
                  nby.model <- naiveBayes(popNames_vector ~ ., data=trainSetData_PC)
                  nby.pred.class <- predict(nby.model,testSetData_PC,type="class")
                  nby.pred.prob <- predict(nby.model,testSetData_PC,type="raw")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(nby.pred.class), as.data.frame(nby.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="tree"){
                  tree.model <- tree(popNames_vector ~ ., data=trainSetData_PC)
                  tree.pred.class <- predict(tree.model,testSetData_PC,type="class")
                  tree.pred.prob <- predict(tree.model,testSetData_PC,type="vector")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(tree.pred.class), as.data.frame(tree.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  tree_node <- as.character(summary(tree.model)$used)#Get loci used at tree node
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  cat("Loci used at tree node",tree_node,file=paste0(dir,"Loci_treenode_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="randomForest"){
                  rf.model <- randomForest(popNames_vector ~ ., data=trainSetData_PC, ntree=ntree, importance=T)
                  rf.pred.class <- predict(rf.model,testSetData_PC,type="response")
                  rf.pred.prob <- predict(rf.model,testSetData_PC,type="prob")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(rf.pred.class), as.data.frame(rf.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  #Output "importance"(Gini index) of loci to files;see randomForest package's "importance" function
                  write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Loci_importance_",f,"_K",k,"_",i,".txt"), quote=F, row.names=T)
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }
              }#for(f in train.loci)
            }else if(loci.sample=="random"){
              for(f in train.loci){
                tempLocusIndex <- sort(sample(1:noLocus, round(noLocus*f)))
                trainLocusName <- locusNames[tempLocusIndex]
                trainLocusIndex_genoMatrix <- NULL #create a train locus index for genoMatrix to be extracted
                for(m in 1:length(trainLocusName)){
                  tempAlleleIndex <- grep(pattern=paste0(trainLocusName[m],"_"), alleleName)
                  trainLocusIndex_genoMatrix <- c(trainLocusIndex_genoMatrix,tempAlleleIndex )
                }
                trainLocusIndex_genoMatrix <- sort(trainLocusIndex_genoMatrix)
                #resahpe the genoMatrix data, comprising only training loci
                genoMatrix_trainLoci <- genoMatrix[,c(trainLocusIndex_genoMatrix)]
                #Scale and center the data if "scaled=T"
                if(scaled){
                  genoMatrix_trainLoci <- as.data.frame(scale(genoMatrix_trainLoci))
                }
                #Add the pop.Names to the last column
                genoMatrix_trainLoci <- cbind(genoMatrix_trainLoci, genoMatrix[,ncol(genoMatrix)]);colnames(genoMatrix_trainLoci)[ncol(genoMatrix_trainLoci)] <- "popNames_vector"
                ##
                trainSetData <- genoMatrix_trainLoci[-fold_index[[i]],]#Get the training set data (training individuals/loci)
                testSetData <- genoMatrix_trainLoci[fold_index[[i]],]#Get the test set data (test individuals/loci)
                ##
                #Peform PCA on training
                PCA_results <- perform.PCA(trainSetData[,1:ncol(trainSetData)-1], method=pca.PCs) #Run PCA without label column
                loadings <- PCA_results[[1]] #loadings (coefficience) of variables and PCs; apply this to test data
                trainSetData_PC <- as.data.frame(PCA_results[[2]])
                trainSetData_PC <- cbind(trainSetData_PC, trainSetData$popNames_vector) ##Will be used for building predicting models
                colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popNames_vector"
                #Convert test data to PC variables based on training's loadings
                testSetData_matrix <- as.matrix(testSetData[,1:ncol(testSetData)-1])
                testSetData_PC <- as.data.frame(testSetData_matrix %*% loadings)
                testSetData_PC <- cbind(testSetData_PC, testSetData$popNames_vector)
                colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popNames_vector"
                #Use training to build models and test on test individuals
                if(model=="svm"){
                  svm.fit <- svm(popNames_vector ~ ., data=trainSetData_PC, kernel=svm.kernel, cost=svm.cost, prob=T, scale=F)
                  svm.pred <- predict(svm.fit, testSetData_PC, type="class",prob=T)
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(svm.pred), attr(svm.pred,"probabilities"))#combine output to data frame
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="knn"){
                  
                }else if (model=="lda"){
                  lda.fit <- lda(popNames_vector ~ ., data=trainSetData_PC)
                  lda.pred <- predict(lda.fit, testSetData_PC)
                  lda.pred.class <- lda.pred$class
                  lda.pred.prob <- lda.pred$posterior
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(lda.pred.class), as.data.frame(lda.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="naiveBayes"){
                  nby.model <- naiveBayes(popNames_vector ~ ., data=trainSetData_PC)
                  nby.pred.class <- predict(nby.model,testSetData_PC,type="class")
                  nby.pred.prob <- predict(nby.model,testSetData_PC,type="raw")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(nby.pred.class), as.data.frame(nby.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="tree"){
                  tree.model <- tree(popNames_vector ~ ., data=trainSetData_PC)
                  tree.pred.class <- predict(tree.model,testSetData_PC,type="class")
                  tree.pred.prob <- predict(tree.model,testSetData_PC,type="vector")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(tree.pred.class), as.data.frame(tree.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  tree_node <- as.character(summary(tree.model)$used)#Get loci used at tree node
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  cat("Loci used at tree node",tree_node,file=paste0(dir,"Loci_treenode_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="randomForest"){
                  rf.model <- randomForest(popNames_vector ~ ., data=trainSetData_PC, ntree=ntree, importance=T)
                  rf.pred.class <- predict(rf.model,testSetData_PC,type="response")
                  rf.pred.prob <- predict(rf.model,testSetData_PC,type="prob")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(rf.pred.class), as.data.frame(rf.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  #Output "importance"(Gini index) of loci to files;see randomForest package's "importance" function
                  write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Loci_importance_",f,"_K",k,"_",i,".txt"), quote=F, row.names=T)
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }
              }# for(f in train.loci)
            }#else if(loci.sample="random")
          }# for(i in 1:k)...
        }#close for if(multiprocess){...}eslse{}  
      }# for(k in k.fold)
      
    }else if(length(x)==5){
      datatype <- "genetics + non-genetics"
      #Split geneitc and non-genetic data first for estimating Fst
      otherVarNames <- x[[4]]#Get non-genetic variable name
      noVars <- length(otherVarNames)
      firstOtherVars_pos <- length(genoMatrix) - noVars #column position of first non-geneitc var
      lastOtherVars_pos <- length(genoMatrix) - 1 #column position of last non-genetic var
      #loop each k fold
      for(k in k.fold){
        #create fold index based on label (genoMatrix$popNames_vector)
        fold_index <- createFolds(genoMatrix$popNames_vector, k=k) #this is a list containing k folds of index, e.g., fold_index$Fold1, $Fold2...
        #take the K-th fold as test set; remaing folds as training set
        if(multiprocess){
          if (processors <= maxCores & processors > 0){
            cl <- makeCluster(processors)
            registerDoParallel(cl,cores=processors)
          }else {
            cl <- makeCluster(maxCores)
            registerDoParallel(cl,cores=maxCores)
          }
          foreach(i=1:k, .export=c("Fsts","perform.PCA"), .packages=c("e1071","klaR","MASS","tree","randomForest")) %dopar% {
            trainSetMatrix <- genoMatrix[-fold_index[[i]],]
            testIndID <- x[[2]][fold_index[[i]]]
            trainSetMatrix_genetic <- trainSetMatrix[,-(firstOtherVars_pos:lastOtherVars_pos)] #test inds. with only genetic variables and popNames
            #check if fst prior for sampling loci
            if(loci.sample=="fst"){
              train_X <- list(trainSetMatrix_genetic, row.names(trainSetMatrix_genetic), x[[3]])#create a new x list for Fsts function; x[[3]] is locus name
              fstTable <- Fsts(train_X)#Estimate locus Fst for training data
              orderFst <- order(fstTable$Fst, decreasing=T)#Order the loci index by Fst value(from highest)
              for(f in train.loci){
                trainLocusIndex_fstTable <- orderFst[1:round(noLocus*f)]#Get training locus index (from fstTable)
                trainLocusName <- as.character(fstTable[trainLocusIndex_fstTable,]$Locus)#Get training locus name
                trainLocusIndex_genoMatrix <- NULL #create a train locus index for genoMatrix to be extracted
                for(m in 1:length(trainLocusName)){
                  tempAlleleIndex <- grep(pattern=paste0(trainLocusName[m],"_"), alleleName)#alleleName is colnames(genoMatrix)
                  trainLocusIndex_genoMatrix <- c(trainLocusIndex_genoMatrix,tempAlleleIndex )
                }
                trainLocusIndex_genoMatrix <- sort(trainLocusIndex_genoMatrix)
                genoMatrix_trainVars <- genoMatrix[,c(trainLocusIndex_genoMatrix,firstOtherVars_pos:lastOtherVars_pos)]#resahpe the genoMatrix data, comprising training loci + non-genetic
                #Scale and center the entire dataset if scaled=T
                if(scaled){
                  genoMatrix_trainVars <- as.data.frame(scale(genoMatrix_trainVars))
                }
                #Add popNames_vector to the last column
                genoMatrix_trainVars <- cbind(genoMatrix_trainVars, genoMatrix$popNames_vector);colnames(genoMatrix_trainVars)[ncol(genoMatrix_trainVars)] <- "popNames_vector"
                ##Split entire data into training and test sets
                trainSetData <- genoMatrix_trainVars[-fold_index[[i]],]#Get the training set data (training individuals/loci)
                testSetData <- genoMatrix_trainVars[fold_index[[i]],]#Get the test set data (test individuals/loci)
                ##
                
                #Determine PCA method for non-genetic data
                if(pca.method=="mixed"){
                  #Perform PCA on training (genetic+non-genetic)
                  PCA_results <- perform.PCA(trainSetData[,1:ncol(trainSetData)-1], method=pca.PCs) #Run PCA without label column
                  loadings <- PCA_results[[1]] #loadings (coefficience) of variables and PCs; apply this to test data
                  trainSetData_PC <- as.data.frame(PCA_results[[2]])
                  trainSetData_PC <- cbind(trainSetData_PC, trainSetData$popNames_vector) ##Will be used for building predictive models
                  colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popNames_vector"
                  #Convert test data to PC variables based on training's loadings
                  testSetData_matrix <- as.matrix(testSetData[,1:ncol(testSetData)-1])
                  testSetData_PC <- as.data.frame(testSetData_matrix %*% loadings)
                  testSetData_PC <- cbind(testSetData_PC, testSetData$popNames_vector)
                  colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popNames_vector"
                }else if(pca.method=="independent"){
                  #Perform PCA on training genetic and non-genetic independently
                  PCA_result_genetics <- perform.PCA(trainSetData[,1:length(trainLocusIndex_genoMatrix)], method=pca.PCs)#Perform PCA on genetic data
                  PCA_result_nongenetics <- perform.PCA(trainSetData[,(length(trainLocusIndex_genoMatrix)+1):(ncol(trainSetData)-1)], method=pca.PCs)#perform PCA on non-genetic data
                  loadings_genetics <- PCA_result_genetics[[1]] #loadings of genetic PCs; apply this to test data
                  loadings_nongenetics <- PCA_result_nongenetics[[1]] #loadings of non-genetic PCs; apply this to test data
                  trainSetData_genetic_PC <- as.data.frame(PCA_result_genetics[[2]]);colnames(trainSetData_genetic_PC) <- sub("PC", "genPC", colnames(trainSetData_genetic_PC))
                  trainSetData_nongenetic_PC <- as.data.frame(PCA_result_nongenetics[[2]]);colnames(trainSetData_nongenetic_PC) <- sub("PC", "nPC", colnames(trainSetData_nongenetic_PC))
                  trainSetData_PC <- cbind(trainSetData_genetic_PC, trainSetData_nongenetic_PC, trainSetData$popNames_vector)#Will use for building predictive models 
                  colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popNames_vector"
                  #Convert test data to PC variables based on training's loadings
                  testSetData_genetic_matrix <- as.matrix(testSetData[,1:length(trainLocusIndex_genoMatrix)]) #make genetic data matrix
                  testSetData_genetic_PC <- as.data.frame(testSetData_genetic_matrix %*% loadings_genetics);colnames(testSetData_genetic_PC)<-sub("PC", "genPC", colnames(testSetData_genetic_PC))  
                  testSetData_nongenetic_matrix <- as.matrix(testSetData[,(length(trainLocusIndex_genoMatrix)+1):(ncol(testSetData)-1)])
                  testSetData_nongenetic_PC <- as.data.frame(testSetData_nongenetic_matrix %*% loadings_nongenetics);colnames(testSetData_nongenetic_PC) <- sub("PC","nPC",colnames(testSetData_nongenetic_PC))
                  testSetData_PC <- cbind(testSetData_genetic_PC, testSetData_nongenetic_PC, testSetData$popNames_vector)
                  colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popNames_vector"
                }else if(pca.method=="original"){
                  #Perform PCA on only genetic data
                  PCA_result_genetics <- perform.PCA(trainSetData[,1:length(trainLocusIndex_genoMatrix)], method=pca.PCs)#Perform PCA on genetic data
                  loadings_genetics <- PCA_result_genetics[[1]] #loadings of genetic PCs; apply this to test data
                  trainSetData_genetic_PC <- as.data.frame(PCA_result_genetics[[2]]);colnames(trainSetData_genetic_PC) <- sub("PC", "genPC", colnames(trainSetData_genetic_PC))
                  #concatenate genetic PCs and original non-genetic data and popNames_vector; note that dataset has PCs only from genetics
                  trainSetData_PC <- cbind(trainSetData_genetic_PC, trainSetData[,(length(trainLocusIndex_genoMatrix)+1):(ncol(trainSetData)-1)], trainSetData$popNames_vector)
                  colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popNames_vector"
                  #Convert test data (genetic part) to PC variables (scores) based on training 
                  testSetData_genetic_matrix <- as.matrix(testSetData[,1:length(trainLocusIndex_genoMatrix)])
                  testSetData_genetic_PC <- as.data.frame(testSetData_genetic_matrix %*% loadings_genetics);colnames(testSetData_genetic_PC)<-sub("PC", "genPC", colnames(testSetData_genetic_PC)) 
                  testSetData_PC <- cbind(testSetData_genetic_PC, testSetData[,(length(trainLocusIndex_genoMatrix)+1):(ncol(trainSetData)-1)], testSetData$popNames_vector)
                  colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popNames_vector"
                }
                
                #Use training to build models and test on test individuals
                if(model=="svm"){
                  svm.fit <- svm(popNames_vector ~ ., data=trainSetData_PC, kernel=svm.kernel, cost=svm.cost, prob=T, scale=F)
                  svm.pred <- predict(svm.fit, testSetData_PC, type="class",prob=T)
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(svm.pred), attr(svm.pred,"probabilities"))#combine output to data frame
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="knn"){
                  
                }else if (model=="lda"){
                  lda.fit <- lda(popNames_vector ~ ., data=trainSetData_PC)
                  lda.pred <- predict(lda.fit, testSetData_PC)
                  lda.pred.class <- lda.pred$class
                  lda.pred.prob <- lda.pred$posterior
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(lda.pred.class), as.data.frame(lda.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="naiveBayes"){
                  nby.model <- naiveBayes(popNames_vector ~ ., data=trainSetData_PC)
                  nby.pred.class <- predict(nby.model,testSetData_PC,type="class")
                  nby.pred.prob <- predict(nby.model,testSetData_PC,type="raw")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(nby.pred.class), as.data.frame(nby.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="tree"){
                  tree.model <- tree(popNames_vector ~ ., data=trainSetData_PC)
                  tree.pred.class <- predict(tree.model,testSetData_PC,type="class")
                  tree.pred.prob <- predict(tree.model,testSetData_PC,type="vector")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(tree.pred.class), as.data.frame(tree.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  tree_node <- as.character(summary(tree.model)$used)#Get loci used at tree node
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  cat("Loci used at tree node",tree_node,file=paste0(dir,"Loci_treenode_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="randomForest"){
                  rf.model <- randomForest(popNames_vector ~ ., data=trainSetData_PC, ntree=ntree, importance=T)
                  rf.pred.class <- predict(rf.model,testSetData_PC,type="response")
                  rf.pred.prob <- predict(rf.model,testSetData_PC,type="prob")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(rf.pred.class), as.data.frame(rf.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  #Output "importance"(Gini index) of loci to files;see randomForest package's "importance" function
                  write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Loci_importance_",f,"_K",k,"_",i,".txt"), quote=F, row.names=T)
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }
              }#for(f in train.loci)
            }else if(loci.sample=="random"){
              for(f in train.loci){
                tempLocusIndex <- sort(sample(1:noLocus, round(noLocus*f)))
                trainLocusName <- locusNames[tempLocusIndex]
                trainLocusIndex_genoMatrix <- NULL #create a train locus index for genoMatrix to be extracted
                for(m in 1:length(trainLocusName)){
                  tempAlleleIndex <- grep(pattern=paste0(trainLocusName[m],"_"), alleleName)
                  trainLocusIndex_genoMatrix <- c(trainLocusIndex_genoMatrix,tempAlleleIndex )
                }
                trainLocusIndex_genoMatrix <- sort(trainLocusIndex_genoMatrix)
                genoMatrix_trainVars <- genoMatrix[,c(trainLocusIndex_genoMatrix,firstOtherVars_pos:lastOtherVars_pos)]#resahpe the genoMatrix data, comprising training loci + non-genetic
                #Scale and center the entire dataset if scaled=T
                if(scaled){
                  genoMatrix_trainVars <- as.data.frame(scale(genoMatrix_trainVars))
                }
                #Add popNames_vector to the last column
                genoMatrix_trainVars <- cbind(genoMatrix_trainVars, genoMatrix$popNames_vector);colnames(genoMatrix_trainVars)[ncol(genoMatrix_trainVars)] <- "popNames_vector"
                ##Split entire data into training and test sets
                trainSetData <- genoMatrix_trainVars[-fold_index[[i]],]#Get the training set data (training individuals/loci)
                testSetData <- genoMatrix_trainVars[fold_index[[i]],]#Get the test set data (test individuals/loci)
                ##
                #Determine PCA method for non-genetic data
                if(pca.method=="mixed"){
                  #Perform PCA on training (genetic+non-genetic)
                  PCA_results <- perform.PCA(trainSetData[,1:ncol(trainSetData)-1], method=pca.PCs) #Run PCA without label column
                  loadings <- PCA_results[[1]] #loadings (coefficience) of variables and PCs; apply this to test data
                  trainSetData_PC <- as.data.frame(PCA_results[[2]])
                  trainSetData_PC <- cbind(trainSetData_PC, trainSetData$popNames_vector) ##Will be used for building predictive models
                  colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popNames_vector"
                  #Convert test data to PC variables based on training's loadings
                  testSetData_matrix <- as.matrix(testSetData[,1:ncol(testSetData)-1])
                  testSetData_PC <- as.data.frame(testSetData_matrix %*% loadings)
                  testSetData_PC <- cbind(testSetData_PC, testSetData$popNames_vector)
                  colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popNames_vector"
                }else if(pca.method=="independent"){
                  #Perform PCA on training genetic and non-genetic independently
                  PCA_result_genetics <- perform.PCA(trainSetData[,1:length(trainLocusIndex_genoMatrix)], method=pca.PCs)#Perform PCA on genetic data
                  PCA_result_nongenetics <- perform.PCA(trainSetData[,(length(trainLocusIndex_genoMatrix)+1):(ncol(trainSetData)-1)], method=pca.PCs)#perform PCA on non-genetic data
                  loadings_genetics <- PCA_result_genetics[[1]] #loadings of genetic PCs; apply this to test data
                  loadings_nongenetics <- PCA_result_nongenetics[[1]] #loadings of non-genetic PCs; apply this to test data
                  trainSetData_genetic_PC <- as.data.frame(PCA_result_genetics[[2]]);colnames(trainSetData_genetic_PC) <- sub("PC", "genPC", colnames(trainSetData_genetic_PC))
                  trainSetData_nongenetic_PC <- as.data.frame(PCA_result_nongenetics[[2]]);colnames(trainSetData_nongenetic_PC) <- sub("PC", "nPC", colnames(trainSetData_nongenetic_PC))
                  trainSetData_PC <- cbind(trainSetData_genetic_PC, trainSetData_nongenetic_PC, trainSetData$popNames_vector)#Will use for building predictive models 
                  colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popNames_vector"
                  #Convert test data to PC variables based on training's loadings
                  testSetData_genetic_matrix <- as.matrix(testSetData[,1:length(trainLocusIndex_genoMatrix)]) #make genetic data matrix
                  testSetData_genetic_PC <- as.data.frame(testSetData_genetic_matrix %*% loadings_genetics);colnames(testSetData_genetic_PC)<-sub("PC", "genPC", colnames(testSetData_genetic_PC))  
                  testSetData_nongenetic_matrix <- as.matrix(testSetData[,(length(trainLocusIndex_genoMatrix)+1):(ncol(testSetData)-1)])
                  testSetData_nongenetic_PC <- as.data.frame(testSetData_nongenetic_matrix %*% loadings_nongenetics);colnames(testSetData_nongenetic_PC) <- sub("PC","nPC",colnames(testSetData_nongenetic_PC))
                  testSetData_PC <- cbind(testSetData_genetic_PC, testSetData_nongenetic_PC, testSetData$popNames_vector)
                  colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popNames_vector"
                }else if(pca.method=="original"){
                  #Perform PCA on only genetic data
                  PCA_result_genetics <- perform.PCA(trainSetData[,1:length(trainLocusIndex_genoMatrix)], method=pca.PCs)#Perform PCA on genetic data
                  loadings_genetics <- PCA_result_genetics[[1]] #loadings of genetic PCs; apply this to test data
                  trainSetData_genetic_PC <- as.data.frame(PCA_result_genetics[[2]]);colnames(trainSetData_genetic_PC) <- sub("PC", "genPC", colnames(trainSetData_genetic_PC))
                  #concatenate genetic PCs and original non-genetic data and popNames_vector; note that dataset has PCs only from genetics
                  trainSetData_PC <- cbind(trainSetData_genetic_PC, trainSetData[,(length(trainLocusIndex_genoMatrix)+1):(ncol(trainSetData)-1)], trainSetData$popNames_vector)
                  colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popNames_vector"
                  #Convert test data (genetic part) to PC variables (scores) based on training 
                  testSetData_genetic_matrix <- as.matrix(testSetData[,1:length(trainLocusIndex_genoMatrix)])
                  testSetData_genetic_PC <- as.data.frame(testSetData_genetic_matrix %*% loadings_genetics);colnames(testSetData_genetic_PC)<-sub("PC", "genPC", colnames(testSetData_genetic_PC)) 
                  testSetData_PC <- cbind(testSetData_genetic_PC, testSetData[,(length(trainLocusIndex_genoMatrix)+1):(ncol(trainSetData)-1)], testSetData$popNames_vector)
                  colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popNames_vector"
                }
                
                #Use training to build models and test on test individuals
                if(model=="svm"){
                  svm.fit <- svm(popNames_vector ~ ., data=trainSetData_PC, kernel=svm.kernel, cost=svm.cost, prob=T, scale=F)
                  svm.pred <- predict(svm.fit, testSetData_PC, type="class",prob=T)
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(svm.pred), attr(svm.pred,"probabilities"))#combine output to data frame
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="knn"){
                  
                }else if (model=="lda"){
                  lda.fit <- lda(popNames_vector ~ ., data=trainSetData_PC)
                  lda.pred <- predict(lda.fit, testSetData_PC)
                  lda.pred.class <- lda.pred$class
                  lda.pred.prob <- lda.pred$posterior
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(lda.pred.class), as.data.frame(lda.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="naiveBayes"){
                  nby.model <- naiveBayes(popNames_vector ~ ., data=trainSetData_PC)
                  nby.pred.class <- predict(nby.model,testSetData_PC,type="class")
                  nby.pred.prob <- predict(nby.model,testSetData_PC,type="raw")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(nby.pred.class), as.data.frame(nby.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="tree"){
                  tree.model <- tree(popNames_vector ~ ., data=trainSetData_PC)
                  tree.pred.class <- predict(tree.model,testSetData_PC,type="class")
                  tree.pred.prob <- predict(tree.model,testSetData_PC,type="vector")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(tree.pred.class), as.data.frame(tree.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  tree_node <- as.character(summary(tree.model)$used)#Get loci used at tree node
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  cat("Loci used at tree node",tree_node,file=paste0(dir,"Loci_treenode_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="randomForest"){
                  rf.model <- randomForest(popNames_vector ~ ., data=trainSetData_PC, ntree=ntree, importance=T)
                  rf.pred.class <- predict(rf.model,testSetData_PC,type="response")
                  rf.pred.prob <- predict(rf.model,testSetData_PC,type="prob")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(rf.pred.class), as.data.frame(rf.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  #Output "importance"(Gini index) of loci to files;see randomForest package's "importance" function
                  write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Loci_importance_",f,"_K",k,"_",i,".txt"), quote=F, row.names=T)
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }
              }#for(f in train.loci)
            }#else if(loci.sample==...
          }#foreach(i=1:k,...
          stopCluster(cl)
        }else{
          for(i in 1:k){
            trainSetMatrix <- genoMatrix[-fold_index[[i]],]
            testIndID <- x[[2]][fold_index[[i]]]
            trainSetMatrix_genetic <- trainSetMatrix[,-(firstOtherVars_pos:lastOtherVars_pos)] #test inds. with only genetic variables and popNames
            #check if fst prior for sampling loci
            if(loci.sample=="fst"){
              train_X <- list(trainSetMatrix_genetic, row.names(trainSetMatrix_genetic), x[[3]])#create a new x list for Fsts function; x[[3]] is locus name
              fstTable <- Fsts(train_X)#Estimate locus Fst for training data
              orderFst <- order(fstTable$Fst, decreasing=T)#Order the loci index by Fst value(from highest)
              for(f in train.loci){
                trainLocusIndex_fstTable <- orderFst[1:round(noLocus*f)]#Get training locus index (from fstTable)
                trainLocusName <- as.character(fstTable[trainLocusIndex_fstTable,]$Locus)#Get training locus name
                trainLocusIndex_genoMatrix <- NULL #create a train locus index for genoMatrix to be extracted
                for(m in 1:length(trainLocusName)){
                  tempAlleleIndex <- grep(pattern=paste0(trainLocusName[m],"_"), alleleName)#alleleName is colnames(genoMatrix)
                  trainLocusIndex_genoMatrix <- c(trainLocusIndex_genoMatrix,tempAlleleIndex )
                }
                trainLocusIndex_genoMatrix <- sort(trainLocusIndex_genoMatrix)
                genoMatrix_trainVars <- genoMatrix[,c(trainLocusIndex_genoMatrix,firstOtherVars_pos:lastOtherVars_pos)]#resahpe the genoMatrix data, comprising training loci + non-genetic
                #Scale and center the entire dataset if scaled=T
                if(scaled){
                  genoMatrix_trainVars <- as.data.frame(scale(genoMatrix_trainVars))
                }
                #Add popNames_vector to the last column
                genoMatrix_trainVars <- cbind(genoMatrix_trainVars, genoMatrix$popNames_vector);colnames(genoMatrix_trainVars)[ncol(genoMatrix_trainVars)] <- "popNames_vector"
                ##Split entire data into training and test sets
                trainSetData <- genoMatrix_trainVars[-fold_index[[i]],]#Get the training set data (training individuals/loci)
                testSetData <- genoMatrix_trainVars[fold_index[[i]],]#Get the test set data (test individuals/loci)
                ##
                
                #Determine PCA method for non-genetic data
                if(pca.method=="mixed"){
                  #Perform PCA on training (genetic+non-genetic)
                  PCA_results <- perform.PCA(trainSetData[,1:ncol(trainSetData)-1], method=pca.PCs) #Run PCA without label column
                  loadings <- PCA_results[[1]] #loadings (coefficience) of variables and PCs; apply this to test data
                  trainSetData_PC <- as.data.frame(PCA_results[[2]])
                  trainSetData_PC <- cbind(trainSetData_PC, trainSetData$popNames_vector) ##Will be used for building predictive models
                  colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popNames_vector"
                  #Convert test data to PC variables based on training's loadings
                  testSetData_matrix <- as.matrix(testSetData[,1:ncol(testSetData)-1])
                  testSetData_PC <- as.data.frame(testSetData_matrix %*% loadings)
                  testSetData_PC <- cbind(testSetData_PC, testSetData$popNames_vector)
                  colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popNames_vector"
                }else if(pca.method=="independent"){
                  #Perform PCA on training genetic and non-genetic independently
                  PCA_result_genetics <- perform.PCA(trainSetData[,1:length(trainLocusIndex_genoMatrix)], method=pca.PCs)#Perform PCA on genetic data
                  PCA_result_nongenetics <- perform.PCA(trainSetData[,(length(trainLocusIndex_genoMatrix)+1):(ncol(trainSetData)-1)], method=pca.PCs)#perform PCA on non-genetic data
                  loadings_genetics <- PCA_result_genetics[[1]] #loadings of genetic PCs; apply this to test data
                  loadings_nongenetics <- PCA_result_nongenetics[[1]] #loadings of non-genetic PCs; apply this to test data
                  trainSetData_genetic_PC <- as.data.frame(PCA_result_genetics[[2]]);colnames(trainSetData_genetic_PC) <- sub("PC", "genPC", colnames(trainSetData_genetic_PC))
                  trainSetData_nongenetic_PC <- as.data.frame(PCA_result_nongenetics[[2]]);colnames(trainSetData_nongenetic_PC) <- sub("PC", "nPC", colnames(trainSetData_nongenetic_PC))
                  trainSetData_PC <- cbind(trainSetData_genetic_PC, trainSetData_nongenetic_PC, trainSetData$popNames_vector)#Will use for building predictive models 
                  colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popNames_vector"
                  #Convert test data to PC variables based on training's loadings
                  testSetData_genetic_matrix <- as.matrix(testSetData[,1:length(trainLocusIndex_genoMatrix)]) #make genetic data matrix
                  testSetData_genetic_PC <- as.data.frame(testSetData_genetic_matrix %*% loadings_genetics);colnames(testSetData_genetic_PC)<-sub("PC", "genPC", colnames(testSetData_genetic_PC))  
                  testSetData_nongenetic_matrix <- as.matrix(testSetData[,(length(trainLocusIndex_genoMatrix)+1):(ncol(testSetData)-1)])
                  testSetData_nongenetic_PC <- as.data.frame(testSetData_nongenetic_matrix %*% loadings_nongenetics);colnames(testSetData_nongenetic_PC) <- sub("PC","nPC",colnames(testSetData_nongenetic_PC))
                  testSetData_PC <- cbind(testSetData_genetic_PC, testSetData_nongenetic_PC, testSetData$popNames_vector)
                  colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popNames_vector"
                }else if(pca.method=="original"){
                  #Perform PCA on only genetic data
                  PCA_result_genetics <- perform.PCA(trainSetData[,1:length(trainLocusIndex_genoMatrix)], method=pca.PCs)#Perform PCA on genetic data
                  loadings_genetics <- PCA_result_genetics[[1]] #loadings of genetic PCs; apply this to test data
                  trainSetData_genetic_PC <- as.data.frame(PCA_result_genetics[[2]]);colnames(trainSetData_genetic_PC) <- sub("PC", "genPC", colnames(trainSetData_genetic_PC))
                  #concatenate genetic PCs and original non-genetic data and popNames_vector; note that dataset has PCs only from genetics
                  trainSetData_PC <- cbind(trainSetData_genetic_PC, trainSetData[,(length(trainLocusIndex_genoMatrix)+1):(ncol(trainSetData)-1)], trainSetData$popNames_vector)
                  colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popNames_vector"
                  #Convert test data (genetic part) to PC variables (scores) based on training 
                  testSetData_genetic_matrix <- as.matrix(testSetData[,1:length(trainLocusIndex_genoMatrix)])
                  testSetData_genetic_PC <- as.data.frame(testSetData_genetic_matrix %*% loadings_genetics);colnames(testSetData_genetic_PC)<-sub("PC", "genPC", colnames(testSetData_genetic_PC)) 
                  testSetData_PC <- cbind(testSetData_genetic_PC, testSetData[,(length(trainLocusIndex_genoMatrix)+1):(ncol(trainSetData)-1)], testSetData$popNames_vector)
                  colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popNames_vector"
                }
                
                #Use training to build models and test on test individuals
                if(model=="svm"){
                  svm.fit <- svm(popNames_vector ~ ., data=trainSetData_PC, kernel=svm.kernel, cost=svm.cost, prob=T, scale=F)
                  svm.pred <- predict(svm.fit, testSetData_PC, type="class",prob=T)
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(svm.pred), attr(svm.pred,"probabilities"))#combine output to data frame
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="knn"){
                  
                }else if (model=="lda"){
                  lda.fit <- lda(popNames_vector ~ ., data=trainSetData_PC)
                  lda.pred <- predict(lda.fit, testSetData_PC)
                  lda.pred.class <- lda.pred$class
                  lda.pred.prob <- lda.pred$posterior
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(lda.pred.class), as.data.frame(lda.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="naiveBayes"){
                  nby.model <- naiveBayes(popNames_vector ~ ., data=trainSetData_PC)
                  nby.pred.class <- predict(nby.model,testSetData_PC,type="class")
                  nby.pred.prob <- predict(nby.model,testSetData_PC,type="raw")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(nby.pred.class), as.data.frame(nby.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="tree"){
                  tree.model <- tree(popNames_vector ~ ., data=trainSetData_PC)
                  tree.pred.class <- predict(tree.model,testSetData_PC,type="class")
                  tree.pred.prob <- predict(tree.model,testSetData_PC,type="vector")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(tree.pred.class), as.data.frame(tree.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  tree_node <- as.character(summary(tree.model)$used)#Get loci used at tree node
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  cat("Loci used at tree node",tree_node,file=paste0(dir,"Loci_treenode_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="randomForest"){
                  rf.model <- randomForest(popNames_vector ~ ., data=trainSetData_PC, ntree=ntree, importance=T)
                  rf.pred.class <- predict(rf.model,testSetData_PC,type="response")
                  rf.pred.prob <- predict(rf.model,testSetData_PC,type="prob")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(rf.pred.class), as.data.frame(rf.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  #Output "importance"(Gini index) of loci to files;see randomForest package's "importance" function
                  write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Loci_importance_",f,"_K",k,"_",i,".txt"), quote=F, row.names=T)
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }
              }#for(f in train.loci)
            }else if(loci.sample=="random"){
              for(f in train.loci){
                tempLocusIndex <- sort(sample(1:noLocus, round(noLocus*f)))
                trainLocusName <- locusNames[tempLocusIndex]
                trainLocusIndex_genoMatrix <- NULL #create a train locus index for genoMatrix to be extracted
                for(m in 1:length(trainLocusName)){
                  tempAlleleIndex <- grep(pattern=paste0(trainLocusName[m],"_"), alleleName)
                  trainLocusIndex_genoMatrix <- c(trainLocusIndex_genoMatrix,tempAlleleIndex )
                }
                trainLocusIndex_genoMatrix <- sort(trainLocusIndex_genoMatrix)
                genoMatrix_trainVars <- genoMatrix[,c(trainLocusIndex_genoMatrix,firstOtherVars_pos:lastOtherVars_pos)]#resahpe the genoMatrix data, comprising training loci + non-genetic
                #Scale and center the entire dataset if scaled=T
                if(scaled){
                  genoMatrix_trainVars <- as.data.frame(scale(genoMatrix_trainVars))
                }
                #Add popNames_vector to the last column
                genoMatrix_trainVars <- cbind(genoMatrix_trainVars, genoMatrix$popNames_vector);colnames(genoMatrix_trainVars)[ncol(genoMatrix_trainVars)] <- "popNames_vector"
                ##Split entire data into training and test sets
                trainSetData <- genoMatrix_trainVars[-fold_index[[i]],]#Get the training set data (training individuals/loci)
                testSetData <- genoMatrix_trainVars[fold_index[[i]],]#Get the test set data (test individuals/loci)
                ##
                #Determine PCA method for non-genetic data
                if(pca.method=="mixed"){
                  #Perform PCA on training (genetic+non-genetic)
                  PCA_results <- perform.PCA(trainSetData[,1:ncol(trainSetData)-1], method=pca.PCs) #Run PCA without label column
                  loadings <- PCA_results[[1]] #loadings (coefficience) of variables and PCs; apply this to test data
                  trainSetData_PC <- as.data.frame(PCA_results[[2]])
                  trainSetData_PC <- cbind(trainSetData_PC, trainSetData$popNames_vector) ##Will be used for building predictive models
                  colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popNames_vector"
                  #Convert test data to PC variables based on training's loadings
                  testSetData_matrix <- as.matrix(testSetData[,1:ncol(testSetData)-1])
                  testSetData_PC <- as.data.frame(testSetData_matrix %*% loadings)
                  testSetData_PC <- cbind(testSetData_PC, testSetData$popNames_vector)
                  colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popNames_vector"
                }else if(pca.method=="independent"){
                  #Perform PCA on training genetic and non-genetic independently
                  PCA_result_genetics <- perform.PCA(trainSetData[,1:length(trainLocusIndex_genoMatrix)], method=pca.PCs)#Perform PCA on genetic data
                  PCA_result_nongenetics <- perform.PCA(trainSetData[,(length(trainLocusIndex_genoMatrix)+1):(ncol(trainSetData)-1)], method=pca.PCs)#perform PCA on non-genetic data
                  loadings_genetics <- PCA_result_genetics[[1]] #loadings of genetic PCs; apply this to test data
                  loadings_nongenetics <- PCA_result_nongenetics[[1]] #loadings of non-genetic PCs; apply this to test data
                  trainSetData_genetic_PC <- as.data.frame(PCA_result_genetics[[2]]);colnames(trainSetData_genetic_PC) <- sub("PC", "genPC", colnames(trainSetData_genetic_PC))
                  trainSetData_nongenetic_PC <- as.data.frame(PCA_result_nongenetics[[2]]);colnames(trainSetData_nongenetic_PC) <- sub("PC", "nPC", colnames(trainSetData_nongenetic_PC))
                  trainSetData_PC <- cbind(trainSetData_genetic_PC, trainSetData_nongenetic_PC, trainSetData$popNames_vector)#Will use for building predictive models 
                  colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popNames_vector"
                  #Convert test data to PC variables based on training's loadings
                  testSetData_genetic_matrix <- as.matrix(testSetData[,1:length(trainLocusIndex_genoMatrix)]) #make genetic data matrix
                  testSetData_genetic_PC <- as.data.frame(testSetData_genetic_matrix %*% loadings_genetics);colnames(testSetData_genetic_PC)<-sub("PC", "genPC", colnames(testSetData_genetic_PC))  
                  testSetData_nongenetic_matrix <- as.matrix(testSetData[,(length(trainLocusIndex_genoMatrix)+1):(ncol(testSetData)-1)])
                  testSetData_nongenetic_PC <- as.data.frame(testSetData_nongenetic_matrix %*% loadings_nongenetics);colnames(testSetData_nongenetic_PC) <- sub("PC","nPC",colnames(testSetData_nongenetic_PC))
                  testSetData_PC <- cbind(testSetData_genetic_PC, testSetData_nongenetic_PC, testSetData$popNames_vector)
                  colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popNames_vector"
                }else if(pca.method=="original"){
                  #Perform PCA on only genetic data
                  PCA_result_genetics <- perform.PCA(trainSetData[,1:length(trainLocusIndex_genoMatrix)], method=pca.PCs)#Perform PCA on genetic data
                  loadings_genetics <- PCA_result_genetics[[1]] #loadings of genetic PCs; apply this to test data
                  trainSetData_genetic_PC <- as.data.frame(PCA_result_genetics[[2]]);colnames(trainSetData_genetic_PC) <- sub("PC", "genPC", colnames(trainSetData_genetic_PC))
                  #concatenate genetic PCs and original non-genetic data and popNames_vector; note that dataset has PCs only from genetics
                  trainSetData_PC <- cbind(trainSetData_genetic_PC, trainSetData[,(length(trainLocusIndex_genoMatrix)+1):(ncol(trainSetData)-1)], trainSetData$popNames_vector)
                  colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popNames_vector"
                  #Convert test data (genetic part) to PC variables (scores) based on training 
                  testSetData_genetic_matrix <- as.matrix(testSetData[,1:length(trainLocusIndex_genoMatrix)])
                  testSetData_genetic_PC <- as.data.frame(testSetData_genetic_matrix %*% loadings_genetics);colnames(testSetData_genetic_PC)<-sub("PC", "genPC", colnames(testSetData_genetic_PC)) 
                  testSetData_PC <- cbind(testSetData_genetic_PC, testSetData[,(length(trainLocusIndex_genoMatrix)+1):(ncol(trainSetData)-1)], testSetData$popNames_vector)
                  colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popNames_vector"
                }
                
                #Use training to build models and test on test individuals
                if(model=="svm"){
                  svm.fit <- svm(popNames_vector ~ ., data=trainSetData_PC, kernel=svm.kernel, cost=svm.cost, prob=T, scale=F)
                  svm.pred <- predict(svm.fit, testSetData_PC, type="class",prob=T)
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(svm.pred), attr(svm.pred,"probabilities"))#combine output to data frame
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="knn"){
                  
                }else if (model=="lda"){
                  lda.fit <- lda(popNames_vector ~ ., data=trainSetData_PC)
                  lda.pred <- predict(lda.fit, testSetData_PC)
                  lda.pred.class <- lda.pred$class
                  lda.pred.prob <- lda.pred$posterior
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(lda.pred.class), as.data.frame(lda.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="naiveBayes"){
                  nby.model <- naiveBayes(popNames_vector ~ ., data=trainSetData_PC)
                  nby.pred.class <- predict(nby.model,testSetData_PC,type="class")
                  nby.pred.prob <- predict(nby.model,testSetData_PC,type="raw")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(nby.pred.class), as.data.frame(nby.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="tree"){
                  tree.model <- tree(popNames_vector ~ ., data=trainSetData_PC)
                  tree.pred.class <- predict(tree.model,testSetData_PC,type="class")
                  tree.pred.prob <- predict(tree.model,testSetData_PC,type="vector")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(tree.pred.class), as.data.frame(tree.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  tree_node <- as.character(summary(tree.model)$used)#Get loci used at tree node
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  cat("Loci used at tree node",tree_node,file=paste0(dir,"Loci_treenode_",f,"_K",k,"_",i,".txt"), sep="\n")
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }else if(model=="randomForest"){
                  rf.model <- randomForest(popNames_vector ~ ., data=trainSetData_PC, ntree=ntree, importance=T)
                  rf.pred.class <- predict(rf.model,testSetData_PC,type="response")
                  rf.pred.prob <- predict(rf.model,testSetData_PC,type="prob")
                  outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(rf.pred.class), as.data.frame(rf.pred.prob))
                  colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                  ##Output assignment results (prob.& used train loci) to files
                  write.table(outcome_matrix, file=paste0(dir,"Out_",f,"_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
                  cat(trainLocusName, file=paste0(dir,"Loci_",f,"_K",k,"_",i,".txt"), sep="\n")
                  #Output "importance"(Gini index) of loci to files;see randomForest package's "importance" function
                  write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Loci_importance_",f,"_K",k,"_",i,".txt"), quote=F, row.names=T)
                  if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",f,"_K",k,"_",i,".txt"), quote=F)}
                  #
                }
              }#for(f in train.loci)
            }#else if(loci.sample==...
          }#for(i in 1:k)...
        }#close for if(multiprocess){}else{}
      }#for(k in k.fold)
    }#else if(length(x)==5)
    
    #Output a metadata file
    version <- as.character(packageVersion("assignPOP"))
    cat(" Analysis Description ( R - assignPOP ver.",version,")\n",
        "Perform assign.kfold() @", format(Sys.time()),"\n\n",
        "k.fold =",k.fold,"\n",
        "train.loci =",train.loci,"\n",
        "Total assignment tests =",sum(k.fold*length(train.loci)),"\n",
        "Fst locus sample method:",loci.sample,"\n",
        "Data scaled and centerd:",scaled,"\n",
        "PC retaining criteria:",pca.PCs,"\n",
        "PCA for non-genetic data:",pca.method,"\n",
        "Machine learning model:",model,"\n\n",
        "Input Data (",datatype,")\n",
        "Number of individuals:",sum(popSizes),"\n",
        "Number of loci:",noLocus,"\n",
        "Number of non-genetic variables:",noVars,"\n",
        "Number of populations:",noPops,"\n",
        names(popSizes),"\n",popSizes,"\n",
        file=paste0(dir,"AnalysisInfo.txt"))
    #Print some message to R console
    cat("\n  K-fold cross-validation done!!")
    cat(paste0("\n  ",sum(k.fold*length(train.loci))," assignment tests completed!!"))
    
  }else if(is.data.frame(x)){
    #Analyze non-genetic data
    #checking pca.method
    if(is.character(pca.method)){
      if(skipQ){
        anspca <- "N"
      }else{
        anspca <- readline("  Perform PCA on dataset for dimensionality reduction? (enter Y/N): ")
      }
      if(grepl(pattern="Y", toupper(anspca))){
        pca.method <- TRUE
      }else if(grepl(pattern="N", toupper(anspca))){
        pca.method <- FALSE
      }else{
        stop("Please enter either Yes or No.")
      }
    }
    popName <- NULL; i <- NULL
    IndID <- x[ ,1] #get individual ID
    colnames(x)[ncol(x)] <- "popName"
    dataMatrix <- x[-1] #exclude first ID column  
    popSizes <- table(x[ncol(x)])
    pops <- names(popSizes)
    noPops <- length(popSizes)
    noVars <- ncol(x)-2
    varNames <- names(dataMatrix[1:ncol(dataMatrix)-1])#get variable name
    noFactorVar <- 0
    ##Check if pop size and k fold value can work
    if(max(k.fold) > min(popSizes)){
      stop("Max. K value is greater than a small pop. Please adjust your k.fold setting or increase sample size.")
    }
    #Check data type (numeric or categorical)
    for(n in 1:noVars){
      var_name <- varNames[n]
      var_type <- class(dataMatrix[,n])
      cat(paste0("  ",var_name,"(",var_type,")"))
    }
    if(skipQ){
      ans0 <- "Y"
    }else{
      ans0 <- readline("  Are they correct? (enter Y/N): ")
    }
    if(grepl(pattern="N", toupper(ans0))){
      cat("  please enter variable names for changing data type (separate names by a whitespace if multiple)\n")
      ans1 <- readline("  enter here: ")
      ans1 <- str_trim(ans1, side="both")
      ans1 <- unlist(strsplit(ans1,split=" "))#check out variable name to be processed
      noChangeVar <- length(ans1)
      #Check if entry is correct
      if(!all(ans1 %in% varNames)){ #if any of entry not in varNames is true 
        stop("Please enter correct feature names.")
      }
      #Process variables and convert factor data to dummy variable (binary data)
      for(name in ans1){
        ans2 <- readline(paste0("  Which data type should '",name,"' be? (enter numeric or factor): "))
        if(grepl(pattern="N",toupper(ans2))){
          dataMatrix[,name] <- as.numeric(as.character(dataMatrix[,name]))
        }else if(grepl(pattern="F",toupper(ans2))){
          dataMatrix[,name] <- as.factor(dataMatrix[,name])
        }
      }
      #Convert factor data to dummy data
      for(name in varNames){
        if(is.factor(dataMatrix[,name])){
          noFactorVar <- noFactorVar + 1 #count number of categorical varibales
          #Convert factor variable to numeric binary variable (dummy variable)
          dummyData <- as.data.frame(model.matrix( ~ dataMatrix[,name]-1, data=dataMatrix))#get dummy variable data frame
          names(dummyData) <- substring(names(dummyData), 19, 1000L)#extract meaningful wording, or remove some funny wording
          names(dummyData) <- sub("\\b", paste0(name,"."), names(dummyData))#append original variabel name at the beginning
          dataMatrix[,name] <- NULL #remove original factor data column
          dataMatrix <- cbind(dummyData, dataMatrix) #column bind dummy data while making popName last column
        }
      }
    }else if(grepl(pattern="Y",toupper(ans0))){
      #check through data and covert factor to dummy
      for(name in varNames){
        if(is.factor(dataMatrix[,name])){
          noFactorVar <- noFactorVar + 1 #count number of categorical varibales
          #Convert factor variable to numeric binary variable (dummy variable)
          dummyData <- as.data.frame(model.matrix( ~ dataMatrix[,name]-1, data=dataMatrix))#get dummy variable data frame
          names(dummyData) <- substring(names(dummyData), 19, 1000L)#extract meaningful wording, or remove some funny wording
          names(dummyData) <- sub("\\b", paste0(name,"."), names(dummyData))#append original variabel name at the beginning
          dataMatrix[,name] <- NULL #remove original factor data column
          dataMatrix <- cbind(dummyData, dataMatrix) #column bind dummy data while making popName last column
        }
      }
    }
    #create a folder to save outfiles
    dir.create(file.path(dir))
    #Detect CPU core/thread numbers and use n-1 threads for parallel ananlysis
    if(multiprocess){
      maxCores <- detectCores()-1
      if (processors <= maxCores & processors > 0){
        #cl <- makeCluster(processors)
        #registerDoParallel(cl,cores=processors)
        cat("\n  Parallel computing is on. Analyzing data using",processors,"cores/threads of CPU...\n")
      }else {
        #cl <- makeCluster(maxCores)
        #registerDoParallel(cl,cores=maxCores)
        cat("\n  Parallel computing is on. Analyzing data using",maxCores,"cores/threads of CPU...\n")
      }
    }else {
      cat("\n  Parallele computing is off. Analyzing data using 1 CPU core...\n")
    }
    #Scale and center the entire features if set scaled=T
    if(scaled){
      cat("\n  Scaling and centering entire data set...")
      dataMatrix_scaled <- as.data.frame(scale(dataMatrix[,1:(ncol(dataMatrix)-1)]))
      dataMatrix <- cbind(dataMatrix_scaled, dataMatrix$popName);colnames(dataMatrix)[ncol(dataMatrix)] <- "popName"
    }
    #Loop each fold
    for(k in k.fold){
      #create fold index based on label (genoMatrix$popNames_vector)
      fold_index <- createFolds(dataMatrix$popName, k=k) #this is a list containing k folds of index, e.g., fold_index$Fold1, $Fold2...
      #take the K-th fold as test set; remaing folds as training set
      if(multiprocess){
        if (processors <= maxCores & processors > 0){
          cl <- makeCluster(processors)
          registerDoParallel(cl,cores=processors)
        }else {
          cl <- makeCluster(maxCores)
          registerDoParallel(cl,cores=maxCores)
        }
        foreach(i=1:k, .export="perform.PCA", .packages=c("e1071","klaR","MASS","tree","randomForest")) %dopar% {
          trainSetData <- dataMatrix[-fold_index[[i]],]
          testIndID <- IndID[fold_index[[i]]]; testIndID <- droplevels(testIndID)
          #Create test dataset
          testSetData <- dataMatrix[fold_index[[i]],]
          #
          #Perform PCA on training data (pca.method=T)
          if(pca.method){
            PCA_results <- perform.PCA(trainSetData[,1:ncol(trainSetData)-1], method=pca.PCs)#Run PCA on training without label
            loadings <- PCA_results[[1]] #loadings (coefficience) of variables and PCs; apply this to test data
            trainSetData_PC <- as.data.frame(PCA_results[[2]])
            trainSetData_PC <- cbind(trainSetData_PC, trainSetData$popName) ##Will be used for building predicting models
            colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popName"
            #Convert test data to PC variables based on training's loadings
            testSetData_matrix <- as.matrix(testSetData[,1:ncol(testSetData)-1])
            testSetData_PC <- as.data.frame(testSetData_matrix %*% loadings)
            testSetData_PC <- cbind(testSetData_PC, testSetData$popName)
            colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popName"
          }else if(pca.method==F){
            trainSetData_PC <- trainSetData
            testSetData_PC <- testSetData
          }
          #Use training to build models and test on test individuals
          if(model=="svm"){
            svm.fit <- svm(popName ~ ., data=trainSetData_PC, kernel=svm.kernel, cost=svm.cost, prob=T, scale=F)
            svm.pred <- predict(svm.fit, testSetData_PC, type="class",prob=T)
            outcome_matrix <- cbind(testIndID, testSetData_PC$popName, as.data.frame(svm.pred), attr(svm.pred,"probabilities"))#combine output to data frame
            colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
            ##Output assignment results (prob.& used train loci) to files
            write.table(outcome_matrix, file=paste0(dir,"Out_1_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
            if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_1_K",k,"_",i,".txt"), quote=F)}
            #
          }else if(model=="knn"){
            
          }else if(model=="lda"){
            lda.fit <- lda(popName ~ ., data=trainSetData_PC)
            lda.pred <- predict(lda.fit, testSetData_PC)
            lda.pred.class <- lda.pred$class
            lda.pred.prob <- lda.pred$posterior
            outcome_matrix <- cbind(testIndID, testSetData_PC$popName, as.data.frame(lda.pred.class), as.data.frame(lda.pred.prob))
            colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
            ##Output assignment results (prob.& used train loci) to files
            write.table(outcome_matrix, file=paste0(dir,"Out_1_K",k,"_",i,".txt"), quote=F, row.names=F )
            if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_1_K",k,"_",i,".txt"), quote=F)}
            #
          }else if(model=="naiveBayes"){
            nby.model <- naiveBayes(popName ~ ., data=trainSetData_PC)
            nby.pred.class <- predict(nby.model,testSetData_PC,type="class")
            nby.pred.prob <- predict(nby.model,testSetData_PC,type="raw")
            outcome_matrix <- cbind(testIndID, testSetData_PC$popName, as.data.frame(nby.pred.class), as.data.frame(nby.pred.prob))
            colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
            ##Output assignment results (prob.& used train loci) to files
            write.table(outcome_matrix, file=paste0(dir,"Out_1_K",k,"_",i,".txt"), quote=F, row.names=F )
            if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_1_K",k,"_",i,".txt"), quote=F)}
            #
          }else if(model=="tree"){
            tree.model <- tree(popName ~ ., data=trainSetData_PC)
            tree.pred.class <- predict(tree.model,testSetData_PC,type="class")
            tree.pred.prob <- predict(tree.model,testSetData_PC,type="vector")
            outcome_matrix <- cbind(testIndID, testSetData_PC$popName, as.data.frame(tree.pred.class), as.data.frame(tree.pred.prob))
            colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
            tree_node <- as.character(summary(tree.model)$used)#Get loci used at tree node
            ##Output assignment results (prob.& used train loci) to files
            write.table(outcome_matrix, file=paste0(dir,"Out_1_K",k,"_",i,".txt"), quote=F, row.names=F )
            cat("Variables used at tree node",tree_node,file=paste0(dir,"Loci_treenode_1_K",k,"_",i,".txt"), sep="\n")
            if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_1_K",k,"_",i,".txt"), quote=F)}
            #
          }else if(model=="randomForest"){
            rf.model <- randomForest(popName ~ ., data=trainSetData_PC, ntree=ntree, importance=T)
            rf.pred.class <- predict(rf.model,testSetData_PC,type="response")
            rf.pred.prob <- predict(rf.model,testSetData_PC,type="prob")
            outcome_matrix <- cbind(testIndID, testSetData_PC$popName, as.data.frame(rf.pred.class), as.data.frame(rf.pred.prob))
            colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
            ##Output assignment results (prob.& used train loci) to files
            write.table(outcome_matrix, file=paste0(dir,"Out_1_K",k,"_",i,".txt"), quote=F, row.names=F )
            write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Var_importance_1_K",k,"_",i,".txt"), quote=F, row.names=T)
            if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_1_K",k,"_",i,".txt"), quote=F)}
          }
        }#foreach(i=1:k,...)
        stopCluster(cl)
      }else{
        for(i in 1:k){
          trainSetData <- dataMatrix[-fold_index[[i]],]
          testIndID <- IndID[fold_index[[i]]]; testIndID <- droplevels(testIndID)
          #Create test dataset
          testSetData <- dataMatrix[fold_index[[i]],]
          #
          #Perform PCA on training data (pca.method=T)
          if(pca.method){
            PCA_results <- perform.PCA(trainSetData[,1:ncol(trainSetData)-1], method=pca.PCs)#Run PCA on training without label
            loadings <- PCA_results[[1]] #loadings (coefficience) of variables and PCs; apply this to test data
            trainSetData_PC <- as.data.frame(PCA_results[[2]])
            trainSetData_PC <- cbind(trainSetData_PC, trainSetData$popName) ##Will be used for building predicting models
            colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popName"
            #Convert test data to PC variables based on training's loadings
            testSetData_matrix <- as.matrix(testSetData[,1:ncol(testSetData)-1])
            testSetData_PC <- as.data.frame(testSetData_matrix %*% loadings)
            testSetData_PC <- cbind(testSetData_PC, testSetData$popName)
            colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popName"
          }else if(pca.method==F){
            trainSetData_PC <- trainSetData
            testSetData_PC <- testSetData
          }
          #Use training to build models and test on test individuals
          if(model=="svm"){
            svm.fit <- svm(popName ~ ., data=trainSetData_PC, kernel=svm.kernel, cost=svm.cost, prob=T, scale=F)
            svm.pred <- predict(svm.fit, testSetData_PC, type="class",prob=T)
            outcome_matrix <- cbind(testIndID, testSetData_PC$popName, as.data.frame(svm.pred), attr(svm.pred,"probabilities"))#combine output to data frame
            colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
            ##Output assignment results (prob.& used train loci) to files
            write.table(outcome_matrix, file=paste0(dir,"Out_1_K",k,"_",i,".txt"), quote=F, row.names=F )#File annotation: Out_(prop.loci)_(K)_(k-th Fold).txt
            if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_1_K",k,"_",i,".txt"), quote=F)}
            #
          }else if(model=="knn"){
            
          }else if(model=="lda"){
            lda.fit <- lda(popName ~ ., data=trainSetData_PC)
            lda.pred <- predict(lda.fit, testSetData_PC)
            lda.pred.class <- lda.pred$class
            lda.pred.prob <- lda.pred$posterior
            outcome_matrix <- cbind(testIndID, testSetData_PC$popName, as.data.frame(lda.pred.class), as.data.frame(lda.pred.prob))
            colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
            ##Output assignment results (prob.& used train loci) to files
            write.table(outcome_matrix, file=paste0(dir,"Out_1_K",k,"_",i,".txt"), quote=F, row.names=F )
            if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_1_K",k,"_",i,".txt"), quote=F)}
            #
          }else if(model=="naiveBayes"){
            nby.model <- naiveBayes(popName ~ ., data=trainSetData_PC)
            nby.pred.class <- predict(nby.model,testSetData_PC,type="class")
            nby.pred.prob <- predict(nby.model,testSetData_PC,type="raw")
            outcome_matrix <- cbind(testIndID, testSetData_PC$popName, as.data.frame(nby.pred.class), as.data.frame(nby.pred.prob))
            colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
            ##Output assignment results (prob.& used train loci) to files
            write.table(outcome_matrix, file=paste0(dir,"Out_1_K",k,"_",i,".txt"), quote=F, row.names=F )
            if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_1_K",k,"_",i,".txt"), quote=F)}
            #
          }else if(model=="tree"){
            tree.model <- tree(popName ~ ., data=trainSetData_PC)
            tree.pred.class <- predict(tree.model,testSetData_PC,type="class")
            tree.pred.prob <- predict(tree.model,testSetData_PC,type="vector")
            outcome_matrix <- cbind(testIndID, testSetData_PC$popName, as.data.frame(tree.pred.class), as.data.frame(tree.pred.prob))
            colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
            tree_node <- as.character(summary(tree.model)$used)#Get loci used at tree node
            ##Output assignment results (prob.& used train loci) to files
            write.table(outcome_matrix, file=paste0(dir,"Out_1_K",k,"_",i,".txt"), quote=F, row.names=F )
            cat("Variables used at tree node",tree_node,file=paste0(dir,"Loci_treenode_1_K",k,"_",i,".txt"), sep="\n")
            if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_1_K",k,"_",i,".txt"), quote=F)}
            #
          }else if(model=="randomForest"){
            rf.model <- randomForest(popName ~ ., data=trainSetData_PC, ntree=ntree, importance=T)
            rf.pred.class <- predict(rf.model,testSetData_PC,type="response")
            rf.pred.prob <- predict(rf.model,testSetData_PC,type="prob")
            outcome_matrix <- cbind(testIndID, testSetData_PC$popName, as.data.frame(rf.pred.class), as.data.frame(rf.pred.prob))
            colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
            ##Output assignment results (prob.& used train loci) to files
            write.table(outcome_matrix, file=paste0(dir,"Out_1_K",k,"_",i,".txt"), quote=F, row.names=F )
            write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Var_importance_1_K",k,"_",i,".txt"), quote=F, row.names=T)
            if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_1_K",k,"_",i,".txt"), quote=F)}
          }
        }#for(i in 1:k)...
      }#close for if(multiprocess){...}esle{...}
    }#for(k in k.fold)
    
    #
    version <- as.character(packageVersion("assignPOP"))
    cat(" Analysis Description ( R - assignPOP ver.",version,")\n",
        "Perform assign.kfold() @", format(Sys.time()),"\n\n",
        "k.fold =",k.fold,"\n",
        "Total assignment tests =",sum(k.fold),"\n",
        "Data scaled and centerd:",scaled,"\n",
        "PCA for dimensionality reduction:",pca.method,"\n",
        "PC retaining criteria:",pca.PCs,"\n",
        "Machine learning model:",model,"\n\n",
        "Input Data (Non-genetics)\n",
        "Number of individuals:",sum(popSizes),"\n",
        "Number of non-genetic variables:",noVars,"\n",
        "Number of categorical variables:",noFactorVar,"\n",
        "Number of numeric variable:",noVars-noFactorVar,"\n",
        "Number of populations:",noPops,"\n",
        names(popSizes),"\n",popSizes,"\n",
        file=paste0(dir,"AnalysisInfo.txt"))
    #Print some message to R console
    cat("\n  K-fold cross-validation done!!")
    cat(paste0("\n  ",sum(k.fold)," assignment tests completed!!"))
  }#else if(is.data.frame(x))
  
}# End