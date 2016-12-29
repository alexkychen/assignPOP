#' Population assignment test using Monte-Carlo cross-validation
#'
#' This function employs Monte-Carlo cross-validation for assignment tests. The results help evaluate if known data set has sufficient discriminatory power. It accepts genetic-only [object returned from read.genpop() or reducel.allele()], integrated [object returned from compile.data()], or non-genetic [R data frame with header] data as input, and outputs results to text files. Several built-in options are provided. See below for more details.
#' @param x An input object which should be the object (list) returned from the function read.genpop(), reduce.allele(), or compile.data(). It could also be a data frame (with column name) returned from read.csv() or read.table() if you're analyzing non-genetic data, such as morphormetrics, chemistry data. The non-genetic data frame should have sample ID in the first column and population label in the last column. 
#' @param train.inds The number (integer greater than 1) or proportion (float between 0 and 1) of individuals (observations) from each population to be used as training data. Use a numeric vector to specify multiple sets of training individuals. No mixture of integer and float in a vector.
#' @param train.loci The proportion (float between 0 and 1) of loci to be used as training data. Use a numeric vector to specify multiple sets of training loci. This argument will be ignored if you're analyzing non-genetic data.
#' @param loci.sample Locus sampling method, "fst" or "random". If loci.sample="fst" (default) and train.loci=0.1, it means that top 10 percent of high Fst loci will be sampled as training loci. On the other hand, if loci.sample="random", then random 10 percent of loci will be sampled as training loci. This argument will be ignored if you're analyzing non-genetic data.
#' @param iterations Resampling times (an integer) for each combination of training individuals and loci.
#' @param dir A character string to specify the folder name for saving output files. A slash at the end must be included (e.g., dir="YourFolderName/"). Otherwise, the files will be saved under your working directory.
#' @param scaled A logical variable (TRUE or FALSE) to specify whether to center (make mean of each feature to 0) and scale (make standard deviation of each feature to 1) the entire dataset before performing PCA and cross-validation. Default is FALSE. As genetic data has converted to numeric data between 0 and 1, to scale or not to scale the genetic data should not be critical. However, it is recommended to set scaled=TRUE when integrated data contains various scales of features.  
#' @param pca.method Either a character string ("mixed", "independent", or "original") or logical variable (TRUE or FALSE) to specify how to perform PCA on non-genetic data (PCA is always performed on genetic data). The character strings are used when analyzing integrated (genetic plus non-genetic) data. If using "mixed" (default), PCA is perfromed across the genetic and non-genetic data, resulting in each PC summarizing mixed variations of genetic and non-genetic data. If using "independent", PCA is independently performed on non-genetic data. Genetic PCs and non-genetic PCs are then used as new features. If using "original", original non-genetic data and genetic PCs are used as features. The logical variable is used when analyzing non-genetic data.If TRUE, it performs PCA on the training data and applys the loadings to the test data. Scores of training and test data will be used as new features. 
#' @param pca.PCs A criterion to retain number of PCs. By default, it uses Kaiser-Guttman criterion that any PC has the eigenvalue greater than 1 will be retained as the new variable/feature. Users can set an integer to specify the number of PCs to be retained.
#' @param pca.loadings A logical variable (TRUE or FALSE) to determine whether to output the loadings of training data to text files. Default is FALSE. Just a heads-up, the output files could take some storage space, if set TRUE.
#' @param model A character string to specify which classifier to use for creating predictive models. The current options include "lda", "svm", "naiveBayes", "tree", and "randomForest". Default is "svm"(support vector machine).
#' @param svm.kernel A character string to specify which kernel to be used when using "svm" classifier.
#' @param svm.cost A number to specify the cost for "svm" method.
#' @param ntree A integer to specify how many trees to build when using "randomForest" method.
#' @param processors The number of processors to be used for parallel running. By default, it uses N-1 processors in your computer.
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
assign.MC <- function(x, train.inds=c(0.5,0.7,0.9), train.loci=c(0.1,0.25,0.5, 1), loci.sample="fst", iterations=20, dir=NULL,
                      scaled=FALSE, pca.method="mixed", pca.PCs="kaiser-guttman", pca.loadings=F, model="svm", svm.kernel="linear", svm.cost=1, ntree=50, processors=999, ...){
  #check if dir is correctly entered
  if(is.null(dir)){
    stop("Please provide a folder name ending with '/' in argument 'dir' ")
  }else if(substr(dir, start=nchar(dir), stop=nchar(dir))!="/"){
    stop("Please put a forward slash '/' in the end of your folder name (in argument 'dir'). ")
  }
  ##Check data types
  if(!is.data.frame(x)){ #check if input x is a list returned from read.genpop(), reduce.allele(), or compile.data()
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
    #Create a folder to save outfiles
    dir.create(file.path(dir))
    # Detect CPU core numbers and use n-1 cores to run the program below
    maxCores <- detectCores()-1
    if (processors <= maxCores & processors > 0){
      cl <- makeCluster(processors)
      registerDoParallel(cl,cores=processors)
    }else {
      cl <- makeCluster(maxCores)
      registerDoParallel(cl,cores=maxCores)
    }
    cat("  ",maxCores,"cores/threads of CPU will be used for analysis...\n")
    #Determine if object x includes non-genetic data. If so, length of x (list) will be 5. If not (only genetic data), the length will be 3.
    if(length(x)==3){ #Process genetics-only data
      datatype <- "genetics";noVars <- 0; pca.method <- "NA"
      # Run program in parallel using multi-core if it will
      foreach(i=1:iterations, .export=c("Fsts","perform.PCA"), .packages=c("e1071","klaR","MASS","tree","randomForest")) %dopar% {
        for(j in 1:length(train.inds)){
          trainSet_index <- NULL
          if(all(train.inds < 1)){ #sample training individuals per pop by proportion
            for(p in pops){
              onePopRowName <- row.names(subset(genoMatrix,popNames_vector==p))
              sampleOnePopRow <- sample(onePopRowName, round(length(onePopRowName)*train.inds[j]))
              trainSet_index <- c(trainSet_index,sampleOnePopRow)
            }
          } else if(all(train.inds > 1)){ #sample training individuals per pop by fixed number
            for(p in pops){
              onePopRowName <- row.names(subset(genoMatrix,popNames_vector==p))
              sampleOnePopRow <- sample(onePopRowName, train.inds[j])
              trainSet_index <- c(trainSet_index,sampleOnePopRow)
            }
          }
          #
          trainSet_index <- sort(as.numeric(trainSet_index))
          trainSetMatrix <- genoMatrix[trainSet_index,]
          testIndID <- x[[2]][-trainSet_index]#Get test individual ID for later print out
          
          #Check how loci should be sampled, if "prior"(by default), calculate and sort Fst loci
          if(loci.sample=="fst"){
            train_X <- list(trainSetMatrix, row.names(trainSetMatrix), x[[3]])#create a new x list for Fsts function; x[[3]] is locus name
            fstTable <- Fsts(train_X)#Estimate locus Fst for training data
            orderFst <- order(fstTable$Fst, decreasing=T)#Order the loci index by Fst value(from highest)
            #Loop through each fst level
            for(k in train.loci){
              trainLocusIndex_fstTable <- orderFst[1:round(noLocus*k)]#Get training locus index (from fstTable)
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
              trainSetData <- genoMatrix_trainLoci[trainSet_index,]#Get the training set data (training individuals/loci)
              testSetData <- genoMatrix_trainLoci[-trainSet_index,]#Get the test set data (test individuals/loci)
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
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
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
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
                #
              }else if(model=="naiveBayes"){
                nby.model <- naiveBayes(popNames_vector ~ ., data=trainSetData_PC)
                nby.pred.class <- predict(nby.model,testSetData_PC,type="class")
                nby.pred.prob <- predict(nby.model,testSetData_PC,type="raw")
                outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(nby.pred.class), as.data.frame(nby.pred.prob))
                colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                ##Output assignment results (prob.& used train loci) to files
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
                #
              }else if(model=="tree"){
                tree.model <- tree(popNames_vector ~ ., data=trainSetData_PC)
                tree.pred.class <- predict(tree.model,testSetData_PC,type="class")
                tree.pred.prob <- predict(tree.model,testSetData_PC,type="vector")
                outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(tree.pred.class), as.data.frame(tree.pred.prob))
                colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                tree_node <- as.character(summary(tree.model)$used)#Get loci used at tree node
                ##Output assignment results (prob.& used train loci) to files
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                cat("Features used at tree node",tree_node,file=paste0(dir,"Feature_treenode_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
                #
              }else if(model=="randomForest"){
                rf.model <- randomForest(popNames_vector ~ ., data=trainSetData_PC, ntree=ntree, importance=T)
                rf.pred.class <- predict(rf.model,testSetData_PC,type="response")
                rf.pred.prob <- predict(rf.model,testSetData_PC,type="prob")
                outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(rf.pred.class), as.data.frame(rf.pred.prob))
                colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                ##Output assignment results (prob.& used train loci) to files
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                #Output "importance"(Gini index) of loci to files;see randomForest package's "importance" function
                write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Feature_importance_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=T)
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
                #
              }
            }#for(k in train.loci)
            
          }else if(loci.sample=="random"){
            # Random select a proption of loci (ignore Fst)
            for(k in train.loci){
              tempLocusIndex <- sort(sample(1:noLocus, round(noLocus*k)))
              trainLocusName <- locusNames[tempLocusIndex]
              trainLocusIndex_genoMatrix <- NULL #create a train locus index for genoMatrix to be extracted
              for(m in 1:length(trainLocusName)){
                tempAlleleIndex <- grep(pattern=paste0(trainLocusName[m],"_"), alleleName)
                trainLocusIndex_genoMatrix <- c(trainLocusIndex_genoMatrix,tempAlleleIndex )
              }
              trainLocusIndex_genoMatrix <- sort(trainLocusIndex_genoMatrix)
              genoMatrix_trainLoci <- genoMatrix[,c(trainLocusIndex_genoMatrix)]#resahpe the genoMatrix data, comprising only training loci
              #Scale and center the data if "scaled=T"
              if(scaled){
                genoMatrix_trainLoci <- as.data.frame(scale(genoMatrix_trainLoci))
              }
              #Add the pop.Names to the last column
              genoMatrix_trainLoci <- cbind(genoMatrix_trainLoci, genoMatrix[,ncol(genoMatrix)]);colnames(genoMatrix_trainLoci)[ncol(genoMatrix_trainLoci)] <- "popNames_vector"
              trainSetData <- genoMatrix_trainLoci[trainSet_index,]#Get the training set data (training individuals/loci)
              testSetData <- genoMatrix_trainLoci[-trainSet_index,]#Get the test set data (test individuals/loci)
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
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
                #
              }else if(model=="knn"){
                #
              }else if(model=="lda"){
                lda.fit <- lda(popNames_vector ~ ., data=trainSetData_PC)
                lda.pred <- predict(lda.fit, testSetData_PC)
                lda.pred.class <- lda.pred$class
                lda.pred.prob <- lda.pred$posterior
                outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(lda.pred.class), as.data.frame(lda.pred.prob))
                colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                ##Output assignment results (prob.& used train loci) to files
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
                #
              }else if(model=="naiveBayes"){
                nby.model <- naiveBayes(popNames_vector ~ ., data=trainSetData_PC)
                nby.pred.class <- predict(nby.model,testSetData_PC,type="class")
                nby.pred.prob <- predict(nby.model,testSetData_PC,type="raw")
                outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(nby.pred.class), as.data.frame(nby.pred.prob))
                colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                ##Output assignment results (prob.& used train loci) to files
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
                #
              }else if(model=="tree"){
                tree.model <- tree(popNames_vector ~ ., data=trainSetData_PC)
                tree.pred.class <- predict(tree.model,testSetData_PC,type="class")
                tree.pred.prob <- predict(tree.model,testSetData_PC,type="vector")
                outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(tree.pred.class), as.data.frame(tree.pred.prob))
                colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                tree_node <- as.character(summary(tree.model)$used)#Get loci used at tree node
                ##Output assignment results (prob.& used train loci) to files
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                cat("Features used at tree node",tree_node,file=paste0(dir,"Feature_treenode_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
                #
              }else if(model=="randomForest"){
                rf.model <- randomForest(popNames_vector ~ ., data=trainSetData_PC, ntree=ntree, importance=T)
                rf.pred.class <- predict(rf.model,testSetData_PC,type="response")
                rf.pred.prob <- predict(rf.model,testSetData_PC,type="prob")
                outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(rf.pred.class), as.data.frame(rf.pred.prob))
                colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                ##Output assignment results (prob.& used train loci) to files
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                #Output "importance"(Gini index) of loci to files;see randomForest package's "importance" function
                write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Feature_importance_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=T)
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
                #
              }
            }#for(k in train.loci)
          }#else if(loci.sample=="random"
        } # for(j in 1:length(train.inds))
      } # foreach(i=1:iterations, ...
      
    }else if(length(x)==5){ #Process genetic and non-genetic data
      datatype <- "genetics + non-genetics"
      #Split geneitc and non-genetic data first for estimating Fst
      otherVarNames <- x[[4]]#Get non-genetic variable name
      noVars <- length(otherVarNames)
      #
      foreach(i=1:iterations, .export=c("Fsts","perform.PCA"), .packages=c("e1071","klaR","MASS","tree","randomForest")) %dopar% {
        for(j in 1:length(train.inds)){
          trainSet_index <- NULL
          if(all(train.inds < 1)){ #sample training individuals per pop by proportion
            for(p in pops){
              onePopRowName <- row.names(subset(genoMatrix,popNames_vector==p))
              sampleOnePopRow <- sample(onePopRowName, round(length(onePopRowName)*train.inds[j]))
              trainSet_index <- c(trainSet_index,sampleOnePopRow)
            }
          } else if(all(train.inds > 1)){ #sample training individuals per pop by fixed number
            for(p in pops){
              onePopRowName <- row.names(subset(genoMatrix,popNames_vector==p))
              sampleOnePopRow <- sample(onePopRowName, train.inds[j])
              trainSet_index <- c(trainSet_index,sampleOnePopRow)
            }
          }
          #
          trainSet_index <- sort(as.numeric(trainSet_index))
          trainSetMatrix <- genoMatrix[trainSet_index,] #test inds. with all variables and popNames_vector
          testIndID <- x[[2]][-trainSet_index]#Get test individual ID for later print out
          #Separate genetic and non-genetic data
          firstOtherVars_pos <- length(genoMatrix) - noVars #column position of first non-geneitc var
          lastOtherVars_pos <- length(genoMatrix) - 1 #column position of last non-genetic var
          trainSetMatrix_genetic <- trainSetMatrix[,-(firstOtherVars_pos:lastOtherVars_pos)] #test inds. with only genetic variables and popNames
          #trainSetMatrix_otherVars <- trainSetMatrix[,ncol(genoMatrix)-noVars:ncol(genoMatrix)]#test inds. with non-genetic and popNames
          
          if(loci.sample=="fst"){
            train_X <- list(trainSetMatrix_genetic, row.names(trainSetMatrix_genetic), x[[3]])#create a new x list for Fsts function; x[[3]] is locus name
            fstTable <- Fsts(train_X)#Estimate locus Fst for training data
            orderFst <- order(fstTable$Fst, decreasing=T)#Order the loci index by Fst value(from highest)
            #Loop through each fst level
            for(k in train.loci){
              trainLocusIndex_fstTable <- orderFst[1:round(noLocus*k)]#Get training locus index (from fstTable)
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
              trainSetData <- genoMatrix_trainVars[trainSet_index,] #Get the training set data (training individuals/loci)
              testSetData <- genoMatrix_trainVars[-trainSet_index,] #Get the test set data (test individuals/loci)
              
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
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
                #
              }else if(model=="knn"){
                #
              }else if(model=="lda"){
                lda.fit <- lda(popNames_vector ~ ., data=trainSetData_PC)
                lda.pred <- predict(lda.fit, testSetData_PC)
                lda.pred.class <- lda.pred$class
                lda.pred.prob <- lda.pred$posterior
                outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(lda.pred.class), as.data.frame(lda.pred.prob))
                colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                ##Output assignment results (prob.& used train loci) to files
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
                #
              }else if(model=="naiveBayes"){
                nby.model <- naiveBayes(popNames_vector ~ ., data=trainSetData_PC)
                nby.pred.class <- predict(nby.model,testSetData_PC,type="class")
                nby.pred.prob <- predict(nby.model,testSetData_PC,type="raw")
                outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(nby.pred.class), as.data.frame(nby.pred.prob))
                colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                ##Output assignment results (prob.& used train loci) to files
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
                #
              }else if(model=="tree"){
                tree.model <- tree(popNames_vector ~ ., data=trainSetData_PC)
                tree.pred.class <- predict(tree.model,testSetData_PC,type="class")
                tree.pred.prob <- predict(tree.model,testSetData_PC,type="vector")
                outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(tree.pred.class), as.data.frame(tree.pred.prob))
                colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                tree_node <- as.character(summary(tree.model)$used)#Get loci used at tree node
                ##Output assignment results (prob.& used train loci) to files
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                cat("Features used at tree node",tree_node,file=paste0(dir,"Feature_treenode_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
                #
              }else if(model=="randomForest"){
                rf.model <- randomForest(popNames_vector ~ ., data=trainSetData_PC, ntree=ntree, importance=T)
                rf.pred.class <- predict(rf.model,testSetData_PC,type="response")
                rf.pred.prob <- predict(rf.model,testSetData_PC,type="prob")
                outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(rf.pred.class), as.data.frame(rf.pred.prob))
                colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                ##Output assignment results (prob.& used train loci) to files
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                #Output "importance"(Gini index) of loci to files;see randomForest package's "importance" function
                write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Feature_importance_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=T)
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
                #
              }
            }#for(k in train.loci)
          }else if(loci.sample=="random"){
            # Random select a proption of loci (ignore Fst)
            for(k in train.loci){
              tempLocusIndex <- sort(sample(1:noLocus, round(noLocus*k)))
              trainLocusName <- locusNames[tempLocusIndex]
              trainLocusIndex_genoMatrix <- NULL
              for(m in 1:length(trainLocusName)){
                tempAlleleIndex <- grep(pattern=paste0(trainLocusName[m],"_"), alleleName)
                trainLocusIndex_genoMatrix <- c(trainLocusIndex_genoMatrix,tempAlleleIndex )
              }
              trainLocusIndex_genoMatrix <- sort(trainLocusIndex_genoMatrix)
              genoMatrix_trainVars <- genoMatrix[,c(trainLocusIndex_genoMatrix,firstOtherVars_pos:lastOtherVars_pos)]#resahpe the genoMatrix data, comprising training loci + non-genetic + popNames
              #Scale and center the entire dataset if scaled=T
              if(scaled){
                genoMatrix_trainVars <- as.data.frame(scale(genoMatrix_trainVars))
              }
              #Add popNames_vector to the last column
              genoMatrix_trainVars <- cbind(genoMatrix_trainVars, genoMatrix$popNames_vector);colnames(genoMatrix_trainVars)[ncol(genoMatrix_trainVars)] <- "popNames_vector"
              ##Split entire data into training and test sets
              trainSetData <- genoMatrix_trainVars[trainSet_index,]#Get the training set data (training individuals/loci)
              testSetData <- genoMatrix_trainVars[-trainSet_index,]#Get the test set data (test individuals/loci)
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
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
                #
              }else if(model=="knn"){
                
              }else if(model=="lda"){
                lda.fit <- lda(popNames_vector ~ ., data=trainSetData_PC)
                lda.pred <- predict(lda.fit, testSetData_PC)
                lda.pred.class <- lda.pred$class
                lda.pred.prob <- lda.pred$posterior
                outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(lda.pred.class), as.data.frame(lda.pred.prob))
                colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                ##Output assignment results (prob.& used train loci) to files
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
                #
              }else if(model=="naiveBayes"){
                nby.model <- naiveBayes(popNames_vector ~ ., data=trainSetData_PC)
                nby.pred.class <- predict(nby.model,testSetData_PC,type="class")
                nby.pred.prob <- predict(nby.model,testSetData_PC,type="raw")
                outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(nby.pred.class), as.data.frame(nby.pred.prob))
                colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                ##Output assignment results (prob.& used train loci) to files
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
                #
              }else if(model=="tree"){
                tree.model <- tree(popNames_vector ~ ., data=trainSetData_PC)
                tree.pred.class <- predict(tree.model,testSetData_PC,type="class")
                tree.pred.prob <- predict(tree.model,testSetData_PC,type="vector")
                outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(tree.pred.class), as.data.frame(tree.pred.prob))
                colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                tree_node <- as.character(summary(tree.model)$used)#Get loci used at tree node
                ##Output assignment results (prob.& used train loci) to files
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                cat("Features used at tree node",tree_node,file=paste0(dir,"Feature_treenode_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
                #
              }else if(model=="randomForest"){
                rf.model <- randomForest(popNames_vector ~ ., data=trainSetData_PC, ntree=ntree, importance=T)
                rf.pred.class <- predict(rf.model,testSetData_PC,type="response")
                rf.pred.prob <- predict(rf.model,testSetData_PC,type="prob")
                outcome_matrix <- cbind(testIndID, testSetData_PC$popNames_vector, as.data.frame(rf.pred.class), as.data.frame(rf.pred.prob))
                colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
                ##Output assignment results (prob.& used train loci) to files
                write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=F )
                cat(trainLocusName, file=paste0(dir,"Loci_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
                #Output "importance"(Gini index) of loci to files;see randomForest package's "importance" function
                write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Feature_importance_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=T)
                if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
              }
            }#for(k in train.loci)
          }#else if(loci.sample=="random")
        }#for(j in 1:length(train.inds))
      }#foreach(i=1:iterations,...
    }#else if(length(x)==4)
    stopCluster(cl)
    #Output a metadata file
    cat(" Analysis Description (R - assignPOP ver.",packageVersion("assignPOP"),")\n",
        "Perform assign.MC() @", format(Sys.time()),"\n\n",
        "train.inds =",train.inds,"\n",
        "train.loci =",train.loci,"\n",
        "iterations =",iterations,"\n",
        "Total assignment tests =",length(train.inds)*length(train.loci)*iterations,"\n",
        "fst locus sample method:",loci.sample,"\n",
        "Data scaled and centerd:",scaled,"\n",
        "PC retaining criteria:",pca.PCs,"\n",
        "PCA for non-genetic data:",pca.method,"\n",
        "Machine-learning model:",model,"\n\n",
        "Input Data (",datatype,")\n",
        "Number of individuals:",sum(popSizes),"\n",
        "Number of loci:",noLocus,"\n",
        "Number of non-genetic variables:",noVars,"\n",
        "Number of populations:",noPops,"\n",
        names(popSizes),"\n",popSizes,"\n",
        file=paste0(dir,"AnalysisInfo.txt"))
    #Print some message to R console
    cat("\n  Monte-Carlo cross-validation done!!")
    cat(paste0("\n  ",length(train.inds)*length(train.loci)*iterations," assignment tests completed!!"))
    #Print a warning message if only one PC is used (when using Kaiser-Guttman criteria), meaning none of PC's eigenvalue greater than 1.
    #if(!is.null(PCA_results[[3]])){
    #  message(PCA_results[[3]])
    #}
    
  }else if(is.data.frame(x)){ #check if input x is a data frame returned from read.csv() or read.table()
    #Analyze non-genetic data
    #checking pca.method
    if(is.character(pca.method)){
      anspca <- readline("  Perform PCA on dataset for dimensionality reduction? (enter Y/N): ")
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
    #Create a folder to save outfiles
    dir.create(file.path(dir))
    # Detect CPU core numbers and use n-1 cores to run the program below
    maxCores <- detectCores()-1
    if (processors <= maxCores & processors > 0){
      cl <- makeCluster(processors)
      registerDoParallel(cl,cores=processors)
    }else {
      cl <- makeCluster(maxCores)
      registerDoParallel(cl,cores=maxCores)
    }
    cat("  ",maxCores,"cores/threads of CPU will be used for analysis...\n")
    #Check data type (numeric or categorical)
    for(n in 1:noVars){
      var_name <- varNames[n]
      var_type <- class(dataMatrix[,n])
      cat(paste0("  ",var_name,"(",var_type,")"))
    }
    ans0 <- readline("  Are they correct? (enter Y/N): ")
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
    #Scale and center the entire features if set scaled=T
    if(scaled){
      cat("\n  Scaling and centering entire data set...")
      dataMatrix_scaled <- as.data.frame(scale(dataMatrix[,1:(ncol(dataMatrix)-1)]))
      dataMatrix <- cbind(dataMatrix_scaled, dataMatrix$popName);colnames(dataMatrix)[ncol(dataMatrix)] <- "popName"
    }
    #Start resampling cross-validation using multi cores/threads if desired
    foreach(i=1:iterations, .export="perform.PCA", .packages=c("e1071","klaR","MASS","tree","randomForest")) %dopar% {
      for( j in 1:length(train.inds)){
        trainSet_index <- NULL
        if(all(train.inds < 1)){
          for(p in pops){
            onePopRowName <- row.names(subset(dataMatrix, popName==p))
            sampleOnePopRow <- sample(onePopRowName, round(length(onePopRowName)*train.inds[j]))
            trainSet_index <- c(trainSet_index,sampleOnePopRow)
          }#for(p in pops)
        }else if(all(train.inds > 1)){
          for(p in pops){
            onePopRowName <- row.names(subset(dataMatrix, popName==p))
            sampleOnePopRow <- sample(onePopRowName, train.inds[j])
            trainSet_index <- c(trainSet_index, sampleOnePopRow)
          }
        }#if(all(train.inds < 1))
        #Create training dataset and get test ind id
        trainSet_index <- sort(as.numeric(trainSet_index))
        trainSetMatrix <- dataMatrix[trainSet_index,]
        testIndID <- IndID[-trainSet_index]; testIndID <- droplevels(testIndID)
        #Create test dataset
        testSetData <- dataMatrix[-trainSet_index,]
        
        #Determine whether perform PCA on training
        if(pca.method==T){
          PCA_results <- perform.PCA(trainSetMatrix[,1:ncol(trainSetMatrix)-1] ,method=pca.PCs)#Run PCA exclude last column ID
          loadings <- PCA_results[[1]] #loadings (coefficience) of variables and PCs; apply this to test data
          trainSetData_PC <- as.data.frame(PCA_results[[2]])
          trainSetData_PC <- cbind(trainSetData_PC, trainSetMatrix$popName) ##Will be used for building predicting models
          colnames(trainSetData_PC)[ncol(trainSetData_PC)] <- "popName"
          #Convert test data to PC variables based on training's loadings
          testSetData_matrix <- as.matrix(testSetData[,1:ncol(testSetData)-1])
          testSetData_PC <- as.data.frame(testSetData_matrix %*% loadings)
          testSetData_PC <- cbind(testSetData_PC, testSetData$popName)
          colnames(testSetData_PC)[ncol(testSetData_PC)] <- "popName"
        }else if(pca.method==F){
          trainSetData_PC <- trainSetMatrix
          testSetData_PC <- testSetData
        }
        #Use training to build predictive models and test on test individuals
        if(model=="svm"){
          svm.fit <- svm(popName ~ ., data=trainSetData_PC, kernel=svm.kernel, cost=svm.cost, prob=T, scale=F)
          svm.pred <- predict(svm.fit, testSetData_PC, type="class",prob=T)
          outcome_matrix <- cbind(testIndID, testSetData_PC$popName, as.data.frame(svm.pred), attr(svm.pred,"probabilities"))#combine output to data frame
          colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
          ##Output assignment results (prob.& used train loci) to files
          write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_1_",i,".txt"), quote=F, row.names=F )
          if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_1_",i,".txt"), quote=F)}
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
          write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_1_",i,".txt"), quote=F, row.names=F )
          if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_1_",i,".txt"), quote=F)}
          #
        }else if(model=="naiveBayes"){
          nby.model <- naiveBayes(popName ~ ., data=trainSetData_PC)
          nby.pred.class <- predict(nby.model,testSetData_PC,type="class")
          nby.pred.prob <- predict(nby.model,testSetData_PC,type="raw")
          outcome_matrix <- cbind(testIndID, testSetData_PC$popName, as.data.frame(nby.pred.class), as.data.frame(nby.pred.prob))
          colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
          ##Output assignment results (prob.& used train loci) to files
          write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_1_",i,".txt"), quote=F, row.names=F )
          if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_1_",i,".txt"), quote=F)}
          #
        }else if(model=="tree"){
          tree.model <- tree(popName ~ ., data=trainSetData_PC)
          tree.pred.class <- predict(tree.model,testSetData_PC,type="class")
          tree.pred.prob <- predict(tree.model,testSetData_PC,type="vector")
          outcome_matrix <- cbind(testIndID, testSetData_PC$popName, as.data.frame(tree.pred.class), as.data.frame(tree.pred.prob))
          colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
          tree_node <- as.character(summary(tree.model)$used)#Get loci used at tree node
          ##Output assignment results (prob.& used train loci) to files
          write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_1_",i,".txt"), quote=F, row.names=F )
          cat("Variables used at tree node",tree_node,file=paste0(dir,"Var_treenode_",train.inds[j],"_1_",i,".txt"), sep="\n")
          if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_1_",i,".txt"), quote=F)}
          #
        }else if(model=="randomForest"){
          rf.model <- randomForest(popName ~ ., data=trainSetData_PC, ntree=ntree, importance=T)
          rf.pred.class <- predict(rf.model,testSetData_PC,type="response")
          rf.pred.prob <- predict(rf.model,testSetData_PC,type="prob")
          outcome_matrix <- cbind(testIndID, testSetData_PC$popName, as.data.frame(rf.pred.class), as.data.frame(rf.pred.prob))
          colnames(outcome_matrix)[1:3] <- c("Ind.ID","origin.pop","pred.pop")
          ##Output assignment results (prob.& used train loci) to files
          write.table(outcome_matrix, file=paste0(dir,"Out_",train.inds[j],"_1_",i,".txt"), quote=F, row.names=F )
          write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Var_importance_",train.inds[j],"_1_",i,".txt"), quote=F, row.names=T)
          if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_1_",i,".txt"), quote=F)}
        }
      }#for( j in 1:length(train.inds))
    }#foreach
    stopCluster(cl)
    #Output a metadata file
    cat(" Analysis Description (R - assignPOP ver.",packageVersion("assignPOP"),")\n",
        "Perform assign.MC() @", format(Sys.time()),"\n\n",
        "train.inds =",train.inds,"\n",
        "iterations =",iterations,"\n",
        "Total assignment tests =",length(train.inds)*iterations,"\n",
        "Data scaled and centerd:",scaled,"\n",
        "PCA for dimensionality reduction:",pca.method,"\n",
        "PC retaining criteria:",pca.PCs,"\n",
        "Machine-learning model:",model,"\n\n",
        "Input Data (Non-genetics)\n",
        "Number of individuals:",sum(popSizes),"\n",
        "Number of non-genetic variables:",noVars,"\n",
        "Number of categorical variables:",noFactorVar,"\n",
        "Number of numeric variable:",noVars-noFactorVar,"\n",
        "Number of populations:",noPops,"\n",
        names(popSizes),"\n",popSizes,"\n",
        file=paste0(dir,"AnalysisInfo.txt"))
    #Print some message to R console
    cat("\n  Monte-Carlo cross-validation done!!")
    cat(paste0("\n  ",length(train.inds)*iterations," assignment tests completed!!"))

  }#else if(is.data.frame(x))
  
}#END


### Function to estimate locus Fst (following Nei 1973)
Fsts <- function(x){
  popNames_vector <- NULL
  genoMatrix <- x[[1]] #Get genoMatrix data
  genoNames <- names(genoMatrix)#Get column names of genoMatrix
  locusNames <- x[[3]]#Get locus name
  noLocus <- length(locusNames)
  pops <- names(table(genoMatrix$popNames_vector))#Get what pops in data
  locusFst <- as.data.frame(matrix(ncol=2, nrow=0))#Create a data frame to save locus Fst
  colnames(locusFst) <- c("Locus","Fst")
  #Start to process for each locus
  for(i in 1:noLocus){
    getLocusName <- locusNames[i]#Get a locus name as pattern for recognizing genoMatrix columns
    getAlleleIndex <- grep(pattern=getLocusName, genoNames)#Get allele index for extracing locus data from genoMatrix
    oneLocusMatrix <- genoMatrix[,c(getAlleleIndex,ncol(genoMatrix))]#Extract one locus data plus last column (pop name)
    tempTable <- data.frame(matrix(ncol=0,nrow=1)) #Empty table to save pop heterozygosity
    #subset data for each pop
    for(p in pops){
      onePopMatrix <- subset(oneLocusMatrix,popNames_vector==p)
      alleleSum <- colSums(onePopMatrix[1:ncol(onePopMatrix)-1])#sum the value for each allele
      alleleFreq <- (alleleSum)/sum(alleleSum)#estimate allele freq in a pop
      heterozygosity <- 1 - sum(alleleFreq^2)
      heteroCell <- as.data.frame(heterozygosity)
      colnames(heteroCell) <- p
      tempTable <- cbind(tempTable, heteroCell)
    }
    #Estimate mean of heterozygosity across pops (observed, Hs)
    Hs <- mean(as.numeric(tempTable[1,]))
    #Estimate allele frequency across all pops/samples
    alleleSums <- colSums(oneLocusMatrix[1:ncol(oneLocusMatrix)-1])
    alleleFreqs <- (alleleSums)/sum(alleleSums)#estimate allele freq
    #Estimate mean of expected heterozygosity across pops (expected, Ht)
    Ht <- 1 - sum(alleleFreqs^2)
    #Calculate Fst
    Fst <- (Ht - Hs)/Ht
    #result <- c(getLocusName, Fst)
    locusFst <- rbind(locusFst, data.frame(Locus=getLocusName, Fst=Fst))
  }

  return(locusFst)
}

#### Function to convert independent variables to PC variables ####
perform.PCA <- function(matrix, method){
  #matrix <- scale(matrix) #Scaling data (center to zero)
  PCA_results <- prcomp(matrix, scale=F, center=F)
  warnMessage <- NULL
  if(method=="kaiser-guttman"){ #Select PCs which eigenvalue greater than 1
    noPC_eig_g1 <- length(PCA_results$sdev[PCA_results$sdev > 1])#count how many PC's eigenvalue greater than 1
    if(noPC_eig_g1==0){
      warnMessage <- "None of PC's eigenvalue greater than 1. Only the first PC is used. Consider to assign a number to 'pca.PCs'."
    }
    #Extract first 'noPC_eig_g1' PCs of variables ($rotation, aka loadings)
    var_mx <- PCA_results$rotation[,1:noPC_eig_g1] #class::matrix
    #Estimate scores for observation based on the PCs
    obs_mx <- as.matrix(matrix) %*% var_mx #This should be the same matrix as PCA_results$x (score matrix)
    return(list(var_mx, obs_mx))

  }else if(method=="broken-stick"){ #Select PC which eigenvalue greater than broken-stick value
    #future development
  }else if(is.numeric(method)){
    noPCs <- method
    #Extract first number of PCs ($rotation, aka loadings)
    var_mx <- PCA_results$rotation[,1:noPCs] #class::matrix
    #Estimate scores for observation based on the PCs
    obs_mx <- as.matrix(matrix) %*% var_mx #This should be the same matrix as PCA_results$x (score matrix)
  }
  return(list(var_mx, obs_mx, warnMessage))
}

