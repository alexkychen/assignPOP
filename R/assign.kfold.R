#' Population assignment test using K-fold cross-validation
#'
#' This function employs K-fold cross-validation for assignment tests. It uses the object returned from the function read.genpop(), reduce.allele(), or compile.data() as input, and outputs results to text files. Typing a folder name in the argument "dir" is recommended. See below for more details.
#' @param x An input object which should be the object (list) returned from the function read.genpop(), reduce.allele(), or compile.data().
#' @param k.fold The number of groups to be divided for each population. Use a numeric vector to specify multiple sets of k-folds.
#' @param train.loci The proportion (float between 0 and 1) of loci to be used as training data. Use a numeric vector to specify multiple sets of training loci.
#' @param loci.sample Locus sampling method, "fst" or "random". If loci.sample="fst" (default) and train.loci=0.1, it means that top 10 percent of high Fst loci will be sampled as training loci. On the other hand, if loci.sample="random", then random 10 percent of loci will be sampled as training loci.
#' @param dir A character string to specify the folder name for saving output files. A slash at the end must be included (e.g., dir="YourFolderName/"). Otherwise, the files will be saved under your working directory.
#' @param pca.method A criterion to retain number of PCs. By default, it uses Kaiser-Guttman criterion that any PC has the eigenvalue greater than 1 will be retained as the new variable/feature. Users can set an integer to specify the number of PCs to be retained.
#' @param pca.loadings A logical variable (False or True) to determine whether it prints the loadings of training data to output text files. Default is False, if set True, the overall output files could be large.
#' @param model A character string to specify which classifier to use for creating predictive models. The current options include "lda", "svm", "naiveBayes", "tree", and "randomForest".
#' @param svm.kernel A character string to specify which kernel to be used when using "svm" classifier.
#' @param svm.cost A number to specify the cost for "svm" method.
#' @param ntree A integer to specify how many trees to build when using "randomForest" method.
#' @param processors The number of processors to be used for parallel running. By default, it uses N-1 processors in your computer.
#' @return You don't need to specify a name for the returned object when using this function. It automatically output results in text files to your designated folder.
#' @export
#'
assign.kfold <- function(x, k.fold = c(3,4,5), train.loci=c(0.1,0.25,0.5, 1), loci.sample="fst", dir=NULL,
                         pca.method="kaiser-guttman", pca.loadings=F, model="svm", svm.kernel="linear", svm.cost=1, ntree=50, processors=999, ...){
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
  #Determine if object x includes non-genetic data. If so, length of x (list) will be 5. If not (only genetic data), the length will be 3.
  if(length(x)==3){#Process genetics-only data
    datatype <- "genetics";noVars <- 0
    #loop each k fold
    for(k in k.fold){
      #create fold index based on label (genoMatrix$popNames_vector)
      fold_index <- createFolds(genoMatrix$popNames_vector, k=k) #this is a list containing k folds of index, e.g., fold_index$Fold1, $Fold2...
      #take the K-th fold as test set; remaing folds as training set
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
              tempAlleleIndex <- grep(pat=paste0(trainLocusName[m],"_"), alleleName)#alleleName is colnames(genoMatrix)
              trainLocusIndex_genoMatrix <- c(trainLocusIndex_genoMatrix,tempAlleleIndex )
            }
            trainLocusIndex_genoMatrix <- sort(trainLocusIndex_genoMatrix)
            genoMatrix_trainLoci <- genoMatrix[,c(trainLocusIndex_genoMatrix, ncol(genoMatrix))]#resahpe the genoMatrix data, comprising only training loci
            ##
            trainSetData <- genoMatrix_trainLoci[-fold_index[[i]],]#Get the training set data (training individuals/loci)
            testSetData <- genoMatrix_trainLoci[fold_index[[i]],]#Get the test set data (test individuals/loci)
            ##
            #Peform PCA on training
            PCA_results <- perform.PCA(trainSetData[,1:ncol(trainSetData)-1], method=pca.method) #Run PCA without label column
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
              svm.fit <- svm(popNames_vector ~ ., data=trainSetData_PC, kernel=svm.kernel, cost=svm.cost, prob=T, scale=T)
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
              tempAlleleIndex <- grep(pat=paste0(trainLocusName[m],"_"), alleleName)
              trainLocusIndex_genoMatrix <- c(trainLocusIndex_genoMatrix,tempAlleleIndex )
            }
            trainLocusIndex_genoMatrix <- sort(trainLocusIndex_genoMatrix)
            genoMatrix_trainLoci <- genoMatrix[,c(trainLocusIndex_genoMatrix, ncol(genoMatrix))]#resahpe the genoMatrix data, comprising only training loci
            ##
            trainSetData <- genoMatrix_trainLoci[-fold_index[[i]],]#Get the training set data (training individuals/loci)
            testSetData <- genoMatrix_trainLoci[fold_index[[i]],]#Get the test set data (test individuals/loci)
            ##
            #Peform PCA on training
            PCA_results <- perform.PCA(trainSetData[,1:ncol(trainSetData)-1], method=pca.method) #Run PCA without label column
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
              svm.fit <- svm(popNames_vector ~ ., data=trainSetData_PC, kernel=svm.kernel, cost=svm.cost, prob=T, scale=T)
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
              tempAlleleIndex <- grep(pat=paste0(trainLocusName[m],"_"), alleleName)#alleleName is colnames(genoMatrix)
              trainLocusIndex_genoMatrix <- c(trainLocusIndex_genoMatrix,tempAlleleIndex )
            }
            trainLocusIndex_genoMatrix <- sort(trainLocusIndex_genoMatrix)
            genoMatrix_trainVars <- genoMatrix[,c(trainLocusIndex_genoMatrix,firstOtherVars_pos:lastOtherVars_pos,ncol(genoMatrix))]#resahpe the genoMatrix data, comprising training loci + non-genetic + popNames
            ##Split entire data into training and test sets
            trainSetData <- genoMatrix_trainVars[-fold_index[[i]],]#Get the training set data (training individuals/loci)
            testSetData <- genoMatrix_trainVars[fold_index[[i]],]#Get the test set data (test individuals/loci)
            ##
            #Peform PCA on training
            PCA_results <- perform.PCA(trainSetData[,1:ncol(trainSetData)-1], method=pca.method) #Run PCA without label column
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
              svm.fit <- svm(popNames_vector ~ ., data=trainSetData_PC, kernel=svm.kernel, cost=svm.cost, prob=T, scale=T)
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
              tempAlleleIndex <- grep(pat=paste0(trainLocusName[m],"_"), alleleName)
              trainLocusIndex_genoMatrix <- c(trainLocusIndex_genoMatrix,tempAlleleIndex )
            }
            trainLocusIndex_genoMatrix <- sort(trainLocusIndex_genoMatrix)
            genoMatrix_trainVars <- genoMatrix[,c(trainLocusIndex_genoMatrix,firstOtherVars_pos:lastOtherVars_pos,ncol(genoMatrix))]#resahpe the genoMatrix data, comprising training loci + non-genetic + popNames
            ##Split entire data into training and test sets
            trainSetData <- genoMatrix_trainVars[-fold_index[[i]],]#Get the training set data (training individuals/loci)
            testSetData <- genoMatrix_trainVars[fold_index[[i]],]#Get the test set data (test individuals/loci)
            ##
            #Peform PCA on training
            PCA_results <- perform.PCA(trainSetData[,1:ncol(trainSetData)-1], method=pca.method) #Run PCA without label column
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
              svm.fit <- svm(popNames_vector ~ ., data=trainSetData_PC, kernel=svm.kernel, cost=svm.cost, prob=T, scale=T)
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
    }#for(k in k.fold)
  }#else if(length(x)==5)
  stopCluster(cl)
  #
  #Output a metadata file
  cat(" Analysis Description (R - assignPOP ver.1.0)\n",
      "Perform assign.kfold() @", format(Sys.time()),"\n\n",
      "k.fold =",k.fold,"\n",
      "train.loci =",train.loci,"\n",
      "Total assignment tests =",sum(k.fold*length(train.loci)),"\n",
      "fst locus sample method:",loci.sample,"\n",
      "PC retaining criteria:",pca.method,"\n",
      "Machine-learning model:",model,"\n\n",
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

}# End
