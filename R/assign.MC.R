#' Population assignment test using Monte-Carlo cross-validation
#'
#' This function employs Monte-Carlo cross-validation for assignment tests. It uses the object returned from the function read.genpop(), reduce.allele(), or compile.data() as input, and outputs results to text files. Typing a folder name in the argument "dir" is recommended. See below for more details.
#' @param x An input object which should be the object (list) returned from the function read.genpop(), reduce.allele(), or compile.data().
#' @param train.inds The number (integer greater than 1) or proportion (float between 0 and 1) of individuals (observations) from each population to be used as training data. Use a numeric vector to specify multiple sets of training individuals. No mixture of integer and float in a vector.
#' @param train.loci The proportion (float between 0 and 1) of loci to be used as training data. Use a numeric vector to specify multiple sets of training loci.
#' @param loci.sample Locus sampling method, "fst" or "random". If loci.sample="fst" (default) and train.loci=0.1, it means that top 10 percent of high Fst loci will be sampled as training loci. On the other hand, if loci.sample="random", then random 10 percent of loci will be sampled as training loci.
#' @param iterations Resampling times (an integer) for each combination of training individuals and loci.
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
assign.MC <- function(x, train.inds=c(0.5,0.7,0.9), train.loci=c(0.1,0.25,0.5, 1), loci.sample="fst", iterations=20, dir=NULL,
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
  if(length(x)==3){ #Process genetics-only data
    datatype <- "genetics";noVars <- 0
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
              tempAlleleIndex <- grep(pat=paste0(trainLocusName[m],"_"), alleleName)#alleleName is colnames(genoMatrix)
              trainLocusIndex_genoMatrix <- c(trainLocusIndex_genoMatrix,tempAlleleIndex )
            }
            trainLocusIndex_genoMatrix <- sort(trainLocusIndex_genoMatrix)
            genoMatrix_trainLoci <- genoMatrix[,c(trainLocusIndex_genoMatrix, ncol(genoMatrix))]#resahpe the genoMatrix data, comprising only training loci
            ##
            trainSetData <- genoMatrix_trainLoci[trainSet_index,]#Get the training set data (training individuals/loci)
            testSetData <- genoMatrix_trainLoci[-trainSet_index,]#Get the test set data (test individuals/loci)
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
              cat("Loci used at tree node",tree_node,file=paste0(dir,"Loci_treenode_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
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
              write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Loci_importance_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=T)
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
              tempAlleleIndex <- grep(pat=paste0(trainLocusName[m],"_"), alleleName)
              trainLocusIndex_genoMatrix <- c(trainLocusIndex_genoMatrix,tempAlleleIndex )
            }
            trainLocusIndex_genoMatrix <- sort(trainLocusIndex_genoMatrix)
            genoMatrix_trainLoci <- genoMatrix[,c(trainLocusIndex_genoMatrix, ncol(genoMatrix))]#resahpe the genoMatrix data, comprising only training loci
            ##
            trainSetData <- genoMatrix_trainLoci[trainSet_index,]#Get the training set data (training individuals/loci)
            testSetData <- genoMatrix_trainLoci[-trainSet_index,]#Get the test set data (test individuals/loci)
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
              cat("Loci used at tree node",tree_node,file=paste0(dir,"Loci_treenode_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
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
              write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Loci_importance_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=T)
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
              tempAlleleIndex <- grep(pat=paste0(trainLocusName[m],"_"), alleleName)#alleleName is colnames(genoMatrix)
              trainLocusIndex_genoMatrix <- c(trainLocusIndex_genoMatrix,tempAlleleIndex )
            }
            trainLocusIndex_genoMatrix <- sort(trainLocusIndex_genoMatrix)
            genoMatrix_trainVars <- genoMatrix[,c(trainLocusIndex_genoMatrix,firstOtherVars_pos:lastOtherVars_pos,ncol(genoMatrix))]#resahpe the genoMatrix data, comprising training loci + non-genetic + popNames
            ##Split entire data into training and test sets
            trainSetData <- genoMatrix_trainVars[trainSet_index,]#Get the training set data (training individuals/loci)
            testSetData <- genoMatrix_trainVars[-trainSet_index,]#Get the test set data (test individuals/loci)
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
              cat("Loci used at tree node",tree_node,file=paste0(dir,"Loci_treenode_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
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
              write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Loci_importance_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=T)
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
              tempAlleleIndex <- grep(pat=paste0(trainLocusName[m],"_"), alleleName)
              trainLocusIndex_genoMatrix <- c(trainLocusIndex_genoMatrix,tempAlleleIndex )
            }
            trainLocusIndex_genoMatrix <- sort(trainLocusIndex_genoMatrix)
            genoMatrix_trainVars <- genoMatrix[,c(trainLocusIndex_genoMatrix,firstOtherVars_pos:lastOtherVars_pos,ncol(genoMatrix))]#resahpe the genoMatrix data, comprising training loci + non-genetic + popNames
            ##Split entire data into training and test sets
            trainSetData <- genoMatrix_trainVars[trainSet_index,]#Get the training set data (training individuals/loci)
            testSetData <- genoMatrix_trainVars[-trainSet_index,]#Get the test set data (test individuals/loci)
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
              cat("Loci used at tree node",tree_node,file=paste0(dir,"Loci_treenode_",train.inds[j],"_",k,"_",i,".txt"), sep="\n")
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
              write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Loci_importance_",train.inds[j],"_",k,"_",i,".txt"), quote=F, row.names=T)
              if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"Loadings_",train.inds[j],"_",k,"_",i,".txt"), quote=F)}
              #
            }
          }#for(k in train.loci)
        }#else if(loci.sample=="random")
      }#for(j in 1:length(train.inds))
    }#foreach(i=1:iterations,...
  }#else if(length(x)==4)
  stopCluster(cl)
  #Output a metadata file
  cat(" Analysis Description (R - assignPOP ver.1.0)\n",
      "Perform assign.MC() @", format(Sys.time()),"\n\n",
      "train.inds =",train.inds,"\n",
      "train.loci =",train.loci,"\n",
      "iterations =",iterations,"\n",
      "Total assignment tests =",length(train.inds)*length(train.loci)*iterations,"\n",
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
  cat("\n  Monte-Carlo cross-validation done!!")
  cat(paste0("\n  ",length(train.inds)*length(train.loci)*iterations," assignment tests completed!!"))
  #Print a warning message if only one PC is used (when using Kaiser-Guttman criteria), meaning none of PC's eigenvalue greater than 1.
  #if(!is.null(PCA_results[[3]])){
  #  message(PCA_results[[3]])
  #}
}#END


### Function to estimate locus Fst (following Nei 1973)
Fsts <- function(x){
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
    getAlleleIndex <- grep(pat=getLocusName, genoNames)#Get allele index for extracing locus data from genoMatrix
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
  PCA_results <- prcomp(matrix, cor=TRUE)
  warnMessage <- NULL
  if(method=="kaiser-guttman"){ #Select PCs which eigenvalue greater than 1
    noPC_eig_g1 <- length(PCA_results$sdev[PCA_results$sdev > 1])#count how many PC's eigenvalue greater than 1
    if(noPC_eig_g1==0){
      warnMessage <- "None of PC's eigenvalue greater than 1. Only the first PC is used. Consider to assign a number to 'pca.method'."
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

