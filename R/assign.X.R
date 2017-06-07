#' Perform a population assignment test on unknown individuals using known data
#' 
#' This function assigns unknown individuals to possible source populations based on known individuals and genetic or non-genetic or integrated data.
#' @param x1 An input object containing data from known individuals for building predictive models. It could be a list object returned from the function read.genpop(), reduce.allele() or compile.data(). Or, it could be a data frame containing non-genetic data returned from read.csv() or read.table().
#' @param x2 An input object containing data from unknown individuals to be predicted. It could be a list object returned from read.genpop(), reduce.allele(), or compile.data(). Or, it could be a data frame containing non-genetic data returned from read.csv() or read.table(). The x1 and x2 should be the same type (both are either lists or data frames). 
#' @param dir A character string to specify the folder name for saving output files. A slash at the end must be included (e.g., dir="YourFolderName/"). Otherwise, the files will be saved under your working directory.
#' @param scaled A logical variable (TRUE or FALSE) to specify whether to center (make mean of each feature to 0) and scale (make standard deviation of each feature to 1) the dataset before performing PCA and cross-validation. Default is FALSE. As genetic data has converted to numeric data between 0 and 1, to scale or not to scale the genetic data should not be critical. However, it is recommended to set scaled=TRUE when integrated data contains various scales of features.  
#' @param pca.method Either a character string ("mixed", "independent", or "original") or logical variable (TRUE or FALSE) to specify how to perform PCA on non-genetic data (PCA is always performed on genetic data). The character strings are used when analyzing integrated (genetic plus non-genetic) data. If using "mixed" (default), PCA is perfromed across the genetic and non-genetic data, resulting in each PC summarizing mixed variations of genetic and non-genetic data. If using "independent", PCA is independently performed on non-genetic data. Genetic PCs and non-genetic PCs are then used as new features. If using "original", original non-genetic data and genetic PCs are used as features. The logical variable is used when analyzing non-genetic data.If TRUE, it performs PCA on the training data and applys the loadings to the test data. Scores of training and test data will be used as new features. 
#' @param pca.PCs A criterion to retain number of PCs. By default, it uses Kaiser-Guttman criterion that any PC has the eigenvalue greater than 1 will be retained as the new variable/feature. Users can set an integer to specify the number of PCs to be retained.
#' @param pca.loadings A logical variable (TRUE or FALSE) to determine whether to output the loadings of training data to text files. Default is FALSE. Just a heads-up, the output files could take some storage space, if set TRUE.
#' @param model A character string to specify which classifier to use for creating predictive models. The current options include "lda", "svm", "naiveBayes", "tree", and "randomForest". Default is "svm"(support vector machine).
#' @param svm.kernel A character string to specify which kernel to be used when using "svm" classifier.
#' @param svm.cost A number to specify the cost for "svm" method.
#' @param ntree A integer to specify how many trees to build when using "randomForest" method.
#' @param mplot A logical variable to specify whether making a membership probability plot right after the assignment test is done. Set TRUE to make the plot. Otherwise, it will prompt a question.
#' @param ... Other arguments that could be potentially used for various models 
#' @return This function outputs assignment results and other analytical information in text files that will be saved under your designated folder. It also outputs a membership probability plot, if permitted.
#' @import stringr
#' @import foreach
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom caret createFolds
#' @importFrom MASS lda
#' @importFrom e1071 svm naiveBayes
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom tree tree
#' @importFrom randomForest randomForest importance
#' @importFrom utils packageVersion
#' @export
#' 
assign.X <- function(x1, x2, dir=NULL, scaled=F, pca.method="mixed", pca.PCs="kaiser-guttman", pca.loadings=F, 
                     model="svm", svm.kernel="linear", svm.cost=1, ntree=50, mplot=FALSE, ...){
  #check if x1 and x2 are the same type
  if(!class(x1)==class(x2)){
    stop("Input data sets are not the same type. Enter '?assign.X' to see description of x1 and x2 arguments")
  }
  #check if dir is correctly entered
  if(is.null(dir)){
    stop("Please provide a folder name ending with '/' in argument 'dir' ")
  }else if(substr(dir, start=nchar(dir), stop=nchar(dir))!="/"){
    stop("Please put a forward slash '/' in the end of your folder name (in argument 'dir'). ")
  }
  ##check data type
  if(!is.data.frame(x1)){#check if input x is a list returned from read.genpop(), reduce.allele(), or compile.data()
    #Analyze genetic or integrated data
    #checking pca.method
    if(!is.character(pca.method)){ #if pca.method is not character string, print message and stop analyzing
      stop("Please specify a correct parameter, 'mixed' or 'independent' or 'original', for argument 'pca.method' ")
    }
    #claim variables
    trainMatrix <- x1[[1]] #get training dataset
    trainVarName <- colnames(trainMatrix)[1:ncol(trainMatrix)-1] 
    testMatrix <- x2[[1]] #get test dataset
    testVarName <- colnames(testMatrix)[1:ncol(testMatrix)-1]
    testIndID <- x2[[2]]
    popSizes <- table(trainMatrix$popNames_vector)#get number of individual for each pop in table
    pops <- names(popSizes)#Get what pops in data
    noPops <- length(popSizes)#count number of pops
    #Create a folder to save outfiles
    dir.create(file.path(dir))
    #check if features are identical between x1 and x2 datasets
    if(isTRUE(all.equal(trainVarName, testVarName))){
      cat("\n  Known and unknown datasets have identical features. ")
      common_VarName <- trainVarName
      trainDataSet <- trainMatrix
      testDataSet <- testMatrix[,1:ncol(testMatrix)-1]
    }else {
      cat("\n  Known and unknown datasets have unequal features.")
      ans0 <- readline(" Continue and use common features for assignment? (enter Y/N): ")
      if(grepl(pattern="N", toupper(ans0))){
        warning("Program stops. Please revise your feature names accordingly.")
        on.exit()
      }else if(grepl(pattern="Y",toupper(ans0))){
        cat("\n  Automatically identify common features between datasets...")
        #Identify common features and modify datasets
        common_VarName <- intersect(trainVarName, testVarName)#Identify common feature names
        keepTrainVars <- c(common_VarName, "popNames_vector")
        trainDataSet <- trainMatrix[keepTrainVars]
        testDataSet <- testMatrix[common_VarName]
        cat("\n  ",length(common_VarName),"features are used for assignment.")
      }
    } 
    #center and scale the training data if scaled=T
    if(scaled){
      cat("\n  Scaling and centering data sets...")
      trainDataSet <- scale(trainDataSet[,1:ncol(trainDataSet)-1]) #Remove last column (popNames_vector) and scale/center data
      centering <- attr(trainDataSet,"scaled:center") #Get the mean of each feature
      scaling <- attr(trainDataSet,"scaled:scale") #Get the standard deviation of each feature
      trainDataSet <- as.data.frame(trainDataSet) #Convert matrix to data frame
      trainDataSet <- cbind(trainDataSet, trainMatrix$popNames_vector);colnames(trainDataSet)[ncol(trainDataSet)] <- "popNames_vector"
      #apply center(mean) and scale(sd) to test data
      testDataSet <- sweep(testDataSet, 2, centering, FUN="-") #substract mean by column 
      testDataSet <- sweep(testDataSet, 2, scaling, FUN="/") #divide sd by column
    } 
    #Check whether the dataset is genetic-only or integrated data
    if(length(x1)==3){ #Process genetic-only data
      datatype <- "genetics";noVars <- 0; pca.method <- "NA"
      #Peform PCA 
      cat("\n  Performing PCA on genetic data for dimensionality reduction...")
      PCA_results <- perform.PCA(trainDataSet[,1:ncol(trainDataSet)-1], method=pca.PCs) #Run PCA without label column
      loadings <- PCA_results[[1]] #loadings (coefficience) of variables and PCs; apply this to test data
      trainDataSet_PC <- as.data.frame(PCA_results[[2]])
      trainDataSet_PC <- cbind(trainDataSet_PC, trainDataSet$popNames_vector) ##Will be used for building predictive models
      colnames(trainDataSet_PC)[ncol(trainDataSet_PC)] <- "popNames_vector"
      #Convert test data to PC variables based on training's loadings
      testDataSet_matrix <- as.matrix(testDataSet)
      testDataSet_PC <- as.data.frame(testDataSet_matrix %*% loadings)
      #
      #Peform assignment using one of the following models
      if(model=="svm"){
        svm.fit <- svm(popNames_vector ~ ., data=trainDataSet_PC, kernel=svm.kernel, cost=svm.cost, prob=T, scale=F)
        svm.pred <- predict(svm.fit, testDataSet_PC, type="class",prob=T)
        outcome_matrix <- cbind(testIndID, as.data.frame(svm.pred), attr(svm.pred,"probabilities"))#combine output to data frame
        colnames(outcome_matrix)[1:2] <- c("Ind.ID","pred.pop")
        ##Output assignment results to files
        write.table(outcome_matrix, file=paste0(dir,"AssignmentResult.txt"), quote=F, row.names=F )
        cat(common_VarName, file=paste0(dir,"UsedLoci.txt"), sep="\n")
        if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"PC_Loadings.txt"), quote=F)}
        #
      }else if(model=="lda"){
        lda.fit <- lda(popNames_vector ~ ., data=trainDataSet_PC)
        lda.pred <- predict(lda.fit, testDataSet_PC)
        lda.pred.class <- lda.pred$class
        lda.pred.prob <- lda.pred$posterior
        outcome_matrix <- cbind(testIndID, as.data.frame(lda.pred.class), as.data.frame(lda.pred.prob))
        colnames(outcome_matrix)[1:2] <- c("Ind.ID","pred.pop")
        ##Output assignment results to files
        write.table(outcome_matrix, file=paste0(dir,"AssignmentResult.txt"), quote=F, row.names=F )
        cat(common_VarName, file=paste0(dir,"UsedLoci.txt"), sep="\n")
        if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"PC_Loadings.txt"), quote=F)}
        #
      }else if(model=="naiveBayes"){
        nby.model <- naiveBayes(popNames_vector ~ ., data=trainDataSet_PC)
        nby.pred.class <- predict(nby.model,testDataSet_PC,type="class")
        nby.pred.prob <- predict(nby.model,testDataSet_PC,type="raw")
        outcome_matrix <- cbind(testIndID, as.data.frame(nby.pred.class), as.data.frame(nby.pred.prob))
        colnames(outcome_matrix)[1:2] <- c("Ind.ID","pred.pop")
        ##Output assignment results to files
        write.table(outcome_matrix, file=paste0(dir,"AssignmentResult.txt"), quote=F, row.names=F )
        cat(common_VarName, file=paste0(dir,"UsedLoci.txt"), sep="\n")
        if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"PC_Loadings.txt"), quote=F)}
        #
      }else if(model=="tree"){
        tree.model <- tree(popNames_vector ~ ., data=trainDataSet_PC)
        tree.pred.class <- predict(tree.model,testDataSet_PC,type="class")
        tree.pred.prob <- predict(tree.model,testDataSet_PC,type="vector")
        outcome_matrix <- cbind(testIndID, as.data.frame(tree.pred.class), as.data.frame(tree.pred.prob))
        colnames(outcome_matrix)[1:2] <- c("Ind.ID","pred.pop")
        tree_node <- as.character(summary(tree.model)$used)#Get loci used at tree node
        ##Output assignment results to files
        write.table(outcome_matrix, file=paste0(dir,"AssignmentResult.txt"), quote=F, row.names=F )
        cat(common_VarName, file=paste0(dir,"UsedLoci.txt"), sep="\n")
        cat("Features used at tree node",tree_node,file=paste0(dir,"Feature_treenode.txt"), sep="\n")
        if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"PC_Loadings.txt"), quote=F)}
        #
      }else if(model=="randomForest"){
        rf.model <- randomForest(popNames_vector ~ ., data=trainDataSet_PC, ntree=ntree, importance=T)
        rf.pred.class <- predict(rf.model,testDataSet_PC,type="response")
        rf.pred.prob <- predict(rf.model,testDataSet_PC,type="prob")
        outcome_matrix <- cbind(testIndID, as.data.frame(rf.pred.class), as.data.frame(rf.pred.prob))
        colnames(outcome_matrix)[1:2] <- c("Ind.ID","pred.pop")
        ##Output assignment results to files
        write.table(outcome_matrix, file=paste0(dir,"AssignmentResult.txt"), quote=F, row.names=F )
        cat(common_VarName, file=paste0(dir,"UsedLoci.txt"), sep="\n")
        #Output "importance"(Gini index) of loci to files;see randomForest package's "importance" function
        write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Feature_importance.txt"), quote=F, row.names=T)
        if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"PC_Loadings.txt"), quote=F)}
        #
      }
      
    }else if(length(x1)==5){ #Process integrated data
      datatype <- "genetics + non-genetics"
      otherVarName <- x1[[4]]#Get non-genetic variable name
      #Identify the column position of non-genetic variable in dataset
      temp_vector <- NULL
      for(var in otherVarName){
        pos <- grep(pattern=var, common_VarName)
        temp_vector <- c(temp_vector, pos)
      }
      otherVar_startPos <- min(temp_vector) #Get the first column position of non-genetic data
      
      #Determine how to perform PCA on non-genetic data
      if(pca.method=="mixed"){
        cat("\n  Performing PCA on concatenated genetic and non-genetic data...")
        PCA_results <- perform.PCA(trainDataSet[,1:ncol(trainDataSet)-1], method=pca.PCs) #Run PCA without label column
        loadings <- PCA_results[[1]] #loadings (coefficience) of variables and PCs; apply this to test data
        trainDataSet_PC <- as.data.frame(PCA_results[[2]])
        trainDataSet_PC <- cbind(trainDataSet_PC, trainDataSet$popNames_vector) ##Will be used for building predictive models
        colnames(trainDataSet_PC)[ncol(trainDataSet_PC)] <- "popNames_vector"
        #Convert test data to PC variables based on training's loadings
        testDataSet_matrix <- as.matrix(testDataSet)
        testDataSet_PC <- as.data.frame(testDataSet_matrix %*% loadings)
        #
      }else if(pca.method=="independent"){
        #Perform PCA on training genetic and non-genetic independently
        cat("\n  Performing PCA on genetic and non-genetic data independently...")
        PCA_result_genetics <- perform.PCA(trainDataSet[,1:otherVar_startPos-1], method=pca.PCs)#Perform PCA on genetic data
        PCA_result_nongenetics <- perform.PCA(trainDataSet[,otherVar_startPos:(ncol(trainDataSet)-1)], method=pca.PCs)#perform PCA on non-genetic data
        loadings_genetics <- PCA_result_genetics[[1]] #loadings of genetic PCs; apply this to test data
        loadings_nongenetics <- PCA_result_nongenetics[[1]] #loadings of non-genetic PCs; apply this to test data
        trainDataSet_genetic_PC <- as.data.frame(PCA_result_genetics[[2]]); colnames(trainDataSet_genetic_PC) <- sub("PC", "genPC", colnames(trainDataSet_genetic_PC))
        trainDataSet_nongenetic_PC <- as.data.frame(PCA_result_nongenetics[[2]]);colnames(trainDataSet_nongenetic_PC) <- sub("PC", "nPC", colnames(trainDataSet_nongenetic_PC))
        trainDataSet_PC <- cbind(trainDataSet_genetic_PC, trainDataSet_nongenetic_PC, trainDataSet$popNames_vector)#Will use for building predictive models 
        colnames(trainDataSet_PC)[ncol(trainDataSet_PC)] <- "popNames_vector"
        #Convert test data to PC variables based on training's loadings
        testDataSet_genetic_matrix <- as.matrix(testDataSet[,1:otherVar_startPos-1]) #make genetic data matrix
        testDataSet_genetic_PC <- as.data.frame(testDataSet_genetic_matrix %*% loadings_genetics);colnames(testDataSet_genetic_PC)<-sub("PC", "genPC", colnames(testDataSet_genetic_PC))  
        testDataSet_nongenetic_matrix <- as.matrix(testDataSet[,otherVar_startPos:ncol(testDataSet)])
        testDataSet_nongenetic_PC <- as.data.frame(testDataSet_nongenetic_matrix %*% loadings_nongenetics);colnames(testDataSet_nongenetic_PC) <- sub("PC","nPC",colnames(testDataSet_nongenetic_PC))
        testDataSet_PC <- cbind(testDataSet_genetic_PC, testDataSet_nongenetic_PC )
        #
      }else if(pca.method=="original"){
        #Perform PCA on only genetic data
        cat("\n  Performing PCA on only genetic data...")
        PCA_result_genetics <- perform.PCA(trainDataSet[,1:otherVar_startPos-1], method=pca.PCs)#Perform PCA on genetic data
        loadings_genetics <- PCA_result_genetics[[1]] #loadings of genetic PCs; apply this to test data
        trainDataSet_genetic_PC <- as.data.frame(PCA_result_genetics[[2]]); colnames(trainDataSet_genetic_PC) <- sub("PC", "genPC", colnames(trainDataSet_genetic_PC))
        #concatenate genetic PCs and original non-genetic data and popNames_vector; note that dataset has PCs only from genetics
        trainDataSet_PC <- cbind(trainDataSet_genetic_PC, trainDataSet[,otherVar_startPos:ncol(trainDataSet)])
        colnames(trainDataSet_PC)[ncol(trainDataSet_PC)] <- "popNames_vector"
        #Convert test data (genetic part) to PC variables (scores) based on training 
        testDataSet_genetic_matrix <- as.matrix(testDataSet[,1:otherVar_startPos-1]) #make genetic data matrix
        testDataSet_genetic_PC <- as.data.frame(testDataSet_genetic_matrix %*% loadings_genetics);colnames(testDataSet_genetic_PC)<-sub("PC", "genPC", colnames(testDataSet_genetic_PC)) 
        testDataSet_PC <- cbind(testDataSet_genetic_PC, testDataSet[,otherVar_startPos:ncol(testDataSet)])
        #
      }
      ##Perform assignment test using one of the following classifiers
      if(model=="svm"){
        svm.fit <- svm(popNames_vector ~ ., data=trainDataSet_PC, kernel=svm.kernel, cost=svm.cost, prob=T, scale=F)
        svm.pred <- predict(svm.fit, testDataSet_PC, type="class",prob=T)
        outcome_matrix <- cbind(testIndID, as.data.frame(svm.pred), attr(svm.pred,"probabilities"))#combine output to data frame
        colnames(outcome_matrix)[1:2] <- c("Ind.ID","pred.pop")
        ##Output assignment results to files
        write.table(outcome_matrix, file=paste0(dir,"AssignmentResult.txt"), quote=F, row.names=F )
        cat(common_VarName, file=paste0(dir,"UsedFeatures.txt"), sep="\n")
        if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"PC_Loadings.txt"), quote=F)}
        #
      }else if(model=="lda"){
        lda.fit <- lda(popNames_vector ~ ., data=trainDataSet_PC)
        lda.pred <- predict(lda.fit, testDataSet_PC)
        lda.pred.class <- lda.pred$class
        lda.pred.prob <- lda.pred$posterior
        outcome_matrix <- cbind(testIndID, as.data.frame(lda.pred.class), as.data.frame(lda.pred.prob))
        colnames(outcome_matrix)[1:2] <- c("Ind.ID","pred.pop")
        ##Output assignment results to files
        write.table(outcome_matrix, file=paste0(dir,"AssignmentResult.txt"), quote=F, row.names=F )
        cat(common_VarName, file=paste0(dir,"UsedFeatures.txt"), sep="\n")
        if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"PC_Loadings.txt"), quote=F)}
        #
      }else if(model=="naiveBayes"){
        nby.model <- naiveBayes(popNames_vector ~ ., data=trainDataSet_PC)
        nby.pred.class <- predict(nby.model,testDataSet_PC,type="class")
        nby.pred.prob <- predict(nby.model,testDataSet_PC,type="raw")
        outcome_matrix <- cbind(testIndID, as.data.frame(nby.pred.class), as.data.frame(nby.pred.prob))
        colnames(outcome_matrix)[1:2] <- c("Ind.ID","pred.pop")
        ##Output assignment results to files
        write.table(outcome_matrix, file=paste0(dir,"AssignmentResult.txt"), quote=F, row.names=F )
        cat(common_VarName, file=paste0(dir,"UsedFeatures.txt"), sep="\n")
        if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"PC_Loadings.txt"), quote=F)}
        #
      }else if(model=="tree"){
        tree.model <- tree(popNames_vector ~ ., data=trainDataSet_PC)
        tree.pred.class <- predict(tree.model,testDataSet_PC,type="class")
        tree.pred.prob <- predict(tree.model,testDataSet_PC,type="vector")
        outcome_matrix <- cbind(testIndID, as.data.frame(tree.pred.class), as.data.frame(tree.pred.prob))
        colnames(outcome_matrix)[1:2] <- c("Ind.ID","pred.pop")
        tree_node <- as.character(summary(tree.model)$used)#Get loci used at tree node
        ##Output assignment results to files
        write.table(outcome_matrix, file=paste0(dir,"AssignmentResult.txt"), quote=F, row.names=F )
        cat(common_VarName, file=paste0(dir,"UsedFeatures.txt"), sep="\n")
        cat("Features used at tree node",tree_node,file=paste0(dir,"Feature_treenode.txt"), sep="\n")
        if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"PC_Loadings.txt"), quote=F)}
        #
      }else if(model=="randomForest"){
        rf.model <- randomForest(popNames_vector ~ ., data=trainDataSet_PC, ntree=ntree, importance=T)
        rf.pred.class <- predict(rf.model,testDataSet_PC,type="response")
        rf.pred.prob <- predict(rf.model,testDataSet_PC,type="prob")
        outcome_matrix <- cbind(testIndID, as.data.frame(rf.pred.class), as.data.frame(rf.pred.prob))
        colnames(outcome_matrix)[1:2] <- c("Ind.ID","pred.pop")
        ##Output assignment results to files
        write.table(outcome_matrix, file=paste0(dir,"AssignmentResult.txt"), quote=F, row.names=F )
        cat(common_VarName, file=paste0(dir,"UsedFeatures.txt"), sep="\n")
        #Output "importance"(Gini index) of loci to files;see randomForest package's "importance" function
        write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Feature_importance.txt"), quote=F, row.names=T)
        if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"PC_Loadings.txt"), quote=F)}
        #
      }
    }#else if(length(x1)==5)
    #Count number of unknown inds assigned to pops
    res_popSizes <- table(outcome_matrix$pred.pop)
    
    #Output a metadata file
    version <- as.character(packageVersion("assignPOP"))
    cat(" Analysis Description ( R - assignPOP ver.",version,")\n",
        "Perform assign.X() @", format(Sys.time()),"\n\n",
        "Data scaled and centerd:",scaled,"\n",
        "PC retaining criteria:",pca.PCs,"\n",
        "PCA for non-genetic data:",pca.method,"\n",
        "Machine learning model:",model,"\n\n",
        "Input Data (",datatype,")\n",
        "Number of known individuals:",sum(popSizes),"\n",
        "Number of known populations", noPops,"\n",
        names(popSizes),"\n",popSizes,"\n\n",
        "Number of unknown individuals:",nrow(testDataSet),"\n",
        "Number of unknown individuals assigned to populations:\n",
        names(res_popSizes),"\n",res_popSizes,"\n",
        file=paste0(dir,"AnalysisInfo.txt"))
    #Print some message to R console
    cat("\n  Assignment test is done! See results in your designated folder.")
    cat("\n  Predicted populations and probabilities are saved in [AssignmentResult.txt]")
    #Make a membership probability plot if answer is "Y"
    if(mplot){
      ans1 <- "Y"
    }else{
      ans1 <- readline("  Do you want to make a membership probability plot now? (enter Y/N): ")
    }
    if(grepl(pattern ="N", toupper(ans1))){
      on.exit()
    }else if(grepl(pattern ="Y", toupper(ans1))){
      value <- NULL; variable <- NULL; Ind.ID <- NULL
      ndf <- melt(outcome_matrix, id.vars=c("Ind.ID","pred.pop")) #Reshape the data, making probabilities in one single column (var name="value")
      stackplot <- ggplot(ndf, aes(x=Ind.ID, y=value, fill=variable))+
        geom_bar(stat="identity", width=1)+ # width=1 allows no space between bars
        #scale_fill_grey()+ # Make the bar color in grey scale
        ylab("Probability")+
        guides(fill=guide_legend(title=NULL))+ #Hiding title of legend
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),#hiding grid of the panel
              strip.background = element_rect(colour="black", fill="white", linetype="solid"),#change facet title background color
              plot.title = element_text(size=16, vjust=0.8),
              legend.text = element_text(size=14),
              strip.text.x = element_text(size=16),
              axis.title.y = element_text(size=16), axis.text.y = element_text(size=14, colour="black"),
              axis.title.x = element_blank(), axis.text.x = element_text(angle=90, size=7) )
      return(stackplot)
    }#else if(grepl(pattern ="Y", toupper(ans1)))
    
  }else if(is.data.frame(x1)){
    #Analyze non-genetic data
    datatype <- "non-genetic"
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
    #claim variables
    trainDataSet <- x1[,-1]; colnames(trainDataSet)[ncol(trainDataSet)] <- "popName" #Remove ID column and assign popName to pop column name
    testDataSet <- x2[,-1] #Remove ID column for unknown dataset
    trainVarNames <- colnames(trainDataSet) #This includes feature and target (popName) names
    testVarNames <- colnames(testDataSet) #This includes feature names
    testIndID <- x2[,1] #Get unknown individual IDs
    popName <- x1[,ncol(x1)]
    popSizes <- table(popName)
    noFactorVar <- 0 #Variable to count number of factor variables
    common_VarNames <- intersect(trainVarNames, testVarNames) #Identify common features in train and test datasets
    cat("\n  ",length(common_VarNames),"features are found and will be used for prediction.")
    #Create a folder to save outfiles
    dir.create(file.path(dir))
    #Check data type
    cat("\n  checking data type...")
    for(var in common_VarNames){
      varType <- class(x1[[var]]) #use double square brackets to turn one-column data frame to a vector 
      cat(paste0("  ",var,"(",varType,")"))
    }
    ans0 <- readline("  Are they correct? (enter Y/N): ")
    if(grepl(pattern="N", toupper(ans0))){
      cat("  please enter variable names for changing data type (separate names by a whitespace if multiple)\n")
      ans1 <- readline("  enter here: ")
      ans1 <- str_trim(ans1, side="both")
      ans1 <- unlist(strsplit(ans1,split=" "))#check out variable name to be processed
      noChangeVar <- length(ans1)
      #Check if entry is correct
      if(!all(ans1 %in% common_VarNames)){ #if any of entry not in common_VarNames is true 
        stop("Please enter correct feature names.")
      }
      #Process variables and convert factor data to dummy variable (binary data)
      for(name in ans1){
        ans2 <- readline(paste0("  Which data type should '",name,"' be? (enter numeric or factor): "))
        if(grepl(pattern="N",toupper(ans2))){
          trainDataSet[,name] <- as.numeric(as.character(trainDataSet[,name]))
          testDataSet[,name] <- as.numeric(as.character(testDataSet[,name]))
        }else if(grepl(pattern="F",toupper(ans2))){
          trainDataSet[,name] <- as.factor(trainDataSet[,name])
          testDataSet[,name] <- as.factor(testDataSet[,name])
        }
      }
      #Convert factor data to dummy data
      for(name in common_VarNames){
        if(is.factor(trainDataSet[,name])){
          noFactorVar <- noFactorVar + 1 #count number of categorical varibales
          #Convert factor variable to numeric binary variable (dummy variable)
           #Process for known dataset (x1)
          dummyData <- as.data.frame(model.matrix( ~ trainDataSet[,name]-1, data=trainDataSet))#get dummy variable data frame
          names(dummyData) <- substring(names(dummyData), 21, 1000L)#extract meaningful wording, or remove some funny wording
          names(dummyData) <- sub("\\b", paste0(name,"."), names(dummyData))#append original variabel name at the beginning
          trainDataSet[,name] <- NULL #remove original factor data column
          trainDataSet <- cbind(dummyData, trainDataSet) #column bind dummy data while making popName last column
           #Process for unknown dataset (x2)
          dummyDataU <- as.data.frame(model.matrix( ~ testDataSet[,name]-1, data=testDataSet))
          names(dummyDataU) <- substring(names(dummyDataU), 20, 1000L)
          names(dummyDataU) <- sub("\\b", paste0(name,"."), names(dummyDataU))
          testDataSet[,name] <- NULL
          testDataSet <- cbind(dummyDataU, testDataSet)
        }#if(is.factor(trainMatrix[,name]))
      }# for(name in common_VarNames)
    }else if(grepl(pattern="Y",toupper(ans0))){
      #check through data and covert factor to dummy
      for(name in common_VarNames){
        if(is.factor(trainDataSet[,name])){
          noFactorVar <- noFactorVar + 1 #count number of categorical varibales
          #Convert factor variable to numeric binary variable (dummy variable)
           #Process for known dataset (x1)
          dummyData <- as.data.frame(model.matrix( ~ trainDataSet[,name]-1, data=trainDataSet))#get dummy variable data frame
          names(dummyData) <- substring(names(dummyData), 21, 1000L)#extract meaningful wording, or remove some funny wording
          names(dummyData) <- sub("\\b", paste0(name,"."), names(dummyData))#append original variabel name at the beginning
          trainDataSet[,name] <- NULL #remove original factor data column
          trainDataSet <- cbind(dummyData, trainDataSet) #column bind dummy data while making popName last column
           #Process for unknown dataset (x2)
          dummyDataU <- as.data.frame(model.matrix( ~ testDataSet[,name]-1, data=testDataSet))
          names(dummyDataU) <- substring(names(dummyDataU), 20, 1000L)
          names(dummyDataU) <- sub("\\b", paste0(name,"."), names(dummyDataU))
          testDataSet[,name] <- NULL
          testDataSet <- cbind(dummyDataU, testDataSet)
        }#if(is.factor(x1[,name]))
      }#for(name in common_VarNames)
    }#else if(grepl(pattern="Y",toupper(ans0)))
    #Scale and center data set if scaled=T
    if(scaled){
      cat("\n  Scaling and centering data sets...")
      trainDataSet <- scale(trainDataSet[,1:ncol(trainDataSet)-1]) #Remove last column (pop target) and scale/center data
      centering <- attr(trainDataSet,"scaled:center") #Get the mean of each feature
      scaling <- attr(trainDataSet,"scaled:scale") #Get the standard deviation of each feature
      trainDataSet <- as.data.frame(trainDataSet) #Convert matrix to data frame
      trainDataSet <- cbind(trainDataSet, popName)
      #apply center(mean) and scale(sd) to test data
      testDataSet <- sweep(testDataSet, 2, centering, FUN="-") #substract mean by column 
      testDataSet <- sweep(testDataSet, 2, scaling, FUN="/") #divide sd by column
    }#if(scaled)
    #Determine if performing PCA on training for dimensionality reduction
    if(pca.method==T){
      cat("\n  Performing PCA on data sets...")
      PCA_results <- perform.PCA(trainDataSet[,1:ncol(trainDataSet)-1] ,method=pca.PCs)#Run PCA exclude last column
      loadings <- PCA_results[[1]] #loadings (coefficience) of variables and PCs; apply this to test data
      trainDataSet_PC <- as.data.frame(PCA_results[[2]])
      trainDataSet_PC <- cbind(trainDataSet_PC, popName) ##Will be used for building predicting models
      #Convert test data to PC variables based on training's loadings
      testDataSet_matrix <- as.matrix(testDataSet)
      testDataSet_PC <- as.data.frame(testDataSet_matrix %*% loadings)
    }else if(pca.method==F){
      trainDataSet_PC <- trainDataSet
      testDataSet_PC <- testDataSet
    }
    #
    #Perform assignment using one of the following models
    if(model=="svm"){
      svm.fit <- svm(popName ~ ., data=trainDataSet_PC, kernel=svm.kernel, cost=svm.cost, prob=T, scale=F)
      svm.pred <- predict(svm.fit, testDataSet_PC, type="class",prob=T)
      outcome_matrix <- cbind(testIndID, as.data.frame(svm.pred), attr(svm.pred,"probabilities"))#combine output to data frame
      colnames(outcome_matrix)[1:2] <- c("Ind.ID","pred.pop")
      ##Output assignment results to files
      write.table(outcome_matrix, file=paste0(dir,"AssignmentResult.txt"), quote=F, row.names=F )
      cat(common_VarNames, file=paste0(dir,"UsedFeatures.txt"), sep="\n")
      if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"PC_Loadings.txt"), quote=F)}
      #
    }else if(model=="lda"){
      lda.fit <- lda(popName ~ ., data=trainDataSet_PC)
      lda.pred <- predict(lda.fit, testDataSet_PC)
      lda.pred.class <- lda.pred$class
      lda.pred.prob <- lda.pred$posterior
      outcome_matrix <- cbind(testIndID, as.data.frame(lda.pred.class), as.data.frame(lda.pred.prob))
      colnames(outcome_matrix)[1:2] <- c("Ind.ID","pred.pop")
      ##Output assignment results to files
      write.table(outcome_matrix, file=paste0(dir,"AssignmentResult.txt"), quote=F, row.names=F )
      cat(common_VarNames, file=paste0(dir,"UsedFeatures.txt"), sep="\n")
      if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"PC_Loadings.txt"), quote=F)}
      #
    }else if(model=="naiveBayes"){
      nby.model <- naiveBayes(popName ~ ., data=trainDataSet_PC)
      nby.pred.class <- predict(nby.model,testDataSet_PC,type="class")
      nby.pred.prob <- predict(nby.model,testDataSet_PC,type="raw")
      outcome_matrix <- cbind(testIndID, as.data.frame(nby.pred.class), as.data.frame(nby.pred.prob))
      colnames(outcome_matrix)[1:2] <- c("Ind.ID","pred.pop")
      ##Output assignment results to files
      write.table(outcome_matrix, file=paste0(dir,"AssignmentResult.txt"), quote=F, row.names=F )
      cat(common_VarNames, file=paste0(dir,"UsedFeatures.txt"), sep="\n")
      if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"PC_Loadings.txt"), quote=F)}
      #
    }else if(model=="tree"){
      tree.model <- tree(popName ~ ., data=trainDataSet_PC)
      tree.pred.class <- predict(tree.model,testDataSet_PC,type="class")
      tree.pred.prob <- predict(tree.model,testDataSet_PC,type="vector")
      outcome_matrix <- cbind(testIndID, as.data.frame(tree.pred.class), as.data.frame(tree.pred.prob))
      colnames(outcome_matrix)[1:2] <- c("Ind.ID","pred.pop")
      tree_node <- as.character(summary(tree.model)$used)#Get loci used at tree node
      ##Output assignment results to files
      write.table(outcome_matrix, file=paste0(dir,"AssignmentResult.txt"), quote=F, row.names=F )
      cat(common_VarNames, file=paste0(dir,"UsedFeatures.txt"), sep="\n")
      cat("Features used at tree node",tree_node,file=paste0(dir,"Feature_treenode.txt"), sep="\n")
      if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"PC_Loadings.txt"), quote=F)}
      #
    }else if(model=="randomForest"){
      rf.model <- randomForest(popName ~ ., data=trainDataSet_PC, ntree=ntree, importance=T)
      rf.pred.class <- predict(rf.model,testDataSet_PC,type="response")
      rf.pred.prob <- predict(rf.model,testDataSet_PC,type="prob")
      outcome_matrix <- cbind(testIndID, as.data.frame(rf.pred.class), as.data.frame(rf.pred.prob))
      colnames(outcome_matrix)[1:2] <- c("Ind.ID","pred.pop")
      ##Output assignment results to files
      write.table(outcome_matrix, file=paste0(dir,"AssignmentResult.txt"), quote=F, row.names=F )
      cat(common_VarNames, file=paste0(dir,"UsedFeatures.txt"), sep="\n")
      #Output "importance"(Gini index) of loci to files;see randomForest package's "importance" function
      write.table(as.data.frame(importance(rf.model,type=2)),file=paste0(dir,"Feature_importance.txt"), quote=F, row.names=T)
      if(pca.loadings){ write.table(as.data.frame(loadings), file = paste0(dir,"PC_Loadings.txt"), quote=F)}
      #
    }
    #Count number of unknown inds assigned to pops
    res_popSizes <- table(outcome_matrix$pred.pop)
    #Output a metadata file
    version <- as.character(packageVersion("assignPOP"))
    cat(" Analysis Description ( R - assignPOP ver.",version,")\n",
        "Perform assign.X() @", format(Sys.time()),"\n\n",
        "Data scaled and centerd:",scaled,"\n",
        "PC retaining criteria:",pca.PCs,"\n",
        "PCA for non-genetic data:",pca.method,"\n",
        "Machine learning model:",model,"\n\n",
        "Input Data (",datatype,")\n",
        "Number of known individuals:",sum(popSizes),"\n",
        "Number of non-genetic variables:",length(common_VarNames),"\n",
        "Number of categorical variables:",noFactorVar,"\n",
        "Number of numeric variable:",length(common_VarNames)-noFactorVar,"\n",
        "Number of known populations", length(popSizes),"\n",
        names(popSizes),"\n",popSizes,"\n\n",
        "Number of unknown individuals:",nrow(testDataSet),"\n",
        "Number of unknown individuals assigned to populations:\n",
        names(res_popSizes),"\n",res_popSizes,"\n",
        file=paste0(dir,"AnalysisInfo.txt"))
    #Print some message to R console
    cat("\n  Assignment test is done! See results in your designated folder.")
    cat("\n  Predicted populations and probabilities are saved in [AssignmentResult.txt]")
    #Make a membership probability plot if answer is "Y"
    if(mplot){
      ans1 <- "Y"
    }else{
      ans1 <- readline("  Do you want to make a membership probability plot now? (enter Y/N): ")
    }
    if(grepl(pattern="N", toupper(ans1))){
      on.exit()
    }else if(grepl(pattern="Y", toupper(ans1))){
      value <- NULL; variable <- NULL; Ind.ID <- NULL
      ndf <- melt(outcome_matrix, id.vars=c("Ind.ID","pred.pop")) #Reshape the data, making probabilities in one single column (var name="value")
      stackplot <- ggplot(ndf, aes(x=Ind.ID, y=value, fill=variable))+
        geom_bar(stat="identity", width=1)+ # width=1 allows no space between bars
        #scale_fill_grey()+ # Make the bar color in grey scale
        ylab("Probability")+
        guides(fill=guide_legend(title=NULL))+ #Hiding title of legend
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),#hiding grid of the panel
              strip.background = element_rect(colour="black", fill="white", linetype="solid"),#change facet title background color
              plot.title = element_text(size=16, vjust=0.8),
              legend.text = element_text(size=14),
              strip.text.x = element_text(size=16),
              axis.title.y = element_text(size=16), axis.text.y = element_text(size=14, colour="black"),
              axis.title.x = element_blank(), axis.text.x = element_text(angle=90, size=7) )
      return(stackplot)
    }#else if(grepl(pattern="Y", toupper(ans1)))
  }#else if(is.data.frame(x1))
}#End