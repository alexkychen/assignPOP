#' Make a boxplot (ggplot2 style) of assignment accuracy from cross-validation results
#'
#' This functions allows you to make a boxplot of assignment accuracies estimated from Monte-Carlo or K-fold cross-validation results.
#' @param df A dataframe of your assignment accuracy results. It could be the object returned from the function accuracy.MC() or accuracy.kfold() or a data frame imported to R via other functions (e.g., read.table(...)).
#' @param pop Population names (one or multiple string characters) for making the plot. By default, it uses "all", meaning overall assignment accuracies. It creates faceted plot with one population per panel, if multiple population names are given. The specified population name should match what you entered in read.genpop() earlier.
#' @return This function returns a boxplot plot using the ggplot2 library. Users can modified the plot (e.g., change color, text, etc.) using functions provided by ggplot2 library.
#' @import ggplot2
#' @importFrom reshape2 melt
#' @examples
#' Your_df <- read.table(system.file("extdata/Rate.txt", package="assignPOP"), header=TRUE)
#' accuracy.plot(Your_df, pop="all")
#' @export
#'
accuracy.plot <- function(df, pop="all"){
  #claim variables
  train.inds <- NULL; train.loci <- NULL; value <- NULL; KF <- NULL 
  #validate specified pop names
  df_popName <- substring(colnames(df)[4:ncol(df)], 13, 1000L)
  if(!all(pop %in% df_popName)){ #if specified pop name not in df 
    stop(paste0("Pop name not found. Please use one or more of the following names [ ",toString(df_popName)," ] in argument 'pop'. "))
  }
  #check if assignment results of either Monte-Carlo or K-fold
  firstColname <- names(df)[1]
  if(firstColname=="train.inds"){
    #check if training inds is in proportion or fixed number
    if(is.factor(df$train.inds)){ #this is used when df is returned object from accuracy.MC() or accuracy.kfold()
      checkTrainInds <- as.numeric(levels(df$train.inds))[df$train.inds]
    }else if(is.numeric(df$train.inds)){ #this is used when df is read from read.table()
      checkTrainInds <- as.numeric(unique(df$train.inds))
    }
    #claim x label
    if(all(checkTrainInds > 1)){
      x_label <- "Number of individuals used in training set"
    }else if(all(checkTrainInds < 1)){
      x_label <- "Proportion of individuals used in training set"
    }
    #Convert training.inds & train.loci to factors
    df$train.inds <- factor(df$train.inds)
    df$train.loci <- factor(df$train.loci)
    #
    #check pop names (one or multiple)
    if(length(pop)==1){
      col <- paste0("assign.rate.",pop)
      #see if multiple levels of train loci used.(e.g.,10%, 20%...of loci)
      if(length(unique(df$train.loci)) > 1 ){
        boxplot <- ggplot(df, aes_string(y=col, x="train.inds", fill="train.loci"))+
          geom_boxplot()+
          xlab(x_label) + ylab("Assignment accuracy")+
          scale_fill_discrete(name="Prop. of\ntrain loci",guide=guide_legend(reverse=TRUE))+
          theme_bw()+
          theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())
                #strip.text.y = element_text(size=16, face="bold"),
                #legend.title = element_text(size=17),
                #legend.text = element_text(size=16),
                #axis.text=element_text(size=16, colour="black"),
                #axis.title.y=element_text(size=20, vjust=1.5),
                #axis.title.x=element_text(size=20, vjust=0.1))
        return(boxplot)
      #see if only one level of train loci used (e.g.,used all loci)
      }else if(length(unique(df$train.loci))==1){
        boxplot <- ggplot(df, aes_string(y=col, x="train.inds"))+
          geom_boxplot()+
          xlab(x_label) + ylab("Assignment accuracy")+
          theme_bw()+
          theme(legend.position="none",
                panel.grid.major = element_blank(), panel.grid.minor=element_blank())
                #axis.text=element_text(size=16, colour="black"),
                #axis.title.y=element_text(size=20, vjust=1.5),
                #axis.title.x=element_text(size=20, vjust=0.1))
        return(boxplot)
      }#else if(length(unique...))
    #If more than one pop names are specified, make faceted plots  
    }else if(length(pop)>1){
      #Get the first three column data (train.inds, train.loci, iter)
      dfn <- df[, 1:3]
      #Get the pop assign.rate column by arg. 'pop'
      for(pname in pop){
        coln <- paste0("assign.rate.",pname) 
        colnrate <- df[,coln]
        dfn <- cbind(dfn,colnrate)
        colnames(dfn)[ncol(dfn)] <- pname
      }
      #reshape data frame
      dfre <- melt(dfn, id=c(1,2,3)) #keep the first three column
      dfre$train.inds <- as.factor(dfre$train.inds)
      dfre$train.loci <- as.factor(dfre$train.loci)
      if("all" %in% pop){
        levels(dfre$variable) <- sub("all", "Overall", levels(dfre$variable)) #change "all" to "Overall" if exists
      }
      #check if train.loci has multiple levels
      if(length(unique(dfre$train.loci)) > 1){
        boxplot <- ggplot(dfre, aes(x=train.inds, y=value, fill=train.loci))+
          geom_boxplot()+
          facet_grid(. ~ variable)+
          xlab(x_label) + ylab("Assignment accuracy") +
          scale_fill_discrete(name="Prop. of\ntrain loci",guide=guide_legend(reverse=TRUE))+
          theme_bw()+
          theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())
                #strip.text.x = element_text(size=16, face="bold"),
                #legend.title = element_text(size=18),
                #legend.text = element_text(size=16),
                #axis.text=element_text(size=16, colour="black"),
                #axis.title.y=element_text(size=20, vjust=1.5),
                #axis.title.x=element_text(size=20, vjust=0.1))
        return(boxplot)
        
      #check if train.loci is one level
      }else if(length(unique(dfre$train.loci))==1){
        boxplot <- ggplot(dfre, aes(x=train.inds, y=value))+
          geom_boxplot()+
          facet_grid(. ~ variable)+
          xlab(x_label) + ylab("Assignment accuracy") +
          theme_bw()+
          theme(legend.position="none",
                panel.grid.major = element_blank(), panel.grid.minor=element_blank())
                #strip.text.x = element_text(size=16, face="bold"),
                #axis.text=element_text(size=16, colour="black"),
                #axis.title.y=element_text(size=20, vjust=1.5),
                #axis.title.x=element_text(size=20, vjust=0.1))
        return(boxplot)
        
      }#else if(length(unique(dfre$train.loci))==1)
      #
    }#else if(length(pop)>1)
    #
  }else if(firstColname=="KF"){
    #check pop names
    if(length(pop)==1){
      #Convert training.inds & train.loci to factors
      df$KF <- factor(df$KF)
      df$fold <- factor(df$fold)
      df$train.loci <- factor(df$train.loci)
      col <- paste0("assign.rate.",pop)
      #
      if(length(unique(df$train.loci)) > 1 ){ #see if multiple levels of train loci used.(e.g.,10%, 20%...of loci)
        boxplot <- ggplot(df, aes_string(y=col, x="KF" ,fill="train.loci"))+
          geom_boxplot()+
          #geom_point(size=5, position=dodge)+
          xlab("K") + ylab("Assignment accuracy")+
          scale_fill_discrete(name="Prop. of\ntrain loci",guide=guide_legend(reverse=TRUE))+
          theme_bw()+
          theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())
                #strip.text.y = element_text(size=16, face="bold"),
                #legend.title = element_text(size=17),
                #legend.text = element_text(size=16),
                #axis.text=element_text(size=16, colour="black"),
                #axis.title.y=element_text(size=20, vjust=1.5),
                #axis.title.x=element_text(size=20, vjust=0.1))
        return(boxplot)
        
      }else if(length(unique(df$train.loci))==1){ #see if only one level of train loci used (e.g.,used all loci)
        boxplot <- ggplot(df, aes_string(y=col, x="KF"))+
          geom_boxplot()+
          xlab("K") + ylab("Assignment accuracy")+
          theme_bw()+
          theme(legend.position="none",
                panel.grid.major = element_blank(), panel.grid.minor=element_blank())
                #axis.text=element_text(size=16, colour="black"),
                #axis.title.y=element_text(size=20, vjust=1.5),
                #axis.title.x=element_text(size=20, vjust=0.1))
        return(boxplot)
      } #else if(length(...
      
    }else if(length(pop)>1){
      #Get the first three column data (train.inds, train.loci, iter)
      dfn <- df[, 1:3]
      #Get the pop assign.rate column by arg. 'pop'
      for(pname in pop){
        coln <- paste0("assign.rate.",pname) 
        colnrate <- df[,coln]
        dfn <- cbind(dfn,colnrate)
        colnames(dfn)[ncol(dfn)] <- pname
      }
      #reshape data frame
      dfre <- melt(dfn, id=c(1,2,3)) #keep the first three column
      dfre$KF <- as.factor(dfre$KF)
      dfre$fold <- as.factor(dfre$fold)
      dfre$train.loci <- as.factor(dfre$train.loci)
      if("all" %in% pop){
        levels(dfre$variable) <- sub("all", "Overall", levels(dfre$variable)) #change "all" to "Overall" if exists
      }
      #check if train.loci have multiple levels
      if(length(unique(dfre$train.loci)) > 1 ){
        boxplot <- ggplot(dfre, aes(x=KF, y=value, fill=train.loci))+
          geom_boxplot()+
          facet_grid(. ~ variable)+
          xlab("K") + ylab("Assignment accuracy") +
          scale_fill_discrete(name="Prop. of\ntrain loci",guide=guide_legend(reverse=TRUE))+ #Reverse box order in legend
          theme_bw()+
          theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())
                #strip.text.x = element_text(size=16, face="bold"),
                #legend.title = element_text(size=18),
                #legend.text = element_text(size=16),
                #axis.text=element_text(size=16, colour="black"),
                #axis.title.y=element_text(size=20, vjust=1.5),
                #axis.title.x=element_text(size=20, vjust=0.1))
        return(boxplot)
      #check if train.loci is one level  
      }else if(length(unique(dfre$train.loci))==1 ){
        boxplot <- ggplot(dfre, aes(x=KF, y=value))+
          geom_boxplot()+
          facet_grid(. ~ variable)+
          xlab("K") + ylab("Assignment accuracy") +
          theme_bw()+
          theme(legend.position="none",
                panel.grid.major = element_blank(), panel.grid.minor=element_blank())
                #strip.text.x = element_text(size=16, face="bold"),
                #axis.text=element_text(size=16, colour="black"),
                #axis.title.y=element_text(size=20, vjust=1.5),
                #axis.title.x=element_text(size=20, vjust=0.1))
        return(boxplot)
      }
    }#else if(length(pop)>1)
    
  }#else if(firstColname=="KF")
  
}#End
