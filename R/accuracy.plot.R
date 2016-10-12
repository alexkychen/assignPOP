#' Make a boxplot (ggplot2 style) of assignment accuracy from cross-validation results
#'
#' This functions allows you to make a boxplot of assignment accuracies estimated from Monte-Carlo or K-fold cross-validation results.
#' @param df A dataframe of your assignment accuracy results. It could be the object returned from the function accuracy.MC() or accuracy.kfold() or a data frame imported to R via other functions (e.g., read.table(...)).
#' @param pop A population name for making the plot. By default, it uses "all", meaning overall assignment accuracies. If a population name is specified, it only creates the plot for that population. The specified population name should match what you entered in read.genpop() earlier.
#' @return This function returns a boxplot plot using the ggplot2 library. Users can modified the plot (e.g., change color, text, etc.) using functions provided by ggplot2 library.
#' @examples Your_df <- read.table("YourFolderName/Rate_of_N_tests_M_pops.txt", header=T)
#' accuracy.plot(Your_df, pop="all")
#' @export
#'
accuracy.plot <- function(df, pop="all"){
  #check if assignment results of either Monte-Carlo or K-fold
  firstColname <- names(df)[1]
  if(firstColname=="train.inds"){
    #check if training inds is in proportion or fixed number
    checkTrainInds <- as.numeric(levels(df$train.inds))[df$train.inds]
    if(all(checkTrainInds > 1)){
      x_label <- "Number of individuals used in training set"
    }else if(all(checkTrainInds < 1)){
      x_label <- "Proportion of individuals used in training set"
    }
    #Convert training.inds & train.loci to factors
    df$train.inds <- factor(df$train.inds)
    df$train.loci <- factor(df$train.loci)
    col <- paste0("assign.rate.",pop)

    if(length(unique(df$train.loci)) > 1 ){ #see if multiple levels of train loci used.(e.g.,10%, 20%...of loci)
      boxplot <- ggplot(df, aes_string(y=col, x="train.inds", fill="train.loci"))+
        geom_boxplot()+
        xlab(x_label) + ylab("Assignment accuracy")+
        scale_fill_discrete(guide = guide_legend(reverse=TRUE), #Reverse box order in legend
                            name="Prop. of\ntrain loci")+
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),
              #strip.text.y = element_text(size=16, face="bold"),
              legend.title = element_text(size=17),
              legend.text = element_text(size=16),
              axis.text=element_text(size=16, colour="black"),
              axis.title.y=element_text(size=20, vjust=1.5),
              axis.title.x=element_text(size=20, vjust=0.1))
      return(boxplot)

    }else if(length(unique(df$fst.level))==1){ #see if only one level of train loci used (e.g.,used all loci)
      boxplot <- ggplot(df, aes_string(y=col, x="train.inds"))+
        geom_boxplot()+
        xlab(x_label) + ylab("Assignment accuracy")+
        theme_bw()+
        theme(legend.position="none",
              panel.grid.major = element_blank(), panel.grid.minor=element_blank(),
              #strip.text.y = element_text(size=16, face="bold"),
              legend.title = element_text(size=17),
              legend.text = element_text(size=16),
              axis.text=element_text(size=16, colour="black"),
              axis.title.y=element_text(size=20, vjust=1.5),
              axis.title.x=element_text(size=20, vjust=0.1))
      return(boxplot)
    } #else if(length(...

  }else if(firstColname=="KF"){
    #Convert training.inds & train.loci to factors
    df$KF <- factor(df$KF)
    df$fold <- factor(df$fold)
    df$train.loci <- factor(df$train.loci)
    col <- paste0("assign.rate.",pop)
    dodge <- position_dodge(width=0.5)

    if(length(unique(df$train.loci)) > 1 ){ #see if multiple levels of train loci used.(e.g.,10%, 20%...of loci)
      boxplot <- ggplot(df, aes_string(y=col, x="KF" ,fill="train.loci"))+
        geom_boxplot()+
        #geom_point(size=5, position=dodge)+
        xlab("K") + ylab("Assignment accuracy")+
        scale_fill_discrete(guide = guide_legend(reverse=TRUE), #Reverse box order in legend
                            name="Prop. of\ntrain loci")+
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),
              #strip.text.y = element_text(size=16, face="bold"),
              legend.title = element_text(size=17),
              legend.text = element_text(size=16),
              axis.text=element_text(size=16, colour="black"),
              axis.title.y=element_text(size=20, vjust=1.5),
              axis.title.x=element_text(size=20, vjust=0.1))
      return(boxplot)

    }else if(length(unique(df$fst.level))==1){ #see if only one level of train loci used (e.g.,used all loci)
      boxplot <- ggplot(df, aes_string(y=col, x="KF"))+
        geom_boxplot()+
        xlab("K") + ylab("Assignment accuracy")+
        theme_bw()+
        theme(legend.position="none",
              panel.grid.major = element_blank(), panel.grid.minor=element_blank(),
              #strip.text.y = element_text(size=16, face="bold"),
              legend.title = element_text(size=17),
              legend.text = element_text(size=16),
              axis.text=element_text(size=16, colour="black"),
              axis.title.y=element_text(size=20, vjust=1.5),
              axis.title.x=element_text(size=20, vjust=0.1))
      return(boxplot)
    } #else if(length(...
  }#else if(firstColname=="KF")

}
