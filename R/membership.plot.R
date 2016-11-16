#' Make a membership probability plot using results from K-fold cross-validation (ggplot2 style)
#'
#' This function allows you to make a membership probability plot (stacked-bar plot) using results estimated from K-fold cross-validation.
#' @param dir A character string to specify the folder that has your K-fold cross-validation assignment results. A slash should be entered at the end.
#' @param style An option for output style. If style=1, it creates the plot which individuals on the x-axis are in random order. If style=2, individuals are sorted by probabilities within each population. If style=3, individuals of different folds are in seperate plots. If style=4, individuals are separated by fold and sorted by probability.
#' @return his function returns a stacked-bar plot using the ggplot2 library. Users can modified the plot (e.g., change color, text, etc.) using functions provided by ggplot2 library.
#' @examples # membership.plot(dir="YourFolderName/")
#' #If style is not specified, the program will ask you to enter after executing the function.
#' @import ggplot2
#' @export
membership.plot <- function(dir=NULL, style=NULL){
  Ind.ID <- NULL; value <- NULL; variable <- NULL #some NULL variable to handle R CMD check
  #Read all "Out_*" file names in a specified directory
  fileName_vec <- list.files(path=dir, pattern="Out_*")
  fileName_vec <- sort(fileName_vec)
  noFiles <- length(fileName_vec)#count number of files
  #Read one of files and get pop names
  result01 <- read.table(paste0(dir,fileName_vec[1]), header=T)
  pops <- names(result01)[4:length(names(result01))] #read column name and get the pop names between 4th to last column
  noPops <- length(pops)#Number of pops
  #create vectors for saving data
  k_fold_vec <- NULL
  train_loci_vec <- NULL
  #Read through file name and collect k.fold and train.loci
  for(i in 1:noFiles){
    oneFileName <- unlist(strsplit(fileName_vec[i], split="_")) #split file name to 4 elements (e.g.,"Out"  "0.1"  "K3"  "1.txt");"0.1"is train.loci,"K3"is k fold,"1"is fold
    train_loci_vec[i] <- oneFileName[2]
    k_fold_vec[i] <- gsub("K","",oneFileName[3])#Remove "k" character in the "K3" string
  }
  k.fold <- unique(k_fold_vec)#identify unique levels of K
  train.loci <- unique(train_loci_vec)#identify unique levels of training loci

  cat("\n  K = ");cat(paste0(k.fold," "));cat(" are found.") #print out detected "K" on console
  ans_k <- readline("  Please enter one of the K numbers: ") #ask user to enter one of the K numbers
  ans_k <- str_trim(ans_k, side="both") #clean any space
  if(!ans_k %in% k.fold){
    cat("\n  Entry is not correct.")
    break
  }

  cat(paste0("\n  ",length(train.loci)," levels of training loci are found."))
  cat("\n  Levels[train.loci]: ");cat(paste0(train.loci," "))
  ans_t <- readline("  Please enter one of the levels: ")
  ans_t <- str_trim(ans_t, side="both")
  if(!ans_t %in% train.loci){
    cat("\n  Entry is not correct.")
    break
  }

  #Read selected files
  df_mas <- data.frame(matrix(ncol=0,nrow=0))
  for(i in 1:ans_k){
    oneFile <- read.table(paste0(dir,"Out_",ans_t,"_K",ans_k,"_",i,".txt"), header=T)
    sampleSize <- nrow(oneFile)
    fold_n <- rep(paste0("fold_",i),sampleSize)
    oneFile <- cbind(oneFile, fold_n)
    df_mas <- rbind(df_mas, oneFile)
  }

  if(is.null(style)){
    cat("\n  Finally, select one of the output styles.")
    cat("\n  [1] Random order (Individuals on x-axis are in random order)")
    cat("\n  [2] Sorted by probability (Individuals are sorted by probabilities within each group)")
    cat("\n  [3] Separated by fold (Individuals of different folds are in separate plots)")
    cat("\n  [4] Separated and Sorted (Individuals are separated by fold and sorted by probability)")
    style <- readline("  Please enter 1, 2, 3, or 4: ")
    style <- str_trim(style, side="both")
    if(!style %in% c(1,2,3,4)){
      cat("\n  Entry is not correct.")
      break
    }
  }

  if(style==2){ # Individuals are sorted based on the probability of their own populations
    #Separate inds among pops and sort
    df_mas_2 <- data.frame(matrix(ncol=0,nrow=0))
    for(p in pops){
      df_pop <- df_mas[which(df_mas$origin.pop==p),]#Subset samples for each pop based on the origin.pop variable
      df_pop$Ind.ID <- as.factor(df_pop$Ind.ID) #convert ind id charater to factor data
      df_pop$Ind.ID <- droplevels(df_pop$Ind.ID) #Drop unexisting ind id levels
      df_pop$Ind.ID <- reorder(df_pop$Ind.ID, -df_pop[ ,p])#Reorder ind id levels based on the prob. of its pop
      df_mas_2 <- rbind(df_mas_2, df_pop)
    }
    ndf <- melt(df_mas_2, id.vars=c("Ind.ID","origin.pop","pred.pop","fold_n"))#Reshape the data, making probabilities in one single column (var name="value")
    stackplot <- ggplot(ndf, aes(x=Ind.ID, y=value, fill=variable))+
      geom_bar(stat="identity", width=1)+ # width=1 allows no space between bars
      #scale_fill_grey()+ # Make the bar color in grey scale
      facet_grid( . ~ origin.pop, scales="free_x", space="free_x")+ #scales="free" allows each facet includes the data that exist; space="free" allows facet size being proportionally adjusted
      ylab("Probability")+
      labs(title=paste0("K = ",ans_k,"  ,  ","train locus level = ",ans_t))+
      coord_cartesian(ylim=c(0, 1.005))+ #add 0.005 on y to give tiny space between panel and facet strip
      guides(fill=guide_legend(title=NULL, reverse=T))+ #Hiding title of legend
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),#hiding grid of the panel
            strip.background = element_rect(colour="black", fill="white", linetype="solid"),#change facet title background color
            plot.title = element_text(size=16, vjust=0.8),
            legend.text = element_text(size=14),
            strip.text.x = element_text(size=16),
            axis.title.y = element_text(size=16), axis.text.y = element_text(size=14, colour="black"),
            axis.title.x = element_blank(), axis.text.x = element_text(angle=90, size=7) )
    return(stackplot)

  }else if(style==3){ #Individuals are separated by each fold
    ndf <- melt(df_mas, id.vars=c("Ind.ID","origin.pop","pred.pop","fold_n"))#Reshape the data, making probabilities in one single column (var name="value")
    stackplot <- ggplot(ndf, aes(x=Ind.ID, y=value, fill=variable))+
      geom_bar(stat="identity", width=1)+ # width=1 allows no space between bars
      #scale_fill_grey()+ # Make the bar color in grey scale
      facet_grid( fold_n ~ origin.pop, scales="free_x", space="free_x")+ #scales="free" allows each facet includes the data that exist; space="free" allows facet size being proportionally adjusted
      ylab("Probability")+
      labs(title=paste0("K = ",ans_k,"  ,  ","train locus level = ",ans_t))+
      coord_cartesian(ylim=c(0, 1.005))+ #add 0.005 on y to give tiny space between panel and facet strip
      guides(fill=guide_legend(title=NULL, reverse=T))+ #Hiding title of legend
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),#hiding grid of the panel
            strip.background = element_rect(colour="black", fill="white", linetype="solid"),#change facet title background color
            plot.title = element_text(size=16, vjust=0.8),
            legend.text = element_text(size=14),
            strip.text.x = element_text(size=16),
            axis.title.y = element_text(size=16), axis.text.y = element_text(size=14, colour="black"),
            axis.title.x = element_blank(), axis.text.x = element_text(angle=90, size=7) )
    return(stackplot)

  }else if(style==4){ #Individuals are separated by fold and sorted by probability
    df_mas_4 <- data.frame(matrix(ncol=0,nrow=0))
    for(p in pops){
      df_pop <- df_mas[which(df_mas$origin.pop==p),]#Subset samples for each pop based on the origin.pop variable
      df_pop$Ind.ID <- as.factor(df_pop$Ind.ID) #convert ind id charater to factor data
      df_pop$Ind.ID <- droplevels(df_pop$Ind.ID) #Drop unexisting ind id levels
      df_pop$Ind.ID <- factor(df_pop$Ind.ID, levels=df_pop$Ind.ID[order(df_pop$fold_n, -df_pop[,p])], ordered=T)##Reorder samples by fold.k and then probability values of that pop
      df_mas_4 <- rbind(df_mas_4, df_pop)
    }
    ndf <- melt(df_mas_4, id.vars=c("Ind.ID","origin.pop","pred.pop","fold_n"))#Reshape the data, making probabilities in one single column (var name="value")
    stackplot <- ggplot(ndf, aes(x=Ind.ID, y=value, fill=variable))+
      geom_bar(stat="identity", width=1)+ # width=1 allows no space between bars
      #scale_fill_grey()+ # Make the bar color in grey scale
      facet_grid( fold_n ~ origin.pop, scales="free_x", space="free_x")+ #scales="free" allows each facet includes the data that exist; space="free" allows facet size being proportionally adjusted
      ylab("Probability")+
      labs(title=paste0("K = ",ans_k,"  ,  ","train locus level = ",ans_t))+
      coord_cartesian(ylim=c(0, 1.005))+ #add 0.005 on y to give tiny space between panel and facet strip
      guides(fill=guide_legend(title=NULL, reverse=T))+ #Hiding title of legend
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),#hiding grid of the panel
            strip.background = element_rect(colour="black", fill="white", linetype="solid"),#change facet title background color
            plot.title = element_text(size=16, vjust=0.8),
            legend.text = element_text(size=14),
            strip.text.x = element_text(size=16),
            axis.title.y = element_text(size=16), axis.text.y = element_text(size=14, colour="black"),
            axis.title.x = element_blank(), axis.text.x = element_text(angle=90, size=7) )
    return(stackplot)

  }else {
    ndf <- melt(df_mas, id.vars=c("Ind.ID","origin.pop","pred.pop","fold_n"))#Reshape the data, making probabilities in one single column (var name="value")
    stackplot <- ggplot(ndf, aes(x=Ind.ID, y=value, fill=variable))+
      geom_bar(stat="identity", width=1)+ # width=1 allows no space between bars
      #scale_fill_grey()+ # Make the bar color in grey scale
      facet_grid( . ~ origin.pop, scales="free_x", space="free_x")+ #scales="free" allows each facet includes the data that exist; space="free" allows facet size being proportionally adjusted
      ylab("Probability")+
      labs(title=paste0("K = ",ans_k,"  ,  ","train locus level = ",ans_t))+
      coord_cartesian(ylim=c(0, 1.005))+ #add 0.005 on y to give tiny space between panel and facet strip
      guides(fill=guide_legend(title=NULL, reverse=T))+ #Hiding title of legend
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),#hiding grid of the panel
            strip.background = element_rect(colour="black", fill="white", linetype="solid"),#change facet title background color
            plot.title = element_text(size=16, vjust=0.8),
            legend.text = element_text(size=14),
            strip.text.x = element_text(size=16),
            axis.title.y = element_text(size=16), axis.text.y = element_text(size=14, colour="black"),
            axis.title.x = element_blank(), axis.text.x = element_text(angle=90, size=7) )
    return(stackplot)
  }

} #End
