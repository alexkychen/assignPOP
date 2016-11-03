#' Compile genetic and other non-genetic data
#'
#' This function allows you to combine genetic and other non-genetic classifier data, such as morphometrics, of the observations for assignment tests.
#' @param x A returned object (list) from the function read.genpop() or reduce.allele().
#' @param add.x A filename of additional non-genetic data that has sample ID in the first column of the data frame. The sample ID of an individual in your GENEPOP and this additional data must be identifcal.
#' @param method A method to match sample ID between genetic and non-genetic data. The "common" method only concatenate the data that has sample ID in both files. If an individual only exists in one of the files, this individual will be discarded.
#' @return This function returns a new object (list) that comprises 5 items. [[1]] data matrix including genetic and non-genetic data, [[2]] a sample ID vector, [[3]] a locus name vector, [[4]] a vector of non-genetic variable names, and [[5]] the number of non-genetic variables.
#' @examples # infile_com <- compile.data(x, "YourAddFile.txt")
#' @import stringr
#' @importFrom reshape2 melt
#' @export
#'
compile.data <- function(x, add.x, method="common"){
  #Read genetic and non-genetic data
  genoMatrix <- x[[1]]
  if(grepl(pattern=".csv", add.x)){
    cat("  Import a .CSV file.")
    add.df <- read.csv(add.x, header=T)

  }else {
    cat("  Import a table which elements separated by space.")
    add.df <- read.table(add.x, header=T)
  }
  #Analyze non-genetic data
  varNames <- names(add.df[2:ncol(add.df)])#get variable names (exclude ind ID column)
  noVars <- length(varNames)#count number of variables
  cat(paste0("\n  ",noVars," additional variables detected."))
  cat("\n  Checking variable data type...\n")
  for(i in 1:noVars){
    var_name <- names(add.df)[i+1]
    var_type <- class(add.df[,i+1])
    cat(paste0("   ",var_name,"(",var_type,")"))
  }
  ans_0 <- readline("  Are they correct? (enter Y/N): ")
  if(grepl(pattern="N",toupper(ans_0))){
    cat("  please enter variable names for changing data type (separate names by a whitespace if multiple)\n")
    ans_1 <- readline("  enter here: ")
    ans_1 <- str_trim(ans_1, side="both")
    ans_1 <- unlist(strsplit(ans_1,split=" "))#check out variable name to be processed
    noChangeVar <- length(ans_1)
    #Process variables and convert factor data to dummy variable (binary data)
    for(name in ans_1){
      ans_2 <- readline(paste0("  Which data type should '",name,"' be? (enter numeric or factor): "))
      if(grepl(pattern="N",toupper(ans_2))){
        add.df[,name] <- as.numeric(as.character(add.df[,name]))
      }else if(grepl(pattern="F",toupper(ans_2))){
        add.df[,name] <- as.factor(add.df[,name])
      }
    }
    #Look through non-genetic variables and convert to dummy if is factor
    for(name in varNames){
      if(is.factor(add.df[,name])){
        #Convert factor variable to numeric binary variable (dummy variable)
        dummyData <- as.data.frame(model.matrix( ~ add.df[,name]-1, data=add.df))#get dummy variable data frame
        names(dummyData) <- substring(names(dummyData), 15, 1000L)#extract meaningful wording, or remove some funny wording
        names(dummyData) <- sub("\\b", paste0(name,"."), names(dummyData))#appedn original variabel name at the beginning
        add.df[,name] <- NULL #remove original factor data column
        add.df <- cbind(add.df, dummyData) #append new dummy variable column
      }
    }
  }else if(grepl(pattern="Y",toupper(ans_0))){
    #Look through non-genetic variables and convert to dummy if is factor
    for(name in varNames){
      if(is.factor(add.df[,name])){
        #Convert factor variable to numeric binary variable (dummy variable)
        dummyData <- as.data.frame(model.matrix( ~ add.df[,name]-1, data=add.df))#get dummy variable data frame
        names(dummyData) <- substring(names(dummyData), 15, 1000L)#extract meaningful wording, or remove some funny wording
        names(dummyData) <- sub("\\b", paste0(name,"."), names(dummyData))#appedn original variabel name at the beginning
        add.df[,name] <- NULL #remove original factor data column
        add.df <- cbind(add.df, dummyData) #append new dummy variable column
      }
    }
  }
  #Concatenate genetic and non-geneitc data
  if(method=="common"){
    #Identify individual IDs
    geneData_indID <- x[[2]] #individual ID from genetic data
    addData_indID <- as.character(add.df[,1]) #individual ID from non-genetic data
    common_indID <- intersect(geneData_indID, addData_indID) #get the common ind ID of two datasets
    #Subset data by ind ID
    genoMatrix_wID <- cbind(geneData_indID, genoMatrix)#Concatenate ind ID back to genoMatrix
    genoMatrix_com <- genoMatrix_wID[(genoMatrix_wID$geneData_indID %in% common_indID),]#subset genoMatrix data based on common ind ID
    add.df_com <- add.df[(add.df[,1] %in% common_indID),]#subset non-genetic data based on common ind ID
    #Reorder non-genetic data rows by genoMatrix_com ID column
    add.df_com <- add.df_com[match(genoMatrix_com[,1], add.df_com[,1]), ]
    #Insert non-genetic data to genoMatrix, if two data sets have exact same ind ID
    if(identical(as.character(add.df_com[,1]),as.character(genoMatrix_com[,1]))){
      #Concatenate three items: new genoMatrix data, new non-geneitc data, and pop name column
      comMatrix <- cbind(genoMatrix_com[,2:(ncol(genoMatrix_com)-1)], add.df_com[,2:ncol(add.df_com)], genoMatrix_com$popNames_vector)
      colnames(comMatrix)[ncol(comMatrix)] <- "popNames_vector"#rename the last column
      rownames(comMatrix) <- NULL
    }else {
      stop("  Error: Individual ID are not identical between two data sets")
    }
    #Count number of non-genetic columns (new variables)
    noVarCols <- ncol(add.df_com)-1 #total number of columns minus first ID column
    #Print some message
    cat("\n  New data set created!!")
    cat(paste0("\n  It has ",nrow(comMatrix)," observations by ",ncol(comMatrix)," variables"))
    cat(paste0("\n  including ",length(x[[3]])," loci(",ncol(genoMatrix)-1," alleles) plus ",noVars," additional variables(",ncol(add.df)-1," columns)"))

    return(list(comMatrix, common_indID, x[[3]], varNames, noVarCols ))

  }else if(method=="all"){
    ###(reserve for future development)
  }
}
