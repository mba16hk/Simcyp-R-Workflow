#load relevant libraries
library(XLConnect)
library(tools)
library(dplyr)

#import all additional functions
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('additional_functions.R')

#function that formats the input compounds
ProcessInputs <- function (file_directory){
  
  #ensure file exists and extract extensions
  if (file.exists(file_directory)) {
    
    #extract file extensions
    extension <- file_ext(file_directory)
    
  } else {
    
    return(cat("The file does not exist"))
    
  }
  
  
  #load the data based on the extracted extension
  if (extension == 'xlsx'){
    
    import<-loadWorkbook(file_directory, create=FALSE)
    cmpd_data <- readWorksheet(import,  header = TRUE, sheet = 1)
    
  } else if (extension == 'csv'){
    
    cmpd_data <- read.csv(file_directory, header = T, sep = ',')
    
  } else {
    
    return(cat("File type not compatible."))
    
  }
  
  data_headers<-toupper(colnames(cmpd_data))
  colnames(cmpd_data)<-data_headers
  
  #Keep columns of interest to us
  data <- data.frame(cmpd_data$CODE,
                     cmpd_data$COMPOUND,
                     cmpd_data$INCHIKEY,
                     cmpd_data$SMILES)
  
  indicies <- which(colnames(cmpd_data) %in% c('DOSE','DOSEUNITS','VSSMETHOD'))
  additional_headers <- colnames(cmpd_data)[indicies]
  additional_data <- as.data.frame(cmpd_data[,indicies])
  colnames(additional_data)<-additional_headers
  
  if (ncol(additional_data)!=0){
    data<- cbind(data, additional_data)
  }
  
  # Rules for input data #
  # - All capital letters
  # - no special characters or any punctuation
  # - no spaces
  
  #Rename the columns
  X <- c("Code","Compound", "InChiKey", "SMILES", c(additional_headers))
  colnames(data) <- X
  
  #move around additional columns after CAS or before SMILES
  data <- data %>% relocate(c(additional_headers), .before = SMILES)
  
  return(data)
}
