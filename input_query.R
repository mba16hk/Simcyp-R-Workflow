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
    
    cmpd_data <- read.csv(file_directory, header = T, sep = ',',fileEncoding = 'UTF-8-BOM')
    
  } else {
    
    return(cat("File type not compatible."))
    
  }
  
  data_headers<-toupper(colnames(cmpd_data))
  colnames(cmpd_data)<-data_headers
  
  #Keep columns of interest to us
  data <- cmpd_data[,c('CODE','SMILES','COMPOUND','INCHIKEY')]
  
  if ('INCHIKEY' %!in% colnames(data)){
    return(message('Please provide compound inchikeys and ensure the inchikey header is spelled as inchikey'))
  } else if ('CODE' %!in% colnames(data)){
    return(message('Please provide compound codes and ensure the code header is spelled as code'))
  }
  
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
  X <- c(toupper(colnames(data)))
  colnames(data) <- X
  
  #move around additional columns after CAS or before SMILES
  if ('SMILES'%in%colnames(data)){
    data <- data %>% relocate(c(additional_headers), .before = SMILES)
  }
  
  return(data)
}
