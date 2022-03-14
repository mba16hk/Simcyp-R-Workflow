
ExpDataSearch <- function (httk_data, experimental_data_directory,
                           CL_threshold = 0, fu_threshold = 0, BP_threshold = 0,
                           mean_flag = 0) {
  
    
  #extract file extensions
  extension <- file_ext(experimental_data_directory)

  
  #load the data based on the extracted extension
  if (extension == 'xlsx'){
    import<-loadWorkbook(experimental_data_directory, create=FALSE)
    exp_data <- readWorksheet(import,  header = TRUE, sheet = 1)
  } else {
    return(cat("File type not compatible."))
  }
  
  #create a column for the case study each compound belongs to
  CS <- sub("\\-.*", "", httk_data$Code)
  httk_data <- cbind(httk_data, CS)
  
  CS <- sub("\\-.*", "", exp_data$CODE)
  exp_data <- cbind(exp_data, CS)
  
  #inchikeys and case-study IDs
  inchi_CSID <- keep_df_cols(httk_data,c('InChIKey','Code'))
  exp_data <- merge (exp_data, inchi_CSID,
                     by.y = 'Code' , by.x = 'CODE' , all.x = T)
  
  # Extract CLint, fu, BP data with their inchi Keys
  relevant_exp_data <- keep_df_cols(exp_data, c('InChIKey','CLINT',
                                                 'BP','FU'))
  
  #ensure only numbers are kept
  relevant_exp_data$CLINT <- suppressWarnings(as.numeric(relevant_exp_data$CLINT))
  relevant_exp_data$BP <- suppressWarnings(as.numeric(relevant_exp_data$BP))
  relevant_exp_data$FU <- suppressWarnings(as.numeric(relevant_exp_data$FU))
  
  #CL_threshold=BP_threshold=fu_threshold=0
  
  #ensure thresholds are met
  indicies_below_threshold <- which(relevant_exp_data$CLINT< CL_threshold)
  relevant_exp_data$CLINT[indicies_below_threshold] <- NA
  
  indicies_below_threshold <- which(relevant_exp_data$BP< BP_threshold)
  relevant_exp_data$BP[indicies_below_threshold] <- NA
  
  indicies_below_threshold <- which(relevant_exp_data$FU< fu_threshold)
  relevant_exp_data$FU[indicies_below_threshold] <- NA
  
  #merge experimental data with httk data
  httk_exp_data<-merge(httk_data, relevant_exp_data, by ="InChIKey", all.x = T)
  
  #Add column for BP source
  indicies<-which(!is.na(httk_exp_data$BP))
  empty_vector<-rep(NA,nrow(httk_exp_data))
  empty_vector[indicies]<-'Experimental Data'
  httk_exp_data$BP_source<-empty_vector
  
  #first merge two columns together for the source and value
  httk_exp_data$FU_value <-ifelse(is.na(httk_exp_data$FU),
                                  httk_exp_data$fu_value,
                                  httk_exp_data$FU)
  httk_exp_data$fu_source <-ifelse(!is.na(httk_exp_data$FU),
                                  'Experimental data',
                                  httk_exp_data$fu_source) 
  
  httk_exp_data$CLint_value_fin <-ifelse(is.na(httk_exp_data$CLINT),
                                  httk_exp_data$CLint_value,
                                  httk_exp_data$CLINT)
  httk_exp_data$CLint_source <-ifelse(!is.na(httk_exp_data$CLINT),
                                   'Experimental data',
                                   httk_exp_data$CLint_source)
  
  #if mean_flag =0, the provided experimental data is preferred over httk
  if (mean_flag == 1) {
    
    #if there are two fu_values, averge them
    httk_exp_data$fu_value <-ifelse(!is.na(httk_exp_data$fu_value.y) & 
                                      !is.na(httk_exp_data$fu_value.x),
                                    (httk_exp_data$fu_value.y + 
                                       httk_exp_data$fu_value.x)/2,
                                    httk_exp_data$fu_value)
    
    #if there are two clint values, average them
    httk_exp_data$CLint_value_fin <-ifelse(!is.na(httk_exp_data$Clint_value) & 
                                      !is.na(httk_exp_data$CLint_value),
                                    (httk_exp_data$Clint_value + 
                                       httk_exp_data$CLint_value)/2,
                                    httk_exp_data$CLint_value_fin)
    
    #Adjust source columns for fu
    httk_exp_data$fu_source <-ifelse(!is.na(httk_exp_data$fu_value.y) & 
                                       !is.na(httk_exp_data$fu_value.x),
                                    paste('Averaged form:',
                                          httk_exp_data$fu_source,
                                          '& Experimental data', sep =' '),
                                    httk_exp_data$fu_source)
    
    #Adjust source columns for CLint
    httk_exp_data$CLint_source <-ifelse(!is.na(httk_exp_data$Clint_value) & 
                                          !is.na(httk_exp_data$CLint_value),
                                     paste('Averaged form:',
                                           httk_exp_data$CLint_source,
                                           '& Experimental data', sep =' '),
                                     httk_exp_data$CLint_source)
    
  } 
  
  #delete unwanted columns
  rmv<-c("CLint_value","CLINT","fu_value","FU","SMILES")
  httk_exp_data<-rm_df_cols(httk_exp_data,rmv)
  
  #organise df columns
  httk_exp_data <- httk_exp_data %>% relocate(CS, .after = Code)
  httk_exp_data <- httk_exp_data %>% relocate(FU_value, .before = fu_units)
  httk_exp_data <- httk_exp_data %>% relocate(CLint_value_fin, .before = CLint_units)
  
  #rename clint column
  httk_exp_data <- httk_exp_data %>% rename(CLint_value = CLint_value_fin)
  httk_exp_data <- httk_exp_data %>% rename(fu_value = FU_value)
  
  #Further Organisation
  httk_exp_data <- httk_exp_data %>% relocate(InChIKey, .after = COMPOUND.NAME)
  httk_exp_data$COMPOUND.NAME <- toupper(httk_exp_data$COMPOUND.NAME)
  
  #organise in ascending CS number
  httk_exp_data <- httk_exp_data[order(httk_exp_data$Code),]

  return(httk_exp_data)

}

CaseStudyTable <- function (httk_exp_data, experimental_data_directory){
  
  #ensure file exists and extract extensions
  if (file.exists(experimental_data_directory)) {
    
    #extract file extensions
    extension <- file_ext(experimental_data_directory)
    
  } else {
    
    return(cat("The file does not exist"))
    
  }
  
  #load the data based on the extracted extension
  if (extension == 'xlsx'){
    import<-loadWorkbook(experimental_data_directory, create=FALSE)
    exp_data <- readWorksheet(import,  header = TRUE, sheet = 1)
  } else {
    return(cat("File type not compatible."))
  }
  
  #identify case studies in experimental data
  CS <- sub("\\-.*", "", exp_data$CODE)
  exp_data <- cbind(exp_data, CS)
  
  #inchikeys and case-study IDs
  inchi_CSID <- keep_df_cols(httk_exp_data,c('InChIKey','Code'))
  exp_data <- merge (exp_data, inchi_CSID,
                     by.y = 'Code' , by.x = 'CODE' , all.x = T)
  
  #Identify Casestudies and freq of occurrence of experimental data
  case_studies <- data.frame(table(httk_exp_data$CS))
  case_studies$CLint_Data <- 0
  case_studies$fu_Data <- 0
  case_studies$BP_Data <- 0
  
  for (i in 1:nrow(case_studies)){
    
    CS_indicies<-which(httk_exp_data$CS == case_studies[i,1])
    
    clint_data <- which(!is.na(httk_exp_data$CLint_value[CS_indicies]))
    case_studies[i,3] <- length(clint_data)
    
    fu_data <- which(!is.na(httk_exp_data$fu_value[CS_indicies]))
    case_studies[i,4] <- length(fu_data)
    
    bp_data <- which(!is.na(httk_exp_data$BP_value[CS_indicies]))
    case_studies[i,5] <- length(bp_data)
    
  }
  
  return(case_studies)
  
}
