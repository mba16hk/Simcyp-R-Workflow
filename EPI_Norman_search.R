### Epi Suite Search

NormanSearch <- function (info, logP_selection_flag =1 ){
  
  #set the working directory to the norman suspect list database (which uses EPI Suite values)
  sus_dir<-paste0(dirname(dirname(rstudioapi::getSourceEditorContext()$path)),
                  '/data_files/Norman_susdat.csv')
  
  #read the suspect database
  susdat<-read.csv(sus_dir,sep=',', header=T)
  colnames(susdat)[1] <- 'Norman_ID'
  
  #search the suspect database using standard inchikeys and MS ready inchikeys
  found_in_susdat<-filter(susdat, 
                          susdat$StdInChIKey %in%  info$INCHIKEY)
  
  #remove duplicates based on the Norman ID
  found_in_susdat<-found_in_susdat[!duplicated(found_in_susdat$Norman_ID),]
  
  #remove compounds without a molecular weight
  found_in_susdat<-found_in_susdat[!is.na(found_in_susdat$Monoiso_Mass),]
  
  #if one of the inchikeys provided is NA, remove the values
  NA_vals <- which(is.na(c(found_in_susdat$MS_Ready_StdInChIKey,
                           found_in_susdat$StdInChIKey)))
  if (length(NA_vals)>0){
    found_in_susdat <- found_in_susdat[-NA_vals,]
  }
  
  keep<-c('StdInChIKey','SMILES','Monoiso_Mass', 'Molecular_Formula', 'StdInChI',
          'logKow_EPISuite', 'Exp_logKow_EPISuite','PubChem_CID','CAS_RN_Dashboard')
  found_in_susdat<- keep_df_cols(found_in_susdat, keep)
  
  #merge codes
  colnames(info) <- toupper(colnames(info))
  code_idx <- which(info$INCHIKEY %in% found_in_susdat$StdInChIKey)
  data_codes <- info[code_idx,c('INCHIKEY','CODE')]
  found_in_susdat <- merge(found_in_susdat,data_codes,
                           by.x='StdInChIKey',by.y='INCHIKEY',all=T)
  
  #convert to numeric values
  found_in_susdat$Monoiso_Mass <- as.numeric(found_in_susdat$Monoiso_Mass)
  found_in_susdat$logKow_EPISuite <- as.numeric(found_in_susdat$logKow_EPISuite)
  found_in_susdat$Exp_logKow_EPISuite <- as.numeric(found_in_susdat$Exp_logKow_EPISuite)
  
  
  #allow the user to prioritise which logKow data to select from
  if (logP_selection_flag == 1){
    
    #take the experimental logKow if present. Otherwise take the estimated logKow
    found_in_susdat$logKow_EPISuite<-ifelse(!is.na(found_in_susdat$Exp_logKow_EPISuite),
                                            found_in_susdat$Exp_logKow_EPISuite,
                                            found_in_susdat$logKow_EPISuite)
    #organise the data source references
    found_in_susdat$data_source <- ifelse(!is.na(found_in_susdat$Exp_logKow_EPISuite),
                                          'EPI SUITE (EPA) - Experimental Data', 'EPI SUITE (EPA)')
    found_in_susdat$data_source <- ifelse(!is.na(found_in_susdat$logKow_EPISuite) & 
                                            found_in_susdat$data_source!='EPI SUITE (EPA) - Experimental Data',
                                          'EPI SUITE (EPA) - Estimated Data',
                                          found_in_susdat$data_source)
    
    #remove unwanted column
    found_in_susdat<- rm_df_cols(found_in_susdat, c('Exp_logKow_EPISuite'))
    
    #rename Kow column to LogP, detect which column contains Kow
    kow_index<-str_detect(colnames(found_in_susdat),'Kow')
    colnames(found_in_susdat)[kow_index] <- 'LogP'
    
  } else if (logP_selection_flag == 2) {
    
    #take the predicted logKow if present. Otherwise take the experimental logKow
    found_in_susdat$Exp_logKow_EPISuite<-ifelse(!is.na(found_in_susdat$logKow_EPISuite),
                                                found_in_susdat$logKow_EPISuite,
                                                found_in_susdat$Exp_logKow_EPISuite)
    #organise the data source references
    found_in_susdat$data_source <- ifelse(!is.na(found_in_susdat$logKow_EPISuite),
                                          'EPI SUITE (EPA) - Estimated Data', 'EPI SUITE (EPA)')
    found_in_susdat$data_source <- ifelse(!is.na(found_in_susdat$Exp_logKow_EPISuite) &
                                            found_in_susdat$data_source!='EPI SUITE (EPA) - Estimated Data',
                                          'EPI SUITE (EPA) - Experimental Data',
                                          found_in_susdat$data_source)
    
    #remove unwanted column
    found_in_susdat<- rm_df_cols(found_in_susdat, c('logKow_EPISuite'))
    
    #rename Kow column to LogP, detect which column contains Kow
    kow_index<-str_detect(colnames(found_in_susdat),'Kow')
    colnames(found_in_susdat)[kow_index] <- 'LogP'
    
  } else{
    
    return(cat('Set logP_selection_flag to 1 for prioritizing experimental data, 
               or Set logP_selection_flag to 2 for prioritizing predicted data'))
    
  }
  
  
  if (length(keep_df_cols(info,c('DOSE','DOSEUNITS','VSSMETHOD'))) > 0){
    found_in_susdat <- merge(found_in_susdat, 
                             keep_df_cols(info,c('INCHIKEY','DOSE','DOSEUNITS','VSSMETHOD')), 
                             by.x = 'StdInChIKey', by.y = 'INCHIKEY', all.x = T)
  }
  
  rm(susdat, keep)
  

  colnames(found_in_susdat)[1] <- 'INCHIKEY'
  colnames(found_in_susdat)[2] <- 'CAS'
  colnames(found_in_susdat)[7] <- 'MW'
  return(found_in_susdat)

  
}

NotFoundInsusdat <- function (NOT_FOUND, ChEMBL_sus_search){
  
  #Update NOT_FOUND df
  NOT_FOUND<-filter(NOT_FOUND,NOT_FOUND$INCHIKEY %!in% ChEMBL_sus_search$InChIKey)
  
  NOT_FOUND<- keep_df_cols(NOT_FOUND,c('CODE','SMILES'))
  return(NOT_FOUND)
  
}

CAS_and_DTXSID<- function(info){
  
  #set the working directory to the norman suspect list database (which uses EPI Suite values)
  sus_dir<-paste0(dirname(dirname(rstudioapi::getSourceEditorContext()$path)),
                  '/data_files/Norman_susdat.csv')
  
  #read the suspect database
  susdat<-read.csv(sus_dir,sep=',', header=T)
  colnames(susdat)[1] <- 'Norman_ID'
  
  #search the suspect database using standard inchikeys and MS ready inchikeys
  found_in_susdat<-filter(susdat, 
                          susdat$StdInChIKey %in%  info$INCHIKEY)
  
  #remove duplicates based on the Norman ID if the duplicates are present
  if (TRUE %in% duplicated(found_in_susdat$Norman_ID)){
    found_in_susdat<-found_in_susdat[!duplicated(found_in_susdat$Norman_ID),]
  }
  
  #only keep the columns of interest
  keep_cols <- c('CAS_RN','StdInChIKey','DTXSID')
  CAS_DTXSID<- keep_df_cols(found_in_susdat,keep_cols)
  
  #merge the info dataframe with the CAS_DTXSID by inchikey
  CAS_DTXSID <- merge(CAS_DTXSID, info, by.x = "StdInChIKey", by.y = "INCHIKEY", all= T)
  
  #if there are compounds with missing CAS values remove them
  missing_cas <- which(is.na(CAS_DTXSID$CAS_RN))
  if (length(missing_cas) > 0){
    CAS_DTXSID <- CAS_DTXSID[-missing_cas,]
  }
  
  #remove prefix in CAS_RN
  CAS_DTXSID$CAS_RN <- gsub("CAS_RN:\\s","",CAS_DTXSID$CAS_RN)
  
  #rename CAS column
  CAS_DTXSID <- CAS_DTXSID %>% rename(CAS = CAS_RN)
  
  #if there is an OR seprator (|), make multiple entries for the CAS
  multiple_CAS_indicies <- which(str_detect(CAS_DTXSID$CAS, "\\|")==T)
  multiple_CAS <- CAS_DTXSID[multiple_CAS_indicies,]
  
  if (nrow(multiple_CAS>0)){
    new_data <- data.frame()
    for (i in 1:nrow(multiple_CAS)){
      
      Additional_CAS_numbers <- unlist(str_split(multiple_CAS$CAS[i],'\\|'))
      number_of_new_CAS <- length(Additional_CAS_numbers)
      data <- as.data.frame(cbind(rep(multiple_CAS$StdInChIKey[i],number_of_new_CAS),
                                  Additional_CAS_numbers,
                                  rep(multiple_CAS$DTXSID[i],number_of_new_CAS),
                                  rep(multiple_CAS$CODE[i],number_of_new_CAS),
                                  rep(multiple_CAS$COMPOUND[i],number_of_new_CAS),
                                  rep(multiple_CAS$SMILES[i],number_of_new_CAS)))
      colnames(data) <- colnames(CAS_DTXSID)
      #print(CAS_DTXSID)
      new_data <- rbind(new_data,data)
      
    }
    #remove the entries with "|"
    CAS_DTXSID <- CAS_DTXSID[-multiple_CAS_indicies,]
    #add the new CAS data
    CAS_DTXSID <- rbind(CAS_DTXSID,new_data)
  }
  
  return(CAS_DTXSID)
}

