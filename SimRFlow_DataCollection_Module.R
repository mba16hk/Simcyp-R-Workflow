### SimRFlow data collection module incorporated. The idea is that 
#     1) Searcgh ChEMBL API for all needed data
#     2) Missing LogP will be populated by EPI Suite
#     3) Missing PSA and HBD will be poulated by Pubchem

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#import the scripts 
source('input_query.R')
source('ChEMBL_API.R')
source('PubChem_API.R')
source('EPI_Norman_Search.R')

SimRFlow_DataCollection <- function (data_table, PubChem = T, Norman = T){
  
  ## Search ChEMBL for all parameters
  message('Collecting Data from ChEMBL')
  chembl_data <- ChEMBL_API_Search(data_table)
  chembl_data <- as.data.frame(lapply(chembl_data, unlist))
  chembl_data$Source <- 'ChEMBL API'
  
  if (PubChem==F & Norman == F){
    SimRFlow_Data <- merge(data_table,chembl_data,by='INCHIKEY', all = T)
    return(SimRFlow_Data)
  }
  
  #find missing compounds
  missing_comps_inchi <- data_table$INCHIKEY[which(data_table$INCHIKEY %!in% chembl_data$INCHIKEY)]
  
  
  if (PubChem == T){
    #search for missing compounds in PubChem
    
    
    if(length(missing_comps_inchi)>0){
      
      message('Collecting Data from PubChem')
      missing_compounds <- data_table[which(data_table$INCHIKEY %!in% chembl_data$INCHIKEY),]
      pubchem_data <- PubChem_API_Search(missing_compounds)
      
      #### if there are duplicates (everything identical but SMILES), rmv duplicate
      
      # identify duplicated rows
      duplicated_rows1 <- duplicated(pubchem_data$InChIKey)
      duplicated_rows2 <-duplicated(pubchem_data$InChIKey, fromLast = TRUE)
      duplicated_rows1 <- ifelse(duplicated_rows2==T,T,duplicated_rows1)
      
      duplicates <- pubchem_data[duplicated_rows1,]
      
      #remove duplcated rows from pubchem data
      pubchem_data <- pubchem_data[!duplicated_rows1,]
      
      #check if everything is identical but SMILES, and remove duplicates
      rmv_duplicates_idx <- duplicated(rm_df_cols(duplicates,c('CanonicalSMILES')))
      duplicates <- duplicates[!rmv_duplicates_idx,]
      
      ## Since the duplicates have been removed, add them back to pubchem data
      pubchem_data <- rbind(pubchem_data,duplicates)
      pubchem_data$Source <- 'PubChem API'
      
    }
    
    message('Merging PubChem and ChEMBL Data')
    ## merge pubchem and ChEMBL data together
    Pubchem_Chembl <- merge(chembl_data,pubchem_data,
                            by.x = c('SMILES','INCHIKEY','MolFormula','MW','CXLogP','PSA','HBD','HBA','Source'),
                            by.y = c('CanonicalSMILES','InChIKey','MolecularFormula','MolecularWeight',
                                     'XLogP','TPSA','HBondDonorCount','HBondAcceptorCount','Source'),
                            all = T)
    
    SimRFlow_Data <- Pubchem_Chembl
    
    #Organise dataframe columns
    SimRFlow_Data <- SimRFlow_Data %>% relocate(c('CXLogD','Compound_type','Acid_pka','Base_pka','quat_nitrogens'),
                                                .before = Source)
    SimRFlow_Data <- SimRFlow_Data %>% relocate(c('CHEMBL_id'),.before = MolFormula)
    SimRFlow_Data <- SimRFlow_Data %>% relocate(c('CXLogP'),.after = MW)
    
    if (Norman == F){
      
      SimRFlow_Data <- merge(data_table,SimRFlow_Data,by='INCHIKEY', all = T)
      
      return(SimRFlow_Data)
      
    } else{
      
      message('Collecting Data from EPI Suite')
      ### If there are missing LogP values, get them from EPI Suite Norman Database
      missing_logP_idx <- which(is.na(SimRFlow_Data$CXLogP))
      missing_logP_inchi <- SimRFlow_Data$INCHIKEY[missing_logP_idx]
      missing_logP_table <- data_table[which(data_table$INCHIKEY %in% missing_logP_inchi),]
      EPI_data <- keep_df_cols(NormanSearch(missing_logP_table),c("INCHIKEY","LogP"))
      EPI_data$Source <- 'EPI Suite'
      
      ### merge the missing logP with the ChemBL_PubChem dataframe
      SimRFlow_Data <- merge(SimRFlow_Data,EPI_data,by.x=c('INCHIKEY','CXLogP', 'Source'), by.y = c("INCHIKEY","LogP",'Source'), all = T)
      
      #Organise dataframe columns
      SimRFlow_Data <- SimRFlow_Data %>% relocate(c('CXLogP'),.after = MW)
      SimRFlow_Data <- SimRFlow_Data %>% relocate(c('Source'),.after = ALogP)
      SimRFlow_Data <- merge(data_table,SimRFlow_Data,by='INCHIKEY', all = T)
    }
   
  } else{
    
    if(length(missing_comps_inchi)>0){
      if (Norman == T){
        message('Collecting Data from EPI Suite')
        missing_compounds <- data_table[which(data_table$INCHIKEY %!in% chembl_data$INCHIKEY),]
        EPI_data <- keep_df_cols(NormanSearch(missing_compounds),c("INCHIKEY","LogP"))
        EPI_data$Source <- 'EPI Suite'
        ### merge the missing logP with the ChemBL_PubChem dataframe
        SimRFlow_Data <- merge(chembl_data,EPI_data,by.x=c('INCHIKEY','CXLogP', 'Source'), by.y = c("INCHIKEY","LogP", 'Source'), all = T)
        
        #Organise dataframe columns
        SimRFlow_Data <- SimRFlow_Data %>% relocate(c('CXLogP'),.after = MW)
        SimRFlow_Data <- SimRFlow_Data %>% relocate(c('Source'),.after = ALogP)
        SimRFlow_Data <- merge(data_table,SimRFlow_Data,by='INCHIKEY', all = T)
      }
    }
    
  }
  
  
  return(SimRFlow_Data)
  
}



