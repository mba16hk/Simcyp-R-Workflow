
#search the EPI Suite database for the missing compounds
EPISearch <- function (info, NOT_FOUND, ChEMBL_search){
  
  #set the working directory to the EPI suite database
  epi_dir<-paste0(dirname(dirname(rstudioapi::getSourceEditorContext()$path)),
                  '/data_files/EPISUITE_INCHIKEY.csv')
  
  #read the EPI Suite database
  EPI<-read.csv(epi_dir,sep=',', header=T)
  
  #search EPI suite using CAS registration numbers
  found_in_epi<-filter(EPI, 
                       str_extract(EPI$PUBCHEM_inchikey, "[^-]+") %in%  str_extract(NOT_FOUND$InChiKey, "[^-]+"))
  found_in_epi<-found_in_epi[!duplicated(found_in_epi$PUBCHEM_CID),]
  NA_vals <- which(is.na(found_in_epi$InchiKey))
  if (length(NA_vals)>0){
    found_in_epi <- found_in_epi[-NA_vals,]
  }
  found_in_epi <- found_in_epi[!duplicated(found_in_epi$PUBCHEM_inchikey),]
  
  # data_from_epi<-data.frame()
  # for (cmpd_name in NOT_FOUND$Compound){
  #   if (is.na(cmpd_name)){
  #     next
  #   }
  #   found_compound <- str_detect(toupper(found_in_epi$NAMES),toupper(cmpd_name))
  #   data_from_epi<- rbind(data_from_epi, unlist(found_in_epi[found_compound,]))
  # }
  
  #merge based on inchikey
  NOT_FOUND$shortend_inchi <- str_extract(NOT_FOUND$InChiKey, "[^-]+")
  found_in_epi$shortend_inchi <- str_extract(found_in_epi$PUBCHEM_inchikey, "[^-]+")
  
  epi_dat <- merge(found_in_epi,NOT_FOUND, by = 'shortend_inchi',all.x=T)
  epi_dat$Compound <-toupper(epi_dat$Compound)  
  keep<-c('Compound','InChiKey', 'Code',
          'SMILES', 'EPI_MOLWT','EPI_KOW', 'EPI_PKA')
  found_in_epi<- keep_df_cols(epi_dat, keep)
  found_in_epi$data_source <- 'EPI SUITE (EPA) - Experimental Data'
  
  if (length(keep_df_cols(info,c('DOSE','DOSEUNITS','VSSMETHOD'))) > 0){
    found_in_epi <- merge(found_in_epi, 
                          keep_df_cols(info,c('InChiKey','DOSE','DOSEUNITS','VSSMETHOD')), 
                          by = 'InChiKey', all.x = T)
  }
  
  rm(EPI)
  
  #search for logKow and MLWT values from susdat db (which uses EPI Suite values)
  # sus_dir<-paste0(dirname(dirname(rstudioapi::getSourceEditorContext()$path)),
  #                 '/data_files/Additional_dbs/susdat_2020-06-18-123341.csv')
  # susdat<-read.csv(sus_dir, header=T, sep=',')
  # 
  # found_in_susdat<-filter(susdat,susdat$SMILES %in% NOT_FOUND$SMILES)
  # keep<-c('SMILES', 'CAS_RN_Dashboard','StdInChIKey','Monoiso_Mass',
  #         'logKow_EPISuite', 'Exp_logKow_EPISuite')
  # found_in_susdat<- keep_df_cols(found_in_susdat, keep)
  
  #take the experimental logKow if present. Otherwise take the estimated logKow
  # found_in_susdat$logKow_EPISuite<-ifelse(!is.na(found_in_susdat$Exp_logKow_EPISuite),
  #                                         found_in_susdat$Exp_logKow_EPISuite,
  #                                         found_in_susdat$logKow_EPISuite)
  # found_in_susdat<- rm_df_cols(found_in_susdat, c('Exp_logKow_EPISuite'))
  # rm(susdat,i, keep)
  
  #find compounds unique to susdat and identify their names
  # only_in_susdat<-filter(found_in_susdat, found_in_susdat$SMILES %!in% found_in_epi$SMILES)
  # cmpd_names<-NOT_FOUND[which(NOT_FOUND$SMILES %in% only_in_susdat$SMILES),]
  # cmpd_names<-keep_df_cols(cmpd_names,c('Compound','SMILES'))
  # only_in_susdat<-merge(only_in_susdat,cmpd_names, by='SMILES', all=T)
  
  #merge EPI and susdat data
  # found_in_epi<-merge(found_in_epi, only_in_susdat,
  #                     by.y=c('Compound','SMILES','CAS_RN_Dashboard','StdInChIKey',
  #                            'Monoiso_Mass','logKow_EPISuite'),
  #                     by.x=c('Compound','SMILES', 'CAS', 'InChiKey','MOLWT','KOW'),
  #                     all=T)
  # rm(found_in_susdat,only_in_susdat, cmpd_names)
  
  #find available columns and store in a vector
  potential_headers <- c('DOSE','DOSEUNITS','VSSMETHOD')
  additional_names<-keep_df_cols(info,potential_headers)
  if (length(additional_names)==0){
    additional_names <- c()
  } else{
    indicies <- which(potential_headers %in% colnames(info))
    additional_names <- potential_headers[indicies]
  }
  
  #merge the EPI Suite data with the ChEMBL Search
  ChEMBL_EPI_search<-merge(ChEMBL_search, found_in_epi,
                           by.x = c('COMPOUND.NAME', 'SMILES', 'InChIKey', 'Code',
                                    'MW', 'logPow','Acidic..pKa','data_source', additional_names), 
                           by.y = c('Compound','SMILES', 'InChiKey','Code','EPI_MOLWT',
                                    'EPI_KOW','EPI_PKA','data_source',additional_names),
                           all= T)
  #Organise dataframe columns
  ChEMBL_EPI_search <- ChEMBL_EPI_search %>% relocate(MWfreebase, .after = MW)
  ChEMBL_EPI_search <- ChEMBL_EPI_search %>% relocate(Basic.pKa, .after = Acidic..pKa)
  ChEMBL_EPI_search <- ChEMBL_EPI_search %>% relocate(Molregno, .after = COMPOUND.NAME)
  
  ChEMBL_EPI_search <- ChEMBL_EPI_search %>% relocate(Code, .before = COMPOUND.NAME)
  ChEMBL_EPI_search <- ChEMBL_EPI_search %>% relocate(SMILES, .after = COMPOUND.NAME)
  ChEMBL_EPI_search <- ChEMBL_EPI_search %>% relocate(data_source, .after = test)
  
  #repopulate SMILES and compound code
  missing_codes_indicies <- which(is.na(ChEMBL_EPI_search$Code))
  missing_codes_inchi <- ChEMBL_EPI_search$InChIKey[missing_codes_indicies]
  found_codes_indicies <- which(info$InChiKey %in% missing_codes_inchi)
  
  if (length(missing_codes_indicies)>0){
    ChEMBL_EPI_search$Code[missing_codes_indicies]<-info$Code[found_codes_indicies]
    ChEMBL_EPI_search$SMILES[missing_codes_indicies]<-info$SMILES[found_codes_indicies]
  }
  
  return(ChEMBL_EPI_search)
  
}

NotFoundInEPI <- function (NOT_FOUND, Chemble_epi_file){
  
  #Update NOT_FOUND df
  NOT_FOUND<-filter(NOT_FOUND,NOT_FOUND$SMILES %!in% Chemble_epi_file$SMILES)
  
  NOT_FOUND<- keep_df_cols(NOT_FOUND,c('Code','SMILES'))
  return(NOT_FOUND)
  
}

getCASnumbers<- function(data){
  
  #set the working directory to the EPI suite database
  epi_dir<-paste0(dirname(dirname(rstudioapi::getSourceEditorContext()$path)),
                  '/data_files/EPISUITE_INCHIKEY.csv')
  
  #read the EPI Suite database
  EPI<-read.csv(epi_dir,sep=',', header=T)
  
  found_in_epi<-filter(EPI, 
                       str_extract(EPI$PUBCHEM_inchikey, "[^-]+") %in%  str_extract(data$InChiKey, "[^-]+"))
  found_in_epi<-found_in_epi[!duplicated(found_in_epi$PUBCHEM_CID),]
  NA_vals <- which(is.na(found_in_epi$InchiKey))
  if (length(NA_vals)>0){
    found_in_epi <- found_in_epi[-NA_vals,]
  }
  found_in_epi <- found_in_epi[!duplicated(found_in_epi$PUBCHEM_inchikey),]
  
  keep_cols <- c('PUBCHEM_inchikey','EPI_CAS')
  CASnumbers<- keep_df_cols(found_in_epi,keep_cols)
  
  #merge based on inchikey
  CASnumbers$shortend_inchi <- str_extract(CASnumbers$PUBCHEM_inchikey, "[^-]+")
  data$shortend_inchi <- str_extract(data$InChiKey, "[^-]+")
  
  CAS_data <- merge(CASnumbers,data, by = 'shortend_inchi',all.x=T)
  CAS_data <- CAS_data[!duplicated(CAS_data$shortend_inchi),]
  
  #only keep inchikey and CAS
  keep_cols <- c('InChiKey','EPI_CAS')
  CASNO <- keep_df_cols(CAS_data,keep_cols)
  CASNO <- CASNO %>% rename(CAS = EPI_CAS)
  
  #remove leading zeros in CAS numbers to better query httk database
  CASNO$CAS <- sub("^0+", "", CASNO$CAS)    
  
  return(CASNO)
}
