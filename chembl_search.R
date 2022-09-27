#load libraries of interest
library(RSQLite)

CHEMBLSearch <- function (data){
  
  #navigate to where chembl SQL is saved
  newwd<-paste0(dirname(dirname(rstudioapi::getSourceEditorContext()$path)),'/data_files/chembl_29_sqlite')

  filenamedb <- paste0(newwd,"/chembl_29.db")
  col_idx <- which(toupper(colnames(data))=='INCHIKEY')
  query <- data[,col_idx]
  
  ChEMBL_search <- InChIKey_search(filenamedb, query)
  
  if (nrow(ChEMBL_search)==0){
    message('No Data has been found in ChEMBL')
    return(ChEMBL_search)
  }
  
  #Names from chembl are sometimes NA, repopulate from data
  nameless_compounds_inchikeys<-ChEMBL_search$InChIKey[which(is.na(ChEMBL_search$COMPOUND.NAME))]
  ordered_inchi_data <- data[match(nameless_compounds_inchikeys, data[,col_idx]), ]
  ChEMBL_search$COMPOUND.NAME[which(is.na(ChEMBL_search$COMPOUND.NAME))]<-toupper(ordered_inchi_data$COMPOUND)

  ## check search results for duplicates based on molregno, SMILES, ChEMBLID
  ChEMBL_search_CHECK1 <- ChEMBL_search[!duplicated(ChEMBL_search$Molregno),]
  ChEMBL_search_CHECK2 <- ChEMBL_search_CHECK1[!duplicated(ChEMBL_search_CHECK1$SMILES),]
  ChEMBL_search_CHECK3 <- ChEMBL_search_CHECK2[!duplicated(ChEMBL_search_CHECK2$ChEMBL.ID),]

  ## remove duplicates from search results and proceed
  ChEMBL_search <- ChEMBL_search_CHECK3
  
  ##merge it with CAS and Code values
  ChEMBL_search<- merge(ChEMBL_search, data,
             by.x = 'InChIKey', by.y = 'INCHIKEY')
  
  #remove unwanted cols
  if ('SMILES'%in%colnames(data)){
    ChEMBL_search <- rm_df_cols(ChEMBL_search,c('COMPOUND','SMILES.y'))
    ChEMBL_search <- ChEMBL_search %>% rename(SMILES = SMILES.x)
  }
  
  #add column for source of data
  ChEMBL_search$data_source <- 'ChEMBL v29'
  
  ChEMBL_search <- ChEMBL_search %>% relocate(CODE, .before = InChIKey)
  ChEMBL_search <- ChEMBL_search %>% relocate(COMPOUND.NAME, .before = InChIKey)
  
  return(ChEMBL_search)
}

NotInChEMBL<- function (data, chembl_file){
  
  ## report compounds not found based on InChiKey search
  NOT_FOUND <- filter(data, !(data$INCHIKEY %in% chembl_file$InChIKey))
  
  return(NOT_FOUND)
  
}
