#load the httk library
library(httk)

httkSearch <- function (physchem_data, CAS_DTXSID, info,
                        
                        #allow user to specify which of the 3 databases to use (default is to use all of them)
                        obach = T, wambaugh = T, chem_phys_in_vitro = T,
                        
                        #apply an operation on the collected fu and Clint values
                        fu_operation = 'arithmetic mean', CLint_operation = 'arithmetic mean'){
  
  if (nrow(CAS_DTXSID) == 0) {
    return(message('no data from httk can be found for any of the compounds'))
  }
  
  if (obach != T & wambaugh != T & chem_phys_in_vitro != T){
    return(message('To use this function, you must search through at least 1 of the 3 databases. Set one of the database flags to T.'))
  }
  #------------ Obach2008 -------------------------------------
  
  if (obach == T){ #if flag is set to TRUE search obach dataset
    
    #Start by searching Obach for Clint & fu using CAS numbers
    Obach_data<-filter(Obach2008, Obach2008$CAS %in% CAS_DTXSID$CAS)
    if (nrow(Obach_data)>0){
      Obach_data$Obach2008_source <- 'Obach 2008'
    }
    keep<- c('CAS','fu','CL (mL/min/kg)','Obach2008_source')
    Obach_data <- keep_df_cols(Obach_data,keep)
  }
  

  #------------ wambaugh2019 -----------------------------------
  
  if(wambaugh == T){ #if flag is set to TRUE search wambaugh dataset
    
    #next is wambaugh 2019 database which contains clint, fu, logP, pKa using CAS
    wambaugh2019_data<-filter(wambaugh2019, wambaugh2019$CAS %in% CAS_DTXSID$CAS)
    if (nrow(wambaugh2019_data)>0){
      wambaugh2019_data$wambaugh2019_source <- 'Wambaugh 2019'
    }
    keep<- c('CAS','Human.Clint.Point','Human.Funbound.plasma.Point',
             'DSSTox_Substance_Id','wambaugh2019_source')
    wambaugh2019_data <- keep_df_cols(wambaugh2019_data,keep)
    
    #wambaugh2019_data <- wambaugh2019_data %>% rename(DTXSID = DSSTox_Substance_Id)
  }
  
  
  #------------ chem.physical_and_invitro.data ------------------
  
  if (chem_phys_in_vitro == T){ #if flag is set to TRUE search chem.physical_and_invitro.data
    
    #use chem.physical_and_invitro.data  to get more information on CLint, fu, and BP 
    # query using using DTXSID and CAS
    in_vitro_data<-filter(chem.physical_and_invitro.data,
                          chem.physical_and_invitro.data$CAS %in% CAS_DTXSID$CAS)
    
    #only keep human data
    human_in_vitro <- in_vitro_data[str_detect(in_vitro_data$All.Species,'Human'),]
    
    
    #only keep columns of interest
    keep <- c('CAS','Human.Clint','Human.Clint.Reference', 'logP', 'MW',
              'Human.Funbound.plasma','Human.Funbound.plasma.Reference',
              'Human.Rblood2plasma','Human.Rblood2plasma.Reference','DTXSID',
              "pKa_Accept","pKa_Accept.Reference","pKa_Donor","pKa_Donor.Reference")
    human_in_vitro <- keep_df_cols(human_in_vitro, keep)
    colnames(human_in_vitro) <- c('CAS','DTXSID', "logP","MW","pKa_Accept","pKa_Accept.Reference","pKa_Donor",
                                  "pKa_Donor.Reference",'CLint1','Clint_source_iv','fu1','PPB_source_iv','BP','BP_source')
    
    #ensure that the CAS and DTXSID match 
    matched_indicies <- match.df(human_in_vitro,CAS_DTXSID,c('CAS','DTXSID'))
    human_in_vitro <- human_in_vitro[matched_indicies,]
    
    #find the arithmetic means of CLints and Fu Values
    for (i in 1:nrow(human_in_vitro)){
      
      human_in_vitro$CLint_iv[i] <- avg_string_of_nums(human_in_vitro$CLint1[i],
                                                       separator = ',|;')
      
      human_in_vitro$fu[i] <- avg_string_of_nums(human_in_vitro$fu1[i],
                                                 separator = ',|;')
    }
    
    #remove unwanted columns
    rm_cols<- c('CLint1','fu1')
    human_in_vitro <- rm_df_cols(human_in_vitro, rm_cols)
    
    # #remove rows made of NAs by finding indicies
    # ind <- apply(human_in_vitro[,3:ncol(human_in_vitro)], 1, function(x) all(is.na(x)))
    # indicies <- which(ind == TRUE)
    # human_in_vitro <- human_in_vitro[-indicies,]
    
    #-----pKa values processing --------
    
    #if values are None, convert to NA
    human_in_vitro$pKa_Accept<- ifelse(human_in_vitro$pKa_Accept == 'None',
                                       NA, human_in_vitro$pKa_Accept)
    human_in_vitro$pKa_Donor<- ifelse(human_in_vitro$pKa_Donor == 'None',
                                      NA, human_in_vitro$pKa_Donor)
    
    #remove the source if there is no value associated with it
    human_in_vitro$pKa_Accept.Reference<- ifelse(is.na(human_in_vitro$pKa_Accept),
                                                 NA, human_in_vitro$pKa_Accept.Reference)
    human_in_vitro$pKa_Donor.Reference<- ifelse(is.na(human_in_vitro$pKa_Donor),
                                                NA, human_in_vitro$pKa_Donor.Reference)
    
    
    #Choose correct pKa values using the pKa selection function
    human_in_vitro$pKa_basic <- NA
    human_in_vitro$pKa_acidic <- NA
    for (i in 1:nrow(human_in_vitro)){
      
      #take the largest pKa value for bases
      human_in_vitro$pKa_basic[i] <- pKa_selection(human_in_vitro$pKa_Accept[i],
                                                   separator = ',|;',
                                                   selection = 'biggest')
      
      #take the smallest pKa value for acids
      human_in_vitro$pKa_acidic[i] <- pKa_selection(human_in_vitro$pKa_Donor[i],
                                                    separator = ',|;',
                                                    selection = 'smallest')
    }
    
    #remove unwanted columns
    human_in_vitro <- rm_df_cols(human_in_vitro,c('pKa_Accept','pKa_Donor'))
    
    #rename reference columns
    human_in_vitro <- human_in_vitro %>% rename(pKa_acidic_ref = pKa_Donor.Reference)
    human_in_vitro <- human_in_vitro %>% rename(pKa_basic_ref = pKa_Accept.Reference)
  }
  
  #------------ merge the data according to the user-selected databases to use ########
  
  if (obach == T & wambaugh == T & chem_phys_in_vitro == T){
    
    #merge all data from 3 databases into a single dataframe
    httk_data_1 <- merge(Obach_data,wambaugh2019_data,by='CAS',all = T)
    httk_data <- merge(httk_data_1,human_in_vitro, by = 'CAS', all = T)
    
  } else if (obach == T & wambaugh == T & chem_phys_in_vitro != T) {
    
    #merge 2 databases into a single dataframe
    httk_data <- merge(Obach_data, wambaugh2019_data, by = 'CAS', all = T)
    
  } else if (obach == T & wambaugh != T & chem_phys_in_vitro == T) {
    
    #merge 2 databases into a single dataframe
    httk_data <- merge(Obach_data, human_in_vitro, by = 'CAS', all = T)
    
  } else if (obach != T & wambaugh == T & chem_phys_in_vitro == T) {
    
    #merge 2 databases into a single dataframe
    httk_data <- merge(wambaugh2019_data, human_in_vitro, by = 'CAS', all = T)
    
  } else if (obach != T & wambaugh != T & chem_phys_in_vitro == T) {
    
    #merge 2 databases into a single dataframe
    httk_data <- human_in_vitro
    
  } else if (obach != T & wambaugh == T & chem_phys_in_vitro != T) {
    
    #merge 2 databases into a single dataframe
    httk_data <- wambaugh2019_data
    
  } else if (obach == T & wambaugh != T & chem_phys_in_vitro != T) {
    
    #merge 2 databases into a single dataframe
    httk_data <- Obach_data
    
  }
  
  #################################################################################
  
  #-----------------process fu data------------------------------------
  
  #find the columns which contain fraction unbound_data
  fraction_unbound_vector <- str_detect(colnames(httk_data),'fu|Funbound')
  fu_vector_to_remove <- colnames(httk_data)[fraction_unbound_vector]
  fraction_unbound_df <- httk_data[,fraction_unbound_vector]
  
  #convert to dataframe with only numeric entries
  if (is.null(ncol(fraction_unbound_df))){
    fraction_unbound_df <- suppressWarnings(as.numeric(fraction_unbound_df))
  } else{
    fraction_unbound_df <- suppressWarnings(as.data.frame(apply(fraction_unbound_df, 2,function(x) as.numeric(x))))
  }
  

  if (is.null(ncol(fraction_unbound_df))){
    
    #if there is only 1 source, take the data from it only
    httk_data$fu_res <- fraction_unbound_df
    
  } else{
    
    #if there is more than 1 source, apply the different operations
    
    if (fu_operation == 'arithmetic mean'){
      
      httk_data$fu_res <- rowMeans(fraction_unbound_df, na.rm = T)
      
    } else if (fu_operation == 'geometric mean'){
      
      mean_log <- rowMeans(log(fraction_unbound_df), na.rm = T)
      geom_mean <- exp(mean_log)
      httk_data$fu_res <- geom_mean
      
    } else if (fu_operation == 'median'){
      
      httk_data$fu_res <- apply(fraction_unbound_df, 1,function(x) median(as.numeric(x),na.rm = T) )
      
    } else if (fu_operation == 'None'){
      
      httk_data$fu_res <- apply(fraction_unbound_df,1,function(x) paste0(x[!is.na(x)], collapse = ' & '))
      
    } else if (fu_operation == 'minimum'){
      
      httk_data$fu_res <- suppressWarnings(apply(fraction_unbound_df,1,function(x) min(as.numeric(x),na.rm = T)))
      httk_data$fu_res <- ifelse(httk_data$fu_res=='Inf',NA, httk_data$fu_res)
      
    } else if (fu_operation == 'maximum'){
      
      httk_data$fu_res <- suppressWarnings(apply(fraction_unbound_df,1,function(x) max(as.numeric(x),na.rm = T) ))
      httk_data$fu_res <- ifelse(httk_data$fu_res=='Inf',NA, httk_data$fu_res)
      
    } else{
      
      fu_operation <- 'arithmetic mean'
      httk_data$fu_res <- rowMeans(fraction_unbound_df, na.rm = T)
      message('arithmetic mean applied. check the spelling of your input to the fu_operation argument')
      
    }
    
  }
  
  #reorganise dataframe
  httk_data <- httk_data %>% relocate(fu_res, .after = CAS)
  if (chem_phys_in_vitro == T){ #reorganise the DTXSID
    httk_data <- httk_data %>% relocate(DTXSID, .after = CAS)
    
  } else if (wambaugh == T & chem_phys_in_vitro != T){
    httk_data <- httk_data %>% relocate(DSSTox_Substance_Id, .after = CAS)
  }
  
  
  #organise sources for fraction unbound
  potential_sources <- c('Obach2008_source|wambaugh2019_source|PPB_source_iv')
  available_sources <- str_detect(colnames(httk_data), potential_sources)
  available_sources <- httk_data[,available_sources]
    
  #replace empty space sources with NAs
  available_sources <- replace(available_sources,available_sources=='',NA)
    
  #check if sources are duplicated
  if (is.null(ncol(available_sources))){
    fu_sources <- tolower(available_sources)
  } else{
    available_sources_upper <- as.data.frame(apply(available_sources,2,function(x) tolower(x))) #convert to upper case
    final_sources <- apply(available_sources_upper,1,function(x) unique(x[!is.na(x)]))
    
    source_num<- lengths(final_sources)
    
    fu_sources <- c()
    for (k in 1:length(final_sources)){
      
      if(source_num[k]>1){
        
        data <- paste0(final_sources[[k]],collapse = ' & ')
        
      } else if (source_num[k]==0){
        
        data <- NA
        
      } else{
        
        data <- final_sources[[k]]
      }
      
      fu_sources[k]<- data
      
    }
  }
  
  #if there is only 1 fu source no need to apply an operation on it, therefore only report 1 source without operation  
  httk_data$fu_source <- fu_sources
  
  if (fu_operation != 'None'){
    for (i in 1:nrow(httk_data)){
      
      if (is.na(httk_data$fu_source[i])){
        next
      }
      
      if (str_detect(httk_data$fu_source[i],'&')){
        httk_data$fu_source[i] <- paste('httk:', fu_operation, 'of', httk_data$fu_source[i], sep=' ')
        
      } else {
        httk_data$fu_source[i] <- paste('httk:', httk_data$fu_source[i], sep=' ')
        
      }
    }
  }
  
  #remove unwanted columns
  rm_cols <- c(fu_vector_to_remove,'PPB_source_iv')
  httk_data <- rm_df_cols(httk_data,rm_cols)
  
  #set the units for the fraction unbound
  httk_data$fu_units = "fraction"
  
  #calculate the standard deviation across available sources:
  httk_data$fu_SD<-apply(fraction_unbound_df,1,function(x) sd(x, na.rm=T))
  
  #organise the httk_data columns
  httk_data <- httk_data %>% relocate(fu_units, .after = fu_res)
  httk_data <- httk_data %>% relocate(fu_SD, .after = fu_res)
  httk_data <- httk_data %>% relocate(fu_source, .after = fu_units)
  
  #-----------------process Clint data------------------------------------
  
  if (wambaugh == T){ #pre-process the wambaugh data
    
    #rename the Clint column 
    httk_data <- httk_data %>% rename(CLint1 = Human.Clint.Point)
    
    #Remove CLint values of 0 and less
    httk_data$CLint1 <- ifelse(httk_data$CLint1==0 | httk_data$CLint1 < 0,
                               NA,
                               httk_data$CLint1)
    
  }
  
  if (chem_phys_in_vitro == T){ # pre process chem_phys_in_vitro data
    
    httk_data$CLint_iv <- ifelse(httk_data$CLint_iv==0 | httk_data$CLint_iv < 0,
                                 NA,
                                 httk_data$CLint_iv)
  }
  
  #find the columns which contain fraction unbound_data
  Clint_vector <- str_detect(colnames(httk_data),'CLint')
  Clint_vector_to_remove <- colnames(httk_data)[Clint_vector]
  clint_df <- httk_data[,Clint_vector]
  
  #convert to dataframe with only numeric entries
  if (is.null(ncol(clint_df))){
    
    clint_df <- suppressWarnings(as.numeric(clint_df))
  } else if (ncol(clint_df) > 0){
    
    clint_df <- suppressWarnings(as.data.frame(apply(clint_df, 2,function(x) as.numeric(x))))
  }
    
  if (is.null(ncol(clint_df))){
    
    #if there is only 1 source, take the data from it only
    httk_data$CLint <- clint_df
    
  } else if (ncol(clint_df) == 0){
    
    #if there is no CLint data, just populate with NA
    httk_data$CLint <- NA
    
  } else{
    
    #find mean of the 2 CLint values
    if (CLint_operation == 'geometric mean'){
      
      mean_log <- rowMeans(log(clint_df), na.rm = T)
      geom_mean <- exp(mean_log)
      httk_data$CLint <- geom_mean
      
      
    } else if (CLint_operation == 'arithmetic mean'){
      
      httk_data$CLint <- rowMeans(clint_df, na.rm = T)
      
    } else if (CLint_operation == 'None'){
      
      httk_data$CLint <- apply(clint_df,1,function(x) paste0(x[!is.na(x)], collapse = ' & '))
      
    }  else if (CLint_operation == 'minimum'){
      
      httk_data$CLint <- suppressWarnings(apply(clint_df,1,function(x) min(as.numeric(x),na.rm = T)))
      httk_data$CLint <- ifelse(httk_data$CLint=='Inf',NA, httk_data$CLint)
      
    } else if (CLint_operation == 'maximum'){
      
      httk_data$CLint <- suppressWarnings(apply(clint_df,1,function(x) max(as.numeric(x),na.rm = T) ))
      httk_data$CLint <- ifelse(httk_data$CLint=='Inf',NA, httk_data$CLint)
      
    } else{
      
      CLint_operation <- 'arithmetic mean'
      httk_data$CLint <- rowMeans(clint_df, na.rm = T)
      message('arithmetic mean applied. check the spelling of your input to the CLint_operation argument')
    }
    
  }
  
  httk_data$CLint <- suppressWarnings(as.numeric(httk_data$CLint))
  
  #calculate the standard deviation across available sources:
  httk_data$CLint_SD<-apply(clint_df,1,function(x) sd(x, na.rm=T))
  
  #set the units for CLint
  httk_data$CLint_units = "uL/min/million cells"
  
  #organise sources for Clint unbound
  potential_sources <- c('wambaugh2019_source|Clint_source_iv')
  available_sources <- str_detect(colnames(httk_data), potential_sources)
  available_sources <- httk_data[,available_sources]
  
  #replace empty space sources with NAs
  available_sources <- replace(available_sources,available_sources=='',NA)
  
  #check if sources are duplicated
  if (is.null(ncol(available_sources))){
    CLint_sources <- tolower(available_sources)
  } else if (ncol(available_sources) > 0){
    available_sources_upper <- as.data.frame(apply(available_sources,2,function(x) tolower(x))) #convert to upper case
    final_sources <- apply(available_sources_upper,1,function(x) unique(x[!is.na(x)]))
    source_num<- lengths(final_sources)
    
    CLint_sources <- c()
    for (k in 1:length(final_sources)){
      
      if(source_num[k]>1){
        
        data <- paste0(final_sources[[k]],collapse = ' & ')
        
      } else if (source_num[k]==0){
        
        data <- NA
        
      } else{
        
        data <- final_sources[[k]]
      }
      
      CLint_sources[k]<- data
      
    }
    
  } else {
    CLint_sources <- NA
  }

  
  httk_data$CLint_source <- CLint_sources
  for (i in 1:nrow(httk_data)){
    if (is.na(httk_data$CLint_source[i])){
      next
    }
    
    if (str_detect(httk_data$CLint_source[i],'&')){
      httk_data$CLint_source[i] <- paste('httk:', CLint_operation, 'of', httk_data$CLint_source[i], sep=' ')
      
    } else {
      httk_data$CLint_source[i] <- paste('httk:', httk_data$CLint_source[i], sep=' ')
      
    }
  }

  #remove unwanted columns
  if (wambaugh ==T | chem_phys_in_vitro == T){
    rm_cols <- c(Clint_vector_to_remove,'Clint_source_iv','wambaugh2019_source')
    httk_data <- rm_df_cols(httk_data, rm_cols)
  }
  
  #-----------------process sys clearance data------------------------------------
  
  if (obach == T){
    #organise the source of systemic clearance (which comes only from Obach data)
    httk_data <- httk_data %>% rename(Systemic_CL = `CL (mL/min/kg)`)
    if (nrow(Obach_data)>0){
      httk_data <- httk_data %>% rename(Systemic_CL_source = Obach2008_source)
      httk_data$Systemic_CL_source <- ifelse(is.na(httk_data$Systemic_CL),
                                             NA,
                                             paste('httk:',httk_data$Systemic_CL_source,sep=' '))
    }
    
    #convert clearance value from mL/min/kg to L/hr (multiply systemic CL by 0.06 and BW) for later iterations
    
    #set the units for CL
    httk_data$Systemic_CL_units = "mL/min/kg"
  } else{
    
    httk_data$Systemic_CL <- NA
    httk_data$Systemic_CL_source <- NA
    httk_data$Systemic_CL_units = "mL/min/kg"
    
  }
  
 #-----------------------sort out BP ratio if it is missing ---------------------------- 
  if (chem_phys_in_vitro != T){
    httk_data$BP <- NA
    httk_data$BP_source <- NA
  }
  
  
  #organise the httk_data columns
  httk_data <- httk_data %>% relocate(Systemic_CL_units, .after = Systemic_CL)
  httk_data <- httk_data %>% relocate(CLint_units, .after = CLint)
  httk_data <- httk_data %>% relocate(CLint_SD, .after = CLint)
  httk_data <- httk_data %>% relocate(CLint_source, .after = CLint_units)
  httk_data <- httk_data %>% rename(CLint_value = CLint)
  httk_data <- httk_data %>% rename(fu_value = fu_res)
  httk_data <- httk_data %>% rename(BP_value = BP)
  
  #merge the data collected from httk with the CAS rn and InchiKeys
  httk_data_inchi <- merge(httk_data, CAS_DTXSID, by = 'CAS', all.x = T)
  
  #remove the CAS column
  httk_data_inchi<- rm_df_cols(httk_data_inchi, c('CAS'))
  
  #populate additional columns if available
  #find available columns and store in a vector in order to merge with ChEMBL data
  potential_headers <- c('DOSE','DOSEUNITS','VSSMETHOD')
  additional_names<-keep_df_cols(info,potential_headers)
  if (length(additional_names)==0){
    additional_names <- c()
  } else{
    indicies <- which(potential_headers %in% colnames(info))
    additional_names <- potential_headers[indicies]
  }
  
  #merge the httk_data_inchi with the data curated from ChEMBL/sus dat
  physchem_httk <- merge(physchem_data, httk_data_inchi,
                         by.x = c('InChIKey','CODE',additional_names), 
                         by.y = c('StdInChIKey','CODE',additional_names),
                         all = T)
  
  if (chem_phys_in_vitro == T){
  
  #take compound names and smiles if unavailable
  physchem_httk$COMPOUND.NAME <- ifelse(is.na(physchem_httk$COMPOUND.NAME),
                               toupper(physchem_httk$Compound), physchem_httk$COMPOUND.NAME)
  physchem_httk$SMILES.x <- ifelse(is.na(physchem_httk$SMILES.x),
                                        physchem_httk$SMILES.y, physchem_httk$SMILES.x)
  
  
  #NOTE: Chembl httk values are favoured over httk, change these commands if you want otherwise
  
  #take MW from httk if unavilable
  physchem_httk$MW.x <- ifelse(is.na(physchem_httk$MW.x),
                               physchem_httk$MW.y, physchem_httk$MW.x)
  
  #take logP from httk if unavilable
  physchem_httk$logPow <- ifelse(is.na(physchem_httk$logPow),
                               physchem_httk$logP, physchem_httk$logPow)
  
  
  #if the source is not chembl, take acidic and basic pka values form httk
  chembl_mentions <- str_detect(tolower(physchem_httk$data_source),'chembl')
  chembl_mentions <- ifelse(is.na(chembl_mentions),FALSE,chembl_mentions)
  physchem_httk$Acidic..pKa <- ifelse(chembl_mentions | !is.na(physchem_httk$Acidic..pKa),
                                      physchem_httk$Acidic..pKa, physchem_httk$pKa_acidic)
  physchem_httk$Basic.pKa <- ifelse(chembl_mentions | !is.na(physchem_httk$Basic.pKa),
                                      physchem_httk$Basic.pKa, physchem_httk$pKa_basic)
  
  #organise pKa sources from httk
  httk_pKa_sources <- apply(physchem_httk[,c('pKa_acidic_ref','pKa_basic_ref')], 1, 
                      function(x) paste(x[!is.na(x)], collapse = " & "))
  httk_pKa_sources <- ifelse(httk_pKa_sources=='', 'EPI SUITE (EPA)', httk_pKa_sources)
  httk_pKa_sources <- paste('httk:', httk_pKa_sources, sep=' ')
  
  #if there are pKa values but the source is NA, add the httk_pKa_source
  physchem_httk$data_source <- ifelse(is.na(physchem_httk$data_source),
                                      httk_pKa_sources, physchem_httk$data_source)
  
  }
  
  if (wambaugh == T & chem_phys_in_vitro == T){
  
  #check that DTXSID match between invitro dataset and wambaugh
  unmatched <- c()
  for (i in 1:nrow(physchem_httk)){
    if(is.na(physchem_httk$DSSTox_Substance_Id[i])){
      next
    }

    if(physchem_httk$DSSTox_Substance_Id[i] %in% physchem_httk$DTXSID.y[i]){
      next
    } else{
      unmatched <- append(unmatched,i)
    }
  }

  if(!is.null(unmatched)){
    physchem_httk$CLint_value[unmatched] <- NA
    physchem_httk$CLint_source[unmatched] <- NA
    physchem_httk$fu_value[unmatched] <- NA
    physchem_httk$fu_source[unmatched] <- NA

  }
  
  }
  
  #remove unwanted columns
  rmv <- c("DTXSID.x","DSSTox_Substance_Id","logP","MW.y","pKa_basic_ref",
           "pKa_acidic_ref","pKa_basic","pKa_acidic",'SMILES.y','Compound')
  physchem_httk <- rm_df_cols(physchem_httk, rmv)

  #order based on compound codes
  physchem_httk <- physchem_httk[order(physchem_httk$CODE),]

  #oraganise data set
  if ('DOSE' %in% colnames(physchem_httk)){
    physchem_httk <- physchem_httk %>% relocate(DOSE, .after = ChEMBL.ID)
  }
  
  if (chem_phys_in_vitro == T){
  physchem_httk <- physchem_httk %>% relocate(DTXSID.y, .after = ChEMBL.ID)
  #rename columns
  physchem_httk <- physchem_httk %>% rename(SMILES = SMILES.x)
  physchem_httk <- physchem_httk %>% rename(DTXSID = DTXSID.y)
  physchem_httk <- physchem_httk %>% rename(MW = MW.x)
  }

  return(physchem_httk)
  
}

#function to find the average of numeric values separated by a separator
avg_string_of_nums <- function(val_string, separator){ #the input is a string of value(s)
  
  #if the input is not a single number
  if (suppressWarnings(is.na(as.numeric(val_string)))){ 
    
    separated_values <-unlist(str_split(val_string, separator))
    final_value<- mean(as.numeric(separated_values))
    
  }else{
    
    #if there is only a single number, return it
    final_value <- as.numeric(val_string)
  }
  
  return(final_value)
}

pKa_selection <- function(val_string, separator, selection){
  
  if(is.na(val_string)){
    return(NA)
  } else{
    
    separated_values <-unlist(str_split(val_string, separator))
    numeric_values <- as.numeric(separated_values)
    
    if (selection == 'smallest'){
      
      return(min(numeric_values))
      
    } else if (selection == 'biggest'){
      
      return(max(numeric_values))
      
    }
  }
  
}
