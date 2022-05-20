#load the httk library
library(httk)

httkSearch <- function (physchem_data, CAS_DTXSID, info,
                        
                        #apply an operation on the collected fu and Clint values
                        fu_operation = 'arithmetic mean', CLint_operation = 'arithmetic mean'){
  
  if (nrow(CAS_DTXSID) == 0) {
    return(message('no data from httk can be found for any of the compounds'))
  }
  
  #------------ Obach2008 -------------------------------------
  
  #Start by searching Obach for Clint & fu using CAS numbers
  Obach_data<-filter(Obach2008, Obach2008$CAS %in% CAS_DTXSID$CAS)
  if (nrow(Obach_data)>0){
    Obach_data$Obach2008_source <- 'Obach 2008'
  }
  keep<- c('CAS','fu','CL (mL/min/kg)','Obach2008_source')
  Obach_data <- keep_df_cols(Obach_data,keep)
  
  #------------ wambaugh2019 -----------------------------------
  
  #next is wambaugh 2019 database which contains clint, fu, logP, pKa using CAS
  wambaugh2019_data<-filter(wambaugh2019, wambaugh2019$CAS %in% CAS_DTXSID$CAS)
  if (nrow(wambaugh2019_data)>0){
    wambaugh2019_data$wambaugh2019_source <- 'wambaugh2019'
  }
  keep<- c('CAS','Human.Clint.Point','Human.Funbound.plasma.Point',
           'DSSTox_Substance_Id','wambaugh2019_source')
  wambaugh2019_data <- keep_df_cols(wambaugh2019_data,keep)
  
  #------------ chem.physical_and_invitro.data ------------------
  
  #use chem.physical_and_invitro.data  to get more information on CLint, fu, and BP 
  # quary using using DTXSID and CAS
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
                                "pKa_Donor.Reference",'CLint1','CLint_source_iv','fu1','fu_source_iv','BP','BP_source')
  
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
  
  #merge all data into a single dataframe
  httk_data_1 <- merge(Obach_data,wambaugh2019_data,by='CAS',all = T)
  httk_data <- merge(httk_data_1,human_in_vitro, by = 'CAS', all = T)
  
  #-----------------process fu data------------------------------------

  #find mean/median of the 3 fraction unbound values
  if (fu_operation == 'arithmetic mean'){
    
    httk_data$fu_res <- rowMeans(data.frame(as.numeric(httk_data$fu.x),
                                            as.numeric(httk_data$Human.Funbound.plasma.Point),
                                            as.numeric(httk_data$fu.y)),
                                 na.rm = T)
    
  } else if (fu_operation == 'geometric mean'){
    
    df <- data.frame(as.numeric(httk_data$fu.x),
                     as.numeric(httk_data$Human.Funbound.plasma.Point),
                     as.numeric(httk_data$fu.y))
    mean_log <- rowMeans(log(df), na.rm = T)
    geom_mean <- exp(mean_log)
    httk_data$fu_res <- geom_mean
    
  } else if (fu_operation == 'median'){
    
    df <- httk_data[,c('fu.x','fu.y','Human.Funbound.plasma.Point')]
    
    httk_data$fu_res <- apply(df, 1,function(x) median(as.numeric(x),na.rm = T) )
    
  } else if (fu_operation == 'None'){
    
    httk_data$fu_res <- apply(httk_data[,c('fu.x','fu.y','Human.Funbound.plasma.Point')],1,function(x) paste0(x[!is.na(x)], collapse = ' & '))
    
  } else{
    
    httk_data$fu_res <- rowMeans(data.frame(as.numeric(httk_data$fu.x),
                                            as.numeric(httk_data$Human.Funbound.plasma.Point),
                                            as.numeric(httk_data$fu.y)),
                                 na.rm = T)
    message('arithmetic mean applied. check the spelling of your input to the fu_operation argument')
    
  }
  

  #reorganise dataframe
  httk_data <- httk_data %>% relocate(fu_res, .after = CAS)
  httk_data <- httk_data %>% relocate(DTXSID, .after = CAS)
  
  #organise sources for fraction unbound
  if (nrow(Obach_data) >0 & nrow(wambaugh2019_data) >0){
    fu_sources <- apply(httk_data[,c('Obach2008_source','wambaugh2019_source')], 1, 
                        function(x) paste(x[!is.na(x)], collapse = " & "))
    httk_data$fu_source <- fu_sources
    httk_data$fu_source <- ifelse(httk_data$fu_source=='',NA,httk_data$fu_source)
  }
  
  #check if any of the sources haven't been included
  httk_data$fu_source <- ifelse(is.na(httk_data$fu_source) & !is.na(httk_data$fu_source_iv),
                                httk_data$fu_source_iv, httk_data$fu_source)
  httk_data$fu_source <- ifelse(str_detect(tolower(httk_data$fu_source),tolower(str_extract(httk_data$fu_source_iv,'[^\\s]+'))),
                                httk_data$fu_source, paste(httk_data$fu_source,
                                                           '&',httk_data$fu_source_iv,
                                                           sep=' ') )
  httk_data$fu_source <- ifelse(!is.na(httk_data$fu_source), paste('httk:',httk_data$fu_source, sep=' '), httk_data$fu_source)
  
  #remove unwanted columns
  rm_cols <- c('fu.x','Human.Funbound.plasma.Point','fu.y','fu_source_iv')
  httk_data <- rm_df_cols(httk_data,rm_cols)
  
  #set the units for the fraction unbound
  httk_data$fu_units = "fraction"
  
  #organise the httk_data columns
  httk_data <- httk_data %>% relocate(fu_units, .after = fu_res)
  httk_data <- httk_data %>% relocate(fu_source, .after = fu_units)
  
  #-----------------process Clint data------------------------------------
  
  #rename the Clint column 
  httk_data <- httk_data %>% rename(CLint1 = Human.Clint.Point)
  
  #Remove CLint values of 0 and less
  httk_data$CLint1 <- ifelse(httk_data$CLint1==0 | httk_data$CLint1 < 0,
                            NA,
                            httk_data$CLint1)
  httk_data$CLint_iv <- ifelse(httk_data$CLint_iv==0 | httk_data$CLint_iv < 0,
                            NA,
                            httk_data$CLint_iv)
  
  #organise the httk_data columns
  httk_data <- httk_data %>% relocate(CLint_iv, .after = CLint1)
  
  #find mean of the 2 CLint values
  if (CLint_operation == 'geometric mean'){
    
    df <- data.frame(as.numeric(httk_data$CLint1),
                     as.numeric(httk_data$CLint_iv))
    mean_log <- rowMeans(log(df), na.rm = T)
    geom_mean <- exp(mean_log)
    httk_data$CLint <- geom_mean

    
  } else if (CLint_operation == 'arithmetic mean'){
    
    
    httk_data$CLint <- rowMeans(data.frame(as.numeric(httk_data$CLint1),
                                           as.numeric(httk_data$CLint_iv)),
                                na.rm = T)
    
  } else if (CLint_operation == 'None'){
    
    httk_data$CLint <- apply(httk_data[,c('CLint1','CLint_iv')],1,function(x) paste0(x[!is.na(x)], collapse = ' & '))
    
  } else{
    
    httk_data$CLint <- rowMeans(data.frame(as.numeric(httk_data$CLint1),
                                           as.numeric(httk_data$CLint_iv)),
                                na.rm = T)
    message('arithmetic mean applied. check the spelling of your input to the CLint_operation argument')
  }
  
  
  httk_data$CLint <- suppressWarnings(as.numeric(httk_data$CLint))
  
  #set the units for CLint
  httk_data$CLint_units = "uL/min/million cells"
  
  #sort out CLint sources
  httk_data$CLint_source <- NA
  httk_data$CLint_source <- ifelse(!is.na(httk_data$CLint1)& is.na(httk_data$CLint_iv),
                                   httk_data$wambaugh2019_source, httk_data$CLint_source)
  httk_data$CLint_source <- ifelse(is.na(httk_data$CLint1)& !is.na(httk_data$CLint_iv),
                                   httk_data$CLint_source_iv, httk_data$CLint_source)
  httk_data$CLint_source <- ifelse(!is.na(httk_data$CLint_source), paste('httk:',httk_data$CLint_source, sep=' '), httk_data$CLint_source)

  
  #remove unwanted columns
  rm_cols <- c('CLint1','CLint_iv','CLint_source_iv','wambaugh2019_source')
  httk_data <- rm_df_cols(httk_data, rm_cols)
  
  #-----------------process sys clearance data------------------------------------
  
  #organise the source of systemic clearance
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
  
  #organise the httk_data columns
  httk_data <- httk_data %>% relocate(Systemic_CL_units, .after = Systemic_CL)
  httk_data <- httk_data %>% relocate(CLint_units, .after = CLint)
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
                         by.x = c('InChIKey','Code',additional_names), 
                         by.y = c('StdInChIKey','Code',additional_names),
                         all = T)
  
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

  #remove unwanted columns
  rmv <- c("DTXSID.x","DSSTox_Substance_Id","logP","MW.y","pKa_basic_ref",
           "pKa_acidic_ref","pKa_basic","pKa_acidic",'SMILES.y','Compound')
  physchem_httk <- rm_df_cols(physchem_httk, rmv)

  #order based on compound codes
  physchem_httk <- physchem_httk[order(physchem_httk$Code),]

  #oraganise data set
  if ('DOSE' %in% colnames(physchem_httk)){
    physchem_httk <- physchem_httk %>% relocate(DOSE, .after = ChEMBL.ID)
  }
  physchem_httk <- physchem_httk %>% relocate(DTXSID.y, .after = ChEMBL.ID)

  #rename columns
  physchem_httk <- physchem_httk %>% rename(SMILES = SMILES.x)
  physchem_httk <- physchem_httk %>% rename(DTXSID = DTXSID.y)
  physchem_httk <- physchem_httk %>% rename(MW = MW.x)
  

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
