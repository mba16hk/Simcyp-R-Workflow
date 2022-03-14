library(httk)

httkSearch <- function (ChEMBL_EPI_search, CAS_nums){
  
  #------------ Obach2008 -------------------------------------
  
  #Start by searching Obach for Clint & fu
  Obach_data<-filter(Obach2008, Obach2008$CAS %in% CAS_nums$CAS)
  Obach_data$source <- 'Obach2008'
  keep<- c('CAS','fu','CL (mL/min/kg)','source')
  Obach_data <- keep_df_cols(Obach_data,keep)
  
  #------------ wambaugh2019 -----------------------------------
  
  #next is wambaugh 2019 database which contains clint, fu, logP, pKa
  wambaugh2019_data<-filter(wambaugh2019, wambaugh2019$CAS %in% CAS_nums$CAS)
  wambaugh2019_data$source <- 'wambaugh2019'
  keep<- c('CAS','Human.Clint.Point','Human.Funbound.plasma.Point','source')
  wambaugh2019_data <- keep_df_cols(wambaugh2019_data,keep)
  
  #------------ chem.physical_and_invitro.data ------------------
  
  #use chem.physical_and_invitro.data  to get more information on CLint, fu, and BP
  in_vitro_data<-filter(chem.physical_and_invitro.data,
                        chem.physical_and_invitro.data$CAS %in% CAS_nums$CAS)
  
  #only keep human data
  human_in_vitro <- in_vitro_data[str_detect(in_vitro_data$All.Species,'Human'),]
  rm(in_vitro_data)
  
  #only keep columns of interest
  keep <- c('CAS','Human.Clint','Human.Clint.Reference',
            'Human.Funbound.plasma','Human.Funbound.plasma.Reference',
            'Human.Rblood2plasma','Human.Rblood2plasma.Reference')
  human_in_vitro <- keep_df_cols(human_in_vitro, keep)
  colnames(human_in_vitro) <- c('CAS','CLint1','CLint_source_iv','fu1','fu_source_iv','BP','BP_source')
  
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
  
  #remove rows made of NAs by finding indicies
  ind <- apply(human_in_vitro[,2:ncol(human_in_vitro)], 1, function(x) all(is.na(x)))
  indicies <- which(ind == TRUE)
  human_in_vitro <- human_in_vitro[-indicies,]
  
  # #-----pKa values processing --------
  # 
  # #if values are None, convert to NA
  # human_in_vitro$pKa_Accept<- ifelse(human_in_vitro$pKa_Accept == 'None',
  #                                    NA, human_in_vitro$pKa_Accept)
  # human_in_vitro$pKa_Donor<- ifelse(human_in_vitro$pKa_Donor == 'None',
  #                                    NA, human_in_vitro$pKa_Donor)
  # 
  # #remove the source if there is no value associated with it
  # human_in_vitro$pKa_Accept.Reference<- ifelse(is.na(human_in_vitro$pKa_Accept),
  #                                    NA, human_in_vitro$pKa_Accept.Reference)
  # human_in_vitro$pKa_Donor.Reference<- ifelse(is.na(human_in_vitro$pKa_Donor),
  #                                   NA, human_in_vitro$pKa_Donor.Reference)
  # 
  # 
  # #Choose correct pKa values using the pKa selection function
  # human_in_vitro$pKa_basic <- NA
  # human_in_vitro$pKa_acidic <- NA
  # for (i in 1:nrow(human_in_vitro)){
  #   
  #   human_in_vitro$pKa_basic[i] <- pKa_selection(human_in_vitro$pKa_Accept[i],
  #                                                separator = ',|;',
  #                                                selection = 'biggest')
  #   
  #   human_in_vitro$pKa_acidic[i] <- pKa_selection(human_in_vitro$pKa_Donor[i],
  #                                                separator = ',|;',
  #                                                selection = 'smallest')
  # }
  
  #merge all data into a single dataframe
  httk_data_1 <- merge(Obach_data,wambaugh2019_data,by='CAS',all = T)
  httk_data <- merge(httk_data_1,human_in_vitro, by = 'CAS', all = T)
  
  #-----------------process fu data------------------------------------

  #find mean of the 3 fraction unbound values
  httk_data$fu_avg <- rowMeans(data.frame(as.numeric(httk_data$fu.x),
                                          as.numeric(httk_data$Human.Funbound.plasma.Point),
                                          as.numeric(httk_data$fu.y)),
                               na.rm = T)

  #remove unwanted columns and reorganise dataframe
  rm_cols <- c('fu.x','Human.Funbound.plasma.Point','fu.y')
  httk_data <- rm_df_cols(httk_data,rm_cols)
  httk_data <- httk_data %>% relocate(fu_avg, .after = CAS)
  
  #organise sources for fraction unbound
  httk_data$fu_source <- NA
  for (i in 1:nrow(httk_data)){
    
    if(!is.na(httk_data$fu_avg[i])){
      if (is.na(httk_data$source.x[i]) & is.na(httk_data$source.y[i])){
        httk_data$fu_source[i] <- NA
      } else if (is.na(httk_data$source.x[i]) & !is.na(httk_data$source.y[i])){
        httk_data$fu_source[i] <- paste('httk:',httk_data$source.y[i], sep=' ')
      }else if (!is.na(httk_data$source.x[i]) & is.na(httk_data$source.y[i])){
        httk_data$fu_source[i] <- paste('httk:',httk_data$source.x[i], sep=' ')
      } else if (!is.na(httk_data$source.x[i]) & !is.na(httk_data$source.y[i])){
        httk_data$fu_source[i] <- paste('httk:',
                                        httk_data$source.x[i],'&',
                                        httk_data$source.y[i], sep=' ')
      }
      
    }else{
      httk_data$fu_source[i] <- NA
    }
  }
  
  httk_data$fu_source <- ifelse(is.na(httk_data$fu_source),
                                paste('httk:', httk_data$fu_source_iv,sep = ' '),
                                httk_data$fu_source)
  
  #check if any of the sources haven't been included
  list_of_sources <- tolower(str_extract(httk_data$fu_source_iv,'[^\\s]+'))
  for (i in 1:nrow(httk_data)){
    if (str_detect(tolower(httk_data$fu_source[i]),list_of_sources[i])){
      next
    } else{
      
      httk_data$fu_source[i]<-paste(httk_data$fu_source[i],
                                    '&',httk_data$fu_source_iv[i],
                                    sep=' ')
      
    }
  }
  
  #set the units for the fraction unbound
  httk_data$fu_units = "fraction"
  
  #organise the httk_data columns
  httk_data <- httk_data %>% relocate(fu_units, .after = fu_avg)
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
  
  #find mean of the 3 CLint values
  httk_data$CLint <- rowMeans(data.frame(as.numeric(httk_data$CLint1),
                                          as.numeric(httk_data$CLint_iv)),
                               na.rm = T)
  httk_data$CLint <- suppressWarnings(as.numeric(httk_data$CLint))
  
  #set the units for CLint
  httk_data$CLint_units = "uL/min/million cells"
  
  #remove unwanted columns
  rm_cols <- c('CLint1','CLint_iv','fu_source_iv')
  httk_data <- rm_df_cols(httk_data, rm_cols)
  
  #organise the source of CLint values
  httk_data <- httk_data %>% rename(CLint_source = CLint_source_iv)
  
  #check if any of the sources haven't been included
  list_of_sources <- tolower(gsub("[^a-zA-Z]", "", httk_data$source.y))
  for (i in 1:nrow(httk_data)){
    
    if(is.na(list_of_sources[i])){
      next
    }
    
    if (str_detect(tolower(httk_data$CLint_source[i]),list_of_sources[i])){
      next
    } else{
      
      httk_data$CLint_source[i]<-paste(httk_data$CLint_source[i],
                                    '&',httk_data$source.y[i],
                                    sep=' ')
      
    }
  }

  httk_data$CLint_source <- ifelse(is.na(httk_data$CLint_source),
                            NA,
                            paste('httk:',httk_data$CLint_source,sep=' '))
  
  #organise the source of systemic clearance
  httk_data <- httk_data %>% rename(Systemic_CL = `CL (mL/min/kg)`)
  httk_data <- httk_data %>% rename(Systemic_CL_source = source.x)
  httk_data$Systemic_CL_source <- ifelse(is.na(httk_data$Systemic_CL),
                                   NA,
                                   paste('httk:',httk_data$Systemic_CL_source,sep=' '))
  
  #set the units for CLint
  httk_data$Systemic_CL_units = "mL/min/kg"
  
  #organise the httk_data columns
  httk_data <- httk_data %>% relocate(Systemic_CL_units, .after = Systemic_CL)
  httk_data <- httk_data %>% relocate(CLint_units, .after = CLint)
  httk_data <- httk_data %>% relocate(CLint_source, .after = CLint_units)
  
  #merge the data collected from httk with the CAS rn and InchiKeys
  httk_data_inchi <- merge(httk_data, CAS_nums, by = 'CAS', all.x = T)
  
  #remove the CAS column
  httk_data_inchi<- rm_df_cols(httk_data_inchi, c('CAS','source.y'))
  
  #merge the httk_data_inchi with the data curated from ChEMBL/EPI Suite
  physchem_httk <- merge(ChEMBL_EPI_search, httk_data_inchi,
                         by.x = 'InChIKey', by.y = 'InChiKey',
                         all = T)

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
