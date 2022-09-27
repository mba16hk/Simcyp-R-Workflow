#script to process and incorporate additional data to that queries from ChEMBL and susDat
#users can upload their own data here which can override the already collected physchem data

#function that returns list of compounds with missing info 
MissingInformation <- function(info, not_found_chembl,episuite_data, missing_info = T){
  
  if (missing_info == T){
    
    missing_PSA <- which(is.na(episuite_data$PSA))
    missing_HBD <- which(is.na(episuite_data$HBD))
    both_missing <- unique(c(missing_HBD,missing_PSA))
    
    not_found_compounds<-NotFoundInsusdat(not_found_chembl, episuite_data)
    not_found_compounds$Status <- 'Missing Compound'
    #get inchi of missing compounds
    
    missing_inchi<-episuite_data$InChIKey[both_missing]
    #missing_codes <- episuite_data$Code[both_missing]
    
    compounds_with_missing_info <- filter(info,info$INCHIKEY %in% missing_inchi)
    compounds_with_missing_info <- keep_df_cols(compounds_with_missing_info,
                                                c('CODE', 'SMILES'))
    compounds_with_missing_info$Status <- 'Missing Infomation'
    missing_and_nf_compounds<-rbind(compounds_with_missing_info, not_found_compounds)
    missing_and_nf_compounds <- missing_and_nf_compounds %>% relocate(Status, .before = CODE)
    
    return(missing_and_nf_compounds)
    
  } else {
    
    missing_compounds <- NotFoundInsusdat(not_found_chembl, episuite_data)
    missing_compounds$Status <- 'Missing Compound'
    missing_compounds <- missing_compounds %>% relocate(Status, .before = CODE)
    
    return(missing_compounds)
    
  }
  
}

#extract important data from additional data and incorporate with remaining physchem info
AdditionalData <- function (info, additional_data_directory, sus_data,
                            override_existing_data = F){
  
  
  #import the ACD labs datafile
  import<-loadWorkbook(additional_data_directory, create=FALSE)
  additional_data <- readWorksheet(import,  header = TRUE, sheet = 1)
  
  #make all the headers in uppercase
  headers<- toupper(colnames(additional_data))
  colnames(additional_data)<-headers
  
  #keep columns of interest
  keep <- 'CODE|LOGP|PKA(ACID)|PKA(BASE)|MW|HBD|PSA|TYPE'
  cols_to_keep <- str_detect(colnames(additional_data),keep)
  additional_data <- additional_data[,cols_to_keep]
  
  #if the compound code is not available, inform user
  if ('CODE' %!in% colnames(additional_data)){
    return(warning('Ensure a column header "CODE" is present. Additional data not incorporated.'))
  }
  
  #remove rows made of NAs entirely
  if (TRUE %in% apply(additional_data,1, function(x) all(is.na(x)))){
    rmv <- which(apply(additional_data,1, function(x) all(is.na(x)))==T)
    additional_data <- additional_data[-rmv,]
  }
  
  #incorporate data to the physiochemical information from chembl and sus_data
  x <- merge(sus_data, additional_data,
             by= 'CODE', all= T)
  
  #Do not overwrite values extracted by susDat/chembl and copy over data
  if (override_existing_data == F){
    
    if ('MW' %in% colnames(additional_data)){
      x$MW.x <- ifelse(!is.na(x$MW.y)&is.na(x$MW.x),
                       x$MW.y, x$MW.x)
    }
    
    if ('LOGP' %in% colnames(additional_data)){
      x$logPow <- ifelse(!is.na(x$LOGP)&is.na(x$logPow),
                         x$LOGP, x$logPow)
    }
    
    if('HBD' %in% colnames(additional_data)){
      x$HBD.x <- ifelse(!is.na(x$HBD.y)&is.na(x$HBD.x),
                        x$HBD.y, x$HBD.x)
    }
    
    if('PSA' %in% colnames(additional_data)){
      x$PSA.x <- ifelse(!is.na(x$PSA.y)&is.na(x$PSA.x),
                        x$PSA.y, x$PSA.x)
    }
    
    if('PKA(ACID)' %in% colnames(additional_data)){
      x$Acidic..pKa <- ifelse(!is.na(x$PKA(ACID))&is.na(x$Acidic..pKa),
                              x$PKA(ACID), x$Acidic..pKa)
    }
    
    if('PKA(BASE)' %in% colnames(additional_data)){
      x$Basic.pKa <- ifelse(!is.na(x$PKA(BASE))&is.na(x$Basic.pKa),
                            x$PKA(BASE), x$Basic.pKa)
    }
    
    if('TYPE' %in% colnames(additional_data)){
      x$Compound.type <- ifelse(!is.na(x$TYPE)&is.na(x$Compound.type),
                                x$TYPE, x$Compound.type)
    }
  } else if (override_existing_data == T){
    
    if ('MW' %in% colnames(additional_data)){
      x$MW.x <- ifelse(!is.na(x$MW.y), x$MW.y, x$MW.x)
    }
    
    if ('LOGP' %in% colnames(additional_data)){
      x$logPow <- ifelse(!is.na(x$LOGP), x$LOGP, x$logPow)
    }
    
    if('HBD' %in% colnames(additional_data)){
      x$HBD.x <- ifelse(!is.na(x$HBD.y), x$HBD.y, x$HBD.x)
    }
    
    if('PSA' %in% colnames(additional_data)){
      x$PSA.x <- ifelse(!is.na(x$PSA.y),x$PSA.y, x$PSA.x)
    }
    
    if('PKA(ACID)' %in% colnames(additional_data)){
      x$Acidic..pKa <- ifelse(!is.na(x$PKA(ACID)), x$PKA(ACID), x$Acidic..pKa)
    }
    
    if('PKA(BASE)' %in% colnames(additional_data)){
      x$Basic.pKa <- ifelse(!is.na(x$PKA(BASE)),x$PKA(BASE), x$Basic.pKa)
    }
    
    if('TYPE' %in% colnames(additional_data)){
      x$Compound.type <- ifelse(!is.na(x$TYPE),x$TYPE, x$Compound.type)
    }
  }

  
  rmv_cols <- c('MW.y','LOGP','HBD.y','PSA.y','TYPE','PKA(BASE)','PKA(ACID)')
  x<- rm_df_cols(x,rmv_cols)
  
  #rename columns in x
  x <- x %>% rename(MW = MW.x)
  x <- x %>% rename(PSA = PSA.x)
  x <- x %>% rename(HBD = HBD.x)
  
  #import the identifier values for newly found info
  missing_identifier1<-which(is.na(x$INCHIKEY))
  missing_identifier2<-which(is.na(x$SMILES))
  missing_identifier3<-which(is.na(x$COMPOUND.NAME))
  
  missing_identifiers <- unique(c(missing_identifier1, missing_identifier2,missing_identifier3))
  missing_compounds<-x$CODE[missing_identifiers]
  indicies <- which(info$CODE %in% missing_compounds)
  missing_information <- info[indicies,]
  
  y <- merge(x,missing_information, by=c('CODE'),all=T)
  
  y$SMILES.x <- ifelse(!is.na(y$SMILES.y)&is.na(y$SMILES.x),
                       y$SMILES.y, y$SMILES.x)
  y$COMPOUND.NAME <- ifelse(!is.na(y$Compound)&is.na(y$COMPOUND.NAME),
                            y$Compound, y$COMPOUND.NAME)
  y$InChIKey <- ifelse(!is.na(y$InChiKey)&is.na(y$InChIKey),
                       y$InChiKey, y$InChIKey)
  
  #handle the optional columns
  potential_headers <- c('DOSE','DOSEUNITS','VSSMETHOD')
  additional_names<-keep_df_cols(info,potential_headers)
  if (length(additional_names)==0){
    additional_names <- c()
  } else{
    indicies <- which(potential_headers %in% colnames(info))
    additional_names <- potential_headers[indicies]
  }
  
  #each additional name will have two suffixes (.x and .y), merge and rename
  
  new_names.x<-c();new_names.y<-c() #create separate names for .x and .y suffixes
  if(!is.null(additional_names)){
    for (i in 1:length(additional_names)){
      new_names.x <-c(new_names.x,paste0(additional_names[i],'.x'))
      new_names.y <- c(new_names.y,paste0(additional_names[i],'.y'))
    }
  }
  
  #carry over data from .y to .x for each of the additional names
  for (i in 1:length(additional_names)){
    x_data <- y[,new_names.x[i]]
    y_data <- y[,new_names.y[i]]
    xy_data<-data.frame(x_data,y_data)
    xy_data$x_data<-ifelse(is.na(xy_data$x_data),xy_data$y_data, xy_data$x_data)
    y[,new_names.x[i]]<-xy_data$x_data
    y<-rm_df_cols(y,new_names.y[i])
    new_name <- word('dose.x',1,sep = "\\.")
    index<-which(colnames(y)=='DOSE.x')
    colnames(y)[index]<-toupper(new_name)
  }
  
  rmv_cols<- c('InChiKey','SMILES.y','Compound')
  y<-rm_df_cols(y,rmv_cols)
  
  #rename columns in x
  y <- y %>% rename(SMILES = SMILES.x)
  
  #add the source as ACD Labs for missing PSA/HBD
  y$data_source <- ifelse(is.na(y$data_source),"User's Additional Data",y$data_source)
  
  return(y)
  
}