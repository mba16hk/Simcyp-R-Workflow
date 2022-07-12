# Script containing functions related to ACD labs

#extract important data from ACD Labs outputs and incorporate with remaining physiochemical information
ACD_outputs <- function (info, ACD_data_directory, epi_data){
  
  #import the ACD labs datafile
  import<-loadWorkbook(ACD_data_directory, create=FALSE)
  acd_data <- readWorksheet(import,  header = TRUE, sheet = 1)
  
  #remove the first row
  if (NA %in% acd_data[1,]){
    acd_data <- acd_data[-1,]
  }
  
  #make all the headers in uppercase
  headers<- toupper(colnames(acd_data))
  colnames(acd_data)<-headers
  
  #keep columns of interest
  keep <- 'CODE|LOGP|PKA|MOLECULAR.WEIGHT|HYDROGEN.BOND.DONORS|TPSA'
  cols_to_keep <- str_detect(colnames(acd_data),keep)
  acd_data <- acd_data[,cols_to_keep]
  
  #rename ACD labs data cols (to merge with EPI Suite)
  headers<- colnames(acd_data)
  headers[str_detect(headers,'LOGP')]<-'logPow'
  headers[str_detect(headers,'MOLECULAR.WEIGHT')]<-'MW'
  headers[str_detect(headers,'HYDROGEN.BOND.DONORS')]<-'HBD'
  headers[str_detect(headers,'TPSA')]<-'PSA'
  headers[str_detect(headers,'CODE')]<-'Code'
  headers[str_detect(headers,'PKA')]<-'Acidic..pKa'
  #headers[str_detect(headers,'INCHI')]<-'InChIKey'
  colnames(acd_data) <- headers
  
  #step to process the pKas. If base take highest value, if acid take lowest
  # a base is determined if all pKa values are > 7 and acid is if all pKa < 7
  solid_rows<- which(!is.na(acd_data$Code))
  empty_rows_positions <- diff(solid_rows)
  
  cut_start<- which(empty_rows_positions!=1) #positions of where the cut in seq starts
  cut_length<- empty_rows_positions[cut_start]#length of the cut in seq
  cut_end <- solid_rows[cut_start] + (cut_length-1) #the number at the end of the cut
  
  #extracts the pka numbers between the NAs
  multiple_pKas<-list()
  for (i in 1:length(cut_start)){
    multiple_pKas[[i]]<-as.numeric(acd_data$Acidic..pKa[solid_rows[cut_start[i]]:cut_end[i]])
    acd_data$Acidic..pKa[solid_rows[cut_start[i]]:cut_end[i]] <- NA
  }
  
  pKa <- c()
  #determine only 1 pKa value to select
  for (i in 1:length(multiple_pKas)){
    
    
    if (all(multiple_pKas[[i]]>7)){
      
      #if basic, the pKa is the highest value
      pKa[i] <- max(multiple_pKas[[i]])
      
    } else if (all(multiple_pKas[[i]]<7)){
      
      #if all values are acidic, take the lowest
      pKa[i] <- min(multiple_pKas[[i]])
      
    } else {
      
      #if there are both basic and acidic values, 
      # take the one with the greatest difference from 7
      
      diff_from_7<-abs(7-multiple_pKas[[i]])
      pos_of_max <- which(diff_from_7 == max(diff_from_7))
      pKa[i]<-multiple_pKas[[i]][pos_of_max]
      
    }
    
  }
  
  #only keep the calculated pKas
  acd_data$Acidic..pKa[solid_rows[cut_start]]<-pKa
  
  #remove NA rows
  NA_rows <- apply(acd_data, 1, function(x) all(is.na(x)))
  acd_data <- acd_data[!NA_rows,]

  #incorporate data to the physiochemical information from chembl and epi suite
  x <- merge(epi_data, acd_data,
             by = c('Code'), all= T)
  
  #Do not overwrite values extracted by episuite/chembl and copy over data
  x$MW.x <- ifelse(!is.na(x$MW.y)&is.na(x$MW.x),
                   x$MW.y, x$MW.x)
  x$logPow.x <- ifelse(!is.na(x$logPow.y)&is.na(x$logPow.x),
                   x$logPow.y, x$logPow.x)
  x$HBD.x <- ifelse(!is.na(x$HBD.y)&is.na(x$HBD.x),
                       x$HBD.y, x$HBD.x)
  x$PSA.x <- ifelse(!is.na(x$PSA.y)&is.na(x$PSA.x),
                    x$PSA.y, x$PSA.x)
  x$Acidic..pKa.x <- ifelse(!is.na(x$Acidic..pKa.y)&is.na(x$Acidic..pKa.x),
                    x$Acidic..pKa.y, x$Acidic..pKa.x)
  rmv_cols <- c('MW.y','logPow.y','HBD.y','PSA.y','Acidic..pKa.y')
  x<- rm_df_cols(x,rmv_cols)
  
  #rename columns in x
  x <- x %>% rename(MW = MW.x)
  x <- x %>% rename(PSA = PSA.x)
  x <- x %>% rename(HBD = HBD.x)
  x <- x %>% rename(logPow = logPow.x)
  x <- x %>% rename(Acidic..pKa = Acidic..pKa.x)

  #import the identifier values for newly found info
  missing_identifier1<-which(is.na(x$InChIKey))
  missing_identifier2<-which(is.na(x$SMILES))
  
  missing_identifiers <- unique(c(missing_identifier1, missing_identifier2))
  missing_compounds<-x$Code[missing_identifiers]
  indicies <- which(info$Code %in% missing_compounds)
  missing_information <- info[indicies,]
  
  y <- merge(x,missing_information, by=c('Code'),all=T)
  
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
  y$data_source <- ifelse(is.na(y$data_source),'ACD/Labs Percepta',y$data_source)
  
  return(y)

}
