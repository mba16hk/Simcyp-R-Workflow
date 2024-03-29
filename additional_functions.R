#Additional Functions
library(stringr)
library(readxl)
library(tools)
library(dplyr)
library(RSQLite)

#Remove Columns by name
rm_df_cols <-function(df, column_names){
  
  #user inserts a vector of column names, ensure all string
  column_names <- paste(column_names)
  
  #remove the column indicies based on corresponding names
  df<-df[ , -which(names(df) %in% column_names)]
  
  #return the cleaned up dataframe
  return(df)
  
}

# The 'not in' function
'%!in%' <- function(x,y)!('%in%'(x,y))

#function that keeps certain columns of a df, specified by column names
keep_df_cols <-function(df, column_names){
  
  #user inserts a vector of column names, ensure all string
  column_names <- paste(column_names)
  
  #remove the column indicies based on corresponding names
  df<-df[ , -which(names(df) %!in% column_names)]
  
  #return the cleaned up dataframe
  return(df)
}

#Functions for Kilford Equation Calculations
calculate.fu_inc<- function(logPorD){
  
  #Use LogP when calculating for a base or neutral
  #Use logD when calculating for an acid or ampholyte
  #note that logP = logD for neutral compounds
  
  V_R<- 0.005 #taken from Kilford et al
  exponent<-(0.072*logPorD^2)+(0.067*logPorD)-1.126
  fu<- 1/(1+(125*V_R*10^(exponent)))
  return(fu)
}

calculate.logD<-function(logP, pKa){
  
  #calculate logD of a compound at pH 7.4
  #only needed for acids and ampholytes
  fraction<-1/(1+10^(7.4-pKa))
  logD<-log10(10^(logP)*fraction)
  return(logD)
}


#function to check if all columns are present
check_columns<- function(data){
  available_colnames<-colnames(data)
  target_colnames<- c("Code","Compound", "InChiKey", "CAS RegNo", 
                      "SMILES",'DOSE','DOSEUNITS','VSSMETHOD')
  
  if (all(target_colnames %in% available_colnames)){
    return(FALSE)
  }else{
    return(TRUE)
  }
  
}

#function to check if some columns are present
check_some_columns<- function(data){
  
  target_colnames<- c('DOSE','DOSEUNITS','VSSMETHOD')
  indicies <- which(colnames(data) %in% target_colnames)
  available_colnames <- colnames(data)[indicies]
  additional_data <- as.data.frame(data[,indicies])
  colnames(additional_data)<-available_colnames
  
  if (length(available_colnames)==0){
    return(0)
  }
  
  if (length(available_colnames)==1){
    if (available_colnames == 'DOSE'){
      return(1)
    } else if (available_colnames == 'DOSEUNITS'){
      return(2)
    } else if (available_colnames == 'VSSMETHOD'){
      return(3)
    }
  }
  
  if (length(available_colnames)==2){
    if (c('DOSE','DOSEUNITS') %in% available_colnames){
      return(4)
    } else if (c('DOSEUNITS','VSSMETHOD') %in% available_colnames){
      return(5)
    } else if (c('DOSE','VSSMETHOD') %in% available_colnames){
      return(6)
    }
  }
  
}

#function to determine the HTML line code
determine_line <- function(longstrings){
  if (word(longstrings,1,sep=" ")=="header2"){
    longstrings<- sub('^\\w+\\s', '', longstrings)
    return(h2(longstrings))
  } else if(word(longstrings,1,sep=" ")=="header3"){
    longstrings<- sub('^\\w+\\s', '', longstrings)
    return(h3(longstrings))
  } else {
    return(p(longstrings))
  }
}

#function to determine the units of the Simcyp output parameters
determine_units <- function(parameter){
  
  require(stringr)
  
  if (parameter == 'Age'){
    return('(years)')
  } else if (parameter == 'BW'){
    return('(kg)')
  } else if (parameter == 'BSA'){
    return('(m^2)')
  } else if (str_detect(parameter,'Cmax')){
    return('(mg/L)')
  } else if (parameter == 'Tmax'){
    return('(hr)')
  } else if (parameter == 'GFR'){
    return('(mL/min/1.73m^2)')
  } else if (parameter == 'Vss'){
    return('(L/kg)')
  } else if (parameter == 'Ka'){
    return('(1/hr)')
  } else if (parameter == 'AUC' | parameter == 'AUCinf'){
    return('(mg*hr/L)')
  } else if (parameter == 'HalfLife'){
    return('(hr)')
  } else if (str_detect(parameter,'CL')){
    return('(L/hr)')
  } else {
    return(' ')
  }
}

#functions to round up and down logarithmically (base 10)
roundUp <- function(x) 10^ceiling(log10(x))
roundDown <- function(x) 10^floor(log10(x))

#function to extract column within the GUI
extract_column <- function (df, column_index){
  return(df[,column_index])
}

#function to extract tissue types
extract_tissue <- function (df){
  
  #extract everything after the first underscore
  names <- unique(sub("^[^_]*_", "", colnames(df)))
  names <- names[!names %in% c('Group','Times')]
  return(names)
}

#function to convert wide matrices (dataframes) into long form
long_form <- function (df){
  
  #initialise empty long dataframe
  long_data<-data.frame()
  
  #change into long form
  for (i in 1:(ncol(df)-1)){
    sub_ID<-rep(i,nrow(df))
    df2<-cbind(sub_ID,df[,1],df[,(i+1)])
    long_data<-rbind(long_data,df2)
  }
  
  return(long_data)
}


#function to determine if cols in df match with cols of a different dataframe
match.df <- function(df1, df2, colnames_to_match){
  
  #keep the columns we want to compare
  df1 <- keep_df_cols(df1, colnames_to_match)
  df2 <- keep_df_cols(df2, colnames_to_match)
  
  #collapse each df into a single string vector (rows are pasted together)
  df_args <- c(df1, sep=" ")
  df1_collapsed <- do.call(paste, df_args)
  
  df_args <- c(df2, sep=" ")
  df2_collapsed <- do.call(paste, df_args)
  
  #identify matches
  matches <- which(df1_collapsed %in% df2_collapsed)
  
  #identify mismatches
  mismatches <- which(df1_collapsed %!in% df2_collapsed)
  
  #create a empty vector
  matched_mismatched <- rep(NA, nrow(df1))
  
  #populate with binary values (T for matched, F for mismatch)
  matched_mismatched[matches]<- TRUE
  matched_mismatched[mismatches]<- FALSE
  
  return(matched_mismatched)
  
}

#identify the out-of range PSA and HBD values
OutOfRange_Parameter <- function (collected_data){
  
  columns_of_interest <- keep_df_cols(collected_data,c('Code','CODE','MW','PSA','HBD','CXLogP'))
  
  # The mandatory out of ranges for PSA and HBD for Peff
  OutOfRange_PSA_idx <- which(columns_of_interest$PSA>154.4 | columns_of_interest$PSA<16.2)
  out_of_range_PSA <- ifelse(columns_of_interest$PSA>154.4 | columns_of_interest$PSA<16.2, 'PSA', NA)
  OutOfRange_HBD_idx <- which(columns_of_interest$HBD>5 | columns_of_interest$HBD<0)
  out_of_range_HBD <- ifelse(columns_of_interest$HBD>5 | columns_of_interest$HBD<0, 'HBD', NA)
  
  # The out of ranges for MW and LogP for Peff
  OutOfRange_MW_idx <- which(columns_of_interest$MW>455 | columns_of_interest$MW<60)
  out_of_range_MW <- ifelse(columns_of_interest$MW>455 | columns_of_interest$MW<60, 'MW', NA)
  OutOfRange_logP_idx <- which(columns_of_interest$CXLogP>4 | columns_of_interest$CXLogP<c(-3))
  out_of_range_logP <- ifelse(columns_of_interest$CXLogP>4 | columns_of_interest$CXLogP<c(-3), 'logP', NA)
  
  # The out of ranges for MW for fu prediction
  OutofRange_MW_fu_idx <- which(columns_of_interest$MW>825 | columns_of_interest$MW<75)
  out_of_range_MW_fu <- ifelse(columns_of_interest$MW>825 | columns_of_interest$MW<75, 'MW', NA)
  
  columns_of_interest <- cbind(columns_of_interest,out_of_range_PSA,out_of_range_HBD,
                                                     out_of_range_MW,out_of_range_logP,out_of_range_MW_fu)
  
  # out of range data frame
  out_of_range <- columns_of_interest[unique(c(OutOfRange_PSA_idx,OutOfRange_HBD_idx,
                                        OutOfRange_MW_idx,OutOfRange_logP_idx,OutofRange_MW_fu_idx)),]
  
  #out of range parameter
  out_of_range$Affected_Parameter <- NA
  out_of_range$Affected_Parameter <- ifelse(!is.na(out_of_range$out_of_range_PSA)|!is.na(out_of_range$out_of_range_HBD)|
                                            !is.na(out_of_range$out_of_range_MW)| !is.na(out_of_range$out_of_range_logP), 'Peff', out_of_range$Affected_Parameter)
  out_of_range$Affected_Parameter <- ifelse(!is.na(out_of_range$Affected_Parameter)& !is.na(out_of_range$out_of_range_MW_fu), 'Peff & fu', out_of_range$Affected_Parameter)
  out_of_range$Affected_Parameter <- ifelse(is.na(out_of_range$Affected_Parameter)& !is.na(out_of_range$out_of_range_MW_fu), 'fu', out_of_range$Affected_Parameter)
  
  
  # collapse additional columns to provide an indicator of reason of our of range
  out_of_range_indicator <- apply(out_of_range[,c('out_of_range_PSA','out_of_range_HBD',
                                                  'out_of_range_MW','out_of_range_logP','out_of_range_MW_fu')],
                                  1,function(x) unique(x[!is.na(x)]))
  
  out_of_range_indicator_names <- unlist(lapply(out_of_range_indicator,
                                         function(x) paste0(x,collapse = ' & ')))
  
  #add an indicator column
  out_of_range$Out_of_range_Parameter <- out_of_range_indicator_names
  
  #keep only columns of interest
  out_of_range <- rm_df_cols(out_of_range,c('out_of_range_PSA','out_of_range_HBD',
                                            'out_of_range_MW','out_of_range_logP', 'out_of_range_MW_fu'))
  
  return(out_of_range)
  
}
