#Additional Functions
library(stringr)
library(XLConnect)
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

#function to query Chembl using InChiKey
InChIKey_search  <- function(filenamedb,query)
{
  
  #dbListTables(ChEMBLDB)
  
  ## connect to ChEMBL SQL database
  ChEMBLDB <- dbConnect(dbDriver("SQLite"), dbname = filenamedb)
  
  ## identify molregno and canonical SMILES string in ChEMBLdb based on standard InChIkey
  InChiKey <- dbSendQuery(conn = ChEMBLDB, "SELECT molregno, standard_Inchi_key, canonical_smiles FROM compound_structures WHERE standard_InchI_key = ?")
  dbBind(InChiKey, list(query))
  STRUCTURES <- data.frame(dbFetch(InChiKey))
  dbClearResult(InChiKey)
  
  ## identify preferred name and CHEMBLID based on molregno 
  molregno <- dbSendQuery(conn = ChEMBLDB, "SELECT molregno, pref_name, chembl_id FROM molecule_dictionary WHERE molregno = ?")
  dbBind(molregno, list(STRUCTURES$molregno))
  NAMES <- data.frame(dbFetch(molregno))
  dbClearResult(molregno)
  
  ## identify MW (freebase), MW (full), cxlogp, molecular species, cx most acidic pka, cx most basic pka, psa, hbd  
  ## based on molregno 
  physchem <- dbSendQuery(conn = ChEMBLDB, "SELECT molregno, mw_freebase, full_mwt, cx_logp, cx_logd, molecular_species, cx_most_apka, cx_most_bpka,
                          psa, hbd, full_molformula FROM compound_properties WHERE molregno = ?")
  dbBind(physchem, list(STRUCTURES$molregno))
  PROPERTIES <- data.frame(dbFetch(physchem))
  dbClearResult(physchem)
  
  tmp <- merge(NAMES, STRUCTURES, by = "molregno", all = TRUE)
  results <- merge(tmp, PROPERTIES, by.x = "molregno", by.y = "molregno", all = TRUE)
  
  ## use the full vs freebase MW to check what is returned from search
  
  X <- c("Molregno","COMPOUND NAME", "ChEMBL ID", "InChIKey", "SMILES", "MWfreebase", "MW", 
         "logPow", 'logD',"Compound type","Acidic  pKa", "Basic pKa", "PSA", "HBD", "Molecular_Formula")
  colnames(results) <- X
  
  results <- transform(results, test = ifelse(MWfreebase == MW, "OK", "CHECK"))
  RSQLite::dbDisconnect(ChEMBLDB)
  return(results)
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

#Extract useful information from compound databases
extract_info <- function(info_to_extract){
  
  require(Simcyp)
  
  #find the number of simulated individuals
  individuals<- length(unique(GetAllIndividualValues_DB(IndividualValueID$IndividualNo,info_to_extract)))
  
  #extract AUC parameters for all individuals
  AUC_data<-c()
  for (i in 1:individuals){
    dat<-GetAUCFrom_DB(ProfileID$Csys,CompoundID$Substrate,individual = i,info_to_extract)
    AUC_data <- rbind(AUC_data,dat)
  }
  
  #Extract information from dB
  BSA<- GetAllIndividualValues_DB(IndividualValueID$BSA,info_to_extract)# BSA (m^2)
  Age<- GetAllIndividualValues_DB(IndividualValueID$Age,info_to_extract)# Age (years)
  BW<- GetAllIndividualValues_DB(IndividualValueID$BW,info_to_extract)# BW (kg)
  GFR<- GetAllIndividualValues_DB(IndividualValueID$GFR,info_to_extract) #mL/min/1.73m^2
  
  #predicted values
  Vss<- GetAllCompoundResults_DB('idPredictedVss',compound = CompoundID$Substrate ,info_to_extract) #L/kg
  Fg<- GetAllCompoundResults_DB('idfGut', compound = CompoundID$Substrate,info_to_extract)
  Fh<- GetAllCompoundResults_DB('idfLiver', compound = CompoundID$Substrate,info_to_extract)
  Fa<- GetAllCompoundResults_DB('idfaAdj',compound = CompoundID$Substrate, info_to_extract)
  Fu_plasma <- GetAllCompoundResults_DB('idfuAdj', compound = CompoundID$Substrate, info_to_extract)
  Ka<- GetAllCompoundResults_DB('idkaAdj',compound = CompoundID$Substrate, info_to_extract) #absorption rate constant (1/h)
  BP <- GetAllCompoundResults_DB('idbpAdj', compound = CompoundID$Substrate, info_to_extract) # BP ratio
  
  data_from_dB<-cbind(AUC_data,BSA,Age,BW,GFR,Vss,Fg,Fh,Fa,Fu_plasma,Ka,BP)
  
  rmv <- c('ProfileIndex','Inhibition','DiffStoreIndex','Dose',
           'StartTime','EndTime','Tmin','Cmin', 'AUCt_full',
           'Cfirst', 'Clast','LambdaZ')
  
  data_from_dB<- rm_df_cols(data_from_dB,rmv)
  
  return(data_from_dB)
}

# Function from plotting conc-time profiles for each compound
plot_profile <- function(casestudy_ID, Output, units = 'ng/mL', curated_data, logy=F){
  
  require(ggplot2)
  require(scales)
  
  #extract a df for the CSID
  indicies<-which(Output$Group == casestudy_ID)
  df <- Output[indicies,2:ncol(Output)]
  
  #if more than 30 individuals are simulated, only select 30 at random
  if (ncol(df)>31){
    randomly_selected_cols <- sample(2:ncol(df), 2, replace=F)
    df<-df[,c(1,randomly_selected_cols)]
  }
  
  #convert dataframe into long form
  df<- long_form(df)
  colnames(df)<-c('Sub','Times','Conc')
  
  #Create Title for plot
  header<-paste(casestudy_ID,'Concentration-Time Profiles',sep=' ')
  
  if (units == 'ng/mL'){
    
    df <- df
    #mean_profile$Concentration <- mean_profile$Concentration
    ytitle <- paste('Concentration',paste0('(',units,')'), sep=' ')
    
  } else if (units == 'uM'){
    
    #Extract MW from curated data
    df_curated_data <- curated_data[which(curated_data$Code == casestudy_ID),]
    MW <- df_curated_data$MW
    
    ## convert units of cmax from mg/L -> uM #check again
    df$Conc <- (df$Conc/1000)/MW*1e6
    ytitle <- paste('Concentration',paste0('(',units,')'), sep=' ')
    
  } else{
    
    message('Please select either ng/mL or uM. ng/mL used as default')
    ytitle <- paste('Concentration (ng/mL)', sep=' ')
  }
  
  if (logy == T){
    
    #Plot the mean profile (y axis normal scale)
    g1 <- ggplot(df, aes(x = Times, y = Conc, group = factor(Sub))) + geom_line(aes(color=factor(Sub)), size = 1.5, linetype = 1)
    g1 <- g1 + geom_line() 
    g1 <- g1 + xlab("Time (hours)") + ylab(ytitle)
    g1 <- g1 + ggtitle(header) + 
      labs(color  = "Subjects", linetype = 1)+ theme(
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey"),
        axis.text = element_text(colour = "black",size=12),
        axis.title = element_text(colour = "black",size=14,face="bold")
    )
    
    g1 + scale_y_continuous(trans='log10')
    
    #determine the breaks
    
    #start by determining the maximum value of the y axis
    # max_val <- roundUp(max(df$Conc))
    # 
    # #is the maximum value is more than 2x larger than the 
    # #actual maximum value, reduce the maximum value by half
    # if (max(df$Conc)<(0.5*max_val)){
    #   max_val <- 0.5* max_val
    #   max_val2<-roundDown(max_val)
    #   break_set1<-c(max_val2,max_val2/10,max_val2/100,max_val2/1000)
    #   break_set2<-c(max_val, break_set1[1:2]/2)
    # } else{
    #   break_set1<-c(max_val,max_val/10,max_val/100,max_val/1000)
    #   break_set2<-break_set1[1:3]/2
    # }
    # g1 +  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
    #                          breaks=c(0, break_set1, break_set2))
    
    
   
  } else {
    
    #Plot the mean profile (y axis normal scale)
    g1 <- ggplot(df, aes(x = Times, y = Conc, group = factor(Sub))) + geom_line(aes(color=factor(Sub)), size = 1, linetype = 1) 
    g1 <- g1 + xlab("Time (hours)") + ylab(ytitle)
    g1 + ggtitle(header) + 
      labs(color  = "Subjects", linetype = 1)+ theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line.x = element_line(color="black", size = 1),
      axis.line.y = element_line(color="black", size = 1),
      panel.grid.major = element_line(color = "grey"),
      panel.grid.minor = element_line(color = "grey"),
      axis.text = element_text(colour = "black",size=12),
      axis.title = element_text(colour = "black",size=14,face="bold")
    )
  }
  
}

plot_parameters<-function(simcyp_outputs, Compound, 
                          plot_type = 'Distribution', x_variable, y_variable,
                          chosen_col){
  require(ggplot2)
  #Extract data relationg to compound of interest
  compound_IDs<- names(simcyp_outputs)
  index <- which(Compound==compound_IDs)
  outputs<- simcyp_outputs[[index]]
  
  #determine x and y variables
  x_col<- which(colnames(outputs)==x_variable)
  y_col <- which(colnames(outputs)== y_variable)
  
  x_axis_label <- paste(x_variable,determine_units(x_variable), sep=' ')
  
  if (plot_type == 'Distribution'){
    #Create Title for plot
    header<-paste(Compound,'Distribution of', x_variable, sep=' ')
    
    # Distributions are plotted as density plots (only x variable is considered)
    ggplot(outputs, aes(x=outputs[,x_col], colour= chosen_col, fill = chosen_col)) +
      geom_density(alpha=0.4)+ ggtitle(header)+ xlab(x_axis_label)+
      scale_color_manual(values= chosen_col) +
      scale_fill_manual(values= chosen_col) +
      theme(legend.position = "None",
            plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA),
            axis.line.x = element_line(color="black", size = 1),
            axis.line.y = element_line(color="black", size = 1),
            panel.grid.major = element_line(color = "grey"),
            panel.grid.minor = element_line(color = "grey"),
            axis.text = element_text(colour = "black",size=12),
            axis.title = element_text(colour = "black",size=14,face="bold"))
    
  } else if (plot_type == 'Relationship'){
    
    #create y axis label
    y_axis_label <- paste(y_variable, determine_units(y_variable), sep=' ')
    
    #find the correlation coefficient
    cor_coeff<- round(cor(outputs[,x_col],outputs[,y_col]),2)
    
    #Create Title for plot
    header<-paste(Compound,' Relationship between ', x_variable,' and ', y_variable,' (r = ',cor_coeff,')' ,sep='')
    
    #relationships are scatter diagrams where x and y variables are needed
    g1 <- ggplot(outputs, aes(x = outputs[,x_col], y = outputs[,y_col]))
    g1 <- g1 + geom_point(colour = chosen_col) + geom_smooth(method = "lm", se = FALSE, colour = chosen_col)  
    g1 <- g1 + xlab(x_axis_label) + ylab(y_axis_label)
    g1 + ggtitle(header) + theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line.x = element_line(color="black", size = 1),
      axis.line.y = element_line(color="black", size = 1),
      panel.grid.major = element_line(color = "grey"),
      panel.grid.minor = element_line(color = "grey"),
      axis.text = element_text(colour = "black",size=12),
      axis.title = element_text(colour = "black",size=14,face="bold")
    )
  }
  
}

#cross-compound comparison plot
compare_simulated_compound <- function(summary_simcyp, parameter, bar_order, bar_col){
  
  compound_codes<-summary_simcyp$CS_code
  parameters_to_plot <- keep_df_cols(summary_simcyp,parameter)
  
  df<-data.frame(compound_codes,parameters_to_plot)
  colnames(df)<-c('Compounds','parameter')
  
  if(bar_order=='ascending'){
    df <- df[order(df$parameter),]
    
  } else if (bar_order == 'compound_code'){
    df <- df[order(df$Compounds),]
  }
  
  header<-paste('Mean',parameter,'Variation Across Simulated Compounds',sep=' ')
  
  ggplot(df, aes(x=Compounds, y= parameter)) +
    ylab(parameter)+
    geom_bar(stat="identity", colour = bar_col, fill = bar_col)+
    ggtitle(header)+ 
    scale_x_discrete(limits = df$Compounds) +
                       theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line.x = element_line(color="black", size = 1),
      axis.line.y = element_line(color="black", size = 1),
      panel.grid.major = element_line(color = "grey"),
      panel.grid.minor = element_line(color = "grey"),
      axis.text = element_text(colour = "black",size=12),
      axis.title = element_text(colour = "black",size=14,face="bold"))
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
  if (parameter == 'Age'){
    return('(years)')
  } else if (parameter == 'BW'){
    return('(kg)')
  } else if (parameter == 'BSA'){
    return('(m^2)')
  } else if (parameter == 'Cmax'){
    return('(mg/L)')
  } else if (parameter == 'Tmax'){
    return('(hours)')
  } else if (parameter == 'GFR'){
    return('(mL/min/1.73m^2)')
  } else if (parameter == 'Vss'){
    return('(L/kg)')
  } else if (parameter == 'Ka'){
    return('(1/h)')
  } else if (parameter == 'AUC' | parameter == 'AUCinf'){
    return('(mg*h/L)')
  } else if (parameter == 'HalfLife'){
    return('(hours)')
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

Set_SimDuration <- function(Time) {   # Note this function works only for duration entered as 24 hour increments. 
  
  # Set End Day 
  EndDay<- (Time/24)+1
  SetParameter(SimulationParameterID$EndDay,CategoryID$SimulationData, CompoundID$Substrate, as.integer(EndDay))
  #message(paste('Setting End day to ',EndDay  , sep= ' '))
  
  # Set number of doses
  NumDose<- (EndDay-1)*2
  SetParameter(SimulationParameterID$CmpNumDoses1,CategoryID$SimulationData, CompoundID$Substrate, as.integer(NumDose))  # Number dose 2 idNumberDoses2
  SetParameter(SimulationParameterID$CmpNumDoses2,CategoryID$SimulationData, CompoundID$Substrate, as.integer(NumDose))  # Number dose 3 idNumberDoses3
  SetParameter(SimulationParameterID$CmpNumDoses3,CategoryID$SimulationData, CompoundID$Substrate, as.integer(NumDose))  # Number dose 4 idNumberDoses
  #message(paste('Setting number doses to ', NumDose , sep= ' '))
  
  Button<-  if (EndDay<= 2){
    
    SetParameter(SimulationParameterID$VariablePop,CategoryID$SimulationData, CompoundID$Substrate, FALSE)
    #message(paste('Setting Size button to False', sep= ' '))  
  } else {
    
    SetParameter(SimulationParameterID$VariablePop,CategoryID$SimulationData, CompoundID$Substrate, TRUE)
    #message(paste('Setting Size button to True', sep= ' ')) 
  }
  
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
  
  #populate with binary values (T for natched, F for mismatch)
  matched_mismatched[matches]<- TRUE
  matched_mismatched[mismatches]<- FALSE
  
  return(matched_mismatched)
  
}
