
OrganiseInputData <- function (httk_exp_data, info,
                               
                               #things the user can chnage
                               Vss_method = 3, Input_Dose = 100, UNITS = 'mg/kg',
                               
                               #only oral administration allowed for now
                               admin_route = 'Oral', 
                               
                               #can set the BP value for acids if they are missing
                               BP_acids = 0.55,
                               
                               #can change the threshold of base pka to assume they bind to AGP
                               AGP_pKa_threshold = 7){
  
  # -------------- Moleclar Weight -----------------------#
  
  #ensure MW and MWfreebase are equal when one of them is NA
  httk_exp_data$MWfreebase <- ifelse(is.na(httk_exp_data$MWfreebase),
                                     httk_exp_data$MW,
                                     httk_exp_data$MWfreebase)
  
  # -------------- Compound characterisation and pKa values -----------------------#
  
  #Compounds with logP==logD are set to NEUTRAL only if they have no 'compound type' label
  httk_exp_data$Compound.type <- ifelse(is.na(httk_exp_data$Compound.type) & 
                                        (httk_exp_data$logPow == httk_exp_data$logD),
                                        'NEUTRAL',
                                        httk_exp_data$Compound.type)
  
  #If acidic pKa and Basic pKa are both NA ... set the compound to NEUTRAL
  httk_exp_data$Compound.type <- ifelse(is.na(httk_exp_data$Compound.type) & 
                                          is.na(httk_exp_data$Acidic..pKa) & is.na(httk_exp_data$Basic.pKa),
                                        'NEUTRAL',
                                        httk_exp_data$Compound.type)
  
  #if any NA compound types, set them to NEUTRAL
  httk_exp_data$Compound.type <- ifelse(is.na(httk_exp_data$Compound.type),'NEUTRAL',httk_exp_data$Compound.type)
  
  # If the compound is from EPI Suite and it has a pKa (in acidic pKa which we merged there),
  # then if pKa <7 compound is acidic, if >7 compound is basic
  httk_exp_data$Compound.type <- ifelse(is.na(httk_exp_data$Compound.type) & 
                                          httk_exp_data$data_source == 'EPI SUITE (EPA) - Experimental Data' & 
                                          !is.na(httk_exp_data$Acidic..pKa) & httk_exp_data$Acidic..pKa <7,
                                        'ACID', httk_exp_data$Compound.type)
  
  httk_exp_data$Compound.type <- ifelse(is.na(httk_exp_data$Compound.type) & 
                                          httk_exp_data$data_source == 'EPI SUITE (EPA) - Experimental Data' & 
                                          !is.na(httk_exp_data$Acidic..pKa) & httk_exp_data$Acidic..pKa >7,
                                        'BASE', httk_exp_data$Compound.type)
  
  #if acidic pKa or basic pKa are negative, set them to 0
  httk_exp_data$Acidic..pKa <- ifelse(httk_exp_data$Acidic..pKa<0,
                                      0, httk_exp_data$Acidic..pKa)
  
  httk_exp_data$Basic.pKa <- ifelse(httk_exp_data$Basic.pKa<0,
                                      0, httk_exp_data$Basic.pKa)
  
  #for Simcyp, only 1 pKa value is used for acids/bases, and no pKa for NEUTRALS
  
  #remove all pKas for NEUTRALS
  httk_exp_data$Acidic..pKa <- ifelse(httk_exp_data$Compound.type == 'NEUTRAL',
                                      NA, httk_exp_data$Acidic..pKa)
  httk_exp_data$Basic.pKa <- ifelse(httk_exp_data$Compound.type == 'NEUTRAL',
                                      NA, httk_exp_data$Basic.pKa)
  
  #choose highest pKa for bases and lowest pKa for acids
  pKas<-data.frame(as.numeric(httk_exp_data$Acidic..pKa),
                   as.numeric(httk_exp_data$Basic.pKa))
  
  for (i in 1:nrow(httk_exp_data)){
    if(httk_exp_data$Compound.type[i]=='BASE'){
      
      httk_exp_data$Acidic..pKa[i]<- max(pKas[i,],na.rm=T)
      
    } else if(httk_exp_data$Compound.type[i]=='ACID'){
      
      httk_exp_data$Acidic..pKa[i]<- min(pKas[i,],na.rm=T)
      
    }
  }
  
  #-------------create final organised file --------------------------#

  ## build Simcyp compound import dataframe
  CMPD_IMPORT <- data.frame( "Compound_Name" = toupper(as.character(httk_exp_data$COMPOUND.NAME)))
  CMPD_IMPORT$CS_code <- httk_exp_data$Code
  
  #If the compound has no name, use its CS code
  CMPD_IMPORT$Compound_Name <- ifelse(is.na(CMPD_IMPORT$Compound_Name), 
                                      as.character(CMPD_IMPORT$CS_code), 
                                      as.character(CMPD_IMPORT$Compound_Name))
  
  #keep inchikey
  CMPD_IMPORT$InChiKey <- httk_exp_data$InChIKey
  
  #create a compound ID for each compound
  CMPD_IMPORT$Compound_ID <- paste(CMPD_IMPORT$Compound_Name, " ",
                                   "(",CMPD_IMPORT$ChEMBL_ID, "; ",
                                   CMPD_IMPORT$InChiKey,")", sep = "")
  #set route of administration
  CMPD_IMPORT$Route <- admin_route
  
  if ('DOSE' %!in% colnames(httk_exp_data)){
    CMPD_IMPORT$Dose <- Input_Dose
  } else {
    CMPD_IMPORT$Dose <- httk_exp_data$DOSE
  }
  
  if ('DOSEUNITS' %!in% colnames(httk_exp_data)){
    CMPD_IMPORT$Dose_Units <- UNITS
  } else {
    CMPD_IMPORT$Dose_Units <- httk_exp_data$DOSEUNITS
  }
  
  CMPD_IMPORT$MW <- httk_exp_data$MWfreebase
  CMPD_IMPORT$logPow <- httk_exp_data$logPow
  CMPD_IMPORT$cmpd_type <- httk_exp_data$Compound.type
  CMPD_IMPORT$Compound_type <- NA
  
  #convert to names as in Simcyp
  CMPD_IMPORT <- within(CMPD_IMPORT, {
    n <- cmpd_type == "NEUTRAL"
    b <- cmpd_type == "BASE"
    a <- cmpd_type == "ACID"
    z <- cmpd_type == "ZWITTERION"
    
    Compound_type[n] <- "Neutral"
    Compound_type[b] <- "Monoprotic base"
    Compound_type[a] <- "Monoprotic acid"
    Compound_type[z] <- "Ampholyte"
  })
  
  #consider the addition of diprotic acids and bases based on the distinct pKas
  
  CMPD_IMPORT <- rm_df_cols(CMPD_IMPORT,c('cmpd_type','z','a','b'))
  
  CMPD_IMPORT$acidicpKa <- httk_exp_data$Acidic..pKa
  CMPD_IMPORT$basicpKa <- httk_exp_data$Basic.pKa
  CMPD_IMPORT$pKa1 <- 0
  CMPD_IMPORT$pKa2 <- 0
  
  ## assume monoprotic in all cases EXCEPT ampholytes
  CMPD_IMPORT <- transform(CMPD_IMPORT, pKa1 = ifelse(Compound_type == "Monoprotic base", basicpKa, pKa1))
  CMPD_IMPORT <- transform(CMPD_IMPORT, pKa1 = ifelse(Compound_type == "Monoprotic acid", acidicpKa, pKa1))
  CMPD_IMPORT <- transform(CMPD_IMPORT, pKa1 = ifelse(Compound_type == "Ampholyte", acidicpKa, pKa1))
  CMPD_IMPORT <- transform(CMPD_IMPORT, pKa2 = ifelse(Compound_type == "Ampholyte", basicpKa, pKa2))
  
  CMPD_IMPORT <- rm_df_cols(CMPD_IMPORT,c('n','acidicpKa','basicpKa'))
  
  if (!is.null(httk_exp_data$BP_value)){
    CMPD_IMPORT$BP_value <- httk_exp_data$BP_value
    CMPD_IMPORT$BP_value <- ifelse(CMPD_IMPORT$Compound_type == "Monoprotic acid" &
                                     is.na(CMPD_IMPORT$BP_value), BP_acids, CMPD_IMPORT$BP_value)
    
    CMPD_IMPORT$BP_type<- ifelse(CMPD_IMPORT$BP_value != 0 &
                                   !is.na(CMPD_IMPORT$BP_value),"User input", "Predicted")
  }
  
  if (!is.null(httk_exp_data$fu_value)){
    
    CMPD_IMPORT$fu_value<-httk_exp_data$fu_value
    CMPD_IMPORT$fu_option <- ""
    CMPD_IMPORT <- transform(CMPD_IMPORT, 
                             fu_option = ifelse(!is.na(fu_value), "User input", "Predicted"))
    
  }
  
  CMPD_IMPORT$quatN <- "No" ## ASSUMED! NEEDS A SOLUTION OR ALTERNATIVE MODEL! 
  
  #HSA or AGP determination
  CMPD_IMPORT$HSA_AGP <- ifelse(CMPD_IMPORT$Compound_type == "Monoprotic base" 
                                & CMPD_IMPORT$pKa1 >= AGP_pKa_threshold, "AGP", "HSA")
  
  CMPD_IMPORT$Absorption_Model <- "First Order Model" 
  CMPD_IMPORT$Permeability_system <- "PSA/HBD"
  CMPD_IMPORT$PSA <- httk_exp_data$PSA 
  CMPD_IMPORT$HBD <- httk_exp_data$HBD
  
  CMPD_IMPORT$Distribution_model <- "Full PBPK"
  
  if ('VSSMETHOD' %!in% colnames(httk_exp_data)){
    CMPD_IMPORT$Prediction_method <- Vss_method
  } else {
    CMPD_IMPORT$Prediction_method <- httk_exp_data$VSSMETHOD
  }
 
  CMPD_IMPORT$CLR <- "0"
  CMPD_IMPORT$Systemic_CL <- httk_exp_data$Systemic_CL
  CMPD_IMPORT$molregno <- httk_exp_data$Molregno
  CMPD_IMPORT$CLint_sys <- "Hep"
  CMPD_IMPORT$CLint_hep <- httk_exp_data$CLint_value
  CMPD_IMPORT$CLint_hep<-as.numeric(CMPD_IMPORT$CLint_hep)
  CMPD_IMPORT <- CMPD_IMPORT %>% rename(CLint_value = CLint_hep)
  
  #ensure these are numeric
  CMPD_IMPORT$pKa1 <- as.numeric(CMPD_IMPORT$pKa1)
  CMPD_IMPORT$pKa2 <- as.numeric(CMPD_IMPORT$pKa2)
  CMPD_IMPORT$logPow <- as.numeric(CMPD_IMPORT$logPow)
  
  #conform to rules about min and max of pKa1 and pKa2 in Simcyp
  CMPD_IMPORT$pKa1 <- ifelse(CMPD_IMPORT$pKa1 < 0, 0,  CMPD_IMPORT$pKa1)
  CMPD_IMPORT$pKa2 <- ifelse(CMPD_IMPORT$pKa2 < 0, 0,  CMPD_IMPORT$pKa2)
  CMPD_IMPORT$pKa1 <- ifelse(CMPD_IMPORT$pKa1 > 14, 14,  CMPD_IMPORT$pKa1)
  CMPD_IMPORT$pKa2 <- ifelse(CMPD_IMPORT$pKa2 > 14, 14,  CMPD_IMPORT$pKa2)
  
  # --------------- Conduct LogD7.4 and fu_hep calculations using Kilford equations ---------
  
  #logD values from ChEMBL are NOT used, we will rely on Kilford's calculations
  CMPD_IMPORT$logD_7.4<-''
  CMPD_IMPORT <- transform(CMPD_IMPORT, 
                           logD_7.4= ifelse(CMPD_IMPORT$Compound_type != 'Neutral' &
                                              CMPD_IMPORT$Compound_type != 'Ampholyte',
                                            calculate.logD(CMPD_IMPORT$logPow, CMPD_IMPORT$pKa1),
                                            NA))
  
  #calculate logD for ampholytes using the lowest pKa values
  ampholytes<-which(CMPD_IMPORT$Compound_type=='Ampholyte')
  pkas<-keep_df_cols(CMPD_IMPORT,c('pKa1','pKa2'))
  lowest_pka_vals<- apply(pkas[ampholytes,], 1, FUN = min)
  logD_ampholytes<-calculate.logD(CMPD_IMPORT$logPow[ampholytes],lowest_pka_vals)
  CMPD_IMPORT$logD_7.4[ampholytes]<-logD_ampholytes
  
  #calculate fraction unbound in hepatocytes
  CMPD_IMPORT$fu_inc<-''
  CMPD_IMPORT <- transform(CMPD_IMPORT, 
                           fu_inc = ifelse(Compound_type == "Monoprotic base" | Compound_type == "Neutral" | is.na(Compound_type), 
                                           calculate.fu_inc(CMPD_IMPORT$logPow), 
                                           calculate.fu_inc(CMPD_IMPORT$logD_7.4)))
  #organise in ascending CS number
  CMPD_IMPORT <- CMPD_IMPORT[order(CMPD_IMPORT$CS_code),]
  
  return(CMPD_IMPORT)
  
}
