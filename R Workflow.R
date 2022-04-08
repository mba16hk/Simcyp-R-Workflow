library(Simcyp)
library(tidyverse) 

Set_parameters <- function(Compound, trials = 1, subjects, Time) {
  
  info <- Compound
  
  #specify compoind CS ID
  message( "--------------------------------------------")
  message( paste("Compound" , info$CS_code) )
  
  # Set Parameters 
  # ---------------
  
  #Molecular weight
  Simcyp::SetCompoundParameter(CompoundParameterID$Mole_Wh , CompoundID$Substrate, as.numeric(info$MW))
  #Logpow
  Simcyp::SetCompoundParameter(CompoundParameterID$LogP ,CompoundID$Substrate, as.numeric(info$logPow))
  
  #See if it's a PSA/HBD System
  if (info$Permeability_system == 'PSA/HBD'){
    
    Simcyp::SetCompoundParameter(CompoundParameterID$PSARadio ,CompoundID$Substrate, 1L)
   
    if (is.na(info$PSA) & is.na(info$HBD)){
      #PSA
      Simcyp::SetCompoundParameter(CompoundParameterID$PSAValue ,CompoundID$Substrate, as.numeric(0))
      #HBD
      Simcyp::SetCompoundParameter(CompoundParameterID$HBDValue ,CompoundID$Substrate, as.numeric(0))
      message('No HBD or PSA values specified. 0 values set for both.')
    } else{
      
      Simcyp::SetCompoundParameter(CompoundParameterID$PSAValue ,CompoundID$Substrate, as.numeric(info$PSA))
      Simcyp::SetCompoundParameter(CompoundParameterID$HBDValue ,CompoundID$Substrate, as.numeric(info$HBD))
      message(paste('Setting PSA to',as.numeric(info$PSA), 'and HBD to', as.numeric(info$HBD), sep = ' '))
      
    }
    
    
  }
  
  #HSA/AGP
  if (info$HSA_AGP == 'HSA' | is.na(info$HSA_AGP)){
    SetCompoundParameter(Simcyp::CompoundParameterID$PlasmaBindingProtein,CompoundID$Substrate, as.integer(0))
    message(paste("Plasma Binding = HSA. Protein concentration reference value is 45 g/L ", sep="") )
    SetCompoundParameter(Simcyp::CompoundParameterID$ProteinReference,CompoundID$Substrate, as.numeric(45))
    
  } else if (info$HSA_AGP == 'AGP'){
    
    SetCompoundParameter(Simcyp::CompoundParameterID$PlasmaBindingProtein,CompoundID$Substrate, as.integer(1))
    message(paste("Plasma Binding = AGP. Protein concentration reference value is 0.811 g/L ", sep="") )
    SetCompoundParameter(Simcyp::CompoundParameterID$ProteinReference,CompoundID$Substrate, as.numeric(0.811))
  }
  
  #Vss Prediction Method
  #---------------------
  
  #Method1= 0L, Method2=1L, Method 3= 2L
  Method= as.numeric(info$Prediction_method)-1 # R setting is one less.
  
  if (as.numeric(info$Prediction_method)==3){
    Simcyp::SetCompoundParameter(CompoundParameterID$KpPredictionMethod,CompoundID$Substrate, as.integer(Method))
    Simcyp::SetCompoundParameter(CompoundParameterID$UseOliveOilWaterSurrogateForNLP ,CompoundID$Substrate, TRUE)
    message(paste('Vss prediction method is set to Method', info$Prediction_method, sep= ' '))
  } else {
    Simcyp::SetCompoundParameter(CompoundParameterID$KpPredictionMethod,CompoundID$Substrate, as.integer(Method))
    Simcyp::SetCompoundParameter(CompoundParameterID$UseOliveOilWaterSurrogateForNLP ,CompoundID$Substrate, FALSE)
    message(paste('Vss prediction method is set to Method', info$Prediction_method, sep= ' '))
  }
  
  # Fraction unbound (Fu) Settings.
  # -------------------------------
  
  if (is.na(info$fu_value) | info$fu_value==0) {
    
    message("Fu value not specified. Predited Fu used instead.")
    
    # switch to pred , use 0L to switch back to user input fu
    Simcyp::SetCompoundParameter("fu",CompoundID$Substrate, 1L ) 
    # predicted fu from Simcyp
    Pred_fu<-Simcyp::GetCompoundParameter("idFu2",CompoundID$Substrate)  
    
  } else{
    message(paste("Setting User-input Fu =", info$fu_value, sep=""))
    #user input fu  
    Simcyp::SetCompoundParameter("idFu1",CompoundID$Substrate, as.numeric(info$fu_value) )  
  }
  
  # Blood Plasms (BP) Settings.
  # -------------------------------
  
  if (!is.null(info$BP_value)){
    if (info$BP_type=='Predicted'){
      #switch to pred BP
      Simcyp::SetCompoundParameter("BP",CompoundID$Substrate, 1L )  
      message("BP not specified by user. Predicted BP used instead.")
      Pred_bp<-Simcyp::GetCompoundParameter("BP",CompoundID$Substrate)

    } else {
      #insert user-input BP
      Simcyp::SetCompoundParameter(CompoundParameterID$bp,CompoundID$Substrate, as.numeric(info$BP_value))
      message(paste("Setting User-input BP =", info$BP_value, sep=' '))
    }
  }
  
  
  
  # RULES FOR: Compound type and pKa settings. 
  # -----------------------------------------------
  ##Diprotic_Acid 0L, Diprotic_Base 1L, Monoprotic_Acid 2L, Monoprotic_Base 3L, Neutral 4L, Ampholyte 5L
  
  if (info$Compound_type== "Monoprotic acid") { 
    
    Simcyp::SetCompoundParameter(CompoundParameterID$CompoundType,CompoundID$Substrate, 2L )
    Simcyp::SetCompoundParameter(CompoundParameterID$pKa1 ,CompoundID$Substrate, as.numeric(info$pKa1)  )
    message("ACIDIC. THEREFORE, pKa VALUE IS  SET. ")
    
  } else if (info$Compound_type== "Neutral") {
    Simcyp::SetCompoundParameter(CompoundParameterID$CompoundType,CompoundID$Substrate, 4L )
    message("NEUTRAL. No pKa value.")
    
  } else if (info$Compound_type== "Monoprotic base"){
    Simcyp::SetCompoundParameter(CompoundParameterID$CompoundType,CompoundID$Substrate, 3L )
    # AND
    Simcyp::SetCompoundParameter(CompoundParameterID$pKa1 ,CompoundID$Substrate, as.numeric(info$pKa1))
    message("BASIC. THEREFORE, pKa  SET.")
    
  } else { #Compound could be ampholyte or a diprotic acid/base, therefore 2 pKas required
    
    #in this case we assume compound is ampholytic
    Simcyp::SetCompoundParameter(CompoundParameterID$CompoundType,CompoundID$Substrate, 5L )
    # AND
    Simcyp::SetCompoundParameter(CompoundParameterID$pKa1 ,CompoundID$Substrate, as.numeric(info$pKa1))
    Simcyp::SetCompoundParameter(CompoundParameterID$pKa2 ,CompoundID$Substrate, as.numeric(info$pKa2))
    message("AMPHOLYTE. THEREFORE, 2 pKa values SET.")
    
  }
  
  # Set CLint Value if it is not NA
  # -------------------------------
  if (!is.na(info$CLint_value)){
    Simcyp::SetCompoundParameter("idEK_HEP_Clint",CompoundID$Substrate, as.numeric(info$CLint_value))
    message("CLint value provided by the user.")
  }
  
  # Input the calculated fu_inc value
  # ----------------------------------
  Simcyp::SetCompoundParameter(CompoundParameterID$WOMC_HepatocyteFUinc1, 
                       CompoundID$Substrate, as.numeric(info$fu_inc))

  
  #Set the input dose and their respective Units
  #---------------------------------------------
  Simcyp::SetCompoundParameter(CompoundParameterID$Dose, CompoundID$Substrate, as.numeric(info$Dose))
  
  if (info$Dose_Units == "mg/kg"){
    Simcyp::SetCompoundParameter(CompoundParameterID$DoseType, CompoundID$Substrate, 2L)
    Units <- GetCompoundParameter(CompoundParameterID$DoseType,CompoundID$Substrate)
  } else if (info$Dose_Units == "mg"){
    Simcyp::SetCompoundParameter(CompoundParameterID$DoseType, CompoundID$Substrate, 1L)
    Units <- GetCompoundParameter(CompoundParameterID$DoseType,CompoundID$Substrate)
  } else if (info$Dose_Units == "mg/m^2"){
    Simcyp::SetCompoundParameter(CompoundParameterID$DoseType, CompoundID$Substrate, 0L)
    Units <- GetCompoundParameter(CompoundParameterID$DoseType,CompoundID$Substrate)
  }
  
  unit_names<-c( "mg/m^2", "mg", "mg/kg")
  message(paste('Setting dose to', as.numeric(info$Dose),unit_names[Units+1], sep=' '))
  
  # 2.1 Set the Number of Trials
  Simcyp::SetParameter(SimulationParameterID$Group,CategoryID$SimulationData, CompoundID$Substrate, as.integer(trials)) # Trial num
  #print(paste('trials:',GetParameter(SimulationParameterID$Group,CategoryID$SimulationData, CompoundID$Substrate), sep = ' '))

  # 2.2 Set the Number of Subjects for each trial
  Simcyp::SetParameter(Simcyp::SimulationParameterID$Mempgroup,CategoryID$SimulationData, CompoundID$Substrate,as.integer(subjects))
  message(paste('Simulating',GetParameter(Simcyp::SimulationParameterID$Mempgroup,CategoryID$SimulationData, CompoundID$Substrate), 'subjects.', sep = ' '))
  
  GetParameter(SimulationParameterID$Pop,CategoryID$SimulationData, CompoundID$Substrate)   # Population Size. This has not changed yet, one needs to change this also to => trial*subject 
  SetParameter(SimulationParameterID$Pop,CategoryID$SimulationData, CompoundID$Substrate,as.integer(trials*subjects))
  #print(paste('product:',as.integer(trials*subjects), sep = ' '))
  
  #set the simulation duration
  Set_SimDuration(Time)
  
}

SimulateWorkspace <- function (data, workspace, path_user, trials = 1, subjects, Time){
  
  # Initialize data frames to store results 
  Output<- NULL
  
  # START workflow using Predicted Clint (MechKiM Model)
  #simulated_data<-c()
  
  for (i in 1:nrow(data)) { #nrow(Mechkim)
    
    #Empty vectors to save predictions
    Concs<-vector(); Times<-vector()
    Group<-vector(); Subject_ID<- vector()
    
    #Select 1 compound at a time
    Data<- data[i,]
    
    # Set workspace
    capture.output(Simcyp::SetWorkspace(workspace), file='NUL')
    
    # Set the parameters of the compound in the loaded workspace
    Set_parameters(Data, trials, subjects, Time)
    
    #Create a separate db for each compound
    DBfilename <- paste( Data$CS_code,".db",sep="")  # Create the Database file name, with DB Extension
    
    DBfilepath <- file.path(path_user, DBfilename) # File path to save the database results to
    
    # Run simulation, suppress all console output 
    SetParameter("idSeedVariable", 3, 0, 4) #seed 0
    capture.output(Simcyp::Simulate(database=DBfilepath), file='NUL')  
    message( paste("Simulation Complete for ", Data$CS_code) )
    
    conn <- RSQLite::dbConnect(SQLite(),DBfilename)   #DBfilepath
    #simulated_data<-rbind(simulated_data,extract_info(conn))
    
    # Get the subject number 
    Ind_num<- GetAllIndividualValues_DB("idIndividualNo",conn)
    #Group_num<- GetAllIndividualValues_DB("idGroupNo",conn)# Get the Trial number for the population
    
    No_sub<- length(Ind_num)
    
    #Extract Concentration time profile 
    for (j in 1:No_sub ){
      Current_C<- GetProfile_DB("281474959933446", individual=j, 
                                conn, inhibition = FALSE)  #Csys  Data
      Current_T<- GetProfile_DB("281474976645120", individual=j, 
                                conn, inhibition = FALSE)   #Time data
      
      Times  <- Current_T   
      Concs  <- cbind(Concs, Current_C) 
      Subject_ID<- c(Subject_ID, paste( "Subject" , j,  sep=""))
      
    }
    
    RSQLite::dbDisconnect(conn) #Detach database connection 
    
    colnames(Concs) <- c( Subject_ID )  # renaming columns 
    Group  <- rep(Data$CS_code,length(Times)) 
    Output<- rbind(Output, Output<- data.frame(Group, Times, Concs) )
    
  }
  
  return(Output)
  
}

SimcypSimulation <- function (organised_data, trials = 1, subjects, Time){
  
  #Intialise the system files path
  Simcyp::Initialise("C:\\Program Files\\Simcyp Simulator V21\\Screens\\SystemFiles",
                     21,species = SpeciesID$Human, verbose = FALSE)
  
  #Set script to source file location
  path_user <-Simcyp::ScriptLocation()
  setwd(path_user)
  
  #create a list of workspaces in the directory
  SimcypWksz<-unlist(list.files(path_user, pattern="\\.wksz$",full.names=FALSE , recursive=F ))
  
  if (organised_data$Route[1] == 'Oral'){
    
    oral_wkspace_indicies <- str_detect(SimcypWksz,'Oral')
    SimcypWksz = SimcypWksz[oral_wkspace_indicies]
    
  } else if (organised_data$Route[1] == 'IV Bolus'){
    
    #infused for a duration of 30 seconds
    bolus_wkspace_indicies <- str_detect(SimcypWksz,'Bolus')
    SimcypWksz = SimcypWksz[bolus_wkspace_indicies]
    
  } else if (organised_data$Route[1] == 'Dermal'){
    
    #applied to an area of 60 cm^2
    dermal_wkspace_indicies <- str_detect(SimcypWksz,'Dermal')
    SimcypWksz = SimcypWksz[dermal_wkspace_indicies]
    
  }
  
  # Separating the data based on Clint value. 
  Mechkim <- organised_data %>% filter(is.na(CLint_value) )  
  NonMechkim <- organised_data %>% filter(!is.na(CLint_value) ) 
  
  if (nrow(Mechkim != 0)){
    simulated_data_MechKim <- SimulateWorkspace(Mechkim, SimcypWksz[1], path_user, trials, subjects, Time)
  }else{
    simulated_data_MechKim <- data.frame()
  }
  
  if (nrow(NonMechkim != 0)){
    simulated_data_NonMechKim <- SimulateWorkspace(NonMechkim, SimcypWksz[2], path_user, trials, subjects, Time)
  }else{
    simulated_data_NonMechKim <- data.frame()
  }
  
  #Finished with the engine
  Simcyp::Uninitialise()
  
  #combine the outputs from mechkim and nonmechkim
  Outputs<-rbind(simulated_data_MechKim, simulated_data_NonMechKim)
  
  return(Outputs)
  
}
 
AdditionalOutputs <- function (organised_data){
  
  codes<-organised_data$CS_code
  simulated_data <-list()
  for (i in 1:length(codes)){
    info_to_extract <- RSQLite::dbConnect(SQLite(),paste(codes[i],'.db',sep=''))   #DBfilepath
    simulated_data[[i]]<-extract_info(info_to_extract)
    RSQLite::dbDisconnect(info_to_extract) #Detach database connection 
  }
  
  names(simulated_data)<- codes
  
  return(simulated_data)
  
}

SummaryOutputs <- function (simulated_data){
  #function that finds the means of all parameters for each compound
  
  codes <- names(simulated_data)
  summary_data<-c()
  for (i in 1:length(simulated_data)){
    
    #find the means
    mean_data <- colMeans(simulated_data[[i]])
    mean_data<- mean_data[-1] #remove individual number
    
    #find the standard deviations
    sd_data<-apply(simulated_data[[1]],2,sd) #calculate the sd of each column
    sd_data<- sd_data[-1] #remove individual number
    
    #add the term "SD" to each of the column headers
    sd_headers<- names(sd_data)
    sd_headers <- sub("","\\1SD ", sd_headers)
    names(sd_data)<- sd_headers
    
    #order the columns of the means and the sd
    neworder <- order(c(2*(seq_along(mean_data) - 1) + 1,
                        2*seq_along(sd_data)))
    
    mean_sd_data<-as.numeric(append(mean_data, sd_data)[neworder])
    mean_sd_data <- c(codes[i],round(mean_sd_data,3))
    
    summary_data<-as.data.frame(rbind(summary_data, mean_sd_data))
  }
  
  #adjust headers
  all_names<-c(names(mean_data),names(sd_data))
  headers<- all_names[neworder]
  colnames(summary_data)<- c('CS_code',headers)
  
  return(summary_data)
  
}
