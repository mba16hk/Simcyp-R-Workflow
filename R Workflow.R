library(Simcyp)
library(tidyverse) 

Set_parameters <- function(Compound, trials = 1, subjects, Time, oliveoil_water_for_NLP, fu_GFR = T
                           ) {
  
  info <- Compound
  
  #specify compoind CS ID
  message( "--------------------------------------------")
  message( paste("Compound" , info$Code) )
  
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
      #message('No HBD or PSA values specified. 0 values set for both.')
    } else{

      Simcyp::SetCompoundParameter(CompoundParameterID$PSAValue ,CompoundID$Substrate, as.numeric(info$PSA))
      Simcyp::SetCompoundParameter(CompoundParameterID$HBDValue ,CompoundID$Substrate, as.numeric(info$HBD))
      #message(paste('Setting PSA to',as.numeric(info$PSA), 'and HBD to', as.numeric(info$HBD), sep = ' '))

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
    message(paste('Vss prediction method is set to Method', info$Prediction_method, sep= ' '))
    Simcyp::SetCompoundParameter(CompoundParameterID$UseOliveOilWaterSurrogateForNLP ,CompoundID$Substrate, TRUE)
  } else {
    Simcyp::SetCompoundParameter(CompoundParameterID$KpPredictionMethod,CompoundID$Substrate, as.integer(Method))
    Simcyp::SetCompoundParameter(CompoundParameterID$UseOliveOilWaterSurrogateForNLP ,CompoundID$Substrate, FALSE)
    message(paste('Vss prediction method is set to Method', info$Prediction_method, sep= ' '))
  }
  
  # Olive:oil water partition
  # --------------------------
  
  if (as.numeric(info$Prediction_method)==3 & oliveoil_water_for_NLP == T){
    Simcyp::SetCompoundParameter(CompoundParameterID$UseOliveOilWaterSurrogateForNLP ,CompoundID$Substrate, TRUE)
  } else if (as.numeric(info$Prediction_method)==3 & oliveoil_water_for_NLP == F){
    Simcyp::SetCompoundParameter(CompoundParameterID$UseOliveOilWaterSurrogateForNLP ,CompoundID$Substrate, FALSE)
  }
  
  # Fraction unbound (Fu) Settings.
  # -------------------------------
  
  if (is.na(info$fu_value) | info$fu_value==0) {
    
    #message("Fu value not specified. Predited Fu used instead.")
    
    # switch to pred , use 0L to switch back to user input fu
    Simcyp::SetCompoundParameter("fu",CompoundID$Substrate, 1L ) 
    #set the quaternary nitrogen flag
    Simcyp::SetCompoundParameter(CompoundParameterID$QuatNSwitch,CompoundID$Substrate, info$quat_N ) 
    # predicted fu from Simcyp
    fu<-Simcyp::GetCompoundParameter("idFu2",CompoundID$Substrate)  
    print(paste('pred fu',fu,sep=''))
    
  } else{
    message(paste("Setting User-input Fu =", info$fu_value, sep=""))
    #user input fu  
    Simcyp::SetCompoundParameter("idFu1",CompoundID$Substrate, as.numeric(info$fu_value) )  
    fu <- Simcyp::GetCompoundParameter("idFu1",CompoundID$Substrate)  
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
      #message(paste("Setting User-input BP =", info$BP_value, sep=' '))
    }
  }
  
  
  
  # RULES FOR: Compound type and pKa settings. 
  # -----------------------------------------------
  ##Diprotic_Acid 0L, Diprotic_Base 1L, Monoprotic_Acid 2L, Monoprotic_Base 3L, Neutral 4L, Ampholyte 5L
  
  if (info$Compound_type== "Monoprotic acid") { 
    
    Simcyp::SetCompoundParameter(CompoundParameterID$CompoundType,CompoundID$Substrate, 2L )
    Simcyp::SetCompoundParameter(CompoundParameterID$pKa1 ,CompoundID$Substrate, as.numeric(info$pKa1)  )
    #message("ACIDIC. THEREFORE, pKa VALUE IS  SET. ")
    
  } else if (info$Compound_type== "Neutral") {
    Simcyp::SetCompoundParameter(CompoundParameterID$CompoundType,CompoundID$Substrate, 4L )
    message("NEUTRAL. No pKa value.")
    
  } else if (info$Compound_type== "Monoprotic base"){
    Simcyp::SetCompoundParameter(CompoundParameterID$CompoundType,CompoundID$Substrate, 3L )
    # AND
    Simcyp::SetCompoundParameter(CompoundParameterID$pKa1 ,CompoundID$Substrate, as.numeric(info$pKa1))
    #message("BASIC. THEREFORE, pKa  SET.")
    
  } else if (info$Compound_type== "Diprotic base"){
    Simcyp::SetCompoundParameter(CompoundParameterID$CompoundType,CompoundID$Substrate, 1L )
    # AND
    Simcyp::SetCompoundParameter(CompoundParameterID$pKa1 ,CompoundID$Substrate, as.numeric(info$pKa1))
    Simcyp::SetCompoundParameter(CompoundParameterID$pKa2 ,CompoundID$Substrate, as.numeric(info$pKa2))
    #message("Diprotic Base THEREFORE, 2 pKa values SET.")
    
  } else if (info$Compound_type== "Diprotic acid"){
    Simcyp::SetCompoundParameter(CompoundParameterID$CompoundType,CompoundID$Substrate, 0L )
    # AND
    Simcyp::SetCompoundParameter(CompoundParameterID$pKa1 ,CompoundID$Substrate, as.numeric(info$pKa1))
    Simcyp::SetCompoundParameter(CompoundParameterID$pKa2 ,CompoundID$Substrate, as.numeric(info$pKa2))
    #message("Diprotic Acid THEREFORE, 2 pKa values SET.")
    
  }else { #Compound is an ampholyte, therefore 2 pKas required
    
    #in this case we assume compound is ampholytic
    Simcyp::SetCompoundParameter(CompoundParameterID$CompoundType,CompoundID$Substrate, 5L )
    # AND
    Simcyp::SetCompoundParameter(CompoundParameterID$pKa1 ,CompoundID$Substrate, as.numeric(info$pKa1))
    Simcyp::SetCompoundParameter(CompoundParameterID$pKa2 ,CompoundID$Substrate, as.numeric(info$pKa2))
    #message("AMPHOLYTE. THEREFORE, 2 pKa values SET.")
    
  }
  
  # Set CLint Value if it is not NA
  # -------------------------------
  if (!is.na(info$CLint_value)){
    Simcyp::SetCompoundParameter("idWOMC_HEP_Clint",CompoundID$Substrate, as.numeric(info$CLint_value))
    message("CLint value provided by the user.")
  } 
  else if(is.na(info$CLint_value) & fu_GFR == T){
    ## Calculate fu*GFR
    GFR = 9.06 #L/h
    fu_GFR <- as.numeric(fu)*GFR

    Simcyp::SetCompoundParameter("idCLRbase",CompoundID$Substrate, as.numeric(fu_GFR))
    message("fu*GFR calculated and used as Renal Clearance.")
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
  
  if (subjects == 'pop rep'){

    #SetParameter(SimulationParameterID$VirtualPopulationRadio,CategoryID$SimulationData, CompoundID$Substrate,0L)
    #SetParameter(SimulationParameterID$PopRepRadio,CategoryID$SimulationData, CompoundID$Substrate,1L)
    SetParameter(SimulationParameterID$UseAverageMan,CategoryID$SimulationData, CompoundID$Substrate,1L)
    message('Simulating Population Representative')

  } else{
  
  if (is.numeric(subjects)){
    SetParameter(SimulationParameterID$VirtualPopulationRadio,CategoryID$SimulationData, CompoundID$Substrate,1L)
    SetParameter(SimulationParameterID$PopRepRadio,CategoryID$SimulationData, CompoundID$Substrate,0L)
    # 2.1 Set the Number of Trials
    Simcyp::SetParameter(SimulationParameterID$Group,CategoryID$SimulationData, CompoundID$Substrate, as.numeric(trials)) # Trial num
    #print(paste('trials:',GetParameter(SimulationParameterID$Group,CategoryID$SimulationData, CompoundID$Substrate), sep = ' '))
    
    # 2.2 Set the Number of Subjects for each trial
    print(subjects)
    Simcyp::SetParameter(Simcyp::SimulationParameterID$Mempgroup,CategoryID$SimulationData, CompoundID$Substrate,as.integer(subjects))
    message(paste('Simulating',GetParameter(Simcyp::SimulationParameterID$Mempgroup,CategoryID$SimulationData, CompoundID$Substrate), 'subjects.', sep = ' '))
    
    GetParameter(SimulationParameterID$Pop,CategoryID$SimulationData, CompoundID$Substrate)   # Population Size. This has not changed yet, one needs to change this also to => trial*subject 
    SetParameter(SimulationParameterID$Pop,CategoryID$SimulationData, CompoundID$Substrate,as.integer(trials*subjects))
    #print(paste('product:',as.integer(trials*subjects), sep = ' '))
  }
    
  }
  
  #set the simulation duration
  SetParameter(SimulationParameterID$StudyDuration,CategoryID$SimulationData, CompoundID$Substrate, Time)  # idStudyDuration
  
  message('Simcyp Model Parametrisation Completed!')
  
}

SimulateWorkspace <- function (data, workspace, path_user, trials = 1, subjects, Time, seed, oliveoil_water_for_NLP, fu_GFR,
                               
                               #dermal parameters
                               dermal_area, formulation_thickness, formulation_density,
                               
                               #subject design parameters
                               MinAge, MaxAge, Prop_females,
                               
                               #multiple dosin parameters
                               multiple_dosing, Num_doses, dose_interval){
  
  # Initialize data frames to store results 
  Output<- NULL
  
  for (i in 1:nrow(data)) {
    
    #Empty vectors to save predictions
    Concs<-vector(); Times<-vector()
    Group<-vector(); Subject_ID<- vector()
    
    #Select 1 compound at a time
    Data<- data[i,]
    
    # Set workspace
    capture.output(Simcyp::SetWorkspace(workspace), file='NUL')
    
    # Set the parameters of the compound in the loaded workspace
    Set_parameters(Data, trials, subjects, Time, oliveoil_water_for_NLP, fu_GFR)
    
    # set additional parameters if workspace is dermal
    if (str_detect(workspace,'Dermal')){
      
      
      if (!is.na(dermal_area)){
        #allow user to specify dermal area in cm2
        Simcyp::SetCompoundParameter(CompoundParameterID$DermalArea, CompoundID$Substrate, dermal_area)
        user_dermal_area <- GetCompoundParameter(CompoundParameterID$DermalArea,CompoundID$Substrate)
      } else{
        user_dermal_area <- GetCompoundParameter(CompoundParameterID$DermalArea,CompoundID$Substrate)
      }
      
      message(paste('Dermal area of application is',user_dermal_area,'cm^2',sep=' '))
      
      if (!is.na(formulation_thickness)){
        #allow user to specify dermal formulation thickness in cm
        Simcyp::SetCompoundParameter(CompoundParameterID$DermalApplicationLayerThickness, CompoundID$Substrate, formulation_thickness)
        user_formulation_thickness <- GetCompoundParameter(CompoundParameterID$DermalApplicationLayerThickness,CompoundID$Substrate)
      }else{
        user_formulation_thickness <- GetCompoundParameter(CompoundParameterID$DermalApplicationLayerThickness,CompoundID$Substrate)
      }
      
      message(paste('Formulation thickness is',user_formulation_thickness,'cm',sep=' '))
      
      if (!is.na(formulation_density)){
        #allow user to specify dermal formulation density in g/mL
        Simcyp::SetCompoundParameter(CompoundParameterID$DermalApplicationLayerDensity, CompoundID$Substrate, formulation_density)
        user_formulation_density <- GetCompoundParameter(CompoundParameterID$DermalApplicationLayerDensity,CompoundID$Substrate)
      } else{
        user_formulation_density <- GetCompoundParameter(CompoundParameterID$DermalApplicationLayerDensity,CompoundID$Substrate)
      }
      
      message(paste('Formulation density is',user_formulation_density,'g/L',sep=' '))
 
    }
    
    #message('Not Dermal')
    
    # Set additional parameters if there are multiple doses
    if (multiple_dosing == T){
      
      #setting to 0L will switch to single dosing
      Simcyp::SetParameter(SimulationParameterID$CmpSingleOrMultiple0,CategoryID$SimulationData, CompoundID$Substrate,1L)
      
      if (!is.na(dose_interval)){
        #allow user to specify dosing interval
        Simcyp::SetCompoundParameter(CompoundParameterID$Dose_Interval, CompoundID$Substrate, dose_interval)
        interval_of_doses <- GetCompoundParameter(CompoundParameterID$Dose_Interval,CompoundID$Substrate)
      } else{
        interval_of_doses <- GetCompoundParameter(CompoundParameterID$Dose_Interval,CompoundID$Substrate)
      }
      
      message(paste('The dose interval is',interval_of_doses,'hours',sep=' '))
      
      if (!is.na(Num_doses)){
        #allow user to specify the number of doses, otherwise use the default
        Simcyp::SetParameter(SimulationParameterID$CmpNumDosesSwitch0,CategoryID$SimulationData, CompoundID$Substrate,1L)
        Simcyp::SetParameter(SimulationParameterID$CmpNumDoses0,CategoryID$SimulationData, CompoundID$Substrate,Num_doses)
        Number_of_doses <- GetParameter(SimulationParameterID$CmpNumDoses0,CategoryID$SimulationData, CompoundID$Substrate)
      } else{
        Number_of_doses <-GetParameter(SimulationParameterID$CmpNumDoses0,CategoryID$SimulationData, CompoundID$Substrate)
      }
      #message(paste('Number of doses is',Number_of_doses,sep=' '))
      
    }
    
    #message('Single dose')
    
    #set a seed unless specified otherwise
    if (seed == T){
      SetParameter("idSeedVariable", 3, 0, 4) #seed 0
    }
    
    #Set additional subject design parameters
    if (!is.na(MinAge)){
      #allow user to specify the minium age of the population
      Simcyp::SetParameter(SimulationParameterID$Age_Min0,CategoryID$SimulationData, CompoundID$Substrate, MinAge)
      Minimum_Age <- GetParameter(SimulationParameterID$Age_Min0,CategoryID$SimulationData, CompoundID$Substrate)
    } else{
      Minimum_Age <-GetParameter(SimulationParameterID$Age_Min0,CategoryID$SimulationData, CompoundID$Substrate)
    }

    if (!is.na(MaxAge)){
      #allow user to specify the maximum age of the population
      Simcyp::SetParameter(SimulationParameterID$Age_Max0,CategoryID$SimulationData, CompoundID$Substrate, MaxAge)
      Maximum_Age <- GetParameter(SimulationParameterID$Age_Max0,CategoryID$SimulationData, CompoundID$Substrate)
    } else{
      Maximum_Age <-GetParameter(SimulationParameterID$Age_Max0,CategoryID$SimulationData, CompoundID$Substrate)
    }
    message(paste('Population Age Range is from',Minimum_Age,'to', Maximum_Age,sep=' '))

    if (!is.na(Prop_females)){
      #allow user to specify theproportion of females in the population
      Simcyp::SetParameter(SimulationParameterID$ProbF0,CategoryID$SimulationData, CompoundID$Substrate, Prop_females)
      proportion_of_females <- GetParameter(SimulationParameterID$ProbF0,CategoryID$SimulationData, CompoundID$Substrate)
    } else{
      proportion_of_females <-GetParameter(SimulationParameterID$ProbF0,CategoryID$SimulationData, CompoundID$Substrate)
    }
    
    #message(paste('Proportion of females in population is',proportion_of_females,sep=' '))
    
    # File path to save the database results to
    DBfilepath <- file.path(path_user, paste(Data$Code,".db",sep="")) 
    print(DBfilepath)
    
    # Run simulation, suppress all console output 
    capture.output(Simcyp::Simulate(database=DBfilepath), file='NUL')  
    message( paste("Simulation Complete for ", Data$Code) )
    
    Output <- rbind(Output, TissueProfiles(Data$Code))
  }

  return(Output)
  
}

SimcypSimulation <- function (organised_data, trials = 1, subjects, Time, seed = T, fu_GFR,
                              
                              #dermal parameters
                              dermal_area = NA , formulation_thickness = NA, formulation_density = NA,
                              
                              #multiple dosin parameters
                              multiple_dosing = F, Num_doses = 2, dose_interval = 12,
                              
                              #minimum and maximum age and proportion of females
                              MinAge = NA , MaxAge = NA, Prop_females = NA,
                              
                              #additional params
                              oliveoil_water_for_NLP = T){
  
  #Intialise the system files path
  Simcyp::Initialise("C:\\Program Files\\Simcyp Simulator V23\\Screens\\SystemFiles",
                   23,species = SpeciesID$Human, verbose = FALSE)
  
  #Set script to source file location
  path_user <-Simcyp::ScriptLocation()
  setwd(path_user)
  
  #create a list of workspaces in the directory
  SimcypWksz<-unlist(list.files(path_user, pattern="\\.wksz$",full.names=FALSE , recursive=F ))
  
  
  # Set the workspace depending on administration route
  if (organised_data$Route[1] == 'Oral'){
    
    message('Using an Oral Workspace')
    
    oral_wkspace_indicies <- str_detect(SimcypWksz,'Oral')
    SimcypWksz = SimcypWksz[oral_wkspace_indicies]
    
  } else if (organised_data$Route[1] == 'IV Bolus'){
    
    message('Using an IV Bolus Workspace')
    
    #infused for a duration of 30 seconds
    bolus_wkspace_indicies <- str_detect(SimcypWksz,'Bolus')
    SimcypWksz = SimcypWksz[bolus_wkspace_indicies]
    
  } else if (organised_data$Route[1] == 'Dermal'){
    
    message('Using a Dermal Workspace')
    
    #applied to an area of 60 cm^2
    dermal_wkspace_indicies <- str_detect(SimcypWksz,'Dermal')
    SimcypWksz = SimcypWksz[dermal_wkspace_indicies]
    
  } else if (organised_data$Route[1] == 'WP8'){
    
    WP8_wkspace_indicies <- str_detect(SimcypWksz,'WP8')
    SimcypWksz = SimcypWksz[WP8_wkspace_indicies]
  }
  
  # Separating the data based on Clint value. 
  if (fu_GFR == T){
    Mechkim <- data.frame()
    NonMechkim <- organised_data
  } else{
    Mechkim <- organised_data %>% filter(is.na(CLint_value) )  
    NonMechkim <- organised_data %>% filter(!is.na(CLint_value) ) 
  }
 
  
  if (nrow(Mechkim != 0)){
    
    #print(paste('Simcyp Workspace:',SimcypWksz[1],sep=''))
    message('Activating MechKiM')
    simulated_data_MechKim <- SimulateWorkspace(Mechkim, SimcypWksz[1], 
                                                path_user, trials, subjects, Time, seed, fu_GFR = F,
                                                oliveoil_water_for_NLP,
                                                dermal_area, formulation_thickness, formulation_density,
                                                MinAge, MaxAge, Prop_females,
                                                multiple_dosing, Num_doses, dose_interval)
  }else{
    
    simulated_data_MechKim <- data.frame()
  }
  
  if (nrow(NonMechkim != 0)){
   # print(paste('Simcyp Workspace:',SimcypWksz[2],sep=''))
    message('Activating Non-MechKiM')
    print(SimcypWksz[2])
    simulated_data_NonMechKim <- SimulateWorkspace(NonMechkim, SimcypWksz[2], 
                                                   path_user, trials, subjects, Time, seed, fu_GFR,
                                                   oliveoil_water_for_NLP,
                                                   dermal_area, formulation_thickness, formulation_density,
                                                   MinAge, MaxAge, Prop_females,
                                                   multiple_dosing, Num_doses, dose_interval)
  }else{
    simulated_data_NonMechKim <- data.frame()
  }
  
  #Finished with the engine
  Simcyp::Uninitialise()
  
  #combine the outputs from mechkim and nonmechkim
  Outputs<-rbind(simulated_data_MechKim, simulated_data_NonMechKim)
  
  return(Outputs)
  
}
 
AdditionalOutputs <- function (organised_data,directory){
  
  require(RSQLite)
  
  if (!is.null(organised_data) | is.null(directory)){
    codes<-organised_data$Code
    simulated_data <-list()
    for (i in 1:length(codes)){
      info_to_extract <- RSQLite::dbConnect(SQLite(),paste(codes[i],'.db',sep=''))   #DBfilepath
      
      #extract main data
      simulated_data[[i]]<- extract_info(info_to_extract)
      
      #Detach database connection
      RSQLite::dbDisconnect(info_to_extract)  
    }
    
    names(simulated_data)<- codes
    
    return(simulated_data)
  } else{
    
    codes<- unlist(list.files(directory, pattern="\\.db$",full.names=TRUE , recursive=F ))
    simulated_data <-list()
    for (i in 1:length(codes)){
      info_to_extract <- RSQLite::dbConnect(SQLite(),codes[i])   #DBfilepath
      
      #extract main data
      simulated_data[[i]]<- extract_info(info_to_extract)
      
      #Detach database connection
      RSQLite::dbDisconnect(info_to_extract)  
    }
    
    names(simulated_data)<- codes
    
    return(simulated_data)
    
    
  }
  
  
  
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

#Extract useful information from compound databases
extract_info <- function(info_to_extract){
  
  tissues <- c('Brain','Heart','Gut','Lung','Kidney','Pancreas',
               'Liver','Spleen','Skin','Muscle','Adipose','Plasma')
  
  tissues <- toupper(tissues)
  profile_indicies <- c("281470715297801","281470748852233","281470732075017",
                   "281470799183881","281470765629449","281470933401609", 
                   "281470782406665", "281470849515529", "281470832738313",
                   "281470815961097","281470681743369","281474959933446")
  
  #find the number of simulated individuals
  individuals<- length(unique(GetAllIndividualValues_DB(IndividualValueID$IndividualNo,info_to_extract)))
  
  #extract AUC parameters for all individuals
  AUC_data<-c()
  for (i in 1:individuals){
    
    #get plasma cmax from the database files
    dat<-GetAUCFrom_DB(ProfileID$Csys,CompoundID$Substrate,individual = i,info_to_extract)
    
    ########## extract tissue-related Cmax values #######################
    
    #formulate SQL query
    query <- paste("SELECT Cmax FROM AUCData WHERE (Dose = 1 AND Individual =",i,"AND ProfileIndex = ?)", sep = ' ')
    Cmax_val <- dbSendQuery(conn = info_to_extract, query)
    dbBind(Cmax_val, list(profile_indicies))
    
    #extract the Cmax at the different tissue types
    Cmax_tissue <- data.frame(dbFetch(Cmax_val))
    dbClearResult(Cmax_val)
    
    #bind the cmax tissue values with the dat structure 
    dat <- unlist(append(dat,unlist(Cmax_tissue)))
    
    ######################################################################
    
    #get cmax in other tissues
    AUC_data <- rbind(AUC_data,dat)
  }
  
  colnames(AUC_data)[19:ncol(AUC_data)] <- paste0('Cmax_',tissues)
  AUC_data <- as.data.frame(AUC_data)
  
  #Extract information from dB
  BSA<- GetAllIndividualValues_DB(IndividualValueID$BSA,info_to_extract)# BSA (m^2)
  Age<- GetAllIndividualValues_DB(IndividualValueID$Age,info_to_extract)# Age (years)
  BW<- GetAllIndividualValues_DB(IndividualValueID$BW,info_to_extract)# BW (kg)
  GFR<- GetAllIndividualValues_DB(IndividualValueID$GFR,info_to_extract) #mL/min/1.73m^2
  
  #predicted values
  Fg<- GetAllCompoundResults_DB('idfGut', compound = CompoundID$Substrate,info_to_extract)
  Fh<- GetAllCompoundResults_DB('idfLiver', compound = CompoundID$Substrate,info_to_extract)
  Ka<- GetAllCompoundResults_DB('idkaAdj',compound = CompoundID$Substrate, info_to_extract) #absorption rate constant (1/h)
  BP <- GetAllCompoundResults_DB('idbpAdj', compound = CompoundID$Substrate, info_to_extract) # BP ratio
  CLtot <- GetAllCompoundResults_DB('idCLtot', compound = CompoundID$Substrate, info_to_extract) # Systemic Blood clearance (L/h)
  CLH <- GetAllCompoundResults_DB('idCLintH', compound = CompoundID$Substrate, info_to_extract) # Total hepatic Clint (L/h)
  CLR <- GetAllCompoundResults_DB('idCLintR', compound = CompoundID$Substrate, info_to_extract) # Total Renal Clint (L/h)
  
  data_from_dB<-cbind(AUC_data,BSA,Age,BW,GFR,Fg,Fh,Ka,BP,CLtot,CLH,CLR)
  
  
  rmv <- c('ProfileIndex','Inhibition','DiffStoreIndex','Dose',
           'StartTime','EndTime','Tmin','Cmin','Cmax', 'AUCt_full',
           'Cfirst', 'Clast','LambdaZ')
  
  data_from_dB<- rm_df_cols(data_from_dB,rmv)
  
  return(data_from_dB)
}

#return parameters which are not population-dependant
StaticPredictedParameters <- function(summary_simcyp){
  
  Predicted_params<-data.frame()
  kd_vals <- data.frame()
  
  #extract additional data values
  for (j in 1:length(summary_simcyp$CS_code)){
    DBfilename <- paste(summary_simcyp$CS_code[j],".db",sep="")
    
    ## connect to ChEMBL SQL database
    CompoundDB <- dbConnect(dbDriver("SQLite"), dbname = DBfilename)
    
    ## extract parameters of interest
    KD <- dbSendQuery(conn = CompoundDB, "SELECT PredictedVss, PredictedFu, PredictedFa, PredictedHSA_KD, PredictedAGP_KD FROM CompPredictedValues")
    static_output_data <- data.frame(dbFetch(KD))
    dbClearResult(KD)
    RSQLite::dbDisconnect(CompoundDB)
    kd_vals <- rbind(kd_vals,static_output_data)
  }
  
  
  Predicted_params<- as.data.frame(cbind(summary_simcyp$CS_code,kd_vals))
  colnames(Predicted_params)<-c("Code",
                                "Predicted Vss",
                                "Predicted Fu",
                                #"Predicted BP",
                                "Predicted Fa",
                                "Predicted HSA Kd",
                                "Predicted AGP Kd")
  
  Predicted_params$`Predicted HSA Kd`[Predicted_params$`Predicted HSA Kd`=="Inf"]<- '-'
  Predicted_params$`Predicted AGP Kd`[Predicted_params$`Predicted AGP Kd`=="Inf"]<- '-'
 
  
  return(Predicted_params)
  
}

TissueProfiles <- function(compound_code){
  
  # Initialize empty output list where each list element is a dataframe of profile for a tissue type
  Tissue_profiles<- list()
  
  tissues <- c('Brain','Heart','Gut','Lung','Kidney','Pancreas',
               'Liver','Spleen','Skin','Muscle','Adipose','Plasma')
  
  tissues <- toupper(tissues)
  profile_ids <- c("281470715297801","281470748852233","281470732075017",
                   "281470799183881","281470765629449","281470933401609", 
                   "281470782406665", "281470849515529", "281470832738313",
                   "281470815961097","281470681743369","281474959933446")
  
  Output = NULL
  
  for (i in 1:length(tissues)){
    
    #Empty vectors to save predictions
    Concs<-vector(); Times<-vector()
    Group<-vector(); Subject_ID<- vector()
    
    conn <- RSQLite::dbConnect(SQLite(),paste(compound_code,'.db',sep=''))
    # Get the subject number 
    Ind_num<- GetAllIndividualValues_DB("idIndividualNo",conn)
    
    for (j in 1:length(Ind_num)){
      #tissue conc profile
      tissue_conc <- GetProfile_DB(profile_ids[i], individual = j, conn, 
                                   inhibition = FALSE)
      if (i==1){
        #time range
        Times <- GetProfile_DB("281474976645120", individual=j, 
                               conn, inhibition = FALSE)
      }
      
      Concs  <- cbind(Concs, tissue_conc) 
      Subject_ID<- c(Subject_ID, paste( "Subject",j,'_',tissues[i],  sep=""))
    }
    
    colnames(Concs) <- c(Subject_ID)  # renaming columns
    
    if (i==1){
      Group  <- rep(compound_code,length(Times))
      Output<- data.frame(Group, Times, Concs)
    } else{
      Output<- cbind(Output, data.frame(Concs))
    }
    
    RSQLite::dbDisconnect(conn) #Detach database connection
  }
  
  return(Output)
    
}
