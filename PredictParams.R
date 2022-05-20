#File for making predictions for fu, BP and Vss for compounds based on physchem properties
library(Simcyp)

PredictParameters <- function(Compound_dataframe){
  
  #even if fu and BP values are provided, they will be ignored and predicted instead
  Compound_dataframe$fu_value <- NA
  Compound_dataframe$fu_option <- 'Predicted'
  Compound_dataframe$BP_value <- NA
  Compound_dataframe$BP_type <- 'Predicted'
  Compound_dataframe$CLint_value <- NA
  
  #Intialise the system files path
  Simcyp::Initialise("C:\\Program Files\\Simcyp Simulator V21\\Screens\\SystemFiles",
                     21,species = SpeciesID$Human, verbose = FALSE)
  
  #Set script to source file location
  path_user <-Simcyp::ScriptLocation()
  setwd(path_user)
  
  #create a list of workspaces in the directory
  SimcypWksz<-unlist(list.files(path_user, pattern="\\.wksz$",full.names=FALSE , recursive=F ))
  
  # Set workspace for non MechKiM model (faster)
  capture.output(Simcyp::SetWorkspace(SimcypWksz[6]), file='NUL')
  
  Predicted_params<-data.frame()
  kd_vals <- data.frame()
  
  for (i in 1:nrow(Compound_dataframe)){
    Compound <- Compound_dataframe[i,]
    Set_parameters(Compound, trials = 1, subjects=1, Time=3)
    
    # Create the Database file name, with DB Extension
    DBfilename <- paste( Compound$CS_code,".db",sep="")
    #print('DB file created')
    
    # File path to save the database results to
    DBfilepath <- file.path(path_user, DBfilename) 
    #print('DB file saved')
    
    # Run simulation, suppress all console output 
    SetParameter("idSeedVariable", 3, 0, 4) #seed 0
    #print('Seed is set')
    
    capture.output(Simcyp::Simulate(database=DBfilepath), file='NUL')
    #print('Simulation Complete')
    
   
  }
  
  #Finished with the engine
  Simcyp::Uninitialise()
  
  simcyp_outputs <- AdditionalOutputs(Compound_dataframe)
  summary_simcyp <-SummaryOutputs(simcyp_outputs)
  
  #extract Kd values
  for (j in 1:length(summary_simcyp$CS_code)){
    DBfilename <- paste(summary_simcyp$CS_code[j],".db",sep="")
    kd_vals <- rbind(kd_vals,ExtractKD(DBfilename))
  }
  
  
  Predicted_params<- as.data.frame(cbind(summary_simcyp$CS_code,
                           summary_simcyp$Fu_plasma,
                           summary_simcyp$Vss,
                           summary_simcyp$BP,
                           kd_vals))
  colnames(Predicted_params)<-c("Code",
                                "Predicted Fu",
                                "Predicted Vss",
                                "Predicted BP",
                                "Predicted HSA Kd",
                                "Predicted AGP Kd")
  
  return(Predicted_params)
}

#find predicted kd values (regardless if AGP or HSA)
ExtractKD <- function(DBfilename){
  
  ## connect to ChEMBL SQL database
  CompoundDB <- dbConnect(dbDriver("SQLite"), dbname = DBfilename)
  
  ## extract KD values
  KD <- dbSendQuery(conn = CompoundDB, "SELECT PredictedHSA_KD, PredictedAGP_KD FROM CompPredictedValues")
  Kd_vals <- data.frame(dbFetch(KD))
  dbClearResult(KD)
  RSQLite::dbDisconnect(CompoundDB)
  
  return(Kd_vals)
  
}
