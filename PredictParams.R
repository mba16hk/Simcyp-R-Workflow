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
  BP <- c()
  
  for (i in 1:nrow(Compound_dataframe)){
    Compound <- Compound_dataframe[i,]
    Set_parameters(Compound, trials = 1, subjects=1, Time=3, oliveoil_water_for_NLP = T)
    
    # Create the Database file name, with DB Extension
    DBfilename <- paste( Compound$CS_code,".db",sep="")
    #print('DB file created')
    
    # File path to save the database results to
    DBfilepath <- file.path(path_user, DBfilename) 
    #print('DB file saved')
    
    # Run simulation, suppress all console output 
    SetParameter("idSeedVariable", 3, 0, 4) #seed 0
    capture.output(Simcyp::Simulate(database=DBfilepath), file='NUL')
    
    #extract BP data
    info_to_extract <- RSQLite::dbConnect(SQLite(),DBfilepath)   #DBfilepath
    BP[i] <- GetAllCompoundResults_DB('idbpAdj', compound = CompoundID$Substrate, info_to_extract) # BP ratio
    
    #Detach database connection
    RSQLite::dbDisconnect(info_to_extract)
  
  }
  
  #Finished with the engine
  Simcyp::Uninitialise()
  
  simcyp_outputs <- AdditionalOutputs(Compound_dataframe)
  summary_simcyp <-SummaryOutputs(simcyp_outputs)
  
  Predicted_params <- StaticPredictedParameters(summary_simcyp)
  Predicted_params$`Predicted BP ratio` <- BP
  
  return(Predicted_params)
}

