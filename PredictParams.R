#File for making predictions for fu, BP and Vss for compounds based on physchem properties
library(Simcyp)

PredictParameters <- function(Compound_dataframe){
  
  #Intialise the system files path
  Simcyp::Initialise("C:\\Program Files\\Simcyp Simulator V21\\Screens\\SystemFiles",
                     21,species = SpeciesID$Human, verbose = FALSE)
  
  #Set script to source file location
  path_user <-Simcyp::ScriptLocation()
  setwd(path_user)
  
  #create a list of workspaces in the directory
  SimcypWksz<-unlist(list.files(path_user, pattern="\\.wksz$",full.names=FALSE , recursive=F ))
  
  # Set workspace for MechKiM model
  capture.output(Simcyp::SetWorkspace(SimcypWksz[2]), file='NUL')
  
  Predicted_params<-data.frame()
  for (i in 1:nrow(Compound_dataframe)){
    Compound <- Compound_dataframe[i,]
    Set_parameters(Compound, trials = 1, subjects=1, Time=24)
    
    # Create the Database file name, with DB Extension
    DBfilename <- paste( Compound$CS_code,".db",sep="")
    # File path to save the database results to
    DBfilepath <- file.path(path_user, DBfilename) 
    
    # Run simulation, suppress all console output 
    SetParameter("idSeedVariable", 3, 0, 4) #seed 0
    capture.output(Simcyp::Simulate(database=DBfilepath), file='NUL')
   
  }
  
  #Finished with the engine
  Simcyp::Uninitialise()
  
  simcyp_outputs <- AdditionalOutputs(Compound_dataframe)
  summary_simcyp <-SummaryOutputs(simcyp_outputs)
  
  Predicted_params<- as.data.frame(cbind(summary_simcyp$CS_code,
                           summary_simcyp$Fu_plasma,
                           summary_simcyp$Vss))
  colnames(Predicted_params)<-c("Code",
                                "Predicted Fu",
                                "Predicted Vss")
  
  return(Predicted_params)
}