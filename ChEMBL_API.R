# Install and load required library
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(httr)
library(stringr)
library(jsonlite)
library(dplyr)

# Function to search for a compound by InChIKey
ChEMBL_API_Search <- function(data_table) {
  
  inchikey<-data_table$INCHIKEY
  
  #split data into chunks of 20
  chunk_size <- 20
  
  # Calculate the number of chunks
  num_chunks <- ceiling(length(inchikey) / chunk_size)
  
  #Save output in large df
  output_dataframe<- data.frame()
  
  # Process the vector in chunks
  for (i in 1:num_chunks) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, length(inchikey))
    
    chunk <- inchikey[start_idx:end_idx]
    
    # Processing logic 
    
    # Construct the API request URL
    url <- paste("https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_structures__standard_inchi_key__in=",paste0(chunk,collapse=','),"&only=molecule_properties,molecule_structures",sep='')
    
    # Make the GET request
    response <- GET(url)
    
    # Check if the request was successful (status code 200)
    if (response$status_code == 200) {
      # Parse the JSON response
      output <- content(response, encoding = "UTF-8")
      output_molecules <- output$molecules
      first_list <- lapply(output_molecules, `[[`, 1)
      second_list <- lapply(output_molecules, `[[`, 2)
      
      # Combine the sublists into a single dataframe
      result_df_prop <- data.frame(do.call(rbind, first_list))
      result_df_structures <- data.frame(do.call(rbind, second_list))
      result_df <- cbind(result_df_structures,result_df_prop)
    } 
    
    #Get the ChEMBL IDs
    ##get ChEMBLID
    url <- paste("https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_structures__standard_inchi_key__in=",paste0(chunk,collapse=','),'&only=molecule_chembl_id',sep='')
    response <- GET(url)
    
    # Check if the request was successful (status code 200)
    if (response$status_code == 200) {
      # Parse the JSON response
      output <- content(response, encoding = "UTF-8")
      output <- output$molecules
      first_list <- lapply(output, `[[`, 1)
      
      # Combine the sublists into a single dataframe
      chembl_id_df <- data.frame(do.call(rbind, first_list))
      output_df <- cbind(chembl_id_df,result_df)
      colnames(output_df)[1] <- 'CHEMBL_id'
    } 
    
    output_dataframe <- rbind(output_dataframe,output_df)
  }
  
  #remove unwanted columns
  output_dataframe <- rm_df_cols(output_dataframe,c("molfile","standard_inchi",
                                                    "aromatic_rings","hba_lipinski","hbd_lipinski","heavy_atoms",
                                                    "mw_freebase","mw_monoisotopic","np_likeness_score",
                                                    "num_lipinski_ro5_violations","num_ro5_violations",
                                                    "qed_weighted","ro3_pass","rtb"))
  #convert NULL Molecular species into NA
  output_dataframe <- output_dataframe %>% replace(.=="NULL", NA) # replace with NA
  
  #check for structural alerts
  chembl_ids <- paste0(output_dataframe$CHEMBL_id,collapse = ',')
  url <- paste("https://www.ebi.ac.uk/chembl/api/data/compound_structural_alert.json?molecule_chembl_id__in=",chembl_ids,sep='')
  response <- GET(url)
  
  # Check if the request was successful (status code 200)
  if (response$status_code == 200) {
    # Parse the JSON response
    output <- content(response,'text', encoding = "UTF-8")
    output <- fromJSON(output)
    quat_N <- str_detect(output$compound_structural_alerts$alert$alert_name,'quaternary nitrogen')
    output_dataframe$quaternary_nitrogens <- NA
    if (T %in% quat_N){
      with_quat_N <- output$compound_structural_alerts$molecule_chembl_id[quat_N]
      output_dataframe$quaternary_nitrogens <- ifelse(output_dataframe$CHEMBL_id %in% with_quat_N,TRUE,FALSE)
    } else {
      output_dataframe$quaternary_nitrogens <- FALSE
    }
    
  } 
  
  colnames(output_dataframe) <- c("CHEMBL_id","SMILES",'INCHIKEY','ALogP','CXLogD','CXLogP','Acid_pka','Base_pka','MolFormula','MW','HBA','HBD','Compound_type','PSA','quat_nitrogens')
  return(output_dataframe)

}

# file_dir <- 'EXp_Barira_output.csv'
# Desinathon_data <- read.csv(file_dir,header = T,fileEncoding="UTF-8")
# inchis <- Desinathon_data$InChIKey
# output2 <- search_compound_by_inchikey(inchis)



