# Install httr package if not already installed
# install.packages("httr")

# Load the httr package
library(httr)
library(jsonlite)

# Function to search compounds by name
PubChem_API_Search <- function(data_table) {
  
  inchikey<-data_table$INCHIKEY
  
  # Define the base URL for PubChem API
  base_url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
  
  #split data into chunks of 10
  chunk_size <- 10
  
  # Calculate the number of chunks
  num_chunks <- ceiling(length(inchikey) / chunk_size)
  
  #Save output in large df
  output_dataframe<- data.frame()
  
  # Process the vector in chunks
  for (i in 1:num_chunks) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, length(inchikey))
    
    chunk <- inchikey[start_idx:end_idx]
  
  # Define the API endpoint for compound name search
  endpoint <- paste0(base_url, "/compound/InChIkey/", URLencode(paste0(chunk,collapse = ',')), "/property/MolecularFormula,MolecularWeight,CanonicalSMILES,InChIKey,XLogP,HBondDonorCount,HBondAcceptorCount,TPSA/JSON")
  
  # Make the API request
  response <- GET(endpoint)
  
  # Check if the request was successful (status code 200)
  if (response$status == 200) {
    # Parse the JSON response
    json_data <- content(response, "text", encoding = "UTF-8")
    result <- fromJSON(json_data)
    result <- result$PropertyTable$Properties
    ## Check if there are any missing columns
    mandatory_columns <- c('CID',"MolecularFormula","MolecularWeight","CanonicalSMILES","InChIKey","TPSA",             
                           "HBondDonorCount","HBondAcceptorCount","XLogP")
    missing_headers <-which(mandatory_columns %!in% colnames(result))
    if (length(missing_headers)>0){
      missing_header_names <- mandatory_columns[missing_headers]
      padding <- matrix(nrow = nrow(result), ncol = length(missing_header_names),NA)
      orig_ncol_results <- ncol(result)
      result <- cbind(result,padding)
      if(length(missing_header_names)==1){
        colnames(result)[ncol(result)] <- missing_header_names
      }else{
        colnames(result)[orig_ncol_results+1:ncol(result)]<-missing_header_names
      }
      
      
    }
    
    output_dataframe <- rbind(output_dataframe,result)
  } else {
    # Print an error message if the request was not successful
    stop("Error in API request: ", http_status(response)$reason)
  }
  }
  
  output_dataframe <- unique(rm_df_cols(output_dataframe,c('CID')))
  
  return(output_dataframe)
}

