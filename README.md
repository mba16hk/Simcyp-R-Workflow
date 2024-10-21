# SimRFlow
A high throughput workflow comprising data collection and simulation of compounds using Certara's Simcyp Simulator. To use the full workflow, please contact Certara for a Simcyp License. The data collection modules of SimRFlow require no license.

# The following is an example use of SimRFlow functions

```bash
#set the working directory to the working directory of the scripts
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

```bash
#import the scripts 
source('input_query.R')
source('SimRFlow_DataCollection_Module.R')
source('httk_search.R')
source('experimental_data_search.R')
source('organise_simulation_data.R')
source('R Workflow.R')
source('PredictParams.R')
source('Additional_data.R')
source('plotting_functions.R')
```
```bash
### Import and process the data before data collection

#set the file directory of the compound file
file_dir <- 'data_files/compound_file.xlsx'

#preprocess the compound file, extract only columns of interest
data <- ProcessInputs(file_dir)
```

```bash
### Collect physicochemical data and build an understanding of missing parameters

#query the chembl, Pubchem and EPI Suite (Norman) database. If you do not wish to query
# pubchem or Norman, please set them to F instead of T
Physicochemical_data <- SimRFlow_DataCollection(data, PubChem = T, Norman = T)

#determine the compounds not found in your database, or which have missing data
Not_found <- MissingInformation(data,Physicochemical_data, missing_info = T)

#determine which compounds can have parameters which are out of the applicablility domain of Simcyp
out_of_range <- OutOfRange_Parameter(Physicochemical_data)

### OPTIONAL STEP: This step is for users who wish to upload additional data
### This may be the case if the user wants to upload their own physicochemical or PK data 
additional_data_directory<- 'data_files/acd_output.xls'
Physicochemical_data<-AdditionalData(data,additional_data_directory,Physicochemical_data, override_existing_data = F)
```

```bash
### This step can be skipped if users do not wish to query httk as part of their workflow
#extract CAS numbers and DTXSID values to query the httk database
CAS_DTXSID <- CAS_and_DTXSID(data)

#search httk library for experimental data using CAS and DTXSIDs
httk_data <- httkSearch(Physicochemical_data, CAS_DTXSID, data,
                        fu_operation = 'arithmetic mean',CLint_operation = 'median')
```

```bash
### OPTIONAL STEP
### For users who have additional experimental data
experimental_data_directory<-'data_files/experimental_data.xlsx'

#incorporate experimental data with data from httk
httk_exp_data <- ExpDataSearch(httk_data, experimental_data_directory, CL_threshold = 3.8)
```

```bash
#organise the data in preparation for running it through Simcyp
organised_data <- OrganiseInputData(httk_data,info = data, admin_route = 'IV Bolus')
```

```bash
## Using Simcyp

#use the prediction module of SimRFlow to predict fu, BP, Vss and Kd
predictions<- PredictParameters(organised_data)

#run the simulations for each compound in bulk
output_profiles<-SimcypSimulation(organised_data, trials = 1, subjects = 5, Time = 24)
```

```bash
## Process Simcyp outputs

#extract additional outputs (which will be used for plotting)
simcyp_outputs <- AdditionalOutputs(organised_data)

#summarise outputs into human-readable table containing simulated compounds
summary_simcyp <-SummaryOutputs(simcyp_outputs)
```

```bash
###Generating plots

#plot concentration-time profiles. Units can be either 'uM' or 'ng/mL'.
plot_profile('CompoundID', output_profiles, curated_data = httk_data,
              units = 'uM', logy=F, tissue_type = 'PLASMA', CI = F)
plot_profile('CompoundID', output_profiles, curated_data = httk_data,
              units = 'ng/ml', logy=T, tissue_type = 'BRAIN', CI = T)
plot_profile('CompoundID', output_profiles, curated_data = httk_data,
              units = 'ng/ml', logy=F, tissue_type = 'ALL', CI = T)

#plot a scatter plot to see the trend between two simulated parameters
plot_parameters(simcyp_outputs, 'CompoundID', 
             plot_type = 'Relationship', 'BW', 'Tmax',
             'blue')

#plot a distribution plot to see variation of a parameter across a populationf for a given compound            
plot_parameters(simcyp_outputs, 'CompoundID', 
                plot_type = 'Distribution', 'Age', 'Tmax',
                'blue')

#create a plot to compare the simulated compounds in high throughput
compare_simulated_compound(summary_simcyp, parameter = 'Vss', bar_order = 'ascending' , bar_col = 'salmon')

#Create a plot to compare the collected physicochemical data against the simcyp predictions
PhyschemvsPredictedParamsPlot(Physicochemical_data,summary_simcyp,physchem_param = 'LogP',predicted_param = 'Ka',plot_colour = 'blue')
```              

