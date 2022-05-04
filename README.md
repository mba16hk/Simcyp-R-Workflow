# SimRFlow
A high throughput workflow comprising data collection and simulation of compounds using Certara's Simcyp Simulator.

# The following is an example use of SimRFlow functions

```bash
#set the working directory to the working directory of the scripts
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

```bash
#import the scripts 
source('input_query.R')
source('chembl_search.R')
source('susdat_search.R')
source('ACD_Labs.R')
source('httk_search.R')
source('experimental_data_search.R')
source('organise_simulation_data.R')
source('R Workflow.R')
source('PredictParams.R')
```
```bash
#set the file directory of the compound file
file_dir <- 'data_files/compound_file.xlsx'
```

```bash
#preprocess the compound file, extract only columns of interest
data <- ProcessInputs(file_dir)
```

```bash
#query the chembl database
chembl_data <- CHEMBLSearch(data)
```

```bash
#determine the compounds not found in chembl
nf_in_chembl <- CompoundsNotFound(data, chembl_data)
```

```bash
#query the Norman suspect database for the compounds not found in Chembl
sus_data <- SusdatSearch(data, nf_in_chembl, chembl_data)
```

```bash
#determine the compounds not found in the Norman suspect list
NOT_FOUND <- NotFoundInsusdat(nf_in_chembl, sus_data)
```

```bash
### This step can be skipped if users do not wish to query httk as part of their workflow
#extract CAS numbers and DTXSID values to query the httk database
CAS_DTXSID <- CAS_and_DTXSID(data)
```

```bash
### This step is for users who wish to upload additional data
ACD_data_directory<- 'data_files/acd_output.xls'
acd_data<-ACD_inputs(data,nf_in_chembl,sus_data,missing_info=T)
physchem_data <- ACD_outputs(data,ACD_data_directory,sus_data)
```

```bash
### For users who do not have additional data
#physchem_data <- sus_data
```

```bash
#search httk library for experimental data using CAS and DTXSIDs
httk_data <- httkSearch(physchem_data, CAS_DTXSID, data)
```

```bash
### For users who have additional experimental data
experimental_data_directory<-'data_files/experimental_data.xlsx'
```

```bash
#incorporate experimental data with data from httk
httk_exp_data <- ExpDataSearch(httk_data, experimental_data_directory, CL_threshold = 3.8)
```

```bash
#organise the data in preparation for running it through Simcyp
organised_data <- OrganiseInputData(httk_data,info = info, admin_route = 'IV Bolus')

```bash
#use the prediction module of SimRFlow to predict fu, BP, Vss and Kd
predictions<- PredictParameters(organised_data)
```

```bash
#run the simulations for each compound in bulk
output_profiles<-SimcypSimulation(organised_data, trials = 1, subjects = 5, Time = 24)
```

```bash
#extract additional outputs (which will be used for plotting)
simcyp_outputs <- AdditionalOutputs(organised_data)
```

```bash
#summarise outputs into human-readable table containing simulated compounds
summary_simcyp <-SummaryOutputs(simcyp_outputs)
```

```bash
#plot concentration-time profiles. Units can be either 'uM' or 'ng/mL'.
plot_profile('CompoundID', output_profiles, curated_data = httk_data,
              units = 'uM', logy=F)
```

```bash
#plot a scatter plot to see the trend between two simulated parameters
plot_parameters(simcyp_outputs, 'CompoundID', 
             plot_type = 'Relationship', 'BW', 'Tmax',
             'blue')
```

```bash
#plot a distribution plot to see variation of a parameter across a populationf for a given compound            
plot_parameters(simcyp_outputs, 'CompoundID', 
                plot_type = 'Distribution', 'Age', 'Tmax',
                'blue')
```               

```bash
#create a chart to compare a given parameter across simulated compounds
compare_simulated_compound(summary_simcyp,'Vss','compound_code','salmon')
```
