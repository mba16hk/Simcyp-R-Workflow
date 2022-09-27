library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(DT)
library(shinyBS)
library(shinyWidgets)


setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#import the scripts 
source('input_query.R')
source('chembl_search.R')
source('susdat_search.R')
source('Additional_data.R')
source('ACD_Labs.R')
source('httk_search.R')
source('experimental_data_search.R')
source('organise_simulation_data.R')
source('Additional_data.R')
source('R Workflow.R')
source('PredictParams.R')
source('plotting_functions.R')
source('SIVA_invitro_distribution.R')


ui <- dashboardPage( skin = 'black',
    dashboardHeader(title = "SimRFlow"),
    dashboardSidebar(
      sidebarMenu(menuItem('Data Collection', 
                           tabName = 'data_collection',
                           icon = icon('file-invoice')),
                  menuItem('Simcyp Modules',
                           tabName = 'simcyp_modules',
                           icon = icon('laptop-code'),
                           menuSubItem("Predictions", tabName = "get_predictions",icon = icon('bolt')),
                           menuSubItem("Human Simulation", tabName = "run_simulation", icon = icon('male'))),
                  menuItem('Plot Outputs',
                           tabName = 'plotting_outputs',
                           icon = icon('chart-area')),
                  # menuItem('In Vitro Distribution',
                  #          tabName = 'SIVA_module',
                  #          icon = icon('flask'), badgeLabel = "new", badgeColor = "green"),
                  menuItem('Help',
                           tabName = 'help_section',
                           icon = icon('question-circle')))
    ), #end of dashboard sidebar
    
    dashboardBody(
      tabItems(
        tabItem('data_collection',
                #first row
                fluidRow(column(width = 3,
                  #Ask the user to input a file of compounds & Experimental Data
                  box(title = 'Physicochemical Data Curation',
                      
                      tipify(fileInput("file1", "CSV/XLSX Compound File",
                                buttonLabel=list(icon("folder"),"Browse"),
                                multiple = FALSE, #doesn't allow multi-file upload
                                accept = c(".csv",".xlsx"),
                                placeholder = 'test_compounds.xlsx'),
                             "Ensure headers are: SMILES, INCHIKEY, CODE, and COMPOUND",
                             placement="bottom", trigger = "hover"),

                      #button to collect physiochemical information
                      uiOutput('ui.action'), width = 12),
                      
                      #br(),
                      
                      #Inputting the Experimental Data
                      conditionalPanel(
                        condition = "input.search_physchem > 0",
                        
                        box(title = 'Experimental Data Curation',
                         
                        p('Select httk databases for experimenmtal data:'),
                               
                        fluidRow(
                          
                        column(width = 4, tipify(checkboxInput("obach_data", "Obach", TRUE),
                               "Obach 2008: Contains fu & systemic clearance data",
                               placement="bottom", trigger = "hover")),
                        
                        column(width = 4, tipify(checkboxInput("wambaugh_data", "Wambaugh", TRUE),
                               "Wambuagh 2019: Contains fu & CLint data",
                               placement="bottom", trigger = "hover")),
                        
                        column(width = 4, tipify(checkboxInput("chemphys_data", "IV DB", TRUE),
                               "In Vitro database: Contains fu, CLint & BP ratio data",
                               placement="bottom", trigger = "hover"))),
                        
                        fluidRow(
                          
                          column(width = 6, selectInput('fu_op', 'Fu Operation',
                                    choices = list('Arithmetic Mean'= 'arithmetic mean',
                                                   'Geometric Mean'='geometric mean',
                                                   'Median'= 'median',
                                                   'Minimum'= 'minimum',
                                                   'Maximum'= 'maximum'))),
                        
                          column(width = 6, selectInput('clint_op', 'CLint Operation',
                                    choices = list('Arithmetic Mean'= 'arithmetic mean',
                                                   'Geometric Mean'='geometric mean',
                                                   'Minimum'= 'minimum',
                                                   'Maximum'= 'maximum')))),
                        
                        conditionalPanel(
                          condition = "input.include_exp_data == 0",
                          #button to collect experimental information from httk
                          tipify(actionButton("search_httk_db",
                                              label = "Search HTTK",
                                              icon = icon('search'),
                                              style = 'color: #fff; 
                           background-color: #a44f2e; border-color: #8a2b07'),
                                 "Search httk databases for any experimental data.",
                                 placement="bottom", trigger = "hover")
                          
                          
                        ),
                        
                        #see if the user wants to provide experimental data
                        tipify(checkboxInput("include_exp_data", "Provide Experimental Data", FALSE),
                               "Checking the box would allow upload of experimental data.",
                               placement="bottom", trigger = "hover"),
                        
                        conditionalPanel(
                          condition = "input.include_exp_data > 0",
                          tipify(fileInput("file2", "CSV/XLSX Experimental Data File",
                                           multiple = FALSE,
                                           buttonLabel=list(icon("folder"),"Browse"),
                                           accept = c(".xlsx",".csv"),
                                           placeholder = 'exp_data.xlsx'),
                                 "Ensure headers are: CODE, BP, FU and CLINT",
                                 placement="bottom", trigger = "hover"),
                          
                          #numeric input for different thresholds
                          conditionalPanel(
                            condition="output.expdummy",
                            p('Optional thresholding of experimental data:'),
                            fluidRow(
                              column(numericInput('CL_thresh', 
                                                  'CLint', 
                                                  0, min = 0, max = NA), 
                                     width = 4
                                     
                              ),
                              
                              column(
                                numericInput('fu_thresh', 
                                             'Fu', 
                                             0, min = 0, max = 1),
                                width = 4
                              ),
                              
                              column(
                                numericInput('BP_thresh', 
                                             'BP Ratio', 
                                             0, min = 0, max = NA),
                                width = 4
                              )
                            ),
                            
                            #see if the user wants to include compounds with missing information
                            tipify(checkboxInput("exp_mean_flag", "Average inputs with httk", FALSE),
                                   "Computed arithmetic mean of your input data and those from httk.",
                                   placement="bottom", trigger = "hover"),
                            
                            #button to collect physiochemical information
                            tipify(actionButton("search_exp",
                                                label = "Search+Organise Experimental Data",
                                                icon = icon('flask'),
                                                style = 'color: #fff; 
                           background-color: #a44f2e; border-color: #8a2b07'),
                                   "Searches HTTK database and organises user-input experimental data.",
                                   placement="bottom", trigger = "hover"),
                          )), #end of conditional panel
                        width = 12,
                        height = "20em")
                        
                        ),#end of box
                  
                  conditionalPanel(
                    #only create this box if there are compounds which have not been found
                    condition="output.ACDLabs_box",
                    box(title = 'Compounds Not Found',
                        
                        #see if the user wants to include compounds with missing information
                        tipify(checkboxInput("include_compounds", "Include compounds with some missing data", TRUE),
                               "Unchecking would return a CSV of compounds not found at all.",
                               placement="bottom", trigger = "hover"),
                        
                        #output a table of missing compounds based on user preferences
                        column(dataTableOutput("TBL6"),
                               height = "300px",
                               style = "height:300px; overflow-y: scroll;overflow-x: scroll;",
                               width = 12),
                        
                        #some information for the user
                        p("Idea: The CSV file can be optionally imported into ACD/Labs
                        for additional physiochemical data. You can also manuallt collect data of interest. 
                          You may proceed without uploading any additional data."),
                        
                        #check if the user would like to upload information they manually collected from ACD labs
                        checkboxInput("upload_ACD_labs", "Upload additional physchem data to supplement/replace existing data", FALSE),
                        
                        #allow user to upload the file form ACD labs if they choose to.
                        conditionalPanel(
                          condition = "input.upload_ACD_labs > 0",
                          fileInput('file3', 
                                    'Upload Additional Physchem Data',
                                    buttonLabel=list(icon("folder"),"Browse"),
                                    multiple = F),
                          actionButton('refresh_physchem',
                                       'Incorporate Additional Data',
                                       icon = icon('table'),
                                       style = 'color: #fff; 
                           background-color: #a44f2e; border-color: #8a2b07')),
                        #height = "20em",
                        status = 'danger',
                        solidHeader = T,
                        width = 12,
                        collapsible = T))

                  ), #end of fluid row
                  
                  #################################################################
                  #             Data Collection Box and Panels                    #
                  #################################################################
                
                #tabBox height  
                tags$head(tags$style("#data_collection_boxes{height:800px !important;}")),
                
                column(width = 9, box( id = 'data_collection_boxes',
                  title = 'Collected Data Tables',
                  status = 'primary',
                  width = 12,
                  
                  ####### different conditionals for different tab boxes ###########
                  
                  
                  
                  ##### Normal Tab Box #############
                  
                  conditionalPanel(condition = "input.search_physchem == 0 && (input.search_httk_db == 0 && input.search_exp == 0) && input.upload_ACD_labs == 0",
                                   tabBox(
                                     
                                     # Title can include an icon
                                     title = tagList(shiny::icon("database")),
                                     height = "730px",
                                     
                                     #Visualise the Input Data
                                     tabPanel("Input Data", height = "650px", width = 12,
                                              column(dataTableOutput("TBL1"),
                                                     height = "650px",
                                                     style = "height:650px; overflow-y: scroll;overflow-x: scroll;",
                                                     width = 12)),
                                     width = 12
                                   ) #end of tabbox
                                   ),

                  ##### Tab Box with out of range PSA and HBD values and physchem data, condition is physchem search button click #############
                  
                  
                  conditionalPanel(
                    #only create this tab if additional data is uploaded
                    condition="input.search_physchem == 1 && (input.search_httk_db == 0 && input.search_exp == 0)&& input.upload_ACD_labs == 0",
                  tabBox(
                    
                    # Title can include an icon
                    title = tagList(shiny::icon("database")),
                    height = "730px",
                    
                    #Visualise the Input Data
                    tabPanel("Input Data", height = "650px", width = 12,
                             column(dataTableOutput("input_tbl_2"),
                                    height = "650px",
                                    style = "height:650px; overflow-y: scroll;overflow-x: scroll;",
                                    width = 12)),
                    
                    #Visualise the physiochemical Data
                    tabPanel("Physchem Data",  height = "650px", width = 12,
                             column(dataTableOutput("TBL2"), 
                                    height = "650px",
                                    style = "height:650px; overflow-y: scroll;overflow-x: scroll;",
                                    width = 12)),
                    
                    #Visualise the out or range Peff values
                    tabPanel("Out of Range Peff",  height = "640px", width = 12,
                            p("The Peff prediction method relies on PSA/HBS values. The models were calibrated for compounds with 60<MW<455, -3<logP<4, HBD<6 and 16.2<PSA<154.4. The following compounds have out of range parameters and may result in inaccurate predictions."),
                          column(dataTableOutput("TBL_out_of_range"),
                                height = "640px",
                                style = "height:640px; overflow-y: scroll",
                                width = 12)),
                    
                    width = 12
                  )), #end of tabbox
                  
                  ##### Tab Box with experimental data, condition is exp data search button click #############
                  
                  conditionalPanel(
                    #only create this tab if additional data is uploaded
                    condition="input.search_physchem == 1 && (input.search_httk_db == 1 || input.search_exp == 1) && input.upload_ACD_labs == 0",
                    tabBox(
                      
                      # Title can include an icon
                      title = tagList(shiny::icon("database")),
                      height = "730px",
                      
                      #Visualise the Input Data
                      tabPanel("Input Data", height = "650px", width = 12,
                               column(dataTableOutput("input_tbl_3"),
                                      height = "650px",
                                      style = "height:650px; overflow-y: scroll;overflow-x: scroll;",
                                      width = 12)),
                      
                      #Visualise the physiochemical Data
                      tabPanel("Physchem Data",  height = "650px", width = 12,
                               column(dataTableOutput("physchem_tbl_2"), 
                                      height = "650px",
                                      style = "height:650px; overflow-y: scroll;overflow-x: scroll;",
                                      width = 12)),
                      
                      #Visualise the out or range Peff values
                      tabPanel("Out of Range Peff",  height = "640px", width = 12,
                               p("The Peff prediction method relies on PSA/HBS values. The models were calibrated for compounds with 60<MW<455, -3<logP<4, HBD<6 and 16.2<PSA<154.4. The following compounds have out of range parameters and may result in inaccurate predictions."),
                               column(dataTableOutput("TBL_out_of_range2"),
                                      height = "640px",
                                      style = "height:640px; overflow-y: scroll",
                                      width = 12)),
                      
                      #Visualise the experimental Data
                      tabPanel("Experimental Data",  height = "650px", width = 12,
                               column(dataTableOutput("TBL3"),
                                      height = "650px",
                                      style = "height:650px; overflow-y: scroll;overflow-x: scroll;",
                                      width = 12)),
                      
                      width = 12
                    )), #end of tabbox
                  
                  ##### Tab Box with experimental data, condition is exp data search button click and additional physchem data #############
                  
                  conditionalPanel(
                    #only create this tab if additional data is uploaded
                    condition="input.search_physchem == 1 && (input.search_httk_db == 1 || input.search_exp == 1) && input.upload_ACD_labs > 0",
                    tabBox(
                      
                      # Title can include an icon
                      title = tagList(shiny::icon("database")),
                      height = "730px",
                      
                      #Visualise the Input Data
                      tabPanel("Input Data", height = "650px", width = 12,
                               column(dataTableOutput("input_tbl_4"),
                                      height = "650px",
                                      style = "height:650px; overflow-y: scroll;overflow-x: scroll;",
                                      width = 12)),
                      
                      #Visualise the physiochemical Data
                      tabPanel("Physchem Data",  height = "650px", width = 12,
                               column(dataTableOutput("physchem_tbl_3"), 
                                      height = "650px",
                                      style = "height:650px; overflow-y: scroll;overflow-x: scroll;",
                                      width = 12)),
                      
                      #Visualise the out or range Peff values
                      tabPanel("Out of Range Peff",  height = "640px", width = 12,
                               p("The Peff prediction method relies on PSA/HBS values. The models were calibrated for compounds with 60<MW<455, -3<logP<4, HBD<6 and 16.2<PSA<154.4. The following compounds have out of range parameters and may result in inaccurate predictions."),
                               column(dataTableOutput("TBL_out_of_range3"),
                                      height = "640px",
                                      style = "height:640px; overflow-y: scroll",
                                      width = 12)),
                      
                      #Visualise the experimental Data
                      tabPanel("Experimental Data",  height = "650px", width = 12,
                               column(dataTableOutput("exp_dat_tbl_2"),
                                      height = "650px",
                                      style = "height:650px; overflow-y: scroll;overflow-x: scroll;",
                                      width = 12)),
                      
                      #only create this tab if additional data is uploaded
                      tabPanel("Physchem Data + Further Data",  height = "650px", width = 12,
                              column(dataTableOutput("TBL7"),
                                     height = "650px",
                                     style = "height:650px; overflow-y: scroll;overflow-x: scroll;",
                                     width = 12)),
                      
                      width = 12
                    )), #end of tabbox
                  
                  collapsible = TRUE))) #end of box
                
              
                ), #end of tabItem 1
        
        tabItem('get_predictions',
                
                #Set some of simcyp's simulation parameters
                tipify(actionButton("predict_params_button",
                                    label = "Predict Parameters",
                                    icon = icon('laptop-code'),
                                    style = 'color: #fff;
                                     background-color: #a44f2e; 
                                     border-color: #8a2b07'),#end of action button
                       "Predict Fu, Vss, BP and Kd using Simcyp",
                       placement="top", trigger = "hover"),
                
                box(
                  title = 'Predicted Outputs',
                  status = 'primary',
                  width = 12,
                  tabBox(
                    
                    #visualise simcyp outputs table
                    tabPanel('Predicted Outputs',
                             column(dataTableOutput("pred_output_tbl"),
                                    width = 12)), width = 12), collapsible = T)
                
                ),
        
        tabItem('run_simulation',
                fluidRow(
                    #Parameters for simulation
                    box(title = 'Simulation Parameters',
                        status = 'success',
                        height = 600,
                          
                        #Select single or multiple dosing
                        radioButtons('dosing_options','Dosing Options', 
                                     choices = c('Single Dosing',
                                                 'Multiple Dosing'),
                                     inline = T,
                                     selected = F),
                        
                            #conditional panel for dose value
                            conditionalPanel(
                            condition = "output.determine_inputs == 0 || output.determine_inputs == 2 ||output.determine_inputs == 3 || output.determine_inputs == 5",
                            numericInput('dose_value',
                                         'Dose',
                                         100, min = 0, max = NA)),
                        
                        
                            # conditional panel for dose units
                            conditionalPanel(
                              condition = "output.determine_inputs == 0 || output.determine_inputs == 1 ||output.determine_inputs == 3 || output.determine_inputs == 6",
                              selectInput('dose_units',
                                          'Dose Units',
                                          list("mg/kg" = "mg/kg",
                                               "mg" = "mg",
                                               "mg/m^2" = "mg/m^2"))),
                        
                        #,
                        p('Assumptions:'),
                        fluidRow(
                          
                          column(
                            width = 4,
                            tipify(numericInput('acid_bp_ratio', 'Acid BP ratio', 0.55,
                                         min = 0.55, max = NA),
                                   "Acids without any BP ratio value are assumed to have this BP ratio",
                                   placement="bottom", trigger = "hover")
                          ),
                          
                          
                          column(
                            width = 4,
                            tipify(numericInput('agp_threshold', 'pKa threshold for AGP binding to bases', 7,
                                         min = 5.5, max = 8.5),
                                   "Bases that have a pKa value higher than this threshold are assumed to bing to AGP",
                                   placement="bottom", trigger = "hover")
                          ),
                          
                          
                          column(
                            width = 4, conditionalPanel(
                            condition = "output.determine_inputs == 0 || output.determine_inputs == 1 ||output.determine_inputs == 2 || output.determine_inputs == 4",
                            numericInput('pred_method',
                                         'Vss Prediction Method',
                                         3, min = 1, max = 3, step = 1))),
                          
                        ),
                        
                        p('Trial Design:'),
                        fluidRow(
                          
                          column(
                            width = 4,
                            numericInput('subjects', ' Subjects', 10,
                                         min = 1, max = NA, step = 1),
                            tipify(checkboxInput("set_additional_subject_design", "Specify additional subject parameters", FALSE),
                                   "Checking would allow setting more subject design parameters.",
                                   placement="bottom", trigger = "hover")
                          ),
                    
                          
                          column(
                            width = 4,
                            numericInput('sim_time', 'Duration', 24,
                                         min = 4, max = 480, step = 0.5)
                          ),
                          
                          column(
                            width = 4,
                            selectInput('administration_route',
                                        'Administration Route',
                                        list("Oral" = "Oral",
                                             "IV Bolus" = "IV Bolus",
                                             "Dermal" = "Dermal"))
                           
                            )
                        ),
                        
                        fluidRow(
                          
                          column(
                            width = 4,
                            #see if the user wants to include compounds with missing information
                            tipify(checkboxInput("set_seed", "Set seed", TRUE),
                                   "Unchecking would not set a seed for the simulations.",
                                   placement="bottom", trigger = "hover"),
                          ),
                          
                          conditionalPanel(
                            condition = "output.simulate_check",
                            #Set some of simcyp's simulation parameters
                            tipify(actionButton("simulate_button",
                                                label = "Simulate",
                                                icon = icon('laptop-code'),
                                                style = 'color: #fff;
                                     background-color: #a44f2e; 
                                     border-color: #8a2b07'),#end of action button
                                   "Simulate all compounds in the Simcyp Simulator.",
                                   placement="top", trigger = "hover")
                          )
                          
                        ),

                        width = 8), #end of box
                    
                    conditionalPanel(
                      condition = "input.set_additional_subject_design == 1",
                      box(
                        title = 'Subject Design',
                        status = 'primary',
                        width = 4,
                        fluidRow(
                          
                          #inputs for minimum age
                          column(
                            width = 12,
                            sliderInput('age_range', 'Age Range', 
                                        min = 18, max = 65, value = c(20,50), step = 1)
                          ), #end of column
                          
                          #proportion of females
                          column(
                            width = 12,
                            numericInput('female_prop', 'Female Proportion',
                                         value = 0.5, min = 0, max = 1)
                          ) #end of column
                        ) #end of fluidrow
                      ) #end of box
                    ),#end of conditional panel
                    
                    conditionalPanel(
                      condition = "input.administration_route == 'Dermal'",
                      box(
                        title = 'Dermal Parameters',
                        status = 'primary',
                        width = 4,
                        fluidRow(
                          
                          #inputs for dermal area
                          column(
                            width = 12,
                            numericInput('dermal_area_val', 'Area of Application (cm^2)', 
                                        min = 1, max = 500, value = 100, step = 0.5)
                          ), #end of column
                          
                          # inputs for formulation thickness
                          column(
                            width = 12,
                            numericInput('dermal_formulation_thickness', 'Formulation Thickness (cm)',
                                         value = 0.005, min = 0.005, max = 0.1)
                          ), #end of column
                          
                          # inputs for formulation density
                          column(
                            width = 12,
                            numericInput('dermal_formulation_density', 'Formulation Density (g/mL)',
                                         value = 1, min = 0.05, max = 4)
                          ) #end of column
                          
                        ) #end of fluidrow
                      ) #end of box
                    ),#end of conditional panel
          
                  ),# end of fluidrow

                
                conditionalPanel(
                  condition = "input.dosing_options == 'Multiple Dosing'",
                  box(
                    title = 'Multiple Dosing Options',
                    status = 'primary',
                    width = 8,
                    fluidRow(
                      
                      #inputs for number of doses
                      column(
                        width = 6,
                        numericInput('num_doses_val', 'Number of doses', 
                                     min = 1, max = 1000, value = 2, step = 1)
                      ), #end of column
                      
                      # inputs for dosing interval
                      column(
                        width = 6,
                        numericInput('dosing_interval_val', 'Dosing Interval (hrs)',
                                     value = 12, min = 0.001, max = 48)
                      ) #end of column
                      
                    ) #end of fluidrow
                  ) #end of box
                ),#end of conditional panel
                
                box(
                  title = 'Simcyp Output Tables',
                  status = 'primary',
                  width = 12,
                  tabBox(
                    
                    #visualise simcyp outputs table
                    tabPanel('Concentration-Time Profiles',
                             'The concentration time profiles of all simulated compounds and their corresponding subjects in different tissues:',
                             column(dataTableOutput("TBL4"),
                                    width = 12)),
                    
                    
                    #visualise simcyp outputs table
                    tabPanel('Summary of Simulated Outputs',
                             column(dataTableOutput("TBL5"),
                                    width = 12)),
                    
                  width = 12), collapsible = TRUE) #end of box
        ), #end of tabitem 2
        
        tabItem('plotting_outputs',
                
                #Custom CSS
                tags$head(tags$style("#conc_time_box{height:1000px !important;}")),
                
                #plotting
                conditionalPanel(
                  condition = "input.simulate_button > 0",
                box( id = "conc_time_box",
                  
                  title = 'Plot Simulated Outputs',
                  
                  tabBox( height = "900px",
                          
                          ##########################################################################
                          ######################## conc-time-profile ###############################
                          ##########################################################################
                          
                          tabPanel("Concentration-Time Profiles", height = "760px",
                                   column(width = 12,
                                          
                                          #Allow user to enter compound code and tissue-type
                                          fluidRow(column(uiOutput("compound_lists"), width = 6),column(uiOutput("tissue_lists"), width = 6)),
                                          
                                          #select y axis units
                                          selectInput('unit', 'Concentration Units',
                                                      choices = list('ng/mL'= 'ng/mL','uM'='uM')),
                                          
                                          #Select Log Scale or Normal Scale
                                          radioButtons('log_scale','y-axis scaling', 
                                                       choices = c('Logarithmic'= T,
                                                                   'Natural'= F),
                                                       selected = F), 

                                          #plot the concentration-time profile
                                          plotOutput("conc_time_plot", height = '560px')), 
                                          #download the conc-time plot
                                          downloadButton(outputId = "DownloadConcPlot", label = "Download Plot"),
                                   
                                   width = 12),
                          
                          ##########################################################################
                          ##                      Distribution+rln Plot                        #####
                          ##########################################################################
                          
                          tabPanel("Parameter Distributions and Relationships", height = "760px",
                                   column(width = 6, 
                                          
                                          #Enter the compound code  
                                          uiOutput("compound_lists2"),
                                          
                                          #select x axis
                                          selectInput('x_axis_var1', 'x Variable',
                                                      choices = list('Plasma Tmax'= 'Tmax',
                                                                     'Plasma Cmax'='Cmax_PLASMA',
                                                                     'AUC'= 'AUC',
                                                                     'AUCinf'= 'AUCinf',
                                                                     'Half Life' = 'HalfLife',
                                                                     'Accumulation Index'= 'AccumulationIndex',
                                                                     'BSA'='BSA',
                                                                     'Age'='Age',
                                                                     'BW'='BW',
                                                                     'GFR'='GFR',
                                                                     # 'Volume of Distribution'='Vss',
                                                                     # 'Fraction Absorbed'='Fa',
                                                                     # 'Fraction unbound in plasma'= 'Fu_plasma',
                                                                     'Ka'='Ka',
                                                                     'Skin Cmax' = 'Cmax_SKIN',
                                                                     'Kidney Cmax' =  'Cmax_KIDNEY',
                                                                     'Brain Cmax' = 'Cmax_BRAIN',
                                                                     'Heart Cmax' = 'Cmax_HEART',
                                                                     'Gut Cmax' = 'Cmax_GUT',
                                                                     'Lung Cmax' = 'Cmax_LUNG',
                                                                     'Liver Cmax' = 'Cmax_LIVER',
                                                                     'Pancreas Cmax' = 'Cmax_PANCREAS',
                                                                     'Spleen Cmax' = 'Cmax_SPLEEN',
                                                                     'Muscle Cmax' = 'Cmax_MUSCLE',
                                                                     'Adipose Cmax' = 'Cmax_ADIPOSE')),
                                          
                                            #select plot colour
                                            selectInput('plot_col', 'Plot Colour',
                                                        choices = list('Green'= 'darkolivegreen3',
                                                                       'Burgundy' = 'indianred',
                                                                       'Red' = 'firebrick1',
                                                                       'Blue'='royalblue',
                                                                       'Yellow'='goldenrod1',
                                                                       'Pink'='hotpink3',
                                                                       'Purple'='mediumpurple2',
                                                                       'Orange'='orange',
                                                                       'Teal'='lightseagreen',
                                                                       'Black' = 'gray0')),
                                          
                                          #Distribution Plot
                                          plotOutput("distribution_plot",  height = '560px'),
                                          #download the distribution_plot plot
                                          downloadButton(outputId = "DownloaddistPlot", label = "Download Plot")
                                         ), #end of column
                                   
                                   column(width = 6,
                                          
                                          #Enter the compound code  
                                          uiOutput("compound_lists3"),
                                          
                                          fluidRow(column(width = 6,
                                                          #select x axis
                                                          selectInput('x_axis_var2', 'x Variable',
                                                                      choices = list('Plasma Tmax'= 'Tmax',
                                                                                     'Plasma Cmax'='Cmax_PLASMA',
                                                                                     'AUC'= 'AUC',
                                                                                     'AUCinf'= 'AUCinf',
                                                                                     'Half Life' = 'HalfLife',
                                                                                     'Accumulation Index'= 'AccumulationIndex',
                                                                                     'BSA'='BSA',
                                                                                     'Age'='Age',
                                                                                     'BW'='BW',
                                                                                     'GFR'='GFR',
                                                                                     # 'Volume of Distribution'='Vss',
                                                                                     # 'Fraction Absorbed'='Fa',
                                                                                     # 'Fraction unbound in plasma'= 'Fu_plasma',
                                                                                     'Ka'='Ka',
                                                                                     'Skin Cmax' = 'Cmax_SKIN',
                                                                                     'Kidney Cmax' =  'Cmax_KIDNEY',
                                                                                     'Brain Cmax' = 'Cmax_BRAIN',
                                                                                     'Heart Cmax' = 'Cmax_HEART',
                                                                                     'Gut Cmax' = 'Cmax_GUT',
                                                                                     'Lung Cmax' = 'Cmax_LUNG',
                                                                                     'Liver Cmax' = 'Cmax_LIVER',
                                                                                     'Pancreas Cmax' = 'Cmax_PANCREAS',
                                                                                     'Spleen Cmax' = 'Cmax_SPLEEN',
                                                                                     'Muscle Cmax' = 'Cmax_MUSCLE',
                                                                                     'Adipose Cmax' = 'Cmax_ADIPOSE'))),
                                                          
                                                          column(width = 6,
                                                                 #select y axis
                                                                 selectInput('y_axis_var', 'y Variable',
                                                                             choices = list('Plasma Tmax'= 'Tmax',
                                                                                            'Plasma Cmax'='Cmax_PLASMA',
                                                                                            'AUC'= 'AUC',
                                                                                            'AUCinf'= 'AUCinf',
                                                                                            'Half Life' = 'HalfLife',
                                                                                            'Accumulation Index'= 'AccumulationIndex',
                                                                                            'BSA'='BSA',
                                                                                            'Age'='Age',
                                                                                            'BW'='BW',
                                                                                            'GFR'='GFR',
                                                                                            # 'Volume of Distribution'='Vss',
                                                                                            # 'Fraction Absorbed'='Fa',
                                                                                            # 'Fraction unbound in plasma'= 'Fu_plasma',
                                                                                            'Ka'='Ka',
                                                                                            'Skin Cmax' = 'Cmax_SKIN',
                                                                                            'Kidney Cmax' =  'Cmax_KIDNEY',
                                                                                            'Brain Cmax' = 'Cmax_BRAIN',
                                                                                            'Heart Cmax' = 'Cmax_HEART',
                                                                                            'Gut Cmax' = 'Cmax_GUT',
                                                                                            'Lung Cmax' = 'Cmax_LUNG',
                                                                                            'Liver Cmax' = 'Cmax_LIVER',
                                                                                            'Pancreas Cmax' = 'Cmax_PANCREAS',
                                                                                            'Spleen Cmax' = 'Cmax_SPLEEN',
                                                                                            'Muscle Cmax' = 'Cmax_MUSCLE',
                                                                                            'Adipose Cmax' = 'Cmax_ADIPOSE')))
                                                          ),
                                          
   
                                          #select plot colour
                                          selectInput('plot_col2', 'Plot Colour',
                                                      choices = list('Green'= 'darkolivegreen3',
                                                                     'Burgundy' = 'indianred',
                                                                     'Red' = 'firebrick1',
                                                                     'Blue'='royalblue',
                                                                     'Yellow'='goldenrod1',
                                                                     'Pink'='hotpink3',
                                                                     'Purple'='mediumpurple2',
                                                                     'Orange'='orange',
                                                                     'Teal'='lightseagreen',
                                                                     'Black' = 'gray0')),
                                          
                                          plotOutput("relationship_plot",  height = '560px'),
                                          #download the relationship_plot plot
                                          downloadButton(outputId = "DownloadrlnPlot", label = "Download Plot")),#end of column

                                   width = 12), #end of tab Panel
                          
                          tabPanel("Cross-Compound Comparison Charts", height = "760px",
                                   column(width = 12,
                                          
                                          #Allow user to enter compound code and tissue-type
                                          fluidRow(column(# Parameters
                                            selectInput('sim_parameters', 'Simulated Parameters',
                                                        choices = list('Plasma Tmax'= 'Tmax',
                                                                       'Plasma Cmax'='Cmax_PLASMA',
                                                                       'AUC'= 'AUC',
                                                                       'AUCinf'= 'AUCinf',
                                                                       'Half Life' = 'HalfLife',
                                                                       'Accumulation Index'= 'AccumulationIndex',
                                                                       'Volume of Distribution'='Predicted Vss',
                                                                       'Fraction Absorbed'='Predicted Fa',
                                                                       'Fraction unbound in plasma'= 'Predicted Fu',
                                                                       'Ka'='Ka',
                                                                       'Skin Cmax' = 'Cmax_SKIN',
                                                                       'Kidney Cmax' =  'Cmax_KIDNEY',
                                                                       'Brain Cmax' = 'Cmax_BRAIN',
                                                                       'Heart Cmax' = 'Cmax_HEART',
                                                                       'Gut Cmax' = 'Cmax_GUT',
                                                                       'Lung Cmax' = 'Cmax_LUNG',
                                                                       'Liver Cmax' = 'Cmax_LIVER',
                                                                       'Pancreas Cmax' = 'Cmax_PANCREAS',
                                                                       'Spleen Cmax' = 'Cmax_SPLEEN',
                                                                       'Muscle Cmax' = 'Cmax_MUSCLE',
                                                                       'Adipose Cmax' = 'Cmax_ADIPOSE')), width = 6),
                                            
                                            column(selectInput('plot_col_2', 'Bar Colour',
                                                               choices = list('Green'= 'darkolivegreen3',
                                                                              'Burgundy' = 'indianred',
                                                                              'Red' = 'firebrick1',
                                                                              'Blue'='royalblue',
                                                                              'Yellow'='goldenrod1',
                                                                              'Pink'='hotpink3',
                                                                              'Purple'='mediumpurple2',
                                                                              'Orange'='orange',
                                                                              'Teal'='lightseagreen',
                                                                              'Black' = 'gray0')), width = 6)),

                                          #Select in ascending order or compound alphanumeric order
                                          radioButtons('param_order','x-axis ordering', 
                                                       choices = c('Ascending Value'= 'ascending',
                                                                   'Compound Code'= 'compound_code'),
                                                       selected = 'compound_code'),
                                          
                                          #plot the concentration-time profile
                                          plotOutput("comp_comparison_plot", height = '600px')), 
                                   #download the conc-time plot
                                   downloadButton(outputId = "DownloadCompCmparPlot", label = "Download Plot"),
                                   
                                   width = 12),
                          
                          tabPanel("Physchem Properties vs. Predicted Parameters", height = "760px",
                                   column(width = 12,
                                          
                                          #Allow user to enter compound code and tissue-type
                                          fluidRow(column(# Parameters
                                            selectInput('pred_parameters', 'Predicted Parameters',
                                                        choices = list('Plasma Tmax'= 'Tmax',
                                                                       'Plasma Cmax'='Cmax_PLASMA',
                                                                       'AUC'= 'AUC',
                                                                       'AUCinf'= 'AUCinf',
                                                                       'Half Life' = 'HalfLife',
                                                                       'Accumulation Index'= 'AccumulationIndex',
                                                                       'Volume of Distribution'='Predicted Vss',
                                                                       'Fraction Absorbed'='Predicted Fa',
                                                                       'Fraction unbound in plasma'= 'Predicted Fu',
                                                                       'Ka'='Ka',
                                                                       'Skin Cmax' = 'Cmax_SKIN',
                                                                       'Kidney Cmax' =  'Cmax_KIDNEY',
                                                                       'Brain Cmax' = 'Cmax_BRAIN',
                                                                       'Heart Cmax' = 'Cmax_HEART',
                                                                       'Gut Cmax' = 'Cmax_GUT',
                                                                       'Lung Cmax' = 'Cmax_LUNG',
                                                                       'Liver Cmax' = 'Cmax_LIVER',
                                                                       'Pancreas Cmax' = 'Cmax_PANCREAS',
                                                                       'Spleen Cmax' = 'Cmax_SPLEEN',
                                                                       'Muscle Cmax' = 'Cmax_MUSCLE',
                                                                       'Adipose Cmax' = 'Cmax_ADIPOSE')), width = 4),
                                            
                                            column(# Parameters
                                              selectInput('physchem_parameters', 'Physchem Parameters',
                                                          choices = list('Compound Characterisation' = 'Compound_type',
                                                                         'Molecular Weight'= 'MW',
                                                                         'LogP'='logP',
                                                                         'Polar Surface Area'= 'PSA',
                                                                         'Hydrogen Bond Donor Count'= 'HBD')), width = 4),
                                            
                                            column(selectInput('plot_col3', 'Bar Colour',
                                                               choices = list('Green'= 'darkolivegreen3',
                                                                              'Burgundy' = 'indianred',
                                                                              'Red' = 'firebrick1',
                                                                              'Blue'='royalblue',
                                                                              'Yellow'='goldenrod1',
                                                                              'Pink'='hotpink3',
                                                                              'Purple'='mediumpurple2',
                                                                              'Orange'='orange',
                                                                              'Teal'='lightseagreen',
                                                                              'Black' = 'gray0')), width = 4)),
                                          

                                          #plot the concentration-time profile
                                          plotOutput("physchem_prediction_comparison_plot", height = '600px')), 
                                   #download the conc-time plot
                                   downloadButton(outputId = "DownloadphyschemCompCmparPlot", label = "Download Plot"),
                                   
                                   width = 12),
                          
                          width = 12
                  ),
                  
                  status = "primary", 
                  solidHeader = FALSE,
                  width = 12)), #end of box
                
        ), #end of tabItem
        
        tabItem('SIVA_module',
                
                fluidRow(
                  
                  tags$head(tags$style("#SIVA_boxes{height:360px !important;}")),
                  
                       #Ask the user to input a file of compounds & Experimental Data
                       box(id = "SIVA_boxes",title = 'In Vitro Distribution Inputs',
                           
                           fluidRow(column(width = 9, 
                                           tipify(fileInput("SIVA_file", "CSV/TXT Input File",
                                                            buttonLabel=list(icon("folder"),"Browse"),
                                                            multiple = FALSE, #doesn't allow multi-file upload
                                                            accept = c(".txt",".csv"),
                                                            placeholder = 'cmpnds.csv'),
                                                  "Ensure headers 'CAS, INCHIKEY, NAME, FU' are present",
                                                  placement="bottom", trigger = "hover")),
                                    
                                    column(width = 3, style = "margin-top: +25px;",
                                           #button to collect information for SIVA
                                           tipify(actionButton("search_SIVA_data",
                                                               label = "Curate Data",
                                                               icon = icon('search'),
                                                               style = 'color: #fff; background-color: #a44f2e; border-color: #8a2b07'),
                                                  "Search for physchem data and Henry's law constant.",
                                                  placement="bottom", trigger = "hover"))
                                    
                                    ),

                           conditionalPanel(
                             condition = "input.search_SIVA_data > 0",
                             fluidRow(
                               
                               column(width = 4,
                                      boxPad(
                                        color = "red",
                                        style = "height: 75px",
                                        descriptionBlock(
                                          header = textOutput("missing_comp_num"), 
                                          text = "Missing Compounds", 
                                          rightBorder = FALSE,
                                          marginBottom = TRUE
                                        ))),
                               
                               column(
                                 width = 4,
                                 boxPad(
                                   color = "yellow",
                                   style = "height: 75px",
                                   descriptionBlock(
                                     header = textOutput("duplicated_comp_num"),  
                                     text = "Duplicated Compounds", 
                                     rightBorder = FALSE,
                                     marginBottom = TRUE
                                   ))),
                               
                               column(width = 4,
                                      boxPad(
                                        color = "blue",
                                        style = "height: 75px",
                                        descriptionBlock(
                                          header = textOutput("missing_hlc_num"),  
                                          text = "Missing HLC", 
                                          rightBorder = FALSE,
                                          marginBottom = TRUE
                                        )))
                             ),
                             
                             br(),
                             
                             conditionalPanel(condition = "output.missing_comp_num==0 && output.duplicated_comp_num==0 && output.missing_hlc_num ==0",
                                              p("You may proceed to predicting the in vitro distribution.")),
                             
                             conditionalPanel(condition = "output.missing_comp_num!=0 || output.duplicated_comp_num!=0 || output.missing_hlc_num !=0",
                                              p("Optional: You may download the curated data file, populate the missing entries by manullay searching for them, and reupload here:"),
                                              
                                              fluidRow(column(width = 9,tipify(fileInput("SIVA_file_reupload", "CSV/EXCEL Input File",
                                                               buttonLabel=list(icon("folder"),"Browse"),
                                                               multiple = FALSE, #doesn't allow multi-file upload
                                                               accept = c(".xlsx",".csv"),
                                                               placeholder = 'cmpnds_updated.csv'),
                                                     "Optionally re-upload the collected data file containing missing entries.",
                                                     placement="bottom", trigger = "hover")),
                                              
                                              column(width = 3, style = "margin-top: +25px;",
                                                     #button to collect information for SIVA
                                                     tipify(actionButton("update_SIVA_data",
                                                                         label = "Update Data",
                                                                         icon = icon('search'),
                                                                         style = 'color: #fff; 
                           background-color: #a44f2e; border-color: #8a2b07'),
                                                            "Update the SIVA curated data.",
                                                            placement="bottom", trigger = "hover")))
                                              
                                              
                                              ),
                             
                             
                             
                           ),

                           status = "primary", 
                           solidHeader = FALSE,
                           width = 6),
                       
                       box(id = "SIVA_boxes", title = 'In Vitro Distribution Parameters',
                           
                           
                           fluidRow(
                             column(
                               width = 6,
                               
                               #set the nominal concentration
                               numericInput('nominal_conc', 'Nominal Concentration (M)',
                                            value = 1e-5, min = 1e-10, max = 1e10),
                               
                               #set the number of cells
                               numericInput('number_cells', 'Cell Number',
                                            value = 1000, min = 100, max = 20000),
                               
                               #set the pH of the medium
                               numericInput('medium_pH', 'Media pH',
                                            value = 7.4, min = 0, max = 14)
                               
                               
                             ), #end of column
                             
                             column(
                               width = 6,
                               
                               #set the volume of the well
                               numericInput('vol_well', 'Well Volume (uL)',
                                            value = 360, min = 50, max = 10000),
                               
                               #set the well diameter
                               numericInput('well_diameter', 'Well Diameter (mm)',
                                            value = 6.4, min = 0.1, max = 25),
                               
                               #set the serum fraction
                               numericInput('serum_fraction', 'Serum fraction',
                                            value = 0.1, min = 0.00000001, max = 1)
                               
                               
                             )), #end of column
                           
                           #button to collect information for SIVA
                           actionButton("predict_iv_distribution",
                                               label = "Predict In Vitro Distribution",
                                               icon = icon('laptop-code'),
                                               style = 'color: #fff; 
                           background-color: #a44f2e; border-color: #8a2b07; padding:8px; font-size:100%; float:center; width:330px'),

                           status = "primary", 
                           solidHeader = FALSE,
                           width = 6)
                       
                       ),
                
                #Custom CSS
                tags$head(tags$style("#my_box{height:700px !important;}")),
                
                conditionalPanel(
                  condition = "input.search_SIVA_data > 0",
                
                box(id = "my_box",title = 'In Vitro Distribution Outputs',
                    
                    tabBox( height = "630px",
                      #Curated data table
                      tabPanel("Curated Data", height = "560px",
                               column(width = 12, dataTableOutput("SIVA_curated_data_tbl"),style = "height:560px; overflow-y: scroll;overflow-x: scroll;"), 
                                      width = 12),
                      
                      #SIVA output table
                      tabPanel("In Vitro Distribution Predictions", height = "560px",
                               column(width = 12, dataTableOutput("SIVA_output_tbl"),style = "height:560px; overflow-y: scroll;overflow-x: scroll;"),  
                                      width = 12),
                      
                      #SIVA output plots
                      tabPanel("In Vitro Distribution Plots", height = "560px",
                               column(width = 12, plotOutput("SIVA_output_plots",height = '500px')), 
                               downloadButton(outputId = "Cat10", label = "Download Plots"),
                               width = 12),
                      
                      width = 12
                    ),
                    
                    status = "primary", 
                    solidHeader = FALSE,
                    width = 12))
                
                
                
                
                
        ), #end of tabItem
        
        tabItem('help_section',
                uiOutput('HelpPageText'),
                #HTML("<ul><li>...text...</li><li>...more text...</li></ul>") #bullet points
                #h1(HTML(paste0("Hello O",tags$sub("2")))) #subscript
                #img(src = 'happy.jpg', height = 72, width = 72)
                )
        ) #end of tabitems
                
      ) #end of dashboardbody
) #end of dashboardpage

server <- function(input, output, session) {
  
    data <- reactive({
        file1 <- input$file1
        req(file1)
        ProcessInputs(file1$datapath)

    })
    
    #Display the collected information
    output$TBL1 <- renderDataTable(data(),
                                   rownames = FALSE,
                                   extensions = 'Buttons',
                                   options = list(
                                     paging = FALSE,
                                     searching = TRUE,
                                     ordering = TRUE,
                                     #scrollCollapse= TRUE,
                                     dom = 'l<"sep">Bfrtip',
                                     buttons = c('csv', 'excel')
                                   ))
    
    output$input_tbl_2 <- renderDataTable(data(),
                                          rownames = FALSE,
                                          extensions = 'Buttons',
                                          options = list(
                                            paging = FALSE,
                                            searching = TRUE,
                                            ordering = TRUE,
                                            #scrollCollapse= TRUE,
                                            dom = 'l<"sep">Bfrtip',
                                            buttons = c('csv', 'excel')
                                          ))
    
    output$input_tbl_3 <- renderDataTable(data(),
                                          rownames = FALSE,
                                          extensions = 'Buttons',
                                          options = list(
                                            paging = FALSE,
                                            searching = TRUE,
                                            ordering = TRUE,
                                            #scrollCollapse= TRUE,
                                            dom = 'l<"sep">Bfrtip',
                                            buttons = c('csv', 'excel')
                                          ))
    
    output$input_tbl_4 <- renderDataTable(data(),
                                          rownames = FALSE,
                                          extensions = 'Buttons',
                                          options = list(
                                            paging = FALSE,
                                            searching = TRUE,
                                            ordering = TRUE,
                                            #scrollCollapse= TRUE,
                                            dom = 'l<"sep">Bfrtip',
                                            buttons = c('csv', 'excel')
                                          ))
    
    
    
    output$ui.action <- renderUI({
      if (is.null(data())) return()
      tipify(actionButton("search_physchem", 
                          "Physchem Search",
                          icon = icon('atom'),
                          style = 'color: #fff; 
                   background-color: #a44f2e; border-color: #8a2b07'),
             'Searches ChEMBL and SusDat for physiochemical data',
             placement="bottom", trigger = "hover")
    })
    
    chembl_data <- eventReactive(input$search_physchem, {
        #query the chembl database
        CHEMBLSearch(data())
    })
    
    not_found_in_chembl <- eventReactive(input$search_physchem,{
        #extract the compounds not found in chembl
        NotInChEMBL(data(), chembl_data())
    })
    
    
    sus_data <- eventReactive(input$search_physchem,{
        #query the epi suite database
      SusdatSearch(data(), not_found_in_chembl(), chembl_data())
    })
    
    
    output$TBL2 <- renderDataTable(sus_data(),
                                   rownames = FALSE,
                                   extensions = 'Buttons',
                                   options = list(
                                     paging = FALSE,
                                     searching = TRUE,
                                     ordering = TRUE,
                                     #scrollCollapse= TRUE,
                                     dom = 'l<"sep">Bfrtip',
                                     buttons = c('csv', 'excel')
                                   ))
    
    output$physchem_tbl_2 <- renderDataTable(sus_data(),
                                             rownames = FALSE,
                                             extensions = 'Buttons',
                                             options = list(
                                               paging = FALSE,
                                               searching = TRUE,
                                               ordering = TRUE,
                                               #scrollCollapse= TRUE,
                                               dom = 'l<"sep">Bfrtip',
                                               buttons = c('csv', 'excel')
                                             ))
    
    output$physchem_tbl_3 <- renderDataTable(sus_data(),
                                             rownames = FALSE,
                                             extensions = 'Buttons',
                                             options = list(
                                               paging = FALSE,
                                               searching = TRUE,
                                               ordering = TRUE,
                                               #scrollCollapse= TRUE,
                                               dom = 'l<"sep">Bfrtip',
                                               buttons = c('csv', 'excel')
                                             ))
    

    nf_compounds <- eventReactive(c(input$search_physchem, input$include_compounds),{
        #extract the compounds not found in susdat
      MissingInformation(data(), not_found_in_chembl(), sus_data(), missing_info = input$include_compounds)
    })
    
    output$ACDLabs_box<-reactive(!is.null(nf_compounds()))
    outputOptions(output, "ACDLabs_box", suspendWhenHidden = FALSE)
    
    output$TBL6 <- renderDataTable(nf_compounds(),
                                   rownames = FALSE,
                                   extensions = 'Buttons',
                                   options = list(
                                     paging = FALSE,
                                     searching = TRUE,
                                     ordering = TRUE,
                                     #scrollCollapse= TRUE,
                                     dom = 'l<"sep">Bfrtip',
                                     buttons = c('csv', 'excel')
                                   ))
    
    acd_data <- eventReactive(c(input$refresh_physchem,input$upload_ACD_labs),{
      file3 <- input$file3
      req(file3)
      ACD_outputs(data(),file3$datapath,sus_data())
    })
    
    output$TBL7 <- renderDataTable(acd_data(),
                                   rownames = FALSE,
                                   extensions = 'Buttons',
                                   options = list(
                                     paging = FALSE,
                                     searching = TRUE,
                                     ordering = TRUE,
                                     #scrollCollapse= TRUE,
                                     dom = 'l<"sep">Bfrtip',
                                     buttons = c('csv', 'excel')
                                   ))
    
    output$expdummy<-reactive(!is.null(input$file2))
    outputOptions(output, "expdummy", suspendWhenHidden = FALSE)
    
    CAS_DTXSID <- eventReactive(c(input$search_exp,input$search_httk_db),{
      
      CAS_and_DTXSID(data())
    })
    
    physchem_data<-eventReactive(c(input$search_exp,input$search_httk_db),{

      file3 <- input$file3
      #print(file3)
      
      if(!is.null(file3)){
        
        return(acd_data())
        

      } else{
      
        return(sus_data())

      }
    })
    
    out_of_range_values <- eventReactive(input$search_physchem,{
      #find out of range Peff values
      OutOfRange_PSA_HBD(physchem_data())
    })
    
    output$TBL_out_of_range <- renderDataTable(out_of_range_values(),
                                               rownames = FALSE,
                                               extensions = 'Buttons',
                                               options = list(
                                                 paging = FALSE,
                                                 searching = TRUE,
                                                 ordering = TRUE,
                                                 #scrollCollapse= TRUE,
                                                 dom = 'l<"sep">Bfrtip',
                                                 buttons = c('csv', 'excel')
                                               ))
    
    output$TBL_out_of_range2 <- renderDataTable(out_of_range_values(),
                                               rownames = FALSE,
                                               extensions = 'Buttons',
                                               options = list(
                                                 paging = FALSE,
                                                 searching = TRUE,
                                                 ordering = TRUE,
                                                 #scrollCollapse= TRUE,
                                                 dom = 'l<"sep">Bfrtip',
                                                 buttons = c('csv', 'excel')
                                               ))
    
    output$TBL_out_of_range3 <- renderDataTable(out_of_range_values(),
                                                rownames = FALSE,
                                                extensions = 'Buttons',
                                                options = list(
                                                  paging = FALSE,
                                                  searching = TRUE,
                                                  ordering = TRUE,
                                                  #scrollCollapse= TRUE,
                                                  dom = 'l<"sep">Bfrtip',
                                                  buttons = c('csv', 'excel')
                                                ))

    httk_only_data <- eventReactive(input$search_httk_db,{
      #query the httk database
      httkSearch(physchem_data(),CAS_DTXSID(), info = data(), fu_operation = input$fu_op, CLint_operation = input$clint_op,
                 obach = input$obach_data, wambaugh = input$wambaugh_data, chem_phys_in_vitro = input$chemphys_data)
    })

    httk_data <- eventReactive(input$search_exp, {
        #query the httk database
        httkSearch(physchem_data(),CAS_DTXSID(), info = data(), fu_operation = input$fu_op, CLint_operation = input$clint_op,
                   obach = input$obach_data, wambaugh = input$wambaugh_data, chem_phys_in_vitro = input$chemphys_data)
    })
    
    
    httk_exp_data <- eventReactive(input$search_exp, {
        #incorporate findings from experimental data
        file2 <- input$file2
        req(file2)
        ExpDataSearch(httk_data = httk_data(), 
                      experimental_data_directory = file2$datapath,
                      CL_threshold = input$CL_thresh,
                      BP_threshold = input$BP_thresh,
                      fu_threshold = input$fu_thresh, 
                      mean_flag = input$exp_mean_flag)
    })
    
    #check if the experimental data is incorporated
    experimental_data<-eventReactive(c(input$search_exp,input$search_httk_db),{
      if(input$include_exp_data){
        return(httk_exp_data())
        
      } else{
        return(httk_only_data())
        
      }
    })
    
    
    #only if the experimental data has been collected, does the simulate button appear
    output$simulate_check<-reactive(!is.null(experimental_data()))
    outputOptions(output, "simulate_check", suspendWhenHidden = FALSE)
    
    #Display the collected information
    output$TBL3 <- renderDataTable(experimental_data(),
                                   rownames = FALSE,
                                   extensions = 'Buttons',
                                   options = list(
                                     paging = FALSE,
                                     searching = TRUE,
                                     ordering = TRUE,
                                     #scrollCollapse= TRUE,
                                     dom = 'l<"sep">Bfrtip',
                                     buttons = c('csv', 'excel')
                                   ))
    
    output$exp_dat_tbl_2 <- renderDataTable(experimental_data(),
                                   rownames = FALSE,
                                   extensions = 'Buttons',
                                   options = list(
                                     paging = FALSE,
                                     searching = TRUE,
                                     ordering = TRUE,
                                     #scrollCollapse= TRUE,
                                     dom = 'l<"sep">Bfrtip',
                                     buttons = c('csv', 'excel')
                                   ))
    
    output$determine_inputs <- reactive(check_some_columns(data()))
    outputOptions(output,"determine_inputs", suspendWhenHidden = FALSE)
    
    organised_data <- eventReactive(c(input$simulate_button,input$predict_params_button), {
      #organise the physiochemical and experimental data
      OrganiseInputData(experimental_data(), Vss_method = input$pred_method,
                        Input_Dose = input$dose_value, UNITS = input$dose_units,
                        info = data(), admin_route = input$administration_route)

    })
    
    observe({
      lol <<- organised_data()
    })
    
    predicted_variables <- eventReactive(input$predict_params_button,{
      #Get predicted for fu, BP, Vss and Kd
      PredictParameters(organised_data())
    })
    
    output$pred_output_tbl <- renderDataTable(predicted_variables(),
                                   rownames = FALSE,
                                   extensions = list('Buttons'= NULL,
                                                     'Scroller' = NULL),
                                   options = list(
                                     scroller = TRUE,
                                     scrollX = TRUE,
                                     scrollY = '300px',
                                     searching = TRUE,
                                     fixedColumns = TRUE,
                                     autoWidth = TRUE,
                                     #ordering = TRUE,
                                     scrollCollapse= TRUE,
                                     dom = 'tB',
                                     buttons = c('csv', 'excel')
                                   ))
    
    multiple_dosing_flag <- eventReactive(input$simulate_button, {
      if (input$dosing_options == 'Multiple Dosing'){
        return(T)
      }
    })
    
    output_profiles <- eventReactive(input$simulate_button, {
      #Run the simulations through Simcyp
      SimcypSimulation(organised_data(), subjects = input$subjects, #,trials = input$trials
                       Time = as.numeric(input$sim_time), seed = input$set_seed,
                       
                       #set subject parameters
                       MinAge = min(as.numeric(input$age_range)), 
                       MaxAge = max(as.numeric(input$age_range)), 
                       Prop_females = input$female_prop,
                       
                       #dermal values
                       dermal_area = input$dermal_area_val, 
                       formulation_thickness = input$dermal_formulation_thickness, 
                       formulation_density = input$dermal_formulation_density,
                       
                       #multiple dosing options
                       multiple_dosing = F,#multiple_dosing_flag(), there is an issue with multiple dosing
                       Num_doses = as.numeric(input$num_doses_val), 
                       dose_interval = as.numeric(input$dosing_interval_val)
                       )
    })
    
    
    #additional params
    #oliveoil_water_for_NLP = T
    
    # observe({
    #   obs_data <<- output_profiles()
    # })
    
    simcyp_outputs <- eventReactive(input$simulate_button, {
      #All data in simcyp
      AdditionalOutputs(organised_data())
    })

    summary_outputs <- eventReactive(input$simulate_button, {
      #Summary data in simcyp
      suppressWarnings(SummaryOutputs(simcyp_outputs()))
    })
    
    static_predicted_params <- eventReactive(input$simulate_button, {
      #Extract non-population dependant parameters
      StaticPredictedParameters(summary_outputs())
    })
    
    
    #these are the summary simcyp outputs and the static parameters in one place
    all_summary_outputs <- eventReactive(input$simulate_button, {
      fin_data <- cbind(summary_outputs(),rm_df_cols(static_predicted_params(),c('Code')))
      fin_data
    })
    
    #generate the compound codes after hitting the simulate button
    compound_codes <- eventReactive(input$simulate_button, {
      
      #identify the compound codes
      codes<- extract_column(summary_outputs(),1)
      
      #convert compound codes to a list
      compound_codes_list<-as.list(codes)
      names(compound_codes_list)<-codes
      
      return(compound_codes_list)
      
    })
    
    #generate the tissue types after hitting the simulate button
    tissue_types <- eventReactive(input$simulate_button, {
      
      #identify the compound codes
      tissues<- extract_tissue(output_profiles())
      
      #convert tissues to a list
      tissue_list<-as.list(tissues)
      names(tissue_list)<- tissues
      
      return(tissue_list)
      
    })
    
    output$TBL4 <- renderDataTable(output_profiles(),
                                   rownames = FALSE,
                                   extensions = list('Buttons'= NULL,
                                                     'Scroller' = NULL),
                                   options = list(
                                     scroller = TRUE,
                                     scrollX = TRUE,
                                     scrollY = '300px',
                                     searching = TRUE,
                                     fixedColumns = TRUE,
                                     autoWidth = TRUE,
                                     #ordering = TRUE,
                                     scrollCollapse= TRUE,
                                     dom = 'tB',
                                     buttons = c('csv', 'excel')
                                   ))
    
    output$TBL5 <- renderDataTable(summary_outputs(),
                                   rownames = FALSE,
                                   extensions = list('Buttons'= NULL,
                                                     'Scroller' = NULL),
                                   options = list(
                                     #deferRender = TRUE,
                                     #paging = TRUE,
                                     scroller = TRUE,
                                     scrollX = TRUE,
                                     scrollY = '300px',
                                     searching = TRUE,
                                     fixedColumns = TRUE,
                                     autoWidth = TRUE,
                                     #ordering = TRUE,
                                     scrollCollapse= TRUE,
                                     dom = 'tB',
                                     buttons = c('csv', 'excel')
                                   ))
    
    output$compound_lists<- renderUI({
      #generate drop_down lists of compound codes
      selectInput('comp_code', 'Compound Code',
                  choices = compound_codes())
    })
    
    output$compound_lists2<- renderUI({
      #generate drop_down lists of compound codes
      selectInput('comp_code2', 'Compound Code',
                  choices = compound_codes())
    })
    
    output$compound_lists3<- renderUI({
      #generate drop_down lists of compound codes
      selectInput('comp_code3', 'Compound Code',
                  choices = compound_codes())
    })
    
    output$tissue_lists<- renderUI({
      #generate drop_down lists of compound codes
      selectInput('tissue', 'Tissue Type',
                  choices = tissue_types(), selected = 'PLASMA')
    })
    
    ############################################################
    #                   Plotting Outputs                       #
    ############################################################
    
    output$conc_time_plot <- renderPlot({
      plot_profile(input$comp_code, output_profiles(), tissue_type = input$tissue,
           units=input$unit, curated_data = experimental_data(), 
           logy = input$log_scale)
    })
    
    #allow user to download outputs
    output$DownloadConcPlot<- downloadHandler(
      #Specify The File Name 
      filename = 'Concentration-Time Profile.png',
      content = function(file){
        # open the format of file which needs to be downloaded ex: pdf, png etc. 
        png(file, res = 300, height = 12, width = 15, units = "in")
        
        print(plot_profile(input$comp_code, output_profiles(), tissue_type = input$tissue,
                           units=input$unit, curated_data = experimental_data(), 
                           logy = input$log_scale))
        
        dev.off()
      }
    )
    
    output$comp_comparison_plot <- renderPlot({
      compare_simulated_compound(all_summary_outputs(),
                   parameter=input$sim_parameters, bar_col = input$plot_col_2, 
                   bar_order = input$param_order)
    })
    
    #allow user to download outputs
    output$CompComparPlot<- downloadHandler(
      #Specify The File Name 
      filename = 'Cross-compound Comparison Plot.png',
      content = function(file){
        # open the format of file which needs to be downloaded ex: pdf, png etc. 
        png(file, res = 300, height = 12, width = 15, units = "in")
        
        output$comp_comparison_plot <- renderPlot({
          compare_simulated_compound(all_summary_outputs(),
                                     parameter=input$sim_parameters, bar_col = input$plot_col_2, 
                                     bar_order = input$param_order)
        })
        
        dev.off()
      }
    )

    output$distribution_plot <- renderPlot({
      plot_parameters(simcyp_outputs(), Compound = input$comp_code2,
                      plot_type = 'Distribution',
                      x_variable = input$x_axis_var1,
                      y_variable = 'Age',
                      chosen_col = input$plot_col)
    })
    
    #allow user to download outputs
    output$DownloaddistPlot<- downloadHandler(
      #Specify The File Name 
      filename = 'Distribution Plot.png',
      content = function(file){
        # open the format of file which needs to be downloaded ex: pdf, png etc. 
        png(file, res = 300, height = 12, width = 15, units = "in")
        
        output$distribution_plot <- renderPlot({
          plot_parameters(simcyp_outputs(), Compound = input$comp_code2,
                          plot_type = 'Distribution',
                          x_variable = input$x_axis_var1,
                          y_variable = 'Age',
                          chosen_col = input$plot_col)
        })
        
        dev.off()
      }
    )
    
    output$relationship_plot <- renderPlot({
      plot_parameters(simcyp_outputs(), Compound = input$comp_code3,
                      plot_type = 'Relationship',
                      x_variable = input$x_axis_var2,
                      y_variable = input$y_axis_var,
                      chosen_col = input$plot_col2)
    })
    
    
    #allow user to download outputs
    output$DownloadrlnPlot<- downloadHandler(
      #Specify The File Name 
      filename = 'Relationship Plot.png',
      content = function(file){
        # open the format of file which needs to be downloaded ex: pdf, png etc. 
        png(file, res = 300, height = 12, width = 15, units = "in")
        
        output$relationship_plot <- renderPlot({
          plot_parameters(simcyp_outputs(), Compound = input$comp_code3,
                          plot_type = 'Relationship',
                          x_variable = input$x_axis_var2,
                          y_variable = input$y_axis_var,
                          chosen_col = input$plot_col2)
        })
        
        dev.off()
      }
    )
    
    
    output$physchem_prediction_comparison_plot <- renderPlot({
      PhyschemvsPredictedParamsPlot(organised_data(), all_summary_outputs(),
                                    predicted_param = input$pred_parameters,
                                    physchem_param = input$physchem_parameters,
                                    plot_colour = input$plot_col3)
    })
    

    #allow user to download outputs
    output$DownloadphyschemCompCmparPlot<- downloadHandler(
      #Specify The File Name 
      filename = 'Physchem vs. Predicted Param Plot.png',
      content = function(file){
        # open the format of file which needs to be downloaded ex: pdf, png etc. 
        png(file, res = 300, height = 12, width = 15, units = "in")
        
        output$physchem_prediction_comparison_plot <- renderPlot({
          PhyschemvsPredictedParamsPlot(organised_data(), all_summary_outputs(),
                                        predicted_param = input$pred_parameters,
                                        physchem_param = input$physchem_parameters,
                                        plot_colour = input$plot_col3)
        })
        
        dev.off()
      }
    )
    
    ############################################################
    #                    SIVA Outputs                          #
    ############################################################
    
    SIVA_inputs <- reactive({
      SIVA_input_file <- input$SIVA_file
      req(SIVA_input_file)
      SIVA_input_file$datapath
    })
    
    SIVA_collected_data <- eventReactive(input$search_SIVA_data, {
      #incorporate findings from experimental data
      SIVA_inputs <- input$SIVA_file
      req(SIVA_inputs)
      SIVA_data_collection(filename = SIVA_inputs$datapath)
    })
    
    #inform user of any missing data in the description blocks
    output$missing_comp_num <- renderText({ Missing_compound_check(SIVA_inputs(),SIVA_collected_data()) })
    output$missing_hlc_num <- renderText({ Missing_HLC_check(SIVA_collected_data()) })
    output$duplicated_comp_num <- renderText({ duplicated_compounds_check(SIVA_collected_data()) })
    
    SIVA_updated_data <- eventReactive(input$update_SIVA_data,{
      SIVA_input_file <- input$SIVA_file_reupload
      req(SIVA_input_file)
      SIVA_updated_file_path <- SIVA_input_file$datapath
      
      if (file.exists(SIVA_updated_file_path)) {
        
        #extract file extensions
        extension <- file_ext(SIVA_updated_file_path)
        
      } else {
        
        return(cat("The file does not exist"))
        
      }
      
      #load the data based on the extracted extension
      if (extension == 'xlsx'){
        
        import<-loadWorkbook(SIVA_updated_file_path, create=FALSE)
        ICD <- readWorksheet(import,  header = TRUE, sheet = 1)
        
      } else if (extension == 'csv'){
        
        ICD <- read.csv(SIVA_updated_file_path, header = T, sep = ',')
        
      } else {
        
        return(cat("File type not compatible."))
        
      }
      ICD
    })
    
    SIVA_curated_data_final<-eventReactive(c(input$search_SIVA_data,input$update_SIVA_data),{
      
      SIVA_reupload <- input$SIVA_file_reupload
      #print(file3)
      
      if(!is.null(SIVA_reupload)){
        
        return(SIVA_updated_data())
        
        
      } else{
        
        return(SIVA_collected_data())
        
      }
    })
    
    #Display the collected information
    output$SIVA_curated_data_tbl <- renderDataTable(SIVA_curated_data_final(),
                                                    rownames = FALSE,
                                                    extensions = 'Buttons',
                                                    options = list(
                                                      paging = FALSE,
                                                      searching = TRUE,
                                                      ordering = TRUE,
                                                      #scrollCollapse= TRUE,
                                                      dom = 'l<"sep">Bfrtip',
                                                      buttons = c('csv', 'excel')
                                                    ))
    
    observe({
      obs_data <<- SIVA_curated_data_final()
    })
    
    
    #Calculate the SIVA predictions
    SIVA_predictions <- eventReactive(input$predict_iv_distribution, {

      InVitroDistribution(ICD = SIVA_curated_data_final(), 
                          Cinitial = input$nominal_conc, ncell = input$number_cells, 
                          Vwell = input$vol_well, pH_medium = input$medium_pH,
                          f_serum =input$serum_fraction, Sysdiameter = input$well_diameter,
                          save.csv = F)
    })
    
    #Display the predicted outputs
    output$SIVA_output_tbl <- renderDataTable(SIVA_predictions(),
                                                    rownames = FALSE,
                                                    extensions = 'Buttons',
                                              options = list(
                                                paging = FALSE,
                                                searching = TRUE,
                                                ordering = TRUE,
                                                #scrollCollapse= TRUE,
                                                dom = 'l<"sep">Bfrtip',
                                                buttons = c('csv', 'excel')
                                              ))
    
    #plot the SIVA output plots
    output$SIVA_output_plots <- renderPlot({
      SIVA_plots(SIVA_collected_data(), SIVA_predictions())
    })
    
    #allow user to download SIVA outputs
    output$Cat10<- downloadHandler(
      #Specify The File Name 
      filename = 'In Vitro Distribution Plots.png',
      content = function(file){
        # open the format of file which needs to be downloaded ex: pdf, png etc. 
        png(file, res = 300, height = 4, width = 15, units = "in")
        
        print(SIVA_plots(SIVA_collected_data(), SIVA_predictions()))
        
        dev.off()
      }
    )
    
    ############################################################
    #                       Help Page                          #
    ############################################################
    
    output$HelpPageText <- renderUI ({

      rawText<-readLines('README.txt')
      replacedText <- lapply(rawText, determine_line)
      
      return(replacedText)
    })
    
}

shinyApp(ui = ui, server = server)
