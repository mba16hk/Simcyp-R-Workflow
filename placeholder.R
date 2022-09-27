#set the working directory to the working directory of the scripts
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

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
source('Additional_data.R')
source('plotting_functions.R')

#set the file directory of the compound file

#without dose
#file_dir <- 'C:/Users/hkhalidi/OneDrive - Certara/Desktop/R Workflow/Paper/compound_file.xlsx'

#with dose
#file_dir <- 'C:/Users/hkhalidi/OneDrive - Certara/Desktop/R Workflow/data_files/Small_Compound_File.xlsx'

#RISKHUNTR Data
#file_dir <- 'C:/Users/hkhalidi/OneDrive - Certara/Desktop/RISKHUNT3R/WP4/WP4_compounds_fu_search.xlsx'

file_dir <- "C:/Users/hkhalidi/OneDrive - Certara/Desktop/RISKHUNT3R/WP4/WP4_comps.xlsx"

#preprocess the compound file
data <- ProcessInputs(file_dir)

#query the chembl database
chembl_data <- CHEMBLSearch(data)

#determine the compounds not found in chembl
nf_in_chembl <- NotInChEMBL(data, chembl_data)

#query the Norman suspect database for the compounds not found in Chembl
sus_data <- SusdatSearch(data, nf_in_chembl, chembl_data)

#determine the compounds not found in the Norman suspect list
NOT_FOUND <- NotFoundInsusdat(nf_in_chembl, sus_data)

#lol
out_of_range_peff <- OutOfRange_PSA_HBD(sus_data)
MissingInformation(data, nf_in_chembl, sus_data)
# ACD_data_directory<- 'C:/Users/hkhalidi/OneDrive - Certara/Desktop/R Workflow/data_files/acd_output.xls'
# acd_data<-ACD_inputs(data,nf_in_chembl,sus_data,missing_info=T)
# physchem_data <- ACD_outputs(data,ACD_data_directory,sus_data)

physchem_data <- sus_data
#extract CAS numbers and DTXSID values to query the httk database
CAS_DTXSID <- CAS_and_DTXSID(data)

httk_data <- httkSearch(physchem_data, CAS_DTXSID, data, fu_operation = 'arithmetic mean', CLint_operation = 'arithmetic mean')

#write.csv(httk_data,'WP4_compounds_fu_collected.csv',row.names = F)
#experimental_data_directory<-'C:/Users/hkhalidi/OneDrive - Certara/Desktop/R Workflow/data_files/Exemplar_Experimental_file.xlsx'
#httk_exp_data <- ExpDataSearch(httk_data, experimental_data_directory, CL_threshold = 3.8)


organised_data <- OrganiseInputData(httk_data,info = data, admin_route = 'Oral', UNITS = 'mg', Input_Dose = 10, AGP_pKa_threshold = 7)
#hello <- organised_data[1:3,]
#preds<- PredictParameters(organised_data)

output_profiles<-SimcypSimulation(organised_data, trials = 1, subjects = 2, Time = 3, seed = T, multiple_dosing = F)
simcyp_outputs <- AdditionalOutputs(organised_data)
summary_simcyp <-SummaryOutputs(simcyp_outputs)
predicted_params <- StaticPredictedParameters(summary_simcyp)

write.csv(summary_simcyp,'SimRFlow_cmax_HT.csv',row.names = F)

#search for number doses variable
simcyp_dbs <- list(AmountTypeID,CategoryID,CompoundID,CompoundParameterID,CompoundSummaryStatID, CompoundTypeID, DistributionTypePDID,DoseID, DosingRoute, IndividualID, IndividualValueID, InputTypePDID, IVPTID, PDID, PDResultID, PopulationID, PopulationParameterID, PrimaryCompartmentID, ProfileID, ResultID, SecondaryCompartmentID, SheetGuidsFor19, SheetID, SimulationParameterID, SolverSettingsID, SolverTypeID, SpeciesID, TertiaryCompartmentID,WorkspaceID)

for (i in 1:length(simcyp_dbs)){
  print(i)
  print(FindEnum(simcyp_dbs[[i]],'probfemale'))
}



dev.off()
tissues <- c('PLASMA', 'BRAIN', 'ADIPOSE','LIVER')
for (i in tissues){
  plot_name <- paste0('conc-time_plot_',i,'.pdf')
  
  pdf(plot_name,         # File name
      width = 6, height = 4)
  
  print(plot_profile('A1',output_profiles, i,curated_data = httk_data, units = 'uM', logy=F))
  
  dev.off()
}


params <- c('Vss', 'Fu_plasma', 'AccumulationIndex','HalfLife')
cols <- c('indianred','royalblue','goldenrod1','darkolivegreen3')
for (i in 1:length(cols)){
  plot_name <- paste0('comparison_plot_',params[i],'.pdf')
  
  pdf(plot_name,         # File name
      width = 6, height = 3)
  
  print(compare_simulated_compound(summary_simcyp , params[i],'compound_code',cols[i]))
  
  dev.off()
}

params1 <- c('Vss','HalfLife')
params2 <- c('Age','BSA')
cols <- c('indianred','royalblue')
for (i in 1:length(cols)){
  plot_name <- paste0('relationship_plot_',params[i],'.pdf')
  
  pdf(plot_name,         # File name
      width = 5, height = 4)
  
  print(plot_parameters(simcyp_outputs, 'A1', 
                  plot_type = 'Relationship', params2[i], params1[i],
                  cols[i]))
  
  dev.off()
}

params2 <- c('AUC','BP')
cols <- c('goldenrod1','darkolivegreen3')
for (i in 1:length(cols)){
  plot_name <- paste0('dist_plot_',params2[i],'.pdf')
  
  pdf(plot_name,         # File name
      width = 5, height = 4)
  
  print(plot_parameters(simcyp_outputs, 'A1', 
                        plot_type = 'Distribution', params2[i], params1[i],
                        cols[i]))
  
  dev.off()
}

### CONC-TIME PROFILE

library(grid)
library(patchwork)
library(gridExtra)
library(ggpubr)


g1 <- plot_profile('A1',output_profiles,'PLASMA',curated_data = httk_data, units = 'ng/mL', logy=F, compound_name = paste(sprintf('\u2191'),'Data (A2)', sep=' '))
g2 <- plot_profile('A1',output_profiles,'GUT',curated_data = httk_data, units = 'ng/mL', logy=F, compound_name = paste(sprintf('\u2191'),'Data (A2)', sep=' '))
g3 <- plot_profile('A1',output_profiles,'BRAIN',curated_data = httk_data, units = 'ng/mL', logy=F, compound_name = paste(sprintf('\u2191'),'Data (A2)', sep=' '))
g13 <- plot_profile('A9',output_profiles,'PLASMA',curated_data = httk_data, units = 'ng/mL', logy=F, compound_name = paste('- Data (A9)', sep=' '))
g14 <- plot_profile('A9',output_profiles,'GUT',curated_data = httk_data, units = 'ng/mL', logy=F, compound_name = paste('- Data (A9)', sep=' '))
g15 <- plot_profile('A9',output_profiles,'BRAIN',curated_data = httk_data, units = 'ng/mL', logy=F, compound_name = paste('- Data (A9)', sep=' '))
g4 <- plot_profile('A4',output_profiles,'PLASMA',curated_data = httk_data, units = 'ng/mL', logy=F, compound_name = paste(sprintf('\u2193'),'Data (A4)', sep=' '))
g5 <- plot_profile('A4',output_profiles,'GUT',curated_data = httk_data, units = 'ng/mL', logy=F, compound_name = paste(sprintf('\u2193'),'Data (A4)', sep=' '))
g6 <- plot_profile('A4',output_profiles,'BRAIN',curated_data = httk_data, units = 'ng/mL', logy=F, compound_name = paste(sprintf('\u2193'),'Data (A4)', sep=' '))


figure <- ggarrange(g1 +rremove("ylab")  +rremove("xlab"), g13 +rremove("ylab")  +rremove("xlab"), g4 +rremove("ylab")  +rremove("xlab"), g2+rremove("ylab")  +rremove("xlab"),g14+rremove("ylab")  +rremove("xlab"),g5+rremove("ylab")  +rremove("xlab"),g3+rremove("ylab")  +rremove("xlab"), g15+rremove("ylab")  +rremove("xlab"), g6+rremove("ylab")  +rremove("xlab"),
                    labels = c("(A)", "(B)", "(C)","(D)","(E)","(F)", "(G)","(H)","(I)"),
                    ncol = 3, nrow = 3,
                    common.legend = TRUE, legend = "none",
                    align = "hv")
                    
image_1 <- annotate_figure(figure, left = textGrob("Concentration (ng/mL)", rot = 90, vjust = 1, gp = gpar(cex = 1)),
                bottom = textGrob("Time (hours)", gp = gpar(cex = 1)))

image_1 

tiff('conc-time-profs.tiff', 
     width = 3000, height = 2500,
     res = 300, pointsize = 2)

print(image_1)

dev.off()
####################################
require(datasets)
require(grDevices)
require(graphics)


#### PLOT RELATIONSHIP PLOTS

g1 <- plot_parameters(simcyp_outputs, 'A2', 
             plot_type = 'Relationship', 'BSA', 'HalfLife','royalblue', compound_name = paste(sprintf('\u2191'),'Data (A2)', sep=' '))
g2 <- plot_parameters(simcyp_outputs, 'A2', 
                plot_type = 'Relationship', 'Age', 'Vss', 'indianred', compound_name = paste(sprintf('\u2191'),'Data (A2)', sep=' '))
g3 <- plot_parameters(simcyp_outputs, 'A2', 
                plot_type = 'Relationship', 'BW', 'Tmax','darkolivegreen3', compound_name = paste(sprintf('\u2191'),'Data (A2)', sep=' '))

g4 <- plot_parameters(simcyp_outputs, 'A9', 
                plot_type = 'Relationship', 'BSA', 'HalfLife','royalblue', compound_name = paste('- Data (A9)', sep=' '))
g5 <- plot_parameters(simcyp_outputs, 'A9', 
                plot_type = 'Relationship', 'Age', 'Vss', 'indianred', compound_name = paste('- Data (A9)', sep=' '))
g6 <- plot_parameters(simcyp_outputs, 'A9', 
                plot_type = 'Relationship', 'BW', 'Tmax','darkolivegreen3', compound_name = paste('- Data (A9)', sep=' '))

g7 <- plot_parameters(simcyp_outputs, 'A4', 
                plot_type = 'Relationship', 'BSA', 'HalfLife','royalblue', compound_name = paste(sprintf('\u2193'),'Data (A4)', sep=' '))
g8 <- plot_parameters(simcyp_outputs, 'A4', 
                plot_type = 'Relationship', 'Age', 'Vss', 'indianred', compound_name = paste(sprintf('\u2193'),'Data (A4)', sep=' '))
g9 <- plot_parameters(simcyp_outputs, 'A4', 
                plot_type = 'Relationship', 'BW', 'Tmax','darkolivegreen3', compound_name = paste(sprintf('\u2193'),'Data (A4)', sep=' '))

figure <- ggarrange(g1, g4 +rremove("ylab"), g7 +rremove("ylab") , g2 , g5+rremove("ylab"), g8+rremove("ylab"), g3, g6+rremove("ylab"), g9+rremove("ylab"),
                    labels = c("(A)", "(B)", "(C)","(D)","(E)","(F)", "(G)","(H)","(I)"),
                    ncol = 3, nrow = 3,
                    common.legend = TRUE, legend = "none",
                    align = "hv")
figure

#cairo_pdf('relationship_plots.pdf',         # File name
#    width = 12, height = 10)

tiff('relationship_plots.tiff', 
     width = 4000, height = 2500,
     res = 300, pointsize = 2)

print(figure)

dev.off()

################################################################

#### PLOT DISTRIBUTION PLOTS

g1 <- plot_parameters(simcyp_outputs, 'A2', 
                      plot_type = 'Distribution', 'HalfLife', 'HalfLife','royalblue', compound_name = paste(sprintf('\u2191'),'Data (A2)', sep=' '))
g2 <- plot_parameters(simcyp_outputs, 'A2', 
                      plot_type = 'Distribution', 'Vss', 'Vss', 'indianred', compound_name = paste(sprintf('\u2191'),'Data (A2)', sep=' '))
g3 <- plot_parameters(simcyp_outputs, 'A2', 
                      plot_type = 'Distribution', 'BP', 'BP','darkolivegreen3', compound_name = paste(sprintf('\u2191'),'Data (A2)', sep=' '))
g10 <- plot_parameters(simcyp_outputs, 'A2', 
                      plot_type = 'Distribution', 'Cmax_LUNG', 'Cmax_LUNG','goldenrod1', compound_name = paste(sprintf('\u2191'),'Data (A2)', sep=' '))

g4 <- plot_parameters(simcyp_outputs, 'A9', 
                      plot_type = 'Distribution', 'HalfLife', 'HalfLife','royalblue', compound_name = paste('- Data (A9)', sep=' '))
g5 <- plot_parameters(simcyp_outputs, 'A9', 
                      plot_type = 'Distribution', 'Vss', 'Vss', 'indianred', compound_name = paste('- Data (A9)', sep=' '))
g6 <- plot_parameters(simcyp_outputs, 'A9', 
                      plot_type = 'Distribution', 'BP', 'BP','darkolivegreen3', compound_name = paste('- Data (A9)', sep=' '))
g11 <- plot_parameters(simcyp_outputs, 'A9', 
                       plot_type = 'Distribution', 'Cmax_LUNG', 'Cmax_LUNG','goldenrod1', compound_name = paste('- Data (A9)', sep=' '))

g7 <- plot_parameters(simcyp_outputs, 'A4', 
                      plot_type = 'Distribution', 'HalfLife', 'HalfLife','royalblue', compound_name = paste(sprintf('\u2193'),'Data (A4)', sep=' '))
g8 <- plot_parameters(simcyp_outputs, 'A4', 
                      plot_type = 'Distribution', 'Vss', 'Vss', 'indianred', compound_name = paste(sprintf('\u2193'),'Data (A4)', sep=' '))
g9 <- plot_parameters(simcyp_outputs, 'A4', 
                      plot_type = 'Distribution', 'BP', 'BP','darkolivegreen3', compound_name = paste(sprintf('\u2193'),'Data (A4)', sep=' '))
g12 <- plot_parameters(simcyp_outputs, 'A4', 
                       plot_type = 'Distribution', 'Cmax_LUNG', 'Cmax_LUNG','goldenrod1', compound_name = paste(sprintf('\u2193'),'Data (A4)', sep=' '))


figure <- ggarrange(g1, g4 +rremove("ylab"), g7 +rremove("ylab") , g2 , g5+rremove("ylab"), g8+rremove("ylab"), g3, g6+rremove("ylab"), g9+rremove("ylab"), g10, g11+rremove("ylab"), g12+rremove("ylab"),
                    labels = c("(A)", "(B)", "(C)","(D)","(E)","(F)", "(G)","(H)","(I)","(J)","(K)","(L)"),
                    ncol = 3, nrow = 4,
                    common.legend = TRUE, legend = "none",
                    align = "hv")
figure

#cairo_pdf('distribution_plots.pdf',         # File name
#          width = 12.5, height = 12)

tiff('distribution_plots.tiff', 
     width = 3650, height = 2800,
     res = 300, pointsize = 2)

print(figure)

dev.off()

#######################################

#### PLOT COMPOUND COMPARISON PLOTS


g1 <- compare_simulated_compound(summary_simcyp ,'Vss','ascending','indianred')
g2 <-compare_simulated_compound(summary_simcyp ,'HalfLife','ascending','royalblue')
g3 <-compare_simulated_compound(summary_simcyp ,'CLtot','ascending','goldenrod1')

g5 <-compare_simulated_compound(summary_simcyp ,'Vss','comp_code','indianred')
g6 <-compare_simulated_compound(summary_simcyp ,'HalfLife','comp_code','royalblue')
g7 <-compare_simulated_compound(summary_simcyp ,'CLtot','comp_code','goldenrod1')

figure <- ggarrange(g5+ rremove("xlab"), g1+ rremove("ylab") +rremove("xlab"), g6+ rremove("xlab"), g2+ rremove("ylab") +rremove("xlab"), g7+ rremove("xlab"), g3+rremove("ylab") + rremove("xlab"), # remove axis labels from plots
                    labels = c("(A)", "(B)", "(C)","(D)","(E)","(F)"),
                    ncol = 2, nrow = 3,
                    common.legend = TRUE, legend = "none",
                    align = "hv")
#font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "bottom"))
image_1 <- annotate_figure(figure, bottom = textGrob("Compounds", gp = gpar(cex = 7)))

#pdf('comparison-plots.pdf',         # File name
#    width = 12, height = 9)

tiff('comparison-plots.tiff', 
     width = 3500, height = 2500,
     res = 300, pointsize = 2)

print(image_1)

dev.off()


########################################################################
