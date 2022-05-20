library(ggplot2)
library(scales)

# Function from plotting conc-time profiles for each compound
plot_profile <- function(casestudy_ID, Output, tissue_type = 'PLASMA' ,units = 'ng/mL', curated_data, logy=F){
  
  require(ggplot2)
  require(scales)
  
  #extract a df for the CSID
  indicies<-which(Output$Group == casestudy_ID)
  df <- Output[indicies,2:ncol(Output)]
  
  #determine the correct tissue type
  tissue_type <- toupper(tissue_type)
  tissue_index <- str_detect(colnames(df),tissue_type)
  time_index <- str_detect(colnames(df),'time|Time|times|Times')
  set1 <- which(tissue_index==T)
  set2 <- which(time_index == T)
  df <- df[,c(set2,set1)]
  
  #adjust column names
  colnames(df) <- sub("_[^_]+$", "", colnames(df))
  
  #if more than 30 individuals are simulated, only select 30 at random
  if (ncol(df)>31){
    randomly_selected_cols <- sample(2:ncol(df), 2, replace=F)
    df<-df[,c(1,randomly_selected_cols)]
  }
  
  #convert dataframe into long form
  df<- long_form(df)
  colnames(df)<-c('Sub','Times','Conc')
  
  #Create Title for plot
  header<-paste(casestudy_ID, tissue_type,'Concentration-Time Profiles',sep=' ')
  
  if (units == 'ng/mL'){
    
    df <- df
    #mean_profile$Concentration <- mean_profile$Concentration
    ytitle <- paste('Concentration',paste0('(',units,')'), sep=' ')
    
  } else if (units == 'uM'){
    
    #Extract MW from curated data
    df_curated_data <- curated_data[which(curated_data$Code == casestudy_ID),]
    MW <- df_curated_data$MW
    
    ## convert units of cmax from mg/L -> uM #check again
    df$Conc <- (df$Conc/1000)/MW*1e6
    ytitle <- paste('Concentration',paste0('(',units,')'), sep=' ')
    
  } else{
    
    message('Please select either ng/mL or uM. ng/mL used as default')
    ytitle <- paste('Concentration (ng/mL)', sep=' ')
  }
  
  if (logy == T){
    
    #Plot the mean profile (y axis normal scale)
    g1 <- ggplot(df, aes(x = Times, y = Conc, group = factor(Sub))) + geom_line(aes(color=factor(Sub)), size = 1.5, linetype = 1)
    g1 <- g1 + geom_line() 
    g1 <- g1 + xlab("Time (hours)") + ylab(ytitle)
    g1 <- g1 + ggtitle(header) + 
      labs(color  = "Subjects", linetype = 1)+ theme(
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey"),
        axis.text = element_text(colour = "black",size=12),
        axis.title = element_text(colour = "black",size=14,face="bold")
      )
    
    g1 + scale_y_continuous(trans='log10')
    
  } else {
    
    #Plot the mean profile (y axis normal scale)
    g1 <- ggplot(df, aes(x = Times, y = Conc, group = factor(Sub))) + geom_line(aes(color=factor(Sub)), size = 1, linetype = 1) 
    g1 <- g1 + xlab("Time (hours)") + ylab(ytitle)
    g1 + ggtitle(header) + 
      labs(color  = "Subjects", linetype = 1)+ theme(
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        panel.grid.major = element_line(color = "grey"),
        panel.grid.minor = element_line(color = "grey"),
        axis.text = element_text(colour = "black",size=12),
        axis.title = element_text(colour = "black",size=14,face="bold")
      )
  }
  
}

plot_parameters<-function(simcyp_outputs, Compound, 
                          plot_type = 'Distribution', x_variable, y_variable,
                          chosen_col){
  require(ggplot2)
  #Extract data relationg to compound of interest
  compound_IDs<- names(simcyp_outputs)
  index <- which(Compound==compound_IDs)
  outputs<- simcyp_outputs[[index]]
  
  #determine x and y variables
  x_col<- which(colnames(outputs)==x_variable)
  y_col <- which(colnames(outputs)== y_variable)
  
  x_axis_label <- paste(x_variable,determine_units(x_variable), sep=' ')
  
  if (plot_type == 'Distribution'){
    #Create Title for plot
    header<-paste(Compound,'Distribution of', x_variable, sep=' ')
    
    # Distributions are plotted as density plots (only x variable is considered)
    ggplot(outputs, aes(x=outputs[,x_col], colour= chosen_col, fill = chosen_col)) +
      geom_density(alpha=0.4)+ ggtitle(header)+ xlab(x_axis_label)+
      scale_color_manual(values= chosen_col) +
      scale_fill_manual(values= chosen_col) +
      theme(legend.position = "None",
            plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA),
            axis.line.x = element_line(color="black", size = 1),
            axis.line.y = element_line(color="black", size = 1),
            panel.grid.major = element_line(color = "grey"),
            panel.grid.minor = element_line(color = "grey"),
            axis.text = element_text(colour = "black",size=12),
            axis.title = element_text(colour = "black",size=14,face="bold"))
    
  } else if (plot_type == 'Relationship'){
    
    #create y axis label
    y_axis_label <- paste(y_variable, determine_units(y_variable), sep=' ')
    
    #find the correlation coefficient
    cor_coeff<- round(cor(outputs[,x_col],outputs[,y_col]),2)
    
    #Create Title for plot
    header<-paste(Compound,' Relationship between ', x_variable,' and ', y_variable,' (r = ',cor_coeff,')' ,sep='')
    
    #relationships are scatter diagrams where x and y variables are needed
    g1 <- ggplot(outputs, aes(x = outputs[,x_col], y = outputs[,y_col]))
    g1 <- g1 + geom_point(colour = chosen_col) + geom_smooth(method = "lm", se = FALSE, colour = chosen_col)  
    g1 <- g1 + xlab(x_axis_label) + ylab(y_axis_label)
    g1 + ggtitle(header) + theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line.x = element_line(color="black", size = 1),
      axis.line.y = element_line(color="black", size = 1),
      panel.grid.major = element_line(color = "grey"),
      panel.grid.minor = element_line(color = "grey"),
      axis.text = element_text(colour = "black",size=12),
      axis.title = element_text(colour = "black",size=14,face="bold")
    )
  }
  
}

#cross-compound comparison plot
compare_simulated_compound <- function(summary_simcyp, parameter, bar_order, bar_col){
  
  compound_codes<-summary_simcyp$CS_code
  parameters_to_plot <- keep_df_cols(summary_simcyp,parameter)
  
  df<-data.frame(compound_codes,parameters_to_plot)
  colnames(df)<-c('Compounds','parameter')
  
  if(bar_order=='ascending'){
    df <- df[order(df$parameter),]
    
  } else if (bar_order == 'compound_code'){
    df <- df[order(df$Compounds),]
  }
  
  header<-paste('Mean',parameter,'Variation Across Simulated Compounds',sep=' ')
  
  ggplot(df, aes(x=Compounds, y= parameter)) +
    ylab(parameter)+
    geom_bar(stat="identity", colour = bar_col, fill = bar_col)+
    ggtitle(header)+ 
    scale_x_discrete(limits = df$Compounds) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line.x = element_line(color="black", size = 1),
      axis.line.y = element_line(color="black", size = 1),
      panel.grid.major = element_line(color = "grey"),
      panel.grid.minor = element_line(color = "grey"),
      axis.text = element_text(colour = "black",size=12),
      axis.title = element_text(colour = "black",size=14,face="bold"))
}