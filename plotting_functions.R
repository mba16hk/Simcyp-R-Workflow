library(ggplot2)
library(scales)

# Function from plotting conc-time profiles for each compound
plot_profile <- function(casestudy_ID, Output, tissue_type = 'PLASMA' ,units = 'ng/mL',
                         curated_data, logy=F, compound_name = NA, CI = F){
  
  require(ggplot2)
  require(scales)
  
  #extract a df for the CSID
  indicies<-which(Output$Group == casestudy_ID)
  df <- Output[indicies,2:ncol(Output)]
  
  if (tissue_type == 'ALL'){
    
    all_tissues <- c("plasma",
                     "brain",
                     "liver",
                     "kidney",
                     "skin",
                     "pancreas",
                     "gut",
                     "heart",
                     "lung",
                     "spleen",
                     "muscle",
                     "adipose")
    
    list_of_plots <- list()
    
    general_theme <- theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line.x = element_line(color="black", size = 1),
      axis.line.y = element_line(color="black", size = 1),
      panel.grid.major = element_line(color = "grey"),
      panel.grid.minor = element_line(color = "grey"),
      axis.text = element_text(colour = "black",size=9),
      axis.title = element_text(colour = "black",size=12,face="bold"),
      legend.key=element_blank()
    )
    
    for (i in 1:length(all_tissues)){
      
      tissue_type <- toupper(all_tissues[i])
      tissue_index <- str_detect(colnames(df),tissue_type)
      time_index <- str_detect(colnames(df),'time|Time|times|Times')
      set1 <- which(tissue_index==T)
      set2 <- which(time_index == T)
      tmp_df <- df[,c(set2,set1)]
      
      #adjust column names
      colnames(tmp_df) <- sub("_[^_]+$", "", colnames(tmp_df))
      
      
      if (CI == F){
        ## In this case, we will not be plotting confidence intervals
        
        #if more than 30 individuals are simulated, only select 30 at random
        if (ncol(tmp_df)>31){
          randomly_selected_cols <- sample(2:ncol(tmp_df), 2, replace=F)
          tmp_df<-tmp_df[,c(1,randomly_selected_cols)]
        }
        
        #convert dataframe into long form
        tmp_df<- long_form(tmp_df)
        colnames(tmp_df)<-c('Sub','Times','Conc')
        
      } else if (CI == T){
        
        ### we will be plotting confidence intervals
        
        vals <- tmp_df[,2:ncol(tmp_df)]
        confidence_intervals <- apply(vals,1,function (x) {
          c(lower = quantile(x, 0.05), upper = quantile(x, 0.95))}
        )
        confidence_intervals <- as.data.frame(t(confidence_intervals))
        tmp_df <- cbind(tmp_df$Times,rowMeans(vals),confidence_intervals)
        colnames(tmp_df)<-c('Times','Conc','5%','95%')
        
      }
      
      
      if (units == 'ng/mL'){
        
        tmp_df <- tmp_df
        #mean_profile$Concentration <- mean_profile$Concentration
        ytitle <- paste('Concentration',paste0('(',units,')'), sep=' ')
        
      } else if (units == 'uM'){
        
        #Extract MW from curated data
        df_curated_data <- curated_data[which(curated_data$CODE == casestudy_ID),]
        MW <- df_curated_data$MW
        
        ## convert units of cmax from mg/L -> uM #check again
        tmp_df$Conc <- (tmp_df$Conc/1000)/as.numeric(MW)*1e6
        
        if (CI == T){
          tmp_df$`5%`<- (tmp_df$`5%`/1000)/as.numeric(MW)*1e6
          tmp_df$`95%`<- (tmp_df$`95%`/1000)/as.numeric(MW)*1e6
        }
        
        ytitle <- paste('Concentration',paste0('(',units,')'), sep=' ')
        
      } else{
        
        message('Please select either ng/mL or uM. ng/mL used as default')
        ytitle <- paste('Concentration (ng/mL)', sep=' ')
      }
      
      
      if (CI == F){
        g1 <- ggplot(tmp_df, aes(x = Times, y = Conc, group = factor(Sub))) + 
          geom_line(aes(color=factor(Sub)), size = 1.5, linetype = 1) +  
          xlab(NULL) + ylab(NULL) + ggtitle(toupper(all_tissues[i])) + labs(color  = "Subjects", linetype = 1)
      } else if (CI == T){
        g1 <- ggplot(data = tmp_df) + 
          geom_line(aes(x = Times, y = Conc), size = 1.5, linetype = 1, colour = 'green3') + 
          geom_ribbon(aes(x = Times, ymin= `5%`, ymax= `95%`), alpha=0.2) +  
          xlab(NULL) + ylab(NULL) + ggtitle(toupper(all_tissues[i])) 
      }
      
      
      if (logy == T){
      
        g1 <- g1 + general_theme +  scale_y_continuous(trans='log10')
        
        list_of_plots[[i]] <- g1
        
      } else {

        #Plot the mean profile (y axis normal scale)
        g1 <- g1 +  general_theme
        
        list_of_plots[[i]] <- g1
      }
      
    }
    
    require(ggpubr)
    require(grid)
    
    all_tissues_fig <- ggarrange(plotlist = list_of_plots, ncol = 6, nrow = 2,common.legend = T, legend = 'top')
    annotate_figure(all_tissues_fig, left = textGrob(ytitle, rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                    bottom = textGrob("Time (hrs)", gp = gpar(cex = 1.3)))
    
  } else {
    
    
    general_theme <- theme(
          plot.background = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          axis.line.x = element_line(color="black", size = 1),
          axis.line.y = element_line(color="black", size = 1),
          panel.grid.major = element_line(color = "grey"),
          panel.grid.minor = element_line(color = "grey"),
          axis.text = element_text(colour = "black",size=12),
          axis.title = element_text(colour = "black",size=14,face="bold"),
          legend.key=element_blank()
        )
    
    #determine the correct tissue type
    tissue_type <- toupper(tissue_type)
    tissue_index <- str_detect(colnames(df),tissue_type)
    time_index <- str_detect(colnames(df),'time|Time|times|Times')
    set1 <- which(tissue_index==T)
    set2 <- which(time_index == T)
    df <- df[,c(set2,set1)]
    
    #adjust column names
    colnames(df) <- sub("_[^_]+$", "", colnames(df))
    
    if (CI == F){
      ## In this case, we will not be plotting confidence intervals
      
      #if more than 30 individuals are simulated, only select 30 at random
      if (ncol(df)>31){
        randomly_selected_cols <- sample(2:ncol(df), 2, replace=F)
        df<-df[,c(1,randomly_selected_cols)]
      }
      
      #convert dataframe into long form
      df<- long_form(df)
      colnames(df)<-c('Sub','Times','Conc')
      
    } else if (CI == T){
      
      ### we will be plotting confidence intervals
      
      vals <- df[,2:ncol(df)]
      confidence_intervals <- apply(vals,1,function (x) {
        c(lower = quantile(x, 0.05), upper = quantile(x, 0.95))}
      )
      confidence_intervals <- as.data.frame(t(confidence_intervals))
      df <- cbind(df$Times,rowMeans(vals),confidence_intervals)
      colnames(df)<-c('Times','Conc','5%','95%')
      
    }
    
    #Create Title for plot
    if (!is.na(compound_name)){
      header<-paste(compound_name, tissue_type,sep=' ')
    }else{
      header<-paste(casestudy_ID, tissue_type,'Concentration-Time Profiles',sep=' ')
    }
    
    if (units == 'ng/mL'){
      
      df <- df
      #mean_profile$Concentration <- mean_profile$Concentration
      ytitle <- paste('Concentration',paste0('(',units,')'), sep=' ')
      
    } else if (units == 'uM'){
      
      #Extract MW from curated data
      df_curated_data <- curated_data[which(curated_data$CODE == casestudy_ID),]
      MW <- df_curated_data$MW
      
      ## convert units of cmax from mg/L -> uM #check again
      df$Conc <- (df$Conc/1000)/as.numeric(MW)*1e6
      
      if (CI == T){
        df$`5%`<- (df$`5%`/1000)/as.numeric(MW)*1e6
        df$`95%`<- (df$`95%`/1000)/as.numeric(MW)*1e6
      }
      ytitle <- paste('Concentration',paste0('(',units,')'), sep=' ')
      
    } else{
      
      message('Please select either ng/mL or uM. ng/mL used as default')
      ytitle <- paste('Concentration (ng/mL)', sep=' ')
    }
    
    
    if (CI == F){
      g1 <- ggplot(df, aes(x = Times, y = Conc, group = factor(Sub))) + 
        geom_line(aes(color=factor(Sub)), size = 1.5, linetype = 1) +  
        xlab("Time (hours)") + ylab(ytitle) + ggtitle(header) + labs(color  = "Subjects", linetype = 1)
    } else if (CI == T){
      g1 <- ggplot(data = df) + 
        geom_line(aes(x = Times, y = Conc), size = 1.5, linetype = 1, colour = 'green3') + 
        geom_ribbon(aes(x = Times, ymin= `5%`, ymax= `95%`), alpha=0.2) +  
        xlab("Time (hours)") + ylab(ytitle) + ggtitle(header) 
    }
    
    
    if (logy == T){
      
      g1 <- g1 + general_theme +  scale_y_continuous(trans='log10')
      
      g1
      
    } else {
      
      #Plot the mean profile (y axis normal scale)
      g1 <- g1 +  general_theme
      
      g1
    }
    
  }
  
}

plot_parameters<-function(simcyp_outputs, Compound, 
                          plot_type = 'Distribution', x_variable, y_variable,
                          chosen_col, compound_name = NA){
  require(ggplot2)
  #Extract data relationg to compound of interest
  compound_IDs<- names(simcyp_outputs)
  index <- which(Compound==compound_IDs)
  outputs<- simcyp_outputs[[index]]
  
  #determine x and y variables
  x_col<- which(colnames(outputs)==x_variable)
  y_col <- which(colnames(outputs)== y_variable)
  
  x_axis_label <- paste(x_variable,determine_units(x_variable), sep=' ')
  
  if (!is.na(compound_name)){
    header_comp_name <- compound_name
  } else{
    header_comp_name <- Compound
  }
  
  if (plot_type == 'Distribution'){
    #Create Title for plot
    header<-paste(header_comp_name,'Distribution of', x_variable, sep=' ')
    
    #identify NA variables
    NA_idx <- which(is.na(outputs[,x_col]))
    if (length(NA_idx)>0){
      distribution_data <- outputs[-NA_idx,x_col]
    } else{
      distribution_data <- outputs[,x_col]
    }
    
    if (length(distribution_data) == 0){
      distribution_data <- rep(0, length(outputs[,x_col]))
    }
    
    # Distributions are plotted as density plots (only x variable is considered)
    ggplot(outputs, aes(x=distribution_data, colour= chosen_col, fill = chosen_col)) +
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
    header<-paste(header_comp_name, x_variable,'vs.', y_variable,'(r =',cor_coeff,')' ,sep=' ')
    
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
  
  df$parameter<-as.numeric(df$parameter)
  
  if(bar_order=='ascending'){
    df <- df[order(df$parameter),]
    
  } else if (bar_order == 'compound_code'){
    df <- df[order(df$Compounds),]
  }
  
  header<-paste('Mean',parameter,'Variation Across Simulated Compounds',sep=' ')
  
  units <- determine_units(parameter)
  ytitle <- paste(parameter,units,sep=' ')
  
  ggplot(df, aes(x=Compounds, y= parameter)) +
    ylab(parameter)+
    geom_bar(stat="identity", colour = bar_col, fill = bar_col)+
    ggtitle(header)+ ylab(ytitle)+
    scale_x_discrete(limits = df$Compounds) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line.x = element_line(color="black", size = 1),
      axis.line.y = element_line(color="black", size = 1),
      panel.grid.major = element_line(color = "grey"),
      panel.grid.minor = element_line(color = "grey"),
      axis.text = element_text(colour = "black",size=12),
      axis.text.x = element_text(angle = 90),
      axis.title = element_text(colour = "black",size=14,face="bold"))
}

PhyschemvsPredictedParamsPlot <- function(physchem_data,summary_simcyp,physchem_param,predicted_param,plot_colour){
  
  x_variable_idx <- which(colnames(physchem_data)==physchem_param)
  order_by_code_x <- physchem_data[order(physchem_data$Code),]
  x_variable <- order_by_code_x[,x_variable_idx]
  
  y_variable_idx <- which(colnames(summary_simcyp)== predicted_param)
  order_by_code_y <- summary_simcyp[order(summary_simcyp$CS_code),]
  y_variable <- as.numeric(order_by_code_y[,y_variable_idx])
  
  #collate all data
  df <- as.data.frame(cbind(order_by_code_x$Code,x_variable,y_variable))
  
  if(physchem_param=='Compound_type'|physchem_param=='HBD'){
    
    #create a tab to count observations in each bin
    data_count <- count(df,x_variable)
    
    #Create Title for plot
    header<-paste(predicted_param,'vs.', physchem_param,sep=' ')
    
    g1 <- ggplot(df, aes(x = x_variable, y = as.numeric(y_variable))) + 
      geom_boxplot(outlier.colour="red", outlier.shape=8,
                   outlier.size=2, fill=plot_colour, color="black")+
      geom_jitter(shape=16, position=position_jitter(0))
      #geom_text(data = data_count, aes(y = 0, label = n))
    g1 <- g1 + xlab(physchem_param) + ylab(paste(predicted_param,determine_units(predicted_param),sep = ' ')) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
    g1 + ggtitle(header)
    
  } else{
    
    #find the correlation coefficient
    cor_coeff<- round(cor(as.numeric(x_variable),y_variable),2)
    
    #Create Title for plot
    header<-paste(predicted_param,'vs.', physchem_param,'(r =',cor_coeff,')' ,sep=' ')
    
    #relationships are scatter diagrams where x and y variables are needed
    g1 <- ggplot(df, aes(x = as.numeric(x_variable), y = as.numeric(y_variable)))
    g1 <- g1 + geom_point(colour = plot_colour, size = 3) + geom_smooth(method = "lm", se = FALSE, colour = plot_colour)  
    g1 <- g1 + xlab(physchem_param) + ylab(paste(predicted_param,determine_units(predicted_param),sep = ' ')) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
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
