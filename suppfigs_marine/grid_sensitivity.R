## Assesses how the number of grid cells used to calculate delta fish biomass
## for each community affects the number of FishMIP models that agree on
## the direction of fish biomass change

rm(list=ls())
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)

world_data <- map_data("world")


times <- c("2046-2056", "2069-2099")
ssps <- c("ssp126", "ssp585")
num_cells <- c("20")
plots <- list()

store_mat <- matrix(0, nrow = length(times)*length(ssps), ncol = length(num_cells))
store_mat2 <- matrix(0, nrow = 72, ncol = length(times)*length(ssps))

colnames(store_mat) <- num_cells


for(i in 1:length(times)){
  curr_time <- times[i]
  
  for(j in 1:length(ssps)){
    curr_ssp <- ssps[j]  
    store_mat_large <- matrix(0, nrow = 72, ncol = length(num_cells))
    
    for(k in 1:length(num_cells)){
      
      curr_file <- paste("./fishmip_output_raw/fishmip_tcb_response_71_sites_", curr_time, "_", curr_ssp, "_nearest", num_cells[k], "_raw.csv", sep = "")
      curr_dat <- read.csv(curr_file)
      curr_models <- curr_dat[,-c(1:6)]
      curr_models[curr_models > 0] <- TRUE
      curr_models[curr_models <= 0] <- FALSE
        
      kk <- rowSums(curr_models)
      kk[kk < 8] <- 16 - kk[kk < 8] # Correct agreement number for when more models say negative than positive
    
      store_mat2[,(i-1)*j+i] <- kk 
        
      store_mat_large[,k] <- kk
      store_mat[(i-1)*length(ssps)+j,k] <- mean(kk) # Number of models that say positive change in this region
    }
    
    ## MAPS OF NUMBER OF MODELS IN AGREEMENT
    colnames(store_mat_large) <- num_cells
    gg_dat <- as.data.frame(store_mat_large)
    gg_dat$Lat <- curr_dat$latDD
    gg_dat$Lon <- curr_dat$lonDD
   
    cell_plots <- c("1","3","5","10", "20", "50","100")
    # 
    theme_opts1 <- list(theme(panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.background = element_blank(),
                              plot.background = element_rect(fill="white"),
                              legend.position = "bottom",
                              axis.line = element_blank(),
                              axis.text.x =  element_text(size=16, face = "bold"),
                              axis.text.y = element_text(size=16, face = "bold"),
                              axis.ticks = element_blank(),
                              axis.title.x = element_text(size=16, face = "bold"),
                              axis.title.y = element_text(size=16, face = "bold"),
                              plot.title = element_text(size=18, hjust = 0.5, face = "bold")))
    
    plot_list <- list()
    plot_name <- paste(curr_time, curr_ssp, sep = " ")
    for(m in 1:length(cell_plots)){
      curr_gg_dat <- gg_dat[,c("Lon","Lat",cell_plots[m])]
      colnames(curr_gg_dat) <- c("Lon","Lat","value")
    
      plot_list[[(m-1)*2+1]] <-  ggplot(curr_gg_dat) + 
        geom_point(aes(x=Lon, y=Lat, colour = value, size = 3)) +   geom_polygon(data=world_data, mapping=aes(x=long, y=lat, group=group), alpha = 0, colour ='black')+
        coord_cartesian(ylim = c(-30, 20), xlim = c(25,165)) + theme_opts1 + scale_colour_gradient2(low = 'blue', mid = "green", high = 'red', limits = c(8, 16),
                                                                                                    midpoint = 12, name = "No. Models\n in Agreement ")+
        labs(title = paste("No. cells = ", cell_plots[m], sep = ""))
      
      plot_list[[(m-1)*2+2]] <- ggplot(curr_gg_dat) + 
        geom_histogram(aes(value), bins = 8, binwidth = 1,colour = "black", fill = "royalblue") + coord_cartesian(ylim = c(0,30), xlim = c(8,16)) + 
        xlab("No. Models in Agreement") + ylab("No. Communities") + labs(title = paste("No. cells = ", cell_plots[m], sep = "")) +theme_opts1 +
        annotate("text", x = 10, y = 25, label = paste("Mean = ", round(mean(curr_gg_dat$value),1), "\n SD = ", round(sd(curr_gg_dat$value),1)), size = 6, fontface = "bold")
    }
    
    figure <- ggpubr::ggarrange(plotlist = plot_list, ncol = 2, nrow = length(cell_plots), common.legend = TRUE, align = c("hv"))
    ggpubr::ggexport(figure, filename = paste("./cinner_figures/", plot_name,".png", sep = ""), width = 900, height = 1600)
    
    
  }
  
}
