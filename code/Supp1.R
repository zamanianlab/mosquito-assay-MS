# data wrangling/plotting
library(tidyverse)
library(janitor)

# other plotting
library(reprex)
library(cowplot)
library(ggtext)
library(ggbeeswarm)
library(wesanderson)
library(here)

# misc
library(conflicted)

# conflict resolution
conflict_prefer("filter", "dplyr")

# Import data -------------------------------------------------------------------

motility <- readRDS(here("../data/RDS/motility_data.rds"))

# plot theme ------ ------------------------------------------------------------

theme_zlab_white = function(base_size = 20, base_family = "Helvetica") {
  
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Specify axis options
      axis.line = element_blank(),
      #axis.line = element_line(color = "black"),
      axis.text.x = element_text(size = base_size*0.4, color = "black", lineheight = 0.9, angle = 90),  
      axis.text.y = element_text(size = base_size*0.4, color = "black", lineheight = 0.9),  
      axis.ticks = element_blank(),  
      axis.title.x = element_text(size = base_size*0.6, color = "black", margin = margin(10, 20, 0, 0)),  
      axis.title.y = element_text(size = base_size*0.6, color = "black", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "white"),  
      legend.key = element_rect(color = "white",  fill = "white"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.2, color = "black"),  
      legend.title = element_text(size = base_size*0.2, face = "bold", hjust = 0, color = "black"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "white", color  =  NA),  
      panel.border = element_blank(),  
      #panel.border = element_blank(),
      #panel.grid.major = element_line(color = "black"),  
      panel.grid.major = element_blank(), 
      #panel.grid.minor = element_line(color = "black"),  
      panel.grid.minor = element_blank(),  
      panel.spacing = unit(0.5, "lines"),   
      # Specify faceting options
      strip.background = element_rect(fill = "white", color = "white"),  
      strip.text.x = element_text(size = base_size*0.4, color = "black"),  
      strip.text.y = element_text(size = base_size*0.4, color = "black",angle = -90),  
      # Specify plot options
      plot.background = element_rect(colour = NA, fill = NA),  
      plot.title = element_text(size = base_size*1.2, color = "black"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
}


# organize datasets ------------------------------------------------------------

# temp <- cxpimaging %>%
#   select(acquisition) %>%
#   unique(.)

cxpimaging <- motility %>%
  filter(plate_ID %in% c('20250603-p01-KTR_35104', '20250603-p02-KTR_35105',
                         '20250603-p03-KTR_35106', '20250603-p04-KTR_35107',
                         '20250604-p01-KTR_35113', '20250604-p02-KTR_35114',
                         '20250604-p03-KTR_35115', '20250604-p04-KTR_35116',
                         '20250605-p01-KTR_35130', '20250605-p02-KTR_35132',
                         '20250605-p03-KTR_35133', '20250605-p04-KTR_35134',
                         '20250606-p05-KTR_35151', '20250606-p06-KTR_35152')) %>%
  #filter(group == 'cxp_imaging') %>%
  separate(preparation, into = c('preparation', 'acquisition'), sep = '_') %>%
  filter(acquisition %in% c('+GFP', '+TexasRed', 'DAPI10', '20', 'Tlwell', 'Tlplate', 
                            'Tlopenclose'))

(plot <- cxpimaging %>%
  filter(treatment == 'DMSO') %>%
  mutate(
    acquisition = factor(acquisition, 
                         levels = c('Tlwell', 'Tlopenclose', '20', '+GFP', 
                                    '+TexasRed', 'DAPI10', 'Tlplate'), 
                         labels = c('Protocol', 'Shutter Flutter', '20 Frames', '+ GFP WL', 
                                    '+ Texas Red WL', '+ DAPI WL', 'Cross Plate'),
                         ordered = TRUE)
  ) %>%
  ggplot(aes(x = acquisition, y = optical_flow)) +
  geom_boxplot() +
  geom_quasirandom(
    aes(y = optical_flow, group = conc), fill = 'white', shape = 21, size = 1,
    alpha = 0.75, dodge.width = 0.7, width = 0.1, show.legend = TRUE
  ) +
  theme_zlab_white() +
  theme(axis.line = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor = element_line(color = 'black', linewidth = 0.2),
        axis.text.x = element_text(angle = 0)) +
  ylab("Motility") +
  xlab("Imaging Approach"))

# save Figure 

ggsave(here('../Supp1/S1C.pdf'), plot, width = 7.5, height = 2.5, units = 'in')

