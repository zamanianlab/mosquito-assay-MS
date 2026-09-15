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

plate <- read.csv(here("../data/20250519-p08-KTR_tidy.csv"))

# plot theme ------ ------------------------------------------------------------

theme_zlab_white = function(base_size = 28, base_family = "Helvetica") {
  
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Specify axis options
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(size = base_size*0.4, color = "black", lineheight = 0.9, angle = 90),  
      axis.text.y = element_text(size = base_size*0.4, color = "black", lineheight = 0.9),  
      axis.ticks = element_blank(),  
      axis.title.x = element_text(size = base_size*0.6, color = "black", margin = margin(10, 10, 0, 0)),  
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

# calculate r^2 value 
lm_fit <- lm(size ~ manual, data = plate)
r2_val <- summary(lm_fit)$r.squared

(fiji <- ggplot(data = plate, aes(x = manual, y = size, color)) +
  geom_smooth(method = "lm", se = FALSE, color = "gray") +
  geom_point() +
  #geom_text(aes(label = well), vjust = -0.7, size = 2) +
  annotate("text",
            x = Inf, y = Inf,
            label = paste0("R² = ", round(r2_val, 3)),
            hjust = 3, vjust = 1.5,
            size = 4) +
   xlab("Manual Measurement") +
   ylab("Model Object Size") +
   theme_zlab_white() +
   theme(axis.line = element_line(size = 0.5),
         axis.text.x = element_text(angle = 0))
   )

# alternative 
lm_fit <- lm(size ~ length_px, data = plate)
r2_val <- summary(lm_fit)$r.squared

(model <- ggplot(data = plate, aes(x = length_px, y = size, color)) +
    geom_smooth(method = "lm", se = FALSE, color = "gray") +
    geom_point() +
    geom_text(aes(label = well), vjust = -0.7, size = 2) +
    # annotate("text",
    #          x = Inf, y = Inf,
    #          label = paste0("R² = ", round(r2_val, 3)),
    #          hjust = 3, vjust = 1.5,
    #          size = 4) +
    xlab("Length") +
    ylab("Area") +
    theme_zlab_white() +
    theme(axis.line = element_line(size = 0.5),
          axis.text.x = element_text(angle = 0))
)

# save Figure 

ggsave(here('../Figure1/Fig1C.pdf'), fiji, width = 4, height = 2.5, units = 'in')

ggsave(here('../Figure1/Fig1C-alternative.pdf'), model, width = 4, height = 2.5, units = 'in')
