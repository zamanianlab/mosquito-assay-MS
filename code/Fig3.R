# read in libraries ------------------------------------------------------------

# working directory
library(here)

# data wrangling/plotting
library(tidyverse)
library(janitor)
library(imputeTS)

#library(dplyr)

# other plotting
library(cowplot)
library(ggtext)
library(ggbeeswarm)
library(ggpubr)
library(gghighlight)
library(scales)
library(cowplot)
library(patchwork)
library(ggh4x)
library(tibble)
library(drc)
library(jpeg)

# image read in 
library(magick)

# conflict resolution
library(conflicted)
conflict_prefer("filter", "dplyr")
conflicts_prefer(dplyr::select)
conflicts_prefer(cowplot::get_legend)


# plot theme ------ ------------------------------------------------------------

theme_zlab_white = function(base_size = 15, base_family = "Helvetica") {
  
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Specify axis options
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(size = base_size*0.4, color = "black", lineheight = 0.9, angle = 90),  
      axis.text.y = element_text(size = base_size*0.4, color = "black", lineheight = 0.9),  
      axis.ticks = element_blank(),  
      axis.title.x = element_text(size = base_size*0.6, color = "black", margin = margin(0, 10, 0, 0)),  
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

# Import data -------------------------------------------------------------------

motility_df <- readRDS(here('../data/RDS/motility_data.rds')) %>%
  filter(batch %in% c('insecticide')) %>%
  select(-site, -food)

loopbio_df <- readRDS(here('../data/RDS/loopbio_data.rds')) %>%
  filter(batch != 'resistant1') %>%
  rename(., LBflow = optical_flow) %>%
  select(-assay_type, -batch, -assay_date, -site)

merged <- left_join(motility_df, loopbio_df, by = c('species', 'other', 'preparation',
                                                    'treatment', 'conc', 'well',
                                                    'row', 'col', 'plate_ID'))


# temp <- loopbio_df %>%
#   select(plate_ID) %>%
#   unique(.)
# 
# temp2 <- merged %>%
#   select(plate_ID, species, other) %>%
#   unique(.)

# organize datasets ------------------------------------------------------------

# calculate r^2 value 
lm_fit <- lm(optical_flow ~ LBflow, data = merged)
r2_val <- summary(lm_fit)$r.squared

(comparison <- merged %>%
  #filter(assay_date == '20251201') %>%
  ggplot(aes(x = optical_flow, y = LBflow, color = other)) +
  geom_smooth(method = "lm", se = FALSE, color = "gray") +
  geom_point() +
  annotate("text",
            x = Inf, y = Inf,
            label = paste0("R² = ", round(r2_val, 3)),
            hjust = 3, vjust = 1.5,
            size = 4) +
   xlab("IX Plate") +
   ylab("Loopbio Plate") +
   theme_zlab_white() +
   theme(axis.line = element_line(size = 0.5),
         axis.text.x = element_text(angle = 0),
         legend.text = element_text(size = 12))
   )

# random data plot 

(sample <- merged %>%
    filter(assay_date == '20251206',
           species == 'AlboMO') %>%
    mutate(conc = factor(conc, levels = c('1p', '0uM', '0.001ppm', '0.01ppm', '0.1ppm', '1ppm', 
                                          '10ppm', '100ppm', '0.01mg/ml', '0.1mg/ml',
                                          '1mg/ml', '10mg/ml', '100mg/ml', '1000mg/ml',
                                          '500pM', '1nM',
                                          '5nM', '10nM', '50nM', '100nM', '500nM', '1uM',
                                          '5uM', '10uM', '50uM', '100uM'),
                         ordered = TRUE)) %>%
    ggplot(aes(x = treatment, y = optical_flow, color = conc)) +
    geom_boxplot(show.legend = FALSE) +
    geom_quasirandom(fill = 'white', shape = 21, size = 1,
                     alpha = 0.75, dodge.width = 0.7, width = 0.1, show.legend = FALSE
    ) +
    facet_grid(cols = vars(preparation), 
               scales = 'free', space = 'free_x') +
    theme_zlab_white() +
    theme(axis.line = element_line(color = "black", linewidth = 0.5),
          panel.grid.minor = element_line(color = 'black', linewidth = 0.2),
          axis.text.x = element_text(angle = 0)) +
    ylab("IX Motility") + 
    xlab(""))

(sampleLB <- merged %>%
    filter(assay_date == '20251206',
           species == 'AlboMO') %>%
    mutate(conc = factor(conc, levels = c('1p', '0uM', '0.001ppm', '0.01ppm', '0.1ppm', '1ppm', 
                                          '10ppm', '100ppm', '0.01mg/ml', '0.1mg/ml',
                                          '1mg/ml', '10mg/ml', '100mg/ml', '1000mg/ml',
                                          '500pM', '1nM',
                                          '5nM', '10nM', '50nM', '100nM', '500nM', '1uM',
                                          '5uM', '10uM', '50uM', '100uM'),
                         ordered = TRUE)) %>%
    ggplot(aes(x = treatment, y = LBflow, color = conc)) +
    geom_boxplot(show.legend = FALSE) +
    geom_quasirandom(fill = 'white', shape = 21, size = 1,
                     alpha = 0.75, dodge.width = 0.7, width = 0.1, show.legend = FALSE
    ) +
    facet_grid(cols = vars(preparation), 
               scales = 'free', space = 'free_x') +
    theme_zlab_white() +
    theme(axis.line = element_line(color = "black", linewidth = 0.5),
          panel.grid.minor = element_line(color = 'black', linewidth = 0.2),
          axis.text.x = element_text(angle = 0)) +
    ylab("LB Motility") + 
    xlab(""))

(comparison <- plot_grid(sample, sampleLB, ncol = 1))

# check timing situations ----------

P1_30 <- read.csv(here("../data/loopbio/20251206-p07-KTR_tidy.csv")) %>%
  mutate(time = '30sec')
P2_30 <- read.csv(here("../data/loopbio/20251206-p08-KTR_tidy.csv")) %>%
  mutate(time = '30sec')
P1_2  <- read.csv(here("../data/loopbio/20251206-p11-KTR_tidy.csv")) %>%
  mutate(time = '2min')
P1_5  <- read.csv(here("../data/loopbio/20251206-p12-KTR_tidy.csv")) %>%
  mutate(time = '5min')
P2_2  <- read.csv(here("../data/loopbio/20251206-p13-KTR_tidy.csv")) %>%
  mutate(time = '2min')
P2_5  <- read.csv(here("../data/loopbio/20251206-p14-KTR_tidy.csv")) %>%
  mutate(time = '5min')

all <- rbind(P1_30, P2_30, P1_2, P2_2, P1_5, P2_5) %>%
  mutate(plate_ID = case_when(
         plate == '20251206-p07-KTR' ~ '20251206-p07-KTR_35958',
         plate == '20251206-p11-KTR' ~ '20251206-p07-KTR_35958',
         plate == '20251206-p12-KTR' ~ '20251206-p07-KTR_35958',
         plate == '20251206-p08-KTR' ~ '20251206-p08-KTR_35959',
         plate == '20251206-p13-KTR' ~ '20251206-p08-KTR_35959',
         plate == '20251206-p14-KTR' ~ '20251206-p08-KTR_35959')) 

adjusted <- all %>%
  select(-plate) %>%
  pivot_wider(names_from = time, values_from = optical_flow)  %>%
  mutate(col = as.character(col))

IX_data <- readRDS(here('../data/RDS/motility_data.rds')) %>%
  filter(plate_ID %in% c('20251206-p07-KTR_35958', '20251206-p08-KTR_35959')) %>%
  mutate(col = as.character(col))


time_merge <- left_join(adjusted, IX_data, by = c('plate_ID', 'well', 'other', 
                                                  'conc', 'treatment'))

# compare 30 sec to IX ---- 

# calculate r^2 value 
lm_fit <- lm(optical_flow ~ `30sec`, data = time_merge)
r2_val <- summary(lm_fit)$r.squared

(A <- time_merge %>%
    #filter(assay_date == '20251201') %>%
    ggplot(aes(x = optical_flow, y = `30sec`)) +
    geom_smooth(method = "lm", se = FALSE, color = "gray") +
    geom_point(size = 1) +
    annotate("text",
             x = Inf, y = Inf,
             label = paste0("R² = ", round(r2_val, 3)),
             hjust = 2, vjust = 2,
             size = 3) +
    xlab("IX Plate") +
    ylab("Loopbio 30 sec") +
    theme_zlab_white() +
    theme(axis.line = element_line(size = 0.5),
          axis.text.x = element_text(angle = 0, size = 4),
          axis.text.y = element_text(size = 4),
          legend.text = element_text(size = 8))
)

# compare 2 min to IX ---- 

# calculate r^2 value 
lm_fit <- lm(optical_flow ~ `2min`, data = time_merge)
r2_val <- summary(lm_fit)$r.squared

(B <- time_merge %>%
    #filter(assay_date == '20251201') %>%
    ggplot(aes(x = optical_flow, y = `2min`)) +
    geom_smooth(method = "lm", se = FALSE, color = "gray") +
    geom_point(size = 1) +
    annotate("text",
             x = Inf, y = Inf,
             label = paste0("R² = ", round(r2_val, 3)),
             hjust = 2, vjust = 2,
             size = 3) +
    xlab("IX Plate") +
    ylab("Loopbio 2 min") +
    theme_zlab_white() +
    theme(axis.line = element_line(size = 0.5),
          axis.text.x = element_text(angle = 0, size = 4),
          axis.text.y = element_text(size = 4),
          legend.text = element_text(size = 8))
)

# compare 2 min to IX ---- 

# calculate r^2 value 
lm_fit <- lm(optical_flow ~ `5min`, data = time_merge)
r2_val <- summary(lm_fit)$r.squared

(C <- time_merge %>%
    #filter(assay_date == '20251201') %>%
    ggplot(aes(x = optical_flow, y = `5min`)) +
    geom_smooth(method = "lm", se = FALSE, color = "gray") +
    geom_point(size = 1) +
    annotate("text",
             x = Inf, y = Inf,
             label = paste0("R² = ", round(r2_val, 3)),
             hjust = 2, vjust = 2,
             size = 3) +
    xlab("IX Plate") +
    ylab("Loopbio 5 min") +
    theme_zlab_white() +
    theme(axis.line = element_line(size = 0.5),
          axis.text.x = element_text(angle = 0, size = 4),
          axis.text.y = element_text(size = 4),
          legend.text = element_text(size = 8))
)

# make Figure 3 ----- 

(comparison <- plot_grid(A, B, C, nrow = 1, labels = c('B', '', '')))

(schematic <- image_read(here("../images/mos-Fig3.jpeg")))

(schematic_f <- ggdraw() + draw_image(schematic))

(figure3 <- plot_grid(schematic_f, comparison, nrow = 2, labels = c("A", ""),
                      align = 'hv'))

ggsave(here('../Figure3/Fig3.pdf'), figure3, width = 7, height = 4, units = 'in')

# plot w/ time differences 

flipped <- time_merge %>%
  pivot_longer(cols = c('optical_flow', '30sec', '2min', '5min'), names_to = 'time', values_to = 'flow')

(sample <- flipped %>%
    mutate(conc = factor(conc, levels = c('1p', '0uM', '0.001ppm', '0.01ppm', '0.1ppm', '1ppm', 
                                          '10ppm', '100ppm', '0.01mg/ml', '0.1mg/ml',
                                          '1mg/ml', '10mg/ml', '100mg/ml', '1000mg/ml',
                                          '500pM', '1nM',
                                          '5nM', '10nM', '50nM', '100nM', '500nM', '1uM',
                                          '5uM', '10uM', '50uM', '100uM'),
                         ordered = TRUE),
           time = factor(time, levels = c('optical_flow', '30sec', '2min', '5min'),
                         labels = c('IX', 'LB 30s', 'LB 2min', 'LB 5min'),
                         ordered = TRUE)) %>%
    filter(treatment %in% c('Chlorfenapyr', 'Spinosad')) %>%
    ggplot(aes(x = treatment, y = flow, color = conc)) +
    geom_boxplot(show.legend = FALSE) +
    geom_quasirandom(fill = 'white', shape = 21, size = 1,
                     alpha = 0.75, dodge.width = 0.7, width = 0.1, show.legend = FALSE
    ) +
    facet_grid(rows = vars(time), 
               scales = 'free', space = 'free_x') +
    theme_zlab_white() +
    theme(axis.line = element_line(color = "black", linewidth = 0.5),
          panel.grid.minor = element_line(color = 'black', linewidth = 0.2),
          axis.text.x = element_text(angle = 0)) +
    ylab("Motility") + 
    xlab(""))


# save Figure 

ggsave(here('../Figure1/Fig1C.pdf'), fiji, width = 4, height = 2.5, units = 'in')

ggsave(here('../Figure1/Fig1C-alternative.pdf'), model, width = 4, height = 2.5, units = 'in')
