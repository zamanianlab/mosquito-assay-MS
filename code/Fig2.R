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
library(ggh4x)
library(ggpubr)
library(patchwork)

# image read in 
library(magick)

# misc
library(conflicted)
library(drc)

# conflict resolution
conflict_prefer("filter", "dplyr")
conflicts_prefer(dplyr::select)

# plot theme ------ ------------------------------------------------------------

theme_zlab_white = function(base_size = 14, base_family = "Helvetica") {
  
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
      strip.text.x = element_text(size = base_size*0.6, color = "black"),  
      strip.text.y = element_text(size = base_size*0.6, color = "black",angle = -90),  
      # Specify plot options
      plot.background = element_rect(colour = NA, fill = NA),  
      plot.title = element_text(size = base_size*1.2, color = "black"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
}

# Import data and correct metadata----------------------------------------------------

motility_raw <- readRDS(here('../data/RDS/motility_data.rds')) %>%
  filter(assay_date %in% c('20251201', '20251202', '20251203', 
                           '20251204', '20251205', '20251206',
                           '20251207', '20251208'),
         batch == 'insecticide',
         !(treatment == 'Altosid' & conc == '0uM')) %>%
  mutate(solvent = case_when(
    treatment == 'Water' ~ 'Water',
    treatment == 'DMSO' ~ 'DMSO',
    treatment == 'Vectolex' ~ 'Water',
    treatment == 'Spinosad' ~ 'Water',
    treatment == 'Altosid' ~ 'Water',
    treatment == 'Chlorfenapyr' ~ 'DMSO',
    treatment == 'IVM' ~ 'DMSO'
  ),
  other = case_when(
    plate_ID == '20251205-p04-KTR_35939' ~ '48hr',
    plate_ID == '20251204-p06-KTR_35928' ~ '24hr',
    plate_ID == '20251205-p06-KTR_35941' ~ '48hr',
    TRUE ~ other
  ),
  conc = case_when(
    conc == '0.01mg/ml' ~ '0.01mg/L',
    conc == '0.1mg/ml' ~ '0.1mg/L',
    conc == '100mg/ml' ~ '100mg/L',
    conc == '10mg/ml' ~ '10mg/L',
    conc == '1mg/ml' ~ '1mg/L',
    TRUE ~ conc
  ))

development_raw <- readRDS(here('../data/RDS/development_data.rds')) %>%
  filter(assay_date %in% c('20251201', '20251202', '20251203', 
                           '20251204', '20251205', '20251206',
                           '20251207', '20251208'),
         batch == 'insecticide',
         !is.na(object_number)) %>%
  mutate(solvent = case_when(
    treatment == 'Water' ~ 'Water',
    treatment == 'DMSO' ~ 'DMSO',
    treatment == 'Vectolex' ~ 'Water',
    treatment == 'Spinosad' ~ 'Water',
    treatment == 'Altosid' ~ 'Water',
    treatment == 'Chlorfenapyr' ~ 'DMSO',
    treatment == 'IVM' ~ 'DMSO'
  ))

# normalize data ------------------------------------------------------------

# normalize to negative controls 
NegCon <- motility_raw %>%
  filter(treatment %in% c('DMSO', 'Water')) %>%
  group_by(assay_date, species, preparation, other, solvent) %>%
  mutate(negCon = mean(optical_flow)) %>%
  select(negCon, assay_date, species, preparation, other, solvent) %>%
  unique()

# find combined positive control for each assay date 
PosCon <- motility_raw %>%
  filter(treatment == 'IVM') %>%
  group_by(assay_date, species, preparation, other) %>%
  mutate(posCon = mean(optical_flow)) %>%
  select(posCon, assay_date, species, preparation, other) %>%
  unique()


# normalize to negative and positive controls for each assay date 
normalized <- motility_raw %>%
  left_join(., NegCon, by = c('assay_date', 'species', 'preparation', 'other', 'solvent')) %>%
  left_join(., PosCon, by = c('assay_date', 'species', 'preparation', 'other')) %>%
  mutate(norm_mot = (optical_flow - posCon) / (negCon - posCon)) %>%
  mutate(species = factor(species, 
                          levels = c('AeLVP', 'AlboMO', 'Cxp'),
                          labels = c('*Ae. aegypti* LVP', '*Ae. albopictus* MO', '*Cx. pipiens*')))

norm_means <- normalized %>%
  group_by(species, conc, treatment, assay_date, preparation, other) %>%
  summarise(mean_mot = median(norm_mot, na.rm = TRUE)) %>%
  ungroup()

normalized_mot <- normalized %>%
  left_join(., norm_means) %>%
  filter(treatment != 'missed')

# prune outliers 

normalized_mot <- normalized_mot %>%
  # adding assay_date here doesn't change anything 
  group_by(preparation, treatment, conc, other, species) %>%
  mutate(IQR = IQR(norm_mot),
         O_upper = quantile(norm_mot, probs=c( .75), na.rm = FALSE)+1.5*IQR,  
         O_lower = quantile(norm_mot, probs=c( .25), na.rm = FALSE)-1.5*IQR  
  ) %>% 
  filter(O_lower <= norm_mot & norm_mot <= O_upper) %>%
  ungroup()

# motility plot (Fig2B) -------------------------------------------------------

controls5per <- normalized_mot %>%
  filter(preparation == '5/well',
         treatment %in% c('DMSO', 'Water', 'IVM')) %>%
  mutate(treatment = case_when(
    treatment == 'DMSO' ~ 'Negative Control',
    treatment == 'IVM' ~ 'Positive Control',
    treatment == 'Water' ~ 'Negative Control',
    TRUE ~ treatment
  ))

(Fig2B_neg <- controls5per %>%
    filter(treatment == 'Negative Control',
           other != '1hr') %>%
    mutate(other = factor(other, levels = c('24hr', '48hr', '72hr'),
                          ordered = TRUE)) %>%
    ggplot(aes(x = treatment, y = norm_mot, color = species)) +
    geom_hline(yintercept = 1, linetype="dashed", color = "black", size=0.5) +
    geom_hline(yintercept = 0, linetype="dashed", color = "black", size=0.5) +
    geom_point(aes(colour = species), size = 0.9, alpha = 0.3, position=position_dodge(width=0.3)) +
    #geom_jitter(size = 0.75, alpha = 0.75) +
    facet_grid(cols = vars(treatment), rows = vars(other), scales = 'free') +
    facetted_pos_scales(
      y = list(
        other == "1hr" ~ scale_y_continuous(limits = c(-1, 3)),
        other != "1hr" ~ scale_y_continuous(limits = c(-0.5, 2)))) +
    scale_color_manual(name = "Species",
      values = c('*Ae. aegypti* LVP' = 'chartreuse3', '*Ae. albopictus* MO' = 'goldenrod3', '*Cx. pipiens*' = 'darkorchid')) +
    theme_zlab_white() +
    theme(axis.line = element_line(color = "black", linewidth = 0.5),
          legend.position = "none",
          axis.text.x = element_blank(),
          strip.text.y = element_blank()) +
    ylab("Motility") +
    xlab(""))

normalized_mot <- normalized_mot %>%
  mutate(
    treatment_lab = case_when(
      treatment == "Altosid" ~ "Altosid (mg/L)",
      treatment == "Chlorfenapyr" ~ "Chlorfenapyr (µM)",
      treatment == "Spinosad" ~ "Spinosad (ppm)",
      treatment == "Vectolex" ~ "Vectolex (ppm)",
      TRUE ~ treatment
    )
  )


# stats data prep 

# temp <- normalized_mot %>%
#   select(treatment) %>%
#   unique(.)

normalized_mot <- normalized_mot %>%
  filter(preparation == '5/well')

altosid <- normalized_mot %>%
  filter(treatment %in% c("Water", "Altosid"))
chlorfenapyr <- normalized_mot %>%
  filter(treatment %in% c("DMSO", "Chlorfenapyr"))
spinosad <- normalized_mot %>%
  filter(treatment %in% c("Water", "Spinosad"))
vectolex <- normalized_mot %>%
  filter(treatment %in% c("Water", "Vectolex"))

# stats calculations 

stats.altosid <- compare_means(norm_mot~conc, data = altosid, group.by= c("species", "other"),
                                 method = "t.test", ref.group = "0uM") %>%
  mutate(p.signif = factor(p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***")),
         treatment_lab = 'Altosid (mg/L)')
stats.chlorfenapyr <- compare_means(norm_mot~conc, data = chlorfenapyr, group.by= c("species", "other"),
                               method = "t.test", ref.group = "1p") %>%
  mutate(p.signif = factor(p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***")),
         treatment_lab = 'Chlorfenapyr (µM)')
stats.spinosad <- compare_means(norm_mot~conc, data = spinosad, group.by= c("species", "other"),
                               method = "t.test", ref.group = "0uM") %>%
  mutate(p.signif = factor(p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***")),
         treatment_lab = 'Spinosad (ppm)')
stats.vectolex <- compare_means(norm_mot~conc, data = vectolex, group.by= c("species", "other"),
                               method = "t.test", ref.group = "0uM") %>%
  mutate(p.signif = factor(p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***")),
         treatment_lab = 'Vectolex (ppm)')

# stats merge 
# final stats: order conc, other; set y values 

stats <- rbind(stats.altosid, stats.chlorfenapyr, stats.spinosad, stats.vectolex) %>%
  mutate(conc = group2,
         y = case_when(
           other %in% c('24hr', '48hr', '72hr') & species == '*Cx. pipiens*' ~ 1.7,
           other %in% c('24hr', '48hr', '72hr') & species == '*Ae. albopictus* MO' ~ 1.8,
           other %in% c('24hr', '48hr', '72hr') & species == '*Ae. aegypti* LVP' ~ 1.9,
           other == '1hr' & species == '*Cx. pipiens*' ~ 2.4,
           other == '1hr' & species == '*Ae. albopictus* MO' ~ 2.6,
           other == '1hr' & species == '*Ae. aegypti* LVP' ~ 2.8
           )) %>%
  mutate(conc = factor(conc, levels = c('1p', '0uM', '0.001ppm', '0.01ppm', '0.1ppm', '1ppm',
                                        '10ppm', '100ppm',
                                        '0.01mg/L', '0.1mg/L',
                                        '1mg/L', '10mg/L', '100mg/L', '1000mg/L',
                                        '500pM', '1nM',
                                        '5nM', '10nM', '50nM', '100nM', '500nM', '1uM',
                                        '5uM', '10uM', '50uM', '100uM'),
                       labels = c('1%', '0µM', '0.001', '0.01', '0.1', '1',
                                  '10', '100', 
                                  '0.01', '0.1',
                                  '1', '10', '100', '1000',
                                  '0.0005', '0.001',
                                  '0.005', '0.01', '0.05', '0.1', '0.5', '1',
                                  '5', '10', '50', '100'),
                       ordered = TRUE),
         other = factor(other, levels = c('1hr', '24hr', '48hr', '72hr'),
                        ordered = TRUE)) %>%
  filter(other != '1hr')

(Fig2B_main <- normalized_mot %>%
   filter(!treatment %in% c('DMSO', 'Water', 'IVM'),
          other != '1hr') %>%
   mutate(conc = factor(conc, levels = c('1p', '0uM', '0.001ppm', '0.01ppm', '0.1ppm', '1ppm',
                                         '10ppm', '100ppm',
                                         '0.01mg/L', '0.1mg/L',
                                         '1mg/L', '10mg/L', '100mg/L', '1000mg/L',
                                         '500pM', '1nM',
                                         '5nM', '10nM', '50nM', '100nM', '500nM', '1uM',
                                         '5uM', '10uM', '50uM', '100uM'),
                        labels = c('1%', '0µM', '0.001', '0.01', '0.1', '1',
                                   '10', '100', 
                                   '0.01', '0.1',
                                   '1', '10', '100', '1000',
                                   '0.0005', '0.001',
                                   '0.005', '0.01', '0.05', '0.1', '0.5', '1',
                                   '5', '10', '50', '100'),
                        ordered = TRUE),
          other = factor(other, levels = c('24hr', '48hr', '72hr'),
                         ordered = TRUE)) %>%
   ggplot(aes(x = conc, y = norm_mot, color = species)) +
   geom_hline(yintercept = 1, linetype="dashed", color = "black", size=0.5) +
   geom_hline(yintercept = 0, linetype="dashed", color = "black", size=0.5) +
   geom_point(aes(colour = species), size = 0.75, alpha = 0.5, position=position_dodge(width=0.4)) +
   geom_line(aes(group = species, colour = species), size = 0.1, alpha = 0.05) +
   geom_smooth(aes(group = species, colour = species), size = 0.75, alpha = 1, se=F) +
   geom_text(
      data = stats, aes(label = p.signif, y=y, color = species),
      size = 3, show.legend = FALSE) +
   facet_grid(cols = vars(treatment_lab), rows = vars(other), scales = 'free') +
   scale_color_manual(name = "Species",
     values = c('*Ae. aegypti* LVP' = 'chartreuse3', '*Ae. albopictus* MO' = 'goldenrod3', '*Cx. pipiens*' = 'darkorchid')) +
   facetted_pos_scales(
     y = list(
       other == "1hr" ~ scale_y_continuous(limits = c(-1, 3)),
       other != "1hr" ~ scale_y_continuous(limits = c(-0.5, 2)))) +
   theme_zlab_white() +
   theme(axis.line = element_line(color = "black", linewidth = 0.5),
         axis.text.x = element_text(angle = 0),
         legend.position = 'none',
         #legend.direction = 'horizontal',
         strip.text.y = element_blank()) +
   ylab("") +
   xlab("Concentration"))

(Fig2B_pos <- controls5per %>%
    filter(treatment == 'Positive Control',
           other != '1hr') %>%
    mutate(other = factor(other, levels = c('24hr', '48hr', '72hr'),
                          ordered = TRUE)) %>%
    ggplot(aes(x = treatment, y = norm_mot, color = species)) +
    geom_hline(yintercept = 1, linetype="dashed", color = "black", size=0.5) +
    geom_hline(yintercept = 0, linetype="dashed", color = "black", size=0.5) +
    geom_point(aes(colour = species), size = 0.75, alpha = 0.5, position=position_dodge(width=0.3)) +
    facet_grid(cols = vars(treatment), rows = vars(other), scales = 'free') +
    facetted_pos_scales(
      y = list(
        other == "1hr" ~ scale_y_continuous(limits = c(-1, 3)),
        other != "1hr" ~ scale_y_continuous(limits = c(-0.5, 2)))) +
    scale_color_manual(name = "Species",
      values = c('*Ae. aegypti* LVP' = 'chartreuse3', '*Ae. albopictus* MO' = 'goldenrod3', '*Cx. pipiens*' = 'darkorchid')) +
    theme_zlab_white() +
    theme(axis.line = element_line(color = "black", linewidth = 0.5),
          legend.position = "none",
          axis.text.x = element_blank()) +
    ylab("") +
    xlab(""))

(motility_final <- plot_grid(NULL, Fig2B_neg, NULL, Fig2B_main, NULL, Fig2B_pos, nrow = 1,
                             rel_widths = c(-0.05, 0.8, -0.15, 2.2, -0.15, 0.85))) 



# development plot --------------------


# normalize to negative controls 
NegCon <- development_raw %>%
  filter(treatment %in% c('DMSO', 'Water')) %>%
  group_by(assay_date, species, preparation, other, solvent) %>%
  mutate(negCon = mean(size)) %>%
  select(negCon, assay_date, species, preparation, other, solvent) %>%
  unique()

# find combined positive control for each assay date 
PosCon <- development_raw %>%
  filter(treatment == 'IVM') %>%
  group_by(assay_date, species, preparation, other) %>%
  mutate(posCon = mean(size)) %>%
  select(posCon, assay_date, species, preparation, other) %>%
  unique()


# normalize to negative and positive controls for each assay date 
normalized <- development_raw %>%
  left_join(., NegCon, by = c('assay_date', 'species', 'preparation', 'other', 'solvent')) %>%
  left_join(., PosCon, by = c('assay_date', 'species', 'preparation', 'other')) %>%
  mutate(norm_dev = (size - posCon) / (negCon - posCon)) %>%
  mutate(species = factor(species, 
                          levels = c('AeLVP', 'AlboMO', 'Cxp'),
                          labels = c('*Ae. aegypti* LVP', '*Ae. albopictus* MO', '*Cx. pipiens*')))

norm_dev <- normalized %>%
  group_by(species, conc, treatment, assay_date, preparation, other) %>%
  summarise(mean_dev = median(norm_dev, na.rm = TRUE)) %>%
  ungroup()

normalized_dev <- normalized %>%
  left_join(., norm_dev) %>%
  filter(treatment != 'missed')

# prune outliers 

normalized_dev <- normalized_dev %>%
  # adding assay_date here doesn't change anything 
  group_by(preparation, treatment, conc, other, species) %>%
  mutate(IQR = IQR(norm_dev),
         O_upper = quantile(norm_dev, probs=c( .75), na.rm = FALSE)+1.5*IQR,  
         O_lower = quantile(norm_dev, probs=c( .25), na.rm = FALSE)-1.5*IQR  
  ) %>% 
  filter(O_lower <= norm_dev & norm_dev <= O_upper) %>%
  ungroup()

# development plotting ---- 

controls <- normalized_dev %>%
  filter(preparation == '1/well',
         treatment %in% c('DMSO', 'Water', 'IVM')) %>%
  mutate(treatment = case_when(
    treatment == 'DMSO' ~ 'Negative Control',
    treatment == 'IVM' ~ 'Positive Control',
    treatment == 'Water' ~ 'Negative Control',
    TRUE ~ treatment
  ))

(Fig2C_neg <- controls %>%
    filter(treatment == 'Negative Control') %>%
    ggplot(aes(x = treatment, y = norm_dev, color = species)) +
    geom_hline(yintercept = 1, linetype="dashed", color = "black", size=0.5) +
    geom_hline(yintercept = 0, linetype="dashed", color = "black", size=0.5) +
    geom_point(aes(colour = species), size = 0.75, alpha = 0.5, position=position_dodge(width=0.3)) +
    facet_grid(cols = vars(treatment), rows = vars(other), scales = 'free') +
    scale_color_manual(name = "Species",
                       values = c('*Ae. aegypti* LVP' = 'chartreuse3', '*Ae. albopictus* MO' = 'goldenrod3', '*Cx. pipiens*' = 'darkorchid')) +
    ylim(-0.5,2) +
    theme_zlab_white() +
    theme(axis.line = element_line(color = "black", linewidth = 0.5),
          legend.position = "none",
          axis.text.x = element_blank(),
          strip.text.y = element_blank()) +
    ylab("Size") +
    xlab(""))

normalized_dev <- normalized_dev %>%
  mutate(
    treatment_lab = case_when(
      treatment == "Altosid" ~ "Altosid (mg/L)",
      treatment == "Chlorfenapyr" ~ "Chlorfenapyr (µM)",
      treatment == "Spinosad" ~ "Spinosad (ppm)",
      treatment == "Vectolex" ~ "Vectolex (ppm)",
      TRUE ~ treatment
    )
  )


# stats data prep 

# temp <- normalized_mot %>%
#   select(treatment) %>%
#   unique(.)

normalized_dev <- normalized_dev %>%
  filter(preparation == '1/well',
         !(treatment == "Altosid" & conc == '100mg/L')) 

altosid <- normalized_dev %>%
  filter(treatment %in% c("Water", "Altosid"))
chlorfenapyr <- normalized_dev %>%
  filter(treatment %in% c("DMSO", "Chlorfenapyr"))
spinosad <- normalized_dev %>%
  filter(treatment %in% c("Water", "Spinosad"))
vectolex <- normalized_dev %>%
  filter(treatment %in% c("Water", "Vectolex"))

# stats calculations 

stats.altosid <- compare_means(norm_dev~conc, data = altosid, group.by= c("species"),
                               method = "t.test", ref.group = "0uM") %>%
  mutate(p.signif = factor(p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***")),
         treatment_lab = 'Altosid (mg/L)')
stats.chlorfenapyrC <- compare_means(norm_dev~conc, data = chlorfenapyr, group.by= c("species"),
                                    method = "t.test", ref.group = "0.5p") %>%
  mutate(p.signif = factor(p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***")),
         treatment_lab = 'Chlorfenapyr (µM)')
stats.chlorfenapyrA <- compare_means(norm_dev~conc, data = chlorfenapyr, group.by= c("species"),
                                     method = "t.test", ref.group = "1p") %>%
  mutate(p.signif = factor(p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***")),
         treatment_lab = 'Chlorfenapyr (µM)')
stats.spinosad <- compare_means(norm_dev~conc, data = spinosad, group.by= c("species"),
                                method = "t.test", ref.group = "0uM") %>%
  mutate(p.signif = factor(p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***")),
         treatment_lab = 'Spinosad (ppm)')
stats.vectolex <- compare_means(norm_dev~conc, data = vectolex, group.by= c("species"),
                                method = "t.test", ref.group = "0uM") %>%
  mutate(p.signif = factor(p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***")),
         treatment_lab = 'Vectolex (ppm)')

# stats merge 
# final stats: order conc, other; set y values 

stats <- rbind(stats.altosid, stats.chlorfenapyrC, stats.chlorfenapyrA, stats.spinosad, stats.vectolex) %>%
  mutate(conc = group2,
         y = case_when(
           species == '*Cx. pipiens*' ~ 1.7,
           species == '*Ae. albopictus* MO' ~ 1.8,
           species == '*Ae. aegypti* LVP' ~ 1.9
         )) %>%
  mutate(conc = factor(conc, levels = c('1p', '0uM', '0.001ppm', '0.01ppm', '0.1ppm', '1ppm',
                                        '10ppm', '100ppm',
                                        '0.01mg/L', '0.1mg/L',
                                        '1mg/L', '10mg/L', '100mg/L', '1000mg/L',
                                        '500pM', '1nM',
                                        '5nM', '10nM', '50nM', '100nM', '500nM', '1uM',
                                        '5uM', '10uM', '50uM', '100uM'),
                       labels = c('1%', '0µM', '0.001', '0.01', '0.1', '1',
                                  '10', '100', 
                                  '0.01', '0.1',
                                  '1', '10', '100', '1000',
                                  '0.0005', '0.001',
                                  '0.005', '0.01', '0.05', '0.1', '0.5', '1',
                                  '5', '10', '50', '100'),
                       ordered = TRUE),
         conc= factor(conc, levels = c('1%', '0µM', '0.0005', '0.001',
                                       '0.005', '0.01', '0.05', '0.1',
                                       '0.5', '1', '5', '10', '50', '100',
                                       '1000'),
                      labels = c(2, 0, 0.0005, 0.001, 0.005, 0.01, 0.05,
                                 0.1, 0.5, 1, 5, 10, 50, 100, 1000),
                      ordered = TRUE)) 

(Fig2C_main <- normalized_dev %>%
    filter(!treatment %in% c('DMSO', 'Water', 'IVM')) %>%
    mutate(conc = factor(conc, levels = c('1p', '0uM', '0.001ppm', '0.01ppm', '0.1ppm', '1ppm',
                                          '10ppm', '100ppm',
                                          '0.01mg/L', '0.1mg/L',
                                          '1mg/L', '10mg/L', '100mg/L', '1000mg/L',
                                          '500pM', '1nM',
                                          '5nM', '10nM', '50nM', '100nM', '500nM', '1uM',
                                          '5uM', '10uM', '50uM', '100uM'),
                         labels = c('1%', '0µM', '0.001', '0.01', '0.1', '1',
                                    '10', '100', 
                                    '0.01', '0.1',
                                    '1', '10', '100', '1000',
                                    '0.0005', '0.001',
                                    '0.005', '0.01', '0.05', '0.1', '0.5', '1',
                                    '5', '10', '50', '100')),
           other = factor(other, levels = c('1hr', '24hr', '48hr', '72hr'),
                          ordered = TRUE),
           conc= factor(conc, levels = c('1%', '0µM', '0.0005', '0.001',
                                         '0.005', '0.01', '0.05', '0.1',
                                         '0.5', '1', '5', '10', '50', '100',
                                         '1000'),
                        labels = c(2, 0, 0.0005, 0.001, 0.005, 0.01, 0.05,
                                   0.1, 0.5, 1, 5, 10, 50, 100, 1000),
                        ordered = TRUE)) %>%
    ggplot(aes(x = conc, y = norm_dev, color = species)) +
    geom_hline(yintercept = 1, linetype="dashed", color = "black", size=0.5) +
    geom_hline(yintercept = 0, linetype="dashed", color = "black", size=0.5) +
    geom_point(aes(colour = species), size = 0.75, alpha = 0.5, position=position_dodge(width=0.3)) +
    geom_line(aes(group = species, colour = species), size = 0.1, alpha = 0.05) +
    geom_smooth(aes(group = species, colour = species), size = 0.75, alpha = 1, se=F) +
    geom_text(
      data = stats, aes(label = p.signif, y=y, color = species),
      size = 3, show.legend = FALSE) +
    facet_grid(cols = vars(treatment_lab), rows = vars(other), scales = 'free') +
    scale_x_discrete(
      breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000)
    ) +
    scale_color_manual(name = "Species",
                       values = c('*Ae. aegypti* LVP' = 'chartreuse3', '*Ae. albopictus* MO' = 'goldenrod3', '*Cx. pipiens*' = 'darkorchid')) +
    ylim(-0.5,2) +
    theme_zlab_white() +
    theme(axis.line = element_line(color = "black", linewidth = 0.5),
          axis.text.x = element_text(angle = 0),
          legend.position = 'top',
          legend.direction = 'horizontal',
          strip.text.y = element_blank()) +
    ylab("") +
    xlab("Concentration"))

(Fig2C_pos <- controls %>%
    filter(treatment == 'Positive Control') %>%
    ggplot(aes(x = treatment, y = norm_dev, color = species)) +
    geom_hline(yintercept = 1, linetype="dashed", color = "black", size=0.5) +
    geom_hline(yintercept = 0, linetype="dashed", color = "black", size=0.5) +
    geom_point(aes(colour = species), size = 0.75, alpha = 0.5, position=position_dodge(width=0.3)) +
    facet_grid(cols = vars(treatment), rows = vars(other), scales = 'free') +
    scale_color_manual(name = "Species",
                       values = c('*Ae. aegypti* LVP' = 'chartreuse3', '*Ae. albopictus* MO' = 'goldenrod3', '*Cx. pipiens*' = 'darkorchid')) +
    theme_zlab_white() +
    ylim(-0.5,2) +
    theme(axis.line = element_line(color = "black", linewidth = 0.5),
          legend.position = "none",
          axis.text.x = element_blank()) +
    ylab("") +
    xlab(""))


# rearrange plots and get legend 

(Fig2C_neg_clean <- Fig2C_neg +
  guides(color = "none") +
  theme(
    panel.spacing = unit(0.1, "lines"),
    strip.text.y = element_blank(),
    plot.margin = margin(0,0,0,0)
  ))

(Fig2C_main_clean <- Fig2C_main +
  guides(
    color = guide_legend(override.aes = list(size = 1.5, linetype = 1, shape = NA))
  ) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    plot.margin = margin(0,0,0,0)
  ))

(Fig2C_pos_clean <- Fig2C_pos +
  guides(color = "none") +
  theme(
    panel.spacing = unit(0.1, "lines"),
    plot.margin = margin(0,0,0,0)
  ))

# --- Combine plots in a row with a single shared legend ---
(development_final <- Fig2C_neg_clean + Fig2C_main_clean + Fig2C_pos_clean +
    plot_layout(nrow = 1, widths = c(0.6, 2.5, 0.6), guides = "collect") &
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_markdown(size = 8),
      legend.box.spacing = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0)
    ))

# (development_final <- plot_grid(Fig2C_neg_clean, NULL, Fig2C_main_clean, NULL, Fig2C_pos_clean, nrow = 1,
#                              rel_widths = c(0.35, -0.1, 1, -0.1, 0.4))) 

# merge and save Figure ---------------

(schematic <- image_read(here("../images/mos-F2.jpeg")))
(Fig2A <- ggdraw() + draw_image(schematic))

(Figure2 <- plot_grid(Fig2A, motility_final, NULL, development_final, ncol = 1, 
                      rel_heights = c(1.5, 5, -0.1, 2.25), 
                      rel_widths = c(1.5, 1, 1, 1),
                      labels = c('A', 'B', '', 'C')))

ggsave(here('../Figure2/Fig2.pdf'), Figure2, width = 9, height = 11, units = 'in')


# motility plot (Supplmental motility) -------------------------------------------------------

# renormalize motility data (overwrite database names)

# normalize to negative controls 
NegCon <- motility_raw %>%
  filter(treatment %in% c('DMSO', 'Water')) %>%
  group_by(assay_date, species, preparation, other, solvent) %>%
  mutate(negCon = mean(optical_flow)) %>%
  select(negCon, assay_date, species, preparation, other, solvent) %>%
  unique()

# find combined positive control for each assay date 
PosCon <- motility_raw %>%
  filter(treatment == 'IVM') %>%
  group_by(assay_date, species, preparation, other) %>%
  mutate(posCon = mean(optical_flow)) %>%
  select(posCon, assay_date, species, preparation, other) %>%
  unique()


# normalize to negative and positive controls for each assay date 
normalized <- motility_raw %>%
  left_join(., NegCon, by = c('assay_date', 'species', 'preparation', 'other', 'solvent')) %>%
  left_join(., PosCon, by = c('assay_date', 'species', 'preparation', 'other')) %>%
  mutate(norm_mot = (optical_flow - posCon) / (negCon - posCon)) %>%
  mutate(species = factor(species, 
                          levels = c('AeLVP', 'AlboMO', 'Cxp'),
                          labels = c('*Ae. aegypti* LVP', '*Ae. albopictus* MO', '*Cx. pipiens*')))

norm_means <- normalized %>%
  group_by(species, conc, treatment, assay_date, preparation, other) %>%
  summarise(mean_mot = median(norm_mot, na.rm = TRUE)) %>%
  ungroup()

normalized_mot <- normalized %>%
  left_join(., norm_means) %>%
  filter(treatment != 'missed')

# prune outliers 

normalized_mot <- normalized_mot %>%
  # adding assay_date here doesn't change anything 
  group_by(preparation, treatment, conc, other, species) %>%
  mutate(IQR = IQR(norm_mot),
         O_upper = quantile(norm_mot, probs=c( .75), na.rm = FALSE)+1.5*IQR,  
         O_lower = quantile(norm_mot, probs=c( .25), na.rm = FALSE)-1.5*IQR  
  ) %>% 
  filter(O_lower <= norm_mot & norm_mot <= O_upper) %>%
  ungroup()

controls1per <- normalized_mot %>%
  filter(preparation == '1/well',
         treatment %in% c('DMSO', 'Water', 'IVM')) %>%
  mutate(treatment = case_when(
    treatment == 'DMSO' ~ 'Negative Control',
    treatment == 'IVM' ~ 'Positive Control',
    treatment == 'Water' ~ 'Negative Control',
    TRUE ~ treatment
  ))

(FigS3_neg <- controls1per %>%
    filter(treatment == 'Negative Control',
           other != 'other') %>%
    ggplot(aes(x = treatment, y = norm_mot, color = species)) +
    geom_hline(yintercept = 1, linetype="dashed", color = "black", size=0.5) +
    geom_hline(yintercept = 0, linetype="dashed", color = "black", size=0.5) +
    geom_point(aes(colour = species), size = 0.75, alpha = 0.5, position=position_dodge(width=0.3)) +
    facet_grid(cols = vars(treatment), rows = vars(other), scales = 'free') +
    scale_color_manual(name = "Species",
                       values = c('*Ae. aegypti* LVP' = 'chartreuse3', '*Ae. albopictus* MO' = 'goldenrod3', '*Cx. pipiens*' = 'darkorchid')) +
    ylim(-0.5,2) +
    theme_zlab_white() +
    theme(axis.line = element_line(color = "black", linewidth = 0.5),
          legend.position = "none",
          axis.text.x = element_blank(),
          strip.text.y = element_blank()) +
    ylab("Motility") +
    xlab(""))

normalized_mot <- normalized_mot %>%
  mutate(
    treatment_lab = case_when(
      treatment == "Altosid" ~ "Altosid (mg/L)",
      treatment == "Chlorfenapyr" ~ "Chlorfenapyr (µM)",
      treatment == "Spinosad" ~ "Spinosad (ppm)",
      treatment == "Vectolex" ~ "Vectolex (ppm)",
      TRUE ~ treatment
    )
  )


# stats data prep 

# temp <- normalized_mot %>%
#   select(treatment) %>%
#   unique(.)

normalized_mot <- normalized_mot %>%
  # remove altosid concentration where larvae weren't visible
  filter(preparation == '1/well',
         !(treatment == "Altosid" & conc == '100mg/L'),
         !(treatment == "Altosid" & conc == '1000mg/L'),
         !(treatment == "Altosid" & conc == '1000mg/ml')) 

altosid <- normalized_mot %>%
  filter(treatment %in% c("Water", "Altosid"))
chlorfenapyr <- normalized_mot %>%
  filter(treatment %in% c("DMSO", "Chlorfenapyr"))
spinosad <- normalized_mot %>%
  filter(treatment %in% c("Water", "Spinosad"))
vectolex <- normalized_mot %>%
  filter(treatment %in% c("Water", "Vectolex"))

# stats calculations 

stats.altosid <- compare_means(norm_mot~conc, data = altosid, group.by= c("species", "other"),
                               method = "t.test", ref.group = "0uM") %>%
  mutate(p.signif = factor(p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***")),
         treatment_lab = 'Altosid (mg/L)')
stats.chlorfenapyrC <- compare_means(norm_mot~conc, data = chlorfenapyr, group.by= c("species", "other"),
                                     method = "t.test", ref.group = "0.5p") %>%
  mutate(p.signif = factor(p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***")),
         treatment_lab = 'Chlorfenapyr (µM)')
stats.chlorfenapyrA <- compare_means(norm_mot~conc, data = chlorfenapyr, group.by= c("species", "other"),
                                     method = "t.test", ref.group = "1p") %>%
  mutate(p.signif = factor(p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***")),
         treatment_lab = 'Chlorfenapyr (µM)')
stats.spinosad <- compare_means(norm_mot~conc, data = spinosad, group.by= c("species", "other"),
                                method = "t.test", ref.group = "0uM") %>%
  mutate(p.signif = factor(p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***")),
         treatment_lab = 'Spinosad (ppm)')
stats.vectolex <- compare_means(norm_mot~conc, data = vectolex, group.by= c("species", "other"),
                                method = "t.test", ref.group = "0uM") %>%
  mutate(p.signif = factor(p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***")),
         treatment_lab = 'Vectolex (ppm)')

# stats merge 
# final stats: order conc, other; set y values 

stats <- rbind(stats.altosid, stats.chlorfenapyrC, stats.chlorfenapyrA, stats.spinosad, stats.vectolex) %>%
  mutate(conc = group2,
         y = case_when(
           species == '*Cx. pipiens*' ~ 1.7,
           species == '*Ae. albopictus* MO' ~ 1.8,
           species == '*Ae. aegypti* LVP' ~ 1.9
         )) %>%
  mutate(conc = factor(conc, levels = c('1p', '0uM', '0.001ppm', '0.01ppm', '0.1ppm', '1ppm',
                                        '10ppm', '100ppm',
                                        '0.01mg/L', '0.1mg/L',
                                        '1mg/L', '10mg/L', '100mg/L', '1000mg/L',
                                        '500pM', '1nM',
                                        '5nM', '10nM', '50nM', '100nM', '500nM', '1uM',
                                        '5uM', '10uM', '50uM', '100uM'),
                       labels = c('1%', '0µM', '0.001', '0.01', '0.1', '1',
                                  '10', '100', 
                                  '0.01', '0.1',
                                  '1', '10', '100', '1000',
                                  '0.0005', '0.001',
                                  '0.005', '0.01', '0.05', '0.1', '0.5', '1',
                                  '5', '10', '50', '100'),
                       ordered = TRUE),
         conc= factor(conc, levels = c('1%', '0µM', '0.0005', '0.001',
                                       '0.005', '0.01', '0.05', '0.1',
                                       '0.5', '1', '5', '10', '50', '100',
                                       '1000'),
                      labels = c(2, 0, 0.0005, 0.001, 0.005, 0.01, 0.05,
                                 0.1, 0.5, 1, 5, 10, 50, 100, 1000),
                      ordered = TRUE)) %>%
  filter(other != '1hr')

(FigS3_main <- normalized_mot %>%
    filter(!treatment %in% c('DMSO', 'Water', 'IVM'),
           other != '1hr') %>%
    mutate(conc = factor(conc, levels = c('1p', '0uM', '0.001ppm', '0.01ppm', '0.1ppm', '1ppm',
                                          '10ppm', '100ppm',
                                          '0.01mg/L', '0.1mg/L',
                                          '1mg/L', '10mg/L', '100mg/L', '1000mg/L',
                                          '500pM', '1nM',
                                          '5nM', '10nM', '50nM', '100nM', '500nM', '1uM',
                                          '5uM', '10uM', '50uM', '100uM'),
                         labels = c('1%', '0µM', '0.001', '0.01', '0.1', '1',
                                    '10', '100', 
                                    '0.01', '0.1',
                                    '1', '10', '100', '1000',
                                    '0.0005', '0.001',
                                    '0.005', '0.01', '0.05', '0.1', '0.5', '1',
                                    '5', '10', '50', '100')),
           other = factor(other, levels = c('1hr', '24hr', '48hr', '72hr'),
                          ordered = TRUE),
           conc= factor(conc, levels = c('1%', '0µM', '0.0005', '0.001',
                                         '0.005', '0.01', '0.05', '0.1',
                                         '0.5', '1', '5', '10', '50', '100',
                                         '1000'),
                        labels = c(2, 0, 0.0005, 0.001, 0.005, 0.01, 0.05,
                                   0.1, 0.5, 1, 5, 10, 50, 100, 1000),
                        ordered = TRUE)) %>%
    ggplot(aes(x = conc, y = norm_mot, color = species)) +
    geom_hline(yintercept = 1, linetype="dashed", color = "black", size=0.5) +
    geom_hline(yintercept = 0, linetype="dashed", color = "black", size=0.5) +
    geom_point(aes(colour = species), size = 0.75, alpha = 0.5, position=position_dodge(width=0.3)) +
    geom_line(aes(group = species, colour = species), size = 0.1, alpha = 0.05) +
    geom_smooth(aes(group = species, colour = species), size = 0.75, alpha = 1, se=F) +
    geom_text(
      data = stats, aes(label = p.signif, y=y, color = species),
      size = 3, show.legend = FALSE) +
    facet_grid(cols = vars(treatment_lab), rows = vars(other), scales = 'free') +
    scale_x_discrete(
      breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000)
    ) +
    scale_color_manual(name = "Species",
                       values = c('*Ae. aegypti* LVP' = 'chartreuse3', '*Ae. albopictus* MO' = 'goldenrod3', '*Cx. pipiens*' = 'darkorchid')) +
    ylim(-0.5,2) +
    theme_zlab_white() +
    theme(axis.line = element_line(color = "black", linewidth = 0.5),
          axis.text.x = element_text(angle = 0),
          legend.position = 'none',
          #legend.direction = 'horizontal',
          strip.text.y = element_blank()) +
    ylab("") +
    xlab("Concentration"))

(FigS3_pos <- controls1per %>%
    filter(treatment == 'Positive Control',
           other != '1hr') %>%
    ggplot(aes(x = treatment, y = norm_mot, color = species)) +
    geom_hline(yintercept = 1, linetype="dashed", color = "black", size=0.5) +
    geom_hline(yintercept = 0, linetype="dashed", color = "black", size=0.5) +
    geom_point(aes(colour = species), size = 0.75, alpha = 0.5, position=position_dodge(width=0.3)) +
    facet_grid(cols = vars(treatment), rows = vars(other), scales = 'free') +
    scale_color_manual(name = "Species",
                       values = c('*Ae. aegypti* LVP' = 'chartreuse3', '*Ae. albopictus* MO' = 'goldenrod3', '*Cx. pipiens*' = 'darkorchid')) +
    theme_zlab_white() +
    ylim(-0.5,2) +
    theme(axis.line = element_line(color = "black", linewidth = 0.5),
          legend.position = "none",
          axis.text.x = element_blank()) +
    ylab("") +
    xlab(""))

# rearrange plots 

FigS3_neg_clean <- FigS3_neg +
  guides(color = "none") +
  theme(
    panel.spacing = unit(0.1, "lines"),
    strip.text.y = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

FigS3_main_clean <- FigS3_main +
  guides(
    color = guide_legend(override.aes = list(size = 1.5, linetype = 1, shape = NA))
  ) +
  theme(
    panel.spacing = unit(0.1, "lines"),
    plot.margin = margin(0,0,0,0)
  )

FigS3_pos_clean <- FigS3_pos +
  guides(color = "none") +
  theme(
    panel.spacing = unit(0.1, "lines"),
    plot.margin = margin(0,0,0,0)
  )

# --- Combine plots in a row with a single shared legend ---
(motility_supp <- FigS3_neg_clean + FigS3_main_clean + FigS3_pos_clean +
    plot_layout(nrow = 1, widths = c(0.8, 4, 0.8), guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_markdown(size = 8),
      legend.box.spacing = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0)
    ))


ggsave(here('../Supp3/Supp3.pdf'), motility_supp, width = 10, height = 7, units = 'in')

