################################################################################
#   Diatoms vs dinoflagellates: visualising the phytoplankton community with   #
#                          interactive plots | R scripts                       #
################################################################################

# Script to investigate the differences between regions #
# Last modified: 2026/07/30 by V. POCHIC
# Based on an original script by J-Y. DIAS.

# /!\ WARNING /!\
# This script is a supplementary data visualisation associated with the 
# peer-reviewedarticle Dias et al. 2026 (https://doi.org/10.1093/ismeco/ycag174)
# published in ISME Communications
# It uses the same database than the one analysed in the article:
# REPHY database (1987-2022) - https://doi.org/10.17882/47248
# To make it function properly, you first need to execute scripts 0 and 1 in the
# public github repository associated with the article 
# (https://github.com/J-YDi/Diatoms-vs-Dinoflagellates)
# /!\ THANK YOU FOR YOUR ATTENTION /!\

## Description ####

# The goal of this script is to provide an interactive data visualisation for
# Figure 2 in Dias et al. (2026) which represents the evolution of the phyto-
# plankton community in 3 regions of the French coastline, between March 2007
# and August 2022.

# The results are shown and briefly discussed in a dedicated blog post on
# https://phycoplankton.fr/2026/07/30/visualising-16-years-of-phytoplankton-community-dynamics/

# Please note that this work, contrary to the published article, has not been
# peer-reviewed. I (V. Pochic) take full responsibility for any errors in these 
# lines of code or in the blog post, and my co-authors couldn't be blamed for 
# them.

## Required packages ####
# (maybe not all of them are really necessary but eh)
library(ggplot2)
library(plotly)
library(ggthemes)
library(readr)
library(dplyr)
library(tidyr)
library(FactoMineR)
library(factoextra)
library(cowplot)
library(DescTools)
library(deeptime)
library(htmlwidgets)

####------------------------------------------------------------------------####
## Import data ####
# You need to 

data <- read.csv2("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_final.csv", 
                   header = TRUE, fileEncoding = 'ISO-8859-1')
# Adding season info
data <- data |>
  mutate(season = case_when(Month %in% c(12, 01, 02) ~ "Winter",
                            Month %in% c(03, 04, 05) ~ "Spring",
                            Month %in% c(06, 07, 08) ~ "Summer",
                            Month %in% c(09, 10, 11) ~ "Fall", TRUE ~ NA_character_))


data$region <- as.factor(data$region)

# Don't consider the Pertuis Sea region
data <- filter(data, region != "4-Pertuis Sea")

### Graph to show phytoplankton composition by season and cluster #####
# By Genus
# Compute the mean of all genus by season and cluster
data_graph <- data[,c(3,353,24:328)] %>%
  group_by(season, region) %>%
  summarise_all(~mean(., na.rm = TRUE))

# Make it relative 
data_graph$Abdtot <- rowSums(data_graph[,c(3:307)],na.rm=T)
data_graph[,c(3:307)] <- data_graph[,c(3:307)]/data_graph$Abdtot

datag <- pivot_longer(data = data_graph,cols = Actinoptychus:Coscinodiscophycidae,names_to = "Taxon")
datag$Abdtot <- NULL

# Group by season and region, and select the 4 most abundant phyla for each region
# Determine the most abundant phyla for each region and season

data_graph <- datag %>%
  group_by(season, region, Taxon) %>%
  summarise(abondance = sum(value,na.rm=T)) %>%
  ungroup() %>%
  group_by(season, region) %>%
  mutate(rank = rank(desc(abondance)))

# Save that table if you wish:
# write.csv2(data_graph,
# file="output/tableaux/taxa_mostabundant_region_season.csv", row.names = FALSE,
# fileEncoding = 'ISO-8859-1')

# Select those most abundant taxa (16 in total):
selection_taxon <- c("Skeletonema","Pseudo.nitzschia","Chaetoceros","Chaetocerotaceae","Nitzschia","Cryptophyceae",
                      "Cryptomonadales","Cylindrotheca","Leptocylindrus","Akashiwo","Phaeocystis","Akashiwo",
                      "Phaeocystis","Asterionellopsis","Chrysochromulina","Azadinium")
datataxon <- data |>
  select(Month,Year,region,selection_taxon)
# Agglomerate Chaetoceros with Chaetocerotaceae and Cryptophyceae with Cryptomonadales
datataxon$Chaetocerotaceae <- rowSums(datataxon[,c("Chaetoceros","Chaetocerotaceae")],na.rm=T)
datataxon$Cryptophyceae <- rowSums(datataxon[,c("Cryptophyceae","Cryptomonadales")],na.rm=T)
datataxon <- select(datataxon,-c(Chaetoceros,Cryptomonadales))

# Create the "Other" category
dataothers <- data |>
  select(-selection_taxon) |>
  select(region,Actinoptychus:Coscinodiscophycidae)

dataothers$Others <- rowSums(select(dataothers,Actinoptychus:Coscinodiscophycidae),na.rm = T)
dataothers <- select(dataothers,Others)

# Bind them
datagraph <- cbind(datataxon,dataothers)

# Mean by month, region, year
data_graph <- datagraph %>%
  group_by(Month,Year, region) %>%
  summarise_all(~mean(., na.rm = TRUE))

## Cells per L and log10 versions ####
data_graph$Abdtot <- rowSums(data_graph[,c(4:16)],na.rm=T)

# Convert to date format
date_string <- paste(data_graph$Year, data_graph$Month, "01", sep = "-")
data_graph$Date <- as.Date(date_string,format = "%Y-%m-%d")
data_graph$MonthYear <- format(data_graph$Date, "%Y-%m")

datag <- pivot_longer(data = data_graph,cols = Skeletonema:Others,names_to = "Taxon") %>%
  # get rid of NaNs and zeros (better for the interactive plot) and format it 
  # to have the desired number of digits
  mutate(`Cell density` = ifelse(value == 0 | value == 'NaN', NA, 
                                 as.numeric(formatC(value, 
                                                    format = 'E', digits = 2)))) %>%
  # compute log10 of value
  mutate(`log10(C)` = ifelse(value == 0 | value == 'NaN', NA,
                             as.numeric(formatC(log10(value), 
                                                format = 'E', digits = 2)))) %>%
  # and log2 of value
  mutate(`log2(C)` = ifelse(value == 0 | value == 'NaN', NA,
                            as.numeric(formatC(log2(value), 
                                               format = 'E', digits = 2)))) %>%
  # Change the name of Pseudo-nitzschia so it's written correctly
  mutate(Taxon = ifelse(Taxon == 'Pseudo.nitzschia', 'Pseudo-nitzschia',
                        Taxon)) %>%
  # Change the name of the MonthYear variable so it reads better
  mutate(`Year-Month` = MonthYear)

datag$Taxon <- factor(datag$Taxon, levels = c("Asterionellopsis","Chaetocerotaceae","Cylindrotheca","Leptocylindrus"
                                                ,"Nitzschia","Pseudo-nitzschia", "Skeletonema","Akashiwo","Azadinium",
                                                "Chrysochromulina","Phaeocystis","Cryptophyceae","Others"))

df_color <- data.frame(
  name = c('1-Mediterranean sea', '2-Eastern Channel - North Sea',
           '3-Atlantic - Western Channel'),
  color = c('red3', 
            'dodgerblue2',
            'violet'))

# Graph by region (facets), with real cell numbers
abundance <- ggplot(datag) +
  geom_col(aes(x = `Year-Month`, y = `Cell density`, 
               fill = Taxon), position = "stack", na.rm = FALSE, width = 1) +
  facet_wrap_color(~region, scales = "free_y", ncol = 1, colors = df_color) +
  scale_x_discrete(breaks = c('2007-06', '2008-06', '2009-06', '2010-06', 
                              '2011-06', '2012-06', '2013-06', '2014-06',
                              '2015-06', '2016-06', '2017-06', '2018-06',
                              '2019-06', '2020-06', '2021-06', '2022-06'),
                   labels = c('2007', '2008', '2009', '2010', 
                              '2011', '2012', '2013', '2014',
                              '2015', '2016', '2017', '2018',
                              '2019', '2020', '2021', '2022')
                   ) +
  geom_vline(data = subset(datag, format(Date, "%m") == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", linewidth = 0.5,linetype = "dashed",
             show.legend = FALSE)+
  theme_classic() +
  theme(# Axes
        axis.text.y = element_text(),
        axis.text.x = element_text(),
        axis.ticks.x = element_blank(),
        # Legend
        legend.position = "right",
        legend.background = element_rect(color = 'grey15', linewidth = .45,
                                         fill = '#F8F7F5'),
        # Plot background
        panel.background = element_rect(fill = '#F8F7F5'),
        plot.background = element_rect(fill = '#F8F7F5'),
        # strip
        strip.background = element_rect(fill = 'grey40'),
        strip.text = element_text(size = 12, color = '#F8F7F5')
        )+
  # Labels and colors
  labs(x = 'Year', y = NULL, fill = 'Taxon: '
       )+
  scale_fill_manual(values = c(
    "Asterionellopsis"   = "#2B4561",
    "Chaetocerotaceae"   = "#76A7E2",
    "Cylindrotheca"      = "#2E6CD9",
    "Leptocylindrus"     = "#000E53",
    "Nitzschia"          = "#377185",
    "Pseudo-nitzschia"   = "#BBD4F2",
    "Skeletonema"        = "blue",
    "Akashiwo"           = "#B2DF8A",
    "Azadinium"          = "chartreuse4",
    "Chrysochromulina"   = "#FBB646",
    "Phaeocystis"        = "gold1",
    "Cryptophyceae"      = "#FC4D6B",
    "Others"             = "grey"
  ))

# Final graph:
abundance
plotly_abundance <- ggplotly(abundance)
saveWidget(plotly_abundance, file="output/HTML_figures/rephyto_plotly_abundance.html")

# Now with log10 of cell numbers
log10_plot <- ggplot(datag) +
  geom_col(aes(x = `Year-Month`, y = `log10(C)`, 
               fill = Taxon), position = "stack", na.rm = FALSE, width = 1) +
  facet_wrap_color(~region, scales = "free_y", ncol = 1, colors = df_color) +
  scale_x_discrete(breaks = c('2007-06', '2008-06', '2009-06', '2010-06', 
                              '2011-06', '2012-06', '2013-06', '2014-06',
                              '2015-06', '2016-06', '2017-06', '2018-06',
                              '2019-06', '2020-06', '2021-06', '2022-06'),
                   labels = c('2007', '2008', '2009', '2010', 
                              '2011', '2012', '2013', '2014',
                              '2015', '2016', '2017', '2018',
                              '2019', '2020', '2021', '2022')
  ) +
  geom_vline(data = subset(datag, format(Date, "%m") == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", linewidth = 0.5,linetype = "dashed",
             show.legend = FALSE)+
  theme_classic() +
  theme(# Axes
    axis.text.y = element_text(),
    axis.text.x = element_text(),
    axis.ticks.x = element_blank(),
    # Legend
    legend.position = "right",
    legend.background = element_rect(color = 'grey15', linewidth = .45,
                                     fill = '#F8F7F5'),
    # Plot background
    panel.background = element_rect(fill = '#F8F7F5'),
    plot.background = element_rect(fill = '#F8F7F5'),
    # strip
    strip.background = element_rect(fill = 'grey40'),
    strip.text = element_text(size = 12, color = '#F8F7F5')
  )+
  # Labels and colors
  labs(x = 'Year', y = NULL, fill = 'Taxon: '
  )+
  scale_fill_manual(values = c(
    "Asterionellopsis"   = "#2B4561",
    "Chaetocerotaceae"   = "#76A7E2",
    "Cylindrotheca"      = "#2E6CD9",
    "Leptocylindrus"     = "#000E53",
    "Nitzschia"          = "#377185",
    "Pseudo-nitzschia"   = "#BBD4F2",
    "Skeletonema"        = "blue",
    "Akashiwo"           = "#B2DF8A",
    "Azadinium"          = "chartreuse4",
    "Chrysochromulina"   = "#FBB646",
    "Phaeocystis"        = "gold1",
    "Cryptophyceae"      = "#FC4D6B",
    "Others"             = "grey"
  ))

# Final graph:
log10_plot
plotly_log10 <- ggplotly(log10_plot)
saveWidget(plotly_log10, file="output/HTML_figures/rephyto_plotly_log10.html")

## Relative abundance version ####
# When we switch to relative
data_graph$Abdtot <- rowSums(data_graph[,c(4:16)],na.rm=T)
data_graph[,c(4:16)] <- data_graph[,c(4:16)]/data_graph$Abdtot

data_graph$Abdtot <- NULL

# Convert to date format
date_string <- paste(data_graph$Year, data_graph$Month, "01", sep = "-")
data_graph$Date <- as.Date(date_string,format = "%Y-%m-%d")
data_graph$MonthYear <- format(data_graph$Date, "%Y-%m")

datag <- pivot_longer(data = data_graph,cols = Skeletonema:Others,names_to = "Taxon") %>%
  # get rid of NaNs and zeros (better for the interactive plot)
  # and format the numbers (2 decimals)
  mutate(`Proportion` = ifelse(value == 'NaN' | value == 0, NA,
                               as.numeric(formatC((value), 
                                                  format = 'f', digits = 4))
                               )) %>%
  # Change the name of Pseudo-nitzschia so it's written correctly
  mutate(Taxon = ifelse(Taxon == 'Pseudo.nitzschia', 'Pseudo-nitzschia',
                        Taxon)) %>%
  # Change the name of the MonthYear variable so it reads better
  mutate(`Year-Month` = MonthYear)

datag$Taxon <- factor(datag$Taxon, levels = c("Asterionellopsis","Chaetocerotaceae","Cylindrotheca","Leptocylindrus"
                                              ,"Nitzschia","Pseudo-nitzschia", "Skeletonema","Akashiwo","Azadinium",
                                              "Chrysochromulina","Phaeocystis","Cryptophyceae","Others"))

df_color <- data.frame(
  name = c('1-Mediterranean sea', '2-Eastern Channel - North Sea',
           '3-Atlantic - Western Channel'),
  color = c('red3', 
            'dodgerblue2',
            'violet'))

datag$Taxon <- factor(datag$Taxon, levels = c("Asterionellopsis","Chaetocerotaceae","Cylindrotheca","Leptocylindrus"
                                              ,"Nitzschia","Pseudo-nitzschia", "Skeletonema","Akashiwo","Azadinium",
                                              "Chrysochromulina","Phaeocystis","Cryptophyceae","Others"))

# Now with the proportion of each taxon in the community
fraction_plot <- ggplot(datag) +
  geom_col(aes(x = `Year-Month`, y = `Proportion`, 
               fill = Taxon), position = "stack", na.rm = FALSE, width = 1) +
  facet_wrap_color(~region, scales = "free_y", ncol = 1, colors = df_color) +
  scale_x_discrete(breaks = c('2007-06', '2008-06', '2009-06', '2010-06', 
                              '2011-06', '2012-06', '2013-06', '2014-06',
                              '2015-06', '2016-06', '2017-06', '2018-06',
                              '2019-06', '2020-06', '2021-06', '2022-06'),
                   labels = c('2007', '2008', '2009', '2010', 
                              '2011', '2012', '2013', '2014',
                              '2015', '2016', '2017', '2018',
                              '2019', '2020', '2021', '2022')
  ) +
  geom_vline(data = subset(datag, format(Date, "%m") == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", linewidth = 0.5,linetype = "dashed",
             show.legend = FALSE)+
  theme_classic() +
  theme(# Axes
    axis.text.y = element_text(),
    axis.text.x = element_text(),
    axis.ticks.x = element_blank(),
    # Legend
    legend.position = "right",
    legend.background = element_rect(color = 'grey15', linewidth = .45,
                                     fill = '#F8F7F5'),
    # Plot background
    panel.background = element_rect(fill = '#F8F7F5'),
    plot.background = element_rect(fill = '#F8F7F5'),
    # strip
    strip.background = element_rect(fill = 'grey40'),
    strip.text = element_text(size = 12, color = '#F8F7F5')
  )+
  # Labels and colors
  labs(x = 'Year', y = NULL, fill = 'Taxon: '
  )+
  scale_fill_manual(values = c(
    "Asterionellopsis"   = "#2B4561",
    "Chaetocerotaceae"   = "#76A7E2",
    "Cylindrotheca"      = "#2E6CD9",
    "Leptocylindrus"     = "#000E53",
    "Nitzschia"          = "#377185",
    "Pseudo-nitzschia"   = "#BBD4F2",
    "Skeletonema"        = "blue",
    "Akashiwo"           = "#B2DF8A",
    "Azadinium"          = "chartreuse4",
    "Chrysochromulina"   = "#FBB646",
    "Phaeocystis"        = "gold1",
    "Cryptophyceae"      = "#FC4D6B",
    "Others"             = "grey"
  ))

# Final graph:
fraction_plot
plotly_fraction <- ggplotly(fraction_plot)
saveWidget(plotly_fraction, file="output/HTML_figures/rephyto_plotly_fraction.html")

# Good stuff!

####------------------------------End of script-----------------------------####