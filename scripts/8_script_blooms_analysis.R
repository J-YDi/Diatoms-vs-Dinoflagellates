################################################################################
# Diatoms vs dinoflagellates: a network analysis of bloom impacts on diversity #
#                    and phytoplankton associations | R scripts                #
################################################################################

# Script to analyze the effet of a bloom on the different variables #
# 02/28/2025

# Load packages
library(readr)
library(dplyr)
library(ggforce)
library(ggplot2)
library(cowplot)
library(tidyr)
library(DescTools)

# Difference on percentage of dominance (Berger-Parker index) during between bloom and no bloom ####

# Import data
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_blooms_caracterised_final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)
# Adding the season info
data <- data |>
  mutate(season = case_when(Month %in% c(12, 01, 02) ~ "Winter",
                            Month %in% c(03, 04, 05) ~ "Spring",
                            Month %in% c(06, 07, 08) ~ "Summer",
                            Month %in% c(09, 10, 11) ~ "Fall", TRUE ~ NA_character_))

# Pick the same proportion of observations by season and by region as the blooms to compare blooms and non-blooms 
set.seed(123)
nbhiverclus1 <- sample_n(filter(data,season == "Winter" & region == "1-Mediterranean sea" & is.na(Bloom_Phylum)),21)
nbhiverclus2 <- sample_n(filter(data,season == "Winter" & region == "2-Eastern Channel - North Sea" & is.na(Bloom_Phylum)),3)
nbhiverclus3 <- sample_n(filter(data,season == "Winter" & region == "3-Atlantic - Western Channel" & is.na(Bloom_Phylum)),3)

nbautomneclus1 <- sample_n(filter(data,season == "Fall" & region == "1-Mediterranean sea" & is.na(Bloom_Phylum)),28)
nbautomneclus2 <- sample_n(filter(data,season == "Fall" & region == "2-Eastern Channel - North Sea" & is.na(Bloom_Phylum)),9)
nbautomneclus3 <- sample_n(filter(data,season == "Fall" & region == "3-Atlantic - Western Channel" & is.na(Bloom_Phylum)),11)

nbspringclus1 <- sample_n(filter(data,season == "Spring" & region == "1-Mediterranean sea" & is.na(Bloom_Phylum)),29)
nbspringclus2 <- sample_n(filter(data,season == "Spring" & region == "2-Eastern Channel - North Sea" & is.na(Bloom_Phylum)),55)
nbspringclus3 <- sample_n(filter(data,season == "Spring" & region == "3-Atlantic - Western Channel" & is.na(Bloom_Phylum)),62)

nbeteclus1 <- sample_n(filter(data,season == "Summer" & region == "1-Mediterranean sea" & is.na(Bloom_Phylum)),22)
nbeteclus2 <- sample_n(filter(data,season == "Summer" & region == "2-Eastern Channel - North Sea" & is.na(Bloom_Phylum)),39)
nbeteclus3 <- sample_n(filter(data,season == "Summer" & region == "3-Atlantic - Western Channel" & is.na(Bloom_Phylum)),25)

# Bind the random dates
data_hasard <- bind_rows(nbhiverclus1,nbhiverclus2,nbhiverclus3,nbautomneclus1,nbautomneclus2,nbautomneclus3,nbspringclus1,nbspringclus2,nbspringclus3,
                         nbeteclus1,nbeteclus2,nbeteclus3)

# Keep only the blooms and the random dates
data_bloom <- filter(data, !is.na(Bloom_Phylum))

data <- bind_rows(data_bloom,data_hasard)

data$TBloom <- ifelse(is.na(data$Bloom),"Non-bloom","Bloom")

data$region <- as.factor(data$region)

# Color for each region
region_col <- c("1-Mediterranean sea" = "#F8766D","2-Eastern Channel - North Sea" = "#CD9600", 
                "3-Atlantic - Western Channel" = "#00BE67",  "4-Pertuis Sea" = "#00A9FF")

d <- filter(data, region != "4-Pertuis Sea")
d$info <- "All regions"
BKglobal <- ggplot(d)+
  geom_violin(aes(x=TBloom,y=BergerParker,group=TBloom,colour=as.character(region)),size=1,draw_quantiles = c(0.5),show.legend = FALSE,fill = "gray91")+
  geom_sina(aes(x=TBloom,y=BergerParker,group=TBloom,colour = as.character(region)),alpha = .55,show.legend = F)+
  labs(x=NULL,y="Berger-Parker index")+
  scale_y_continuous(limits=c(0,1.10),breaks = c(0,0.33,0.66,1))+
  scale_colour_manual(values = region_col)+
  facet_wrap(~info)+
  theme_bw()

BKregion <- ggplot(filter(data, region !="4-Pertuis Sea"))+
  geom_violin(aes(x=TBloom,y=BergerParker,group=TBloom,colour=as.character(region)),size=1,draw_quantiles = c(0.5),show.legend = FALSE,fill = "gray91")+
  labs(x=NULL,y=NULL)+
  geom_sina(aes(x=TBloom,y=BergerParker,group=TBloom,colour=as.character(region)),alpha = .55,show.legend = F)+
  scale_colour_manual(values=region_col,guide = "none")+
  scale_y_continuous(limits=c(0,1.10),breaks = c(0,0.33,0.66,1))+
  facet_wrap(~region)+
  theme_bw()

plot_grid(BKglobal,BKregion,ncol=2,labels = "AUTO",rel_widths = c(1,2))
#ggsave('BergerParker_bloom_V2.png', path = "output/graphs/bloom", dpi = 600, width = 220, height = 100, units = 'mm')


# Statistical tests
# All blooms all regions
wilcox.test(filter(data,TBloom == "Bloom")$BergerParker,filter(data,TBloom != "Bloom")$BergerParker)

# All blooms by regions
wilcox.test(filter(data,TBloom == "Bloom" & region == "1-Mediterranean sea")$BergerParker,filter(data,TBloom != "Bloom" & region == "1-Mediterranean sea")$BergerParker)
wilcox.test(filter(data,TBloom == "Bloom" & region == "2-Eastern Channel - North Sea")$BergerParker,filter(data,TBloom != "Bloom" & region == "2-Eastern Channel - North Sea")$BergerParker)
wilcox.test(filter(data,TBloom == "Bloom" & region == "3-Atlantic - Western Channel")$BergerParker,filter(data,TBloom != "Bloom" & region == "3-Atlantic - Western Channel")$BergerParker)


# Impact of blooms we consider only dino, diatoms and Eastern Channel's Haptophytes blooms #####

# Import data
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_PCA_coord_moments.final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Keep only the bloom's dates (and before, after)
blooms <- filter(data, Bloom_Phylum == "Bac" | Bloom_Phylum == "Dino" | (region == "2-Eastern Channel - North Sea" & Bloom_Phylum == "Hapt"))
blooms <- filter(blooms, Moment != "Before -1")

# Color for each region
region_col <- c("1-Mediterranean sea" = "#F8766D","2-Eastern Channel - North Sea" = "#CD9600", 
                "3-Atlantic - Western Channel" = "#00BE67",  "4-Pertuis Sea" = "#00A9FF"
                ,"Bac" = "#1f77b4","Dino" = "green")

blooms$Moment <-
  factor(blooms$Moment,
         levels = c("Before","During","After"))

# Visualizations 
# Diversity indices
datag <- pivot_longer(blooms, cols = c(Shannon, Pielou, BergerParker),names_to = "Indice")

datag[datag$Indice == "BergerParker",]$Indice <- "Berger-Parker"

div <-ggplot(datag) +
  
  # Violin plot 
  geom_violin(aes(x = Moment, y = value, group = Moment), 
              size = 0.7, draw_quantiles = c(0.5), show.legend = FALSE, fill = "gray96") +
  
  # Mean and standard error
  stat_summary(aes(x = Moment, y = value, colour = Bloom_Phylum), 
               fun = mean, geom = "point", size = 5, alpha = 0.6, show.legend = FALSE) +
  stat_summary(aes(x = Moment, y = value, colour = Bloom_Phylum), 
               fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.2, size = 0.8, show.legend = FALSE) +
  
  #geom_sina(aes(x = Moment, y = value, group = Moment, colour = Bloom_Phylum), 
  #          alpha = 0, show.legend = TRUE) +
  
  labs(x = "Moment", y = "Value") +
  scale_colour_manual(values = region_col, guide = "none") +
  
  facet_grid(rows = vars(Indice), cols = vars(region), scales = "free_y") +
  theme_bw()
div
#ggsave('OriginalDiversityindex_moments_combine_3.png', path = "output/graphs/bloom", dpi = 600, width = 200, height = 200, units = 'mm')

# ACP dimensions
datag <- pivot_longer(blooms, cols = c(Dim.1,Dim.2,Dim.3),names_to = "Indice")

# Create the graph step by step to avoid a scale issue
Dim1 <- ggplot(filter(datag,Indice == "Dim.1")) +
  
  geom_violin(aes(x = Moment, y = value, group = Moment), 
              size = 0.7, draw_quantiles = c(0.5), show.legend = FALSE, fill = "gray96") +
  
  stat_summary(aes(x = Moment, y = value, colour = Bloom_Phylum), 
               fun = mean, geom = "point", size = 5, alpha = 0.6, show.legend = FALSE) +
  stat_summary(aes(x = Moment, y = value, colour = Bloom_Phylum), 
               fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.2, size = 0.8, show.legend = FALSE) +

  labs(x = NULL, y = NULL) +
  scale_colour_manual(values = region_col, guide = "none") +

  facet_grid(rows = vars(Indice), cols = vars(region), scales = "free_y")+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

Dim2 <- ggplot(filter(datag,Indice == "Dim.2")) +
  
  geom_violin(aes(x = Moment, y = value, group = Moment), 
              size = 0.7, draw_quantiles = c(0.5), show.legend = FALSE, fill = "gray96") +

  stat_summary(aes(x = Moment, y = value, colour = Bloom_Phylum), 
               fun = mean, geom = "point", size = 5, alpha = 0.6, show.legend = FALSE) +
  stat_summary(aes(x = Moment, y = value, colour = Bloom_Phylum), 
               fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.2, size = 0.8, show.legend = FALSE) +
  
  labs(x = NULL, y = NULL) +
  scale_colour_manual(values = region_col, guide = "none") +
  
  facet_grid(rows = vars(Indice), cols = vars(region), scales = "free_y")+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank())

Dim3 <- ggplot(filter(datag,Indice == "Dim.3")) +
  
  # Violin plot 
  geom_violin(aes(x = Moment, y = value, group = Moment), 
              size = 0.7, draw_quantiles = c(0.5), show.legend = FALSE, fill = "gray96") +
  
  stat_summary(aes(x = Moment, y = value, colour = Bloom_Phylum), 
               fun = mean, geom = "point", size = 5, alpha = 0.6, show.legend = FALSE) +
  stat_summary(aes(x = Moment, y = value, colour = Bloom_Phylum), 
               fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.2, size = 0.8, show.legend = FALSE) +
  
  labs(x = "Moment", y =NULL) +
  scale_colour_manual(values = region_col, guide = "none") +

  facet_grid(rows = vars(Indice), cols = vars(region), scales = "free_y")+
 scale_y_continuous(limits = c(-3.5,3.5))+
  theme_bw()+
  theme(strip.background.x = element_blank(),
                   strip.text.x = element_blank())

ACP <- plot_grid(Dim1,Dim2,Dim3,ncol=1)
ACP
#ggsave('ACPdim_moments_combine4.png', path = "output/graphs/bloom", dpi = 600, width = 200, height = 200, units = 'mm')

plot_grid(div,ACP,labels = "AUTO",nrow = 1)
#ggsave('ACPdim_div_moments_combine.png', path = "output/graphs/bloom", dpi = 600, width = 310, height = 190, units = 'mm')

### Statistical tests on the dino and diatoms blooms
# Mediterranean sea
blooms_med <- filter(blooms, region == "1-Mediterranean sea")

kruskal.test(blooms_med$Shannon~blooms_med$Moment)
kruskal.test(blooms_med$Pielou~blooms_med$Moment)
DunnTest(blooms_med$Pielou~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$BergerParker~blooms_med$Moment)

kruskal.test(blooms_med$Dim.1~blooms_med$Moment)
DunnTest(blooms_med$Dim.1~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$Dim.2~blooms_med$Moment)
kruskal.test(blooms_med$Dim.3~blooms_med$Moment)

kruskal.test(blooms_med$P_bac~blooms_med$Moment)
kruskal.test(blooms_med$P_dino~blooms_med$Moment)
kruskal.test(blooms_med$P_autres~blooms_med$Moment)

kruskal.test(blooms_med$P_BacBac~blooms_med$Moment)
kruskal.test(blooms_med$P_BacDino~blooms_med$Moment)
kruskal.test(blooms_med$P_DinoDino~blooms_med$Moment)

# Eastern channel
blooms_manche <- filter(blooms, region == "2-Eastern Channel - North Sea")

kruskal.test(blooms_manche$Shannon~blooms_manche$Moment)
DunnTest(blooms_manche$Shannon~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$Pielou~blooms_manche$Moment)
DunnTest(blooms_manche$Pielou~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$BergerParker~blooms_manche$Moment)
DunnTest(blooms_manche$BergerParker~blooms_manche$Moment,method = "BH")

kruskal.test(blooms_manche$Dim.1~blooms_manche$Moment)
DunnTest(blooms_manche$Dim.1~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$Dim.2~blooms_manche$Moment)
kruskal.test(blooms_manche$Dim.3~blooms_manche$Moment)

kruskal.test(blooms_manche$P_bac~blooms_manche$Moment)
DunnTest(blooms_manche$P_bac~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_dino~blooms_manche$Moment)
DunnTest(blooms_manche$P_dino~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_autres~blooms_manche$Moment)

kruskal.test(blooms_manche$P_BacBac~blooms_manche$Moment)
DunnTest(blooms_manche$P_BacBac~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_BacDino~blooms_manche$Moment)
DunnTest(blooms_manche$P_BacDino~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_DinoDino~blooms_manche$Moment)
DunnTest(blooms_manche$P_DinoDino~blooms_manche$Moment,method = "BH")


# Atlantic ocean
blooms_atlantique <- filter(blooms, region == "3-Atlantic - Western Channel")


kruskal.test(blooms_atlantique$Shannon~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$Shannon~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$Pielou~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$Pielou~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$BergerParker~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$BergerParker~blooms_atlantique$Moment,method = "BH")

kruskal.test(blooms_atlantique$Dim.1~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$Dim.1~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$Dim.2~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$Dim.3~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$Dim.3~blooms_atlantique$Moment,method = "BH")

kruskal.test(blooms_atlantique$P_bac~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$P_dino~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$P_autres~blooms_atlantique$Moment)

kruskal.test(blooms_atlantique$P_BacBac~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$P_BacDino~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$P_DinoDino~blooms_atlantique$Moment)



### Statistical tests on the dinoflagellates blooms
# Mediterranean sea
blooms_med <- filter(blooms, region == "1-Mediterranean sea" & Bloom_Phylum == "Dino")

kruskal.test(blooms_med$Shannon~blooms_med$Moment)
kruskal.test(blooms_med$Pielou~blooms_med$Moment)
kruskal.test(blooms_med$BergerParker~blooms_med$Moment)

kruskal.test(blooms_med$Dim.1~blooms_med$Moment)
kruskal.test(blooms_med$Dim.2~blooms_med$Moment)
kruskal.test(blooms_med$Dim.3~blooms_med$Moment)

kruskal.test(blooms_med$P_bac~blooms_med$Moment)
kruskal.test(blooms_med$P_dino~blooms_med$Moment)
kruskal.test(blooms_med$P_autres~blooms_med$Moment)

kruskal.test(blooms_med$P_BacBac~blooms_med$Moment)
kruskal.test(blooms_med$P_BacDino~blooms_med$Moment)
kruskal.test(blooms_med$P_DinoDino~blooms_med$Moment)
kruskal.test(blooms_med$P_AAutres~blooms_med$Moment)


# Eastern channel
blooms_manche <- filter(blooms, region == "2-Eastern Channel - North Sea" & Bloom_Phylum == "Dino")

kruskal.test(blooms_manche$Shannon~blooms_manche$Moment)
kruskal.test(blooms_manche$Pielou~blooms_manche$Moment)
kruskal.test(blooms_manche$BergerParker~blooms_manche$Moment)

kruskal.test(blooms_manche$Dim.1~blooms_manche$Moment)
kruskal.test(blooms_manche$Dim.2~blooms_manche$Moment)
kruskal.test(blooms_manche$Dim.3~blooms_manche$Moment)

kruskal.test(blooms_manche$P_bac~blooms_manche$Moment)
kruskal.test(blooms_manche$P_dino~blooms_manche$Moment)
kruskal.test(blooms_manche$P_autres~blooms_manche$Moment)

kruskal.test(blooms_manche$P_BacBac~blooms_manche$Moment)
kruskal.test(blooms_manche$P_BacDino~blooms_manche$Moment)
DunnTest(blooms_manche$P_BacDino~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_DinoDino~blooms_manche$Moment)
kruskal.test(blooms_manche$P_AAutres~blooms_manche$Moment)
DunnTest(blooms_manche$P_AAutres~blooms_manche$Moment,method = "BH")


# Atlantic ocean
blooms_atlantique <- filter(blooms, region == "3-Atlantic - Western Channel" & Bloom_Phylum == "Dino")

kruskal.test(blooms_atlantique$Shannon~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$Pielou~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$BergerParker~blooms_atlantique$Moment)

kruskal.test(blooms_atlantique$Dim.1~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$Dim.2~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$Dim.3~blooms_atlantique$Moment)

kruskal.test(blooms_atlantique$P_bac~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$P_dino~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$P_autres~blooms_atlantique$Moment)

kruskal.test(blooms_atlantique$P_BacBac~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$P_BacDino~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$P_DinoDino~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$P_AAutres~blooms_atlantique$Moment)


### Statistical tests on the diatoms blooms
# Mediterranean sea
blooms_med <- filter(blooms, region == "1-Mediterranean sea" & Bloom_Phylum == "Bac")

kruskal.test(blooms_med$Shannon~blooms_med$Moment)
kruskal.test(blooms_med$Pielou~blooms_med$Moment)
DunnTest(blooms_med$Pielou~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$BergerParker~blooms_med$Moment)

kruskal.test(blooms_med$Dim.1~blooms_med$Moment)
DunnTest(blooms_med$Dim.1~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$Dim.2~blooms_med$Moment)
kruskal.test(blooms_med$Dim.3~blooms_med$Moment)

kruskal.test(blooms_med$P_bac~blooms_med$Moment)
kruskal.test(blooms_med$P_dino~blooms_med$Moment)
kruskal.test(blooms_med$P_autres~blooms_med$Moment)

kruskal.test(blooms_med$P_BacBac~blooms_med$Moment)
kruskal.test(blooms_med$P_BacDino~blooms_med$Moment)
kruskal.test(blooms_med$P_DinoDino~blooms_med$Moment)
kruskal.test(blooms_med$P_AAutres~blooms_med$Moment)

# Eastern channel
blooms_manche <- filter(blooms, region == "2-Eastern Channel - North Sea" & Bloom_Phylum == "Bac")

kruskal.test(blooms_manche$Shannon~blooms_manche$Moment)
DunnTest(blooms_manche$Shannon~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$Pielou~blooms_manche$Moment)
DunnTest(blooms_manche$Pielou~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$BergerParker~blooms_manche$Moment)
DunnTest(blooms_manche$BergerParker~blooms_manche$Moment,method = "BH")

kruskal.test(blooms_manche$Dim.1~blooms_manche$Moment)
DunnTest(blooms_manche$Dim.1~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$Dim.2~blooms_manche$Moment)
kruskal.test(blooms_manche$Dim.3~blooms_manche$Moment)

kruskal.test(blooms_manche$P_bac~blooms_manche$Moment)
DunnTest(blooms_manche$P_bac~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_dino~blooms_manche$Moment)
DunnTest(blooms_manche$P_dino~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_autres~blooms_manche$Moment)

kruskal.test(blooms_manche$P_BacBac~blooms_manche$Moment)
DunnTest(blooms_manche$P_BacBac~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_BacDino~blooms_manche$Moment)
kruskal.test(blooms_manche$P_DinoDino~blooms_manche$Moment)
DunnTest(blooms_manche$P_DinoDino~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_AAutres~blooms_manche$Moment)

# Atlantic ocean
blooms_atlantique <- filter(blooms, region == "3-Atlantic - Western Channel" & Bloom_Phylum == "Bac")

kruskal.test(blooms_atlantique$Shannon~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$Shannon~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$Pielou~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$Pielou~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$BergerParker~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$BergerParker~blooms_atlantique$Moment,method = "BH")

kruskal.test(blooms_atlantique$Dim.1~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$Dim.1~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$Dim.2~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$Dim.3~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$Dim.3~blooms_atlantique$Moment,method = "BH")

kruskal.test(blooms_atlantique$P_bac~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$P_dino~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$P_autres~blooms_atlantique$Moment)

kruskal.test(blooms_atlantique$P_BacBac~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$P_BacDino~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$P_BacDino~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$P_DinoDino~blooms_atlantique$Moment)
kruskal.test(blooms_atlantique$P_AAutres~blooms_atlantique$Moment)

# Eastern channel and Haptophyte's blooms
blooms_manche <- filter(blooms, region == "2-Eastern Channel - North Sea" & Bloom_Phylum == "Hapt")

kruskal.test(blooms_manche$Shannon~blooms_manche$Moment)
kruskal.test(blooms_manche$Pielou~blooms_manche$Moment)
kruskal.test(blooms_manche$BergerParker~blooms_manche$Moment)

kruskal.test(blooms_manche$Dim.1~blooms_manche$Moment)
kruskal.test(blooms_manche$Dim.2~blooms_manche$Moment)
kruskal.test(blooms_manche$Dim.3~blooms_manche$Moment)

kruskal.test(blooms_manche$P_bac~blooms_manche$Moment)
DunnTest(blooms_manche$P_bac~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_dino~blooms_manche$Moment)
kruskal.test(blooms_manche$P_autres~blooms_manche$Moment)

kruskal.test(blooms_manche$P_BacBac~blooms_manche$Moment)
kruskal.test(blooms_manche$P_DinoDino~blooms_manche$Moment)
kruskal.test(blooms_manche$P_BacDino~blooms_manche$Moment)
kruskal.test(blooms_manche$P_AAutres~blooms_manche$Moment)


### Calculate the mean for each variable 
### Diatoms and dinoflagellates blooms

# Import data
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_PCA_coord_moments.final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Keep only the bloom's dates (and before, after)
blooms <- filter(data, Bloom_Phylum == "Bac" | Bloom_Phylum == "Dino")

blooms <- select(blooms, Date, Code_point_Libelle, region, Moment,Bloom_Phylum,RShannon,RPielou, RBergerParker, 
                 Shannon, Pielou, BergerParker,Dim.1,Dim.2,Dim.3,P_bac,P_dino,P_autres,
                 P_AAutres,P_BacBac,P_BacDino,P_DinoDino,
                 PP_AAutres,PP_BacBac,PP_BacDino,PP_DinoDino)

blooms_longer <- pivot_longer(blooms, cols = c(RShannon:PP_DinoDino))

blooms_summary <- summarise(group_by(blooms_longer,region,Moment,name),value=mean(value,na.rm=T))

blooms_summary <- pivot_wider(blooms_summary,names_from = name,values_fn = mean)
# Save it
# write.csv2(blooms_summary,file="output/tableaux/Blooms/mean_variable_blooms.csv", row.names = FALSE,dec = ".")


# Figure Impact of blooms we consider only dino, diatoms and Eastern Channel's Haptophytes blooms on association types and composition #####
# Import data
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_PCA_coord_moments.final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Keep only the bloom's dates (and before, after)
blooms <- filter(data, Bloom_Phylum == "Bac" | Bloom_Phylum == "Dino" | (region == "2-Eastern Channel - North Sea" & Bloom_Phylum == "Hapt"))
blooms <- filter(blooms, Moment != "Before -1")
blooms[blooms$Bloom_Phylum == "Bac",]$Bloom_Phylum <- "Bacillariophyceae"
blooms[blooms$Bloom_Phylum == "Dino",]$Bloom_Phylum <- "Dinophyceae"
blooms[blooms$Bloom_Phylum == "Hapt",]$Bloom_Phylum <- "Haptophyta"

# Color for each region
region_col <- c("1-Mediterranean sea" = "#F8766D","2-Eastern Channel - North Sea" = "#CD9600", 
                "3-Atlantic - Western Channel" = "#00BE67",  "4-Pertuis Sea" = "#00A9FF"
                ,"Bac" = "#1f77b4","Dino" = "green")

blooms$Moment <-
  factor(blooms$Moment,
         levels = c("Before","During","After"))
blooms <- select(blooms, region,P_bac,P_dino,P_autres,P_BacBac,P_BacDino,P_DinoDino,P_AAutres,Moment,Bloom_Phylum)

colnames(blooms) <- c("region","Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso.","Moment","Bloom_Phylum")

blooms_longer <- pivot_longer(blooms, cols = c(Bac:`Other asso.`))

blooms_summary <- summarise(group_by(blooms_longer,region,Moment,name,Bloom_Phylum),Val=mean(value,na.rm=T),sd=sd(value,na.rm = T))

blooms_summary_compo <- blooms_summary %>%
  filter(name %in% c("Bac", "Dino", "Other taxa")) %>%
  group_by(Moment, region,Bloom_Phylum) %>%
  mutate(cumsum_Val = cumsum(Val)) 

blooms_summary_asso <- blooms_summary %>%
  filter(name %in% c("Bac-Bac", "Bac-Dino","Dino-Dino", "Other asso.")) %>%
  group_by(Moment, region,Bloom_Phylum) %>%
  mutate(cumsum_Val = cumsum(Val))



compo_bac <- ggplot(filter(blooms_summary_compo, name %in% c("Bac", "Dino","Other taxa") & Bloom_Phylum=="Bacillariophyceae"), aes(Moment, Val)) +
  geom_col(aes(fill = name), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = Val, ymax = Val+sd, group = name,colour=name),
    width = 0.2, position = position_dodge(0.8)
  )+
  scale_fill_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac" = "#56B4E9","Bac-Dino"   = "chocolate1"
    ,"Other taxa"      = "grey" ,"Dino-Dino" = "#009E73","Other asso." = "grey30"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL) +
  scale_colour_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac" = "#56B4E9","Bac-Dino"   = "chocolate1"
    ,"Other taxa"      = "grey" ,"Dino-Dino" = "#009E73","Other asso." = "grey30"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL)+
  facet_wrap(Bloom_Phylum~region)+
  labs(x=NULL,y="Proportion")+
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank())+
  scale_y_continuous(limits = c(0,0.8))

asso_bac <-  ggplot(filter(blooms_summary_asso, name %in% c("Bac-Bac","Bac-Dino", "Dino-Dino", "Other asso.") & Bloom_Phylum=="Bacillariophyceae"), aes(Moment, Val)) +
  geom_col(aes(fill = name), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = Val, ymax = Val+sd, group = name,colour=name),
    width = 0.2, position = position_dodge(0.8)
  )+
  scale_fill_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac" = "#56B4E9","Bac-Dino"   = "chocolate1"
    ,"Other taxa"      = "grey" ,"Dino-Dino" = "#009E73","Other asso." = "grey30"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL) +
  scale_colour_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac" = "#56B4E9","Bac-Dino"   = "chocolate1"
    ,"Other taxa"      = "grey" ,"Dino-Dino" = "#009E73","Other asso." = "grey30"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL)+
  facet_wrap(~region)+
  labs(x=NULL,y="Proportion")+
  theme(legend.position = "none")+
  scale_y_continuous(limits = c(0,0.65))+
  theme(strip.text = element_blank())

bac <- plot_grid(compo_bac,asso_bac,nrow = 2)


compo_dino <- ggplot(filter(blooms_summary_compo, name %in% c("Bac", "Dino", "Other taxa") & Bloom_Phylum=="Dinophyceae"), aes(Moment, Val)) +
  geom_col(aes(fill = name), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = Val, ymax = Val+sd, group = name,colour=name),
    width = 0.2, position = position_dodge(0.8)
  )+
  scale_fill_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac" = "#56B4E9","Bac-Dino"   = "chocolate1"
    ,"Other taxa"      = "grey" ,"Dino-Dino" = "#009E73","Other asso." = "grey30"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL) +
  scale_colour_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac" = "#56B4E9","Bac-Dino"   = "chocolate1"
    ,"Other taxa"      = "grey" ,"Dino-Dino" = "#009E73","Other asso." = "grey30"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL)+
  facet_wrap(Bloom_Phylum~region)+
  labs(x=NULL,y=NULL)+
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())+
  scale_y_continuous(limits = c(0,0.8))

asso_dino <-  ggplot(filter(blooms_summary_asso, name %in% c("Bac-Bac", "Dino-Dino", "Bac-Dino", "Other asso.") & Bloom_Phylum=="Dinophyceae"), aes(Moment, Val)) +
  geom_col(aes(fill = name), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = Val, ymax = Val+sd, group = name,colour=name),
    width = 0.2, position = position_dodge(0.8)
  )+
  scale_fill_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac" = "#56B4E9","Bac-Dino"   = "chocolate1"
    ,"Other taxa"      = "grey" ,"Dino-Dino" = "#009E73","Other asso." = "grey30"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL) +
  scale_colour_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac" = "#56B4E9","Bac-Dino"   = "chocolate1"
    ,"Other taxa"      = "grey" ,"Dino-Dino" = "#009E73","Other asso." = "grey30"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL)+
  facet_wrap(~region)+
  labs(x=NULL,y=NULL)+
  theme(legend.position = "none",axis.text.y = element_blank())+
  scale_y_continuous(limits = c(0,0.65))+
  theme(strip.text = element_blank())

dino <- plot_grid(compo_dino,asso_dino,nrow = 2)

plot_grid(bac,dino,nrow = 1)

compo_hapt <- ggplot(filter(blooms_summary_compo, name %in% c("Bac", "Dino", "Other taxa") & Bloom_Phylum=="Haptophyta"), aes(Moment, Val)) +
  geom_col(aes(fill = name), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = Val, ymax = Val+sd, group = name,colour=name),
    width = 0.2, position = position_dodge(0.8)
  )+
  scale_fill_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac" = "#56B4E9","Bac-Dino"   = "chocolate1"
    ,"Other taxa"      = "grey" ,"Dino-Dino" = "#009E73","Other asso." = "grey30"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL) +
  scale_colour_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac" = "#56B4E9","Bac-Dino"   = "chocolate1"
    ,"Other taxa"      = "grey" ,"Dino-Dino" = "#009E73","Other asso." = "grey30"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL)+
  facet_wrap(Bloom_Phylum~region)+
  labs(x=NULL,y=NULL)+
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())+
  scale_y_continuous(limits = c(0,0.8))

asso_hapt <-  ggplot(filter(blooms_summary_asso, name %in% c("Bac-Bac", "Dino-Dino","Bac-Dino", "Other asso.") & Bloom_Phylum=="Haptophyta"), aes(Moment, Val)) +
  geom_col(aes(fill = name), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = Val, ymax = Val+sd, group = name,colour=name),
    width = 0.2, position = position_dodge(0.8)
  )+
  scale_fill_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac" = "#56B4E9","Bac-Dino"   = "chocolate1"
    ,"Other taxa"      = "grey" ,"Dino-Dino" = "#009E73","Other asso." = "grey30"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL) +
  scale_colour_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac" = "#56B4E9","Bac-Dino"   = "chocolate1"
    ,"Other taxa"      = "grey" ,"Dino-Dino" = "#009E73","Other asso." = "grey30"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL)+
  facet_wrap(~region)+
  labs(x=NULL,y=NULL)+
  theme(legend.position = "none")+
  scale_y_continuous(limits = c(0,0.65))+
  theme(strip.text = element_blank(),axis.text.y = element_blank())
hapt <- plot_grid(compo_hapt,asso_hapt,nrow = 2)

plot_grid(bac,dino,hapt,nrow=1,rel_widths = c(1,1,0.3))
#ggsave('asso_compo_dynam_blooms.png', path = "output/graphs/bloom", dpi = 600, width = 360, height = 170, units = 'mm')

