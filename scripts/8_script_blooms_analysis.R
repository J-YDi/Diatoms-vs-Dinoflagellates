################################################################################
# Diatoms vs dinoflagellates: a network analysis of bloom impacts on diversity #
#                    and phytoplankton associations | R scripts                #
################################################################################

# Script to analyze the effet of a bloom on the different variables #
# 02/26/2026

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
  geom_violin(aes(x=TBloom,y=BergerParker,group=TBloom,colour=CHLOROA),size=1,draw_quantiles = c(0.5),show.legend = FALSE,fill = "gray91")+
  geom_sina(aes(x=TBloom,y=BergerParker,group=TBloom,colour = CHLOROA),alpha = .55,show.legend = F)+
  labs(x=NULL,y="Berger-Parker index")+
  scale_y_continuous(limits=c(0,1.10),breaks = c(0,0.33,0.66,1))+
  scale_colour_viridis_c()+
  facet_wrap(~info)+
  theme_bw()

cor.test(filter(d,TBloom == "Bloom")$BergerParker,filter(d,TBloom == "Bloom")$CHLOROA,method = "spearman")
cor.test(filter(d,TBloom == "Non-bloom")$BergerParker,filter(d,TBloom == "Non-bloom")$CHLOROA,method = "spearman")

# plot(filter(d,TBloom == "Bloom")$BergerParker,filter(d,TBloom == "Bloom")$CHLOROA)
# plot(filter(d,TBloom == "Non-bloom")$BergerParker,filter(d,TBloom == "Non-bloom")$CHLOROA)

BKregion <- ggplot(filter(data, region !="4-Pertuis Sea"))+
  geom_violin(aes(x=TBloom,y=BergerParker,group=TBloom,colour=CHLOROA),size=1,draw_quantiles = c(0.5),show.legend = FALSE,fill = "gray91")+
  labs(x=NULL,y=NULL)+
  geom_sina(aes(x=TBloom,y=BergerParker,group=TBloom,colour=CHLOROA),alpha = .55,show.legend = F)+
  scale_colour_viridis_c()+
  scale_y_continuous(limits=c(0,1.10),breaks = c(0,0.33,0.66,1))+
  facet_wrap(~region)+
  theme_bw()

cor.test(filter(d,region == "1-Mediterranean sea" & TBloom == "Bloom")$BergerParker,filter(d,region == "1-Mediterranean sea" & TBloom == "Bloom")$CHLOROA,method = "spearman")
cor.test(filter(d,region == "1-Mediterranean sea" & TBloom == "Non-bloom")$BergerParker,filter(d,region == "1-Mediterranean sea" & TBloom == "Non-bloom")$CHLOROA,method = "spearman")

# plot(filter(d,region == "1-Mediterranean sea" & TBloom == "Bloom")$BergerParker,filter(d,region == "1-Mediterranean sea" & TBloom == "Bloom")$CHLOROA)
# plot(filter(d,region == "1-Mediterranean sea" & TBloom == "Non-bloom")$BergerParker,filter(d,region == "1-Mediterranean sea" & TBloom == "Non-bloom")$CHLOROA)

cor.test(filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Bloom")$BergerParker,filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Bloom")$CHLOROA,method = "spearman")
cor.test(filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Non-bloom")$BergerParker,filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Non-bloom")$CHLOROA,method = "spearman")

# plot(filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Bloom")$BergerParker,filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Bloom")$CHLOROA)
# plot(filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Non-bloom")$BergerParker,filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Non-bloom")$CHLOROA)

cor.test(filter(d,region == "3-Atlantic - Western Channel" & TBloom == "Bloom")$BergerParker,filter(d,region == "3-Atlantic - Western Channel" & TBloom == "Bloom")$CHLOROA,method = "spearman")
cor.test(filter(d,region == "3-Atlantic - Western Channel" & TBloom == "Non-bloom")$BergerParker,filter(d,region == "3-Atlantic - Western Channel" & TBloom == "Non-bloom")$CHLOROA,method = "spearman")

# plot(filter(d,region == "3-Atlantic - Western Channel" & TBloom == "Bloom")$BergerParker,filter(d,region == "3-Atlantic - Western Channel" & TBloom == "Bloom")$CHLOROA)
# plot(filter(d,region == "3-Atlantic - Western Channel" & TBloom == "Non-bloom")$BergerParker,filter(d,region == "3-Atlantic - Western Channel" & TBloom == "Non-bloom")$CHLOROA)


data$Bloom_Phylum[is.na(data$Bloom_Phylum)] <- "Non-bloom"

BKdetail <- ggplot(filter(data, region != "4-Pertuis Sea", Bloom_Phylum == "Bac" | Bloom_Phylum == "Dino" | Bloom_Phylum == "Hapt" | Bloom_Phylum == "Non-bloom"))+
  geom_violin(aes(x=Bloom_Phylum,y=BergerParker,group=Bloom_Phylum,colour=as.character(region)),size=1,draw_quantiles = c(0.5),show.legend = FALSE,fill = "gray91")+
  labs(x="Type of bloom",y="")+
  scale_colour_manual(values=region_col,guide = "none")+
  scale_y_continuous(limits=c(0,1.10),breaks = c(0,0.33,0.66,1))+
  geom_sina(aes(x=Bloom_Phylum,y=BergerParker,group=Bloom_Phylum,colour = as.character(region)),alpha = .55,show.legend = F)+
  facet_wrap(~region,scales = "free_x")+
  theme_bw()

cor.test(filter(d,region == "1-Mediterranean sea" & TBloom == "Bloom" & Bloom_Phylum =="Bac")$BergerParker,filter(d,region == "1-Mediterranean sea" & TBloom == "Bloom" & Bloom_Phylum =="Bac")$CHLOROA,method = "spearman")
# plot(filter(d,region == "1-Mediterranean sea" & TBloom == "Bloom"& Bloom_Phylum =="Bac")$BergerParker,filter(d,region == "1-Mediterranean sea" & TBloom == "Bloom"& Bloom_Phylum =="Bac")$CHLOROA)

cor.test(filter(d,region == "1-Mediterranean sea" & TBloom == "Bloom" & Bloom_Phylum =="Dino")$BergerParker,filter(d,region == "1-Mediterranean sea" & TBloom == "Bloom" & Bloom_Phylum =="Dino")$CHLOROA,method = "spearman")
# plot(filter(d,region == "1-Mediterranean sea" & TBloom == "Bloom"& Bloom_Phylum =="Dino")$BergerParker,filter(d,region == "1-Mediterranean sea" & TBloom == "Bloom"& Bloom_Phylum =="Dino")$CHLOROA)


cor.test(filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Bloom"& Bloom_Phylum =="Bac")$BergerParker,filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Bloom"& Bloom_Phylum =="Bac")$CHLOROA,method = "spearman")
# plot(filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Bloom"& Bloom_Phylum =="Bac")$BergerParker,filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Bloom"& Bloom_Phylum =="Bac")$CHLOROA)

cor.test(filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Bloom"& Bloom_Phylum =="Dino")$BergerParker,filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Bloom"& Bloom_Phylum =="Dino")$CHLOROA,method = "spearman")
# plot(filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Bloom"& Bloom_Phylum =="Dino")$BergerParker,filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Bloom"& Bloom_Phylum =="Dino")$CHLOROA)

cor.test(filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Bloom"& Bloom_Phylum =="Hapt")$BergerParker,filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Bloom"& Bloom_Phylum =="Hapt")$CHLOROA,method = "spearman")
# plot(filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Bloom"& Bloom_Phylum =="Hapt")$BergerParker,filter(d,region == "2-Eastern Channel - North Sea" & TBloom == "Bloom"& Bloom_Phylum =="Hapt")$CHLOROA)

cor.test(filter(d,region == "3-Atlantic - Western Channel" & TBloom == "Bloom" & Bloom_Phylum =="Bac")$BergerParker,filter(d,region == "3-Atlantic - Western Channel" & TBloom == "Bloom" & Bloom_Phylum =="Bac")$CHLOROA,method = "spearman")
# plot(filter(d,region == "3-Atlantic - Western Channel" & TBloom == "Bloom"& Bloom_Phylum =="Bac")$BergerParker,filter(d,region == "3-Atlantic - Western Channel" & TBloom == "Bloom"& Bloom_Phylum =="Bac")$CHLOROA)

cor.test(filter(d,region == "3-Atlantic - Western Channel" & TBloom == "Bloom" & Bloom_Phylum =="Dino")$BergerParker,filter(d,region == "3-Atlantic - Western Channel" & TBloom == "Bloom" & Bloom_Phylum =="Dino")$CHLOROA,method = "spearman")
# plot(filter(d,region == "3-Atlantic - Western Channel" & TBloom == "Bloom"& Bloom_Phylum =="Dino")$BergerParker,filter(d,region == "3-Atlantic - Western Channel" & TBloom == "Bloom"& Bloom_Phylum =="Dino")$CHLOROA)





plot_grid(BKglobal,BKregion,BKdetail,ncol=1,labels = "AUTO")
#ggsave('BergerParker_bloom.png', path = "output/graphs/bloom", dpi = 600, width = 200, height = 200, units = 'mm')

plot_grid(BKglobal,BKregion,ncol=2,labels = "AUTO",rel_widths = c(1,2))
#ggsave('BergerParker_bloom_V2.png', path = "output/graphs/bloom", dpi = 600, width = 220, height = 100, units = 'mm')


# Statistical tests
# All blooms all regions
wilcox.test(filter(data,TBloom == "Bloom")$BergerParker,filter(data,TBloom != "Bloom")$BergerParker)

# All blooms by regions
wilcox.test(filter(data,TBloom == "Bloom" & region == "1-Mediterranean sea")$BergerParker,filter(data,TBloom != "Bloom" & region == "1-Mediterranean sea")$BergerParker)
wilcox.test(filter(data,TBloom == "Bloom" & region == "2-Eastern Channel - North Sea")$BergerParker,filter(data,TBloom != "Bloom" & region == "2-Eastern Channel - North Sea")$BergerParker)
wilcox.test(filter(data,TBloom == "Bloom" & region == "3-Atlantic - Western Channel")$BergerParker,filter(data,TBloom != "Bloom" & region == "3-Atlantic - Western Channel")$BergerParker)

# Dinophyceae blooms, by regions
wilcox.test(filter(data,Bloom_Phylum == "Dino" & region == "1-Mediterranean sea")$BergerParker,filter(data,TBloom != "Bloom" & region == "1-Mediterranean sea")$BergerParker)
wilcox.test(filter(data,Bloom_Phylum == "Dino" & region == "2-Eastern Channel - North Sea")$BergerParker,filter(data,TBloom != "Bloom" & region == "2-Eastern Channel - North Sea")$BergerParker)
wilcox.test(filter(data,Bloom_Phylum == "Dino" & region == "3-Atlantic - Western Channel")$BergerParker,filter(data,TBloom != "Bloom" & region == "3-Atlantic - Western Channel")$BergerParker)
# Diatoms blooms, by region
wilcox.test(filter(data,Bloom_Phylum == "Bac" & region == "1-Mediterranean sea")$BergerParker,filter(data,TBloom != "Bloom" & region == "1-Mediterranean sea")$BergerParker)
wilcox.test(filter(data,Bloom_Phylum == "Bac" & region == "2-Eastern Channel - North Sea")$BergerParker,filter(data,TBloom != "Bloom" & region == "2-Eastern Channel - North Sea")$BergerParker)
wilcox.test(filter(data,Bloom_Phylum == "Bac" & region == "3-Atlantic - Western Channel")$BergerParker,filter(data,TBloom != "Bloom" & region == "3-Atlantic - Western Channel")$BergerParker)
# Haptophytes by region
wilcox.test(filter(data,Bloom_Phylum == "Hapt" & region == "2-Eastern Channel - North Sea")$BergerParker,filter(data,TBloom != "Bloom" & region == "2-Eastern Channel - North Sea")$BergerParker)


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

datag <- pivot_longer(blooms, cols = c(RShannon, RPielou, RBergerParker),names_to = "Indice")
ggplot(datag) +
  
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

ggsave('Diversityindex_moments_combine_review.png', path = "output/graphs/bloom", dpi = 600, width = 200, height = 200, units = 'mm')

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
ggsave('OriginalDiversityindex_moments_combine_3_review.png', path = "output/graphs/bloom", dpi = 600, width = 200, height = 200, units = 'mm')

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

ggsave('ACPdim_moments_combine4_review.png', path = "output/graphs/bloom", dpi = 600, width = 200, height = 200, units = 'mm')

plot_grid(div,ACP,labels = "AUTO",nrow = 1)
ggsave('ACPdim_div_moments_combine_review.png', path = "output/graphs/bloom", dpi = 600, width = 310, height = 190, units = 'mm')

### Statistical tests on the dino and diatoms blooms
# Mediterranean sea
blooms_med <- filter(blooms, region == "1-Mediterranean sea")

kruskal.test(blooms_med$RShannon~blooms_med$Moment)
DunnTest(blooms_med$RShannon~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$RPielou~blooms_med$Moment)
DunnTest(blooms_med$RPielou~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$RBergerParker~blooms_med$Moment)
DunnTest(blooms_med$RBergerParker~blooms_med$Moment,method = "BH")

kruskal.test(blooms_med$Shannon~blooms_med$Moment)
DunnTest(blooms_med$Shannon~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$Pielou~blooms_med$Moment)
DunnTest(blooms_med$Pielou~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$BergerParker~blooms_med$Moment)
DunnTest(blooms_med$BergerParker~blooms_med$Moment,method = "BH")

kruskal.test(blooms_med$Dim.1~blooms_med$Moment)
DunnTest(blooms_med$Dim.1~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$Dim.2~blooms_med$Moment)
DunnTest(blooms_med$Dim.2~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$Dim.3~blooms_med$Moment)
DunnTest(blooms_med$Dim.3~blooms_med$Moment,method = "BH")

kruskal.test(blooms_med$P_bac~blooms_med$Moment)
DunnTest(blooms_med$P_bac~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$P_dino~blooms_med$Moment)
DunnTest(blooms_med$P_dino~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$P_autres~blooms_med$Moment)
DunnTest(blooms_med$P_bac~blooms_med$Moment,method = "BH")

kruskal.test(blooms_med$PP_BacBac~blooms_med$Moment)
DunnTest(blooms_med$PP_BacBac~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$PP_BacDino~blooms_med$Moment)
DunnTest(blooms_med$PP_BacDino~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$PP_DinoDino~blooms_med$Moment)
DunnTest(blooms_med$PP_DinoDino~blooms_med$Moment,method = "BH")

kruskal.test(blooms_med$P_BacBac~blooms_med$Moment)
DunnTest(blooms_med$P_BacBac~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$P_BacDino~blooms_med$Moment)
DunnTest(blooms_med$P_BacDino~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$P_DinoDino~blooms_med$Moment)
DunnTest(blooms_med$P_DinoDino~blooms_med$Moment,method = "BH")

# Eastern channel
blooms_manche <- filter(blooms, region == "2-Eastern Channel - North Sea")

kruskal.test(blooms_manche$RShannon~blooms_manche$Moment)
DunnTest(blooms_manche$RShannon~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$RPielou~blooms_manche$Moment)
DunnTest(blooms_manche$RPielou~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$RBergerParker~blooms_manche$Moment)
DunnTest(blooms_manche$RBergerParker~blooms_manche$Moment,method = "BH")

kruskal.test(blooms_manche$Shannon~blooms_manche$Moment)
DunnTest(blooms_manche$Shannon~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$Pielou~blooms_manche$Moment)
DunnTest(blooms_manche$Pielou~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$BergerParker~blooms_manche$Moment)
DunnTest(blooms_manche$BergerParker~blooms_manche$Moment,method = "BH")

kruskal.test(blooms_manche$Dim.1~blooms_manche$Moment)
DunnTest(blooms_manche$Dim.1~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$Dim.2~blooms_manche$Moment)
DunnTest(blooms_manche$Dim.2~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$Dim.3~blooms_manche$Moment)
DunnTest(blooms_manche$Dim.3~blooms_manche$Moment,method = "BH")

kruskal.test(blooms_manche$P_bac~blooms_manche$Moment)
DunnTest(blooms_manche$P_bac~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_dino~blooms_manche$Moment)
DunnTest(blooms_manche$P_dino~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_autres~blooms_manche$Moment)
DunnTest(blooms_manche$P_bac~blooms_manche$Moment,method = "BH")

kruskal.test(blooms_manche$P_BacBac~blooms_manche$Moment)
DunnTest(blooms_manche$P_BacBac~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_BacDino~blooms_manche$Moment)
DunnTest(blooms_manche$P_BacDino~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_DinoDino~blooms_manche$Moment)
DunnTest(blooms_manche$P_DinoDino~blooms_manche$Moment,method = "BH")

kruskal.test(blooms_manche$PP_BacBac~blooms_manche$Moment)
DunnTest(blooms_manche$PP_BacBac~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$PP_BacDino~blooms_manche$Moment)
DunnTest(blooms_manche$PP_BacDino~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$PP_DinoDino~blooms_manche$Moment)
DunnTest(blooms_manche$PP_DinoDino~blooms_manche$Moment,method = "BH")

# Atlantic ocean
blooms_atlantique <- filter(blooms, region == "3-Atlantic - Western Channel")
kruskal.test(blooms_atlantique$RShannon~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$RShannon~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$RPielou~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$RPielou~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$RBergerParker~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$RBergerParker~blooms_atlantique$Moment,method = "BH")

kruskal.test(blooms_atlantique$Shannon~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$Shannon~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$Pielou~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$Pielou~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$BergerParker~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$BergerParker~blooms_atlantique$Moment,method = "BH")

kruskal.test(blooms_atlantique$Dim.1~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$Dim.1~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$Dim.2~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$Dim.2~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$Dim.3~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$Dim.3~blooms_atlantique$Moment,method = "BH")

kruskal.test(blooms_atlantique$P_bac~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$P_bac~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$P_dino~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$P_dino~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$P_autres~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$P_bac~blooms_atlantique$Moment,method = "BH")

kruskal.test(blooms_atlantique$P_BacBac~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$P_BacBac~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$P_BacDino~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$P_BacDino~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$P_DinoDino~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$P_DinoDino~blooms_atlantique$Moment,method = "BH")

kruskal.test(blooms_atlantique$PP_BacBac~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$PP_BacBac~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$PP_BacDino~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$PP_BacDino~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$PP_DinoDino~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$PP_DinoDino~blooms_atlantique$Moment,method = "BH")

### Statistical tests on the diatoms blooms
# Mediterranean sea
blooms_med <- filter(blooms, region == "1-Mediterranean sea" & Bloom_Phylum == "Dino")
kruskal.test(blooms_med$RShannon~blooms_med$Moment)
DunnTest(blooms_med$RShannon~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$RPielou~blooms_med$Moment)
DunnTest(blooms_med$RPielou~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$RBergerParker~blooms_med$Moment)
DunnTest(blooms_med$RBergerParker~blooms_med$Moment,method = "BH")

kruskal.test(blooms_med$Shannon~blooms_med$Moment)
DunnTest(blooms_med$Shannon~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$Pielou~blooms_med$Moment)
DunnTest(blooms_med$Pielou~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$BergerParker~blooms_med$Moment)
DunnTest(blooms_med$BergerParker~blooms_med$Moment,method = "BH")

kruskal.test(blooms_med$Dim.1~blooms_med$Moment)
DunnTest(blooms_med$Dim.1~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$Dim.2~blooms_med$Moment)
DunnTest(blooms_med$Dim.2~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$Dim.3~blooms_med$Moment)
DunnTest(blooms_med$Dim.3~blooms_med$Moment,method = "BH")

kruskal.test(blooms_med$P_bac~blooms_med$Moment)
DunnTest(blooms_med$P_bac~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$P_dino~blooms_med$Moment)
DunnTest(blooms_med$P_dino~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$P_autres~blooms_med$Moment)
DunnTest(blooms_med$P_bac~blooms_med$Moment,method = "BH")

kruskal.test(blooms_med$P_BacBac~blooms_med$Moment)
DunnTest(blooms_med$P_BacBac~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$P_BacDino~blooms_med$Moment)
DunnTest(blooms_med$P_BacDino~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$P_DinoDino~blooms_med$Moment)
DunnTest(blooms_med$P_DinoDino~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$P_AAutres~blooms_med$Moment)
DunnTest(blooms_med$P_AAutres~blooms_med$Moment,method = "BH")

kruskal.test(blooms_med$PP_BacBac~blooms_med$Moment)
DunnTest(blooms_med$PP_BacBac~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$PP_BacDino~blooms_med$Moment)
DunnTest(blooms_med$PP_BacDino~blooms_med$Moment,method = "BH")
kruskal.test(blooms_med$PP_DinoDino~blooms_med$Moment)
DunnTest(blooms_med$PP_DinoDino~blooms_med$Moment,method = "BH")

# Eastern channel
blooms_manche <- filter(blooms, region == "2-Eastern Channel - North Sea" & Bloom_Phylum == "Dino")
kruskal.test(blooms_manche$RShannon~blooms_manche$Moment)
DunnTest(blooms_manche$RShannon~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$RPielou~blooms_manche$Moment)
DunnTest(blooms_manche$RPielou~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$RBergerParker~blooms_manche$Moment)
DunnTest(blooms_manche$RBergerParker~blooms_manche$Moment,method = "BH")

kruskal.test(blooms_manche$Shannon~blooms_manche$Moment)
DunnTest(blooms_manche$Shannon~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$Pielou~blooms_manche$Moment)
DunnTest(blooms_manche$Pielou~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$BergerParker~blooms_manche$Moment)
DunnTest(blooms_manche$BergerParker~blooms_manche$Moment,method = "BH")

kruskal.test(blooms_manche$Dim.1~blooms_manche$Moment)
DunnTest(blooms_manche$Dim.1~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$Dim.2~blooms_manche$Moment)
DunnTest(blooms_manche$Dim.2~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$Dim.3~blooms_manche$Moment)
DunnTest(blooms_manche$Dim.3~blooms_manche$Moment,method = "BH")

kruskal.test(blooms_manche$P_bac~blooms_manche$Moment)
DunnTest(blooms_manche$P_bac~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_dino~blooms_manche$Moment)
DunnTest(blooms_manche$P_dino~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_autres~blooms_manche$Moment)
DunnTest(blooms_manche$P_autres~blooms_manche$Moment,method = "BH")

kruskal.test(blooms_manche$P_BacBac~blooms_manche$Moment)
DunnTest(blooms_manche$P_BacBac~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_BacDino~blooms_manche$Moment)
DunnTest(blooms_manche$P_BacDino~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_DinoDino~blooms_manche$Moment)
DunnTest(blooms_manche$P_DinoDino~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_AAutres~blooms_manche$Moment)
DunnTest(blooms_manche$P_AAutres~blooms_manche$Moment,method = "BH")

kruskal.test(blooms_manche$PP_BacBac~blooms_manche$Moment)
DunnTest(blooms_manche$PP_BacBac~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$PP_BacDino~blooms_manche$Moment)
DunnTest(blooms_manche$PP_BacDino~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$PP_DinoDino~blooms_manche$Moment)
DunnTest(blooms_manche$PP_DinoDino~blooms_manche$Moment,method = "BH")

# Atlantic ocean
blooms_atlantique <- filter(blooms, region == "3-Atlantic - Western Channel" & Bloom_Phylum == "Dino")
kruskal.test(blooms_atlantique$RShannon~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$RShannon~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$RPielou~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$RPielou~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$RBergerParker~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$RBergerParker~blooms_atlantique$Moment,method = "BH")

kruskal.test(blooms_atlantique$Shannon~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$Shannon~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$Pielou~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$Pielou~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$BergerParker~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$BergerParker~blooms_atlantique$Moment,method = "BH")

kruskal.test(blooms_atlantique$Dim.1~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$Dim.1~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$Dim.2~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$Dim.2~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$Dim.3~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$Dim.3~blooms_atlantique$Moment,method = "BH")

kruskal.test(blooms_atlantique$P_bac~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$P_bac~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$P_dino~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$P_dino~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$P_autres~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$P_autres~blooms_atlantique$Moment,method = "BH")

kruskal.test(blooms_atlantique$P_BacBac~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$P_BacBac~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$P_BacDino~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$P_BacDino~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$P_DinoDino~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$P_DinoDino~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$P_AAutres~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$P_AAutres~blooms_atlantique$Moment,method = "BH")

kruskal.test(blooms_atlantique$PP_BacBac~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$PP_BacBac~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$PP_BacDino~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$PP_BacDino~blooms_atlantique$Moment,method = "BH")
kruskal.test(blooms_atlantique$PP_DinoDino~blooms_atlantique$Moment)
DunnTest(blooms_atlantique$PP_DinoDino~blooms_atlantique$Moment,method = "BH")

# Eastern channel and Haptophyte's blooms
blooms_manche <- filter(blooms, region == "2-Eastern Channel - North Sea" & Bloom_Phylum == "Hapt")

kruskal.test(blooms_manche$Shannon~blooms_manche$Moment)
DunnTest(blooms_manche$Shannon~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$Pielou~blooms_manche$Moment)
DunnTest(blooms_manche$Pielou~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$BergerParker~blooms_manche$Moment)
DunnTest(blooms_manche$BergerParker~blooms_manche$Moment,method = "BH")

kruskal.test(blooms_manche$Dim.1~blooms_manche$Moment)
DunnTest(blooms_manche$Dim.1~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$Dim.2~blooms_manche$Moment)
DunnTest(blooms_manche$Dim.2~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$Dim.3~blooms_manche$Moment)
DunnTest(blooms_manche$Dim.3~blooms_manche$Moment,method = "BH")

kruskal.test(blooms_manche$P_bac~blooms_manche$Moment)
DunnTest(blooms_manche$P_bac~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_dino~blooms_manche$Moment)
DunnTest(blooms_manche$P_dino~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_autres~blooms_manche$Moment)
DunnTest(blooms_manche$P_autres~blooms_manche$Moment,method = "BH")

kruskal.test(blooms_manche$P_BacBac~blooms_manche$Moment)
DunnTest(blooms_manche$P_BacBac~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_DinoDino~blooms_manche$Moment)
DunnTest(blooms_manche$P_DinoDino~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_BacDino~blooms_manche$Moment)
DunnTest(blooms_manche$P_BacDino~blooms_manche$Moment,method = "BH")
kruskal.test(blooms_manche$P_AAutres~blooms_manche$Moment)
DunnTest(blooms_manche$P_AAutres~blooms_manche$Moment,method = "BH")

# NEW VIZ Impact of blooms we consider only dino, diatoms and Eastern Channel's Haptophytes blooms #####

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

# Calculer les positions empilées des barres
blooms_summary_compo <- blooms_summary %>%
  filter(name %in% c("Bac", "Dino", "Other taxa")) %>%
  group_by(Moment, region,Bloom_Phylum) %>%
  mutate(cumsum_Val = cumsum(Val))  # Calcul de la somme cumulée


compo <- ggplot(filter(blooms_summary_compo, name %in% c("Bac", "Dino", "Other taxa") & (Bloom_Phylum=="Bacillariophyceae" | Bloom_Phylum=="Dinophyceae")), aes(Moment, Val)) +
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
  facet_grid(Bloom_Phylum~region)+
  labs(x=NULL,y="")+
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0,0.8))

# Calculer les positions empilées des barres
blooms_summary_asso <- blooms_summary %>%
  filter(name %in% c("Bac-Bac", "Bac-Dino","Dino-Dino", "Other asso.")) %>%
  group_by(Moment, region,Bloom_Phylum) %>%
  mutate(cumsum_Val = cumsum(Val))  # Calcul de la somme cumulée

asso <- ggplot(filter(blooms_summary_asso, name %in% c("Bac-Bac", "Bac-Dino","Dino-Dino", "Other asso.") & (Bloom_Phylum=="Bacillariophyceae" | Bloom_Phylum=="Dinophyceae")), aes(Moment, Val)) +
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
  facet_grid(Bloom_Phylum~region)+
  labs(x=NULL,y="")+
  theme(legend.position = "none")+
  scale_y_continuous(limits = c(0,0.8))

plot_grid(compo,asso)



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
ggsave('asso_compo_dynam_blooms_review.png', path = "output/graphs/bloom", dpi = 600, width = 360, height = 170, units = 'mm')


# Analysis : time span beetween bloom phases 
# Import data
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_PCA_coord_moments.final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)
span <- select(data,Code_point_Libelle, Date, Moment, Bloom_Phylum, region)

# Creating a df to store the results
data_results_ts<- c("","")
data_results_ts <- as.data.frame(data_results_ts)

for (i in 1:nrow(span)) {
  j = 0
  if (span[i,]$Moment == "Before" & span[i+1,]$Moment == "During" & (span[i,]$Code_point_Libelle != span[i+1,]$Code_point_Libelle )){
    j <- j+1
    data_results_ts[j,1] <- span[i,]$Code_point_Libelle
    data_results_ts[j,2] <- span[i,]$region
    data_results_ts[j,3] <- span[i,]$Bloom_Phylum
    data_results_ts[j,4] <- span[i,]$Date - span[i+1,]$Date
    data_results_ts[j,4] <- "BD"
  }
}


res_BD <- span %>%
  arrange(Code_point_Libelle, Date) %>%
  group_by(Code_point_Libelle) %>%
  mutate(
    prev_Moment = lag(Moment),
    prev_Date   = lag(Date),
    prev_region = lag(region),
    prev_phylum = lag(Bloom_Phylum)
  ) %>%
  filter(Moment == "During", prev_Moment == "Before") %>%
  transmute(
    Code_point_Libelle = Code_point_Libelle,
    region             = prev_region,      
    Bloom_Phylum       = prev_phylum,      
    Date_diff          = as.numeric(Date - prev_Date), 
    Phase              = "BD",
    Dates_compared_txt = paste0(
      "Before: ", format(prev_Date, "%Y-%m-%d"),
      " | During: ", format(Date, "%Y-%m-%d")
    )
  ) %>%
  ungroup()

# --- Transitions During -> After (Phase = DA)
res_DA <- span %>%
  arrange(Code_point_Libelle, Date) %>%
  group_by(Code_point_Libelle) %>%
  mutate(
    prev_Moment = lag(Moment),
    prev_Date   = lag(Date),
    prev_region = lag(region),
    prev_phylum = lag(Bloom_Phylum)
  ) %>%
  filter(Moment == "After", prev_Moment == "During") %>%
  transmute(
    Code_point_Libelle = Code_point_Libelle,
    region             = prev_region,      # région de la phase During
    Bloom_Phylum       = prev_phylum,      # phylum de la phase During
    Date_diff          = as.numeric(Date - prev_Date), # jours
    Phase              = "DA",
    Dates_compared_txt = paste0(
      "During: ", format(prev_Date, "%Y-%m-%d"),
      " | After: ", format(Date, "%Y-%m-%d")
    )
  ) %>%
  ungroup()

data_results_ts <- bind_rows(res_BD, res_DA) %>%
  arrange(Code_point_Libelle, Phase)

data_results_ts$Date_diff <- ifelse(data_results_ts$Phase == "BD",data_results_ts$Date_diff * -1,data_results_ts$Date_diff)

data_results_ts <- filter(data_results_ts,Bloom_Phylum %in% c("Bac","Dino","Hapt"))

data_results_ts$Bloom_Phylum <- ifelse(data_results_ts$Bloom_Phylum == "Bac","Bacillariophyceae",
                                       ifelse(data_results_ts$Bloom_Phylum == "Dino","Dinophyceae",
                                              "Haptophytes"))


ggplot(
  filter(data_results_ts, !(Bloom_Phylum == "Haptophytes" & region == "3-Atlantic - Western Channel"))
) +
  geom_histogram(aes(x = Date_diff, fill = Bloom_Phylum), binwidth = 1,colour = "black") +
  geom_text(aes(x = -5, y = 30, label = "Before")) +
  geom_text(aes(x = 4,  y = 30, label = "After")) +
  geom_vline(aes(xintercept = 0)) +
  facet_grid(region ~ Bloom_Phylum) +
  labs(x = "Days between bloom phases", y = "Number of dates", fill = "Phylum") +
  scale_x_continuous(
    breaks = seq(-30, 30, by = 5),
    limits = c(-30, 30),
    minor_breaks = NULL
  )+
  scale_fill_manual(values = c(
    "Bacillariophyceae" = "#1f77b4", "Dinophyceae" = "green",
    "Haptophytes"= "grey"))

ggsave('ts_blooms_review.png', path = "output/graphs/bloom", dpi = 600, width = 300, height = 200, units = 'mm')


