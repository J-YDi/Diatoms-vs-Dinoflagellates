################################################################################
# Diatoms vs dinoflagellates: a network analysis of bloom impacts on diversity #
#                    and phytoplankton associations | R scripts                #
################################################################################

# Script used to create the map #
# 02/26/2026

# Load packages
library(ggplot2)
library(maps)
library(mapdata)
library(readr)
library(ggthemes)
library(dplyr)


# Import data 
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

data$region <- as.factor(data$region)
data$lon <- as.numeric(data$lon)
data$lat <- as.numeric(data$lat)


# Color for each region
region_col <- c("1-Mediterranean sea" = "#F8766D","2-Eastern Channel - North Sea" = "#CD9600", 
                "3-Atlantic - Western Channel" = "#00BE67",  "4-Pertuis Sea" = "gray85")

# Making the map
Worldmap <- map_data('worldHires')

data <- filter(data, region != "4-Pertuis Sea")

ggplot() + geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), fill = "gray", color = 'gray10', size = .25)+
  coord_fixed(xlim=c(-5.3,9.9), ylim=c(41,51.5), ratio=1.4)+
  labs(y = 'Latitude (degrees)', x = 'Longitude (degrees)')+
  theme_gdocs()+
  geom_point(data = data, aes(x = lon, y = lat,colour=region), size =12)+
  geom_text(data = filter(data,Code_point_Libelle == "Teychan bis") ,aes(x = lon + 1.25, y = lat,label = "Teychan bis" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  geom_text(data = filter(data,Code_point_Libelle == "Ouest Loscolo") ,aes(x = lon + 1.7, y = lat,label = "Ouest Loscolo" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  geom_text(data = filter(data,Code_point_Libelle == "Men er Roue") ,aes(x = lon - 1.6, y = lat,label = "Men er Roue" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  #geom_text(data = filter(data,Code_point_Libelle == "Le Cornard") ,aes(x = lon + 1.1, y = lat,label = "Le Cornard" ), stat = 'unique', size = 3,color="grey89") +
  geom_text(data = filter(data,Code_point_Libelle == "Bois de la Chaise large") ,aes(x = lon + 2.65, y = lat,label = "Bois de la Chaise large" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  #geom_text(data = filter(data,Code_point_Libelle == "Auger") ,aes(x = lon + .8, y = lat,label = "Auger" ), stat = 'unique', size = 3,color="grey89") +
  geom_text(data = filter(data,Code_point_Libelle == "Point 1 Boulogne") ,aes(x = lon + 2.1, y = lat,label = "Point 1 Boulogne" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  geom_text(data = filter(data,Code_point_Libelle == "Loguivy") ,aes(x = lon - 1, y = lat+0.24,label = "Loguivy" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  
  geom_text(data = filter(data,Code_point_Libelle == "les Hébihens") ,aes(x = lon + 1.65, y = lat,label = "les Hébihens" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  geom_text(data = filter(data,Code_point_Libelle == "Géfosse") ,aes(x = lon - 1.2, y = lat,label = "Géfosse" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  geom_text(data = filter(data,Code_point_Libelle == "Cabourg") ,aes(x = lon + 1.3, y = lat,label = "Cabourg" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  geom_text(data = filter(data,Code_point_Libelle == "Antifer ponton pétrolier") ,aes(x = lon + 2.7, y = lat,label = "Antifer ponton pétrolier" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  geom_text(data = filter(data,Code_point_Libelle == "Sète mer") ,aes(x = lon - 1.2, y = lat-0.05,label = "Sète mer" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  geom_text(data = filter(data,Code_point_Libelle == "Parc Leucate 2")[1,] ,aes(x = lon - 1.7, y = lat-0.13,label = "Parc Leucate 2" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  geom_text(data = filter(data,Code_point_Libelle == "Diana centre")[1,] ,aes(x = lon - 1.55, y = lat,label = "Diana centre" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  
  geom_text(data = filter(data,Code_point_Libelle == "Calvi") ,aes(x = lon + .9, y = lat,label = "Calvi" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  geom_text(data = filter(data,Code_point_Libelle == "Bouzigues (a)")[1,] ,aes(x = lon - 1.6, y = lat+0.16,label = "Bouzigues (a)" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  geom_text(data = filter(data,Code_point_Libelle == "Barcares") ,aes(x = lon - 1.3, y = lat + 0.2,label = "Barcares" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  geom_text(data = filter(data,Code_point_Libelle == "Anse de Carteau 2") ,aes(x = lon + 2, y = lat,label = "Anse de Carteau 2" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  geom_text(data = filter(data,Code_point_Libelle == "22B - Toulon gde rade") ,aes(x = lon +2.5, y = lat,label = "22B - Toulon gde rade" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  geom_text(data = filter(data,Code_point_Libelle == "At so") ,aes(x = lon + .9, y = lat,label = "At so" ), stat = 'unique', size = 7,color="black",fontface = "bold") +
  geom_text(data = filter(data,Code_point_Libelle == "At so") ,aes(x = 3.5, y = 48,label = "FRANCE" ), stat = 'unique', size = 9,color="burlywood4",fontface = "bold") +
  
  geom_text(data = filter(data,region == "1-Mediterranean sea") ,aes(x = lon, y = lat,label = "1" ), stat = 'unique', size = 6,color="gray75") +
  geom_text(data = filter(data,region == "2-Eastern Channel - North Sea") ,aes(x = lon, y = lat,label = "2" ), stat = 'unique', size = 6,color="gray75") +
  geom_text(data = filter(data,region == "3-Atlantic - Western Channel") ,aes(x = lon, y = lat,label = "3" ), stat = 'unique', size = 6,color="gray75") +
  #geom_text(data = filter(data,region == "4-Pertuis Sea") ,aes(x = lon, y = lat,label = "4" ), stat = 'unique', size = 4,color="gray80") +
  
  geom_text(aes(x = -3.5, y = 46, label = 'Atlantic ocean'), stat = 'unique', size = 5,color = "lightblue4",fontface="italic",angle=-45) +
  geom_text(aes(x = -3.2, y = 50, label = 'English channel'), stat = 'unique', size = 5,color = "lightblue4",fontface="italic",angle=28) +
  geom_text(aes(x = 5, y = 42, label = "Mediterranean sea"), stat = 'unique', size = 5,color = "lightblue4",fontface="italic",angle=8) +
  
  scale_colour_manual(values = region_col,name="Region")+
  theme(
    panel.grid.major = element_line(color = 'gray10', size = .25),
    panel.background = element_rect(fill = 'lightblue1'),
    panel.border = element_blank(),
    legend.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size=18, color = "black"),
    legend.title = element_text(size = 17, color = "black"),
    plot.background = element_rect(color = "white", linewidth = 0.001),
    legend.position = "bottom",
    plot.title = element_text(size=16, color ="black"),
    axis.text = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 7.5),
      direction = "horizontal"
    )
  )

# Save it 
ggsave('maps_station_lgtermregion_review.png', path = "output/graphs/maps", dpi = 600, width = 300, height = 300, units = 'mm')
