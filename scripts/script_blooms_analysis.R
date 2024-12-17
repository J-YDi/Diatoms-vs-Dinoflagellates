# Script Suite Stage JY Dias # 10/12/2024

library(ggforce)

# Difference on percentage of dominance (Berger-Parker index) during between bloom and no bloom ####

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

BKglobal <- ggplot(filter(data, region != "4-Pertuis Sea"))+
  geom_violin(aes(x=TBloom,y=BergerParker,group=TBloom,colour=as.character(region)),size=1,draw_quantiles = c(0.5),show.legend = FALSE,fill = "gray91")+
  geom_sina(aes(x=TBloom,y=BergerParker,group=TBloom,colour = as.character(region)),alpha = .55,show.legend = F)+
  labs(x=NULL,y="")+
  scale_y_continuous(limits=c(0,1.10),breaks = c(0,0.33,0.66,1))+
  theme_bw()

BKregion <- ggplot(filter(data, region !="4-Pertuis Sea"))+
  geom_violin(aes(x=TBloom,y=BergerParker,group=TBloom,colour=as.character(region)),size=1,draw_quantiles = c(0.5),show.legend = FALSE,fill = "gray91")+
  labs(x=NULL,y="Berger-Parker index")+
  geom_sina(aes(x=TBloom,y=BergerParker,group=TBloom,colour=as.character(region)),alpha = .55,show.legend = F)+
  scale_colour_manual(values=region_col,guide = "none")+
  scale_y_continuous(limits=c(0,1.10),breaks = c(0,0.33,0.66,1))+
  facet_wrap(~region)+
  theme_bw()

data$Bloom_Phylum[is.na(data$Bloom_Phylum)] <- "Non-bloom"

BKdetail <- ggplot(filter(data, region != "4-Pertuis Sea", Bloom_Phylum == "Bac" | Bloom_Phylum == "Dino" | Bloom_Phylum == "Hapt" | Bloom_Phylum == "Non-bloom"))+
  geom_violin(aes(x=Bloom_Phylum,y=BergerParker,group=Bloom_Phylum,colour=as.character(region)),size=1,draw_quantiles = c(0.5),show.legend = FALSE,fill = "gray91")+
  labs(x="Type of bloom",y="")+
  scale_colour_manual(values=region_col,guide = "none")+
  scale_y_continuous(limits=c(0,1.10),breaks = c(0,0.33,0.66,1))+
  geom_sina(aes(x=Bloom_Phylum,y=BergerParker,group=Bloom_Phylum,colour = as.character(region)),alpha = .55,show.legend = F)+
  facet_wrap(~region,scales = "free_x")+
  theme_bw()


plot_grid(BKglobal,BKregion,BKdetail,ncol=1,labels = "AUTO")
ggsave('BergerParker_bloom.png', path = "output/graphs/bloom", dpi = 600, width = 200, height = 200, units = 'mm')

# Statistical tests
wilcox.test(filter(data,TBloom == "Bloom")$BergerParker,filter(data,TBloom != "Bloom")$BergerParker)

wilcox.test(filter(data,TBloom == "Bloom" & region == "1-Mediterranean sea")$BergerParker,filter(data,TBloom != "Bloom" & region == "1-Mediterranean sea")$BergerParker)
wilcox.test(filter(data,TBloom == "Bloom" & region == "2-Eastern Channel - North Sea")$BergerParker,filter(data,TBloom != "Bloom" & region == "2-Eastern Channel - North Sea")$BergerParker)
wilcox.test(filter(data,TBloom == "Bloom" & region == "3-Atlantic - Western Channel")$BergerParker,filter(data,TBloom != "Bloom" & region == "3-Atlantic - Western Channel")$BergerParker)



wilcox.test(filter(data,Bloom_Phylum == "Dino" & region == "1-Mediterranean sea")$BergerParker,filter(data,TBloom != "Bloom" & region == "1-Mediterranean sea")$BergerParker)
wilcox.test(filter(data,Bloom_Phylum == "Dino" & region == "2-Eastern Channel - North Sea")$BergerParker,filter(data,TBloom != "Bloom" & region == "2-Eastern Channel - North Sea")$BergerParker)
wilcox.test(filter(data,Bloom_Phylum == "Dino" & region == "3-Atlantic - Western Channel")$BergerParker,filter(data,TBloom != "Bloom" & region == "3-Atlantic - Western Channel")$BergerParker)

wilcox.test(filter(data,Bloom_Phylum == "Bac" & region == "1-Mediterranean sea")$BergerParker,filter(data,TBloom != "Bloom" & region == "1-Mediterranean sea")$BergerParker)
wilcox.test(filter(data,Bloom_Phylum == "Bac" & region == "2-Eastern Channel - North Sea")$BergerParker,filter(data,TBloom != "Bloom" & region == "2-Eastern Channel - North Sea")$BergerParker)
wilcox.test(filter(data,Bloom_Phylum == "Bac" & region == "3-Atlantic - Western Channel")$BergerParker,filter(data,TBloom != "Bloom" & region == "3-Atlantic - Western Channel")$BergerParker)

wilcox.test(filter(data,Bloom_Phylum == "Hapt" & region == "2-Eastern Channel - North Sea")$BergerParker,filter(data,TBloom != "Bloom" & region == "2-Eastern Channel - North Sea")$BergerParker)


#### Impact of blooms we consider only dino and diatoms blooms 
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_PCA_coord_moments.final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

blooms <- filter(data, Bloom_Phylum == "Bac" | Bloom_Phylum == "Dino")

# Color for each region
region_col <- c("1-Mediterranean sea" = "#F8766D","2-Eastern Channel - North Sea" = "#CD9600", 
                "3-Atlantic - Western Channel" = "#00BE67",  "4-Pertuis Sea" = "#00A9FF"
                ,"Bac" = "#1f77b4","Dino" = "green")

blooms$Moment <-
  factor(blooms$Moment,
         levels = c("Before","During","After"))

datag <- pivot_longer(blooms, cols = c(RShannon, RPielou, RBergerParker),names_to = "Indice")

ggplot(datag) +
  
  # Violin plot avec la médiane
  geom_violin(aes(x = Moment, y = value, group = Moment, colour = as.character(region)), 
              size = 0.8, draw_quantiles = c(0.5), show.legend = FALSE, fill = "gray96") +
  
  # Ajout des lignes représentant les moyennes par Bloom_Phylum
  stat_summary(aes(x = Moment, y = value, colour = Bloom_Phylum),
               fun = mean, geom = "pointrange", alpha = 0.4, 
               size = 1, show.legend = FALSE) +
  
  # Ajout des points individuels avec geom_sina
  geom_sina(aes(x = Moment, y = value, group = Moment, colour = Bloom_Phylum), 
            alpha = 0.20, show.legend = TRUE) +
  
  # Labels et configuration des axes
  labs(x = "Moment", y = "Value") +
  scale_colour_manual(values = region_col, guide = "none") +
  
  # Facet_grid pour aligner les échelles par ligne
  facet_grid(rows = vars(Indice), cols = vars(region), scales = "free_y") +
  theme_bw()
ggsave('Diversityindex_moments_combine.png', path = "output/graphs/bloom", dpi = 600, width = 200, height = 200, units = 'mm')

datag <- pivot_longer(blooms, cols = c(Shannon, Pielou, BergerParker),names_to = "Indice")

ggplot(datag) +
  
  # Violin plot avec la médiane
  geom_violin(aes(x = Moment, y = value, group = Moment, colour = as.character(region)), 
              size = 0.8, draw_quantiles = c(0.5), show.legend = FALSE, fill = "gray96") +
  
  # Ajout des lignes représentant les moyennes par Bloom_Phylum
  stat_summary(aes(x = Moment, y = value, colour = Bloom_Phylum),
               fun = mean, geom = "pointrange", alpha = 0.4, 
               size = 1, show.legend = FALSE) +
  
  # Ajout des points individuels avec geom_sina
  geom_sina(aes(x = Moment, y = value, group = Moment, colour = Bloom_Phylum), 
            alpha = 0.20, show.legend = TRUE) +
  
  # Labels et configuration des axes
  labs(x = "Moment", y = "Value") +
  scale_colour_manual(values = region_col, guide = "none") +
  
  # Facet_grid pour aligner les échelles par ligne
  facet_grid(rows = vars(Indice), cols = vars(region), scales = "free_y") +
  theme_bw()
ggsave('OriginalDiversityindex_moments_combine.png', path = "output/graphs/bloom", dpi = 600, width = 200, height = 200, units = 'mm')

datag <- pivot_longer(blooms, cols = c(Dim.1,Dim.2,Dim.3),names_to = "Indice")

ggplot(datag) +
  
  # Violin plot avec la médiane
  geom_violin(aes(x = Moment, y = value, group = Moment, colour = as.character(region)), 
              size = 0.8, draw_quantiles = c(0.5), show.legend = FALSE, fill = "gray96") +
  
  # Ajout des lignes représentant les moyennes par Bloom_Phylum
  stat_summary(aes(x = Moment, y = value, colour = Bloom_Phylum),
               fun = mean, geom = "pointrange", alpha = 0.4, 
               size = 1, show.legend = FALSE) +
  
  # Ajout des points individuels avec geom_sina
  geom_sina(aes(x = Moment, y = value, group = Moment, colour = Bloom_Phylum), 
            alpha = 0.20, show.legend = TRUE) +
  
  # Labels et configuration des axes
  labs(x = "Moment", y = "Value") +
  scale_colour_manual(values = region_col, guide = "none") +
  
  # Facet_grid pour aligner les échelles par ligne
  facet_grid(rows = vars(Indice), cols = vars(region), scales = "free_y") +
  scale_y_continuous(limits = c(-9,9))+
  theme_bw()
ggsave('ACPdim_moments_combine.png', path = "output/graphs/bloom", dpi = 600, width = 200, height = 200, units = 'mm')

### Statistical tests on the dino and diatoms blooms
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

### Statistical tests on the diatoms blooms
blooms_med <- filter(blooms, region == "1-Mediterranean sea" & Bloom_Phylum == "Bac")
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

blooms_manche <- filter(blooms, region == "2-Eastern Channel - North Sea" & Bloom_Phylum == "Bac")
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

blooms_atlantique <- filter(blooms, region == "3-Atlantic - Western Channel" & Bloom_Phylum == "Bac")
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


### Calculate the mean for each variable 
### Diatoms and dinoflagellates blooms

blooms <- select(blooms, Date, Code_point_Libelle, region, Moment,Bloom_Phylum,RShannon,RPielou, RBergerParker, 
                 Shannon, Pielou, BergerParker,Dim.1,Dim.2,Dim.3,P_bac,P_dino,P_autres,
                 P_AAutres,P_BacBac,P_BacDino,P_DinoDino)

blooms_longer <- pivot_longer(blooms, cols = c(RShannon:P_DinoDino))

blooms_summary <- summarise(group_by(blooms_longer,region,Moment,name),value=mean(value,na.rm=T))

blooms_summary <- pivot_wider(blooms_summary,names_from = name,values_fn = mean)
write.csv2(blooms_summary,file="output/tableaux/Blooms/mean_variable_blooms.csv", row.names = FALSE,dec = ".")

Associations <- pivot_longer(blooms_summary,cols=c(P_AAutres,P_BacBac,P_BacDino,P_DinoDino))
Associations$value <- round(Associations$value,digits = 2)

Compositions <- pivot_longer(blooms_summary,cols=c(P_autres,P_bac,P_dino))
Compositions$value <- round(Compositions$value,digits = 2)

# Create all the df needed
# Associations 
Associations_med_before <- filter(Associations, region == "1-Mediterranean sea" & Moment == "Before") 
Associations_med_during <- filter(Associations, region == "1-Mediterranean sea" & Moment == "During") 
Associations_med_after <- filter(Associations, region == "1-Mediterranean sea" & Moment == "After") 

Associations_manche_before <- filter(Associations, region == "2-Eastern Channel - North Sea" & Moment == "Before") 
Associations_manche_during <- filter(Associations, region == "2-Eastern Channel - North Sea" & Moment == "During") 
Associations_manche_after <- filter(Associations, region == "2-Eastern Channel - North Sea" & Moment == "After") 

Associations_atlantic_before <- filter(Associations, region == "3-Atlantic - Western Channel" & Moment == "Before") 
Associations_atlantic_during <- filter(Associations, region == "3-Atlantic - Western Channel" & Moment == "During") 
Associations_atlantic_after <- filter(Associations, region == "3-Atlantic - Western Channel" & Moment == "After") 

# Compositions
Compositions_med_before <- filter(Compositions, region == "1-Mediterranean sea" & Moment == "Before") 
Compositions_med_during <- filter(Compositions, region == "1-Mediterranean sea" & Moment == "During") 
Compositions_med_after <- filter(Compositions, region == "1-Mediterranean sea" & Moment == "After") 

Compositions_manche_before <- filter(Compositions, region == "2-Eastern Channel - North Sea" & Moment == "Before") 
Compositions_manche_during <- filter(Compositions, region == "2-Eastern Channel - North Sea" & Moment == "During") 
Compositions_manche_after <- filter(Compositions, region == "2-Eastern Channel - North Sea" & Moment == "After") 

Compositions_atlantic_before <- filter(Compositions, region == "3-Atlantic - Western Channel" & Moment == "Before") 
Compositions_atlantic_during <- filter(Compositions, region == "3-Atlantic - Western Channel" & Moment == "During") 
Compositions_atlantic_after <- filter(Compositions, region == "3-Atlantic - Western Channel" & Moment == "After")

## Add some info for the graph
# Compute the cumulative percentages (top of each rectangle)
Associations_med_before$ymax <- cumsum(Associations_med_before$value)
# Compute the bottom of each rectangle
Associations_med_before$ymin <- c(0, head(Associations_med_before$ymax, n=-1))
# Compute label position
Associations_med_before$labelPosition <- (Associations_med_before$ymax + Associations_med_before$ymin) / 2
# Compute a good label
Associations_med_before$label <- paste0(Associations_med_before$value)

# Compute the cumulative percentages (top of each rectangle)
Associations_med_during$ymax <- cumsum(Associations_med_during$value)
# Compute the bottom of each rectangle
Associations_med_during$ymin <- c(0, head(Associations_med_during$ymax, n=-1))
# Compute label position
Associations_med_during$labelPosition <- (Associations_med_during$ymax + Associations_med_during$ymin) / 2
# Compute a good label
Associations_med_during$label <- paste0(Associations_med_during$value)

# Compute the cumulative percentages (top of each rectangle)
Associations_med_after$ymax <- cumsum(Associations_med_during$value)
# Compute the bottom of each rectangle
Associations_med_after$ymin <- c(0, head(Associations_med_after$ymax, n=-1))
# Compute label position
Associations_med_after$labelPosition <- (Associations_med_after$ymax + Associations_med_after$ymin) / 2
# Compute a good label
Associations_med_after$label <- paste0(Associations_med_after$value)


## Add some info for the graph
# Compute the cumulative percentages (top of each rectangle)
Associations_manche_before$ymax <- cumsum(Associations_manche_before$value)
# Compute the bottom of each rectangle
Associations_manche_before$ymin <- c(0, head(Associations_manche_before$ymax, n=-1))
# Compute label position
Associations_manche_before$labelPosition <- (Associations_manche_before$ymax + Associations_manche_before$ymin) / 2
# Compute a good label
Associations_manche_before$label <- paste0(Associations_manche_before$value)

# Compute the cumulative percentages (top of each rectangle)
Associations_manche_during$ymax <- cumsum(Associations_manche_during$value)
# Compute the bottom of each rectangle
Associations_manche_during$ymin <- c(0, head(Associations_manche_during$ymax, n=-1))
# Compute label position
Associations_manche_during$labelPosition <- (Associations_manche_during$ymax + Associations_manche_during$ymin) / 2
# Compute a good label
Associations_manche_during$label <- paste0(Associations_manche_during$value)

# Compute the cumulative percentages (top of each rectangle)
Associations_manche_after$ymax <- cumsum(Associations_manche_during$value)
# Compute the bottom of each rectangle
Associations_manche_after$ymin <- c(0, head(Associations_manche_after$ymax, n=-1))
# Compute label position
Associations_manche_after$labelPosition <- (Associations_manche_after$ymax + Associations_manche_after$ymin) / 2
# Compute a good label
Associations_manche_after$label <- paste0(Associations_manche_after$value)


## Add some info for the graph
# Compute the cumulative percentages (top of each rectangle)
Associations_atlantic_before$ymax <- cumsum(Associations_atlantic_before$value)
# Compute the bottom of each rectangle
Associations_atlantic_before$ymin <- c(0, head(Associations_atlantic_before$ymax, n=-1))
# Compute label position
Associations_atlantic_before$labelPosition <- (Associations_atlantic_before$ymax + Associations_atlantic_before$ymin) / 2
# Compute a good label
Associations_atlantic_before$label <- paste0(Associations_atlantic_before$value)

# Compute the cumulative percentages (top of each rectangle)
Associations_atlantic_during$ymax <- cumsum(Associations_atlantic_during$value)
# Compute the bottom of each rectangle
Associations_atlantic_during$ymin <- c(0, head(Associations_atlantic_during$ymax, n=-1))
# Compute label position
Associations_atlantic_during$labelPosition <- (Associations_atlantic_during$ymax + Associations_atlantic_during$ymin) / 2
# Compute a good label
Associations_atlantic_during$label <- paste0(Associations_atlantic_during$value)

# Compute the cumulative percentages (top of each rectangle)
Associations_atlantic_after$ymax <- cumsum(Associations_atlantic_during$value)
# Compute the bottom of each rectangle
Associations_atlantic_after$ymin <- c(0, head(Associations_atlantic_after$ymax, n=-1))
# Compute label position
Associations_atlantic_after$labelPosition <- (Associations_atlantic_after$ymax + Associations_atlantic_after$ymin) / 2
# Compute a good label
Associations_atlantic_after$label <- paste0(Associations_atlantic_after$value)


## Add some info for the graph
# Compute the cumulative percentages (top of each rectangle)
Compositions_med_before$ymax <- cumsum(Compositions_med_before$value)
# Compute the bottom of each rectangle
Compositions_med_before$ymin <- c(0, head(Compositions_med_before$ymax, n=-1))
# Compute label position
Compositions_med_before$labelPosition <- (Compositions_med_before$ymax + Compositions_med_before$ymin) / 2
# Compute a good label
Compositions_med_before$label <- paste0(Compositions_med_before$value)

# Compute the cumulative percentages (top of each rectangle)
Compositions_med_during$ymax <- cumsum(Compositions_med_during$value)
# Compute the bottom of each rectangle
Compositions_med_during$ymin <- c(0, head(Compositions_med_during$ymax, n=-1))
# Compute label position
Compositions_med_during$labelPosition <- (Compositions_med_during$ymax + Compositions_med_during$ymin) / 2
# Compute a good label
Compositions_med_during$label <- paste0(Compositions_med_during$value)

# Compute the cumulative percentages (top of each rectangle)
Compositions_med_after$ymax <- cumsum(Compositions_med_during$value)
# Compute the bottom of each rectangle
Compositions_med_after$ymin <- c(0, head(Compositions_med_after$ymax, n=-1))
# Compute label position
Compositions_med_after$labelPosition <- (Compositions_med_after$ymax + Compositions_med_after$ymin) / 2
# Compute a good label
Compositions_med_after$label <- paste0(Compositions_med_after$value)


## Add some info for the graph
# Compute the cumulative percentages (top of each rectangle)
Compositions_manche_before$ymax <- cumsum(Compositions_manche_before$value)
# Compute the bottom of each rectangle
Compositions_manche_before$ymin <- c(0, head(Compositions_manche_before$ymax, n=-1))
# Compute label position
Compositions_manche_before$labelPosition <- (Compositions_manche_before$ymax + Compositions_manche_before$ymin) / 2
# Compute a good label
Compositions_manche_before$label <- paste0(Compositions_manche_before$value)

# Compute the cumulative percentages (top of each rectangle)
Compositions_manche_during$ymax <- cumsum(Compositions_manche_during$value)
# Compute the bottom of each rectangle
Compositions_manche_during$ymin <- c(0, head(Compositions_manche_during$ymax, n=-1))
# Compute label position
Compositions_manche_during$labelPosition <- (Compositions_manche_during$ymax + Compositions_manche_during$ymin) / 2
# Compute a good label
Compositions_manche_during$label <- paste0(Compositions_manche_during$value)

# Compute the cumulative percentages (top of each rectangle)
Compositions_manche_after$ymax <- cumsum(Compositions_manche_during$value)
# Compute the bottom of each rectangle
Compositions_manche_after$ymin <- c(0, head(Compositions_manche_after$ymax, n=-1))
# Compute label position
Compositions_manche_after$labelPosition <- (Compositions_manche_after$ymax + Compositions_manche_after$ymin) / 2
# Compute a good label
Compositions_manche_after$label <- paste0(Compositions_manche_after$value)


## Add some info for the graph
# Compute the cumulative percentages (top of each rectangle)
Compositions_atlantic_before$ymax <- cumsum(Compositions_atlantic_before$value)
# Compute the bottom of each rectangle
Compositions_atlantic_before$ymin <- c(0, head(Compositions_atlantic_before$ymax, n=-1))
# Compute label position
Compositions_atlantic_before$labelPosition <- (Compositions_atlantic_before$ymax + Compositions_atlantic_before$ymin) / 2
# Compute a good label
Compositions_atlantic_before$label <- paste0(Compositions_atlantic_before$value)

# Compute the cumulative percentages (top of each rectangle)
Compositions_atlantic_during$ymax <- cumsum(Compositions_atlantic_during$value)
# Compute the bottom of each rectangle
Compositions_atlantic_during$ymin <- c(0, head(Compositions_atlantic_during$ymax, n=-1))
# Compute label position
Compositions_atlantic_during$labelPosition <- (Compositions_atlantic_during$ymax + Compositions_atlantic_during$ymin) / 2
# Compute a good label
Compositions_atlantic_during$label <- paste0(Compositions_atlantic_during$value)

# Compute the cumulative percentages (top of each rectangle)
Compositions_atlantic_after$ymax <- cumsum(Compositions_atlantic_during$value)
# Compute the bottom of each rectangle
Compositions_atlantic_after$ymin <- c(0, head(Compositions_atlantic_after$ymax, n=-1))
# Compute label position
Compositions_atlantic_after$labelPosition <- (Compositions_atlantic_after$ymax + Compositions_atlantic_after$ymin) / 2
# Compute a good label
Compositions_atlantic_after$label <- paste0(Compositions_atlantic_after$value)

# Doing the graph
atlantic_before <- ggplot() +
  geom_rect(data = Compositions_atlantic_before, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_atlantic_before, aes(ymax=ymax+0.02, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_atlantic_before,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_atlantic_before,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)


atlantic_during <- ggplot() +
  geom_rect(data = Compositions_atlantic_during, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_atlantic_during, aes(ymax=ymax+0.01, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_atlantic_during,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_atlantic_during,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)




atlantic_after <- ggplot() +
  geom_rect(data = Compositions_atlantic_after, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_atlantic_after, aes(ymax=ymax+0.01, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_atlantic_after,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_atlantic_after,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)

atlantic <- plot_grid(atlantic_before,atlantic_during,atlantic_after,ncol = 3)

# Doing the graph
med_before <- ggplot() +
  geom_rect(data = Compositions_med_before, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_med_before, aes(ymax=ymax, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_med_before,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_med_before,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)


med_during <- ggplot() +
  geom_rect(data = Compositions_med_during, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_med_during, aes(ymax=ymax, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_med_during,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_med_during,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)




med_after <- ggplot() +
  geom_rect(data = Compositions_med_after, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_med_after, aes(ymax=ymax, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_med_after,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_med_after,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)

med <- plot_grid(med_before,med_during,med_after,ncol = 3)

# Doing the graph
manche_before <- ggplot() +
  geom_rect(data = Compositions_manche_before, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_manche_before, aes(ymax=ymax+0.01, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_manche_before,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_manche_before,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)


manche_during <- ggplot() +
  geom_rect(data = Compositions_manche_during, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_manche_during, aes(ymax=ymax+0.01, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_manche_during,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_manche_during,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)




manche_after <- ggplot() +
  geom_rect(data = Compositions_manche_after, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_manche_after, aes(ymax=ymax+0.01, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_manche_after,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_manche_after,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)

manche <- plot_grid(manche_before,manche_during,manche_after,ncol = 3)

plot_grid(med,manche,atlantic,ncol = 1)
ggsave('asso_compo_moments_combine.png', path = "output/graphs/bloom", dpi = 600, width = 200, height = 200, units = 'mm')

#### Diatoms blooms #####
blooms_summary <- summarise(group_by(blooms_longer,region,Moment,name,Bloom_Phylum),value=mean(value,na.rm=T))

blooms_summary <- pivot_wider(blooms_summary,names_from = name,values_fn = mean)
write.csv2(blooms_summary,file="output/tableaux/Blooms/mean_variable_blooms_phylum.csv", row.names = FALSE,dec = ".")

Associations <- pivot_longer(blooms_summary,cols=c(P_AAutres,P_BacBac,P_BacDino,P_DinoDino))
Associations$value <- round(Associations$value,digits = 2)

Compositions <- pivot_longer(blooms_summary,cols=c(P_autres,P_bac,P_dino))
Compositions$value <- round(Compositions$value,digits = 2)

# Create all the df needed
# Associations 
Associations_med_before <- filter(Associations, region == "1-Mediterranean sea" & Moment == "Before" & Bloom_Phylum == "Bac") 
Associations_med_during <- filter(Associations, region == "1-Mediterranean sea" & Moment == "During" & Bloom_Phylum == "Bac") 
Associations_med_after <- filter(Associations, region == "1-Mediterranean sea" & Moment == "After" & Bloom_Phylum == "Bac") 

Associations_manche_before <- filter(Associations, region == "2-Eastern Channel - North Sea" & Moment == "Before" & Bloom_Phylum == "Bac") 
Associations_manche_during <- filter(Associations, region == "2-Eastern Channel - North Sea" & Moment == "During" & Bloom_Phylum == "Bac") 
Associations_manche_after <- filter(Associations, region == "2-Eastern Channel - North Sea" & Moment == "After" & Bloom_Phylum == "Bac") 

Associations_atlantic_before <- filter(Associations, region == "3-Atlantic - Western Channel" & Moment == "Before" & Bloom_Phylum == "Bac") 
Associations_atlantic_during <- filter(Associations, region == "3-Atlantic - Western Channel" & Moment == "During" & Bloom_Phylum == "Bac") 
Associations_atlantic_after <- filter(Associations, region == "3-Atlantic - Western Channel" & Moment == "After" & Bloom_Phylum == "Bac") 

# Compositions
Compositions_med_before <- filter(Compositions, region == "1-Mediterranean sea" & Moment == "Before" & Bloom_Phylum == "Bac") 
Compositions_med_during <- filter(Compositions, region == "1-Mediterranean sea" & Moment == "During" & Bloom_Phylum == "Bac") 
Compositions_med_after <- filter(Compositions, region == "1-Mediterranean sea" & Moment == "After" & Bloom_Phylum == "Bac") 

Compositions_manche_before <- filter(Compositions, region == "2-Eastern Channel - North Sea" & Moment == "Before" & Bloom_Phylum == "Bac") 
Compositions_manche_during <- filter(Compositions, region == "2-Eastern Channel - North Sea" & Moment == "During" & Bloom_Phylum == "Bac") 
Compositions_manche_after <- filter(Compositions, region == "2-Eastern Channel - North Sea" & Moment == "After" & Bloom_Phylum == "Bac") 

Compositions_atlantic_before <- filter(Compositions, region == "3-Atlantic - Western Channel" & Moment == "Before" & Bloom_Phylum == "Bac") 
Compositions_atlantic_during <- filter(Compositions, region == "3-Atlantic - Western Channel" & Moment == "During" & Bloom_Phylum == "Bac") 
Compositions_atlantic_after <- filter(Compositions, region == "3-Atlantic - Western Channel" & Moment == "After" & Bloom_Phylum == "Bac")

## Add some info for the graph
# Compute the cumulative percentages (top of each rectangle)
Associations_med_before$ymax <- cumsum(Associations_med_before$value)
# Compute the bottom of each rectangle
Associations_med_before$ymin <- c(0, head(Associations_med_before$ymax, n=-1))
# Compute label position
Associations_med_before$labelPosition <- (Associations_med_before$ymax + Associations_med_before$ymin) / 2
# Compute a good label
Associations_med_before$label <- paste0(Associations_med_before$value)

# Compute the cumulative percentages (top of each rectangle)
Associations_med_during$ymax <- cumsum(Associations_med_during$value)
# Compute the bottom of each rectangle
Associations_med_during$ymin <- c(0, head(Associations_med_during$ymax, n=-1))
# Compute label position
Associations_med_during$labelPosition <- (Associations_med_during$ymax + Associations_med_during$ymin) / 2
# Compute a good label
Associations_med_during$label <- paste0(Associations_med_during$value)

# Compute the cumulative percentages (top of each rectangle)
Associations_med_after$ymax <- cumsum(Associations_med_during$value)
# Compute the bottom of each rectangle
Associations_med_after$ymin <- c(0, head(Associations_med_after$ymax, n=-1))
# Compute label position
Associations_med_after$labelPosition <- (Associations_med_after$ymax + Associations_med_after$ymin) / 2
# Compute a good label
Associations_med_after$label <- paste0(Associations_med_after$value)


## Add some info for the graph
# Compute the cumulative percentages (top of each rectangle)
Associations_manche_before$ymax <- cumsum(Associations_manche_before$value)
# Compute the bottom of each rectangle
Associations_manche_before$ymin <- c(0, head(Associations_manche_before$ymax, n=-1))
# Compute label position
Associations_manche_before$labelPosition <- (Associations_manche_before$ymax + Associations_manche_before$ymin) / 2
# Compute a good label
Associations_manche_before$label <- paste0(Associations_manche_before$value)

# Compute the cumulative percentages (top of each rectangle)
Associations_manche_during$ymax <- cumsum(Associations_manche_during$value)
# Compute the bottom of each rectangle
Associations_manche_during$ymin <- c(0, head(Associations_manche_during$ymax, n=-1))
# Compute label position
Associations_manche_during$labelPosition <- (Associations_manche_during$ymax + Associations_manche_during$ymin) / 2
# Compute a good label
Associations_manche_during$label <- paste0(Associations_manche_during$value)

# Compute the cumulative percentages (top of each rectangle)
Associations_manche_after$ymax <- cumsum(Associations_manche_during$value)
# Compute the bottom of each rectangle
Associations_manche_after$ymin <- c(0, head(Associations_manche_after$ymax, n=-1))
# Compute label position
Associations_manche_after$labelPosition <- (Associations_manche_after$ymax + Associations_manche_after$ymin) / 2
# Compute a good label
Associations_manche_after$label <- paste0(Associations_manche_after$value)


## Add some info for the graph
# Compute the cumulative percentages (top of each rectangle)
Associations_atlantic_before$ymax <- cumsum(Associations_atlantic_before$value)
# Compute the bottom of each rectangle
Associations_atlantic_before$ymin <- c(0, head(Associations_atlantic_before$ymax, n=-1))
# Compute label position
Associations_atlantic_before$labelPosition <- (Associations_atlantic_before$ymax + Associations_atlantic_before$ymin) / 2
# Compute a good label
Associations_atlantic_before$label <- paste0(Associations_atlantic_before$value)

# Compute the cumulative percentages (top of each rectangle)
Associations_atlantic_during$ymax <- cumsum(Associations_atlantic_during$value)
# Compute the bottom of each rectangle
Associations_atlantic_during$ymin <- c(0, head(Associations_atlantic_during$ymax, n=-1))
# Compute label position
Associations_atlantic_during$labelPosition <- (Associations_atlantic_during$ymax + Associations_atlantic_during$ymin) / 2
# Compute a good label
Associations_atlantic_during$label <- paste0(Associations_atlantic_during$value)

# Compute the cumulative percentages (top of each rectangle)
Associations_atlantic_after$ymax <- cumsum(Associations_atlantic_during$value)
# Compute the bottom of each rectangle
Associations_atlantic_after$ymin <- c(0, head(Associations_atlantic_after$ymax, n=-1))
# Compute label position
Associations_atlantic_after$labelPosition <- (Associations_atlantic_after$ymax + Associations_atlantic_after$ymin) / 2
# Compute a good label
Associations_atlantic_after$label <- paste0(Associations_atlantic_after$value)


## Add some info for the graph
# Compute the cumulative percentages (top of each rectangle)
Compositions_med_before$ymax <- cumsum(Compositions_med_before$value)
# Compute the bottom of each rectangle
Compositions_med_before$ymin <- c(0, head(Compositions_med_before$ymax, n=-1))
# Compute label position
Compositions_med_before$labelPosition <- (Compositions_med_before$ymax + Compositions_med_before$ymin) / 2
# Compute a good label
Compositions_med_before$label <- paste0(Compositions_med_before$value)

# Compute the cumulative percentages (top of each rectangle)
Compositions_med_during$ymax <- cumsum(Compositions_med_during$value)
# Compute the bottom of each rectangle
Compositions_med_during$ymin <- c(0, head(Compositions_med_during$ymax, n=-1))
# Compute label position
Compositions_med_during$labelPosition <- (Compositions_med_during$ymax + Compositions_med_during$ymin) / 2
# Compute a good label
Compositions_med_during$label <- paste0(Compositions_med_during$value)

# Compute the cumulative percentages (top of each rectangle)
Compositions_med_after$ymax <- cumsum(Compositions_med_during$value)
# Compute the bottom of each rectangle
Compositions_med_after$ymin <- c(0, head(Compositions_med_after$ymax, n=-1))
# Compute label position
Compositions_med_after$labelPosition <- (Compositions_med_after$ymax + Compositions_med_after$ymin) / 2
# Compute a good label
Compositions_med_after$label <- paste0(Compositions_med_after$value)


## Add some info for the graph
# Compute the cumulative percentages (top of each rectangle)
Compositions_manche_before$ymax <- cumsum(Compositions_manche_before$value)
# Compute the bottom of each rectangle
Compositions_manche_before$ymin <- c(0, head(Compositions_manche_before$ymax, n=-1))
# Compute label position
Compositions_manche_before$labelPosition <- (Compositions_manche_before$ymax + Compositions_manche_before$ymin) / 2
# Compute a good label
Compositions_manche_before$label <- paste0(Compositions_manche_before$value)

# Compute the cumulative percentages (top of each rectangle)
Compositions_manche_during$ymax <- cumsum(Compositions_manche_during$value)
# Compute the bottom of each rectangle
Compositions_manche_during$ymin <- c(0, head(Compositions_manche_during$ymax, n=-1))
# Compute label position
Compositions_manche_during$labelPosition <- (Compositions_manche_during$ymax + Compositions_manche_during$ymin) / 2
# Compute a good label
Compositions_manche_during$label <- paste0(Compositions_manche_during$value)

# Compute the cumulative percentages (top of each rectangle)
Compositions_manche_after$ymax <- cumsum(Compositions_manche_during$value)
# Compute the bottom of each rectangle
Compositions_manche_after$ymin <- c(0, head(Compositions_manche_after$ymax, n=-1))
# Compute label position
Compositions_manche_after$labelPosition <- (Compositions_manche_after$ymax + Compositions_manche_after$ymin) / 2
# Compute a good label
Compositions_manche_after$label <- paste0(Compositions_manche_after$value)


## Add some info for the graph
# Compute the cumulative percentages (top of each rectangle)
Compositions_atlantic_before$ymax <- cumsum(Compositions_atlantic_before$value)
# Compute the bottom of each rectangle
Compositions_atlantic_before$ymin <- c(0, head(Compositions_atlantic_before$ymax, n=-1))
# Compute label position
Compositions_atlantic_before$labelPosition <- (Compositions_atlantic_before$ymax + Compositions_atlantic_before$ymin) / 2
# Compute a good label
Compositions_atlantic_before$label <- paste0(Compositions_atlantic_before$value)

# Compute the cumulative percentages (top of each rectangle)
Compositions_atlantic_during$ymax <- cumsum(Compositions_atlantic_during$value)
# Compute the bottom of each rectangle
Compositions_atlantic_during$ymin <- c(0, head(Compositions_atlantic_during$ymax, n=-1))
# Compute label position
Compositions_atlantic_during$labelPosition <- (Compositions_atlantic_during$ymax + Compositions_atlantic_during$ymin) / 2
# Compute a good label
Compositions_atlantic_during$label <- paste0(Compositions_atlantic_during$value)

# Compute the cumulative percentages (top of each rectangle)
Compositions_atlantic_after$ymax <- cumsum(Compositions_atlantic_during$value)
# Compute the bottom of each rectangle
Compositions_atlantic_after$ymin <- c(0, head(Compositions_atlantic_after$ymax, n=-1))
# Compute label position
Compositions_atlantic_after$labelPosition <- (Compositions_atlantic_after$ymax + Compositions_atlantic_after$ymin) / 2
# Compute a good label
Compositions_atlantic_after$label <- paste0(Compositions_atlantic_after$value)

# Doing the graph
atlantic_before <- ggplot() +
  geom_rect(data = Compositions_atlantic_before, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_atlantic_before, aes(ymax=ymax, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_atlantic_before,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_atlantic_before,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)


atlantic_during <- ggplot() +
  geom_rect(data = Compositions_atlantic_during, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_atlantic_during, aes(ymax=ymax, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_atlantic_during,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_atlantic_during,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)




atlantic_after <- ggplot() +
  geom_rect(data = Compositions_atlantic_after, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_atlantic_after, aes(ymax=ymax, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_atlantic_after,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_atlantic_after,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)

atlantic <- plot_grid(atlantic_before,atlantic_during,atlantic_after,ncol = 3)

# Doing the graph
med_before <- ggplot() +
  geom_rect(data = Compositions_med_before, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_med_before, aes(ymax=ymax, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_med_before,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_med_before,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)


med_during <- ggplot() +
  geom_rect(data = Compositions_med_during, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_med_during, aes(ymax=ymax, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_med_during,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_med_during,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)




med_after <- ggplot() +
  geom_rect(data = Compositions_med_after, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_med_after, aes(ymax=ymax, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_med_after,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_med_after,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)

med <- plot_grid(med_before,med_during,med_after,ncol = 3)

# Doing the graph
manche_before <- ggplot() +
  geom_rect(data = Compositions_manche_before, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_manche_before, aes(ymax=ymax, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_manche_before,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_manche_before,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)


manche_during <- ggplot() +
  geom_rect(data = Compositions_manche_during, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_manche_during, aes(ymax=ymax+0.01, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_manche_during,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_manche_during,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)




manche_after <- ggplot() +
  geom_rect(data = Compositions_manche_after, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_manche_after, aes(ymax=ymax+0.01, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_manche_after,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_manche_after,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)

manche <- plot_grid(manche_before,manche_during,manche_after,ncol = 3)

plot_grid(med,manche,atlantic,ncol = 1)
ggsave('asso_compo_moments_combine_diatom.png', path = "output/graphs/bloom", dpi = 600, width = 200, height = 200, units = 'mm')

#### Dinoflagellates blooms #####
blooms_summary <- summarise(group_by(blooms_longer,region,Moment,name,Bloom_Phylum),value=mean(value,na.rm=T))

blooms_summary <- pivot_wider(blooms_summary,names_from = name,values_fn = mean)
write.csv2(blooms_summary,file="output/tableaux/Blooms/mean_variable_blooms_phylum.csv", row.names = FALSE,dec = ".")

Associations <- pivot_longer(blooms_summary,cols=c(P_AAutres,P_BacBac,P_BacDino,P_DinoDino))
Associations$value <- round(Associations$value,digits = 2)

Compositions <- pivot_longer(blooms_summary,cols=c(P_autres,P_bac,P_dino))
Compositions$value <- round(Compositions$value,digits = 2)

# Create all the df needed
# Associations 
Associations_med_before <- filter(Associations, region == "1-Mediterranean sea" & Moment == "Before" & Bloom_Phylum == "Dino") 
Associations_med_during <- filter(Associations, region == "1-Mediterranean sea" & Moment == "During" & Bloom_Phylum == "Dino") 
Associations_med_after <- filter(Associations, region == "1-Mediterranean sea" & Moment == "After" & Bloom_Phylum == "Dino") 

Associations_manche_before <- filter(Associations, region == "2-Eastern Channel - North Sea" & Moment == "Before" & Bloom_Phylum == "Dino") 
Associations_manche_during <- filter(Associations, region == "2-Eastern Channel - North Sea" & Moment == "During" & Bloom_Phylum == "Dino") 
Associations_manche_after <- filter(Associations, region == "2-Eastern Channel - North Sea" & Moment == "After" & Bloom_Phylum == "Dino") 

Associations_atlantic_before <- filter(Associations, region == "3-Atlantic - Western Channel" & Moment == "Before" & Bloom_Phylum == "Dino") 
Associations_atlantic_during <- filter(Associations, region == "3-Atlantic - Western Channel" & Moment == "During" & Bloom_Phylum == "Dino") 
Associations_atlantic_after <- filter(Associations, region == "3-Atlantic - Western Channel" & Moment == "After" & Bloom_Phylum == "Dino") 

# Compositions
Compositions_med_before <- filter(Compositions, region == "1-Mediterranean sea" & Moment == "Before" & Bloom_Phylum == "Dino") 
Compositions_med_during <- filter(Compositions, region == "1-Mediterranean sea" & Moment == "During" & Bloom_Phylum == "Dino") 
Compositions_med_after <- filter(Compositions, region == "1-Mediterranean sea" & Moment == "After" & Bloom_Phylum == "Dino") 

Compositions_manche_before <- filter(Compositions, region == "2-Eastern Channel - North Sea" & Moment == "Before" & Bloom_Phylum == "Dino") 
Compositions_manche_during <- filter(Compositions, region == "2-Eastern Channel - North Sea" & Moment == "During" & Bloom_Phylum == "Dino") 
Compositions_manche_after <- filter(Compositions, region == "2-Eastern Channel - North Sea" & Moment == "After" & Bloom_Phylum == "Dino") 

Compositions_atlantic_before <- filter(Compositions, region == "3-Atlantic - Western Channel" & Moment == "Before" & Bloom_Phylum == "Dino") 
Compositions_atlantic_during <- filter(Compositions, region == "3-Atlantic - Western Channel" & Moment == "During" & Bloom_Phylum == "Dino") 
Compositions_atlantic_after <- filter(Compositions, region == "3-Atlantic - Western Channel" & Moment == "After" & Bloom_Phylum == "Dino")

## Add some info for the graph
# Compute the cumulative percentages (top of each rectangle)
Associations_med_before$ymax <- cumsum(Associations_med_before$value)
# Compute the bottom of each rectangle
Associations_med_before$ymin <- c(0, head(Associations_med_before$ymax, n=-1))
# Compute label position
Associations_med_before$labelPosition <- (Associations_med_before$ymax + Associations_med_before$ymin) / 2
# Compute a good label
Associations_med_before$label <- paste0(Associations_med_before$value)

# Compute the cumulative percentages (top of each rectangle)
Associations_med_during$ymax <- cumsum(Associations_med_during$value)
# Compute the bottom of each rectangle
Associations_med_during$ymin <- c(0, head(Associations_med_during$ymax, n=-1))
# Compute label position
Associations_med_during$labelPosition <- (Associations_med_during$ymax + Associations_med_during$ymin) / 2
# Compute a good label
Associations_med_during$label <- paste0(Associations_med_during$value)

# Compute the cumulative percentages (top of each rectangle)
Associations_med_after$ymax <- cumsum(Associations_med_during$value)
# Compute the bottom of each rectangle
Associations_med_after$ymin <- c(0, head(Associations_med_after$ymax, n=-1))
# Compute label position
Associations_med_after$labelPosition <- (Associations_med_after$ymax + Associations_med_after$ymin) / 2
# Compute a good label
Associations_med_after$label <- paste0(Associations_med_after$value)


## Add some info for the graph
# Compute the cumulative percentages (top of each rectangle)
Associations_manche_before$ymax <- cumsum(Associations_manche_before$value)
# Compute the bottom of each rectangle
Associations_manche_before$ymin <- c(0, head(Associations_manche_before$ymax, n=-1))
# Compute label position
Associations_manche_before$labelPosition <- (Associations_manche_before$ymax + Associations_manche_before$ymin) / 2
# Compute a good label
Associations_manche_before$label <- paste0(Associations_manche_before$value)

# Compute the cumulative percentages (top of each rectangle)
Associations_manche_during$ymax <- cumsum(Associations_manche_during$value)
# Compute the bottom of each rectangle
Associations_manche_during$ymin <- c(0, head(Associations_manche_during$ymax, n=-1))
# Compute label position
Associations_manche_during$labelPosition <- (Associations_manche_during$ymax + Associations_manche_during$ymin) / 2
# Compute a good label
Associations_manche_during$label <- paste0(Associations_manche_during$value)

# Compute the cumulative percentages (top of each rectangle)
Associations_manche_after$ymax <- cumsum(Associations_manche_during$value)
# Compute the bottom of each rectangle
Associations_manche_after$ymin <- c(0, head(Associations_manche_after$ymax, n=-1))
# Compute label position
Associations_manche_after$labelPosition <- (Associations_manche_after$ymax + Associations_manche_after$ymin) / 2
# Compute a good label
Associations_manche_after$label <- paste0(Associations_manche_after$value)


## Add some info for the graph
# Compute the cumulative percentages (top of each rectangle)
Associations_atlantic_before$ymax <- cumsum(Associations_atlantic_before$value)
# Compute the bottom of each rectangle
Associations_atlantic_before$ymin <- c(0, head(Associations_atlantic_before$ymax, n=-1))
# Compute label position
Associations_atlantic_before$labelPosition <- (Associations_atlantic_before$ymax + Associations_atlantic_before$ymin) / 2
# Compute a good label
Associations_atlantic_before$label <- paste0(Associations_atlantic_before$value)

# Compute the cumulative percentages (top of each rectangle)
Associations_atlantic_during$ymax <- cumsum(Associations_atlantic_during$value)
# Compute the bottom of each rectangle
Associations_atlantic_during$ymin <- c(0, head(Associations_atlantic_during$ymax, n=-1))
# Compute label position
Associations_atlantic_during$labelPosition <- (Associations_atlantic_during$ymax + Associations_atlantic_during$ymin) / 2
# Compute a good label
Associations_atlantic_during$label <- paste0(Associations_atlantic_during$value)

# Compute the cumulative percentages (top of each rectangle)
Associations_atlantic_after$ymax <- cumsum(Associations_atlantic_during$value)
# Compute the bottom of each rectangle
Associations_atlantic_after$ymin <- c(0, head(Associations_atlantic_after$ymax, n=-1))
# Compute label position
Associations_atlantic_after$labelPosition <- (Associations_atlantic_after$ymax + Associations_atlantic_after$ymin) / 2
# Compute a good label
Associations_atlantic_after$label <- paste0(Associations_atlantic_after$value)


## Add some info for the graph
# Compute the cumulative percentages (top of each rectangle)
Compositions_med_before$ymax <- cumsum(Compositions_med_before$value)
# Compute the bottom of each rectangle
Compositions_med_before$ymin <- c(0, head(Compositions_med_before$ymax, n=-1))
# Compute label position
Compositions_med_before$labelPosition <- (Compositions_med_before$ymax + Compositions_med_before$ymin) / 2
# Compute a good label
Compositions_med_before$label <- paste0(Compositions_med_before$value)

# Compute the cumulative percentages (top of each rectangle)
Compositions_med_during$ymax <- cumsum(Compositions_med_during$value)
# Compute the bottom of each rectangle
Compositions_med_during$ymin <- c(0, head(Compositions_med_during$ymax, n=-1))
# Compute label position
Compositions_med_during$labelPosition <- (Compositions_med_during$ymax + Compositions_med_during$ymin) / 2
# Compute a good label
Compositions_med_during$label <- paste0(Compositions_med_during$value)

# Compute the cumulative percentages (top of each rectangle)
Compositions_med_after$ymax <- cumsum(Compositions_med_during$value)
# Compute the bottom of each rectangle
Compositions_med_after$ymin <- c(0, head(Compositions_med_after$ymax, n=-1))
# Compute label position
Compositions_med_after$labelPosition <- (Compositions_med_after$ymax + Compositions_med_after$ymin) / 2
# Compute a good label
Compositions_med_after$label <- paste0(Compositions_med_after$value)


## Add some info for the graph
# Compute the cumulative percentages (top of each rectangle)
Compositions_manche_before$ymax <- cumsum(Compositions_manche_before$value)
# Compute the bottom of each rectangle
Compositions_manche_before$ymin <- c(0, head(Compositions_manche_before$ymax, n=-1))
# Compute label position
Compositions_manche_before$labelPosition <- (Compositions_manche_before$ymax + Compositions_manche_before$ymin) / 2
# Compute a good label
Compositions_manche_before$label <- paste0(Compositions_manche_before$value)

# Compute the cumulative percentages (top of each rectangle)
Compositions_manche_during$ymax <- cumsum(Compositions_manche_during$value)
# Compute the bottom of each rectangle
Compositions_manche_during$ymin <- c(0, head(Compositions_manche_during$ymax, n=-1))
# Compute label position
Compositions_manche_during$labelPosition <- (Compositions_manche_during$ymax + Compositions_manche_during$ymin) / 2
# Compute a good label
Compositions_manche_during$label <- paste0(Compositions_manche_during$value)

# Compute the cumulative percentages (top of each rectangle)
Compositions_manche_after$ymax <- cumsum(Compositions_manche_during$value)
# Compute the bottom of each rectangle
Compositions_manche_after$ymin <- c(0, head(Compositions_manche_after$ymax, n=-1))
# Compute label position
Compositions_manche_after$labelPosition <- (Compositions_manche_after$ymax + Compositions_manche_after$ymin) / 2
# Compute a good label
Compositions_manche_after$label <- paste0(Compositions_manche_after$value)


## Add some info for the graph
# Compute the cumulative percentages (top of each rectangle)
Compositions_atlantic_before$ymax <- cumsum(Compositions_atlantic_before$value)
# Compute the bottom of each rectangle
Compositions_atlantic_before$ymin <- c(0, head(Compositions_atlantic_before$ymax, n=-1))
# Compute label position
Compositions_atlantic_before$labelPosition <- (Compositions_atlantic_before$ymax + Compositions_atlantic_before$ymin) / 2
# Compute a good label
Compositions_atlantic_before$label <- paste0(Compositions_atlantic_before$value)

# Compute the cumulative percentages (top of each rectangle)
Compositions_atlantic_during$ymax <- cumsum(Compositions_atlantic_during$value)
# Compute the bottom of each rectangle
Compositions_atlantic_during$ymin <- c(0, head(Compositions_atlantic_during$ymax, n=-1))
# Compute label position
Compositions_atlantic_during$labelPosition <- (Compositions_atlantic_during$ymax + Compositions_atlantic_during$ymin) / 2
# Compute a good label
Compositions_atlantic_during$label <- paste0(Compositions_atlantic_during$value)

# Compute the cumulative percentages (top of each rectangle)
Compositions_atlantic_after$ymax <- cumsum(Compositions_atlantic_during$value)
# Compute the bottom of each rectangle
Compositions_atlantic_after$ymin <- c(0, head(Compositions_atlantic_after$ymax, n=-1))
# Compute label position
Compositions_atlantic_after$labelPosition <- (Compositions_atlantic_after$ymax + Compositions_atlantic_after$ymin) / 2
# Compute a good label
Compositions_atlantic_after$label <- paste0(Compositions_atlantic_after$value)

# Doing the graph
atlantic_before <- ggplot() +
  geom_rect(data = Compositions_atlantic_before, aes(ymax=ymax+0.01, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_atlantic_before, aes(ymax=ymax, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_atlantic_before,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_atlantic_before,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)


atlantic_during <- ggplot() +
  geom_rect(data = Compositions_atlantic_during, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_atlantic_during, aes(ymax=ymax, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_atlantic_during,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_atlantic_during,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)




atlantic_after <- ggplot() +
  geom_rect(data = Compositions_atlantic_after, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_atlantic_after, aes(ymax=ymax, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_atlantic_after,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_atlantic_after,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)

atlantic <- plot_grid(atlantic_before,atlantic_during,atlantic_after,ncol = 3)

# Doing the graph
med_before <- ggplot() +
  geom_rect(data = Compositions_med_before, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_med_before, aes(ymax=ymax, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_med_before,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_med_before,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)


med_during <- ggplot() +
  geom_rect(data = Compositions_med_during, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_med_during, aes(ymax=ymax, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_med_during,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_med_during,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)




med_after <- ggplot() +
  geom_rect(data = Compositions_med_after, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_med_after, aes(ymax=ymax, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_med_after,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_med_after,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)

med <- plot_grid(med_before,med_during,med_after,ncol = 3)

# Doing the graph
manche_before <- ggplot() +
  geom_rect(data = Compositions_manche_before, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_manche_before, aes(ymax=ymax+0.01, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_manche_before,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_manche_before,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)


manche_during <- ggplot() +
  geom_rect(data = Compositions_manche_during, aes(ymax=ymax+0.01, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_manche_during, aes(ymax=ymax, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_manche_during,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_manche_during,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)




manche_after <- ggplot() +
  geom_rect(data = Compositions_manche_after, aes(ymax=ymax+0.01, ymin=ymin, xmax=4, xmin=3, fill=name)) +
  geom_rect(data = Associations_manche_after, aes(ymax=ymax, ymin=ymin, xmax=2.9, xmin=2, fill=name)) +
  geom_text( data = Compositions_manche_after,x=4.7, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  geom_text( data = Associations_manche_after,x=1.2, aes(y=labelPosition, label=label, color=name), size=3) + # x here controls label position (inner / outer)
  scale_fill_manual(values = c("P_bac" = "#56B4E9",
                               "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  scale_color_manual(values = c("P_bac" = "#56B4E9",
                                "P_autres"      = "grey" ,"P_dino" = "#009E73","P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                                ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_no_axes() +
  theme(legend.position = "none")+
  facet_wrap(region~Moment)

manche <- plot_grid(manche_before,manche_during,manche_after,ncol = 3)

plot_grid(med,manche,atlantic,ncol = 1)
ggsave('asso_compo_moments_combine_dino.png', path = "output/graphs/bloom", dpi = 600, width = 200, height = 200, units = 'mm')

#### Test stats entre les moments tout blooms pour les proportions ####
Associations_med <- filter(blooms,region == "1-Mediterranean sea")
Associations_med <- dplyr::select(Associations_med,Moment,P_AAutres,P_BacBac,P_BacDino,P_DinoDino)


# MANOVA correcte avec cbind() pour les variables dépendantes
res_manova <- manova(cbind(P_BacDino, P_BacBac, P_DinoDino) ~ Moment, data = Associations_med)

# Résumé des résultats
summary(res_manova)

summary.aov(res_manova)

Compositions_med <- filter(blooms,region == "1-Mediterranean sea")
Compositions_med <- dplyr::select(Compositions_med,Moment,P_autres,P_bac,P_dino)


# MANOVA correcte avec cbind() pour les variables dépendantes
res_manova <- manova(cbind(P_bac, P_dino) ~ Moment, data = Compositions_med)

# Résumé des résultats
summary(res_manova)

summary.aov(res_manova)

# Manche

Associations_manche <- filter(blooms,region == "2-Eastern Channel - North Sea")
Associations_manche <- dplyr::select(Associations_manche,Moment,P_AAutres,P_BacBac,P_BacDino,P_DinoDino)


# MANOVA correcte avec cbind() pour les variables dépendantes
res_manova <- manova(cbind(P_BacDino, P_BacBac, P_DinoDino) ~ Moment, data = Associations_manche)

# Résumé des résultats
summary(res_manova)

summary.aov(res_manova)

summary(aov(P_BacDino ~ Moment, data = Associations_manche))
summary(aov(P_BacBac ~ Moment, data = Associations_manche))
summary(aov(P_DinoDino ~ Moment, data = Associations_manche))
summary(aov(P_AAutres ~ Moment, data = Associations_manche))

# Test de Tukey post-hoc pour chaque variable
TukeyHSD(aov(P_BacDino ~ Moment, data = Associations_manche))
TukeyHSD(aov(P_BacBac ~ Moment, data = Associations_manche))
TukeyHSD(aov(P_DinoDino ~ Moment, data = Associations_manche))
TukeyHSD(aov(P_AAutres ~ Moment, data = Associations_manche))

Compositions_manche <- filter(blooms,region == "2-Eastern Channel - North Sea")
Compositions_manche <- dplyr::select(Compositions_manche,Moment,P_autres,P_bac,P_dino)


# MANOVA correcte avec cbind() pour les variables dépendantes
res_manova <- manova(cbind(P_bac, P_dino) ~ Moment, data = Compositions_manche)

# Résumé des résultats
summary(res_manova)

summary.aov(res_manova)

summary(aov(P_bac ~ Moment, data = Compositions_manche))
summary(aov(P_dino ~ Moment, data = Compositions_manche))
summary(aov(P_autres ~ Moment, data = Compositions_manche))

# Test de Tukey post-hoc pour chaque variable
TukeyHSD(aov(P_bac ~ Moment, data = Compositions_manche))
TukeyHSD(aov(P_dino ~ Moment, data = Compositions_manche))
TukeyHSD(aov(P_autres ~ Moment, data = Compositions_manche))


# Atlantique

Associations_atlantic <- filter(blooms,region == "3-Atlantic - Western Channel")
Associations_atlantic <- dplyr::select(Associations_atlantic,Moment,P_AAutres,P_BacBac,P_BacDino,P_DinoDino)


# MANOVA correcte avec cbind() pour les variables dépendantes
res_manova <- manova(cbind(P_BacDino, P_BacBac, P_DinoDino) ~ Moment, data = Associations_atlantic)

# Résumé des résultats
summary(res_manova)

summary.aov(res_manova)

summary(aov(P_BacDino ~ Moment, data = Associations_manche))
summary(aov(P_BacBac ~ Moment, data = Associations_manche))
summary(aov(P_DinoDino ~ Moment, data = Associations_manche))
summary(aov(P_AAutres ~ Moment, data = Associations_manche))

# Test de Tukey post-hoc pour chaque variable
TukeyHSD(aov(P_BacDino ~ Moment, data = Associations_manche))
TukeyHSD(aov(P_BacBac ~ Moment, data = Associations_manche))
TukeyHSD(aov(P_DinoDino ~ Moment, data = Associations_manche))
TukeyHSD(aov(P_AAutres ~ Moment, data = Associations_manche))

Compositions_atlantic <- filter(blooms,region == "3-Atlantic - Western Channel")
Compositions_atlantic <- dplyr::select(Compositions_atlantic,Moment,P_autres,P_bac,P_dino)


# MANOVA correcte avec cbind() pour les variables dépendantes
res_manova <- manova(cbind(P_bac, P_dino) ~ Moment, data = Compositions_atlantic)

# Résumé des résultats
summary(res_manova)

summary.aov(res_manova)

summary(aov(P_bac ~ Moment, data = Compositions_atlantic))
summary(aov(P_dino ~ Moment, data = Compositions_atlantic))
summary(aov(P_autres ~ Moment, data = Compositions_atlantic))

# Test de Tukey post-hoc pour chaque variable
TukeyHSD(aov(P_bac ~ Moment, data = Compositions_atlantic))
TukeyHSD(aov(P_dino ~ Moment, data = Compositions_atlantic))
TukeyHSD(aov(P_autres ~ Moment, data = Compositions_atlantic))

#### Test stats entre les moments blooms de diatomées pour les proportions ####
Associations_med <- filter(blooms,region == "1-Mediterranean sea" & Bloom_Phylum == "Bac")
Associations_med <- dplyr::select(Associations_med,Moment,P_AAutres,P_BacBac,P_BacDino,P_DinoDino)


# MANOVA correcte avec cbind() pour les variables dépendantes
res_manova <- manova(cbind(P_BacDino, P_BacBac, P_DinoDino) ~ Moment, data = Associations_med)

# Résumé des résultats
summary(res_manova)

summary.aov(res_manova)

Compositions_med <- filter(blooms,region == "1-Mediterranean sea" & Bloom_Phylum == "Bac")
Compositions_med <- dplyr::select(Compositions_med,Moment,P_autres,P_bac,P_dino)


# MANOVA correcte avec cbind() pour les variables dépendantes
res_manova <- manova(cbind(P_bac, P_dino) ~ Moment, data = Compositions_med)

# Résumé des résultats
summary(res_manova)

summary.aov(res_manova)

# Manche

Associations_manche <- filter(blooms,region == "2-Eastern Channel - North Sea" & Bloom_Phylum == "Bac")
Associations_manche <- dplyr::select(Associations_manche,Moment,P_AAutres,P_BacBac,P_BacDino,P_DinoDino)


# MANOVA correcte avec cbind() pour les variables dépendantes
res_manova <- manova(cbind(P_BacDino, P_BacBac, P_DinoDino) ~ Moment, data = Associations_manche)

# Résumé des résultats
summary(res_manova)

summary.aov(res_manova)

summary(aov(P_BacDino ~ Moment, data = Associations_manche))
summary(aov(P_BacBac ~ Moment, data = Associations_manche))
summary(aov(P_DinoDino ~ Moment, data = Associations_manche))
summary(aov(P_AAutres ~ Moment, data = Associations_manche))

# Test de Tukey post-hoc pour chaque variable
TukeyHSD(aov(P_BacDino ~ Moment, data = Associations_manche))
TukeyHSD(aov(P_BacBac ~ Moment, data = Associations_manche))
TukeyHSD(aov(P_DinoDino ~ Moment, data = Associations_manche))
TukeyHSD(aov(P_AAutres ~ Moment, data = Associations_manche))

Compositions_manche <- filter(blooms,region == "2-Eastern Channel - North Sea" & Bloom_Phylum == "Bac")
Compositions_manche <- dplyr::select(Compositions_manche,Moment,P_autres,P_bac,P_dino)


# MANOVA correcte avec cbind() pour les variables dépendantes
res_manova <- manova(cbind(P_bac, P_dino) ~ Moment, data = Compositions_manche)

# Résumé des résultats
summary(res_manova)

summary.aov(res_manova)

summary(aov(P_bac ~ Moment, data = Compositions_manche))
summary(aov(P_dino ~ Moment, data = Compositions_manche))
summary(aov(P_autres ~ Moment, data = Compositions_manche))

# Test de Tukey post-hoc pour chaque variable
TukeyHSD(aov(P_bac ~ Moment, data = Compositions_manche))
TukeyHSD(aov(P_dino ~ Moment, data = Compositions_manche))
TukeyHSD(aov(P_autres ~ Moment, data = Compositions_manche))


# Atlantique

Associations_atlantic <- filter(blooms,region == "3-Atlantic - Western Channel" & Bloom_Phylum == "Bac")
Associations_atlantic <- dplyr::select(Associations_atlantic,Moment,P_AAutres,P_BacBac,P_BacDino,P_DinoDino)


# MANOVA correcte avec cbind() pour les variables dépendantes
res_manova <- manova(cbind(P_BacDino, P_BacBac, P_DinoDino) ~ Moment, data = Associations_atlantic)

# Résumé des résultats
summary(res_manova)

summary.aov(res_manova)

summary(aov(P_BacDino ~ Moment, data = Associations_atlantic))
summary(aov(P_BacBac ~ Moment, data = Associations_atlantic))
summary(aov(P_DinoDino ~ Moment, data = Associations_atlantic))
summary(aov(P_AAutres ~ Moment, data = Associations_atlantic))

# Test de Tukey post-hoc pour chaque variable
TukeyHSD(aov(P_BacDino ~ Moment, data = Associations_atlantic))
TukeyHSD(aov(P_BacBac ~ Moment, data = Associations_atlantic))
TukeyHSD(aov(P_DinoDino ~ Moment, data = Associations_atlantic))
TukeyHSD(aov(P_AAutres ~ Moment, data = Associations_atlantic))

Compositions_atlantic <- filter(blooms,region == "3-Atlantic - Western Channel" & Bloom_Phylum == "Bac")
Compositions_atlantic <- dplyr::select(Compositions_atlantic,Moment,P_autres,P_bac,P_dino)


# MANOVA correcte avec cbind() pour les variables dépendantes
res_manova <- manova(cbind(P_bac, P_dino) ~ Moment, data = Compositions_atlantic)

# Résumé des résultats
summary(res_manova)

summary.aov(res_manova)

summary(aov(P_bac ~ Moment, data = Compositions_atlantic))
summary(aov(P_dino ~ Moment, data = Compositions_atlantic))
summary(aov(P_autres ~ Moment, data = Compositions_atlantic))

# Test de Tukey post-hoc pour chaque variable
TukeyHSD(aov(P_bac ~ Moment, data = Compositions_atlantic))
TukeyHSD(aov(P_dino ~ Moment, data = Compositions_atlantic))
TukeyHSD(aov(P_autres ~ Moment, data = Compositions_atlantic))

#### Test stats entre les moments blooms de dinoflagellés pour les proportions ####
Associations_med <- filter(blooms,region == "1-Mediterranean sea" & Bloom_Phylum == "Dino")
Associations_med <- dplyr::select(Associations_med,Moment,P_AAutres,P_BacBac,P_BacDino,P_DinoDino)


# MANOVA correcte avec cbind() pour les variables dépendantes
res_manova <- manova(cbind(P_BacDino, P_BacBac, P_DinoDino) ~ Moment, data = Associations_med)

# Résumé des résultats
summary(res_manova)

summary.aov(res_manova)

Compositions_med <- filter(blooms,region == "1-Mediterranean sea" & Bloom_Phylum == "Dino")
Compositions_med <- dplyr::select(Compositions_med,Moment,P_autres,P_bac,P_dino)


# MANOVA correcte avec cbind() pour les variables dépendantes
res_manova <- manova(cbind(P_bac, P_dino) ~ Moment, data = Compositions_med)

# Résumé des résultats
summary(res_manova)

summary.aov(res_manova)

# Manche

Associations_manche <- filter(blooms,region == "2-Eastern Channel - North Sea" & Bloom_Phylum == "Dino")
Associations_manche <- dplyr::select(Associations_manche,Moment,P_AAutres,P_BacBac,P_BacDino,P_DinoDino)


# MANOVA correcte avec cbind() pour les variables dépendantes
res_manova <- manova(cbind(P_BacDino, P_BacBac, P_DinoDino) ~ Moment, data = Associations_manche)

# Résumé des résultats
summary(res_manova)

summary.aov(res_manova)

summary(aov(P_BacDino ~ Moment, data = Associations_manche))
summary(aov(P_BacBac ~ Moment, data = Associations_manche))
summary(aov(P_DinoDino ~ Moment, data = Associations_manche))
summary(aov(P_AAutres ~ Moment, data = Associations_manche))

# Test de Tukey post-hoc pour chaque variable
TukeyHSD(aov(P_BacDino ~ Moment, data = Associations_manche))
TukeyHSD(aov(P_BacBac ~ Moment, data = Associations_manche))
TukeyHSD(aov(P_DinoDino ~ Moment, data = Associations_manche))
TukeyHSD(aov(P_AAutres ~ Moment, data = Associations_manche))

Compositions_manche <- filter(blooms,region == "2-Eastern Channel - North Sea" & Bloom_Phylum == "Dino" )
Compositions_manche <- dplyr::select(Compositions_manche,Moment,P_autres,P_bac,P_dino)


# MANOVA correcte avec cbind() pour les variables dépendantes
res_manova <- manova(cbind(P_bac, P_dino) ~ Moment, data = Compositions_manche)

# Résumé des résultats
summary(res_manova)

summary.aov(res_manova)

summary(aov(P_bac ~ Moment, data = Compositions_manche))
summary(aov(P_dino ~ Moment, data = Compositions_manche))
summary(aov(P_autres ~ Moment, data = Compositions_manche))

# Test de Tukey post-hoc pour chaque variable
TukeyHSD(aov(P_bac ~ Moment, data = Compositions_manche))
TukeyHSD(aov(P_dino ~ Moment, data = Compositions_manche))
TukeyHSD(aov(P_autres ~ Moment, data = Compositions_manche))


# Atlantique

Associations_atlantic <- filter(blooms,region == "3-Atlantic - Western Channel" & Bloom_Phylum == "Dino")
Associations_atlantic <- dplyr::select(Associations_atlantic,Moment,P_AAutres,P_BacBac,P_BacDino,P_DinoDino)


# MANOVA correcte avec cbind() pour les variables dépendantes
res_manova <- manova(cbind(P_BacDino, P_BacBac, P_DinoDino) ~ Moment, data = Associations_atlantic)

# Résumé des résultats
summary(res_manova)

summary.aov(res_manova)

summary(aov(P_BacDino ~ Moment, data = Associations_atlantic))
summary(aov(P_BacBac ~ Moment, data = Associations_atlantic))
summary(aov(P_DinoDino ~ Moment, data = Associations_atlantic))
summary(aov(P_AAutres ~ Moment, data = Associations_atlantic))

# Test de Tukey post-hoc pour chaque variable
TukeyHSD(aov(P_BacDino ~ Moment, data = Associations_atlantic))
TukeyHSD(aov(P_BacBac ~ Moment, data = Associations_atlantic))
TukeyHSD(aov(P_DinoDino ~ Moment, data = Associations_atlantic))
TukeyHSD(aov(P_AAutres ~ Moment, data = Associations_atlantic))

Compositions_atlantic <- filter(blooms,region == "3-Atlantic - Western Channel" & Bloom_Phylum == "Dino")
Compositions_atlantic <- dplyr::select(Compositions_atlantic,Moment,P_autres,P_bac,P_dino)


# MANOVA correcte avec cbind() pour les variables dépendantes
res_manova <- manova(cbind(P_bac, P_dino) ~ Moment, data = Compositions_atlantic)

# Résumé des résultats
summary(res_manova)

summary.aov(res_manova)

summary(aov(P_bac ~ Moment, data = Compositions_atlantic))
summary(aov(P_dino ~ Moment, data = Compositions_atlantic))
summary(aov(P_autres ~ Moment, data = Compositions_atlantic))

# Test de Tukey post-hoc pour chaque variable
TukeyHSD(aov(P_bac ~ Moment, data = Compositions_atlantic))
TukeyHSD(aov(P_dino ~ Moment, data = Compositions_atlantic))
TukeyHSD(aov(P_autres ~ Moment, data = Compositions_atlantic))

