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
