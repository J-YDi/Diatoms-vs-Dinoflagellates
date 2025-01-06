# Script Suite Stage JY Dias # 09/12/2024

# Load packages
library(ggplot2)
library(ggthemes)
library(readr)
library(dplyr)
library(tidyr)
library(FactoMineR)
library(factoextra)
library(cowplot)

# Load data
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)
# Adding season info
data <- data |>
  mutate(season = case_when(Month %in% c(12, 01, 02) ~ "Winter",
                            Month %in% c(03, 04, 05) ~ "Spring",
                            Month %in% c(06, 07, 08) ~ "Summer",
                            Month %in% c(09, 10, 11) ~ "Fall", TRUE ~ NA_character_))


data$region <- as.factor(data$region)

# Color for each region
region_col <- c("1-Mediterranean sea" = "#F8766D","2-Eastern Channel - North Sea" = "#CD9600", 
                "3-Atlantic - Western Channel" = "#00BE67",  "4-Pertuis Sea" = "#00A9FF")

data <- filter(data, region != "4-Pertuis Sea")

### Graph to show phytoplankton composition by season and cluster #####
# Class level
# Prepare data for graph
data$Others <- rowSums(select(data,Raphidophyceae:`Autres protistes`,Dinoflagellata,Xanthophyceae),na.rm = T)
datam <- summarise(group_by(data,region,season,Month,Year), Bacillariophyceae=mean(Bacillariophyceae,na.rm=T),
                   Dinophyceae=mean(Dinophyceae,na.rm=T),Cryptophyceae=mean(Cryptophyceae,na.rm=T),
                   Haptophyta=mean(Haptophyta,na.rm=T),Ciliophora=mean(Ciliophora,na.rm=T),Others=mean(Others,na.rm=T))

# Compute total abundances to have relative composition
datam$Abdtot <- rowSums(datam[,c(5:10)],na.rm=T)
datam[,c(5:10)] <- datam[,c(5:10)]/datam$Abdtot

date_string <- paste(datam$Year, datam$Month, "01", sep = "-")

# Convert to date format
datam$Date <- as.Date(date_string,format = "%Y-%m-%d")

# Make it as the graph needs
datag <- pivot_longer(data = datam,cols = Bacillariophyceae:Others,names_to = "Phylum")
datag$Abdtot <- NULL

datag$MonthYear <- format(datag$Date, "%Y-%m")
datag$Phylum <- as.factor(datag$Phylum)

datag$Phylum <-
  factor(datag$Phylum,
         levels = c("Bacillariophyceae","Dinophyceae","Ciliophora","Cryptophyceae","Haptophyta","Others"))

datag$Date <- as.Date(datag$Date,format="%Y-%m-%d")

ggplot(datag) +
  geom_col(aes(x = MonthYear, y = value * 100, fill = Phylum), position = "stack", na.rm = FALSE, width = 1) +
  facet_wrap(~region, scales = "free", ncol = 1) +
  scale_x_discrete(labels = function(x) format(as.Date(paste(x, "01", sep = "-")), "%Y-%m")) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 5))+
  labs(x="Date",y="Relative abundance (%)",fill="Taxon")+
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  scale_fill_manual(values = c("#56B4E9", "#009E73" ,"#F0E442", "#CC79A7","#996136","grey60"))+
   geom_text(data = subset(datag, format(Date, "%m") == "01"), 
            aes(x = MonthYear, y = 90, label = format(Date, "%Y")),
            color = "black", size = 2.5,angle=90, vjust = 0)+
  geom_vline(data = subset(datag, format(Date, "%m") == "01"), 
                      aes(xintercept = MonthYear),
                      color = "grey3", size = 0.5,linetype = "dashed")+
  guides(fill = guide_legend(ncol = 6, byrow = TRUE)) +
  theme(legend.position = "bottom",legend.text = element_text(size = 7),legend.title = element_text(size=7))
ggsave('Relativeabundance_class.png', path = "output/graphs/description_region",dpi = 500, width = 360, height =200, units = 'mm')

# Make it with log abundance
# Prepare data for graph
data$Others <- rowSums(select(data,Raphidophyceae:`Autres protistes`,Dinoflagellata,Xanthophyceae),na.rm = T)
datam <- summarise(group_by(data,region,season,Month,Year), Bacillariophyceae=mean(Bacillariophyceae,na.rm=T),
                   Dinophyceae=mean(Dinophyceae,na.rm=T),Cryptophyceae=mean(Cryptophyceae,na.rm=T),
                   Haptophyta=mean(Haptophyta,na.rm=T),Ciliophora=mean(Ciliophora,na.rm=T),Others=mean(Others,na.rm=T))

date_string <- paste(datam$Year, datam$Month, "01", sep = "-")

# Convert to date format
datam$Date <- as.Date(date_string,format = "%Y-%m-%d")

# Make it as the graph needs
datag <- pivot_longer(data = datam,cols = Bacillariophyceae:Others,names_to = "Phylum")
datag$MonthYear <- format(datag$Date, "%Y-%m")
datag$Phylum <- as.factor(datag$Phylum)

datag$Phylum <-
  factor(datag$Phylum,
         levels = c("Bacillariophyceae","Dinophyceae","Ciliophora","Cryptophyceae","Haptophyta","Others"))

datag$Date <- as.Date(datag$Date,format="%Y-%m-%d")

ggplot(datag) +
  geom_col(aes(x = MonthYear, y = log(value+1), fill = Phylum), position = "stack", na.rm = FALSE, width = 1) +
  facet_wrap(~region, scales = "free", ncol = 1) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 5))+
  labs(x="Date",y="Log Abundance",fill="Taxon")+
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  scale_fill_manual(values = c("#56B4E9", "#009E73" ,"#F0E442", "#CC79A7","#996136","grey60"))+
  geom_text(data = subset(datag, format(Date, "%m") == "01"), 
            aes(x = MonthYear, y = 63, label = format(Date, "%Y")),
            color = "black", size = 2.5,angle=90, vjust = 0)+
  geom_vline(data = subset(datag, format(Date, "%m") == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", size = 0.5,linetype = "dashed")+
  guides(fill = guide_legend(ncol = 6, byrow = TRUE)) +
  scale_y_continuous(limits = c(0,65))+
  theme(legend.position = "bottom",legend.text = element_text(size = 7),legend.title = element_text(size=7))
ggsave('Abundance_class.png', path = "output/graphs/description_region",dpi = 500, width = 360, height =200, units = 'mm')



# Make it by Genus
# Compute the mean of all genus by season and cluster
data_graph <- data[,c(3,354,24:328)] %>%
  group_by(season, region) %>%
  summarise_all(~mean(., na.rm = TRUE))

# Make it relative 
data_graph$Abdtot <- rowSums(data_graph[,c(3:307)],na.rm=T)
data_graph[,c(3:307)] <- data_graph[,c(3:307)]/data_graph$Abdtot

datag <- pivot_longer(data = data_graph,cols = Actinoptychus:Coscinodiscophycidae,names_to = "Phylum")
datag$Abdtot <- NULL

# Group by season and region, and select the 5 most abundant phyla for each region
# Determine the most abundant phyla for each region and season

data_graph <- datag %>%
  group_by(season, region, Phylum) %>%
  summarise(abondance = sum(value,na.rm=T)) %>%
  ungroup() %>%
  group_by(season, region) %>%
  mutate(rank = rank(desc(abondance)))

write.csv(data_graph,file="output/tableaux/phyla_mostabundant_region_season.csv")

# We select this taxa
selection_phylum <- c("Skeletonema","Pseudo-nitzschia","Chaetoceros","Chaetocerotaceae","Nitzschia","Cryptophyceae",
                      "Cryptomonadales","Cylindrotheca","Leptocylindrus","Akashiwo","Phaeocystis","Akashiwo",
                      "Phaeocystis","Asterionellopsis","Chrysochromulina","Azadinium")

datataxon <- data |>
  select(Month,Year,region,selection_phylum)
# Agglomerate Chaetoceros with Chaetocerotaceae and Cryptophyceae with Cryptomonadales
datataxon$Chaetocerotaceae <- rowSums(datataxon[,c("Chaetoceros","Chaetocerotaceae")],na.rm=T)
datataxon$Cryptophyceae <- rowSums(datataxon[,c("Cryptophyceae","Cryptomonadales")],na.rm=T)
datataxon <- select(datataxon,-c(Chaetoceros,Cryptomonadales))

# Select the others
dataothers <- data |>
  select(-selection_phylum) |>
  select(region,Actinoptychus:Coscinodiscophycidae)

dataothers$Others <- rowSums(select(dataothers,Actinoptychus:Coscinodiscophycidae),na.rm = T)
dataothers <- select(dataothers,Others)

datagraph <- cbind(datataxon,dataothers)


data_graph <- datagraph %>%
  group_by(Month,Year, region) %>%
  summarise_all(~mean(., na.rm = TRUE))

# Convert to date format
date_string <- paste(data_graph$Year, data_graph$Month, "01", sep = "-")
data_graph$Date <- as.Date(date_string,format = "%Y-%m-%d")
data_graph$MonthYear <- format(data_graph$Date, "%Y-%m")

datag <- pivot_longer(data = data_graph,cols = Skeletonema:Others,names_to = "Phylum")

datag$Phylum <- factor(datag$Phylum, levels = c("Others", "Skeletonema","Pseudo-nitzschia","Chaetoceros","Nitzschia","Cryptophyceae",
                                                "Cylindrotheca","Leptocylindrus","Akashiwo",
                                                "Phaeocystis","Asterionellopsis","Chrysochromulina","Azadinium"))

ggplot(datag) +
  geom_col(aes(x = MonthYear, y = log(value+1), fill = Phylum), position = "stack", na.rm = FALSE, width = 1) +
  facet_wrap(~region, scales = "free", ncol = 1) +
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  geom_vline(data = subset(datag, format(Date, "%m") == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", size = 0.5,linetype = "dashed")+
  geom_text(data = subset(datag, format(Date, "%m") == "01"), 
            aes(x = MonthYear, y = 93, label = format(Date, "%Y")),
            color = "black", size = 2.5,angle=90, vjust = 0)+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 5))+
  labs(x="Date",y="Log Abundance",fill="Taxon")+
  scale_y_continuous(limits = c(0,100))+
  scale_fill_manual(values = c("grey","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
                               "#FDBF6F", "#FF7F00", "#6A3D9A", "#FFFF99", "#FB9A99"
                               , "#FFBCFF" , "#BC00BC", "#510051","#B15928"))
ggsave('Abundance_genus.png', path = "output/graphs/description_region",dpi = 500, width = 360, height =200, units = 'mm')

# Make it relative
selection_phylum <- c("Skeletonema","Pseudo-nitzschia","Chaetoceros","Chaetocerotaceae","Nitzschia","Cryptophyceae",
                      "Cryptomonadales","Cylindrotheca","Leptocylindrus","Akashiwo","Phaeocystis","Akashiwo",
                      "Phaeocystis","Asterionellopsis","Chrysochromulina","Azadinium")
datataxon <- data |>
  select(Month,Year,region,selection_phylum)
# Agglomerate Chaetoceros with Chaetocerotaceae and Cryptophyceae with Cryptomonadales
datataxon$Chaetocerotaceae <- rowSums(datataxon[,c("Chaetoceros","Chaetocerotaceae")],na.rm=T)
datataxon$Cryptophyceae <- rowSums(datataxon[,c("Cryptophyceae","Cryptomonadales")],na.rm=T)
datataxon <- select(datataxon,-c(Chaetoceros,Cryptomonadales))

dataothers <- data |>
  select(-selection_phylum) |>
  select(region,Actinoptychus:Coscinodiscophycidae)

dataothers$Others <- rowSums(select(dataothers,Actinoptychus:Coscinodiscophycidae),na.rm = T)
dataothers <- select(dataothers,Others)

datagraph <- cbind(datataxon,dataothers)


data_graph <- datagraph %>%
  group_by(Month,Year, region) %>%
  summarise_all(~mean(., na.rm = TRUE))
# When we switch to relative
data_graph$Abdtot <- rowSums(data_graph[,c(4:16)],na.rm=T)
data_graph[,c(4:16)] <- data_graph[,c(4:16)]/data_graph$Abdtot

data_graph$Abdtot <- NULL


# Convert to date format
date_string <- paste(data_graph$Year, data_graph$Month, "01", sep = "-")
data_graph$Date <- as.Date(date_string,format = "%Y-%m-%d")
data_graph$MonthYear <- format(data_graph$Date, "%Y-%m")

datag <- pivot_longer(data = data_graph,cols = Skeletonema:Others,names_to = "Phylum")

datag$Phylum <- factor(datag$Phylum, levels = c("Asterionellopsis","Chaetocerotaceae","Cylindrotheca","Leptocylindrus"
                                                ,"Nitzschia","Pseudo-nitzschia", "Skeletonema","Akashiwo","Azadinium",
                                                "Chrysochromulina","Phaeocystis","Cryptophyceae","Others"))

ggplot(datag) +
  geom_col(aes(x = MonthYear, y = value*100, fill = Phylum), position = "stack", na.rm = FALSE, width = 1) +
  facet_wrap(~region, scales = "free", ncol = 1) +
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  geom_vline(data = subset(datag, format(Date, "%m") == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", size = 0.5,linetype = "dashed")+
  geom_text(data = subset(datag, format(Date, "%m") == "01"), 
            aes(x = MonthYear, y = 104, label = format(Date, "%Y")),
            color = "black", size = 2,angle=90, vjust = 0)+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 5))+
  labs(x="Date",y="Relative Abundance (%)",fill="Taxon")+
  scale_fill_manual(values = c("#FFBCFF", "#B2DF8A","#FF7F00", "#6A3D9A", "#33A02C", "#1F78B4","#A6CEE3",
                               "#FFFF99","#510051", "#BC00BC", "#FB9A99","#FDBF6F","grey"))
ggsave('Relativeabundance_genus.png', path = "output/graphs/description_region",dpi = 500, width = 360, height =200, units = 'mm')


### Test difference between regions for the networks properties ######

data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_PCA_coord_moments.final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

data <- select(data, Date, Code_point_Libelle, region, Dim.1, Dim.2, Dim.3, Shannon, Pielou, BergerParker)

kruskal.test(data$P_bac~data$region)
DunnTest(data$P_bac~data$region,method = "BH")
kruskal.test(data$P_dino~data$region)
DunnTest(data$P_dino~data$region,method = "BH")
kruskal.test(data$P_autres~data$region)
DunnTest(data$P_bac~data$region,method = "BH")

kruskal.test(data$P_BacBac~data$region)
DunnTest(data$P_BacBac~data$region,method = "BH")
kruskal.test(data$P_BacDino~data$region)
DunnTest(data$P_BacDino~data$region,method = "BH")
kruskal.test(data$P_DinoDino~data$region)
DunnTest(data$P_DinoDino~data$region,method = "BH")
kruskal.test(data$P_AAutres~data$region)
DunnTest(data$P_AAutres~data$region,method = "BH")


kruskal.test(data$Dim.1~data$region)
DunnTest(data$Dim.1~data$region,method = "BH")
kruskal.test(data$Dim.2~data$region)
DunnTest(data$Dim.2~data$region,method = "BH")
kruskal.test(data$Dim.3~data$region)
DunnTest(data$Dim.3~data$region,method = "BH")


kruskal.test(data$Pielou~data$region)
DunnTest(data$Pielou~data$region,method = "BH")
kruskal.test(data$BergerParker~data$region)
DunnTest(data$BergerParker~data$region,method = "BH")
kruskal.test(data$Shannon~data$region)
DunnTest(data$Shannon~data$region,method = "BH")


# graph

datal <- pivot_longer(data,names_to = "Var",cols = Dim.1:BergerParker)
BK <- ggplot(filter(datal,Var == "BergerParker"))+
  geom_boxplot(aes(y=value,group=region,fill=region),size = 1)+
  scale_fill_manual(values=region_col,guide = "none")+
  labs(x=NULL,y="")+
  facet_grid(Var~region,scales = "free_y")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

Pielou <- ggplot(filter(datal,Var == "Pielou"))+
  geom_boxplot(aes(y=value,group=region,fill=region),size = 1)+
  scale_fill_manual(values=region_col,guide = "none")+
  labs(x=NULL,y="")+
  facet_grid(Var~region,scales = "free_y")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

Shannon <- ggplot(filter(datal,Var == "Shannon"))+
  geom_boxplot(aes(y=value,group=region,fill=region),size = 1)+
  scale_fill_manual(values=region_col,guide = "none")+
  labs(x=NULL,y="")+
  facet_grid(Var~region,scales = "free_y")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

Dim.1 <- ggplot(filter(datal,Var == "Dim.1"))+
  geom_boxplot(aes(y=value,group=region,fill=region),size = 1)+
  scale_fill_manual(values=region_col,guide = "none")+
  labs(x=NULL,y="")+
  facet_grid(Var~region,scales = "free_y")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

Dim.2 <- ggplot(filter(datal,Var == "Dim.2"))+
  geom_boxplot(aes(y=value,group=region,fill=region),size = 1)+
  scale_fill_manual(values=region_col,guide = "none")+
  labs(x=NULL,y="")+
  facet_grid(Var~region,scales = "free_y")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

Dim.3 <- ggplot(filter(datal,Var == "Dim.3"))+
  geom_boxplot(aes(y=value,group=region,fill=region),size = 1)+
  scale_fill_manual(values=region_col,guide = "none")+
  labs(x=NULL,y="")+
  facet_grid(Var~region,scales = "free_y")+
  scale_y_continuous(limits = c(-5,7.5))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot_grid(Shannon,Pielou,BK,Dim.1,Dim.2,Dim.3, ncol = 3)
ggsave('difference_region_metrics.png', path = "output/graphs/description_region",dpi = 500, width = 460, height =200, units = 'mm')
