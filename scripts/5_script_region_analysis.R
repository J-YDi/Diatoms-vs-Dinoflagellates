################################################################################
# Diatoms vs dinoflagellates: a network analysis of bloom impacts on diversity #
#                    and phytoplankton associations | R scripts                #
################################################################################

# Script to investigate the differences between regions #
# 02/28/2025

# Load packages
library(ggplot2)
library(ggthemes)
library(readr)
library(dplyr)
library(tidyr)
library(FactoMineR)
library(factoextra)
library(cowplot)
library(DescTools)

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

# Don't consider the Pertuis Sea
data <- filter(data, region != "4-Pertuis Sea")

### Graph to show phytoplankton composition by season and cluster #####
# By Genus
# Compute the mean of all genus by season and cluster
data_graph <- data[,c(3,354,24:328)] %>%
  group_by(season, region) %>%
  summarise_all(~mean(., na.rm = TRUE))

# Make it relative 
data_graph$Abdtot <- rowSums(data_graph[,c(3:307)],na.rm=T)
data_graph[,c(3:307)] <- data_graph[,c(3:307)]/data_graph$Abdtot

datag <- pivot_longer(data = data_graph,cols = Actinoptychus:Coscinodiscophycidae,names_to = "Phylum")
datag$Abdtot <- NULL

# Group by season and region, and select the 4 most abundant phyla for each region
# Determine the most abundant phyla for each region and season

data_graph <- datag %>%
  group_by(season, region, Phylum) %>%
  summarise(abondance = sum(value,na.rm=T)) %>%
  ungroup() %>%
  group_by(season, region) %>%
  mutate(rank = rank(desc(abondance)))

#write.csv(data_graph,file="output/tableaux/phyla_mostabundant_region_season.csv")

# Select those 
selection_phylum <- c("Skeletonema","Pseudo-nitzschia","Chaetoceros","Chaetocerotaceae","Nitzschia","Cryptophyceae",
                      "Cryptomonadales","Cylindrotheca","Leptocylindrus","Akashiwo","Phaeocystis","Akashiwo",
                      "Phaeocystis","Asterionellopsis","Chrysochromulina","Azadinium")
datataxon <- data |>
  select(Month,Year,region,selection_phylum)
# Agglomerate Chaetoceros with Chaetocerotaceae and Cryptophyceae with Cryptomonadales
datataxon$Chaetocerotaceae <- rowSums(datataxon[,c("Chaetoceros","Chaetocerotaceae")],na.rm=T)
datataxon$Cryptophyceae <- rowSums(datataxon[,c("Cryptophyceae","Cryptomonadales")],na.rm=T)
datataxon <- select(datataxon,-c(Chaetoceros,Cryptomonadales))

# Create the "Other" category
dataothers <- data |>
  select(-selection_phylum) |>
  select(region,Actinoptychus:Coscinodiscophycidae)

dataothers$Others <- rowSums(select(dataothers,Actinoptychus:Coscinodiscophycidae),na.rm = T)
dataothers <- select(dataothers,Others)

# Bind them
datagraph <- cbind(datataxon,dataothers)

# Mean by month, region, year
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

# Graph by region
med <- ggplot(filter(datag,region=="1-Mediterranean sea")) +
  geom_col(aes(x = MonthYear, y = value*100, fill = Phylum), position = "stack", na.rm = FALSE, width = 1) +
  facet_wrap(~region, scales = "free", ncol = 1) +
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  geom_vline(data = subset(filter(datag,region=="1-Mediterranean sea"), format(Date, "%m") == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", size = 0.5,linetype = "dashed")+
  geom_text(data = subset(filter(datag,region=="1-Mediterranean sea"), format(Date, "%m") == "07"), 
            aes(x = MonthYear, y = 101, label = format(Date, "%Y")),
            color = "black", size = 3,angle=0, vjust = 0)+
  theme(axis.text.x = element_blank(),legend.position = "none")+
  labs(x=NULL,,y="")+
  scale_fill_manual(values = c("#FFBCFF", "#B2DF8A","#FF7F00", "#6A3D9A", "#33A02C", "#1F78B4","#A6CEE3",
                               "#FFFF99","#510051", "#BC00BC", "#FB9A99","#FDBF6F","grey"))

manche <- ggplot(filter(datag,region=="2-Eastern Channel - North Sea")) +
  geom_col(aes(x = MonthYear, y = value*100, fill = Phylum), position = "stack", na.rm = FALSE, width = 1) +
  facet_wrap(~region, scales = "free", ncol = 1) +
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  geom_vline(data = subset(filter(datag,region=="2-Eastern Channel - North Sea"), format(Date, "%m") == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", size = 0.5,linetype = "dashed")+
  #geom_text(data = subset(filter(datag,region=="2-Eastern Channel - North Sea"), format(Date, "%m") == "07"), 
  #          aes(x = MonthYear, y = 101, label = format(Date, "%Y")),
  #          color = "black", size = 3,angle=0, vjust = 0)+
  theme(axis.text.x = element_blank(),legend.position = "none")+
  labs(x=NULL,y="Relative abundance (%)")+
  scale_fill_manual(values = c("#FFBCFF", "#B2DF8A","#FF7F00", "#6A3D9A", "#33A02C", "#1F78B4","#A6CEE3",
                               "#FFFF99","#510051", "#BC00BC", "#FB9A99","#FDBF6F","grey"))

atl <- ggplot(filter(datag,region=="3-Atlantic - Western Channel")) +
  geom_col(aes(x = MonthYear, y = value*100, fill = Phylum), position = "stack", na.rm = FALSE, width = 1) +
  facet_wrap(~region, scales = "free", ncol = 1) +
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  geom_vline(data = subset(filter(datag,region=="3-Atlantic - Western Channel"), format(Date, "%m") == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", size = 0.5,linetype = "dashed")+
  #geom_text(data = subset(filter(datag,region=="3-Atlantic - Western Channel"), format(Date, "%m") == "07"), 
  #          aes(x = MonthYear, y = 101, label = format(Date, "%Y")),
  #          color = "black", size = 3,angle=0, vjust = 0)+
  theme(axis.text.x = element_text(size = 6),legend.position = "none")+
  labs(x="Date",y="")+
  scale_fill_manual(values = c("#FFBCFF", "#B2DF8A","#FF7F00", "#6A3D9A", "#33A02C", "#1F78B4","#A6CEE3",
                               "#FFFF99","#510051", "#BC00BC", "#FB9A99","#FDBF6F","grey"))
# Final graph:
plot_grid(med,manche,atl,nrow = 3,greedy = "FALSE" )
#ggsave('Relativeabundance_genus.png', path = "output/graphs/description_region",dpi = 500, width = 360, height =200, units = 'mm')




### Test difference between regions for the networks properties ######

# Import data
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_PCA_coord_moments.final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Select the parameters
data <- select(data, Date, Code_point_Libelle, region, Dim.1, Dim.2, Dim.3, Shannon, Pielou, BergerParker,
               P_bac,P_dino,P_autres,P_BacBac,P_BacDino,P_DinoDino,P_AAutres)

# Stastistical tests
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


# Viz
# Differences between station for diversity indexes and PCA's dimensions
colnames(data) <- c("Date","Code_point_Libelle","region","Dim.1","Dim.2","Dim.3","Shannon","Pielou","Berger-Parker",
                    "Bacillariophyceae","Dinophyceae","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso.")

datal <- pivot_longer(data,names_to = "Var",cols = Dim.1:`Berger-Parker`)
BK <- ggplot(filter(datal,Var == "Berger-Parker"))+
  geom_boxplot(aes(y=value,group=region,fill=region),size = 1)+
  scale_fill_manual(values=region_col,guide = "none")+
  labs(x=NULL,y=NULL)+
  facet_grid(Var~region,scales = "free_y")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

Pielou <- ggplot(filter(datal,Var == "Pielou"))+
  geom_boxplot(aes(y=value,group=region,fill=region),size = 1)+
  scale_fill_manual(values=region_col,guide = "none")+
  labs(x=NULL,y=NULL)+
  facet_grid(Var~region,scales = "free_y")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

Shannon <- ggplot(filter(datal,Var == "Shannon"))+
  geom_boxplot(aes(y=value,group=region,fill=region),size = 1)+
  scale_fill_manual(values=region_col,guide = "none")+
  labs(x=NULL,y="Value")+
  facet_grid(Var~region,scales = "free_y")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

Dim.1 <- ggplot(filter(datal,Var == "Dim.1"))+
  geom_boxplot(aes(y=value,group=region,fill=region),size = 1)+
  scale_fill_manual(values=region_col,guide = "none")+
  labs(x=NULL,y="Value")+
  facet_grid(Var~region,scales = "free_y")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

Dim.2 <- ggplot(filter(datal,Var == "Dim.2"))+
  geom_boxplot(aes(y=value,group=region,fill=region),size = 1)+
  scale_fill_manual(values=region_col,guide = "none")+
  labs(x=NULL,y=NULL)+
  facet_grid(Var~region,scales = "free_y")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

Dim.3 <- ggplot(filter(datal,Var == "Dim.3"))+
  geom_boxplot(aes(y=value,group=region,fill=region),size = 1)+
  scale_fill_manual(values=region_col,guide = "none")+
  labs(x=NULL,y=NULL)+
  facet_grid(Var~region,scales = "free_y")+
  scale_y_continuous(limits = c(-5,7.5))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot_grid(Shannon,Pielou,BK,Dim.1,Dim.2,Dim.3, ncol = 3)
#ggsave('difference_region_metrics.png', path = "output/graphs/description_region",dpi = 500, width = 460, height =200, units = 'mm')

# For association networks composition
datal <- pivot_longer(data,names_to = "Var",cols = c(Bacillariophyceae,Dinophyceae,`Other taxa`))
compo <- ggplot(filter(datal))+
  geom_boxplot(aes(y=value,group=region,fill=region),size = 1)+
  scale_fill_manual(values=region_col,guide = "none")+
  labs(x=NULL,y="Proportion")+
  facet_grid(Var~region,scales = "free_y")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

datal <- pivot_longer(data,names_to = "Var",cols = c(`Bac-Bac`,`Bac-Dino`,`Dino-Dino`))
asso <- ggplot(filter(datal))+
  geom_boxplot(aes(y=value,group=region,fill=region),size = 1)+
  scale_fill_manual(values=region_col,guide = "none")+
  labs(x=NULL,y=NULL)+
  facet_grid(Var~region,scales = "free_y")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot_grid(compo,asso)
#ggsave('difference_region_compo_asso.png', path = "output/graphs/description_region",dpi = 500, width = 360, height =200, units = 'mm')
