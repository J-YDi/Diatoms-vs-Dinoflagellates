################################################################################
# Diatoms vs dinoflagellates: a network analysis of bloom impacts on diversity #
#                    and phytoplankton associations | R scripts                #
################################################################################

# Script used to manipulate original datasets and create new ones #
# 02/26/2026

# Load packages
library(dplyr)
library(tidyverse)
library(stringr)
library(FactoMineR)
library(factoextra)
library(vegan)
library(paletteer)
library(mgcv)
library(gratia)
library(scatterpie)
library(OTUtable)
library(vegan)
library(readr)

#### Import data ###### 

# Open the original data 
DataREPHY_MA <- read.csv2('data/REPHY_Manche_Atlantique_1987-2022.csv', fileEncoding = "ISO-8859-1")
DataREPHY_Med <- read.csv2('data/REPHY_Med_1987-2022.csv', fileEncoding = "ISO-8859-1")

# Creating a vector to import all columns as characters
classes_char <- rep('character', 56)

# Load 2023 data
DataREPHY_2023 <- read.csv2('data/Extraction SEANOE_REPHY_phyto-Hydro_Manche_Atl-Med 2023 Semestre1 validé 16012024.csv', 
                            fileEncoding = "ISO-8859-1", colClasses = classes_char)

# Merge the first 2 datasets
DataREPHY <- bind_rows(DataREPHY_MA, DataREPHY_Med) %>%
  # This line allows to remove the problematic row that separates the hydro and phyto datasets
  filter(Passage...Mois != 'Passage : Mois')

# Binding 2023 to the other years 
DataREPHY_8723 <- bind_rows(DataREPHY, DataREPHY_2023)


### Load tables used to enrich the dataset
# Load table "Zones_marines"
ZM <- read.csv('data/Zones_marines.csv', sep = ';', header = TRUE)

# Load table "Liste_phylum.classe" to associate species to a class
PhyClasse <- read.csv('data/Liste_phylum.classe_REPHY.csv', sep =';', header = TRUE, fileEncoding = 'ISO-8859-1')

#### Formatting the dataframe with better column names

# Extracting the numeric code for the ZM
Table1 <- DataREPHY_8723 %>%
  mutate(ZM_Quadrige_Numero = as.numeric(str_extract(Lieu.de.surveillance...Entité.de.classement...Libellé, '[:alnum:]+')))

# Change column names (to match ZM)
colnames(Table1)[which(names(Table1) == "Lieu.de.surveillance...Mnémonique")] <- "Code_point_Mnemonique"
colnames(Table1)[which(names(Table1) == "Lieu.de.surveillance...Libellé")] <- "Code_point_Libelle"
colnames(Table1)[which(names(Table1) == "Passage...Date")] <- "Date"
colnames(Table1)[which(names(Table1) == "Coordonnées.passage...Coordonnées.minx")] <- "lon"
colnames(Table1)[which(names(Table1) == "Coordonnées.passage...Coordonnées.miny")] <- "lat"
colnames(Table1)[which(names(Table1) == "Résultat...Libellé.unité.de.mesure.associé.au.quintuplet")] <- "Mesure_Unite"
colnames(Table1)[which(names(Table1) == "Résultat...Symbole.unité.de.mesure.associé.au.quintuplet")] <- "Mesure_Symbole"
colnames(Table1)[which(names(Table1) == "Résultat...Nom.du.taxon.référent")] <- "Taxon"
colnames(Table1)[which(names(Table1) == "Résultat...Libellé.du.groupe.de.taxon")] <- "Groupe_Taxon"
colnames(Table1)[which(names(Table1) == "Résultat...Valeur.de.la.mesure")] <- "Valeur_mesure"
colnames(Table1)[which(names(Table1) == "Prélèvement...Immersion")] <- "Profondeur.metre"
colnames(Table1)[which(names(Table1) == "Prélèvement...Niveau")] <- "Prelevement.niveau"
colnames(Table1)[which(names(Table1) == "Résultat...Code.paramètre")] <- "Code.parametre"
colnames(Table1)[which(names(Table1) == "Résultat...Libellé.paramètre")] <- "Parametre"
colnames(Table1)[which(names(Table1) == "Résultat...Niveau.de.qualité")] <- "Qualite.resultat"
colnames(Table1)[which(names(Table1) == "Prélèvement...Niveau.de.qualité")] <- "Qualite.prelevement"
colnames(Table1)[which(names(Table1) == "Passage...Service.saisisseur...Libellé")] <- "Service.saisie"
colnames(Table1)[which(names(Table1) == "Passage...Heure")] <- "Heure"
colnames(Table1)[which(names(Table1) == "Résultat...Service.analyste...Code")] <- "Service.analyse"
colnames(Table1)[which(names(Table1) == "Prélèvement...Service.préleveur...Code")] <- "Service.prelevement"
colnames(Table1)[which(names(Table1) == "Prélèvement...Identifiant.interne")] <- "ID.interne.prelevement"
colnames(Table1)[which(names(Table1) == "Passage...Identifiant.interne")] <- "ID.interne.passage"

# Convert some oxygene measurements in ml/l to mg/l as recommanded in the REPHY manual
# Replace , by . to make it numeric
Table1$Valeur_mesure <-  str_replace_all(Table1$Valeur_mesure, ',', '.')
Table1 <- Table1 %>%
  mutate(Valeur_mesure = as.numeric(Valeur_mesure))
# Line for the conversion
Table1$Valeur_mesure <- ifelse(Table1$Code.parametre == "OXYGENE" & Table1$Mesure_Symbole == "ml.l-1", Table1$Valeur_mesure * 1.429 , Table1$Valeur_mesure)
# Indicating the change
Table1$Mesure_Symbole <- ifelse(Table1$Code.parametre == "OXYGENE" & Table1$Mesure_Symbole == "ml.l-1", "mg.l-1 converted" , Table1$Mesure_Symbole)


#### Curate table to keep only desired variables
Table1 <- Table1 %>%
  dplyr::select(c('ZM_Quadrige_Numero', 'Code_point_Mnemonique', 'Code_point_Libelle', 'Date', 
                  'Heure', 'lon', 'lat', 'Mesure_Unite', 'Mesure_Symbole', 'Taxon', 'Valeur_mesure', 'Résultat...Libellé.méthode',
                  'Prelevement.niveau', 'Profondeur.metre', 'Code.parametre', 'Parametre', 
                  'Qualite.prelevement', 'Qualite.resultat', 'ID.interne.prelevement', 'ID.interne.passage'))

# Modifying date format so that it gives the year and month, and getting rid of rows with no Year value
Table1 <- Table1 %>%
  mutate(Date = dmy(Date)) %>%
  mutate(Day = day(Date)) %>%
  mutate(Month = month(Date, label = F)) %>%
  mutate(Year = year(Date)) %>%
  filter(!is.na(Year))

#### Tidying table structure 
## Associate a region with each ZM code
Table1 <- left_join(Table1, ZM, by='ZM_Quadrige_Numero', suffix=c('',''))

# Separate the table for hydrology and phytoplankton
Table1_hydro <- Table1 %>%
  filter(Taxon == "")

Table1_phyto <- Table1 %>%
  filter(Taxon != "")

# Hydrology table
# Keep only sample between 0 and 5m depth
Table1_hydro$Profondeur.metre <- as.numeric(Table1_hydro$Profondeur.metre)
Table1_hydro <- Table1_hydro |> filter(Profondeur.metre <= 5 | is.na(Profondeur.metre)) |> #Some depth are not indicated
  filter(Year >= 1995) 
Table1_hydro_b <- Table1_hydro
# Keep Chlorophylle wanted, the sample measures with recommended methods
Table1_hydro <- filter(Table1_hydro, Code.parametre != "CHLOROA")
Table1_hydro_chloro <- Table1_hydro_b |>
  filter(Code.parametre == "CHLOROA") |>
  filter(Résultat...Libellé.méthode == "Chromatographie liquide - pigments phytoplanctoniques (Van Heukelem et Thomas 2001)" |
           Résultat...Libellé.méthode == "Chromatographie liquide - pigments phytoplanctoniques (Zapata et al. 2000)" |
           Résultat...Libellé.méthode == "Fluorimétrie (Aminot A. Kérouel R. 2004 - Chlorophylle)" |
           Résultat...Libellé.méthode == "Fluorimétrie (Neveux J. et Panouse M. 1987  - Chlorophylle)" |
           Résultat...Libellé.méthode == "Spectrophotométrie monochromatique (Aminot A. Kérouel R. 2004 - Chlorophylle)" |
           Résultat...Libellé.méthode == "Spectrophotométrie monochromatique (Aminot et Chaussepied 1983 - Chlorophylle)" )

Table1_hydro_chloro <- dplyr::select(Table1_hydro_chloro, Date, ID.interne.passage,
                                     Prelevement.niveau,Valeur_mesure,Résultat...Libellé.méthode)

colnames(Table1_hydro_chloro)[4:5] <- c("CHLOROA","Methode_CHLOROA")

# Solve the problem of duplicates due to double measurement of chloro with different methods
doublons_chloro <- Table1_hydro_chloro[duplicated(Table1_hydro_chloro$ID.interne.passage) |
                                         duplicated(Table1_hydro_chloro$ID.interne.passage, fromLast = TRUE), ]

# If there are duplicates we want to select one particular method (so measure) according to their global frequency
resultat_filtre_chloro <- doublons_chloro %>%
  filter(Methode_CHLOROA %in% levels(as.factor(doublons_chloro$Methode_CHLOROA))) %>%
  group_by(ID.interne.passage) %>%
  mutate(Ordre = match(Methode_CHLOROA, c("Spectrophotométrie monochromatique (Aminot A. Kérouel R. 2004 - Chlorophylle)",
                                          "Fluorimétrie (Aminot A. Kérouel R. 2004 - Chlorophylle)",
                                          "Chromatographie liquide - pigments phytoplanctoniques (Van Heukelem et Thomas 2001)",
                                          "Spectrophotométrie monochromatique (Aminot et Chaussepied 1983 - Chlorophylle)",
                                          "Fluorimétrie (Neveux J. et Panouse M. 1987  - Chlorophylle)",
                                          "Chromatographie liquide - pigments phytoplanctoniques (Zapata et al. 2000)"))) %>%
  arrange(desc(Ordre)) %>%
  filter(duplicated(ID.interne.passage) | n()==1)

Table1_hydro_chloro_unique <- subset(Table1_hydro_chloro, !(ID.interne.passage %in% unique(doublons_chloro$ID.interne.passage)))

Table1_hydro_chloro <- bind_rows(Table1_hydro_chloro_unique,resultat_filtre_chloro)

# If there are different measures at different depth to a same place and date we want to select the depth-less 
Table1_hydro <- Table1_hydro |>
  filter(Code.Region != 0) %>%
  filter(Prelevement.niveau == "Surface (0-1m)" |Prelevement.niveau == "2 mètres" |Prelevement.niveau == "de 3 à 5 mètres" |Prelevement.niveau == "Mi-profondeur" ) |>
  group_by(Code.Region, Code_point_Libelle, lon, lat, Year, Month, Date, ID.interne.passage, Prelevement.niveau, Profondeur.metre, Code.parametre) %>%
  summarise(Valeur_mesure = mean(Valeur_mesure), .groups = 'keep') %>%
  pivot_wider(names_from = Code.parametre, values_from = Valeur_mesure)

# Select depth
doublons_hydro <- Table1_hydro[duplicated(Table1_hydro$ID.interne.passage) |
                                 duplicated(Table1_hydro$ID.interne.passage, fromLast = TRUE), ]

# Solve duplicates in hydro :
resultat_filtre <- doublons_hydro %>%
  filter(Prelevement.niveau %in% c("Surface (0-1m)", "2 mètres", "de 3 à 5 mètres","Mi-profondeur")) %>%
  group_by(ID.interne.passage) %>%
  mutate(Ordre = match(Prelevement.niveau, c("Surface (0-1m)", "2 mètres", "de 3 à 5 mètres","Mi-profondeur"))) %>%
  arrange(desc(Ordre)) %>%
  filter(duplicated(ID.interne.passage) | n()==1)

Table1_unique <- subset(Table1_hydro, !(ID.interne.passage %in% unique(doublons_hydro$ID.interne.passage)))

Table1_hydro <- bind_rows(Table1_unique,resultat_filtre)

Table1_hydro <- Table1_hydro |>
  group_by(Code.Region, Code_point_Libelle, lon, lat, Year, Month, Date, ID.interne.passage, Prelevement.niveau,Profondeur.metre)


# Phytoplankton data 
# As for hydrological data : keep Depth between 0 et 5m and the depth-less if multiple ones
Table1_phyto$Profondeur.metre <- as.numeric(Table1_phyto$Profondeur.metre)
Table1_phyto <- Table1_phyto |> filter(Profondeur.metre <= 5 | is.na(Profondeur.metre)) |>
  filter(Year >= 1995) |>
  filter(Code.parametre == "FLORTOT")

Table1_phyto <- Table1_phyto |>
  filter(Code.Region != 0) %>%
  filter(Prelevement.niveau == "Surface (0-1m)" |Prelevement.niveau == "2 mètres" |Prelevement.niveau == "de 3 à 5 mètres" |Prelevement.niveau == "Mi-profondeur" ) |>
  group_by(Code.Region, Code_point_Libelle, lon, lat, Year, Month, Date, ID.interne.passage, Prelevement.niveau,Profondeur.metre, Code.parametre)

### Manipulating the phyto table ###
# Associate the Phylum/Class to the taxon
Table1_phyto <- left_join(Table1_phyto, PhyClasse, by='Taxon', suffix=c('',''))

# Compute genus from Taxon
# We find the genus by selecting only the first word in the string 'Taxon'
Table1_phyto <- Table1_phyto %>%
  mutate(Genus = str_extract(Taxon, 'Pseudo-nitzschia|[:alnum:]+'))

# Group the abundances to the genus level
Table1_phyto_genus <- Table1_phyto %>%
  group_by(Code.Region, Code_point_Libelle, lon, lat, Year, Month, Date, ID.interne.passage, Prelevement.niveau, Profondeur.metre,Code.parametre, Genus) %>%
  summarise(Comptage = sum(Valeur_mesure), .groups = 'keep') %>%
  pivot_wider(names_from = Genus, values_from = Comptage)

Table1_phyto_class <- Table1_phyto[Table1_phyto$Phylum.Classe != "", ]

Table1_phyto_classe <- Table1_phyto_class %>%
  group_by(Code.Region, Code_point_Libelle, lon, lat, Year, Month, Date, ID.interne.passage, Prelevement.niveau, Profondeur.metre,Code.parametre, Phylum.Classe) %>%
  summarise(Comptage = sum(Valeur_mesure), .groups = 'keep') %>%
  pivot_wider(names_from = Phylum.Classe, values_from = Comptage)

# Select depth
doublons_phyto_genus <- Table1_phyto_genus[duplicated(Table1_phyto_genus$ID.interne.passage) |
                                             duplicated(Table1_phyto_genus$ID.interne.passage, fromLast = TRUE), ]

resultat_filtre_genus <- doublons_phyto_genus %>%
  filter(Prelevement.niveau %in% c("Surface (0-1m)", "2 mètres", "de 3 à 5 mètres","Mi-profondeur")) %>%
  group_by(ID.interne.passage) %>%
  mutate(Ordre = match(Prelevement.niveau, c("Surface (0-1m)", "2 mètres", "de 3 à 5 mètres","Mi-profondeur"))) %>%
  arrange(desc(Ordre)) %>%
  filter(duplicated(ID.interne.passage) | n()==1)


Table1_phyto_unique <- subset(Table1_phyto_genus, !(ID.interne.passage %in% unique(doublons_phyto_genus$ID.interne.passage)))

Table1_phyto_genus <- bind_rows(Table1_phyto_unique,resultat_filtre_genus)

# Final phyto data
Table1_phyto_genus <- Table1_phyto_genus |>
  group_by(Code.Region, Code_point_Libelle, lon, lat, Year, Month, Date, ID.interne.passage, Prelevement.niveau)

# Prepare merging phyto and hydro data
Table1_hydro <- dplyr::select(ungroup(Table1_hydro),Date:Prelevement.niveau,SALI:PCYAN )

data_hp <- left_join(Table1_phyto_genus,Table1_hydro,by = join_by(Date, ID.interne.passage,
                                                                  Prelevement.niveau))

Table1_phyto_classe <- select(ungroup(Table1_phyto_classe),Date, ID.interne.passage, Prelevement.niveau,Bacillariophyceae:Xanthophyceae)

data_hpc <- left_join(data_hp,Table1_phyto_classe,by = join_by(Date, ID.interne.passage,
                                                               Prelevement.niveau))

data_hpc <- dplyr::select(data_hpc, Code.Region:Code.parametre, SALI:`NO3+NO2`,`TURB-FNU`, `Actinoptychus`:`Coscinodiscophycidae`,Bacillariophyceae.y:Xanthophyceae)

# Select chlorophyll methods
data_hpc2 <- left_join(data_hpc, Table1_hydro_chloro)

data_hpc2 <- dplyr::select(data_hpc2, Code.Region:Code.parametre,Methode_CHLOROA,CHLOROA, SALI:`NO3+NO2`, `TURB-FNU`, `Actinoptychus`:`Coscinodiscophycidae`,Bacillariophyceae.y:Xanthophyceae)

# The missing depths correspond only to the "Surface (0-1m)" so we say that corresponds
# has a depth of 0.5 the Depth in NA
data_hpc2$Profondeur.metre <- ifelse(is.na(data_hpc2$Profondeur.metre),0.5,data_hpc2$Profondeur.metre) 

data_hpc2 <- data_hpc2 |>
  group_by(Code.Region, Code_point_Libelle, lon, lat, Year, Month, Date, ID.interne.passage, Prelevement.niveau)

colnames(data_hpc2)[c(328,329,331,333,336,337,338,339,341,342)] <- c("Bacillariophyceae","Dinophyceae","Raphidophyceae","Dictyochophyceae","Ciliophora",
                                                                     "Cryptophyceae","Dinoflagellata","Chrysophyceae","Chromista","Cyanophyceae")

# Still some duplicates
doublons_final <- data_hpc2[duplicated(data_hpc2$ID.interne.passage) |
                              duplicated(data_hpc2$ID.interne.passage, fromLast = TRUE), ]

resultat_filtre_final <- doublons_final %>%
  filter(Prelevement.niveau %in% c("Surface (0-1m)", "2 mètres", "de 3 à 5 mètres","Mi-profondeur")) %>%
  group_by(ID.interne.passage) %>%
  mutate(Ordre = match(Prelevement.niveau, c("Surface (0-1m)", "2 mètres", "de 3 à 5 mètres","Mi-profondeur"))) %>%
  arrange(desc(Ordre)) %>%
  filter(duplicated(ID.interne.passage) | n()==1)

data_hpc2_unique <- subset(data_hpc2, !(ID.interne.passage %in% unique(doublons_final$ID.interne.passage)))
data_hpc3 <- bind_rows(data_hpc2_unique,resultat_filtre_final)
data_hpc3 <- data_hpc3 |>
  group_by(Code.Region, Code_point_Libelle, lon, lat, Year, Month, Date, ID.interne.passage, Prelevement.niveau)

data_hpc3 <- select(data_hpc3, -Ordre)

# Save it
#write.csv2(data_hpc3,file="output/data_modif/Table_FLORTOT_Surf_9523_hydro_phyto_chloro.csv", row.names = FALSE,dec = ".")

########## STATION SELECTION AFTER SAMPLING ANALYSIS ########
Table <- read_delim("output/data_modif/Table_FLORTOT_Surf_9523_hydro_phyto_chloro.csv", 
                    delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                        grouping_mark = ""), trim_ws = TRUE)

# Selection with period by station of interest on the criterion of continuity of chlorophyll 5 YEARS criteria
# Create the dataframe for spatial clustering

Table_select <- filter(Table, Code_point_Libelle == "Point 1 Dunkerque" & (Date >= as.Date("2002-04-01") & Date <= as.Date("2012-11-30"))|
                         Code_point_Libelle == "Point 1 Boulogne" & (Date >= as.Date("1995-01-01") & Date <= as.Date("2022-12-31"))|
                         Code_point_Libelle == "Bif" & (Date >= as.Date("2009-08-01") & Date <= as.Date("2017-06-30"))|
                         Code_point_Libelle == "At so" & (Date >= as.Date("2003-03-01") & Date <= as.Date("2022-12-31"))|
                         Code_point_Libelle == "Sète mer" & (Date >= as.Date("1995-01-01") & Date <= as.Date("2023-01-31"))|
                         Code_point_Libelle == "Parc Leucate 2" & (Date >= as.Date("1995-02-01") & Date <= as.Date("2022-08-31"))|
                         Code_point_Libelle == "Grand Rhône" & (Date >= as.Date("2008-01-01") & Date <= as.Date("2015-08-31"))|
                         Code_point_Libelle == "Bouzigues (a)" & (Date >= as.Date("1995-01-01") & Date <= as.Date("2023-02-28"))|
                         Code_point_Libelle == "Barcares" & (Date >= as.Date("1995-01-01") & Date <= as.Date("2022-10-31"))|
                         Code_point_Libelle == "Anse de Carteau 2" & (Date >= as.Date("2003-02-01") & Date <= as.Date("2023-06-30"))|
                         Code_point_Libelle == "Géfosse" & (Date >= as.Date("2004-01-01") & Date <= as.Date("2022-12-31"))|
                         Code_point_Libelle == "Cabourg" & (Date >= as.Date("2002-03-01") & Date <= as.Date("2022-12-31"))|
                         Code_point_Libelle == "Antifer ponton pétrolier" & (Date >= as.Date("2002-04-01") & Date <= as.Date("2022-12-31"))|
                         Code_point_Libelle == "Villefranche" & (Date >= as.Date("1998-10-01") & Date <= as.Date("2014-04-30"))| # A VOIR
                         Code_point_Libelle == "Sud Bastia" & (Date >= as.Date("2007-04-01") & Date <= as.Date("2020-02-29"))|
                         Code_point_Libelle == "Lazaret (a)" & (Date >= as.Date("1999-01-01") & Date <= as.Date("2015-09-30"))|
                         Code_point_Libelle == "Diana centre" & (Date >= as.Date("2001-02-01") & Date <= as.Date("2022-11-30"))|
                         Code_point_Libelle == "Calvi" & (Date >= as.Date("2006-02-01") & Date <= as.Date("2022-08-31"))|
                         Code_point_Libelle == "22B - Toulon gde rade" & (Date >= as.Date("2006-02-01") & Date <= as.Date("2023-06-30"))|
                         Code_point_Libelle == "Ouest Loscolo" & (Date >= as.Date("1996-01-01") & Date <= as.Date("2022-12-31"))|
                         Code_point_Libelle == "Nord Saumonards" & (Date >= as.Date("2009-03-01") & Date <= as.Date("2023-04-30"))|
                         Code_point_Libelle == "Le Croisic (a)" & (Date >= as.Date("2000-02-01") & Date <= as.Date("2015-12-31"))|
                         Code_point_Libelle == "Le Cornard" & (Date >= as.Date("1995-01-01") & Date <= as.Date("2023-05-31"))|
                         Code_point_Libelle == "La Palmyre" & (Date >= as.Date("2011-01-01") & Date <= as.Date("2020-03-31"))|
                         Code_point_Libelle == "Filière w" & (Date >= as.Date("2012-03-01") & Date <= as.Date("2023-05-31"))|
                         Code_point_Libelle == "Boyard" & (Date >= as.Date("1995-01-01") & Date <= as.Date("2015-07-30"))|
                         Code_point_Libelle == "Bois de la Chaise large" & (Date >= as.Date("2007-03-01") & Date <= as.Date("2023-06-30"))|
                         Code_point_Libelle == "Basse Michaud" & (Date >= as.Date("2016-01-01") & Date <= as.Date("2023-06-30"))|
                         Code_point_Libelle == "Auger" & (Date >= as.Date("1995-04-01") & Date <= as.Date("2023-03-31"))|
                         Code_point_Libelle == "Teychan bis" & (Date >= as.Date("1995-01-01") & Date <= as.Date("2023-06-30"))|
                         Code_point_Libelle == "Arcachon - Bouée 7" & (Date >= as.Date("2014-10-01") & Date <= as.Date("2023-06-30"))|
                         Code_point_Libelle == "Ouessant - Youc'h korz" & (Date >= as.Date("2016-01-01") & Date <= as.Date("2022-01-31"))|
                         Code_point_Libelle == "Men er Roue" & (Date >= as.Date("1996-01-01") & Date <= as.Date("2022-12-31"))|
                         Code_point_Libelle == "Men Du" & (Date >= as.Date("1995-02-01") & Date <= as.Date("2003-12-31"))|
                         Code_point_Libelle == "Le Passage (a)" & (Date >= as.Date("1995-01-01") & Date <= as.Date("2000-12-31"))|
                         Code_point_Libelle == "Lanvéoc large" & (Date >= as.Date("2010-01-01") & Date <= as.Date("2020-08-31"))|
                         Code_point_Libelle == "Kervel large" & (Date >= as.Date("2013-11-01") & Date <= as.Date("2023-06-30"))|
                         Code_point_Libelle == "Concarneau large" & (Date >= as.Date("2010-01-01") & Date <= as.Date("2023-06-30"))|
                         Code_point_Libelle == "St Pol large" & (Date >= as.Date("2013-01-01") & Date <= as.Date("2023-06-30"))|
                         Code_point_Libelle == "St Cast" & (Date >= as.Date("1995-01-01") & Date <= as.Date("2004-03-31"))|
                         Code_point_Libelle == "Pen al Lann" & (Date >= as.Date("1995-01-01") & Date <= as.Date("2006-02-28"))|
                         Code_point_Libelle == "Loguivy" & (Date >= as.Date("2007-03-01") & Date <= as.Date("2023-06-30"))|
                         Code_point_Libelle == "les Hébihens" & (Date >= as.Date("2007-03-01") & Date <= as.Date("2023-06-30"))|
                         Code_point_Libelle == "Donville" & (Date >= as.Date("2002-04-01") & Date <= as.Date("2014-08-31"))|
                         Code_point_Libelle == "Bréhat" & (Date >= as.Date("1995-01-01") & Date <= as.Date("2003-06-30"))|
                         Code_point_Libelle == "Marseillan (a)" & (Date >= as.Date("2009-03-01") & Date <= as.Date("2023-02-28"))
)

# Save it
#write.csv2(Table_select,file="output/data_modif/Table_FLORTOT_Surf_9523_Stselect_hydro_phyto_chloro_phylum_period5_final.csv", row.names = FALSE,dec = ".")



# Selection with period by station interest on chlorophyll continuity 15 YEARS criteria
# Create the dataframe for analysis
Table_select <- filter(Table, Code_point_Libelle == "Point 1 Boulogne" & (Date >= as.Date("1995-01-01") & Date <= as.Date("2022-12-31"))|
                         Code_point_Libelle == "At so" & (Date >= as.Date("2003-03-01") & Date <= as.Date("2022-12-31"))|
                         Code_point_Libelle == "Sète mer" & (Date >= as.Date("1995-01-01") & Date <= as.Date("2023-01-31"))|
                         Code_point_Libelle == "Parc Leucate 2" & (Date >= as.Date("1995-02-01") & Date <= as.Date("2022-08-31"))|
                         Code_point_Libelle == "Bouzigues (a)" & (Date >= as.Date("1995-01-01") & Date <= as.Date("2023-02-28"))|
                         Code_point_Libelle == "Barcares" & (Date >= as.Date("1995-01-01") & Date <= as.Date("2022-10-31"))|
                         Code_point_Libelle == "Anse de Carteau 2" & (Date >= as.Date("2003-02-01") & Date <= as.Date("2023-06-30"))|
                         Code_point_Libelle == "Géfosse" & (Date >= as.Date("2004-01-01") & Date <= as.Date("2022-12-31"))|
                         Code_point_Libelle == "Cabourg" & (Date >= as.Date("2002-03-01") & Date <= as.Date("2022-12-31"))|
                         Code_point_Libelle == "Antifer ponton pétrolier" & (Date >= as.Date("2002-04-01") & Date <= as.Date("2022-12-31"))|
                         Code_point_Libelle == "Villefranche" & (Date >= as.Date("1998-10-01") & Date <= as.Date("2014-04-30"))| # A VOIR
                         Code_point_Libelle == "Lazaret (a)" & (Date >= as.Date("1999-01-01") & Date <= as.Date("2015-09-30"))|
                         Code_point_Libelle == "Diana centre" & (Date >= as.Date("2001-02-01") & Date <= as.Date("2022-11-30"))|
                         Code_point_Libelle == "Calvi" & (Date >= as.Date("2006-02-01") & Date <= as.Date("2022-08-31"))|
                         Code_point_Libelle == "22B - Toulon gde rade" & (Date >= as.Date("2006-02-01") & Date <= as.Date("2023-06-30"))|
                         Code_point_Libelle == "Ouest Loscolo" & (Date >= as.Date("1996-01-01") & Date <= as.Date("2022-12-31"))|
                         Code_point_Libelle == "Le Croisic (a)" & (Date >= as.Date("2000-02-01") & Date <= as.Date("2015-12-31"))|
                         Code_point_Libelle == "Le Cornard" & (Date >= as.Date("1995-01-01") & Date <= as.Date("2023-05-31"))|
                         Code_point_Libelle == "Boyard" & (Date >= as.Date("1995-01-01") & Date <= as.Date("2015-07-30"))|
                         Code_point_Libelle == "Bois de la Chaise large" & (Date >= as.Date("2007-03-01") & Date <= as.Date("2023-06-30"))|
                         Code_point_Libelle == "Auger" & (Date >= as.Date("1995-04-01") & Date <= as.Date("2023-03-31"))|
                         Code_point_Libelle == "Teychan bis" & (Date >= as.Date("1995-01-01") & Date <= as.Date("2023-06-30"))|
                         Code_point_Libelle == "Men er Roue" & (Date >= as.Date("1996-01-01") & Date <= as.Date("2022-12-31"))|
                         Code_point_Libelle == "Loguivy" & (Date >= as.Date("2007-03-01") & Date <= as.Date("2023-06-30"))|
                         Code_point_Libelle == "les Hébihens" & (Date >= as.Date("2007-03-01") & Date <= as.Date("2023-06-30"))
)

# Keep the stations where the 15 years are common
Table_select <- filter(Table_select, Code_point_Libelle %in% c("Teychan bis" ,"Ouest Loscolo" , "Men er Roue" ,
                        "Le Cornard" , "Bois de la Chaise large" ,"Auger" ,"Point 1 Boulogne" , "Loguivy" ,
                        "les Hébihens" ,"Géfosse" , "Cabourg" ,"At so" ,"Antifer ponton pétrolier" ,"Sète mer" , "Parc Leucate 2",
                        "Diana centre" , "Calvi" , "Bouzigues (a)" , "Barcares" ,"Anse de Carteau 2","22B - Toulon gde rade"))
Table_select <- filter(Table_select, Date >= as.Date("2007-03-01") & Date <= as.Date("2022-08-31"))

########## Add diversity indexes #####

# Store the old data for the manipulation
TableIDiv <- Table_select
# Select all the genus
TableIDiv <- dplyr::select(TableIDiv, ID.interne.passage, Asterionellopsis:Coscinodiscophycidae)

# Create a dataframe to store the diversity indexes by date and station
Data_results <- as.data.frame(c("",""))
Data_results[1,1] <- 0
Data_results[1,2] <- 0

# Calculate Shannon, Simpson, Berger-Parker and Pielou indexes
for (i in 1:7588){
  Table2 <- TableIDiv[i,]
  abondance_sums <- apply(Table2[,2:305],2,sum,na.rm=T)
  #print(i/7588 * 100)
  Data_results[i,1] <- diversity(abondance_sums,index = "shannon")
  Data_results[i,2] <- diversity(abondance_sums,index = "simpson")
  Data_results[i,3] <- max(abondance_sums) / sum(abondance_sums)
  Data_results[i,4] <- pielou(sample = abondance_sums)
  Data_results[i,5] <- TableIDiv[i,1]
}
colnames(Data_results) <- c("Shannon","Simpson","BergerParker","Pielou","ID.interne.passage")

Data_results$Shannon <- as.numeric(Data_results$Shannon)

data_ok <- left_join(Table_select,Data_results)
# Calculate the species richness
data_ok$Rspe <- rowSums(TableIDiv[,-1] != 0,na.rm = T)
# Delete the sampling where phytoplankton was not counted
data_ok <- filter(data_ok, Rspe != 0 )

# Save it
#write.csv2(data_ok,file="output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_div_final.csv", row.names = FALSE,dec = ".")


#####################################################################################
#                                                                                   #
#                 NEED TO RUN 1_Script_clustering_stations.R                        #
#                                                                                   #
#####################################################################################

# Load data
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_div_final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Load clustering data information
data_clust <- read_delim("output/tableaux/clusters_EM_01_k5_final.csv", 
                         delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                             grouping_mark = ""), trim_ws = TRUE)
# Change the clusters' name
data_clust$cluster <- ifelse(data_clust$cluster == 1, "1-Mediterranean sea",
                        ifelse(data_clust$cluster == 2, "2-Eastern Channel - North Sea",
                            ifelse(data_clust$cluster == 3, "3-Atlantic - Western Channel",
                                   ifelse(data_clust$cluster == 4, "4-Pertuis Sea",
                                          ifelse(data_clust$cluster == 5, "5-Rhone Estuary",data_clust$cluster)))))

# Binding the two datasets
data_clusterised <- left_join(data,data_clust)
# Arrange columns as I want
data_ok <- dplyr::select(data_clusterised,Code.Region:Code_point_Libelle,cluster,lon:Rspe)

colnames(data_ok)[3] <- "region"

# Save it
#write.csv2(data_ok,file="output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_final.csv", row.names = FALSE,dec = ".")


#####################################################################################
#                                                                                   #
#                 NEED TO RUN 2_script_build_networks.R                             #
#                                                                                   #
#####################################################################################
# Import the data with the network metrics
# Mediterranean sea
data_med <- read_delim("output/tableaux/Networks/subnetworks/results_metrics_reseaux_cluster1.csv", 
                         delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                             grouping_mark = ""), trim_ws = TRUE)
# Eastern Channel
data_manche <- read_delim("output/tableaux/Networks/subnetworks/results_metrics_reseaux_cluster2.csv", 
                       delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                           grouping_mark = ""), trim_ws = TRUE)
# Atlantic
data_atlantic <- read_delim("output/tableaux/Networks/subnetworks/results_metrics_reseaux_cluster3.csv", 
                       delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                           grouping_mark = ""), trim_ws = TRUE)
# Combine all
data_metrics <- rbind(data_med,data_manche,data_atlantic)

# Import data without this metrics
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_blooms_caracterised_final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)
# Bind them
datawithmetrics <- left_join(data,data_metrics)

# Still some duplicates to delete
doublons_final <- datawithmetrics[duplicated(datawithmetrics$ID.interne.passage) |
                              duplicated(datawithmetrics$ID.interne.passage, fromLast = TRUE), ]

resultat_filtre_final <- doublons_final %>%
  filter(Prelevement.niveau %in% c("Surface (0-1m)", "2 mètres", "de 3 à 5 mètres","Mi-profondeur")) %>%
  group_by(ID.interne.passage) %>%
  mutate(Ordre = match(Prelevement.niveau, c("Surface (0-1m)", "2 mètres", "de 3 à 5 mètres","Mi-profondeur"))) %>%
  arrange(desc(Ordre)) %>%
  filter(duplicated(ID.interne.passage) | n()==1)

datawithmetrics_unique <- subset(datawithmetrics, !(ID.interne.passage %in% unique(doublons_final$ID.interne.passage)))
datawithmetrics <- bind_rows(datawithmetrics_unique,resultat_filtre_final)

# Still some duplicates
doublons_final <- datawithmetrics[duplicated(datawithmetrics$ID.interne.passage) |
                                    duplicated(datawithmetrics$ID.interne.passage, fromLast = TRUE), ]

resultat_filtre_final <- doublons_final %>%
  filter(Prelevement.niveau %in% c("Surface (0-1m)", "2 mètres", "de 3 à 5 mètres","Mi-profondeur")) %>%
  group_by(ID.interne.passage) %>%
  mutate(Ordre = match(Prelevement.niveau, c("Surface (0-1m)", "2 mètres", "de 3 à 5 mètres","Mi-profondeur"))) %>%
  arrange(desc(Ordre)) %>%
  filter(duplicated(ID.interne.passage) | n()==1)

datawithmetrics_unique <- subset(datawithmetrics, !(ID.interne.passage %in% unique(doublons_final$ID.interne.passage)))
datawithmetrics <- bind_rows(datawithmetrics_unique,resultat_filtre_final)
datawithmetrics <- select(datawithmetrics,-Ordre)

# that's ok, save it
#write.csv2(datawithmetrics,file="output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetricsfinal.csv", row.names = FALSE,dec = ".")

### Add the observations coordinates on the PCA to the data
# Import data
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)
# Import PCA results
PCA_coord <- read_delim("output/tableaux/Networks/PCA_coords.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)
# Bind
data <- cbind(data,PCA_coord)
# Save
#write.csv2(data,file="output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_PCA_coord.final.csv", row.names = FALSE,dec = ".")

#### Define the different moments : before, during and after the bloom ####
# Import data
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_PCA_coord.final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Make it chronological
data <- data %>%
  arrange(Code_point_Libelle, Date)

data$Moment <- "Nothing"
data[!is.na(data$Bloom),]$Moment <- "During" 
# Determine the "before" moment
for (i in 1:nrow(data)){
  if (is.na(data[i,]$Bloom) & !is.na(data[i+1,]$Bloom) & (data[i,]$Code_point_Libelle == data[i+1,]$Code_point_Libelle )){
    data[i,"Moment"]$Moment <- "Before"
  }
  if (!is.na(data[i,]$Bloom) & !is.na(data[i+1,]$Bloom) & (data[i,]$Code_point_Libelle == data[i+1,]$Code_point_Libelle )){
    data[i,"Moment"]$Moment <- "During"
  }
  if (!is.na(data[i,]$Bloom) & !is.na(data[i+1,]$Bloom) & (data[i,]$Code_point_Libelle == data[i+1,]$Code_point_Libelle )){
    data[i,"Moment"]$Moment <- "During"
  }
}

# Determine the after moment
for (i in 2:(nrow(data)-1)){
  if (data[i,]$Moment == "Nothing" & data[i-1,]$Moment == "During"  & (data[i,]$Code_point_Libelle == data[i-1,]$Code_point_Libelle )){
    data[i,"Moment"]$Moment <- "After"
  }
}

# Determine the during if we are between to date of bloom so we can't define an after and before for these blooms
for (i in 2:(nrow(data)-1)){
  if (data[i-1,]$Moment == "During" & data[i+1,]$Moment == "During"  & (data[i,]$Code_point_Libelle == data[i-1,]$Code_point_Libelle )){
    data[i,"Moment"]$Moment <- "During"
  }
}

# Indicate the taxa responsible of the bloom for the before and after moments
data <- data %>%
  mutate(
    Bloom = if_else(Moment == "Before", 
                    lead(Bloom),  
                    Bloom),
    Bloom_Phylum = if_else(Moment == "Before", 
                    lead(Bloom_Phylum),  
                    Bloom_Phylum),
    Bloom_Genre = if_else(Moment == "Before", 
                           lead(Bloom_Genre),  
                         Bloom_Genre)
    
  )

data <- data %>%
  mutate(
    Bloom = if_else(Moment == "After", 
                    lag(Bloom),  
                    Bloom) ,
    Bloom_Phylum = if_else(Moment == "After", 
                    lag(Bloom_Phylum),  
                    Bloom_Phylum) ,
    Bloom_Genre = if_else(Moment == "After", 
                    lag(Bloom_Genre),  
                    Bloom_Genre) 
  )

# Save for the before -1 
data <- data %>%
  mutate(
    Bloom = if_else(Moment == "Before -1", 
                    lead(Bloom),  
                    Bloom),
    Bloom_Phylum = if_else(Moment == "Before -1", 
                           lead(Bloom_Phylum),  
                           Bloom_Phylum),
    Bloom_Genre = if_else(Moment == "Before -1", 
                          lead(Bloom_Genre),  
                          Bloom_Genre)
    
  )

# Save it
#write.csv2(data,file="output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_PCA_coord_moments.final.csv", row.names = FALSE,dec = ".")

