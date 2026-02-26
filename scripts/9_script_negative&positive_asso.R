################################################################################
# Diatoms vs dinoflagellates: a network analysis of bloom impacts on diversity #
#                    and phytoplankton associations | R scripts                #
################################################################################

# Script used to work on negative and positive associations #
# 02/26/2026

# Load packages
library(readr)
library(NetCoMi)
library(dplyr)
library(igraph)
library(highcharter)
library(pulsar)
library(lubridate)
library(tidyr)
library(DescTools)
library(FactoMineR)
library(factoextra)
library(missMDA)
library(vegan)
library(corrplot)
library(cowplot)

# Import data
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_blooms_caracterised_final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Compute the global network #####
### Cluster 1-Mediterranean sea ####
# Select region 1
CL1 <- filter(data, region == "1-Mediterranean sea")

# Select all genus
CL1 <- dplyr::select(CL1,Actinoptychus:Coscinodiscophycidae)

# Replace NA's by 0 to make it uniform
CL1[is.na(CL1)]<- 0

CL1 <- as.matrix(CL1)

# Build the global network
net_cluster1 <- netConstruct(data = CL1, dataType = "counts",measure = "spearman", 
                             filtTax = "numbSamp",filtTaxPar = list(numbSamp = 30),
                             filtSamp = "none",sparsMethod = "t-test",alpha = 0.05, zeroMethod = "none",adjust = "adaptBH",
                             normMethod = "none", dissFunc = "signed")

# Build info about the graph
net_props_cluster1 <- netAnalyze(net_cluster1,
                                 clustMethod = "cluster_fast_greedy",
                                 hubPar = c("eigenvector"),
                                 graphlet = F,
                                 connectivity = T)
# Store the associations list
# Keep only the positive associations
edgelist_1_all <- net_cluster1$edgelist1[,c("v1","v2","asso")]
# Save the edge list
write.csv2(edgelist_1_all,file="output/tableaux/Networks/edgelist_med_all.csv", row.names = FALSE,dec = ".")

# Store the association matrix
assoMat_1 <- as.data.frame(net_cluster1$assoMat1)
# Keep all the associations
assoMat_1 <- assoMat_1[, colSums(assoMat_1) != 0]
assoMat_1_all <- assoMat_1[rowSums(assoMat_1) != 0, ]
write.csv2(assoMat_1_all,file="output/tableaux/Networks/assoMat_med_all.csv", row.names = FALSE,dec = ".")
# Save the network composition
taxon_list_1_all <- as.data.frame(rownames(assoMat_1_all))
colnames(taxon_list_1_all) <- "Taxon"
write.csv2(taxon_list_1_all,file="output/tableaux/Networks/taxon_list_med_all.csv", row.names = FALSE,dec = ".")

### Working on the Mediterranean sea #####
# Import data
data_edge <- read_delim("output/tableaux/Networks/edgelist_med_all.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Load table "Liste_phylum.classe" to associate species to a class
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt",header = T,sep = "\t")

# Create a new dataframe that indicate the association type
# Working on the edges
edge1 <- as.data.frame(data_edge$v1)
colnames(edge1) <- "Taxon"
edge1 <- left_join(edge1,phylum_classe)
colnames(edge1) <- c("v1","v1_classe")

edge2 <- as.data.frame(data_edge$v2)
colnames(edge2) <- "Taxon"
edge2 <- left_join(edge2,phylum_classe)
colnames(edge2) <- c("v2","v2_classe")

edgelist_genus_med <- cbind(edge1,edge2,data_edge$asso)
# Associate the class
edgelist_genus_med$link_genus <- paste0(edgelist_genus_med$v1_classe,"-",edgelist_genus_med$v2_classe) 
# Function to normalize the association type, that's make Dinophyceae-Bacillariophyceae => Bacillariophyceae-Dinophyceae
normaliser_paire <- function(x) {
  parties <- sort(unlist(strsplit(x, "-")))
  paste(parties, collapse = "-")
}
# Applying it 
edgelist_genus_med$link_genus <- sapply(edgelist_genus_med$link_genus, normaliser_paire)

edgelist_genus_med$Sign <- ifelse(edgelist_genus_med$`data_edge$asso` > 0, "+","-")
# Save it
write.csv2(edgelist_genus_med,file="output/tableaux/Networks/edgelist_genus_med_all.csv", row.names = FALSE,dec = ".")

### Cluster 2-Eastern Channel - North Sea ####
# Select region 2
CL2 <- filter(data, region == "2-Eastern Channel - North Sea")

# Select all genus
CL2 <- dplyr::select(CL2,Actinoptychus:Coscinodiscophycidae)

# Replace NA's by 0 to make it uniform
CL2[is.na(CL2)]<- 0

CL2 <- as.matrix(CL2)

# Build the global network
net_cluster2 <- netConstruct(data = CL2, dataType = "counts",measure = "spearman", 
                             filtTax = "numbSamp",filtTaxPar = list(numbSamp = 30),
                             filtSamp = "none",sparsMethod = "t-test",alpha = 0.05, zeroMethod = "none",adjust = "adaptBH",
                             normMethod = "none", dissFunc = "signed")

# Build info about the graph
net_props_cluster2 <- netAnalyze(net_cluster2,
                                 clustMethod = "cluster_fast_greedy",
                                 hubPar = c("eigenvector"),
                                 graphlet = F,
                                 connectivity = T)

# Store the associations list
# Keep only the positive associations
edgelist_2_all <- net_cluster2$edgelist1[,c("v1","v2","asso")]
# Save the edge list
write.csv2(edgelist_2_all,file="output/tableaux/Networks/edgelist_manche_all.csv", row.names = FALSE,dec = ".")

# Store the association matrix
assoMat_2 <- as.data.frame(net_cluster2$assoMat1)
# Keep all associations
assoMat_2 <- assoMat_2[, colSums(assoMat_2) != 0]
assoMat_2_all <- assoMat_2[rowSums(assoMat_2) != 0, ]
write.csv2(assoMat_2_all,file="output/tableaux/Networks/assoMat_manche_all.csv", row.names = FALSE,dec = ".")
# Save the network composition
taxon_list_2_all <- as.data.frame(rownames(assoMat_2))
colnames(taxon_list_2_all) <- "Taxon"
write.csv2(taxon_list_2_all,file="output/tableaux/Networks/taxon_list_manche_all.csv", row.names = FALSE,dec = ".")

# Import data
data_edge <- read_delim("output/tableaux/Networks/edgelist_manche_all.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)
# Load table "Liste_phylum.classe" to associate species to a class
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt",header = T,sep = "\t")

# Create a new dataframe that indicate the association type
# Working on the edges
edge1 <- as.data.frame(data_edge$v1)
colnames(edge1) <- "Taxon"
edge1 <- left_join(edge1,phylum_classe)
colnames(edge1) <- c("v1","v1_classe")

edge2 <- as.data.frame(data_edge$v2)
colnames(edge2) <- "Taxon"
edge2 <- left_join(edge2,phylum_classe)
colnames(edge2) <- c("v2","v2_classe")

edgelist_genus_manche <- cbind(edge1,edge2,data_edge$asso)
# Associate the class
edgelist_genus_manche$link_genus <- paste0(edgelist_genus_manche$v1_classe,"-",edgelist_genus_manche$v2_classe) 

# Function to normalize the association type, that's make Dinophyceae-Bacillariophyceae => Bacillariophyceae-Dinophyceae
normaliser_paire <- function(x) {
  parties <- sort(unlist(strsplit(x, "-")))
  paste(parties, collapse = "-")
}
# Applying it 
edgelist_genus_manche$link_genus <- sapply(edgelist_genus_manche$link_genus, normaliser_paire)

edgelist_genus_manche$Sign <- ifelse(edgelist_genus_manche$`data_edge$asso` > 0, "+","-")
# Save it
write.csv2(edgelist_genus_manche,file="output/tableaux/Networks/edgelist_genus_manche_all.csv", row.names = FALSE,dec = ".")


### Cluster 3-Atlantic - Western Channel ####
# Select cluster 
CL3 <- filter(data, region == "3-Atlantic - Western Channel" )

# Select all genus
CL3 <- dplyr::select(CL3,Actinoptychus:Coscinodiscophycidae)

# Replace NA's by 0 to make it uniform
CL3[is.na(CL3)]<- 0

CL3 <- as.matrix(CL3)

# Build the global network
net_cluster3 <- netConstruct(data = CL3, dataType = "counts",measure = "spearman", 
                             filtTax = "numbSamp",filtTaxPar = list(numbSamp = 30),
                             filtSamp = "none",sparsMethod = "t-test",alpha = 0.05, zeroMethod = "none",adjust = "adaptBH",
                             normMethod = "none", dissFunc = "signed")

# Build info about the graph
net_props_cluster3 <- netAnalyze(net_cluster3,
                                 clustMethod = "cluster_fast_greedy",
                                 hubPar = c("eigenvector"),
                                 graphlet = F,
                                 connectivity = T)

# Store the associations list
# Keep all associations
edgelist_3_all <- net_cluster3$edgelist1[,c("v1","v2","asso")]
# Save the edge list
write.csv2(edgelist_3_all,file="output/tableaux/Networks/edgelist_atlantic_all.csv", row.names = FALSE,dec = ".")

# Store the association matrix
assoMat_3 <- as.data.frame(net_cluster3$assoMat1)
# Keep all associations
assoMat_3 <- assoMat_3[, colSums(assoMat_3) != 0]
assoMat_3_all <- assoMat_3[rowSums(assoMat_3) != 0, ]
write.csv2(assoMat_3_all,file="output/tableaux/Networks/assoMat_atlantic_all.csv", row.names = FALSE,dec = ".")
# Save the network composition
taxon_list_3_all <- as.data.frame(rownames(assoMat_3))
colnames(taxon_list_3_all) <- "Taxon"
  write.csv2(taxon_list_3_all,file="output/tableaux/Networks/taxon_list_atlantic_all.csv", row.names = FALSE,dec = ".")

# Import data
data_edge <- read_delim("output/tableaux/Networks/edgelist_Atlantic_all.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)
# Load table "Liste_phylum.classe" to associate species to a class
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt",header = T,sep = "\t")

# Create a new dataframe that indicate the association type
# Working on the edges
edge1 <- as.data.frame(data_edge$v1)
colnames(edge1) <- "Taxon"
edge1 <- left_join(edge1,phylum_classe)
colnames(edge1) <- c("v1","v1_classe")

edge2 <- as.data.frame(data_edge$v2)
colnames(edge2) <- "Taxon"
edge2 <- left_join(edge2,phylum_classe)
colnames(edge2) <- c("v2","v2_classe")

edgelist_genus_atlantic <- cbind(edge1,edge2,data_edge$asso)
# Associate the class
edgelist_genus_atlantic$link_genus <- paste0(edgelist_genus_atlantic$v1_classe,"-",edgelist_genus_atlantic$v2_classe) 

# Function to normalize the association type, that's make Dinophyceae-Bacillariophyceae => Bacillariophyceae-Dinophyceae
normaliser_paire <- function(x) {
  parties <- sort(unlist(strsplit(x, "-")))
  paste(parties, collapse = "-")
}
# Applying it
edgelist_genus_atlantic$link_genus <- sapply(edgelist_genus_atlantic$link_genus, normaliser_paire)
edgelist_genus_atlantic$Sign <- ifelse(edgelist_genus_atlantic$`data_edge$asso` > 0, "+","-")
# Save it
write.csv2(edgelist_genus_atlantic,file="output/tableaux/Networks/edgelist_genus_atlantic_all.csv", row.names = FALSE,dec = ".")

# Compute the temporal networks #####
### Cluster 1-Mediterranean sea ####
# Import data 
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_blooms_caracterised_final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Import dataset to associate species to a class
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt", header = TRUE, sep = "\t")

# Import the association matrix
assoMat_all <- read_delim("output/tableaux/Networks/assoMat_med_all.csv", 
                      delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                          grouping_mark = ""), trim_ws = TRUE)
# Import the edge list
edge_all <- read_delim("output/tableaux/Networks/edgelist_genus_med_all.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Create the graph
cluster1 <- graph_from_adjacency_matrix(as.matrix(assoMat_all), weighted = TRUE, mode = "undirected", diag = FALSE)


# Preparation to temporal network construction 
# Index of species
phyto_index <- as.data.frame(V(cluster1))
phyto_index$Taxon <- rownames(phyto_index)
colnames(phyto_index)[1] <- "Pindex"

# Table preparation to retrieve station/date information
CL1df <- filter(data,region == "1-Mediterranean sea")

# Select all genus
CL1 <- dplyr::select(CL1df,Actinoptychus:Coscinodiscophycidae)

# Replace NA's by 0 to make it uniform
CL1[is.na(CL1)]<- 0

CL1 <- as.matrix(CL1)

# Creating a df to store the results
data_results_reseaux <- c("","")
data_results_reseaux <- as.data.frame(data_results_reseaux)



# Compute the subgraphs
for (i in 1:nrow(CL1)){ 
  spe <- as.data.frame(CL1[i,]) # Pick the abundances of the day
  colnames(spe) <- "Count" 
  spe$Taxon <- rownames(spe) # Pick the phytoplankton name
  spe <- left_join(phyto_index,spe, by = join_by(Taxon))
  spe <- filter(spe,Count>0) # Pick the phytoplankton present at this day
  spe$Pindex <- as.numeric(spe$Pindex) # Index associate to the phyto taxa
  spe <- left_join(spe,phylum_classe)
  station <- CL1df[i,]$Code_point_Libelle # Station at this line
  date <- CL1df[i,]$Date # Date at this line
  
  if (nrow(spe) != 0){ # If there are at least one taxa we can make a graph
    vids <- spe$Pindex # Indicate which phytoplankton are present at this day to keep them
    sub <- igraph::subgraph(cluster1, vids) # Makes the subgraph
    
    sub_edges <- filter(edge_all,v1 %in% spe$Taxon & v2 %in% spe$Taxon) #Indicate the edges in the subgraph
    
    
    data_results_reseaux[i,1] <- station
    data_results_reseaux[i,2] <- date

    data_results_reseaux[i,3] <- nrow(filter(spe,Phylum.Classe == "Bacillariophyceae")) # number of Bacillariophyceae
    data_results_reseaux[i,4] <- nrow(filter(spe,Phylum.Classe == "Bacillariophyceae"))/nrow(spe) # Pourcentage of Bacillariophyceae
    data_results_reseaux[i,5] <- nrow(filter(spe,Phylum.Classe == "Dinophyceae")) # number of Dinophyceae
    data_results_reseaux[i,6] <- nrow(filter(spe,Phylum.Classe == "Dinophyceae"))/nrow(spe) # Pourcentage of Dinophyceae
    data_results_reseaux[i,7] <- nrow(filter(spe,Phylum.Classe != "Bacillariophyceae" & Phylum.Classe != "Dinophyceae")) # number of other taxa
    data_results_reseaux[i,8] <- nrow(filter(spe,Phylum.Classe != "Bacillariophyceae"& Phylum.Classe != "Dinophyceae"))/nrow(spe) # Pourcentage of other taxa
    
    data_results_reseaux[i,9] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Bacillariophyceae-Bacillariophyceae")) # number of Bac-Bac association -
    data_results_reseaux[i,10] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Bacillariophyceae-Bacillariophyceae"))/nrow(sub_edges) # Pourcentage of Bac-Bac association - out of total +/-
    data_results_reseaux[i,11] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Bacillariophyceae-Bacillariophyceae"))/nrow(filter(sub_edges,Sign == "-")) # Pourcentage of Bac-Bac association - out of total -
    data_results_reseaux[i,12] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Bacillariophyceae-Dinophyceae")) # number of Bac-Dino association -
    data_results_reseaux[i,13] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Bacillariophyceae-Dinophyceae"))/nrow(sub_edges) # Pourcentage of Bac-Dino association - out of total +/-
    data_results_reseaux[i,14] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Bacillariophyceae-Dinophyceae"))/nrow(filter(sub_edges,Sign == "-")) # Pourcentage of Bac-Dino association - out of total -
    data_results_reseaux[i,15] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Dinophyceae-Dinophyceae")) # number of Dino-Dino association -
    data_results_reseaux[i,16] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Dinophyceae-Dinophyceae"))/nrow(sub_edges) # Pourcentage of Dino-Dino association - out of total +/-
    data_results_reseaux[i,17] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Dinophyceae-Dinophyceae"))/nrow(filter(sub_edges,Sign == "-")) # Pourcentage of Dino-Dino association - out of total -
    data_results_reseaux[i,18] <- nrow(filter(sub_edges,Sign == "-" & link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" )) # number of other associations -
    data_results_reseaux[i,19] <- nrow(filter(sub_edges,Sign == "-" & link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" ))/nrow(sub_edges) # Pourcentage of other associations - out of total +/-
    data_results_reseaux[i,20] <- nrow(filter(sub_edges,Sign == "-" & link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" ))/nrow(filter(sub_edges,Sign == "-")) # Pourcentage of other associations - out of total -

    
  } 
  else { 
    data_results_reseaux[i,1] <- station
    data_results_reseaux[i,2] <- date
    data_results_reseaux[i,3] <- NA #
    data_results_reseaux[i,4] <- NA # 
    data_results_reseaux[i,5] <- NA # 
    data_results_reseaux[i,6] <- NA # 
    
    data_results_reseaux[i,7] <- NA # 
    data_results_reseaux[i,8] <- NA # 
    data_results_reseaux[i,9] <- NA # 
    
    data_results_reseaux[i,10] <- NA # 
    data_results_reseaux[i,11] <- NA # 
    data_results_reseaux[i,12] <- NA # 
    data_results_reseaux[i,13] <- NA # 
    
    data_results_reseaux[i,14] <- NA # 
    data_results_reseaux[i,15] <- NA # 
    
    data_results_reseaux[i,16] <- NA # 
    
    data_results_reseaux[i,17] <- NA # 
    data_results_reseaux[i,18] <- NA # 
    data_results_reseaux[i,19] <- NA # 
    data_results_reseaux[i,20] <- NA # 
    
  }
  
  
  colnames(data_results_reseaux) <- c("Code_point_Libelle","Date","N_Bac","P_Bac","N_Dino","P_Dino","N_others","P_others","N_BacBac_neg","P_BacBac_neg_tot","P_BacBac_neg_neg",
  "N_BacDino_neg","P_BacDino_neg_tot","P_BacDino_neg_neg","N_DinoDino_neg","P_DinoDino_neg_tot","P_DinoDino_neg_neg","N_OtherA_neg","P_OtherA_neg_tot","P_OtherA_neg_neg"
  )
  print(i/nrow(CL1)*100)
}
# Save it
write.csv2(data_results_reseaux,file="output/tableaux/Networks/subnetworks/results_metrics_reseaux_cluster1_neg.csv", row.names = FALSE,dec = ".")

### Cluster 2-Eastern Channel - North Sea ####
# Import data
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_blooms_caracterised_final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)
# Import dataset to associate species to a class
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt", header = TRUE, sep = "\t")

# Import the association matrix
assoMat_all <- read_delim("output/tableaux/Networks/assoMat_manche_all.csv", 
                      delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                          grouping_mark = ""), trim_ws = TRUE)
# Import the edge list
edge_all <- read_delim("output/tableaux/Networks/edgelist_genus_manche_all.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Create the graph
cluster2 <- graph_from_adjacency_matrix(as.matrix(assoMat_all), weighted = TRUE, mode = "undirected", diag = FALSE)


# Preparation to temporal network construction
# Index of species
phyto_index <- as.data.frame(V(cluster2))
phyto_index$Taxon <- rownames(phyto_index)
colnames(phyto_index)[1] <- "Pindex"

# Table preparation to retrieve station/date information
CL2df <- filter(data,region == "2-Eastern Channel - North Sea")

# Select all genus
CL2 <- dplyr::select(CL2df,Actinoptychus:Coscinodiscophycidae)

# Replace NA's by 0 to make it uniform
CL2[is.na(CL2)]<- 0

CL2 <- as.matrix(CL2)

# Creating a df to store the results
data_results_reseaux <- c("","")
data_results_reseaux <- as.data.frame(data_results_reseaux)



# Compute the subgraphs
for (i in 1:nrow(CL2)){ 
  spe <- as.data.frame(CL2[i,]) # Pick the abundances of the day
  colnames(spe) <- "Count" 
  spe$Taxon <- rownames(spe) # Pick the phytoplankton name
  spe <- left_join(phyto_index,spe, by = join_by(Taxon))
  spe <- filter(spe,Count>0) # Pick the phytoplankton present at this day
  spe$Pindex <- as.numeric(spe$Pindex) # Index associate to the phyto taxa
  spe <- left_join(spe,phylum_classe)
  station <- CL2df[i,]$Code_point_Libelle # Station at this line
  date <- CL2df[i,]$Date  # Date at this line
  
    if (nrow(spe) != 0){ # If there are at least one taxa we can make a graph
      vids <- spe$Pindex # Indicate which phytoplankton are present at this day to keep them
      sub <- igraph::subgraph(cluster2, vids) # Makes the subgraph
      
      sub_edges <- filter(edge_all,v1 %in% spe$Taxon & v2 %in% spe$Taxon) #Indicate the edges in the subgraph
      
      
      data_results_reseaux[i,1] <- station
      data_results_reseaux[i,2] <- date
      
      data_results_reseaux[i,3] <- nrow(filter(spe,Phylum.Classe == "Bacillariophyceae")) # number of Bacillariophyceae
      data_results_reseaux[i,4] <- nrow(filter(spe,Phylum.Classe == "Bacillariophyceae"))/nrow(spe) # Pourcentage of Bacillariophyceae
      data_results_reseaux[i,5] <- nrow(filter(spe,Phylum.Classe == "Dinophyceae")) # number of Dinophyceae
      data_results_reseaux[i,6] <- nrow(filter(spe,Phylum.Classe == "Dinophyceae"))/nrow(spe) # Pourcentage of Dinophyceae
      data_results_reseaux[i,7] <- nrow(filter(spe,Phylum.Classe != "Bacillariophyceae" & Phylum.Classe != "Dinophyceae")) # number of other taxa
      data_results_reseaux[i,8] <- nrow(filter(spe,Phylum.Classe != "Bacillariophyceae"& Phylum.Classe != "Dinophyceae"))/nrow(spe) # Pourcentage of other taxa
      
      data_results_reseaux[i,9] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Bacillariophyceae-Bacillariophyceae")) # number of Bac-Bac association -
      data_results_reseaux[i,10] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Bacillariophyceae-Bacillariophyceae"))/nrow(sub_edges) # Pourcentage of Bac-Bac association - out of total +/-
      data_results_reseaux[i,11] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Bacillariophyceae-Bacillariophyceae"))/nrow(filter(sub_edges,Sign == "-")) # Pourcentage of Bac-Bac association - out of total -
      data_results_reseaux[i,12] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Bacillariophyceae-Dinophyceae")) # number of Bac-Dino association -
      data_results_reseaux[i,13] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Bacillariophyceae-Dinophyceae"))/nrow(sub_edges) # Pourcentage of Bac-Dino association - out of total +/-
      data_results_reseaux[i,14] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Bacillariophyceae-Dinophyceae"))/nrow(filter(sub_edges,Sign == "-")) # Pourcentage of Bac-Dino association - out of total -
      data_results_reseaux[i,15] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Dinophyceae-Dinophyceae")) # number of Dino-Dino association -
      data_results_reseaux[i,16] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Dinophyceae-Dinophyceae"))/nrow(sub_edges) # Pourcentage of Dino-Dino association - out of total +/-
      data_results_reseaux[i,17] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Dinophyceae-Dinophyceae"))/nrow(filter(sub_edges,Sign == "-")) # Pourcentage of Dino-Dino association - out of total -
      data_results_reseaux[i,18] <- nrow(filter(sub_edges,Sign == "-" & link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" )) # number of other associations -
      data_results_reseaux[i,19] <- nrow(filter(sub_edges,Sign == "-" & link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" ))/nrow(sub_edges) # Pourcentage of other associations - out of total +/-
      data_results_reseaux[i,20] <- nrow(filter(sub_edges,Sign == "-" & link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" ))/nrow(filter(sub_edges,Sign == "-")) # Pourcentage of other associations - out of total -
      
      
    } 
    else { 
      data_results_reseaux[i,1] <- station
      data_results_reseaux[i,2] <- date
      data_results_reseaux[i,3] <- NA #
      data_results_reseaux[i,4] <- NA # 
      data_results_reseaux[i,5] <- NA # 
      data_results_reseaux[i,6] <- NA # 
      
      data_results_reseaux[i,7] <- NA # 
      data_results_reseaux[i,8] <- NA # 
      data_results_reseaux[i,9] <- NA # 
      
      data_results_reseaux[i,10] <- NA # 
      data_results_reseaux[i,11] <- NA # 
      data_results_reseaux[i,12] <- NA # 
      data_results_reseaux[i,13] <- NA # 
      
      data_results_reseaux[i,14] <- NA # 
      data_results_reseaux[i,15] <- NA # 
      
      data_results_reseaux[i,16] <- NA # 
      
      data_results_reseaux[i,17] <- NA # 
      data_results_reseaux[i,18] <- NA # 
      data_results_reseaux[i,19] <- NA # 
      data_results_reseaux[i,20] <- NA # 
      
    }
    
    
    colnames(data_results_reseaux) <- c("Code_point_Libelle","Date","N_Bac","P_Bac","N_Dino","P_Dino","N_others","P_others","N_BacBac_neg","P_BacBac_neg_tot","P_BacBac_neg_neg",
                                        "N_BacDino_neg","P_BacDino_neg_tot","P_BacDino_neg_neg","N_DinoDino_neg","P_DinoDino_neg_tot","P_DinoDino_neg_neg","N_OtherA_neg","P_OtherA_neg_tot","P_OtherA_neg_neg"
    )
    print(i/nrow(CL2)*100)
  }
# Save it
write.csv2(data_results_reseaux,file="output/tableaux/Networks/subnetworks/results_metrics_reseaux_cluster2_neg.csv", row.names = FALSE,dec = ".")


### Cluster 3-Atlantic - Western Channel ####
# Import data 
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_blooms_caracterised_final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)
# Import dataset to associate species to a class
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt", header = TRUE, sep = "\t")

# Import the association matrix
assoMat_all <- read_delim("output/tableaux/Networks/assoMat_atlantic_all.csv", 
                      delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                          grouping_mark = ""), trim_ws = TRUE)
# Import the edge list
edge_all <- read_delim("output/tableaux/Networks/edgelist_genus_atlantic_all.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Create the graph
cluster3 <- graph_from_adjacency_matrix(as.matrix(assoMat_all), weighted = TRUE, mode = "undirected", diag = FALSE)


# Preparation to temporal network construction   
# Index of species
phyto_index <- as.data.frame(V(cluster3))
phyto_index$Taxon <- rownames(phyto_index)
colnames(phyto_index)[1] <- "Pindex"

# Table preparation to retrieve station/date information
CL3df <- filter(data,region == "3-Atlantic - Western Channel")

# Select all genus
CL3 <- dplyr::select(CL3df,Actinoptychus:Coscinodiscophycidae)

# Replace NA's by 0 to make it uniform
CL3[is.na(CL3)]<- 0

CL3 <- as.matrix(CL3)

# Creating a df to store the results
data_results_reseaux <- c("","")
data_results_reseaux <- as.data.frame(data_results_reseaux)



# Compute the subgraphs
for (i in 1:nrow(CL3)){ 
  spe <- as.data.frame(CL3[i,]) # Pick the abundances of the day
  colnames(spe) <- "Count" 
  spe$Taxon <- rownames(spe) # Pick the phytoplankton name
  spe <- left_join(phyto_index,spe, by = join_by(Taxon))
  spe <- filter(spe,Count>0) # Pick the phytoplankton present at this day
  spe$Pindex <- as.numeric(spe$Pindex) # Index associate to the phyto taxa
  spe <- left_join(spe,phylum_classe)
  station <- CL3df[i,]$Code_point_Libelle # Station at this line
  date <- CL3df[i,]$Date # Date at this line
  
  if (nrow(spe) != 0){ # If there are at least one taxa we can make a graph
    vids <- spe$Pindex # Indicate which phytoplankton are present at this day to keep them
    sub <- igraph::subgraph(cluster3, vids) # Makes the subgraph
    
    sub_edges <- filter(edge_all,v1 %in% spe$Taxon & v2 %in% spe$Taxon) #Indicate the edges in the subgraph
    
    
    data_results_reseaux[i,1] <- station
    data_results_reseaux[i,2] <- date
    
    data_results_reseaux[i,3] <- nrow(filter(spe,Phylum.Classe == "Bacillariophyceae")) # number of Bacillariophyceae
    data_results_reseaux[i,4] <- nrow(filter(spe,Phylum.Classe == "Bacillariophyceae"))/nrow(spe) # Pourcentage of Bacillariophyceae
    data_results_reseaux[i,5] <- nrow(filter(spe,Phylum.Classe == "Dinophyceae")) # number of Dinophyceae
    data_results_reseaux[i,6] <- nrow(filter(spe,Phylum.Classe == "Dinophyceae"))/nrow(spe) # Pourcentage of Dinophyceae
    data_results_reseaux[i,7] <- nrow(filter(spe,Phylum.Classe != "Bacillariophyceae" & Phylum.Classe != "Dinophyceae")) # number of other taxa
    data_results_reseaux[i,8] <- nrow(filter(spe,Phylum.Classe != "Bacillariophyceae"& Phylum.Classe != "Dinophyceae"))/nrow(spe) # Pourcentage of other taxa
    
    data_results_reseaux[i,9] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Bacillariophyceae-Bacillariophyceae")) # number of Bac-Bac association -
    data_results_reseaux[i,10] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Bacillariophyceae-Bacillariophyceae"))/nrow(sub_edges) # Pourcentage of Bac-Bac association - out of total +/-
    data_results_reseaux[i,11] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Bacillariophyceae-Bacillariophyceae"))/nrow(filter(sub_edges,Sign == "-")) # Pourcentage of Bac-Bac association - out of total -
    data_results_reseaux[i,12] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Bacillariophyceae-Dinophyceae")) # number of Bac-Dino association -
    data_results_reseaux[i,13] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Bacillariophyceae-Dinophyceae"))/nrow(sub_edges) # Pourcentage of Bac-Dino association - out of total +/-
    data_results_reseaux[i,14] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Bacillariophyceae-Dinophyceae"))/nrow(filter(sub_edges,Sign == "-")) # Pourcentage of Bac-Dino association - out of total -
    data_results_reseaux[i,15] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Dinophyceae-Dinophyceae")) # number of Dino-Dino association -
    data_results_reseaux[i,16] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Dinophyceae-Dinophyceae"))/nrow(sub_edges) # Pourcentage of Dino-Dino association - out of total +/-
    data_results_reseaux[i,17] <- nrow(filter(sub_edges,Sign == "-" & link_genus == "Dinophyceae-Dinophyceae"))/nrow(filter(sub_edges,Sign == "-")) # Pourcentage of Dino-Dino association - out of total -
    data_results_reseaux[i,18] <- nrow(filter(sub_edges,Sign == "-" & link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" )) # number of other associations -
    data_results_reseaux[i,19] <- nrow(filter(sub_edges,Sign == "-" & link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" ))/nrow(sub_edges) # Pourcentage of other associations - out of total +/-
    data_results_reseaux[i,20] <- nrow(filter(sub_edges,Sign == "-" & link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" ))/nrow(filter(sub_edges,Sign == "-")) # Pourcentage of other associations - out of total -
    
    
  } 
  else { 
    data_results_reseaux[i,1] <- station
    data_results_reseaux[i,2] <- date
    data_results_reseaux[i,3] <- NA #
    data_results_reseaux[i,4] <- NA # 
    data_results_reseaux[i,5] <- NA # 
    data_results_reseaux[i,6] <- NA # 
    
    data_results_reseaux[i,7] <- NA # 
    data_results_reseaux[i,8] <- NA # 
    data_results_reseaux[i,9] <- NA # 
    
    data_results_reseaux[i,10] <- NA # 
    data_results_reseaux[i,11] <- NA # 
    data_results_reseaux[i,12] <- NA # 
    data_results_reseaux[i,13] <- NA # 
    
    data_results_reseaux[i,14] <- NA # 
    data_results_reseaux[i,15] <- NA # 
    
    data_results_reseaux[i,16] <- NA # 
    
    data_results_reseaux[i,17] <- NA # 
    data_results_reseaux[i,18] <- NA # 
    data_results_reseaux[i,19] <- NA # 
    data_results_reseaux[i,20] <- NA # 
    
  }
  
  
  colnames(data_results_reseaux) <- c("Code_point_Libelle","Date","N_Bac","P_Bac","N_Dino","P_Dino","N_others","P_others","N_BacBac_neg","P_BacBac_neg_tot","P_BacBac_neg_neg",
                                      "N_BacDino_neg","P_BacDino_neg_tot","P_BacDino_neg_neg","N_DinoDino_neg","P_DinoDino_neg_tot","P_DinoDino_neg_neg","N_OtherA_neg","P_OtherA_neg_tot","P_OtherA_neg_neg"
  )
  print(i/nrow(CL3)*100)
}
# Save it
write.csv2(data_results_reseaux,file="output/tableaux/Networks/subnetworks/results_metrics_reseaux_cluster3_neg.csv", row.names = FALSE,dec = ".")

# Import data
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_PCA_coord_moments.final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

med_neg <- read_delim("output/tableaux/Networks/subnetworks/results_metrics_reseaux_cluster1_neg.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)
manche_neg <- read_delim("output/tableaux/Networks/subnetworks/results_metrics_reseaux_cluster2_neg.csv", 
                      delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                          grouping_mark = ""), trim_ws = TRUE)
atl_neg <- read_delim("output/tableaux/Networks/subnetworks/results_metrics_reseaux_cluster3_neg.csv", 
                      delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                          grouping_mark = ""), trim_ws = TRUE)
met_neg <- bind_rows(med_neg,manche_neg,atl_neg)

datawithmetrics <- left_join(data,met_neg)
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
write.csv2(datawithmetrics,file="output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_PCA_coord_moments.final_plusneg.csv", row.names = FALSE,dec = ".")

data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_PCA_coord_moments.final_plusneg.csv", 
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
blooms <- select(blooms, region,P_Bac,P_Dino,P_others,P_BacBac_neg_tot,P_BacDino_neg_tot,P_DinoDino_neg_tot,P_OtherA_neg_tot,Moment,Bloom_Phylum)

colnames(blooms) <- c("region","Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)","Moment","Bloom_Phylum")

blooms_longer <- pivot_longer(blooms, cols = c(Bac:`Other asso. (-)`))

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
    "Bac-Bac (-)" = "blue","Bac-Dino (-)"   = "orange"
    ,"Other taxa"      = "grey" ,"Dino-Dino (-)" = "green3","Other asso. (-)" = "grey"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL) +
  scale_colour_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac (-)" = "blue","Bac-Dino (-)"   = "orange"
    ,"Other taxa"      = "grey" ,"Dino-Dino (-)" = "green3","Other asso. (-)" = "grey"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)"),name = NULL)+
  facet_grid(Bloom_Phylum~region)+
  labs(x=NULL,y="")+
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0,0.8))

# Calculer les positions empilées des barres
blooms_summary_asso <- blooms_summary %>%
  filter(name %in% c("Bac-Bac (-)", "Bac-Dino (-)","Dino-Dino (-)", "Other asso. (-)")) %>%
  group_by(Moment, region,Bloom_Phylum) %>%
  mutate(cumsum_Val = cumsum(Val))  # Calcul de la somme cumulée

asso <- ggplot(filter(blooms_summary_asso, name %in% c("Bac-Bac (-)", "Bac-Dino (-)","Dino-Dino (-)", "Other asso. (-)") & (Bloom_Phylum=="Bacillariophyceae" | Bloom_Phylum=="Dinophyceae")), aes(Moment, Val)) +
  geom_col(aes(fill = name), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = Val, ymax = Val+sd, group = name,colour=name),
    width = 0.2, position = position_dodge(0.8)
  )+
  scale_fill_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac (-)" = "blue","Bac-Dino (-)"   = "orange"
    ,"Other taxa"      = "grey" ,"Dino-Dino (-)" = "green3","Other asso. (-)" = "grey"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL) +
  scale_colour_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac (-)" = "blue","Bac-Dino (-)"   = "orange"
    ,"Other taxa"      = "grey" ,"Dino-Dino (-)" = "green3","Other asso. (-)" = "grey"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)"),name = NULL)+
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

asso_bac <-  ggplot(filter(blooms_summary_asso, name %in% c("Bac-Bac (-)","Bac-Dino (-)", "Dino-Dino (-)", "Other asso. (-)") & Bloom_Phylum=="Bacillariophyceae"), aes(Moment, Val)) +
  geom_col(aes(fill = name), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = Val, ymax = Val+sd, group = name,colour=name),
    width = 0.2, position = position_dodge(0.8)
  )+
  scale_fill_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac (-)" = "blue","Bac-Dino (-)"   = "orange"
    ,"Other taxa"      = "grey" ,"Dino-Dino (-)" = "green3","Other asso. (-)" = "grey"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL) +
  scale_colour_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac (-)" = "blue","Bac-Dino (-)"   = "orange"
    ,"Other taxa"      = "grey" ,"Dino-Dino (-)" = "green3","Other asso. (-)" = "grey"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)"),name = NULL)+
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

asso_dino <-  ggplot(filter(blooms_summary_asso, name %in% c("Bac-Bac (-)", "Dino-Dino (-)", "Bac-Dino (-)", "Other asso. (-)") & Bloom_Phylum=="Dinophyceae"), aes(Moment, Val)) +
  geom_col(aes(fill = name), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = Val, ymax = Val+sd, group = name,colour=name),
    width = 0.2, position = position_dodge(0.8)
  )+
  scale_fill_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac (-)" = "blue","Bac-Dino (-)"   = "orange"
    ,"Other taxa"      = "grey" ,"Dino-Dino (-)" = "green3","Other asso. (-)" = "grey"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL) +
  scale_colour_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac (-)" = "blue","Bac-Dino (-)"   = "orange"
    ,"Other taxa"      = "grey" ,"Dino-Dino (-)" = "green3","Other asso. (-)" = "grey"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)"),name = NULL)+
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

asso_hapt <-  ggplot(filter(blooms_summary_asso, name %in% c("Bac-Bac (-)", "Dino-Dino (-)","Bac-Dino (-)", "Other asso. (-)") & Bloom_Phylum=="Haptophyta"), aes(Moment, Val)) +
  geom_col(aes(fill = name), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = Val, ymax = Val+sd, group = name,colour=name),
    width = 0.2, position = position_dodge(0.8)
  )+
  scale_fill_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac (-)" = "blue","Bac-Dino (-)"   = "orange"
    ,"Other taxa"      = "grey" ,"Dino-Dino (-)" = "green3","Other asso. (-)" = "grey"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL) +
  scale_colour_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac (-)" = "blue","Bac-Dino (-)"   = "orange"
    ,"Other taxa"      = "grey" ,"Dino-Dino (-)" = "green3","Other asso. (-)" = "grey"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)"),name = NULL)+

  facet_wrap(~region)+
  labs(x=NULL,y=NULL)+
  theme(legend.position = "none")+
  scale_y_continuous(limits = c(0,0.65))+
  theme(strip.text = element_blank(),axis.text.y = element_blank())
hapt <- plot_grid(compo_hapt,asso_hapt,nrow = 2)

plot_grid(bac,dino,hapt,nrow=1,rel_widths = c(1,1,0.3))
ggsave('asso_compo_dynam_blooms_review_neg.png', path = "output/graphs/bloom", dpi = 600, width = 360, height = 170, units = 'mm')


#### We take negative and positive associations together
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_PCA_coord_moments.final_plusneg.csv", 
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
blooms <- select(blooms, region,P_Bac,P_Dino,P_others,N_BacBac,N_BacBac_neg,N_BacDino,N_BacDino_neg,N_DinoDino,N_DinoDino_neg,N_AAutres,N_OtherA_neg,Moment,Bloom_Phylum)

colnames(blooms) <- c("region","Bac","Dino","Other taxa","Bac-Bac (+)","Bac-Bac (-)","Bac-Dino (+)","Bac-Dino (-)","Dino-Dino (+)","Dino-Dino (-)","Other asso. (+)","Other asso. (-)","Moment","Bloom_Phylum")

blooms[,5:12] <- blooms[,5:12]/rowSums(blooms[,5:12])


blooms_longer <- pivot_longer(blooms, cols = c(Bac:`Other asso. (-)`))

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
    "Bac-Bac (-)" = "blue","Bac-Dino (-)"   = "orange"
    ,"Other taxa"      = "grey" ,"Dino-Dino (-)" = "green3","Other asso. (-)" = "grey"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL) +
  scale_colour_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac (-)" = "blue","Bac-Dino (-)"   = "orange"
    ,"Other taxa"      = "grey" ,"Dino-Dino (-)" = "green3","Other asso. (-)" = "grey"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)"),name = NULL)+
  facet_grid(Bloom_Phylum~region)+
  labs(x=NULL,y="")+
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0,0.8))

#### We take negative and positive associations together but show only positive
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_PCA_coord_moments.final_plusneg.csv", 
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
blooms <- select(blooms, region,P_Bac,P_Dino,P_others,N_BacBac,N_BacBac_neg,N_BacDino,N_BacDino_neg,N_DinoDino,N_DinoDino_neg,N_AAutres,N_OtherA_neg,Moment,Bloom_Phylum)

colnames(blooms) <- c("region","Bac","Dino","Other taxa","Bac-Bac (+)","Bac-Bac (-)","Bac-Dino (+)","Bac-Dino (-)","Dino-Dino (+)","Dino-Dino (-)","Other asso. (+)","Other asso. (-)","Moment","Bloom_Phylum")

blooms[,5:12] <- blooms[,5:12]/rowSums(blooms[,5:12])


blooms_longer <- pivot_longer(blooms, cols = c(Bac:`Other asso. (-)`))

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
    "Bac-Bac (-)" = "blue","Bac-Dino (-)"   = "orange"
    ,"Other taxa"      = "grey" ,"Dino-Dino (-)" = "green3","Other asso. (-)" = "grey"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac","Bac-Dino","Dino-Dino","Other asso."),name = NULL) +
  scale_colour_manual(values = c(
    "Bac" = "#1f77b4", "Dino" = "green",
    "Bac-Bac (-)" = "blue","Bac-Dino (-)"   = "orange"
    ,"Other taxa"      = "grey" ,"Dino-Dino (-)" = "green3","Other asso. (-)" = "grey"
  ), breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)"),name = NULL)+
  facet_grid(Bloom_Phylum~region)+
  labs(x=NULL,y="")+
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0,0.8))

# Calculer les positions empilées des barres
blooms_summary_asso <- blooms_summary %>%
  filter(name %in% c("Bac-Bac (+)", "Bac-Dino (+)","Dino-Dino (+)", "Other asso. (+)")) %>%
  group_by(Moment, region,Bloom_Phylum) %>%
  mutate(cumsum_Val = cumsum(Val))  # Calcul de la somme cumulée

asso <- ggplot(filter(blooms_summary_asso, name %in% c("Bac-Bac (+)", "Bac-Dino (+)","Dino-Dino (+)", "Other asso. (+)") & (Bloom_Phylum=="Bacillariophyceae" | Bloom_Phylum=="Dinophyceae")), 
               aes(Moment, Val)) + geom_col(aes(fill = name), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar( aes(ymin = Val, ymax = Val+sd, group = name,colour=name), width = 0.2, position = position_dodge(0.8) )+ 
  scale_fill_manual(values = c( "Bac" = "#1f77b4", "Dino" = "green", "Bac-Bac (-)" = "blue","Bac-Dino (-)" = "orange" ,"Other taxa" = "grey" ,
                                "Dino-Dino (-)" = "green3","Other asso. (-)" = "grey", "Bac-Bac (+)" = "#56B4E9","Bac-Dino (+)" = "chocolate1" ,
                                "Other taxa" = "grey" ,"Dino-Dino (+)" = "#009E73","Other asso. (+)" = "grey30" ), 
                    breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)", 
                               "Bac-Bac (+)","Bac-Dino (+)","Dino-Dino (+)","Other asso.(+)"),name = NULL) + 
  scale_colour_manual(values = c( "Bac" = "#1f77b4", "Dino" = "green", "Bac-Bac (-)" = "blue","Bac-Dino (-)" = "orange" ,"Other taxa" = "grey" ,
                                  "Dino-Dino (-)" = "green3","Other asso. (-)" = "grey", "Bac-Bac (+)" = "#56B4E9","Bac-Dino (+)" = "chocolate1" ,"Other taxa" = "grey" ,
                                  "Dino-Dino (+)" = "#009E73","Other asso. (+)" = "grey30" ), breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)",
                                                                                                         "Bac-Bac (+)","Bac-Dino (+)","Dino-Dino (+)","Other asso.(+)"),name = NULL) + 
  facet_grid(Bloom_Phylum~region)+ 
  labs(x=NULL,y="")+ 
  theme(legend.position = "none")+ 
  scale_y_continuous(limits = c(0,0.6))+
  theme(strip.text = element_blank())

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

asso_bac <- ggplot(filter(blooms_summary_asso, name %in% c("Bac-Bac (+)", "Bac-Dino (+)","Dino-Dino (+)", "Other asso. (+)") & (Bloom_Phylum=="Bacillariophyceae")), 
                   aes(Moment, Val)) + geom_col(aes(fill = name), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar( aes(ymin = Val, ymax = Val+sd, group = name,colour=name), width = 0.2, position = position_dodge(0.8) )+ 
  scale_fill_manual(values = c( "Bac" = "#1f77b4", "Dino" = "green", "Bac-Bac (-)" = "blue","Bac-Dino (-)" = "orange" ,"Other taxa" = "grey" ,
                                "Dino-Dino (-)" = "green3","Other asso. (-)" = "grey", "Bac-Bac (+)" = "#56B4E9","Bac-Dino (+)" = "chocolate1" ,
                                "Other taxa" = "grey" ,"Dino-Dino (+)" = "#009E73","Other asso. (+)" = "grey30" ), 
                    breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)", 
                               "Bac-Bac (+)","Bac-Dino (+)","Dino-Dino (+)","Other asso.(+)"),name = NULL) + 
  scale_colour_manual(values = c( "Bac" = "#1f77b4", "Dino" = "green", "Bac-Bac (-)" = "blue","Bac-Dino (-)" = "orange" ,"Other taxa" = "grey" ,
                                  "Dino-Dino (-)" = "green3","Other asso. (-)" = "grey", "Bac-Bac (+)" = "#56B4E9","Bac-Dino (+)" = "chocolate1" ,"Other taxa" = "grey" ,
                                  "Dino-Dino (+)" = "#009E73","Other asso. (+)" = "grey30" ), breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)",
                                                                                                         "Bac-Bac (+)","Bac-Dino (+)","Dino-Dino (+)","Other asso.(+)"),name = NULL) + 
  facet_wrap(Bloom_Phylum~region)+ 
  labs(x=NULL,y="Proportion")+ 
  theme(legend.position = "none")+ 
  scale_y_continuous(limits = c(0,0.6))#+
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

asso_dino <- ggplot(filter(blooms_summary_asso, name %in% c("Bac-Bac (+)", "Bac-Dino (+)","Dino-Dino (+)", "Other asso. (+)") & (Bloom_Phylum=="Dinophyceae")), 
                    aes(Moment, Val)) + geom_col(aes(fill = name), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar( aes(ymin = Val, ymax = Val+sd, group = name,colour=name), width = 0.2, position = position_dodge(0.8) )+ 
  scale_fill_manual(values = c( "Bac" = "#1f77b4", "Dino" = "green", "Bac-Bac (-)" = "blue","Bac-Dino (-)" = "orange" ,"Other taxa" = "grey" ,
                                "Dino-Dino (-)" = "green3","Other asso. (-)" = "grey", "Bac-Bac (+)" = "#56B4E9","Bac-Dino (+)" = "chocolate1" ,
                                "Other taxa" = "grey" ,"Dino-Dino (+)" = "#009E73","Other asso. (+)" = "grey30" ), 
                    breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)", 
                               "Bac-Bac (+)","Bac-Dino (+)","Dino-Dino (+)","Other asso.(+)"),name = NULL) + 
  scale_colour_manual(values = c( "Bac" = "#1f77b4", "Dino" = "green", "Bac-Bac (-)" = "blue","Bac-Dino (-)" = "orange" ,"Other taxa" = "grey" ,
                                  "Dino-Dino (-)" = "green3","Other asso. (-)" = "grey", "Bac-Bac (+)" = "#56B4E9","Bac-Dino (+)" = "chocolate1" ,"Other taxa" = "grey" ,
                                  "Dino-Dino (+)" = "#009E73","Other asso. (+)" = "grey30" ), breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)",
                                                                                                         "Bac-Bac (+)","Bac-Dino (+)","Dino-Dino (+)","Other asso.(+)"),name = NULL) + 
  facet_wrap(Bloom_Phylum~region)+ 
  labs(x=NULL,y=NULL)+ 
  theme(legend.position = "none")+ 
  scale_y_continuous(limits = c(0,0.6))+
  theme(axis.text.y = element_blank())#+
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

asso_hapt <- ggplot(filter(blooms_summary_asso, name %in% c("Bac-Bac (+)", "Bac-Dino (+)","Dino-Dino (+)", "Other asso. (+)") & (Bloom_Phylum=="Haptophyta")), 
                    aes(Moment, Val)) + geom_col(aes(fill = name), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar( aes(ymin = Val, ymax = Val+sd, group = name,colour=name), width = 0.2, position = position_dodge(0.8) )+ 
  scale_fill_manual(values = c( "Bac" = "#1f77b4", "Dino" = "green", "Bac-Bac (-)" = "blue","Bac-Dino (-)" = "orange" ,"Other taxa" = "grey" ,
                                "Dino-Dino (-)" = "green3","Other asso. (-)" = "grey", "Bac-Bac (+)" = "#56B4E9","Bac-Dino (+)" = "chocolate1" ,
                                "Other taxa" = "grey" ,"Dino-Dino (+)" = "#009E73","Other asso. (+)" = "grey30" ), 
                    breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)", 
                               "Bac-Bac (+)","Bac-Dino (+)","Dino-Dino (+)","Other asso.(+)"),name = NULL) + 
  scale_colour_manual(values = c( "Bac" = "#1f77b4", "Dino" = "green", "Bac-Bac (-)" = "blue","Bac-Dino (-)" = "orange" ,"Other taxa" = "grey" ,
                                  "Dino-Dino (-)" = "green3","Other asso. (-)" = "grey", "Bac-Bac (+)" = "#56B4E9","Bac-Dino (+)" = "chocolate1" ,"Other taxa" = "grey" ,
                                  "Dino-Dino (+)" = "#009E73","Other asso. (+)" = "grey30" ), breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)",
                                                                                                         "Bac-Bac (+)","Bac-Dino (+)","Dino-Dino (+)","Other asso.(+)"),name = NULL) + 
  facet_wrap(Bloom_Phylum~region)+ 
  labs(x=NULL,y=NULL)+ 
  theme(legend.position = "none")+ 
  scale_y_continuous(limits = c(0,0.6))+
  theme(axis.text.y = element_blank())#+
  theme(strip.text = element_blank())

hapt <- plot_grid(compo_hapt,asso_hapt,nrow = 2)

# Calculer les positions empilées des barres
blooms_summary_asso <- blooms_summary %>%
  filter(name %in% c("Bac-Bac (-)", "Bac-Dino (-)","Dino-Dino (-)", "Other asso. (-)")) %>%
  group_by(Moment, region,Bloom_Phylum) %>%
  mutate(cumsum_Val = cumsum(Val))  # Calcul de la somme cumulée

asso_neg <- ggplot(filter(blooms_summary_asso, name %in% c("Bac-Bac (-)", "Bac-Dino (-)","Dino-Dino (-)", "Other asso. (-)") & (Bloom_Phylum=="Bacillariophyceae" | Bloom_Phylum=="Dinophyceae")), 
               aes(Moment, Val)) + geom_col(aes(fill = name), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar( aes(ymin = Val, ymax = Val+sd, group = name,colour=name), width = 0.2, position = position_dodge(0.8) )+ 
  scale_fill_manual(values = c( "Bac" = "#1f77b4", "Dino" = "green", "Bac-Bac (-)" = "blue","Bac-Dino (-)" = "orange" ,"Other taxa" = "grey" ,
                                "Dino-Dino (-)" = "green3","Other asso. (-)" = "grey", "Bac-Bac (+)" = "#56B4E9","Bac-Dino (+)" = "chocolate1" ,
                                "Other taxa" = "grey" ,"Dino-Dino (+)" = "#009E73","Other asso. (+)" = "grey30" ), 
                    breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)", 
                               "Bac-Bac (+)","Bac-Dino (+)","Dino-Dino (+)","Other asso.(+)"),name = NULL) + 
  scale_colour_manual(values = c( "Bac" = "#1f77b4", "Dino" = "green", "Bac-Bac (-)" = "blue","Bac-Dino (-)" = "orange" ,"Other taxa" = "grey" ,
                                  "Dino-Dino (-)" = "green3","Other asso. (-)" = "grey", "Bac-Bac (+)" = "#56B4E9","Bac-Dino (+)" = "chocolate1" ,"Other taxa" = "grey" ,
                                  "Dino-Dino (+)" = "#009E73","Other asso. (+)" = "grey30" ), breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)",
                                                                                                         "Bac-Bac (+)","Bac-Dino (+)","Dino-Dino (+)","Other asso.(+)"),name = NULL) + 
  facet_grid(Bloom_Phylum~region)+ 
  labs(x=NULL,y="")+ 
  theme(legend.position = "none")+ 
  scale_y_continuous(limits = c(0,0.6))+
  theme(strip.text = element_blank())

plot_grid(compo,asso_neg)



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

asso_bac_neg <- ggplot(filter(blooms_summary_asso, name %in% c("Bac-Bac (-)", "Bac-Dino (-)","Dino-Dino (-)", "Other asso. (-)") & (Bloom_Phylum=="Bacillariophyceae")), 
                   aes(Moment, Val)) + geom_col(aes(fill = name), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar( aes(ymin = Val, ymax = Val+sd, group = name,colour=name), width = 0.2, position = position_dodge(0.8) )+ 
  scale_fill_manual(values = c( "Bac" = "#1f77b4", "Dino" = "green", "Bac-Bac (-)" = "blue","Bac-Dino (-)" = "orange" ,"Other taxa" = "grey" ,
                                "Dino-Dino (-)" = "green3","Other asso. (-)" = "grey", "Bac-Bac (+)" = "#56B4E9","Bac-Dino (+)" = "chocolate1" ,
                                "Other taxa" = "grey" ,"Dino-Dino (+)" = "#009E73","Other asso. (+)" = "grey30" ), 
                    breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)", 
                               "Bac-Bac (+)","Bac-Dino (+)","Dino-Dino (+)","Other asso.(+)"),name = NULL) + 
  scale_colour_manual(values = c( "Bac" = "#1f77b4", "Dino" = "green", "Bac-Bac (-)" = "blue","Bac-Dino (-)" = "orange" ,"Other taxa" = "grey" ,
                                  "Dino-Dino (-)" = "green3","Other asso. (-)" = "grey", "Bac-Bac (+)" = "#56B4E9","Bac-Dino (+)" = "chocolate1" ,"Other taxa" = "grey" ,
                                  "Dino-Dino (+)" = "#009E73","Other asso. (+)" = "grey30" ), breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)",
                                                                                                         "Bac-Bac (+)","Bac-Dino (+)","Dino-Dino (+)","Other asso.(+)"),name = NULL) + 
  facet_grid(Bloom_Phylum~region)+ 
  labs(x=NULL,y="Proportion")+ 
  theme(legend.position = "none")+ 
  scale_y_continuous(limits = c(0,0.6))+
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

asso_dino_neg <- ggplot(filter(blooms_summary_asso, name %in% c("Bac-Bac (-)", "Bac-Dino (-)","Dino-Dino (-)", "Other asso. (-)") & (Bloom_Phylum=="Dinophyceae")), 
                    aes(Moment, Val)) + geom_col(aes(fill = name), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar( aes(ymin = Val, ymax = Val+sd, group = name,colour=name), width = 0.2, position = position_dodge(0.8) )+ 
  scale_fill_manual(values = c( "Bac" = "#1f77b4", "Dino" = "green", "Bac-Bac (-)" = "blue","Bac-Dino (-)" = "orange" ,"Other taxa" = "grey" ,
                                "Dino-Dino (-)" = "green3","Other asso. (-)" = "grey", "Bac-Bac (+)" = "#56B4E9","Bac-Dino (+)" = "chocolate1" ,
                                "Other taxa" = "grey" ,"Dino-Dino (+)" = "#009E73","Other asso. (+)" = "grey30" ), 
                    breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)", 
                               "Bac-Bac (+)","Bac-Dino (+)","Dino-Dino (+)","Other asso.(+)"),name = NULL) + 
  scale_colour_manual(values = c( "Bac" = "#1f77b4", "Dino" = "green", "Bac-Bac (-)" = "blue","Bac-Dino (-)" = "orange" ,"Other taxa" = "grey" ,
                                  "Dino-Dino (-)" = "green3","Other asso. (-)" = "grey", "Bac-Bac (+)" = "#56B4E9","Bac-Dino (+)" = "chocolate1" ,"Other taxa" = "grey" ,
                                  "Dino-Dino (+)" = "#009E73","Other asso. (+)" = "grey30" ), breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)",
                                                                                                         "Bac-Bac (+)","Bac-Dino (+)","Dino-Dino (+)","Other asso.(+)"),name = NULL) + 
  facet_grid(Bloom_Phylum~region)+ 
  labs(x=NULL,y=NULL)+ 
  theme(legend.position = "none")+ 
  scale_y_continuous(limits = c(0,0.6))+
  theme(axis.text.y = element_blank())+
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

asso_hapt_neg <- ggplot(filter(blooms_summary_asso, name %in% c("Bac-Bac -+)", "Bac-Dino (-)","Dino-Dino (-)", "Other asso. (-)") & (Bloom_Phylum=="Haptophyta")), 
                    aes(Moment, Val)) + geom_col(aes(fill = name), position = position_dodge(0.8), width = 0.7)+
  geom_errorbar( aes(ymin = Val, ymax = Val+sd, group = name,colour=name), width = 0.2, position = position_dodge(0.8) )+ 
  scale_fill_manual(values = c( "Bac" = "#1f77b4", "Dino" = "green", "Bac-Bac (-)" = "blue","Bac-Dino (-)" = "orange" ,"Other taxa" = "grey" ,
                                "Dino-Dino (-)" = "green3","Other asso. (-)" = "grey", "Bac-Bac (+)" = "#56B4E9","Bac-Dino (+)" = "chocolate1" ,
                                "Other taxa" = "grey" ,"Dino-Dino (+)" = "#009E73","Other asso. (+)" = "grey30" ), 
                    breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)", 
                               "Bac-Bac (+)","Bac-Dino (+)","Dino-Dino (+)","Other asso.(+)"),name = NULL) + 
  scale_colour_manual(values = c( "Bac" = "#1f77b4", "Dino" = "green", "Bac-Bac (-)" = "blue","Bac-Dino (-)" = "orange" ,"Other taxa" = "grey" ,
                                  "Dino-Dino (-)" = "green3","Other asso. (-)" = "grey", "Bac-Bac (+)" = "#56B4E9","Bac-Dino (+)" = "chocolate1" ,"Other taxa" = "grey" ,
                                  "Dino-Dino (+)" = "#009E73","Other asso. (+)" = "grey30" ), breaks = c("Bac","Dino","Other taxa","Bac-Bac (-)","Bac-Dino (-)","Dino-Dino (-)","Other asso. (-)",
                                                                                                         "Bac-Bac (+)","Bac-Dino (+)","Dino-Dino (+)","Other asso.(+)"),name = NULL) + 
  facet_grid(Bloom_Phylum~region)+ 
  labs(x=NULL,y=NULL)+ 
  theme(legend.position = "none")+ 
  scale_y_continuous(limits = c(0,0.6))+
  theme(axis.text.y = element_blank())+
  theme(strip.text = element_blank())

hapt <- plot_grid(compo_hapt,asso_hapt,nrow = 2)

compo_tot <- plot_grid(compo_bac,compo_dino,compo_hapt,nrow=1,rel_widths = c(1,1,0.3))
asso_pos <- plot_grid(asso_bac,asso_dino,asso_hapt,nrow=1,rel_widths = c(1,1,0.3))
asso_neg <- plot_grid(asso_bac_neg,asso_dino_neg,asso_hapt_neg,nrow=1,rel_widths = c(1,1,0.3))
plot_grid(asso_pos,asso_neg,nrow = 2)

plot_grid(bac,dino,hapt,nrow=1,rel_widths = c(1,1,0.3))
ggsave('asso_compo_dynam_blooms_review_pos_tot_3.png', path = "output/graphs/bloom", dpi = 600, width = 360, height = 170, units = 'mm')


### Compare composition of each global networks ####
# Import datasets
# Mediterranean sea
med <- read_delim("output/tableaux/Networks/taxon_list_med_all.csv", 
                  delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                      grouping_mark = ""), trim_ws = TRUE)
med$region <- "1-Mediterranean sea"

med$Number <- 1

# English channel
manche <- read_delim("output/tableaux/Networks/taxon_list_manche_all.csv", 
                     delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                         grouping_mark = ""), trim_ws = TRUE)
manche$region <- "2-Eastern Channel - North Sea"

manche$Number <- 1

# Atlantic
atlantic <- read_delim("output/tableaux/Networks/taxon_list_atlantic_all.csv", 
                       delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                           grouping_mark = ""), trim_ws = TRUE)
atlantic$region <- "3-Atlantic - Western Channel"

atlantic$Number <- 1

data <- rbind(med,atlantic,manche)

# Load table "Liste_phylum.classe" to associate species to a class
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt",header = T,sep = "\t")

# Associate class to taxa
data <- left_join(data,phylum_classe)

data$Phylum.Classe <- factor(data$Phylum.Classe, levels = c("Autres protistes","Chromista","Khakista","Chrysophyceae","Cyanophyceae","Dictyochophyceae",
                                                            "Raphidophyceae","Cryptophyceae", "Ebriophyceae","Haptophyta","Ciliophora",
                                                            "Chlorophyta","Euglenozoa","Dinophyceae","Bacillariophyceae"
))
# Composition for all the global networks
tableau_compo <- summarise(group_by(data,region,Phylum.Classe), Number = sum(Number)) |>
  mutate(Freq = Number / sum(Number))
# Save it
write.csv2(tableau_compo,file="output/tableaux/Networks/compo_global_networks_all.csv", row.names = FALSE,dec = ".")

# Frequence of the association analysis
freq_table_med <-  edgelist_genus_med %>%
  rowwise() %>%
  mutate(Combination = paste(sort(c(v1_classe, v2_classe)), collapse = "-")) %>%  
  ungroup() %>%
  group_by(Combination,Sign) %>%
  summarise(Number = n(), .groups = "drop") %>%
  separate(Combination, into = c("v1_classe", "v2_classe"), sep = "-")  

freq_table_med$Frequence <- freq_table_med$Number/sum(freq_table_med$Number)
freq_table_med$link_genus <- paste0(freq_table_med$v1_classe,"-",freq_table_med$v2_classe) 

#Ordering the class
freq_table_med <- freq_table_med[order(freq_table_med$Number,decreasing = T),]
freq_table_med$v1_classe <- factor(freq_table_med$v1_classe, levels = c("Bacillariophyceae", "Dinophyceae","Haptophyta", 
                                                                        "Chlorophyta", "Chrysophyceae", 
                                                                        "Cryptophyceae", "Cyanophyceae", 
                                                                        "Dictyochophyceae", "Ebriophyceae", 
                                                                        "Euglenozoa"))
freq_table_med$region <- "1-Mediterranean sea"
# Save it
write.csv2(freq_table_med,file="output/tableaux/Networks/freq_table_med_all.csv", row.names = FALSE,dec = ".")

# Frequence of the association analysis
freq_table_manche <-  edgelist_genus_manche %>%
  rowwise() %>%
  mutate(Combination = paste(sort(c(v1_classe, v2_classe)), collapse = "-")) %>%  
  ungroup() %>%
  group_by(Combination,Sign) %>%
  summarise(Number = n(), .groups = "drop") %>%
  separate(Combination, into = c("v1_classe", "v2_classe"), sep = "-")  

freq_table_manche$Frequence <- freq_table_manche$Number/sum(freq_table_manche$Number)
freq_table_manche$link_genus <- paste0(freq_table_manche$v1_classe,"-",freq_table_manche$v2_classe) 



#Ordering the class
freq_table_manche <- freq_table_manche[order(freq_table_manche$Number,decreasing = T),]
freq_table_manche$v1_classe <- factor(freq_table_manche$v1_classe, levels = c("Bacillariophyceae", "Dinophyceae","Haptophyta","Ciliophora", 
                                                                              "Chlorophyta", "Chrysophyceae", 
                                                                              "Cryptophyceae", "Cyanophyceae", 
                                                                              "Dictyochophyceae", "Ebriophyceae", 
                                                                              "Euglenozoa","Raphidophyceae"))

freq_table_manche$region <- "2-Eastern Channel - North Sea"

# Save it
write.csv2(freq_table_manche,file="output/tableaux/Networks/freq_table_manche_all.csv", row.names = FALSE,dec = ".")

# Frequence of the association analysis
freq_table_atlantic <-  edgelist_genus_atlantic %>%
  rowwise() %>%
  mutate(Combination = paste(sort(c(v1_classe, v2_classe)), collapse = "-")) %>%  
  ungroup() %>%
  group_by(Combination,Sign) %>%
  summarise(Number = n(), .groups = "drop") %>%
  separate(Combination, into = c("v1_classe", "v2_classe"), sep = "-")  

freq_table_atlantic$Frequence <- freq_table_atlantic$Number/sum(freq_table_atlantic$Number)
freq_table_atlantic$link_genus <- paste0(freq_table_atlantic$v1_classe,"-",freq_table_atlantic$v2_classe) 



#Ordering the class
freq_table_atlantic <- freq_table_atlantic[order(freq_table_atlantic$Number,decreasing = T),]
freq_table_atlantic$v1_classe <- factor(freq_table_atlantic$v1_classe, levels = c("Bacillariophyceae", "Dinophyceae","Haptophyta","Ciliophora", 
                                                                                  "Chlorophyta", "Chrysophyceae", 
                                                                                  "Cryptophyceae", "Cyanophyceae", 
                                                                                  "Dictyochophyceae", "Ebriophyceae", 
                                                                                  "Euglenozoa","Raphidophyceae","Autres protistes","Chromista",
                                                                                  "Khakista"))

freq_table_atlantic$region <- "3-Atlantic - Western Channel"
# Save it
write.csv2(freq_table_atlantic,file="output/tableaux/Networks/freq_table_atlantic_all.csv", row.names = FALSE,dec = ".")


# Same but for the associations types 
# Import the table that contains the association frequencies 
# For the mediterranean sea
med <- read_delim("output/tableaux/Networks/freq_table_med_all.csv", 
                  delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                      grouping_mark = ""), trim_ws = TRUE)
# English channel
manche <- read_delim("output/tableaux/Networks/freq_table_manche_all.csv", 
                     delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                         grouping_mark = ""), trim_ws = TRUE)
# Atlantic
atlantic <- read_delim("output/tableaux/Networks/freq_table_atlantic_all.csv", 
                       delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                           grouping_mark = ""), trim_ws = TRUE)
# Bind all
freq_table_big <- rbind(med,manche,atlantic)

# Name the "other" associations
freq_table_big$link_genus <- ifelse(freq_table_big$link_genus ==  "Bacillariophyceae-Bacillariophyceae", "Bacillariophyceae-Bacillariophyceae",
                                    ifelse(freq_table_big$link_genus ==  "Bacillariophyceae-Dinophyceae","Bacillariophyceae-Dinophyceae",
                                           ifelse(freq_table_big$link_genus ==  "Dinophyceae-Dinophyceae","Dinophyceae-Dinophyceae" ,"Other association"
                                           )))

freq_table_big$link_genus <- factor(freq_table_big$link_genus, levels = c("Other association","Dinophyceae-Dinophyceae","Bacillariophyceae-Dinophyceae","Bacillariophyceae-Bacillariophyceae"))

# Create a new dataset
freq_table_big2 <- freq_table_big |>
  group_by(region, link_genus, Sign) |>
  summarise(Number = sum(Number), .groups = "drop") |>
  group_by(region) |>
  mutate(Freq = Number / sum(Number)) |>
  ungroup()
write.csv2(freq_table_big2,file="output/tableaux/Networks/asso_global_networks_all.csv", row.names = FALSE,dec = ".")


# Making a new graph to synthetize association and composition
# Import the dataset about composition
tableau_compo <- read_delim("output/tableaux/Networks/compo_global_networks_all.csv", 
                            delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                                grouping_mark = ""), trim_ws = TRUE)
# Indicate the "other" phyla
tableau_compo$Phylum.Classe <- ifelse(tableau_compo$Phylum.Classe == "Bacillariophyceae" | tableau_compo$Phylum.Classe == "Dinophyceae", 
                                      tableau_compo$Phylum.Classe,"Other phyla")

tableau_compo <- summarise(group_by(tableau_compo,region,Phylum.Classe), Number = sum(Number)) |>
  mutate(Freq = Number / sum(Number))

# Create the dataframe as needed to make the graph
tableau_compo <- tableau_compo %>%
  group_by(region) %>%
  arrange(desc(Phylum.Classe)) %>%
  mutate(ymin = cumsum(Freq) - Freq,
         ymax = cumsum(Freq),
         labelPosition = (ymax + ymin) / 2,
         label = paste0(round(Freq * 100, 1), "%"))

# Import the frequencies data
freq_table_big <- read_delim("output/tableaux/Networks/asso_global_networks_all.csv", 
                             delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                                 grouping_mark = ""), trim_ws = TRUE)


# Create the dataframe as needed to make the graph
freq_table_big <- freq_table_big %>%
  mutate(
    sign_order = ifelse(Sign == "-", 1, 2),   # 1 = négatif, 2 = positif
    link_genus_sign = paste0(link_genus, " (", Sign, ")")
  )
freq_table_big <- freq_table_big %>%
  mutate(
    alpha_val = ifelse(Sign == "-", 0.4, 1)  # ajuste 0.3–0.6 selon rendu
  )

freq_table_big <- freq_table_big %>%
  group_by(region) %>%
  arrange(sign_order, link_genus_sign) %>%
  mutate(
    ymin = cumsum(Freq) - Freq,
    ymax = cumsum(Freq),
    labelPosition = (ymax + ymin) / 2,
    label = paste0(round(Freq * 100, 1), "%")
  )

# freq_table_big <- freq_table_big %>%
#   group_by(region) %>%
#   arrange(desc(link_genus)) %>%
#   mutate(ymin = cumsum(Freq) - Freq,
#          ymax = cumsum(Freq),
#          labelPosition = (ymax + ymin) / 2,
#          label = paste0(round(Freq * 100, 1), "%"))

freq_table_big$link_genus_sign <- paste0(freq_table_big$link_genus," (",freq_table_big$Sign,")")

# Making the graph
global_compo_asso <- ggplot() +
  geom_rect(aes(ymax = ymax, ymin = ymin, xmax = 3.9, xmin = 3, fill = Phylum.Classe),data = tableau_compo) +
  geom_text(aes(x = 4.5, y = labelPosition, label = label, colour = Phylum.Classe), size = 3,data = tableau_compo) +
  geom_rect(aes(ymax = ymax, ymin = ymin, xmax = 2.9, xmin = 2, fill = link_genus_sign),data = freq_table_big) +
  geom_text(aes(x = 1.4, y = labelPosition, label = label, colour = link_genus_sign), size = 3,data = freq_table_big) +
  scale_fill_manual(values = c(
    "Bacillariophyceae" = "#1f77b4", "Dinophyceae" = "green",
    "Bacillariophyceae-Bacillariophyceae (+)" = "#56B4E9","Bacillariophyceae-Bacillariophyceae (-)" = "blue",
    "Bacillariophyceae-Dinophyceae (+)"   = "chocolate1","Bacillariophyceae-Dinophyceae (-)"   = "orange"
    ,"Other phyla"      = "grey" ,"Dinophyceae-Dinophyceae (+)" = "#009E73","Dinophyceae-Dinophyceae (-)" = "green3",
    "Other association (+)" = "grey30","Other association (-)" = "grey90"
  ), breaks = c(
    "Bacillariophyceae","Dinophyceae","Other phyla", "Bacillariophyceae-Bacillariophyceae (+)", "Bacillariophyceae-Dinophyceae (+)", 
    "Dinophyceae-Dinophyceae (+)", "Other association (+)","Bacillariophyceae-Bacillariophyceae (-)", "Bacillariophyceae-Dinophyceae (-)", 
    "Dinophyceae-Dinophyceae (-)", "Other association (-)"
  ),name = NULL) +
  scale_colour_manual(values = c(
    "Bacillariophyceae" = "#1f77b4", "Dinophyceae" = "green",
    "Bacillariophyceae-Bacillariophyceae (+)" = "#56B4E9","Bacillariophyceae-Bacillariophyceae (-)" = "blue",
    "Bacillariophyceae-Dinophyceae (+)"   = "chocolate1","Bacillariophyceae-Dinophyceae (-)"   = "orange"
    ,"Other phyla"      = "grey" ,"Dinophyceae-Dinophyceae (+)" = "#009E73","Dinophyceae-Dinophyceae (-)" = "green3",
    "Other association (+)" = "grey30","Other association (-)" = "grey"
  ), breaks = c(
    "Bacillariophyceae","Dinophyceae","Other phyla", "Bacillariophyceae-Bacillariophyceae (+)", "Bacillariophyceae-Dinophyceae (+)", 
    "Dinophyceae-Dinophyceae (+)", "Other association (+)","Bacillariophyceae-Bacillariophyceae (-)", "Bacillariophyceae-Dinophyceae (-)", 
    "Dinophyceae-Dinophyceae (-)", "Other association (-)"
  ),name = NULL) +
  coord_polar(theta = "y") +
  xlim(c(-1, 4.6)) +
  theme_no_axes() +
  facet_wrap(~region) +
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(ncol = 7), colour = guide_legend(ncol = 7))
global_compo_asso
ggsave('asso_compo_regions_combine_withneg_2.png', path = "output/graphs/description_region/", dpi = 600, width = 300, height = 150, units = 'mm')

