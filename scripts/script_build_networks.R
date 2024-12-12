# Script Suite Stage JY Dias # 12/12/2024

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
edgelist_1 <- filter(net_cluster1$edgelist1[,c("v1","v2","asso")],asso > 0)

write.csv2(edgelist_1,file="output/tableaux/Networks/edgelist_med.csv", row.names = FALSE,dec = ".")

# Store the association matrix
assoMat_1 <- as.data.frame(net_cluster1$assoMat1)
# Keep only the positive associations
assoMat_1[assoMat_1 < 0] <- 0
assoMat_1 <- assoMat_1[, colSums(assoMat_1) != 0]
assoMat_1 <- assoMat_1[rowSums(assoMat_1) != 0, ]
write.csv2(assoMat_1,file="output/tableaux/Networks/assoMat_med.csv", row.names = FALSE,dec = ".")
taxon_list_1 <- as.data.frame(rownames(assoMat_1))
colnames(taxon_list_1) <- "Taxon"
write.csv2(taxon_list_1,file="output/tableaux/Networks/taxon_list_med.csv", row.names = FALSE,dec = ".")


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
edgelist_2 <- filter(net_cluster2$edgelist1[,c("v1","v2","asso")],asso > 0)

write.csv2(edgelist_2,file="output/tableaux/Networks/edgelist_manche.csv", row.names = FALSE,dec = ".")

# Store the association matrix
assoMat_2 <- as.data.frame(net_cluster2$assoMat1)
# Keep only the positive associations
assoMat_2[assoMat_2 < 0] <- 0
assoMat_2 <- assoMat_2[, colSums(assoMat_2) != 0]
assoMat_2 <- assoMat_2[rowSums(assoMat_2) != 0, ]
write.csv2(assoMat_2,file="output/tableaux/Networks/assoMat_manche.csv", row.names = FALSE,dec = ".")
taxon_list_2 <- as.data.frame(rownames(assoMat_2))
colnames(taxon_list_2) <- "Taxon"
write.csv2(taxon_list_2,file="output/tableaux/Networks/taxon_list_manche.csv", row.names = FALSE,dec = ".")

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
# Keep only the positive associations
edgelist_3 <- filter(net_cluster3$edgelist1[,c("v1","v2","asso")],asso > 0)

write.csv2(edgelist_3,file="output/tableaux/Networks/edgelist_atlantic.csv", row.names = FALSE,dec = ".")

# Store the association matrix
assoMat_3 <- as.data.frame(net_cluster3$assoMat1)
# Keep only the positive associations
assoMat_3[assoMat_3 < 0] <- 0
assoMat_3 <- assoMat_3[, colSums(assoMat_3) != 0]
assoMat_3 <- assoMat_3[rowSums(assoMat_3) != 0, ]
write.csv2(assoMat_3,file="output/tableaux/Networks/assoMat_atlantic.csv", row.names = FALSE,dec = ".")
taxon_list_3 <- as.data.frame(rownames(assoMat_3))
colnames(taxon_list_3) <- "Taxon"
write.csv2(taxon_list_3,file="output/tableaux/Networks/taxon_list_atlantic.csv", row.names = FALSE,dec = ".")

# Compute the temporal networks #####
### Cluster 1-Mediterranean sea ####
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_blooms_caracterised_final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt", header = TRUE, sep = "\t")

# Charger la matrice d'adjacence
assoMat <- read_delim("output/tableaux/Networks/assoMat_med.csv", 
                      delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                          grouping_mark = ""), trim_ws = TRUE)

edge <- read_delim("output/tableaux/Networks/edgelist_genus_med.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Créer le graphe
cluster1 <- graph_from_adjacency_matrix(as.matrix(assoMat), weighted = TRUE, mode = "undirected", diag = FALSE)


# Preparation  
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
  spe <- as.data.frame(CL1[i,])
  colnames(spe) <- "Count" 
  spe$Taxon <- rownames(spe)
  spe <- left_join(phyto_index,spe, by = join_by(Taxon))
  spe <- filter(spe,Count>0)
  spe$Pindex <- as.numeric(spe$Pindex)
  spe <- left_join(spe,phylum_classe)
  station <- CL1df[i,]$Code_point_Libelle
  date <- CL1df[i,]$Date
  
  if (nrow(spe) != 0){
    vids <- spe$Pindex
    sub <- igraph::subgraph(cluster1, vids)
    
    
    # Metrics 
    S_net <-  vcount(sub) # Number of nodes
    L_net <- ecount(sub) # Number of links
    Z_mat <- L_net / S_net # linkage density or average number of links per nodes
    C_net <- edge_density(sub, loops = FALSE) #connectance
    
    # Average path length
    avg_path_length <- mean_distance(
      sub,
      directed = FALSE,
      unconnected = TRUE # if the graphs is disconnected, only the existing paths are considered
    )
    
    Edge_connect <- edge.connectivity(sub) # Edge connectivity = adhesion
    
    wc <- cluster_fast_greedy(sub)
    Modularity <- modularity(sub,membership = membership(wc)) # Modularity
    
    Vert_connect <- vertex.connectivity(sub) # Vertex connectivity = adhesion
    m_degree <- mean(igraph::degree(sub)) # mean number of links
    assort <- assortativity_degree(sub,directed = F) #assortativity
    diss <- mean(1 - E(sub)$weight) # Dissilarity as defined in NetCoMi
    trans <- transitivity(sub,type = "global") #Transitivity
    mean_edge_bet <- mean(edge_betweenness(sub)) # Mean edge betweeness
    
    adj <- as.matrix(as_adjacency_matrix(sub, attr = "weight",)) # OK
    diag(adj) <- 1
    nat_connect <- natural.connectivity(as.matrix(adj)) # natural connectivity
    
    sub_edges <- filter(edge,v1 %in% spe$Taxon & v2 %in% spe$Taxon)
    
    
    data_results_reseaux[i,1] <- station
    data_results_reseaux[i,2] <- date
    data_results_reseaux[i,3] <- S_net # nombre de noeuds
    data_results_reseaux[i,4] <- L_net # nombre de liens
    data_results_reseaux[i,5] <- Z_mat # Densite des liens
    data_results_reseaux[i,6] <- C_net # connectance
    
    data_results_reseaux[i,7] <- avg_path_length # longueur moyen des liens
    data_results_reseaux[i,8] <- Edge_connect # adhesion
    data_results_reseaux[i,9] <- Modularity #modularite
    
    data_results_reseaux[i,10] <- m_degree #Nombre de liens moyens
    data_results_reseaux[i,11] <- assort # Assortativite
    data_results_reseaux[i,12] <- diss # dissimilarite
    data_results_reseaux[i,13] <- trans # transitivite
    
    data_results_reseaux[i,14] <- mean_edge_bet # Nombre moyen de voisins
    data_results_reseaux[i,15] <- nat_connect # Connectivite naturelle
    
    data_results_reseaux[i,16] <- length(wc) # Nombre de cluster
    data_results_reseaux[i,17] <- nrow(filter(spe,Phylum.Classe == "Bacillariophyceae")) # Nombre de Bacillariophyceae
    data_results_reseaux[i,18] <- nrow(filter(spe,Phylum.Classe == "Bacillariophyceae"))/nrow(spe) # Pourcentage de Bacillariophyceae
    data_results_reseaux[i,19] <- nrow(filter(spe,Phylum.Classe == "Dinophyceae")) # Nombre de Dinophyceae
    data_results_reseaux[i,20] <- nrow(filter(spe,Phylum.Classe == "Dinophyceae"))/nrow(spe) # Pourcentage de Dinophyceae
    data_results_reseaux[i,21] <- nrow(filter(spe,Phylum.Classe != "Bacillariophyceae" & Phylum.Classe != "Dinophyceae")) # Nombre d'autres taxons
    data_results_reseaux[i,22] <- nrow(filter(spe,Phylum.Classe != "Bacillariophyceae"& Phylum.Classe != "Dinophyceae"))/nrow(spe) # Pourcentage d'autres taxons
    
    data_results_reseaux[i,23] <- nrow(filter(sub_edges, link_genus == "Bacillariophyceae-Bacillariophyceae")) # Nombre de relation Bac-Bac
    data_results_reseaux[i,24] <- nrow(filter(sub_edges, link_genus == "Bacillariophyceae-Bacillariophyceae"))/nrow(sub_edges) # Pourcentage de relation Bac-Bac
    data_results_reseaux[i,25] <- nrow(filter(sub_edges, link_genus == "Bacillariophyceae-Dinophyceae")) # Nombre de relation Bac-Dino
    data_results_reseaux[i,26] <- nrow(filter(sub_edges, link_genus == "Bacillariophyceae-Dinophyceae"))/nrow(sub_edges) # Pourcentage de relation Bac-Dino
    data_results_reseaux[i,27] <- nrow(filter(sub_edges, link_genus == "Dinophyceae-Dinophyceae")) # Nombre de relation Dino-Dino
    data_results_reseaux[i,28] <- nrow(filter(sub_edges, link_genus == "Dinophyceae-Dinophyceae"))/nrow(sub_edges) # Pourcentage de relation Dino-Dino
    data_results_reseaux[i,29] <- nrow(filter(sub_edges, link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" )) # Nombre d'autres associations
    data_results_reseaux[i,30] <- nrow(filter(sub_edges,link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" ))/nrow(sub_edges) # Pourcentage d'autres relations
    
  } 
  else { 
    data_results_reseaux[i,1] <- station
    data_results_reseaux[i,2] <- date
    data_results_reseaux[i,3] <- NA # nombre de noeuds
    data_results_reseaux[i,4] <- NA # nombre de liens
    data_results_reseaux[i,5] <- NA # Densite des liens
    data_results_reseaux[i,6] <- NA # connectance
    
    data_results_reseaux[i,7] <- NA # longueur moyen des liens
    data_results_reseaux[i,8] <- NA # adhesion
    data_results_reseaux[i,9] <- NA #modularite
    
    data_results_reseaux[i,10] <- NA #Nombre de liens moyens
    data_results_reseaux[i,11] <- NA # Assortativite
    data_results_reseaux[i,12] <- NA # dissimilarite
    data_results_reseaux[i,13] <- NA # transitivite
    
    data_results_reseaux[i,14] <- NA # Nombre moyen de voisins
    data_results_reseaux[i,15] <- NA # Connectivite naturelle
    
    data_results_reseaux[i,16] <- NA # Nombre de cluster
    
    data_results_reseaux[i,17] <- NA # Nombre de bacillariophyceae
    data_results_reseaux[i,18] <- NA # Pourcentage de bacillariophyceae
    data_results_reseaux[i,19] <- NA # Nombre de dinophyceae
    data_results_reseaux[i,20] <- NA # Pourcentage de dinophyceae
    data_results_reseaux[i,21] <- NA # Nombre d'autres taxons
    data_results_reseaux[i,22] <- NA # Pourcentage d'autres taxons
    
    data_results_reseaux[i,23] <- NA # Nombre de Bac_Bac
    data_results_reseaux[i,24] <- NA # % Bac_Bac
    data_results_reseaux[i,25] <- NA # N Bac_Dino
    data_results_reseaux[i,26] <- NA # % Bac_Dino
    data_results_reseaux[i,27] <- NA # N Dino_Dino
    data_results_reseaux[i,28] <- NA # % Dino_Dino
    data_results_reseaux[i,29] <- NA # N Autres
    data_results_reseaux[i,30] <- NA # % Autres
    
  }
  
  colnames(data_results_reseaux) <- c("Code_point_Libelle","Date","N_noeuds","N_liens","D_liens","C_tance",
                                      "Avg_p_length","Adhes","Mod","meanN_liens","Assort","Diss","Trans","meanN_voisins",
                                      "Nat_connect","N_clust","N_bac","P_bac","N_dino","P_dino","N_autres","P_autres",
                                      "N_BacBac","P_BacBac","N_BacDino","P_BacDino","N_DinoDino","P_DinoDino","N_AAutres",
                                      "P_AAutres")
  
  #plot(sub,main = paste0(station,date),layout = layout_with_fr)
  #png(paste0("output/graphs/Reseaux/TS_CLUST1/",station,date,".png"))
  #plot(sub,main = paste0(station,date),layout = layout_with_fr)
  #dev.off()
  #htmlwidgets::saveWidget(viz_sub, paste0("output/graphs/Reseaux/HTML/",station,date,".html"))
  print(i/nrow(CL1)*100)
}
# Save it
write.csv2(data_results_reseaux,file="output/tableaux/Networks/subnetworks/results_metrics_reseaux_cluster1.csv", row.names = FALSE,dec = ".")

### Cluster 2-Eastern Channel - North Sea ####

data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_blooms_caracterised_final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt", header = TRUE, sep = "\t")

# Charger la matrice d'adjacence
assoMat <- read_delim("output/tableaux/Networks/assoMat_manche.csv", 
                      delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                          grouping_mark = ""), trim_ws = TRUE)

edge <- read_delim("output/tableaux/Networks/edgelist_genus_manche.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Créer le graphe
cluster2 <- graph_from_adjacency_matrix(as.matrix(assoMat), weighted = TRUE, mode = "undirected", diag = FALSE)


# Preparation  
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
  spe <- as.data.frame(CL2[i,])
  colnames(spe) <- "Count" 
  spe$Taxon <- rownames(spe)
  spe <- left_join(phyto_index,spe, by = join_by(Taxon))
  spe <- filter(spe,Count>0)
  spe$Pindex <- as.numeric(spe$Pindex)
  spe <- left_join(spe,phylum_classe)
  station <- CL2df[i,]$Code_point_Libelle
  date <- CL2df[i,]$Date
  
  if (nrow(spe) != 0){
    vids <- spe$Pindex
    sub <- igraph::subgraph(cluster2, vids)
    
    
    # Metrics 
    S_net <-  vcount(sub) # Number of nodes
    L_net <- ecount(sub) # Number of links
    Z_mat <- L_net / S_net # linkage density or average number of links per nodes
    C_net <- edge_density(sub, loops = FALSE) #connectance
    
    # Average path length
    avg_path_length <- mean_distance(
      sub,
      directed = FALSE,
      unconnected = TRUE # if the graphs is disconnected, only the existing paths are considered
    )
    
    Edge_connect <- edge.connectivity(sub) # Edge connectivity = adhesion
    
    wc <- cluster_fast_greedy(sub)
    Modularity <- modularity(sub,membership = membership(wc)) # Modularity
    
    Vert_connect <- vertex.connectivity(sub) # Vertex connectivity = adhesion
    m_degree <- mean(igraph::degree(sub)) # mean number of links
    assort <- assortativity_degree(sub,directed = F) #assortativity
    diss <- mean(1 - E(sub)$weight) # Dissilarity as defined in NetCoMi
    trans <- transitivity(sub,type = "global") #Transitivity
    mean_edge_bet <- mean(edge_betweenness(sub)) # Mean edge betweeness
    
    adj <- as.matrix(as_adjacency_matrix(sub, attr = "weight",)) # OK
    diag(adj) <- 1
    nat_connect <- natural.connectivity(as.matrix(adj)) # natural connectivity
    
    sub_edges <- filter(edge,v1 %in% spe$Taxon & v2 %in% spe$Taxon)
    
    
    data_results_reseaux[i,1] <- station
    data_results_reseaux[i,2] <- date
    data_results_reseaux[i,3] <- S_net # nombre de noeuds
    data_results_reseaux[i,4] <- L_net # nombre de liens
    data_results_reseaux[i,5] <- Z_mat # Densite des liens
    data_results_reseaux[i,6] <- C_net # connectance
    
    data_results_reseaux[i,7] <- avg_path_length # longueur moyen des liens
    data_results_reseaux[i,8] <- Edge_connect # adhesion
    data_results_reseaux[i,9] <- Modularity #modularite
    
    data_results_reseaux[i,10] <- m_degree #Nombre de liens moyens
    data_results_reseaux[i,11] <- assort # Assortativite
    data_results_reseaux[i,12] <- diss # dissimilarite
    data_results_reseaux[i,13] <- trans # transitivite
    
    data_results_reseaux[i,14] <- mean_edge_bet # Nombre moyen de voisins
    data_results_reseaux[i,15] <- nat_connect # Connectivite naturelle
    
    data_results_reseaux[i,16] <- length(wc) # Nombre de cluster
    data_results_reseaux[i,17] <- nrow(filter(spe,Phylum.Classe == "Bacillariophyceae")) # Nombre de Bacillariophyceae
    data_results_reseaux[i,18] <- nrow(filter(spe,Phylum.Classe == "Bacillariophyceae"))/nrow(spe) # Pourcentage de Bacillariophyceae
    data_results_reseaux[i,19] <- nrow(filter(spe,Phylum.Classe == "Dinophyceae")) # Nombre de Dinophyceae
    data_results_reseaux[i,20] <- nrow(filter(spe,Phylum.Classe == "Dinophyceae"))/nrow(spe) # Pourcentage de Dinophyceae
    data_results_reseaux[i,21] <- nrow(filter(spe,Phylum.Classe != "Bacillariophyceae" & Phylum.Classe != "Dinophyceae")) # Nombre d'autres taxons
    data_results_reseaux[i,22] <- nrow(filter(spe,Phylum.Classe != "Bacillariophyceae"& Phylum.Classe != "Dinophyceae"))/nrow(spe) # Pourcentage d'autres taxons
    
    data_results_reseaux[i,23] <- nrow(filter(sub_edges, link_genus == "Bacillariophyceae-Bacillariophyceae")) # Nombre de relation Bac-Bac
    data_results_reseaux[i,24] <- nrow(filter(sub_edges, link_genus == "Bacillariophyceae-Bacillariophyceae"))/nrow(sub_edges) # Pourcentage de relation Bac-Bac
    data_results_reseaux[i,25] <- nrow(filter(sub_edges, link_genus == "Bacillariophyceae-Dinophyceae")) # Nombre de relation Bac-Dino
    data_results_reseaux[i,26] <- nrow(filter(sub_edges, link_genus == "Bacillariophyceae-Dinophyceae"))/nrow(sub_edges) # Pourcentage de relation Bac-Dino
    data_results_reseaux[i,27] <- nrow(filter(sub_edges, link_genus == "Dinophyceae-Dinophyceae")) # Nombre de relation Dino-Dino
    data_results_reseaux[i,28] <- nrow(filter(sub_edges, link_genus == "Dinophyceae-Dinophyceae"))/nrow(sub_edges) # Pourcentage de relation Dino-Dino
    data_results_reseaux[i,29] <- nrow(filter(sub_edges, link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" )) # Nombre d'autres associations
    data_results_reseaux[i,30] <- nrow(filter(sub_edges,link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" ))/nrow(sub_edges) # Pourcentage d'autres relations
    
  } 
  else { 
    data_results_reseaux[i,1] <- station
    data_results_reseaux[i,2] <- date
    data_results_reseaux[i,3] <- NA # nombre de noeuds
    data_results_reseaux[i,4] <- NA # nombre de liens
    data_results_reseaux[i,5] <- NA # Densite des liens
    data_results_reseaux[i,6] <- NA # connectance
    
    data_results_reseaux[i,7] <- NA # longueur moyen des liens
    data_results_reseaux[i,8] <- NA # adhesion
    data_results_reseaux[i,9] <- NA #modularite
    
    data_results_reseaux[i,10] <- NA #Nombre de liens moyens
    data_results_reseaux[i,11] <- NA # Assortativite
    data_results_reseaux[i,12] <- NA # dissimilarite
    data_results_reseaux[i,13] <- NA # transitivite
    
    data_results_reseaux[i,14] <- NA # Nombre moyen de voisins
    data_results_reseaux[i,15] <- NA # Connectivite naturelle
    
    data_results_reseaux[i,16] <- NA # Nombre de cluster
    
    data_results_reseaux[i,17] <- NA # Nombre de bacillariophyceae
    data_results_reseaux[i,18] <- NA # Pourcentage de bacillariophyceae
    data_results_reseaux[i,19] <- NA # Nombre de dinophyceae
    data_results_reseaux[i,20] <- NA # Pourcentage de dinophyceae
    data_results_reseaux[i,21] <- NA # Nombre d'autres taxons
    data_results_reseaux[i,22] <- NA # Pourcentage d'autres taxons
    
    data_results_reseaux[i,23] <- NA # Nombre de Bac_Bac
    data_results_reseaux[i,24] <- NA # % Bac_Bac
    data_results_reseaux[i,25] <- NA # N Bac_Dino
    data_results_reseaux[i,26] <- NA # % Bac_Dino
    data_results_reseaux[i,27] <- NA # N Dino_Dino
    data_results_reseaux[i,28] <- NA # % Dino_Dino
    data_results_reseaux[i,29] <- NA # N Autres
    data_results_reseaux[i,30] <- NA # % Autres
    
  }
  
  colnames(data_results_reseaux) <- c("Code_point_Libelle","Date","N_noeuds","N_liens","D_liens","C_tance",
                                      "Avg_p_length","Adhes","Mod","meanN_liens","Assort","Diss","Trans","meanN_voisins",
                                      "Nat_connect","N_clust","N_bac","P_bac","N_dino","P_dino","N_autres","P_autres",
                                      "N_BacBac","P_BacBac","N_BacDino","P_BacDino","N_DinoDino","P_DinoDino","N_AAutres",
                                      "P_AAutres")
  
  #plot(sub,main = paste0(station,date),layout = layout_with_fr)
  #png(paste0("output/graphs/Reseaux/TS_CLUST1/",station,date,".png"))
  #plot(sub,main = paste0(station,date),layout = layout_with_fr)
  #dev.off()
  #htmlwidgets::saveWidget(viz_sub, paste0("output/graphs/Reseaux/HTML/",station,date,".html"))
  print(i/nrow(CL2)*100)
}
# Save it
write.csv2(data_results_reseaux,file="output/tableaux/Networks/subnetworks/results_metrics_reseaux_cluster2.csv", row.names = FALSE,dec = ".")


### Cluster 3-Atlantic - Western Channel ####
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_blooms_caracterised_final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt", header = TRUE, sep = "\t")

# Charger la matrice d'adjacence
assoMat <- read_delim("output/tableaux/Networks/assoMat_atlantic.csv", 
                      delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                          grouping_mark = ""), trim_ws = TRUE)

edge <- read_delim("output/tableaux/Networks/edgelist_genus_atlantic.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Créer le graphe
cluster3 <- graph_from_adjacency_matrix(as.matrix(assoMat), weighted = TRUE, mode = "undirected", diag = FALSE)


# Preparation  
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
  spe <- as.data.frame(CL3[i,])
  colnames(spe) <- "Count" 
  spe$Taxon <- rownames(spe)
  spe <- left_join(phyto_index,spe, by = join_by(Taxon))
  spe <- filter(spe,Count>0)
  spe$Pindex <- as.numeric(spe$Pindex)
  spe <- left_join(spe,phylum_classe)
  station <- CL3df[i,]$Code_point_Libelle
  date <- CL3df[i,]$Date
  
  if (nrow(spe) != 0){
    vids <- spe$Pindex
    sub <- igraph::subgraph(cluster3, vids)
    
    
    # Metrics 
    S_net <-  vcount(sub) # Number of nodes
    L_net <- ecount(sub) # Number of links
    Z_mat <- L_net / S_net # linkage density or average number of links per nodes
    C_net <- edge_density(sub, loops = FALSE) #connectance
    
    # Average path length
    avg_path_length <- mean_distance(
      sub,
      directed = FALSE,
      unconnected = TRUE # if the graphs is disconnected, only the existing paths are considered
    )
    
    Edge_connect <- edge.connectivity(sub) # Edge connectivity = adhesion
    
    wc <- cluster_fast_greedy(sub)
    Modularity <- modularity(sub,membership = membership(wc)) # Modularity
    
    Vert_connect <- vertex.connectivity(sub) # Vertex connectivity = adhesion
    m_degree <- mean(igraph::degree(sub)) # mean number of links
    assort <- assortativity_degree(sub,directed = F) #assortativity
    diss <- mean(1 - E(sub)$weight) # Dissilarity as defined in NetCoMi
    trans <- transitivity(sub,type = "global") #Transitivity
    mean_edge_bet <- mean(edge_betweenness(sub)) # Mean edge betweeness
    
    adj <- as.matrix(as_adjacency_matrix(sub, attr = "weight",)) # OK
    diag(adj) <- 1
    nat_connect <- natural.connectivity(as.matrix(adj)) # natural connectivity
    
    sub_edges <- filter(edge,v1 %in% spe$Taxon & v2 %in% spe$Taxon)
    
    
    data_results_reseaux[i,1] <- station
    data_results_reseaux[i,2] <- date
    data_results_reseaux[i,3] <- S_net # nombre de noeuds
    data_results_reseaux[i,4] <- L_net # nombre de liens
    data_results_reseaux[i,5] <- Z_mat # Densite des liens
    data_results_reseaux[i,6] <- C_net # connectance
    
    data_results_reseaux[i,7] <- avg_path_length # longueur moyen des liens
    data_results_reseaux[i,8] <- Edge_connect # adhesion
    data_results_reseaux[i,9] <- Modularity #modularite
    
    data_results_reseaux[i,10] <- m_degree #Nombre de liens moyens
    data_results_reseaux[i,11] <- assort # Assortativite
    data_results_reseaux[i,12] <- diss # dissimilarite
    data_results_reseaux[i,13] <- trans # transitivite
    
    data_results_reseaux[i,14] <- mean_edge_bet # Nombre moyen de voisins
    data_results_reseaux[i,15] <- nat_connect # Connectivite naturelle
    
    data_results_reseaux[i,16] <- length(wc) # Nombre de cluster
    data_results_reseaux[i,17] <- nrow(filter(spe,Phylum.Classe == "Bacillariophyceae")) # Nombre de Bacillariophyceae
    data_results_reseaux[i,18] <- nrow(filter(spe,Phylum.Classe == "Bacillariophyceae"))/nrow(spe) # Pourcentage de Bacillariophyceae
    data_results_reseaux[i,19] <- nrow(filter(spe,Phylum.Classe == "Dinophyceae")) # Nombre de Dinophyceae
    data_results_reseaux[i,20] <- nrow(filter(spe,Phylum.Classe == "Dinophyceae"))/nrow(spe) # Pourcentage de Dinophyceae
    data_results_reseaux[i,21] <- nrow(filter(spe,Phylum.Classe != "Bacillariophyceae" & Phylum.Classe != "Dinophyceae")) # Nombre d'autres taxons
    data_results_reseaux[i,22] <- nrow(filter(spe,Phylum.Classe != "Bacillariophyceae"& Phylum.Classe != "Dinophyceae"))/nrow(spe) # Pourcentage d'autres taxons
    
    data_results_reseaux[i,23] <- nrow(filter(sub_edges, link_genus == "Bacillariophyceae-Bacillariophyceae")) # Nombre de relation Bac-Bac
    data_results_reseaux[i,24] <- nrow(filter(sub_edges, link_genus == "Bacillariophyceae-Bacillariophyceae"))/nrow(sub_edges) # Pourcentage de relation Bac-Bac
    data_results_reseaux[i,25] <- nrow(filter(sub_edges, link_genus == "Bacillariophyceae-Dinophyceae")) # Nombre de relation Bac-Dino
    data_results_reseaux[i,26] <- nrow(filter(sub_edges, link_genus == "Bacillariophyceae-Dinophyceae"))/nrow(sub_edges) # Pourcentage de relation Bac-Dino
    data_results_reseaux[i,27] <- nrow(filter(sub_edges, link_genus == "Dinophyceae-Dinophyceae")) # Nombre de relation Dino-Dino
    data_results_reseaux[i,28] <- nrow(filter(sub_edges, link_genus == "Dinophyceae-Dinophyceae"))/nrow(sub_edges) # Pourcentage de relation Dino-Dino
    data_results_reseaux[i,29] <- nrow(filter(sub_edges, link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" )) # Nombre d'autres associations
    data_results_reseaux[i,30] <- nrow(filter(sub_edges,link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" ))/nrow(sub_edges) # Pourcentage d'autres relations
    
  } 
  else { 
    data_results_reseaux[i,1] <- station
    data_results_reseaux[i,2] <- date
    data_results_reseaux[i,3] <- NA # nombre de noeuds
    data_results_reseaux[i,4] <- NA # nombre de liens
    data_results_reseaux[i,5] <- NA # Densite des liens
    data_results_reseaux[i,6] <- NA # connectance
    
    data_results_reseaux[i,7] <- NA # longueur moyen des liens
    data_results_reseaux[i,8] <- NA # adhesion
    data_results_reseaux[i,9] <- NA #modularite
    
    data_results_reseaux[i,10] <- NA #Nombre de liens moyens
    data_results_reseaux[i,11] <- NA # Assortativite
    data_results_reseaux[i,12] <- NA # dissimilarite
    data_results_reseaux[i,13] <- NA # transitivite
    
    data_results_reseaux[i,14] <- NA # Nombre moyen de voisins
    data_results_reseaux[i,15] <- NA # Connectivite naturelle
    
    data_results_reseaux[i,16] <- NA # Nombre de cluster
    
    data_results_reseaux[i,17] <- NA # Nombre de bacillariophyceae
    data_results_reseaux[i,18] <- NA # Pourcentage de bacillariophyceae
    data_results_reseaux[i,19] <- NA # Nombre de dinophyceae
    data_results_reseaux[i,20] <- NA # Pourcentage de dinophyceae
    data_results_reseaux[i,21] <- NA # Nombre d'autres taxons
    data_results_reseaux[i,22] <- NA # Pourcentage d'autres taxons
    
    data_results_reseaux[i,23] <- NA # Nombre de Bac_Bac
    data_results_reseaux[i,24] <- NA # % Bac_Bac
    data_results_reseaux[i,25] <- NA # N Bac_Dino
    data_results_reseaux[i,26] <- NA # % Bac_Dino
    data_results_reseaux[i,27] <- NA # N Dino_Dino
    data_results_reseaux[i,28] <- NA # % Dino_Dino
    data_results_reseaux[i,29] <- NA # N Autres
    data_results_reseaux[i,30] <- NA # % Autres
    
  }
  
  colnames(data_results_reseaux) <- c("Code_point_Libelle","Date","N_noeuds","N_liens","D_liens","C_tance",
                                      "Avg_p_length","Adhes","Mod","meanN_liens","Assort","Diss","Trans","meanN_voisins",
                                      "Nat_connect","N_clust","N_bac","P_bac","N_dino","P_dino","N_autres","P_autres",
                                      "N_BacBac","P_BacBac","N_BacDino","P_BacDino","N_DinoDino","P_DinoDino","N_AAutres",
                                      "P_AAutres")
  
  #plot(sub,main = paste0(station,date),layout = layout_with_fr)
  #png(paste0("output/graphs/Reseaux/TS_CLUST1/",station,date,".png"))
  #plot(sub,main = paste0(station,date),layout = layout_with_fr)
  #dev.off()
  #htmlwidgets::saveWidget(viz_sub, paste0("output/graphs/Reseaux/HTML/",station,date,".html"))
  print(i/nrow(CL3)*100)
}
# Save it
write.csv2(data_results_reseaux,file="output/tableaux/Networks/subnetworks/results_metrics_reseaux_cluster3.csv", row.names = FALSE,dec = ".")


