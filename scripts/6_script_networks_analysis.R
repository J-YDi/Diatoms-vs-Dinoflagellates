################################################################################
# Diatoms vs dinoflagellates: a network analysis of bloom impacts on diversity #
#                    and phytoplankton associations | R scripts                #
################################################################################

# Script to analyze composition and associations type in the networks and make 
# a PCA to sum up the different metrics #
# 02/28/2025

# Load packages
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)
library(igraph)
library(ggraph)
library(FactoMineR)
library(factoextra)
library(vegan)
library(corrplot)
library(ggforce)

### Working on the Mediterranean sea #####
# Import data
data <- read_delim("output/tableaux/Networks/edgelist_med.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Load table "Liste_phylum.classe" to associate species to a class
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt",header = T,sep = "\t")

# Create a new dataframe that indicate the association type
# Working on the edges
edge1 <- as.data.frame(data$v1)
colnames(edge1) <- "Taxon"
edge1 <- left_join(edge1,phylum_classe)
colnames(edge1) <- c("v1","v1_classe")

edge2 <- as.data.frame(data$v2)
colnames(edge2) <- "Taxon"
edge2 <- left_join(edge2,phylum_classe)
colnames(edge2) <- c("v2","v2_classe")

edgelist_genus_med <- cbind(edge1,edge2,data$asso)
# Associate the class
edgelist_genus_med$link_genus <- paste0(edgelist_genus_med$v1_classe,"-",edgelist_genus_med$v2_classe) 
# Function to normalize the association type, that's make Dinophyceae-Bacillariophyceae => Bacillariophyceae-Dinophyceae
normaliser_paire <- function(x) {
  parties <- sort(unlist(strsplit(x, "-")))
  paste(parties, collapse = "-")
}
# Applying it 
edgelist_genus_med$link_genus <- sapply(edgelist_genus_med$link_genus, normaliser_paire)
# Save it
#write.csv2(edgelist_genus_med,file="output/tableaux/Networks/edgelist_genus_med.csv", row.names = FALSE,dec = ".")

# Frequence of the association analysis
freq_table_med <-  edgelist_genus_med %>%
  rowwise() %>%
  mutate(Combination = paste(sort(c(v1_classe, v2_classe)), collapse = "-")) %>%  
  ungroup() %>%
  group_by(Combination) %>%
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
#write.csv2(freq_table_med,file="output/tableaux/Networks/freq_table_med.csv", row.names = FALSE,dec = ".")


# Load the association matrix
assoMat <- read_delim("output/tableaux/Networks/assoMat_med.csv", 
                      delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                          grouping_mark = ""), trim_ws = TRUE)

# Create the graph
cluster1 <- graph_from_adjacency_matrix(as.matrix(assoMat), weighted = TRUE, mode = "undirected", diag = FALSE)

# Clustering
wcini <- cluster_fast_greedy(cluster1)

# Load table "Liste_phylum.classe" to associate species to a class
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt", header = TRUE, sep = "\t")

# Associate the taxa to the class
compo_reseau <- data.frame(Taxon = wcini$names, cluster = wcini$membership) %>%
  left_join(phylum_classe, by = c("Taxon" = "Taxon")) %>%
  mutate(color = case_when(
    Phylum.Classe == "Bacillariophyceae" ~ "#1f77b4",
    Phylum.Classe == "Dinophyceae" ~ "green",
    Phylum.Classe == "Ciliophora" ~ "grey",
    Phylum.Classe == "Cryptophyceae" ~ "grey",
    Phylum.Classe == "Haptophyta" ~ "grey",
    TRUE ~ "grey"
  ))

# Import the edges list
edgelist_genus_med <- read_delim("output/tableaux/Networks/edgelist_genus_med.csv", 
                      delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                          grouping_mark = ""), trim_ws = TRUE)
# Adding some customization info for ggraph
V(cluster1)$Phylum <- compo_reseau$Phylum.Classe #class
V(cluster1)$size <- 15 # size of the node
E(cluster1)$weight <- E(cluster1)$weight * 5 #edges size proportional to the correlation
V(cluster1)$label <- V(cluster1)$name # name of the node
V(cluster1)$color <- compo_reseau$color # color of the node 
E(cluster1)$info <- edgelist_genus_med$link_genus #type of association


# making the global network plot
med_global <- ggraph(cluster1, layout = "auto") +  
  geom_edge_link(aes(edge_width = weight,edge_colour = info), alpha = 0.5, lineend = "round", linejoin = "bevel") +
  geom_node_point(aes(size = size, fill = color, stroke = 1), shape = 21,alpha = 0.5) +
  #geom_node_text(aes(label = label, fontface = "bold"),colour = "black", 
  #               repel = TRUE, size = 3, check_overlap = TRUE,
  #               vjust = 0.5, hjust = 1) + 
  scale_edge_width(range = c(0.1, 2)) +
  scale_size_continuous(range = c(3, 10)) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c(
      "#1f77b4" = "#1f77b4",
      "green" = "green",
      "grey" = "grey"
    )
  )+
  scale_edge_colour_manual(values = c(
    "Bacillariophyceae-Bacillariophyceae" = "#56B4E9","Bacillariophyceae-Dinophyceae"   = "chocolate1"    ,"Bacillariophyceae-Euglenozoa" = "grey30"
   ,"Bacillariophyceae-Dictyochophyceae"  = "grey30","Bacillariophyceae-Chlorophyta"   = "grey30"   ,"Bacillariophyceae-Chrysophyceae"    = "grey30"
   ,"Bacillariophyceae-Cyanophyceae"      = "grey30","Bacillariophyceae-Haptophyta"    = "grey30"   ,"Dinophyceae-Dinophyceae" = "#009E73"
   , "Dinophyceae-Euglenozoa"             = "grey30", "Dinophyceae-Ebriophyceae"       = "grey30"   , "Bacillariophyceae-Ebriophyceae"    = "grey30"
   , "Bacillariophyceae-Cryptophyceae"    = "grey30", "Dinophyceae-Haptophyta"         = "grey30"   , "Chrysophyceae-Dinophyceae"         = "grey30"
   , "Chlorophyta-Dinophyceae"            = "grey30", "Cryptophyceae-Dinophyceae"      = "grey30"   , "Cyanophyceae-Dinophyceae"          = "grey30"
   , "Chlorophyta-Euglenozoa"             = "grey30", "Cryptophyceae-Euglenozoa"       = "grey30"   , "Chrysophyceae-Euglenozoa"          = "grey30"
   , "Cyanophyceae-Euglenozoa"            = "grey30", "Chlorophyta-Dictyochophyceae"   = "grey30"   , "Dictyochophyceae-Dinophyceae"      = "grey30"
   , "Cryptophyceae-Dictyochophyceae"     = "grey30", "Dictyochophyceae-Haptophyta"    = "grey30"   , "Dictyochophyceae-Euglenozoa"       = "grey30"
   , "Chrysophyceae-Dictyochophyceae"     = "grey30", "Cyanophyceae-Dictyochophyceae"  = "grey30"   , "Chlorophyta-Cryptophyceae"         = "grey30"
   , "Chlorophyta-Haptophyta"             = "grey30", "Chlorophyta-Chrysophyceae"      = "grey30"   , "Chlorophyta-Cyanophyceae"          = "grey30"
   , "Cryptophyceae-Haptophyta"           = "grey30", "Chrysophyceae-Cryptophyceae"    = "grey30"   , "Cryptophyceae-Cyanophyceae"        = "grey30"
   , "Euglenozoa-Haptophyta"              = "grey30", "Chrysophyceae-Haptophyta"       = "grey30"   , "Ebriophyceae-Haptophyta"           = "grey30"
   , "Cyanophyceae-Haptophyta"            = "grey30", "Ebriophyceae-Euglenozoa"        = "grey30"   , "Euglenozoa-Euglenozoa"             = "grey30"
   , "Chrysophyceae-Cyanophyceae"         = "grey30", "Ebriophyceae-Ebriophyceae" = "grey30"
  ))
med_global
#ggsave('Med_global_network.png', path = "output/graphs/networks/global_networks/", dpi = 600, width = 400, height = 300, units = 'mm')


### Working on the English channel #####
# Import data
data <- read_delim("output/tableaux/Networks/edgelist_manche.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)
# Load table "Liste_phylum.classe" to associate species to a class
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt",header = T,sep = "\t")

# Create a new dataframe that indicate the association type
# Working on the edges
edge1 <- as.data.frame(data$v1)
colnames(edge1) <- "Taxon"
edge1 <- left_join(edge1,phylum_classe)
colnames(edge1) <- c("v1","v1_classe")

edge2 <- as.data.frame(data$v2)
colnames(edge2) <- "Taxon"
edge2 <- left_join(edge2,phylum_classe)
colnames(edge2) <- c("v2","v2_classe")

edgelist_genus_manche <- cbind(edge1,edge2,data$asso)
# Associate the class
edgelist_genus_manche$link_genus <- paste0(edgelist_genus_manche$v1_classe,"-",edgelist_genus_manche$v2_classe) 

# Function to normalize the association type, that's make Dinophyceae-Bacillariophyceae => Bacillariophyceae-Dinophyceae
normaliser_paire <- function(x) {
  parties <- sort(unlist(strsplit(x, "-")))
  paste(parties, collapse = "-")
}
# Applying it 
edgelist_genus_manche$link_genus <- sapply(edgelist_genus_manche$link_genus, normaliser_paire)
# Save it
#write.csv2(edgelist_genus_manche,file="output/tableaux/Networks/edgelist_genus_manche.csv", row.names = FALSE,dec = ".")

# Frequence of the association analysis
freq_table_manche <-  edgelist_genus_manche %>%
  rowwise() %>%
  mutate(Combination = paste(sort(c(v1_classe, v2_classe)), collapse = "-")) %>%  
  ungroup() %>%
  group_by(Combination) %>%
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
#write.csv2(freq_table_manche,file="output/tableaux/Networks/freq_table_manche.csv", row.names = FALSE,dec = ".")

# Load the association matrix
assoMat <- read_delim("output/tableaux/Networks/assoMat_manche.csv", 
                      delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                          grouping_mark = ""), trim_ws = TRUE)

# Create the graph
cluster2 <- graph_from_adjacency_matrix(as.matrix(assoMat), weighted = TRUE, mode = "undirected", diag = FALSE)

# Clustering
wcini <- cluster_fast_greedy(cluster2)

# Load table "Liste_phylum.classe" to associate species to a class
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt", header = TRUE, sep = "\t")

# Associate the taxa to the class
compo_reseau <- data.frame(Taxon = wcini$names, cluster = wcini$membership) %>%
  left_join(phylum_classe, by = c("Taxon" = "Taxon")) %>%
  mutate(color = case_when(
    Phylum.Classe == "Bacillariophyceae" ~ "#1f77b4",
    Phylum.Classe == "Dinophyceae" ~ "green",
    Phylum.Classe == "Ciliophora" ~ "grey",
    Phylum.Classe == "Cryptophyceae" ~ "grey",
    Phylum.Classe == "Haptophyta" ~ "grey",
    TRUE ~ "grey"
  ))

# Import the edges list
edgelist_genus_med <- read_delim("output/tableaux/Networks/edgelist_genus_manche.csv", 
                                 delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                                     grouping_mark = ""), trim_ws = TRUE)
# Adding some customization info for ggraph
V(cluster2)$Phylum <- compo_reseau$Phylum.Classe # classe
V(cluster2)$size <- 15 # size of the node
E(cluster2)$weight <- E(cluster2)$weight * 5  # edges size proportional to the correlation
V(cluster2)$label <- V(cluster2)$name # mode's name
V(cluster2)$color <- compo_reseau$color # color of the node
E(cluster2)$info <- edgelist_genus_med$link_genus # type of association


# making the global network plot
manche_global <- ggraph(cluster2, layout = "auto") +  
  geom_edge_link(aes(edge_width = weight,edge_colour = info), alpha = 0.5, lineend = "round", linejoin = "bevel") +
  geom_node_point(aes(size = size, fill = color, stroke = 1), shape = 21,alpha = 0.5) +
  #geom_node_text(aes(label = label, fontface = "bold"),colour = "black", 
  #               repel = TRUE, size = 3, check_overlap = TRUE,
  #               vjust = 0.5, hjust = 1) + 
  scale_edge_width(range = c(0.1, 2)) +
  scale_size_continuous(range = c(3, 10)) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c(
      "#1f77b4" = "#1f77b4",
      "green" = "green",
      "grey" = "grey"
    )
  )+
  scale_edge_colour_manual(values = c(
    "Bacillariophyceae-Bacillariophyceae" = "#56B4E9","Bacillariophyceae-Dinophyceae"   = "chocolate1"    ,"Bacillariophyceae-Euglenozoa" = "grey30"
    ,"Bacillariophyceae-Dictyochophyceae"  = "grey30","Bacillariophyceae-Chlorophyta"   = "grey30"   ,"Bacillariophyceae-Chrysophyceae"    = "grey30"
    ,"Bacillariophyceae-Cyanophyceae"      = "grey30","Bacillariophyceae-Haptophyta"    = "grey30"   ,"Dinophyceae-Dinophyceae" = "#009E73"
    , "Dinophyceae-Euglenozoa"             = "grey30", "Dinophyceae-Ebriophyceae"       = "grey30"   , "Bacillariophyceae-Ebriophyceae"    = "grey30"
    , "Bacillariophyceae-Cryptophyceae"    = "grey30", "Dinophyceae-Haptophyta"         = "grey30"   , "Chrysophyceae-Dinophyceae"         = "grey30"
    , "Chlorophyta-Dinophyceae"            = "grey30", "Cryptophyceae-Dinophyceae"      = "grey30"   , "Cyanophyceae-Dinophyceae"          = "grey30"
    , "Chlorophyta-Euglenozoa"             = "grey30", "Cryptophyceae-Euglenozoa"       = "grey30"   , "Chrysophyceae-Euglenozoa"          = "grey30"
    , "Cyanophyceae-Euglenozoa"            = "grey30", "Chlorophyta-Dictyochophyceae"   = "grey30"   , "Dictyochophyceae-Dinophyceae"      = "grey30"
    , "Cryptophyceae-Dictyochophyceae"     = "grey30", "Dictyochophyceae-Haptophyta"    = "grey30"   , "Dictyochophyceae-Euglenozoa"       = "grey30"
    , "Chrysophyceae-Dictyochophyceae"     = "grey30", "Cyanophyceae-Dictyochophyceae"  = "grey30"   , "Chlorophyta-Cryptophyceae"         = "grey30"
    , "Chlorophyta-Haptophyta"             = "grey30", "Chlorophyta-Chrysophyceae"      = "grey30"   , "Chlorophyta-Cyanophyceae"          = "grey30"
    , "Cryptophyceae-Haptophyta"           = "grey30", "Chrysophyceae-Cryptophyceae"    = "grey30"   , "Cryptophyceae-Cyanophyceae"        = "grey30"
    , "Euglenozoa-Haptophyta"              = "grey30", "Chrysophyceae-Haptophyta"       = "grey30"   , "Ebriophyceae-Haptophyta"           = "grey30"
    , "Cyanophyceae-Haptophyta"            = "grey30", "Ebriophyceae-Euglenozoa"        = "grey30"   , "Euglenozoa-Euglenozoa"             = "grey30"
    , "Chrysophyceae-Cyanophyceae"         = "grey30", "Ebriophyceae-Ebriophyceae" = "grey30"
  ))

manche_global 
#ggsave('Manche_global_network.png', path = "output/graphs/networks/global_networks/", dpi = 600, width = 400, height = 300, units = 'mm')


### Working on the Atlantic #####
# Import data
data <- read_delim("output/tableaux/Networks/edgelist_Atlantic.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)
# Load table "Liste_phylum.classe" to associate species to a class
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt",header = T,sep = "\t")

# Create a new dataframe that indicate the association type
# Working on the edges
edge1 <- as.data.frame(data$v1)
colnames(edge1) <- "Taxon"
edge1 <- left_join(edge1,phylum_classe)
colnames(edge1) <- c("v1","v1_classe")

edge2 <- as.data.frame(data$v2)
colnames(edge2) <- "Taxon"
edge2 <- left_join(edge2,phylum_classe)
colnames(edge2) <- c("v2","v2_classe")

edgelist_genus_atlantic <- cbind(edge1,edge2,data$asso)
# Associate the class
edgelist_genus_atlantic$link_genus <- paste0(edgelist_genus_atlantic$v1_classe,"-",edgelist_genus_atlantic$v2_classe) 

# Function to normalize the association type, that's make Dinophyceae-Bacillariophyceae => Bacillariophyceae-Dinophyceae
normaliser_paire <- function(x) {
  parties <- sort(unlist(strsplit(x, "-")))
  paste(parties, collapse = "-")
}
# Applying it
edgelist_genus_atlantic$link_genus <- sapply(edgelist_genus_atlantic$link_genus, normaliser_paire)
# Save it
#write.csv2(edgelist_genus_atlantic,file="output/tableaux/Networks/edgelist_genus_atlantic.csv", row.names = FALSE,dec = ".")

# Frequence of the association analysis
freq_table_atlantic <-  edgelist_genus_atlantic %>%
  rowwise() %>%
  mutate(Combination = paste(sort(c(v1_classe, v2_classe)), collapse = "-")) %>%  
  ungroup() %>%
  group_by(Combination) %>%
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
#write.csv2(freq_table_atlantic,file="output/tableaux/Networks/freq_table_atlantic.csv", row.names = FALSE,dec = ".")

# Load the association matrix
assoMat <- read_delim("output/tableaux/Networks/assoMat_atlantic.csv", 
                      delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                          grouping_mark = ""), trim_ws = TRUE)

# Create the graph
cluster3 <- graph_from_adjacency_matrix(as.matrix(assoMat), weighted = TRUE, mode = "undirected", diag = FALSE)

# Clustering 
wcini <- cluster_fast_greedy(cluster3)

# Load table "Liste_phylum.classe" to associate species to a class
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt", header = TRUE, sep = "\t")

# Associate the taxa to the class
compo_reseau <- data.frame(Taxon = wcini$names, cluster = wcini$membership) %>%
  left_join(phylum_classe, by = c("Taxon" = "Taxon")) %>%
  mutate(color = case_when(
    Phylum.Classe == "Bacillariophyceae" ~ "#1f77b4",
    Phylum.Classe == "Dinophyceae" ~ "green",
    Phylum.Classe == "Ciliophora" ~ "grey",
    Phylum.Classe == "Cryptophyceae" ~ "grey",
    Phylum.Classe == "Haptophyta" ~ "grey",
    TRUE ~ "grey"
  ))

# Import the edges list
edgelist_genus_med <- read_delim("output/tableaux/Networks/edgelist_genus_atlantic.csv", 
                                 delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                                     grouping_mark = ""), trim_ws = TRUE)
# Adding some customization info for ggraph
V(cluster3)$Phylum <- compo_reseau$Phylum.Classe # class
V(cluster3)$size <- 15 # node's size
E(cluster3)$weight <- E(cluster3)$weight * 5  #edges size proportional to the correlation 
V(cluster3)$label <- V(cluster3)$name #node's name
V(cluster3)$color <- compo_reseau$color # color of the node
E(cluster3)$info <- edgelist_genus_med$link_genus # type of association


# making the global network plot
atlantic_global <-ggraph(cluster3, layout = "auto") +  
  geom_edge_link(aes(edge_width = weight,edge_colour = info), alpha = 0.5, lineend = "round", linejoin = "bevel") +
  geom_node_point(aes(size = size, fill = color, stroke = 1), shape = 21,alpha = 0.5) +
  #geom_node_text(aes(label = label, fontface = "bold"),colour = "black", 
  #               repel = TRUE, size = 3, check_overlap = TRUE,
  #               vjust = 0.5, hjust = 1) + 
  scale_edge_width(range = c(0.1, 2)) +
  scale_size_continuous(range = c(3, 10)) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c(
      "#1f77b4" = "#1f77b4",
      "green" = "green",
      "grey" = "grey"
    )
  )+
  scale_edge_colour_manual(values = c(
    "Bacillariophyceae-Bacillariophyceae" = "#56B4E9","Bacillariophyceae-Dinophyceae"   = "chocolate1"    ,"Bacillariophyceae-Euglenozoa" = "grey30"
    ,"Bacillariophyceae-Dictyochophyceae"  = "grey30","Bacillariophyceae-Chlorophyta"   = "grey30"   ,"Bacillariophyceae-Chrysophyceae"    = "grey30"
    ,"Bacillariophyceae-Cyanophyceae"      = "grey30","Bacillariophyceae-Haptophyta"    = "grey30"   ,"Dinophyceae-Dinophyceae" = "#009E73"
    , "Dinophyceae-Euglenozoa"             = "grey30", "Dinophyceae-Ebriophyceae"       = "grey30"   , "Bacillariophyceae-Ebriophyceae"    = "grey30"
    , "Bacillariophyceae-Cryptophyceae"    = "grey30", "Dinophyceae-Haptophyta"         = "grey30"   , "Chrysophyceae-Dinophyceae"         = "grey30"
    , "Chlorophyta-Dinophyceae"            = "grey30", "Cryptophyceae-Dinophyceae"      = "grey30"   , "Cyanophyceae-Dinophyceae"          = "grey30"
    , "Chlorophyta-Euglenozoa"             = "grey30", "Cryptophyceae-Euglenozoa"       = "grey30"   , "Chrysophyceae-Euglenozoa"          = "grey30"
    , "Cyanophyceae-Euglenozoa"            = "grey30", "Chlorophyta-Dictyochophyceae"   = "grey30"   , "Dictyochophyceae-Dinophyceae"      = "grey30"
    , "Cryptophyceae-Dictyochophyceae"     = "grey30", "Dictyochophyceae-Haptophyta"    = "grey30"   , "Dictyochophyceae-Euglenozoa"       = "grey30"
    , "Chrysophyceae-Dictyochophyceae"     = "grey30", "Cyanophyceae-Dictyochophyceae"  = "grey30"   , "Chlorophyta-Cryptophyceae"         = "grey30"
    , "Chlorophyta-Haptophyta"             = "grey30", "Chlorophyta-Chrysophyceae"      = "grey30"   , "Chlorophyta-Cyanophyceae"          = "grey30"
    , "Cryptophyceae-Haptophyta"           = "grey30", "Chrysophyceae-Cryptophyceae"    = "grey30"   , "Cryptophyceae-Cyanophyceae"        = "grey30"
    , "Euglenozoa-Haptophyta"              = "grey30", "Chrysophyceae-Haptophyta"       = "grey30"   , "Ebriophyceae-Haptophyta"           = "grey30"
    , "Cyanophyceae-Haptophyta"            = "grey30", "Ebriophyceae-Euglenozoa"        = "grey30"   , "Euglenozoa-Euglenozoa"             = "grey30"
    , "Chrysophyceae-Cyanophyceae"         = "grey30", "Ebriophyceae-Ebriophyceae" = "grey30"
  ))
atlantic_global
#ggsave('Atlantic_global_network.png', path = "output/graphs/networks/global_networks/", dpi = 600, width = 400, height = 300, units = 'mm')

# All networks in one plot
global_networks <- plot_grid(med_global,manche_global,atlantic_global,nrow = 1)
#ggsave("global_networks.png", path = "output/graphs/networks/global_networks", dpi = 600, width = 400, height = 900, units = 'mm')

### Compare composition of each global networks ####
# Import datasets
# Mediterranean sea
med <- read_delim("output/tableaux/Networks/taxon_list_med.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)
med$region <- "1-Mediterranean sea"

med$Number <- 1

# English channel
manche <- read_delim("output/tableaux/Networks/taxon_list_manche.csv", 
                  delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                      grouping_mark = ""), trim_ws = TRUE)
manche$region <- "2-Eastern Channel - North Sea"

manche$Number <- 1

# Atlantic
atlantic <- read_delim("output/tableaux/Networks/taxon_list_atlantic.csv", 
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
#write.csv2(tableau_compo,file="output/tableaux/Networks/compo_global_networks.csv", row.names = FALSE,dec = ".")


# Same but for the associations types 
# Import the table that contains the association frequencies 
# For the mediterranean sea
med <- read_delim("output/tableaux/Networks/freq_table_med.csv", 
                  delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                      grouping_mark = ""), trim_ws = TRUE)
# English channel
manche <- read_delim("output/tableaux/Networks/freq_table_manche.csv", 
                  delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                      grouping_mark = ""), trim_ws = TRUE)
# Atlantic
atlantic <- read_delim("output/tableaux/Networks/freq_table_atlantic.csv", 
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
freq_table_big2 <- summarise(group_by(freq_table_big,region,link_genus), Number = sum(Number)) |>
  mutate(Freq = Number / sum(Number))
#write.csv2(freq_table_big2,file="output/tableaux/Networks/asso_global_networks.csv", row.names = FALSE,dec = ".")


# Making a new graph to synthetize association and composition
# Import the dataset about composition
tableau_compo <- read_delim("output/tableaux/Networks/compo_global_networks.csv", 
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
freq_table_big <- read_delim("output/tableaux/Networks/asso_global_networks.csv", 
                             delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                                 grouping_mark = ""), trim_ws = TRUE)
# Create the dataframe as needed to make the graph
freq_table_big <- freq_table_big %>%
  group_by(region) %>%
  arrange(desc(link_genus)) %>%
  mutate(ymin = cumsum(Freq) - Freq,
         ymax = cumsum(Freq),
         labelPosition = (ymax + ymin) / 2,
         label = paste0(round(Freq * 100, 1), "%"))

# Making the graph
global_compo_asso <- ggplot() +
  geom_rect(aes(ymax = ymax, ymin = ymin, xmax = 3.9, xmin = 3, fill = Phylum.Classe),data = tableau_compo) +
  geom_text(aes(x = 4.5, y = labelPosition, label = label, colour = Phylum.Classe), size = 3,data = tableau_compo) +
  geom_rect(aes(ymax = ymax, ymin = ymin, xmax = 2.9, xmin = 2, fill = link_genus),data = freq_table_big) +
  geom_text(aes(x = 1.4, y = labelPosition, label = label, colour = link_genus), size = 3,data = freq_table_big) +
  scale_fill_manual(values = c(
    "Bacillariophyceae" = "#1f77b4", "Dinophyceae" = "green",
    "Bacillariophyceae-Bacillariophyceae" = "#56B4E9","Bacillariophyceae-Dinophyceae"   = "chocolate1"
    ,"Other phyla"      = "grey" ,"Dinophyceae-Dinophyceae" = "#009E73","Other association" = "grey30"
  ), breaks = c(
    "Bacillariophyceae","Dinophyceae","Other phyla", "Bacillariophyceae-Bacillariophyceae", "Bacillariophyceae-Dinophyceae", 
    "Dinophyceae-Dinophyceae", "Other association"
  ),name = NULL) +
  scale_colour_manual(values = c(
    "Bacillariophyceae" = "#1f77b4", "Dinophyceae" = "green",
    "Bacillariophyceae-Bacillariophyceae" = "#56B4E9","Bacillariophyceae-Dinophyceae"   = "chocolate1"
    ,"Other phyla"      = "grey" ,"Dinophyceae-Dinophyceae" = "#009E73","Other association" = "grey30"
  ), breaks = c(
    "Bacillariophyceae","Dinophyceae","Other phyla", "Bacillariophyceae-Bacillariophyceae", "Bacillariophyceae-Dinophyceae", 
    "Dinophyceae-Dinophyceae", "Other association"
  ),name = NULL)+
  coord_polar(theta = "y") +
  xlim(c(-1, 4.6)) +
  theme_no_axes() +
  facet_wrap(~region) +
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(ncol = 7), colour = guide_legend(ncol = 7))
global_compo_asso
#ggsave('asso_compo_regions_combine.png', path = "output/graphs/description_region/", dpi = 600, width = 300, height = 150, units = 'mm')

### PCA on the networks metrics #####
# Import data
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_final.csv", 
                       delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                           grouping_mark = ""), trim_ws = TRUE)
# Select the variables
data_pca <- select(data,N_noeuds:N_clust)
colnames(data_pca) <- c("N_nodes","N_edges","Dens","C_tance","Avg_p_length","Adhes","Mod","Avg_degree","Assort","Diss","Trans",
                        "Avg_edge_bet","Nat_connect","N_clust")
# Making the PCA
PCA_results <- PCA(data_pca)
fviz_screeplot(PCA_results) # Screeplot
# To check for a spatial dissimilarity we check by see the individuals position by region
fviz_pca_ind(PCA_results,col.ind = data$region,addEllipses = T,alpha.ind = 1,geom = c("point"))
# Axes 2 and 3
fviz_pca_var(PCA_results,axes = c(3,2))
# Store the individuals coordinates
PCA_coord <- PCA_results$ind$coord

fviz_pca_var(PCA_results,axes = c(1,2),repel = T,title="")
#ggsave('PCA_1.png', path = "output/graphs/networks/", dpi = 600, width = 100, height = 100, units = 'mm')
fviz_pca_var(PCA_results,axes = c(2,3),repel = T,title="")
#ggsave('PCA_2.png', path = "output/graphs/networks/", dpi = 600, width = 100, height = 100, units = 'mm')

PCA_coord <- as.data.frame(PCA_coord)
PCA_coord <- select(PCA_coord, Dim.1,Dim.2,Dim.3)
# Save it
#write.csv2(PCA_coord,file="output/tableaux/Networks/PCA_coords.csv", row.names = FALSE,dec = ".")


### Correlate PCA coordinates and biodiv indexes ####
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_PCA_coord_moments.final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Select the variables 
data_cor <- select(data,Dim.1:Dim.3,Shannon,Pielou,BergerParker)

# Consider only the date which have a value for all variables
Table.corr_all.comp <- data_cor[complete.cases(data_cor),]
correlations <- cor(Table.corr_all.comp,method = "spearman")

# Function to make a correlation plot with p-value
cor.mtest <- function(Table.corr_all.comp, ...) {
  mat <- as.matrix(Table.corr_all.comp)
  n <- ncol(Table.corr_all.comp)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Calculate and doing the correlation plot
p.mat <- cor.mtest(Table.corr_all.comp)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(correlations, method="color", col=col(200),  
         type="upper", order="original",tl.cex = 0.6,number.cex = 0.6,
         addCoef.col = "black",
         tl.col="black", tl.srt=45,
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         diag=F, 
         title = ""
)

