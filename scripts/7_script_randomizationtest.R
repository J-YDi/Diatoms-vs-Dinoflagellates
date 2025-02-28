################################################################################
# Diatoms vs dinoflagellates: a network analysis of bloom impacts on diversity #
#                    and phytoplankton associations | R scripts                #
################################################################################

# Script to perform randomization test #
# 02/28/2025

# Load packages  
library(readr)
library(ggplot2)
library(igraph)
library(tidyr)
library(dplyr)

# Randomization test #####


# Function to create null matrices 
generate_null_matrix <- function(mat) {
  n <- nrow(mat)  # Matrix size
  # Extract the indices from the upper triangular half (without diagonal)  upper_indices <- which(upper.tri(mat))
  upper_indices <- which(upper.tri(mat))
  
  # Extract values from the upper triangular half
  upper_values <- mat[upper_indices]
  
  # Swap these values randomly
  permuted_values <- sample(upper_values)
  
  # Create a new null matrix
  null_mat <- matrix(0, n, n)
  
  # Place the swapped values in the upper triangular half
  null_mat[upper_indices] <- permuted_values
  
  # Complete the matrix so that it is symmetrical
  null_mat <- null_mat + t(null_mat)
  
  # Return the null matrix
  return(null_mat)
}

# Function to normalize the association pairs as it Bac-Dino = Dino-Bac
normaliser_paire <- function(x) {
  parties <- sort(unlist(strsplit(x, "-")))
  paste(parties, collapse = "-")
}

#### Mediterranean Sea #####
# Load data
assoMat_med <- read_delim("output/tableaux/Networks/assoMat_med.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Load table "Liste_phylum.classe" to associate species to a class
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt",header = T,sep = "\t")

# Store the columns and lines names
column <- colnames(assoMat_med)
linesname <- rownames(assoMat_med)

# Creating a df to store the results
data_null_reseaux <- c("","")
data_null_reseaux <- as.data.frame(data_null_reseaux)

# Create null matrices
for (i in 1:10000){
# Create a null matrix
null_mat <- generate_null_matrix(as.matrix(assoMat_med))

null_mat <- as.data.frame(null_mat)
# associate the taxonomic info
colnames(null_mat) <- column
rownames(null_mat) <- linesname

# Create null networks
null_net <- graph_from_adjacency_matrix(as.matrix(null_mat), weighted = TRUE, mode = "undirected", diag = FALSE)
# Indicate the associations
edge_list_null <- igraph::as_data_frame(null_net, what = "edges")

# Associate an association type
edge1 <- as.data.frame(edge_list_null$from)
colnames(edge1) <- "Taxon"
edge1 <- left_join(edge1,phylum_classe)
colnames(edge1) <- c("v1","v1_classe")

edge2 <- as.data.frame(edge_list_null$to)
colnames(edge2) <- "Taxon"
edge2 <- left_join(edge2,phylum_classe)
colnames(edge2) <- c("v2","v2_classe")

edgelist_genus_null <- cbind(edge1,edge2,edge_list_null$weight)
edgelist_genus_null$link_genus <- paste0(edgelist_genus_null$v1_classe,"-",edgelist_genus_null$v2_classe) 

edgelist_genus_null$link_genus <- sapply(edgelist_genus_null$link_genus, normaliser_paire)

copy_edge1 <- edge1
colnames(copy_edge1) <- c("Genus","Classe")
copy_edge2 <- edge2
colnames(copy_edge2) <- c("Genus","Classe")
taxons <- rbind(copy_edge1,copy_edge2)
taxons <- unique(taxons)

# Store the information about composition and associations in the null networks

data_null_reseaux[i,1] <- nrow(filter(edgelist_genus_null, link_genus == "Bacillariophyceae-Bacillariophyceae")) # Number of Bac-Bac
data_null_reseaux[i,2] <- nrow(filter(edgelist_genus_null, link_genus == "Bacillariophyceae-Bacillariophyceae"))/nrow(edgelist_genus_null) # Pourcentage of Bac-Bac
data_null_reseaux[i,3] <- nrow(filter(edgelist_genus_null, link_genus == "Bacillariophyceae-Dinophyceae")) # Number of Bac-Dino
data_null_reseaux[i,4] <- nrow(filter(edgelist_genus_null, link_genus == "Bacillariophyceae-Dinophyceae"))/nrow(edgelist_genus_null) # Pourcentage of Bac-Dino
data_null_reseaux[i,5] <- nrow(filter(edgelist_genus_null, link_genus == "Dinophyceae-Dinophyceae")) # Number of Dino-Dino
data_null_reseaux[i,6] <- nrow(filter(edgelist_genus_null, link_genus == "Dinophyceae-Dinophyceae"))/nrow(edgelist_genus_null) # Pourcentage of Dino-Dino
data_null_reseaux[i,7] <- nrow(filter(edgelist_genus_null, link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" )) # Number of other associations
data_null_reseaux[i,8] <- nrow(filter(edgelist_genus_null,link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" ))/nrow(edgelist_genus_null) # Pourcentage of other associations
data_null_reseaux[i,9] <- nrow(edgelist_genus_null) # Number of associations

data_null_reseaux[i,10] <- nrow(filter(taxons,Classe == "Bacillariophyceae")) # Number of Bacillariophyceae
data_null_reseaux[i,11] <- nrow(filter(taxons,Classe == "Bacillariophyceae"))/nrow(taxons) # Pourcentage of Bacillariophyceae
data_null_reseaux[i,12] <- nrow(filter(taxons,Classe == "Dinophyceae")) # Number of Dinophyceae
data_null_reseaux[i,13] <- nrow(filter(taxons,Classe == "Dinophyceae"))/nrow(taxons) # Pourcentage of Dinophyceae
data_null_reseaux[i,14] <- nrow(filter(taxons,Classe != "Bacillariophyceae" & Classe != "Dinophyceae")) # Number of other taxons
data_null_reseaux[i,15] <- nrow(filter(taxons,Classe != "Bacillariophyceae"& Classe != "Dinophyceae"))/nrow(taxons) # Pourcentage of other taxons
data_null_reseaux[i,16] <- nrow(taxons) # "Species" richness


colnames(data_null_reseaux) <- c("N_BacBac","P_BacBac","N_BacDino","P_BacDino","N_DinoDino","P_DinoDino","N_AAutres",
                                    "P_AAutres","N_asso","N_bac","P_bac","N_dino","P_dino","N_autres","P_autres","N_taxons")
#print(i*100/10000)
}
# Save it
#write.csv2(data_null_reseaux,file="output/tableaux/Networks/results_null_med.csv", row.names = FALSE,dec = ".")

## Compare the observed values to the null networks
# Import data of the observed values
freq_asso <- read_delim("output/tableaux/Networks/asso_global_networks.csv", 
                       delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                           grouping_mark = ""), trim_ws = TRUE)
# Import data about the mediterranean null networks
null_med <- read_delim("output/tableaux/Networks/results_null_med.csv", 
                  delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                      grouping_mark = ""), trim_ws = TRUE)
# Calculate p-values
distribution <- null_med$P_BacDino 
valeur <- freq_asso$Freq[freq_asso$region == "1-Mediterranean sea" & freq_asso$link_genus == "Bacillariophyceae-Dinophyceae"]
t_obs <- (valeur - mean(distribution)) / sd(distribution)
df <- length(distribution) - 1  # df

p_value_bacdino <- 2 * (1 - pt(abs(t_obs), df))


ggplot()+
  geom_histogram(aes(distribution, y=..density..),bins = 20,fill = "chocolate1",colour="black")+
  geom_density(aes(distribution),alpha=.2, fill="chocolate1")+
  geom_vline(aes(xintercept = valeur),colour="red",size=2)+
  labs(title = "Bac_Dino Med",subtitle = p_value_bacdino)

distribution <- null_med$P_BacBac 
valeur <- filter(freq_asso,region == "1-Mediterranean sea" & link_genus == "Bacillariophyceae-Bacillariophyceae")$Freq
t_obs <- (valeur - mean(distribution)) / sd(distribution)
df <- length(distribution) - 1  # df

p_value_bacbac <- 2 * (1 - pt(abs(t_obs), df))

ggplot()+
  geom_histogram(aes(distribution, y=..density..),bins = 20,fill = "#56B4E9",colour="black")+
  geom_vline(aes(xintercept = valeur),colour="red",size=2)+
  geom_density(aes(distribution),alpha=.2, fill="#56B4E9")+
  labs(title = "Bac_Bac Med",subtitle = p_value_bacbac)

distribution <- null_med$P_DinoDino
valeur <- filter(freq_asso,region == "1-Mediterranean sea" & link_genus == "Dinophyceae-Dinophyceae")$Freq
t_obs <- (valeur - mean(distribution)) / sd(distribution)
df <- length(distribution) - 1  # df

p_value_dinodino <- 2 * (1 - pt(abs(t_obs), df))
ggplot()+
  geom_histogram(aes(distribution, y=..density..),bins = 20,fill = "#009E73",colour="black")+
  geom_density(aes(distribution),alpha=.2, fill="#009E73")+
  geom_vline(aes(xintercept = valeur),colour="red",size=2)+
  labs(title = "Dino_Dino Med",subtitle = p_value_dinodino)

distribution <- null_med$P_AAutres
valeur <- filter(freq_asso,region == "1-Mediterranean sea" & link_genus == "Other association")$Freq
t_obs <- (valeur - mean(distribution)) / sd(distribution)
df <- length(distribution) - 1  # df

p_value_aautres <- 2 * (1 - pt(abs(t_obs), df))
ggplot()+
  geom_histogram(aes(x= distribution, y=..density..),bins = 20,fill = "grey30",colour="black")+
  geom_density(aes(distribution),alpha=.2, fill="grey30")+
  geom_vline(aes(xintercept = valeur),colour="red",size=2)+
  labs(title = "AAutres Med",subtitle = p_value_aautres)

# Fisher's method to get global p-value
p_values <- c(p_value_aautres, p_value_bacbac, p_value_bacdino, p_value_dinodino)

chi_squared_stat <- -2 * sum(log(p_values))
p_global <- 1 - pchisq(chi_squared_stat, df = 2 * length(p_values))

p_global




#### Eastern Channel - North Sea #####
# Load data
assoMat_manche <- read_delim("output/tableaux/Networks/assoMat_manche.csv", 
                          delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                              grouping_mark = ""), trim_ws = TRUE)
# Load table "Liste_phylum.classe" to associate species to a class
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt",header = T,sep = "\t")

# Store the columns and lines names
column <- colnames(assoMat_manche)
linesname <- rownames(assoMat_manche)

# Creating a df to store the results
data_null_reseaux <- c("","")
data_null_reseaux <- as.data.frame(data_null_reseaux)

# Create null matrices
for (i in 1:10000){
  # Create a null matrix
  null_mat <- generate_null_matrix(as.matrix(assoMat_manche))
  
  null_mat <- as.data.frame(null_mat)
  # associate the taxonomic info
  colnames(null_mat) <- column
  rownames(null_mat) <- linesname
  
  # Create null networks
  null_net <- graph_from_adjacency_matrix(as.matrix(null_mat), weighted = TRUE, mode = "undirected", diag = FALSE)
  # Indicate the associations
  edge_list_null <- igraph::as_data_frame(null_net, what = "edges")
  
  # Associate an association type
  edge1 <- as.data.frame(edge_list_null$from)
  colnames(edge1) <- "Taxon"
  edge1 <- left_join(edge1,phylum_classe)
  colnames(edge1) <- c("v1","v1_classe")
  
  edge2 <- as.data.frame(edge_list_null$to)
  colnames(edge2) <- "Taxon"
  edge2 <- left_join(edge2,phylum_classe)
  colnames(edge2) <- c("v2","v2_classe")
  
  edgelist_genus_null <- cbind(edge1,edge2,edge_list_null$weight)
  edgelist_genus_null$link_genus <- paste0(edgelist_genus_null$v1_classe,"-",edgelist_genus_null$v2_classe) 
  
  edgelist_genus_null$link_genus <- sapply(edgelist_genus_null$link_genus, normaliser_paire)
  
  copy_edge1 <- edge1
  colnames(copy_edge1) <- c("Genus","Classe")
  copy_edge2 <- edge2
  colnames(copy_edge2) <- c("Genus","Classe")
  taxons <- rbind(copy_edge1,copy_edge2)
  taxons <- unique(taxons)
  
  # Store the information about composition and associations in the null networks
  data_null_reseaux[i,1] <- nrow(filter(edgelist_genus_null, link_genus == "Bacillariophyceae-Bacillariophyceae")) # Number of Bac-Bac
  data_null_reseaux[i,2] <- nrow(filter(edgelist_genus_null, link_genus == "Bacillariophyceae-Bacillariophyceae"))/nrow(edgelist_genus_null) # Pourcentage of Bac-Bac
  data_null_reseaux[i,3] <- nrow(filter(edgelist_genus_null, link_genus == "Bacillariophyceae-Dinophyceae")) # Number of Bac-Dino
  data_null_reseaux[i,4] <- nrow(filter(edgelist_genus_null, link_genus == "Bacillariophyceae-Dinophyceae"))/nrow(edgelist_genus_null) # Pourcentage of Bac-Dino
  data_null_reseaux[i,5] <- nrow(filter(edgelist_genus_null, link_genus == "Dinophyceae-Dinophyceae")) # Number of Dino-Dino
  data_null_reseaux[i,6] <- nrow(filter(edgelist_genus_null, link_genus == "Dinophyceae-Dinophyceae"))/nrow(edgelist_genus_null) # Pourcentage of Dino-Dino
  data_null_reseaux[i,7] <- nrow(filter(edgelist_genus_null, link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" )) #  Number of other associations
  data_null_reseaux[i,8] <- nrow(filter(edgelist_genus_null,link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" ))/nrow(edgelist_genus_null) # Pourcentage of other associations
  data_null_reseaux[i,9] <- nrow(edgelist_genus_null) # Number of associations
  
  data_null_reseaux[i,10] <- nrow(filter(taxons,Classe == "Bacillariophyceae")) # Number of Bacillariophyceae
  data_null_reseaux[i,11] <- nrow(filter(taxons,Classe == "Bacillariophyceae"))/nrow(taxons) # Pourcentage of Bacillariophyceae
  data_null_reseaux[i,12] <- nrow(filter(taxons,Classe == "Dinophyceae")) # Number of Dinophyceae
  data_null_reseaux[i,13] <- nrow(filter(taxons,Classe == "Dinophyceae"))/nrow(taxons) # Pourcentage of Dinophyceae
  data_null_reseaux[i,14] <- nrow(filter(taxons,Classe != "Bacillariophyceae" & Classe != "Dinophyceae")) # Number of other taxa
  data_null_reseaux[i,15] <- nrow(filter(taxons,Classe != "Bacillariophyceae"& Classe != "Dinophyceae"))/nrow(taxons) # Pourcentage of other taxa
  data_null_reseaux[i,16] <- nrow(taxons) # "Species" richness
  
  
  colnames(data_null_reseaux) <- c("N_BacBac","P_BacBac","N_BacDino","P_BacDino","N_DinoDino","P_DinoDino","N_AAutres",
                                   "P_AAutres","N_asso","N_bac","P_bac","N_dino","P_dino","N_autres","P_autres","N_taxons")
  #print(i*100/10000)
}
# Save it
# write.csv2(data_null_reseaux,file="output/tableaux/Networks/results_null_manche.csv", row.names = FALSE,dec = ".")

## Compare the observed values to the null networks
# Import data of the observed values
freq_asso <- read_delim("output/tableaux/Networks/asso_global_networks.csv", 
                        delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                            grouping_mark = ""), trim_ws = TRUE)
# Import data about the eastern channel null networks
null_manche <- read_delim("output/tableaux/Networks/results_null_manche.csv", 
                       delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                           grouping_mark = ""), trim_ws = TRUE)
# Calculate p-values
distribution <- null_manche$P_BacDino 
valeur <- filter(freq_asso,region == "2-Eastern Channel - North Sea" & link_genus == "Bacillariophyceae-Dinophyceae")$Freq
t_obs <- (valeur - mean(distribution)) / sd(distribution)
df <- length(distribution) - 1  # df

p_value_bacdino <- 2 * (1 - pt(abs(t_obs), df))

ggplot()+
  geom_histogram(aes(distribution, y=..density..),bins = 20,fill = "chocolate1",colour="black")+
  geom_density(aes(distribution),alpha=.2, fill="chocolate1")+
  geom_vline(aes(xintercept = valeur),colour="red",size=2)+
  labs(title = "Bac_Dino manche",subtitle = p_value_bacdino)

distribution <- null_manche$P_BacBac 
valeur <- filter(freq_asso,region == "2-Eastern Channel - North Sea" & link_genus == "Bacillariophyceae-Bacillariophyceae")$Freq
t_obs <- (valeur - mean(distribution)) / sd(distribution)
df <- length(distribution) - 1  # df


p_value_bacbac <- 2 * (1 - pt(abs(t_obs), df))
ggplot()+
  geom_histogram(aes(distribution, y=..density..),bins = 20,fill = "#56B4E9",colour="black")+
  geom_vline(aes(xintercept = valeur),colour="red",size=2)+
  geom_density(aes(distribution),alpha=.2, fill="#56B4E9")+
  labs(title = "Bac_Bac manche",subtitle = p_value_bacbac)

distribution <- null_manche$P_DinoDino
valeur <- filter(freq_asso,region == "2-Eastern Channel - North Sea" & link_genus == "Dinophyceae-Dinophyceae")$Freq
t_obs <- (valeur - mean(distribution)) / sd(distribution)
df <- length(distribution) - 1  # df

p_value_dinodino <- 2 * (1 - pt(abs(t_obs), df))
ggplot()+
  geom_histogram(aes(distribution, y=..density..),bins = 20,fill = "#009E73",colour="black")+
  geom_density(aes(distribution),alpha=.2, fill="#009E73")+
  geom_vline(aes(xintercept = valeur),colour="red",size=2)+
  labs(title = "Dino_Dino manche",subtitle = p_value_dinodino)

distribution <- null_manche$P_AAutres
valeur <- filter(freq_asso,region == "2-Eastern Channel - North Sea" & link_genus == "Other association")$Freq
t_obs <- (valeur - mean(distribution)) / sd(distribution)
df <- length(distribution) - 1  # df

p_value_aautres <- 2 * (1 - pt(abs(t_obs), df))
ggplot()+
  geom_histogram(aes(x= distribution, y=..density..),bins = 20,fill = "grey30",colour="black")+
  geom_density(aes(distribution),alpha=.2, fill="grey30")+
  geom_vline(aes(xintercept = valeur),colour="red",size=2)+
  labs(title = "AAutres manche",subtitle = p_value_aautres)


# Fisher's method to get global p-value
p_values <- c(p_value_aautres, p_value_bacbac, p_value_bacdino, p_value_dinodino)
chi_squared_stat <- -2 * sum(log(p_values))

p_global <- 1 - pchisq(chi_squared_stat, df = 2 * length(p_values))

p_global



#### Atlantic - Western Channel #####
# Load data
assoMat_atlantique <- read_delim("output/tableaux/Networks/assoMat_atlantic.csv", 
                             delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                                 grouping_mark = ""), trim_ws = TRUE)
# Load table "Liste_phylum.classe" to associate species to a class
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt",header = T,sep = "\t")

# Store the columns and lines names
column <- colnames(assoMat_atlantique)
linesname <- rownames(assoMat_atlantique)

# Creating a df to store the results
data_null_reseaux <- c("","")
data_null_reseaux <- as.data.frame(data_null_reseaux)

# Create null matrices
for (i in 1:10000){
  # Create a null matrix
  null_mat <- generate_null_matrix(as.matrix(assoMat_atlantique))
  
  null_mat <- as.data.frame(null_mat)
  # associate the taxonomic info
  colnames(null_mat) <- column
  rownames(null_mat) <- linesname
  
  # Create null networks
  null_net <- graph_from_adjacency_matrix(as.matrix(null_mat), weighted = TRUE, mode = "undirected", diag = FALSE)
  # Indicate the associations
  edge_list_null <- igraph::as_data_frame(null_net, what = "edges")
  
  # Associate an association type
  edge1 <- as.data.frame(edge_list_null$from)
  colnames(edge1) <- "Taxon"
  edge1 <- left_join(edge1,phylum_classe)
  colnames(edge1) <- c("v1","v1_classe")
  
  edge2 <- as.data.frame(edge_list_null$to)
  colnames(edge2) <- "Taxon"
  edge2 <- left_join(edge2,phylum_classe)
  colnames(edge2) <- c("v2","v2_classe")
  
  edgelist_genus_null <- cbind(edge1,edge2,edge_list_null$weight)
  edgelist_genus_null$link_genus <- paste0(edgelist_genus_null$v1_classe,"-",edgelist_genus_null$v2_classe) 
  
  edgelist_genus_null$link_genus <- sapply(edgelist_genus_null$link_genus, normaliser_paire)
  
  copy_edge1 <- edge1
  colnames(copy_edge1) <- c("Genus","Classe")
  copy_edge2 <- edge2
  colnames(copy_edge2) <- c("Genus","Classe")
  taxons <- rbind(copy_edge1,copy_edge2)
  taxons <- unique(taxons)
  
  # Store the information about composition and associations in the null networks
  
  data_null_reseaux[i,1] <- nrow(filter(edgelist_genus_null, link_genus == "Bacillariophyceae-Bacillariophyceae")) # Number of Bac-Bac
  data_null_reseaux[i,2] <- nrow(filter(edgelist_genus_null, link_genus == "Bacillariophyceae-Bacillariophyceae"))/nrow(edgelist_genus_null) # Pourcentage of Bac-Bac
  data_null_reseaux[i,3] <- nrow(filter(edgelist_genus_null, link_genus == "Bacillariophyceae-Dinophyceae")) # Number of Bac-Dino
  data_null_reseaux[i,4] <- nrow(filter(edgelist_genus_null, link_genus == "Bacillariophyceae-Dinophyceae"))/nrow(edgelist_genus_null) # Pourcentage of Bac-Dino
  data_null_reseaux[i,5] <- nrow(filter(edgelist_genus_null, link_genus == "Dinophyceae-Dinophyceae")) # Number of Dino-Dino
  data_null_reseaux[i,6] <- nrow(filter(edgelist_genus_null, link_genus == "Dinophyceae-Dinophyceae"))/nrow(edgelist_genus_null) # Pourcentage of Dino-Dino
  data_null_reseaux[i,7] <- nrow(filter(edgelist_genus_null, link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" )) # Number of other associations
  data_null_reseaux[i,8] <- nrow(filter(edgelist_genus_null,link_genus != "Dinophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Dinophyceae" & link_genus != "Bacillariophyceae-Bacillariophyceae" ))/nrow(edgelist_genus_null) # Pourcentage of other associations
  data_null_reseaux[i,9] <- nrow(edgelist_genus_null) # Number of associations
  
  data_null_reseaux[i,10] <- nrow(filter(taxons,Classe == "Bacillariophyceae")) # Number of Bacillariophyceae
  data_null_reseaux[i,11] <- nrow(filter(taxons,Classe == "Bacillariophyceae"))/nrow(taxons) # Pourcentage of Bacillariophyceae
  data_null_reseaux[i,12] <- nrow(filter(taxons,Classe == "Dinophyceae")) # Number of Dinophyceae
  data_null_reseaux[i,13] <- nrow(filter(taxons,Classe == "Dinophyceae"))/nrow(taxons) # Pourcentage of Dinophyceae
  data_null_reseaux[i,14] <- nrow(filter(taxons,Classe != "Bacillariophyceae" & Classe != "Dinophyceae")) # Number of other taxons
  data_null_reseaux[i,15] <- nrow(filter(taxons,Classe != "Bacillariophyceae"& Classe != "Dinophyceae"))/nrow(taxons) # Pourcentage of other taxons
  data_null_reseaux[i,16] <- nrow(taxons) # "Species richness"
  
  
  colnames(data_null_reseaux) <- c("N_BacBac","P_BacBac","N_BacDino","P_BacDino","N_DinoDino","P_DinoDino","N_AAutres",
                                   "P_AAutres","N_asso","N_bac","P_bac","N_dino","P_dino","N_autres","P_autres","N_taxons")
  #print(i*100/1000)
}
# Save it
# write.csv2(data_null_reseaux,file="output/tableaux/Networks/results_null_atlantique.csv", row.names = FALSE,dec = ".")

## Compare the observed values to the null networks
# Import data of the observed values
freq_asso <- read_delim("output/tableaux/Networks/asso_global_networks.csv", 
                        delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                            grouping_mark = ""), trim_ws = TRUE)
# Import data about the atlantic null networks
null_atlantique <- read_delim("output/tableaux/Networks/results_null_atlantique.csv", 
                          delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                              grouping_mark = ""), trim_ws = TRUE)
# Calculate p-values
distribution <- null_atlantique$P_BacDino 
valeur <- filter(freq_asso,region == "3-Atlantic - Western Channel" & link_genus == "Bacillariophyceae-Dinophyceae")$Freq
t_obs <- (valeur - mean(distribution)) / sd(distribution)
df <- length(distribution) - 1  # df

p_value_bacdino <- 2 * (1 - pt(abs(t_obs), df))
ggplot()+
  geom_histogram(aes(distribution, y=..density..),bins = 20,fill = "chocolate1",colour="black")+
  geom_density(aes(distribution),alpha=.2, fill="chocolate1")+
  geom_vline(aes(xintercept = valeur),colour="red",size=2)+
  labs(title = "Bac_Dino atlantique",subtitle =p_value_bacdino )

distribution <- null_atlantique$P_BacBac 
valeur <- filter(freq_asso,region == "3-Atlantic - Western Channel" & link_genus == "Bacillariophyceae-Bacillariophyceae")$Freq
t_obs <- (valeur - mean(distribution)) / sd(distribution)
df <- length(distribution) - 1  # df

p_value_bacbac <- 2 * (1 - pt(abs(t_obs), df))
ggplot()+
  geom_histogram(aes(distribution, y=..density..),bins = 20,fill = "#56B4E9",colour="black")+
  geom_vline(aes(xintercept = valeur),colour="red",size=2)+
  geom_density(aes(distribution),alpha=.2, fill="#56B4E9")+
  labs(title = "Bac_Bac atlantique",subtitle =p_value_bacbac)

distribution <- null_atlantique$P_DinoDino
valeur <- filter(freq_asso,region == "3-Atlantic - Western Channel" & link_genus == "Dinophyceae-Dinophyceae")$Freq
t_obs <- (valeur - mean(distribution)) / sd(distribution)
df <- length(distribution) - 1  # df

p_value_dinodino <- 2 * (1 - pt(abs(t_obs), df))
ggplot()+
  geom_histogram(aes(distribution, y=..density..),bins = 20,fill = "#009E73",colour="black")+
  geom_density(aes(distribution),alpha=.2, fill="#009E73")+
  geom_vline(aes(xintercept = valeur),colour="red",size=2)+
  labs(title = "Dino_Dino atlantique",subtitle = p_value_dinodino)

distribution <- null_atlantique$P_AAutres
valeur <- filter(freq_asso,region == "3-Atlantic - Western Channel" & link_genus == "Other association")$Freq
t_obs <- (valeur - mean(distribution)) / sd(distribution)
df <- length(distribution) - 1  # df


p_value_aautres <- 2 * (1 - pt(abs(t_obs), df))
ggplot()+
  geom_histogram(aes(x= distribution, y=..density..),bins = 20,fill = "grey30",colour="black")+
  geom_density(aes(distribution),alpha=.2, fill="grey30")+
  geom_vline(aes(xintercept = valeur),colour="red",size=2)+
  labs(title = "AAutres atlantique",subtitle = p_value_aautres)

# Fisher's method to get global p-value
p_values <- c(p_value_aautres, p_value_bacbac, p_value_bacdino, p_value_dinodino)
chi_squared_stat <- -2 * sum(log(p_values))

p_global <- 1 - pchisq(chi_squared_stat, df = 2 * length(p_values))

p_global
