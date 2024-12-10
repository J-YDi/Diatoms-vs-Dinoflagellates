# Script Suite Stage JY Dias # 10/12/2024

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
