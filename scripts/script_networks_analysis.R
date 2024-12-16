# Script Suite Stage JY Dias # 12/12/2024

library(ggalluvial)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)

### Working on the Mediterranean sea #####
# Import data
data <- read_delim("output/tableaux/Networks/edgelist_med.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt",header = T,sep = "\t")


edge1 <- as.data.frame(data$v1)
colnames(edge1) <- "Taxon"
edge1 <- left_join(edge1,phylum_classe)
colnames(edge1) <- c("v1","v1_classe")

edge2 <- as.data.frame(data$v2)
colnames(edge2) <- "Taxon"
edge2 <- left_join(edge2,phylum_classe)
colnames(edge2) <- c("v2","v2_classe")

edgelist_genus_med <- cbind(edge1,edge2,data$asso)
edgelist_genus_med$link_genus <- paste0(edgelist_genus_med$v1_classe,"-",edgelist_genus_med$v2_classe) 

normaliser_paire <- function(x) {
  parties <- sort(unlist(strsplit(x, "-")))
  paste(parties, collapse = "-")
}

edgelist_genus_med$link_genus <- sapply(edgelist_genus_med$link_genus, normaliser_paire)
write.csv2(edgelist_genus_med,file="output/tableaux/Networks/edgelist_genus_med.csv", row.names = FALSE,dec = ".")

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
write.csv2(freq_table_med,file="output/tableaux/Networks/freq_table_med.csv", row.names = FALSE,dec = ".")


# Create the alluvial graph
ggplot(freq_table_med, aes(axis1 = v1_classe, axis2 = v2_classe, y = Number)) +
  geom_alluvium(aes(fill = v1_classe), alpha = 0.6, curve_type = "quintic") + 
  geom_stratum(aes(fill = after_stat(stratum)),colour="white",width = 0.05) +  
  scale_x_discrete(expand = c(-0.1, 0.1)) +
  scale_fill_manual(
    values = c(
      "Bacillariophyceae" = "#1f77b4", "Chlorophyta" = "#ff7f0e", "Chrysophyceae" = "#2ca02c",
      "Cryptophyceae" = "#d62728", "Cyanophyceae" = "cyan", "Dictyochophyceae" = "#9467bd", 
      "Dinophyceae" = "green", "Ebriophyceae" = "yellow", "Euglenozoa" = "red","Haptophyta" = "maroon"
    )
  ) +
  labs(title = "Mediterranean sea global network associations",
       x = NULL,
       y = NULL, fill = "Phylum") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank(),  
    panel.grid = element_blank()  
  )
ggsave('Mediterranean_alluvial_global.png', path = "output/graphs/networks/global_networks", dpi = 600, width = 200, height = 200, units = 'mm')


ggplot(freq_table_med, aes(y = reorder(link_genus, Frequence), x = Frequence*100)) +
  geom_bar(stat = "identity", fill = "#F8766D") +  
  geom_text(aes(label = round(Frequence,digits = 3)*100), hjust = -0.2, vjust = 0.5, size = 3) +  
  labs(title = "Mediterranean sea global network associations",
       x = "Frequence (%)", 
       y = "Association") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
ggsave('Mediterranean_freq_global.png', path = "output/graphs/networks/global_networks", dpi = 600, width = 200, height = 200, units = 'mm')

sum(filter(freq_table_med,v1_classe=="Bacillariophyceae")$Frequence)
sum(filter(freq_table_med,v1_classe=="Dinophyceae")$Frequence)


# Charger la matrice d'adjacence
assoMat <- read_delim("output/tableaux/Networks/assoMat_med.csv", 
                      delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                          grouping_mark = ""), trim_ws = TRUE)

# Créer le graphe
cluster1 <- graph_from_adjacency_matrix(as.matrix(assoMat), weighted = TRUE, mode = "undirected", diag = FALSE)

# Clustering des nœuds
wcini <- cluster_fast_greedy(cluster1)

# Importer les informations taxonomiques
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt", header = TRUE, sep = "\t")

# Associer les clusters et taxonomies
compo_reseau <- data.frame(Taxon = wcini$names, cluster = wcini$membership) %>%
  left_join(phylum_classe, by = c("Taxon" = "Taxon")) %>%
  mutate(color = case_when(
    Phylum.Classe == "Bacillariophyceae" ~ "#56B4E9",
    Phylum.Classe == "Dinophyceae" ~ "#009E73",
    Phylum.Classe == "Ciliophora" ~ "#F0E442",
    Phylum.Classe == "Cryptophyceae" ~ "#CC79A7",
    Phylum.Classe == "Haptophyta" ~ "#996136",
    TRUE ~ "grey"
  ))


edgelist_genus_med <- read_delim("output/tableaux/Networks/edgelist_genus_med.csv", 
                      delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                          grouping_mark = ""), trim_ws = TRUE)
# Adding some customization info for ggraph
V(cluster1)$Phylum <- compo_reseau$Phylum.Classe
V(cluster1)$size <- 15
E(cluster1)$weight <- E(cluster1)$weight * 5  
V(cluster1)$label <- V(cluster1)$name
V(cluster1)$color <- compo_reseau$color
E(cluster1)$info <- edgelist_genus_med$link_genus


# making the global network plot
med_global <- ggraph(cluster1, layout = "auto") +  
  geom_edge_link(aes(edge_width = weight,edge_colour = info), alpha = 0.5, lineend = "round", linejoin = "bevel") +
  geom_node_point(aes(size = size, fill = color, stroke = 1), shape = 21,alpha = 0.5) +
  geom_node_text(aes(label = label, fontface = "bold"),colour = "black", 
                 repel = TRUE, size = 3, check_overlap = TRUE,
                 vjust = 0.5, hjust = 1) + 
  scale_edge_width(range = c(0.1, 2)) +
  scale_size_continuous(range = c(3, 10)) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c(
      "#56B4E9" = "#56B4E9",
      "#009E73" = "#009E73",
      "#F0E442" = "#F0E442",
      "#CC79A7" = "#CC79A7",
      "#996136" = "#996136",
      "grey" = "grey"
    )
  )+
  scale_edge_colour_manual(values = c(
    "Bacillariophyceae-Bacillariophyceae" = "#56B4E9","Bacillariophyceae-Dinophyceae"   = "chocolate1"    ,"Bacillariophyceae-Euglenozoa" = "grey"
   ,"Bacillariophyceae-Dictyochophyceae"  = "grey","Bacillariophyceae-Chlorophyta"   = "grey"   ,"Bacillariophyceae-Chrysophyceae"    = "grey"
   ,"Bacillariophyceae-Cyanophyceae"      = "grey","Bacillariophyceae-Haptophyta"    = "grey"   ,"Dinophyceae-Dinophyceae" = "#009E73"
   , "Dinophyceae-Euglenozoa"             = "grey", "Dinophyceae-Ebriophyceae"       = "grey"   , "Bacillariophyceae-Ebriophyceae"    = "grey"
   , "Bacillariophyceae-Cryptophyceae"    = "grey", "Dinophyceae-Haptophyta"         = "grey"   , "Chrysophyceae-Dinophyceae"         = "grey"
   , "Chlorophyta-Dinophyceae"            = "grey", "Cryptophyceae-Dinophyceae"      = "grey"   , "Cyanophyceae-Dinophyceae"          = "grey"
   , "Chlorophyta-Euglenozoa"             = "grey", "Cryptophyceae-Euglenozoa"       = "grey"   , "Chrysophyceae-Euglenozoa"          = "grey"
   , "Cyanophyceae-Euglenozoa"            = "grey", "Chlorophyta-Dictyochophyceae"   = "grey"   , "Dictyochophyceae-Dinophyceae"      = "grey"
   , "Cryptophyceae-Dictyochophyceae"     = "grey", "Dictyochophyceae-Haptophyta"    = "grey"   , "Dictyochophyceae-Euglenozoa"       = "grey"
   , "Chrysophyceae-Dictyochophyceae"     = "grey", "Cyanophyceae-Dictyochophyceae"  = "grey"   , "Chlorophyta-Cryptophyceae"         = "grey"
   , "Chlorophyta-Haptophyta"             = "grey", "Chlorophyta-Chrysophyceae"      = "grey"   , "Chlorophyta-Cyanophyceae"          = "grey"
   , "Cryptophyceae-Haptophyta"           = "grey", "Chrysophyceae-Cryptophyceae"    = "grey"   , "Cryptophyceae-Cyanophyceae"        = "grey"
   , "Euglenozoa-Haptophyta"              = "grey", "Chrysophyceae-Haptophyta"       = "grey"   , "Ebriophyceae-Haptophyta"           = "grey"
   , "Cyanophyceae-Haptophyta"            = "grey", "Ebriophyceae-Euglenozoa"        = "grey"   , "Euglenozoa-Euglenozoa"             = "grey"
   , "Chrysophyceae-Cyanophyceae"         = "grey", "Ebriophyceae-Ebriophyceae" = "grey"
  ))

ggsave('Med_global_network.png', path = "output/graphs/networks/global_networks/", dpi = 600, width = 400, height = 300, units = 'mm')


##### Type of association through time ######
data <- read_delim("output/tableaux/Networks/subnetworks/results_metrics_reseaux_cluster1.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

data_association <- pivot_longer(data,cols = c(P_BacBac,P_BacDino,P_DinoDino,P_AAutres),names_to = "Asso_type" )
data_association <- select(data_association,Code_point_Libelle,Date,Asso_type,value)

data_association$MonthYear <- format(data_association$Date, "%Y-%m")
data_association$Month <- format(data_association$Date, "%m")

data_association_graph <- summarise(group_by(data_association,Code_point_Libelle,MonthYear,Asso_type,Month),value=mean(value,na.rm=T))

ggplot(data_association_graph) +
  geom_col(aes(y = value, x = MonthYear, fill = Asso_type), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")+
  geom_vline(data = subset(data_association_graph,Month == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", size = 0.5,linetype = "dashed")+
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  scale_fill_manual(values = c("P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  geom_text(data = subset(data_association_graph, Month == "01"), 
            aes(x = MonthYear, y = 0.95, label = sub("-.*", "", MonthYear)),
            color = "black", size = 2.5,angle=90, vjust = 0)+
  facet_wrap(~Code_point_Libelle,nrow = 4)+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 5))
ggsave('med_typeassociation_station.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 400, height = 300, units = 'mm')

data_association_graph <- summarise(group_by(data_association,Code_point_Libelle,Asso_type,Month),value=mean(value,na.rm=T))

ggplot(data_association_graph) +
  geom_col(aes(y = value, x = Month, fill = Asso_type), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")+
  facet_wrap(~Code_point_Libelle,nrow = 4)
ggsave('med_typeassociation_station_mois.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 400, height = 300, units = 'mm')



data_compo <- pivot_longer(data,cols = c(P_bac,P_dino,P_autres),names_to = "Taxon" )
data_compo <- select(data_compo,Code_point_Libelle,Date,Taxon,value)

data_compo$MonthYear <- format(data_compo$Date, "%Y-%m")
data_compo$Month <- format(data_compo$Date, "%m")

data_compo_graph <- summarise(group_by(data_compo,Code_point_Libelle,MonthYear,Taxon,Month),value=mean(value,na.rm=T))

ggplot(data_compo_graph) +
  geom_col(aes(y = value, x = MonthYear, fill = Taxon), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")+
  geom_vline(data = subset(data_compo_graph,Month == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", size = 0.5,linetype = "dashed")+
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  scale_fill_manual(values = c("P_bac" = "#56B4E9","P_autres"      = "grey" ,"P_dino" = "#009E73"))+
  geom_text(data = subset(data_compo_graph, Month == "01"), 
            aes(x = MonthYear, y = 0.95, label = sub("-.*", "", MonthYear)),
            color = "black", size = 2.5,angle=90, vjust = 0)+
  facet_wrap(~Code_point_Libelle,nrow = 4)+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 5))
ggsave('med_compo_station.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 400, height = 300, units = 'mm')

data_compo_graph <- summarise(group_by(data_compo,Code_point_Libelle,Taxon,Month),value=mean(value,na.rm=T))

ggplot(data_compo_graph) +
  geom_col(aes(y = value, x = Month, fill = Taxon), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")+
  facet_wrap(~Code_point_Libelle,nrow = 4)
ggsave('med_compo_station_mois.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 400, height = 300, units = 'mm')


data_association_graph <- summarise(group_by(data_association,MonthYear,Asso_type,Month),value=mean(value,na.rm=T))

ggplot(data_association_graph) +
  geom_col(aes(y = value, x = MonthYear, fill = Asso_type), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")+
  geom_vline(data = subset(data_association_graph,Month == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", size = 0.5,linetype = "dashed")+
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  scale_fill_manual(values = c("P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  geom_text(data = subset(data_association_graph, Month == "01"), 
            aes(x = MonthYear, y = 0.95, label = sub("-.*", "", MonthYear)),
            color = "black", size = 2.5,angle=90, vjust = 0)+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 5))
ggsave('med_typeassociation.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 260, height = 170, units = 'mm')

data_association_graph <- summarise(group_by(data_association,Asso_type,Month),value=mean(value,na.rm=T))

ggplot(data_association_graph) +
  geom_col(aes(y = value, x = Month, fill = Asso_type), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")
ggsave('med_typeassociation_mois.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 400, height = 300, units = 'mm')



data_compo_graph <- summarise(group_by(data_compo,MonthYear,Taxon,Month),value=mean(value,na.rm=T))

ggplot(data_compo_graph) +
  geom_col(aes(y = value, x = MonthYear, fill = Taxon), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")+
  geom_vline(data = subset(data_compo_graph,Month == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", size = 0.5,linetype = "dashed")+
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  scale_fill_manual(values = c("P_bac" = "#56B4E9","P_autres"      = "grey" ,"P_dino" = "#009E73"))+
  geom_text(data = subset(data_compo_graph, Month == "01"), 
            aes(x = MonthYear, y = 0.95, label = sub("-.*", "", MonthYear)),
            color = "black", size = 2.5,angle=90, vjust = 0)+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 5))
ggsave('med_compo.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 260, height = 170, units = 'mm')

data_compo_graph <- summarise(group_by(data_compo,Taxon,Month),value=mean(value,na.rm=T))

ggplot(data_compo_graph) +
  geom_col(aes(y = value, x = Month, fill = Taxon), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")
ggsave('med_compo_mois.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 400, height = 300, units = 'mm')

data_pca1 <- select(data,N_bac,N_dino,N_autres,N_BacBac,N_BacDino,N_DinoDino,N_AAutres)
PCA(data_pca1)

data_pca2 <- select(data,N_noeuds:N_clust,P_BacBac,P_BacDino,P_DinoDino,P_AAutres)
PCA(data_pca2)

data_pca3 <- select(data,N_noeuds:N_clust)
PCA(data_pca3)

data_pca4 <- select(data,N_noeuds:N_clust,N_BacBac,N_BacDino,N_DinoDino,N_AAutres)
PCA(data_pca4)

Table.corr_all.comp <- data_pca2[complete.cases(data_pca2),]
correlations <- cor(Table.corr_all.comp,method = "spearman")

# ... : Arguments supplémentaire à passer à la fonction cor.test
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

# Matrice de p-value de la corrélation
p.mat <- cor.mtest(Table.corr_all.comp)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(correlations, method="color", col=col(200),  
         type="upper", order="alphabet",tl.cex = 0.7,
         addCoef.col = "black", # Ajout du coefficient de corrélation
         tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
         # Combiner avec le niveau de significativité
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # Cacher les coefficients de corrélation sur la diagonale
         diag=F, 
         title = ""
)

Table.corr_all.comp <- data_pca4[complete.cases(data_pca4),]
correlations <- cor(Table.corr_all.comp,method = "spearman")

# ... : Arguments supplémentaire à passer à la fonction cor.test
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

# Matrice de p-value de la corrélation
p.mat <- cor.mtest(Table.corr_all.comp)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(correlations, method="color", col=col(200),  
         type="upper", order="alphabet",tl.cex = 0.7,
         addCoef.col = "black", # Ajout du coefficient de corrélation
         tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
         # Combiner avec le niveau de significativité
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # Cacher les coefficients de corrélation sur la diagonale
         diag=F, 
         title = ""
)


### Working on the English channel #####
# Import data
data <- read_delim("output/tableaux/Networks/edgelist_manche.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt",header = T,sep = "\t")


edge1 <- as.data.frame(data$v1)
colnames(edge1) <- "Taxon"
edge1 <- left_join(edge1,phylum_classe)
colnames(edge1) <- c("v1","v1_classe")

edge2 <- as.data.frame(data$v2)
colnames(edge2) <- "Taxon"
edge2 <- left_join(edge2,phylum_classe)
colnames(edge2) <- c("v2","v2_classe")

edgelist_genus_manche <- cbind(edge1,edge2,data$asso)
edgelist_genus_manche$link_genus <- paste0(edgelist_genus_manche$v1_classe,"-",edgelist_genus_manche$v2_classe) 

normaliser_paire <- function(x) {
  parties <- sort(unlist(strsplit(x, "-")))
  paste(parties, collapse = "-")
}

edgelist_genus_manche$link_genus <- sapply(edgelist_genus_manche$link_genus, normaliser_paire)
write.csv2(edgelist_genus_manche,file="output/tableaux/Networks/edgelist_genus_manche.csv", row.names = FALSE,dec = ".")

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
write.csv2(freq_table_manche,file="output/tableaux/Networks/freq_table_manche.csv", row.names = FALSE,dec = ".")


# Create the alluvial graph
ggplot(freq_table_manche, aes(axis1 = v1_classe, axis2 = v2_classe, y = Number)) +
  geom_alluvium(aes(fill = v1_classe), alpha = 0.6, curve_type = "quintic") + 
  geom_stratum(aes(fill = after_stat(stratum)),colour="white",width = 0.05) +  
  scale_x_discrete(expand = c(-0.1, 0.1)) +
  scale_fill_manual(
    values = c(
      "Bacillariophyceae" = "#1f77b4", "Chlorophyta" = "#ff7f0e", "Chrysophyceae" = "#2ca02c",
      "Cryptophyceae" = "#d62728", "Cyanophyceae" = "cyan", "Dictyochophyceae" = "#9467bd", 
      "Dinophyceae" = "green", "Ebriophyceae" = "yellow", "Euglenozoa" = "red","Haptophyta" = "maroon",
      "Ciliophora" ="pink2","Raphidophyceae" = "magenta"
    )
  ) +
  labs(title = "Eastern Channel - North Sea global network associations",
       x = NULL,
       y = NULL, fill = "Phylum") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank(),  
    panel.grid = element_blank()  
  )
ggsave('Channel_alluvial_global.png', path = "output/graphs/networks/global_networks", dpi = 600, width = 200, height = 200, units = 'mm')


ggplot(freq_table_manche, aes(y = reorder(link_genus, Frequence), x = Frequence*100)) +
  geom_bar(stat = "identity", fill = "#CD9600") +  
  geom_text(aes(label = round(Frequence,digits = 3)*100), hjust = -0.2, vjust = 0.5, size = 3) +  
  labs(title = "Eastern Channel - North Sea global network associations",
       x = "Frequence (%)", 
       y = "Association") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
ggsave('Channel_freq_global.png', path = "output/graphs/networks/global_networks", dpi = 600, width = 200, height = 200, units = 'mm')

sum(filter(freq_table_manche,v1_classe=="Bacillariophyceae")$Frequence)
sum(filter(freq_table_manche,v1_classe=="Dinophyceae")$Frequence)


# Charger la matrice d'adjacence
assoMat <- read_delim("output/tableaux/Networks/assoMat_manche.csv", 
                      delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                          grouping_mark = ""), trim_ws = TRUE)

# Créer le graphe
cluster2 <- graph_from_adjacency_matrix(as.matrix(assoMat), weighted = TRUE, mode = "undirected", diag = FALSE)

# Clustering des nœuds
wcini <- cluster_fast_greedy(cluster2)

# Importer les informations taxonomiques
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt", header = TRUE, sep = "\t")

# Associer les clusters et taxonomies
compo_reseau <- data.frame(Taxon = wcini$names, cluster = wcini$membership) %>%
  left_join(phylum_classe, by = c("Taxon" = "Taxon")) %>%
  mutate(color = case_when(
    Phylum.Classe == "Bacillariophyceae" ~ "#56B4E9",
    Phylum.Classe == "Dinophyceae" ~ "#009E73",
    Phylum.Classe == "Ciliophora" ~ "#F0E442",
    Phylum.Classe == "Cryptophyceae" ~ "#CC79A7",
    Phylum.Classe == "Haptophyta" ~ "#996136",
    TRUE ~ "grey"
  ))


edgelist_genus_med <- read_delim("output/tableaux/Networks/edgelist_genus_manche.csv", 
                                 delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                                     grouping_mark = ""), trim_ws = TRUE)
# Adding some customization info for ggraph
V(cluster2)$Phylum <- compo_reseau$Phylum.Classe
V(cluster2)$size <- 15
E(cluster2)$weight <- E(cluster2)$weight * 5  
V(cluster2)$label <- V(cluster2)$name
V(cluster2)$color <- compo_reseau$color
E(cluster2)$info <- edgelist_genus_med$link_genus


# making the global network plot
manche_global <- ggraph(cluster2, layout = "auto") +  
  geom_edge_link(aes(edge_width = weight,edge_colour = info), alpha = 0.5, lineend = "round", linejoin = "bevel") +
  geom_node_point(aes(size = size, fill = color, stroke = 1), shape = 21,alpha = 0.5) +
  geom_node_text(aes(label = label, fontface = "bold"),colour = "black", 
                 repel = TRUE, size = 3, check_overlap = TRUE,
                 vjust = 0.5, hjust = 1) + 
  scale_edge_width(range = c(0.1, 2)) +
  scale_size_continuous(range = c(3, 10)) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c(
      "#56B4E9" = "#56B4E9",
      "#009E73" = "#009E73",
      "#F0E442" = "#F0E442",
      "#CC79A7" = "#CC79A7",
      "#996136" = "#996136",
      "grey" = "grey"
    )
  )+
  scale_edge_colour_manual(values = c(
    "Bacillariophyceae-Bacillariophyceae" = "#56B4E9","Bacillariophyceae-Dinophyceae"   = "chocolate1"    ,"Bacillariophyceae-Euglenozoa" = "grey"
    ,"Bacillariophyceae-Dictyochophyceae"  = "grey","Bacillariophyceae-Chlorophyta"   = "grey"   ,"Bacillariophyceae-Chrysophyceae"    = "grey"
    ,"Bacillariophyceae-Cyanophyceae"      = "grey","Bacillariophyceae-Haptophyta"    = "grey"   ,"Dinophyceae-Dinophyceae" = "#009E73"
    , "Dinophyceae-Euglenozoa"             = "grey", "Dinophyceae-Ebriophyceae"       = "grey"   , "Bacillariophyceae-Ebriophyceae"    = "grey"
    , "Bacillariophyceae-Cryptophyceae"    = "grey", "Dinophyceae-Haptophyta"         = "grey"   , "Chrysophyceae-Dinophyceae"         = "grey"
    , "Chlorophyta-Dinophyceae"            = "grey", "Cryptophyceae-Dinophyceae"      = "grey"   , "Cyanophyceae-Dinophyceae"          = "grey"
    , "Chlorophyta-Euglenozoa"             = "grey", "Cryptophyceae-Euglenozoa"       = "grey"   , "Chrysophyceae-Euglenozoa"          = "grey"
    , "Cyanophyceae-Euglenozoa"            = "grey", "Chlorophyta-Dictyochophyceae"   = "grey"   , "Dictyochophyceae-Dinophyceae"      = "grey"
    , "Cryptophyceae-Dictyochophyceae"     = "grey", "Dictyochophyceae-Haptophyta"    = "grey"   , "Dictyochophyceae-Euglenozoa"       = "grey"
    , "Chrysophyceae-Dictyochophyceae"     = "grey", "Cyanophyceae-Dictyochophyceae"  = "grey"   , "Chlorophyta-Cryptophyceae"         = "grey"
    , "Chlorophyta-Haptophyta"             = "grey", "Chlorophyta-Chrysophyceae"      = "grey"   , "Chlorophyta-Cyanophyceae"          = "grey"
    , "Cryptophyceae-Haptophyta"           = "grey", "Chrysophyceae-Cryptophyceae"    = "grey"   , "Cryptophyceae-Cyanophyceae"        = "grey"
    , "Euglenozoa-Haptophyta"              = "grey", "Chrysophyceae-Haptophyta"       = "grey"   , "Ebriophyceae-Haptophyta"           = "grey"
    , "Cyanophyceae-Haptophyta"            = "grey", "Ebriophyceae-Euglenozoa"        = "grey"   , "Euglenozoa-Euglenozoa"             = "grey"
    , "Chrysophyceae-Cyanophyceae"         = "grey", "Ebriophyceae-Ebriophyceae" = "grey"
  ))

ggsave('Manche_global_network.png', path = "output/graphs/networks/global_networks/", dpi = 600, width = 400, height = 300, units = 'mm')


##### Type of association through time ######
data <- read_delim("output/tableaux/Networks/subnetworks/results_metrics_reseaux_cluster2.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

data_association <- pivot_longer(data,cols = c(P_BacBac,P_BacDino,P_DinoDino,P_AAutres),names_to = "Asso_type" )
data_association <- select(data_association,Code_point_Libelle,Date,Asso_type,value)

data_association$MonthYear <- format(data_association$Date, "%Y-%m")
data_association$Month <- format(data_association$Date, "%m")

data_association_graph <- summarise(group_by(data_association,Code_point_Libelle,MonthYear,Asso_type,Month),value=mean(value,na.rm=T))

ggplot(data_association_graph) +
  geom_col(aes(y = value, x = MonthYear, fill = Asso_type), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")+
  geom_vline(data = subset(data_association_graph,Month == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", size = 0.5,linetype = "dashed")+
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  scale_fill_manual(values = c("P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  geom_text(data = subset(data_association_graph, Month == "01"), 
            aes(x = MonthYear, y = 0.95, label = sub("-.*", "", MonthYear)),
            color = "black", size = 2.5,angle=90, vjust = 0)+
  facet_wrap(~Code_point_Libelle,nrow = 4)+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 5))
ggsave('manche_typeassociation_station.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 400, height = 300, units = 'mm')

data_association_graph <- summarise(group_by(data_association,Code_point_Libelle,Asso_type,Month),value=mean(value,na.rm=T))

ggplot(data_association_graph) +
  geom_col(aes(y = value, x = Month, fill = Asso_type), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")+
  facet_wrap(~Code_point_Libelle,nrow = 4)
ggsave('manche_typeassociation_station_mois.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 400, height = 300, units = 'mm')

data_compo <- pivot_longer(data,cols = c(P_bac,P_dino,P_autres),names_to = "Taxon" )
data_compo <- select(data_compo,Code_point_Libelle,Date,Taxon,value)

data_compo$MonthYear <- format(data_compo$Date, "%Y-%m")
data_compo$Month <- format(data_compo$Date, "%m")

data_compo_graph <- summarise(group_by(data_compo,Code_point_Libelle,MonthYear,Taxon,Month),value=mean(value,na.rm=T))

ggplot(data_compo_graph) +
  geom_col(aes(y = value, x = MonthYear, fill = Taxon), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")+
  geom_vline(data = subset(data_compo_graph,Month == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", size = 0.5,linetype = "dashed")+
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  scale_fill_manual(values = c("P_bac" = "#56B4E9","P_autres"      = "grey" ,"P_dino" = "#009E73"))+
  geom_text(data = subset(data_compo_graph, Month == "01"), 
            aes(x = MonthYear, y = 0.95, label = sub("-.*", "", MonthYear)),
            color = "black", size = 2.5,angle=90, vjust = 0)+
  facet_wrap(~Code_point_Libelle,nrow = 4)+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 5))
ggsave('manche_compo_station.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 400, height = 300, units = 'mm')

data_compo_graph <- summarise(group_by(data_compo,Code_point_Libelle,Taxon,Month),value=mean(value,na.rm=T))

ggplot(data_compo_graph) +
  geom_col(aes(y = value, x = Month, fill = Taxon), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")+
  facet_wrap(~Code_point_Libelle,nrow = 4)
ggsave('manche_compo_station_mois.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 400, height = 300, units = 'mm')

data_association_graph <- summarise(group_by(data_association,MonthYear,Asso_type,Month),value=mean(value,na.rm=T))

ggplot(data_association_graph) +
  geom_col(aes(y = value, x = MonthYear, fill = Asso_type), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")+
  geom_vline(data = subset(data_association_graph,Month == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", size = 0.5,linetype = "dashed")+
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  scale_fill_manual(values = c("P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  geom_text(data = subset(data_association_graph, Month == "01"), 
            aes(x = MonthYear, y = 0.95, label = sub("-.*", "", MonthYear)),
            color = "black", size = 2.5,angle=90, vjust = 0)+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 5))
ggsave('manche_typeassociation.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 260, height = 170, units = 'mm')

data_association_graph <- summarise(group_by(data_association,Asso_type,Month),value=mean(value,na.rm=T))

ggplot(data_association_graph) +
  geom_col(aes(y = value, x = Month, fill = Asso_type), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")
ggsave('manche_typeassociation_mois.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 400, height = 300, units = 'mm')


data_compo_graph <- summarise(group_by(data_compo,MonthYear,Taxon,Month),value=mean(value,na.rm=T))

ggplot(data_compo_graph) +
  geom_col(aes(y = value, x = MonthYear, fill = Taxon), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")+
  geom_vline(data = subset(data_compo_graph,Month == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", size = 0.5,linetype = "dashed")+
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  scale_fill_manual(values = c("P_bac" = "#56B4E9","P_autres"      = "grey" ,"P_dino" = "#009E73"))+
  geom_text(data = subset(data_compo_graph, Month == "01"), 
            aes(x = MonthYear, y = 0.95, label = sub("-.*", "", MonthYear)),
            color = "black", size = 2.5,angle=90, vjust = 0)+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 5))
ggsave('manche_compo.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 260, height = 170, units = 'mm')

data_compo_graph <- summarise(group_by(data_compo,Taxon,Month),value=mean(value,na.rm=T))

ggplot(data_compo_graph) +
  geom_col(aes(y = value, x = Month, fill = Taxon), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")
ggsave('manche_compo_mois.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 400, height = 300, units = 'mm')


data_pca1 <- select(data,N_bac,N_dino,N_autres,N_BacBac,N_BacDino,N_DinoDino,N_AAutres)
PCA(data_pca1)

data_pca2 <- select(data,N_noeuds:N_clust,P_BacBac,P_BacDino,P_DinoDino,P_AAutres)
PCA(data_pca2)

data_pca3 <- select(data,N_noeuds:N_clust)
PCA(data_pca3)

data_pca4 <- select(data,N_noeuds:N_clust,N_BacBac,N_BacDino,N_DinoDino,N_AAutres)
PCA(data_pca4)

Table.corr_all.comp <- data_pca2[complete.cases(data_pca2),]
correlations <- cor(Table.corr_all.comp,method = "spearman")

# ... : Arguments supplémentaire à passer à la fonction cor.test
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

# Matrice de p-value de la corrélation
p.mat <- cor.mtest(Table.corr_all.comp)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(correlations, method="color", col=col(200),  
         type="upper", order="alphabet",tl.cex = 0.7,
         addCoef.col = "black", # Ajout du coefficient de corrélation
         tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
         # Combiner avec le niveau de significativité
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # Cacher les coefficients de corrélation sur la diagonale
         diag=F, 
         title = ""
)

Table.corr_all.comp <- data_pca4[complete.cases(data_pca4),]
correlations <- cor(Table.corr_all.comp,method = "spearman")

# ... : Arguments supplémentaire à passer à la fonction cor.test
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

# Matrice de p-value de la corrélation
p.mat <- cor.mtest(Table.corr_all.comp)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(correlations, method="color", col=col(200),  
         type="upper", order="alphabet",tl.cex = 0.7,
         addCoef.col = "black", # Ajout du coefficient de corrélation
         tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
         # Combiner avec le niveau de significativité
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # Cacher les coefficients de corrélation sur la diagonale
         diag=F, 
         title = ""
)



### Working on the Atlantic #####
# Import data
data <- read_delim("output/tableaux/Networks/edgelist_Atlantic.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt",header = T,sep = "\t")


edge1 <- as.data.frame(data$v1)
colnames(edge1) <- "Taxon"
edge1 <- left_join(edge1,phylum_classe)
colnames(edge1) <- c("v1","v1_classe")

edge2 <- as.data.frame(data$v2)
colnames(edge2) <- "Taxon"
edge2 <- left_join(edge2,phylum_classe)
colnames(edge2) <- c("v2","v2_classe")

edgelist_genus_atlantic <- cbind(edge1,edge2,data$asso)
edgelist_genus_atlantic$link_genus <- paste0(edgelist_genus_atlantic$v1_classe,"-",edgelist_genus_atlantic$v2_classe) 

normaliser_paire <- function(x) {
  parties <- sort(unlist(strsplit(x, "-")))
  paste(parties, collapse = "-")
}

edgelist_genus_atlantic$link_genus <- sapply(edgelist_genus_atlantic$link_genus, normaliser_paire)
write.csv2(edgelist_genus_atlantic,file="output/tableaux/Networks/edgelist_genus_atlantic.csv", row.names = FALSE,dec = ".")

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
write.csv2(freq_table_atlantic,file="output/tableaux/Networks/freq_table_atlantic.csv", row.names = FALSE,dec = ".")


# Create the alluvial graph
ggplot(freq_table_atlantic, aes(axis1 = v1_classe, axis2 = v2_classe, y = Number)) +
  geom_alluvium(aes(fill = v1_classe), alpha = 0.6, curve_type = "quintic") + 
  geom_stratum(aes(fill = after_stat(stratum)),colour="white",width = 0.05) +  
  scale_x_discrete(expand = c(-0.1, 0.1)) +
  scale_fill_manual(
    values = c(
      "Bacillariophyceae" = "#1f77b4", "Chlorophyta" = "#ff7f0e", "Chrysophyceae" = "#2ca02c",
      "Cryptophyceae" = "#d62728", "Cyanophyceae" = "cyan", "Dictyochophyceae" = "#9467bd", 
      "Dinophyceae" = "green", "Ebriophyceae" = "yellow", "Euglenozoa" = "red","Haptophyta" = "maroon",
      "Ciliophora" ="pink2","Raphidophyceae" = "magenta","Autres protistes" ="grey","Chromista"="tomato","Khakista"="blue"
    )
  ) +
  labs(title = "Atlantic - Western Channel global network associations",
       x = NULL,
       y = NULL, fill = "Phylum") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank(),  
    panel.grid = element_blank()  
  )
ggsave('Atlantic_alluvial_global.png', path = "output/graphs/networks/global_networks", dpi = 600, width = 200, height = 200, units = 'mm')


ggplot(freq_table_atlantic, aes(y = reorder(link_genus, Frequence), x = Frequence*100)) +
  geom_bar(stat = "identity", fill = "#00BE67") +  
  geom_text(aes(label = round(Frequence,digits = 3)*100), hjust = -0.2, vjust = 0.5, size = 3) +  
  labs(title = "Eastern Channel - North Sea global network associations",
       x = "Frequence (%)", 
       y = "Association") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
ggsave('Atlantic_freq_global.png', path = "output/graphs/networks/global_networks", dpi = 600, width = 200, height = 200, units = 'mm')

sum(filter(freq_table_atlantic,v1_classe=="Bacillariophyceae")$Frequence)
sum(filter(freq_table_atlantic,v1_classe=="Dinophyceae")$Frequence)

# Charger la matrice d'adjacence
assoMat <- read_delim("output/tableaux/Networks/assoMat_atlantic.csv", 
                      delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                          grouping_mark = ""), trim_ws = TRUE)

# Créer le graphe
cluster3 <- graph_from_adjacency_matrix(as.matrix(assoMat), weighted = TRUE, mode = "undirected", diag = FALSE)

# Clustering des nœuds
wcini <- cluster_fast_greedy(cluster3)

# Importer les informations taxonomiques
phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt", header = TRUE, sep = "\t")

# Associer les clusters et taxonomies
compo_reseau <- data.frame(Taxon = wcini$names, cluster = wcini$membership) %>%
  left_join(phylum_classe, by = c("Taxon" = "Taxon")) %>%
  mutate(color = case_when(
    Phylum.Classe == "Bacillariophyceae" ~ "#56B4E9",
    Phylum.Classe == "Dinophyceae" ~ "#009E73",
    Phylum.Classe == "Ciliophora" ~ "#F0E442",
    Phylum.Classe == "Cryptophyceae" ~ "#CC79A7",
    Phylum.Classe == "Haptophyta" ~ "#996136",
    TRUE ~ "grey"
  ))


edgelist_genus_med <- read_delim("output/tableaux/Networks/edgelist_genus_atlantic.csv", 
                                 delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                                     grouping_mark = ""), trim_ws = TRUE)
# Adding some customization info for ggraph
V(cluster3)$Phylum <- compo_reseau$Phylum.Classe
V(cluster3)$size <- 15
E(cluster3)$weight <- E(cluster3)$weight * 5  
V(cluster3)$label <- V(cluster3)$name
V(cluster3)$color <- compo_reseau$color
E(cluster3)$info <- edgelist_genus_med$link_genus


# making the global network plot
atlantic_global <-ggraph(cluster3, layout = "auto") +  
  geom_edge_link(aes(edge_width = weight,edge_colour = info), alpha = 0.5, lineend = "round", linejoin = "bevel") +
  geom_node_point(aes(size = size, fill = color, stroke = 1), shape = 21,alpha = 0.5) +
  geom_node_text(aes(label = label, fontface = "bold"),colour = "black", 
                 repel = TRUE, size = 3, check_overlap = TRUE,
                 vjust = 0.5, hjust = 1) + 
  scale_edge_width(range = c(0.1, 2)) +
  scale_size_continuous(range = c(3, 10)) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c(
      "#56B4E9" = "#56B4E9",
      "#009E73" = "#009E73",
      "#F0E442" = "#F0E442",
      "#CC79A7" = "#CC79A7",
      "#996136" = "#996136",
      "grey" = "grey"
    )
  )+
  scale_edge_colour_manual(values = c(
    "Bacillariophyceae-Bacillariophyceae" = "#56B4E9","Bacillariophyceae-Dinophyceae"   = "chocolate1"    ,"Bacillariophyceae-Euglenozoa" = "grey"
    ,"Bacillariophyceae-Dictyochophyceae"  = "grey","Bacillariophyceae-Chlorophyta"   = "grey"   ,"Bacillariophyceae-Chrysophyceae"    = "grey"
    ,"Bacillariophyceae-Cyanophyceae"      = "grey","Bacillariophyceae-Haptophyta"    = "grey"   ,"Dinophyceae-Dinophyceae" = "#009E73"
    , "Dinophyceae-Euglenozoa"             = "grey", "Dinophyceae-Ebriophyceae"       = "grey"   , "Bacillariophyceae-Ebriophyceae"    = "grey"
    , "Bacillariophyceae-Cryptophyceae"    = "grey", "Dinophyceae-Haptophyta"         = "grey"   , "Chrysophyceae-Dinophyceae"         = "grey"
    , "Chlorophyta-Dinophyceae"            = "grey", "Cryptophyceae-Dinophyceae"      = "grey"   , "Cyanophyceae-Dinophyceae"          = "grey"
    , "Chlorophyta-Euglenozoa"             = "grey", "Cryptophyceae-Euglenozoa"       = "grey"   , "Chrysophyceae-Euglenozoa"          = "grey"
    , "Cyanophyceae-Euglenozoa"            = "grey", "Chlorophyta-Dictyochophyceae"   = "grey"   , "Dictyochophyceae-Dinophyceae"      = "grey"
    , "Cryptophyceae-Dictyochophyceae"     = "grey", "Dictyochophyceae-Haptophyta"    = "grey"   , "Dictyochophyceae-Euglenozoa"       = "grey"
    , "Chrysophyceae-Dictyochophyceae"     = "grey", "Cyanophyceae-Dictyochophyceae"  = "grey"   , "Chlorophyta-Cryptophyceae"         = "grey"
    , "Chlorophyta-Haptophyta"             = "grey", "Chlorophyta-Chrysophyceae"      = "grey"   , "Chlorophyta-Cyanophyceae"          = "grey"
    , "Cryptophyceae-Haptophyta"           = "grey", "Chrysophyceae-Cryptophyceae"    = "grey"   , "Cryptophyceae-Cyanophyceae"        = "grey"
    , "Euglenozoa-Haptophyta"              = "grey", "Chrysophyceae-Haptophyta"       = "grey"   , "Ebriophyceae-Haptophyta"           = "grey"
    , "Cyanophyceae-Haptophyta"            = "grey", "Ebriophyceae-Euglenozoa"        = "grey"   , "Euglenozoa-Euglenozoa"             = "grey"
    , "Chrysophyceae-Cyanophyceae"         = "grey", "Ebriophyceae-Ebriophyceae" = "grey"
  ))

ggsave('Atlantic_global_network.png', path = "output/graphs/networks/global_networks/", dpi = 600, width = 400, height = 300, units = 'mm')

##### Type of association through time ######
data <- read_delim("output/tableaux/Networks/subnetworks/results_metrics_reseaux_cluster3.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

data_association <- pivot_longer(data,cols = c(P_BacBac,P_BacDino,P_DinoDino,P_AAutres),names_to = "Asso_type" )
data_association <- select(data_association,Code_point_Libelle,Date,Asso_type,value)

data_association$MonthYear <- format(data_association$Date, "%Y-%m")
data_association$Month <- format(data_association$Date, "%m")

data_association_graph <- summarise(group_by(data_association,Code_point_Libelle,MonthYear,Asso_type,Month),value=mean(value,na.rm=T))

ggplot(data_association_graph) +
  geom_col(aes(y = value, x = MonthYear, fill = Asso_type), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")+
  geom_vline(data = subset(data_association_graph,Month == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", size = 0.5,linetype = "dashed")+
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  scale_fill_manual(values = c("P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  geom_text(data = subset(data_association_graph, Month == "01"), 
            aes(x = MonthYear, y = 0.95, label = sub("-.*", "", MonthYear)),
            color = "black", size = 2.5,angle=90, vjust = 0)+
  facet_wrap(~Code_point_Libelle,nrow = 4)+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 5))
ggsave('atlantic_typeassociation_station.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 400, height = 300, units = 'mm')

data_association_graph <- summarise(group_by(data_association,Code_point_Libelle,Asso_type,Month),value=mean(value,na.rm=T))

ggplot(data_association_graph) +
  geom_col(aes(y = value, x = Month, fill = Asso_type), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")+
  facet_wrap(~Code_point_Libelle,nrow = 4)
ggsave('atlantic_typeassociation_station_mois.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 400, height = 300, units = 'mm')

data_compo <- pivot_longer(data,cols = c(P_bac,P_dino,P_autres),names_to = "Taxon" )
data_compo <- select(data_compo,Code_point_Libelle,Date,Taxon,value)

data_compo$MonthYear <- format(data_compo$Date, "%Y-%m")
data_compo$Month <- format(data_compo$Date, "%m")

data_compo_graph <- summarise(group_by(data_compo,Code_point_Libelle,MonthYear,Taxon,Month),value=mean(value,na.rm=T))

ggplot(data_compo_graph) +
  geom_col(aes(y = value, x = MonthYear, fill = Taxon), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")+
  geom_vline(data = subset(data_compo_graph,Month == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", size = 0.5,linetype = "dashed")+
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  scale_fill_manual(values = c("P_bac" = "#56B4E9","P_autres"      = "grey" ,"P_dino" = "#009E73"))+
  geom_text(data = subset(data_compo_graph, Month == "01"), 
            aes(x = MonthYear, y = 0.95, label = sub("-.*", "", MonthYear)),
            color = "black", size = 2.5,angle=90, vjust = 0)+
  facet_wrap(~Code_point_Libelle,nrow = 4)+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 5))
ggsave('atlantic_compo_station.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 400, height = 300, units = 'mm')

data_compo_graph <- summarise(group_by(data_compo,Code_point_Libelle,Taxon,Month),value=mean(value,na.rm=T))

ggplot(data_compo_graph) +
  geom_col(aes(y = value, x = Month, fill = Taxon), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")+
  facet_wrap(~Code_point_Libelle,nrow = 4)
ggsave('atlantic_compo_station_mois.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 400, height = 300, units = 'mm')

data_association_graph <- summarise(group_by(data_association,MonthYear,Asso_type,Month),value=mean(value,na.rm=T))

ggplot(data_association_graph) +
  geom_col(aes(y = value, x = MonthYear, fill = Asso_type), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")+
  geom_vline(data = subset(data_association_graph,Month == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", size = 0.5,linetype = "dashed")+
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  scale_fill_manual(values = c("P_BacBac" = "#56B4E9","P_BacDino"   = "chocolate1"
                               ,"P_AAutres"      = "grey" ,"P_DinoDino" = "#009E73"))+
  geom_text(data = subset(data_association_graph, Month == "01"), 
            aes(x = MonthYear, y = 0.95, label = sub("-.*", "", MonthYear)),
            color = "black", size = 2.5,angle=90, vjust = 0)+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 5))
ggsave('atlantic_typeassociation.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 260, height = 170, units = 'mm')

data_association_graph <- summarise(group_by(data_association,Asso_type,Month),value=mean(value,na.rm=T))

ggplot(data_association_graph) +
  geom_col(aes(y = value, x = Month, fill = Asso_type), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")
ggsave('atlantic_typeassociation_mois.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 400, height = 300, units = 'mm')

data_compo_graph <- summarise(group_by(data_compo,MonthYear,Taxon,Month),value=mean(value,na.rm=T))

ggplot(data_compo_graph) +
  geom_col(aes(y = value, x = MonthYear, fill = Taxon), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")+
  geom_vline(data = subset(data_compo_graph,Month == "01"), 
             aes(xintercept = MonthYear),
             color = "grey3", size = 0.5,linetype = "dashed")+
  scale_x_discrete(labels = rep(c("M","A","M","J","J","A","S","O","N","D","J","F"),16))+
  scale_fill_manual(values = c("P_bac" = "#56B4E9","P_autres"      = "grey" ,"P_dino" = "#009E73"))+
  geom_text(data = subset(data_compo_graph, Month == "01"), 
            aes(x = MonthYear, y = 0.95, label = sub("-.*", "", MonthYear)),
            color = "black", size = 2.5,angle=90, vjust = 0)+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 5))
ggsave('atlantic_compo.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 260, height = 170, units = 'mm')

data_compo_graph <- summarise(group_by(data_compo,Taxon,Month),value=mean(value,na.rm=T))

ggplot(data_compo_graph) +
  geom_col(aes(y = value, x = Month, fill = Taxon), position = "stack", na.rm = FALSE)+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")
ggsave('atlantic_compo_mois.png', path = "output/graphs/networks/subnetworks/", dpi = 600, width = 400, height = 300, units = 'mm')


data_pca1 <- select(data,N_bac,N_dino,N_autres,N_BacBac,N_BacDino,N_DinoDino,N_AAutres)
PCA(data_pca1)

data_pca2 <- select(data,N_noeuds:N_clust,P_BacBac,P_BacDino,P_DinoDino,P_AAutres)
PCA(data_pca2)

data_pca3 <- select(data,N_noeuds:N_clust)
PCA(data_pca3)

data_pca4 <- select(data,N_noeuds:N_clust,N_BacBac,N_BacDino,N_DinoDino,N_AAutres)
PCA(data_pca4)

Table.corr_all.comp <- data_pca2[complete.cases(data_pca2),]
correlations <- cor(Table.corr_all.comp,method = "spearman")

# ... : Arguments supplémentaire à passer à la fonction cor.test
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

# Matrice de p-value de la corrélation
p.mat <- cor.mtest(Table.corr_all.comp)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(correlations, method="color", col=col(200),  
         type="upper", order="alphabet",tl.cex = 0.7,
         addCoef.col = "black", # Ajout du coefficient de corrélation
         tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
         # Combiner avec le niveau de significativité
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # Cacher les coefficients de corrélation sur la diagonale
         diag=F, 
         title = ""
)

Table.corr_all.comp <- data_pca4[complete.cases(data_pca4),]
correlations <- cor(Table.corr_all.comp,method = "spearman")

# ... : Arguments supplémentaire à passer à la fonction cor.test
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

# Matrice de p-value de la corrélation
p.mat <- cor.mtest(Table.corr_all.comp)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(correlations, method="color", col=col(200),  
         type="upper", order="alphabet",tl.cex = 0.7,
         addCoef.col = "black", # Ajout du coefficient de corrélation
         tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
         # Combiner avec le niveau de significativité
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # Cacher les coefficients de corrélation sur la diagonale
         diag=F, 
         title = ""
)




### Compare composition of each global networks ####
med <- read_delim("output/tableaux/Networks/taxon_list_med.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)
med$region <- "1-Mediterranean sea"

med$Number <- 1

manche <- read_delim("output/tableaux/Networks/taxon_list_manche.csv", 
                  delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                      grouping_mark = ""), trim_ws = TRUE)
manche$region <- "2-Eastern Channel - North Sea"

manche$Number <- 1

atlantic <- read_delim("output/tableaux/Networks/taxon_list_atlantic.csv", 
                       delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                           grouping_mark = ""), trim_ws = TRUE)
atlantic$region <- "3-Atlantic - Western Channel"

atlantic$Number <- 1

data <- rbind(med,atlantic,manche)

phylum_classe <- read.table("data/Liste_phylum.classe_REPHY_modif.txt",header = T,sep = "\t")

data <- left_join(data,phylum_classe)

data$Phylum.Classe <- factor(data$Phylum.Classe, levels = c("Autres protistes","Chromista","Khakista","Chrysophyceae","Cyanophyceae","Dictyochophyceae",
                                                            "Raphidophyceae","Cryptophyceae", "Ebriophyceae","Haptophyta","Ciliophora",
                                                            "Chlorophyta","Euglenozoa","Dinophyceae","Bacillariophyceae"
                                                                                  ))

tableau_compo <- summarise(group_by(data,region,Phylum.Classe), Number = sum(Number)) |>
  mutate(Freq = Number / sum(Number))
write.csv2(tableau_compo,file="output/tableaux/Networks/compo_global_networks.csv", row.names = FALSE,dec = ".")


ggplot(tableau_compo) +
  geom_col(aes(y = region, x = Number, fill = Phylum.Classe), position = "stack", na.rm = FALSE, width = 0.7)+
  scale_fill_manual(values = c(
    "Bacillariophyceae" = "#1f77b4", "Chlorophyta" = "#ff7f0e", "Chrysophyceae" = "#2ca02c",
    "Cryptophyceae" = "#d62728", "Cyanophyceae" = "cyan", "Dictyochophyceae" = "#9467bd", 
    "Dinophyceae" = "green", "Ebriophyceae" = "yellow", "Euglenozoa" = "red","Haptophyta" = "maroon",
    "Ciliophora" ="pink2","Raphidophyceae" = "magenta","Autres protistes" ="grey","Chromista"="tomato","Khakista"="blue"
  ))+
  labs(x="Number of taxa",y="Global network",fill="Classe")+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")
ggsave('composition_number_global.png', path = "output/graphs/networks/global_networks", dpi = 600, width = 300, height = 200, units = 'mm')


ggplot(tableau_compo) +
  geom_col(aes(y = region, x = Freq*100, fill = Phylum.Classe), position = "stack", na.rm = FALSE, width = 0.7)+
  scale_fill_manual(values = c(
    "Bacillariophyceae" = "#1f77b4", "Chlorophyta" = "#ff7f0e", "Chrysophyceae" = "#2ca02c",
    "Cryptophyceae" = "#d62728", "Cyanophyceae" = "cyan", "Dictyochophyceae" = "#9467bd", 
    "Dinophyceae" = "green", "Ebriophyceae" = "yellow", "Euglenozoa" = "red","Haptophyta" = "maroon",
    "Ciliophora" ="pink2","Raphidophyceae" = "magenta","Autres protistes" ="grey","Chromista"="tomato","Khakista"="blue"
  ))+
  labs(x="Taxa proportion (%)",y="Global network",fill="Classe")+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")
ggsave('composition_freq_global.png', path = "output/graphs/networks/global_networks", dpi = 600, width = 300, height = 200, units = 'mm')


#plot_grid(med_global,manche_global,atlantic_global,nrow = 3)
#ggsave("global_networks.png", path = "output/graphs/networks/global_networks", dpi = 600, width = 300, height = 900, units = 'mm')


# Same but for the associations types 
med <- read_delim("output/tableaux/Networks/freq_table_med.csv", 
                  delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                      grouping_mark = ""), trim_ws = TRUE)
manche <- read_delim("output/tableaux/Networks/freq_table_manche.csv", 
                  delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                      grouping_mark = ""), trim_ws = TRUE)
atlantic <- read_delim("output/tableaux/Networks/freq_table_atlantic.csv", 
                  delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                      grouping_mark = ""), trim_ws = TRUE)

freq_table_big <- rbind(med,manche,atlantic)

freq_table_big$link_genus <- ifelse(freq_table_big$link_genus ==  "Bacillariophyceae-Bacillariophyceae", "Bacillariophyceae-Bacillariophyceae",
        ifelse(freq_table_big$link_genus ==  "Bacillariophyceae-Dinophyceae","Bacillariophyceae-Dinophyceae",
          ifelse(freq_table_big$link_genus ==  "Dinophyceae-Dinophyceae","Dinophyceae-Dinophyceae" ,"Autres"
       )))

freq_table_big$link_genus <- factor(freq_table_big$link_genus, levels = c("Autres","Dinophyceae-Dinophyceae","Bacillariophyceae-Dinophyceae","Bacillariophyceae-Bacillariophyceae"))

# Modifier dans les freq_table pour avoir que certaines et ensuite passer le reste en autre
# Faire le graphe en dessous 

ggplot(freq_table_big) +
  geom_col(aes(y = region, x = Frequence, fill = link_genus), position = "stack", na.rm = FALSE, width = 0.7)+
  labs(x="Frequence",y="Global network",fill="Associations")+
  scale_fill_manual(values = c(
    "Bacillariophyceae-Bacillariophyceae" = "#56B4E9","Bacillariophyceae-Dinophyceae"   = "chocolate1"
    ,"Autres"      = "grey" ,"Dinophyceae-Dinophyceae" = "#009E73"

  ))+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")
ggsave("global_networks_associations.png", path = "output/graphs/networks/global_networks", dpi = 600, width = 400, height = 300, units = 'mm')

ggplot(freq_table_big) +
  geom_col(aes(y = region, x = Number, fill = link_genus), position = "stack", na.rm = FALSE, width = 0.7)+
  labs(x="Number of association",y="Global network",fill="Associations")+
  scale_fill_manual(values = c(
    "Bacillariophyceae-Bacillariophyceae" = "#56B4E9","Bacillariophyceae-Dinophyceae"   = "chocolate1"
    ,"Autres"      = "grey" ,"Dinophyceae-Dinophyceae" = "#009E73"
    
  ))+
  guides(fill = guide_legend(ncol = 7, byrow = F)) +
  theme(legend.position = "bottom")
ggsave("global_networks_associations_number.png", path = "output/graphs/networks/global_networks", dpi = 600, width = 400, height = 300, units = 'mm')


## Correlation between diversity index and associations, compositions of subnetworks ######
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetricsfinal.csv", 
                       delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                           grouping_mark = ""), trim_ws = TRUE)

data_med <- filter(data, region == "1-Mediterranean sea")

data_med_1 <- select(data_med,Shannon:Rspe,P_BacBac,P_DinoDino,P_BacDino,P_AAutres)

Table.corr_all.comp <- data_med_1[complete.cases(data_med_1),]
correlations <- cor(Table.corr_all.comp,method = "spearman")

# ... : Arguments supplémentaire à passer à la fonction cor.test
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

# Matrice de p-value de la corrélation
p.mat <- cor.mtest(Table.corr_all.comp)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(correlations, method="color", col=col(200),  
         type="upper", order="alphabet",tl.cex = 0.7,
         addCoef.col = "black", # Ajout du coefficient de corrélation
         tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
         # Combiner avec le niveau de significativité
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # Cacher les coefficients de corrélation sur la diagonale
         diag=F, 
         title = ""
)

data_med_2 <- select(data_med,Shannon:Rspe,N_BacBac,N_DinoDino,N_BacDino,N_AAutres)

Table.corr_all.comp <- data_med_2[complete.cases(data_med_2),]
correlations <- cor(Table.corr_all.comp,method = "spearman")

# ... : Arguments supplémentaire à passer à la fonction cor.test
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

# Matrice de p-value de la corrélation
p.mat <- cor.mtest(Table.corr_all.comp)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(correlations, method="color", col=col(200),  
         type="upper", order="alphabet",tl.cex = 0.7,
         addCoef.col = "black", # Ajout du coefficient de corrélation
         tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
         # Combiner avec le niveau de significativité
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # Cacher les coefficients de corrélation sur la diagonale
         diag=F, 
         title = ""
)

data_med_3 <- select(data_med,Shannon:Rspe,N_bac,N_dino,N_autres)

Table.corr_all.comp <- data_med_3[complete.cases(data_med_3),]
correlations <- cor(Table.corr_all.comp,method = "spearman")

# ... : Arguments supplémentaire à passer à la fonction cor.test
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

# Matrice de p-value de la corrélation
p.mat <- cor.mtest(Table.corr_all.comp)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(correlations, method="color", col=col(200),  
         type="upper", order="alphabet",tl.cex = 0.7,
         addCoef.col = "black", # Ajout du coefficient de corrélation
         tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
         # Combiner avec le niveau de significativité
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # Cacher les coefficients de corrélation sur la diagonale
         diag=F, 
         title = ""
)

data_med_4 <- select(data_med,Shannon:Rspe,P_bac,P_dino,P_autres)

Table.corr_all.comp <- data_med_4[complete.cases(data_med_4),]
correlations <- cor(Table.corr_all.comp,method = "spearman")

# ... : Arguments supplémentaire à passer à la fonction cor.test
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

# Matrice de p-value de la corrélation
p.mat <- cor.mtest(Table.corr_all.comp)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(correlations, method="color", col=col(200),  
         type="upper", order="alphabet",tl.cex = 0.7,
         addCoef.col = "black", # Ajout du coefficient de corrélation
         tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
         # Combiner avec le niveau de significativité
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # Cacher les coefficients de corrélation sur la diagonale
         diag=F, 
         title = ""
)

### PCA on the metrics
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_final.csv", 
                       delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                           grouping_mark = ""), trim_ws = TRUE)
data_pca <- select(data,N_noeuds:N_clust)
PCA_results <- PCA(data_pca)
fviz_screeplot(PCA_results)
fviz_pca_ind(PCA_results,col.ind = data$region,addEllipses = T,alpha.ind = 1,geom = c("point"))
fviz_pca_var(PCA_results,axes = c(3,2))
PCA_coord <- PCA_results$ind$coord

PCA_coord <- as.data.frame(PCA_coord)
PCA_coord <- select(PCA_coord, Dim.1,Dim.2,Dim.3)
write.csv2(PCA_coord,file="output/tableaux/Networks/PCA_coords.csv", row.names = FALSE,dec = ".")


#### correlation with the axis to summarise the different metrics

data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_withmetrics&networksdiv_PCA_coord.final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)
data_cor <- select(data,N_noeuds:N_clust,Dim.1:Dim.3)

Table.corr_all.comp <- data_cor[complete.cases(data_cor),]
correlations <- cor(Table.corr_all.comp,method = "spearman")

# ... : Arguments supplémentaire à passer à la fonction cor.test
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

# Matrice de p-value de la corrélation
p.mat <- cor.mtest(Table.corr_all.comp)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(correlations, method="color", col=col(200),  
         type="upper", order="original",tl.cex = 0.7,
         addCoef.col = "black", # Ajout du coefficient de corrélation
         tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
         # Combiner avec le niveau de significativité
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # Cacher les coefficients de corrélation sur la diagonale
         diag=F, 
         title = ""
)
