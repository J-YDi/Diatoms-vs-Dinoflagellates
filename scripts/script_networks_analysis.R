# Script Suite Stage JY Dias # 10/12/2024

library(ggalluvial)

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

