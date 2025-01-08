# Script Suite Stage JY Dias # 09/12/2024

# Loading packages
library(cluster)
library(ClustGeo)
library(sp)
library(readr)
library(FactoMineR)
library(missMDA)
library(factoextra)
library(dplyr)
library(maps)
library(mapdata)
library(ggplot2)
library(ggthemes)

# Prerequisites hclust geo function modified to have access to the compromised distance matrix
HCLUSTMODIF <- function (D0, D1 = NULL, alpha = 0, scale = TRUE, wt = NULL) 
{
  if (class(D0) != "dist") 
    stop("DO must be of class dist (use as.dist)", call. = FALSE)
  if (!is.null(D1) && (class(D1) != "dist")) 
    stop("D1 must be of class dist (use as.dist)", call. = FALSE)
  n <- as.integer(attr(D0, "Size"))
  if (is.null(n)) 
    stop("invalid dissimilarities", call. = FALSE)
  if (is.na(n) || n > 65536L) 
    stop("size cannot be NA nor exceed 65536", call. = FALSE)
  if (n < 2) 
    stop("must have n >= 2 objects to cluster", call. = FALSE)
  if (!is.null(D1) && length(D0) != length(D1)) 
    stop("the two dissimilarity structures must have the same size", 
         call. = FALSE)
  if ((max(alpha) > 1) || (max(alpha) < 0)) 
    stop("Values alpha must be in [0,1]", call. = FALSE)
  if ((scale == TRUE) && (!is.null(D1))) {
    D0 <- D0/max(D0)
    D1 <- D1/max(D1)
  }
  delta0 <- wardinit(D0, wt)
  if (!is.null(D1)) 
    delta1 <- wardinit(D1, wt)
  else delta1 <- 0
  delta <- (1 - alpha) * delta0 + alpha * delta1
  res <- delta
  return(res)
}



##### Spatial clustering based on all the stations that had 5 years of continuity for the chlorophyll measures #####
# Import data
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_9523_Stselect_hydro_phyto_chloro_phylum_period5_final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Select abiotic variables 
Table_PCA <- data[, c("SALI","TEMP","TURB","NH4","PO4","SIOH","OXYGENE","NO3+NO2","TURB-FNU")]

# Imputation of missing hydro values
ncomp <- estim_ncpPCA(Table_PCA)
res.imp <- imputePCA(Table_PCA, ncp = ncomp$ncp,method = "EM")

# Doing the PCA
PCA_result <- PCA(res.imp$completeObs, ncp=20)
# We retrieve the coordinates
coord_PCA <- PCA_result$ind$coord

# Let's see the screeplot
fviz_eig(PCA_result,addlabels = T)
# We create a dataset with the coordinates in the dimensions of the PCA
Table_abio_PCA <- cbind(data$Code_point_Libelle,data$ID.interne.passage, coord_PCA)
Table_abio_PCA <- as.data.frame(Table_abio_PCA)
Table_abio_PCA$Dim.1 <- as.numeric(Table_abio_PCA$Dim.1)
Table_abio_PCA$Dim.2 <- as.numeric(Table_abio_PCA$Dim.2)
Table_abio_PCA$Dim.3 <- as.numeric(Table_abio_PCA$Dim.3)
Table_abio_PCA$Dim.4 <- as.numeric(Table_abio_PCA$Dim.4)
Table_abio_PCA$Dim.5 <- as.numeric(Table_abio_PCA$Dim.5)
Table_abio_PCA$Dim.6 <- as.numeric(Table_abio_PCA$Dim.6)
Table_abio_PCA$Dim.7 <- as.numeric(Table_abio_PCA$Dim.7)
Table_abio_PCA$Dim.8 <- as.numeric(Table_abio_PCA$Dim.8)
Table_abio_PCA$Dim.9 <- as.numeric(Table_abio_PCA$Dim.9)

colnames(Table_abio_PCA)[1:2] <- c("Code_point_Libelle","ID.interne.passage")

# We take the median of the positions of each station in each dimension of the PCA
Table_abio_st_med <- summarise(group_by(Table_abio_PCA, Code_point_Libelle), median_dim1=median(Dim.1,na.rm=T)
                               , median_dim2=median(Dim.2,na.rm=T),
                               median_dim3=median(Dim.3,na.rm=T),
                               median_dim4=median(Dim.4,na.rm=T),
                               median_dim5=median(Dim.5,na.rm=T),
                               median_dim6=median(Dim.6,na.rm=T),
                               median_dim7=median(Dim.7,na.rm=T),
                               median_dim8=median(Dim.8,na.rm=T),
                               median_dim9=median(Dim.9,na.rm=T))

# We create a dataset which contains the geographical coordinates of the stations
data_geo <- summarise(group_by(data, Code_point_Libelle,lon,lat))
# Delete some duplicates because of different GPS coordinates 
data_geo <- data_geo[-c(10,12:19,26:28,42:51,58:72),]
data_geo$lat <- as.numeric(data_geo$lat)
data_geo$lon <- as.numeric(data_geo$lon)


# Define the number of clusters
K <- 5

# Calculate distance matrix on hydro data
D0 <- dist(Table_abio_st_med)
tree <- hclustgeo(D0)

# Convert geographic positions to spatial object
firstPoints <- SpatialPoints(coords = cbind(data_geo$lat,data_geo$lon))
# Make a spatial distance between objects
euclidDist <- sp::spDists(firstPoints,longlat = FALSE)

# Distance matrix on geography
D1 <- as.dist(euclidDist)

# Choice of alpha parameter
range.alpha <- seq(0,1,0.1)
cr <- choicealpha(D0, D1, range.alpha, 
                  K, graph = FALSE)
cr$Q # proportion of explained inertia
# Graph associated
plot(cr,xlab="alpha",ylab="inertia Q")
# We choose alpha = 0.1
tree <- hclustgeo(D0,D1,alpha=0.1)
# HCLUSTMODIF function to extract the compromise distance matrix between D0 and D1
diss <- HCLUSTMODIF(D0,D1,alpha=0.1)
# Tree for the clustering HYDROGEO
tree[["labels"]] <- Table_abio_st_med$Code_point_Libelle
# Plot it
plot(tree,hang = -1, 
     xlab = "", sub = "",
     main = "Ward dendrogram with alpha 0.1 and 5 clusters")
# Defining the clusters
cluster <- cutree(tree,K)
data_clust_result <- as.data.frame(cbind(data_geo,cluster))
colnames(data_clust_result)[4] <- "cluster"


# We evaluate the quality of the clustering
silhouette_score <- silhouette(as.numeric(data_clust_result$cluster), diss)
# Result :
cat("Indice de Silhouette global 0.1 :", mean(silhouette_score[, "sil_width"]), "\n")

# Map to show the different cluster
Worldmap <- map_data('worldHires')
ggplot() + geom_polygon(data = Worldmap, aes(x = long, y = lat, group = group), fill = 'gray', color = 'gray10', size = .25)+
  coord_fixed(xlim=c(-5.5,9.5), ylim=c(41,51.5), ratio=1.4)+
  labs(y = 'Latitude (degrés)', x = 'Longitude (degrés)',
       subtitle = "EM-0,1-k5",colour = "Cluster")+
  theme_gdocs()+
  geom_point(data = data_clust_result, aes(x = lon, y = lat,colour=as.character(cluster)), size =5)+
  theme(panel.grid.major = element_line(color = 'gray10', size = .25), panel.grid.minor = NULL, panel.ontop = FALSE,
        panel.background = element_rect(fill = 'lightblue2'))+
  guides(color = guide_legend(override.aes = list(size = 10)))
#ggsave('maps_clust_EM_01_k5.png', path = "output/graphs/Clusterings", dpi = 600, width = 400, height = 200, units = 'mm')

#Save the clusters 
data_clust_result <- select(data_clust_result,Code_point_Libelle,cluster)
write.csv2(data_clust_result,file="output/tableaux/clusters_EM_01_k5_final.csv", row.names = FALSE,dec = ".")
