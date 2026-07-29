# Diatoms vs dinoflagellates: a temporal network analysis of bloom impacts on phytoplankton diversity and  community structure in French coastal waters | R script and generated data
## Dias Jean-Yves $^1,^2$, Pochic Victor $^1$ $^,$ $^3$ , Chaffron Samuel $^4$ $^,$ $^5$, Gernez Pierre $^1$
#### $^1$ Nantes Université, Institut des Substances et Organismes de la Mer, ISOMER, UR 2160, F-44000 Nantes, France ; $^2$ Laboratoire de Biologie des Organismes et des Écosystèmes Aquatiques-BOREA, Muséum national d’Histoire naturelle (MNHN), SU, CNRS, IRD, UA, F-75005 Paris, France ; $^3$ Ifremer, COAST, F-44000 Nantes, France ; $^4$ Nantes Université, École Centrale Nantes, CNRS, LS2N, UMR 6004, F-44000 Nantes, France ; $^5$ Research Federation for the Study of Global Ocean Systems Ecology and Evolution, FR2022/Tara Oceans GOSEE, F-75016 Paris, France
## **Article access : https://doi.org/10.1093/ismeco/ycag174

#### Github repository organization

##### Scripts folder :
+ 0_data_import.R : Script used to manipulate original datasets and create new ones
+ 1_script_clustering_stations.R: Script used to create the different regions
+ 2_script_build_networks.R : Script used to build the global and temporal networks with the different metrics
+ 3_script_maps.R : Script used to create the map
+ 4_script_Chlatimeseries_processing.R : Script to detect blooms
+ 5_script_region_analysis.R : Script to investigate the differences between regions
+ 6_script_networks_analysis.R : Script to analyze composition and associations type in the networks and make a PCA to sum up the different metrics
+ 7_script_randomizationtest.R : Script to investigate the relationship between composition and associations type and test of randomization
+ 8_script_blooms_analysis.R : Script to analyze the effet of a bloom on the different variables.
+ 9_script_negative&positive_asso.R : Script to analyze both negative and positive associations.
+ Supp_script_interactive_figure_2.R : Script to generate interactive versions of Figure 2.

##### data_modif folder : 
+ contains all the original dataset modified for the study, all can be obtained by the different scripts and the original raw data (available at SEANOE, thanks to IFREMER, https://doi.org/10.17882/47248 or upon request)

##### tableaux folder : 
+ contains all dataset created for the study,mainly results and intermediate results, all can be obtained by the different scripts and the files from the data_modif folder **except Table_bloom_R_v3c.csv**



###### Contact : jean-yves.dias@sorbonne-universite.fr
