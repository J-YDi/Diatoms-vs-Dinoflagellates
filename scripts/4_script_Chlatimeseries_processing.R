# Script Suite Stage JY Dias # 10/12/2024

# Loading packages
library(readr)
library(interp)
library("pastecs")
library("stlplus")
library(trend)
library(dplyr)

# Import data
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Create the function to detect outliers
outliers <- function(x) {
  median_dev <- abs(x - median(x,na.rm = T))
  indexes <- which(median_dev > (2.323561*1.486) * mad(x,na.rm = T))
  return(indexes)
}

# Loop to regularize, deseasonalize and detect outliers for each station

for (i in c(1:21)){
  # Determine the station to analyze
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  # Keep it only
  Table <- filter(data, Code_point_Libelle == station)
  #windows(title = "A regulariser")
  #print(ggplot(Table) + geom_path(aes(x=Date, y=CHLOROA)) + geom_rug(aes(x=Date)))
  # Make it by chronological order
  Table <- Table[order(Table$Date), ]
  # Determine the mean interval between sampling date
  interval <- as.numeric(mean(diff(Table$Date)))
  # and the minimal interval
  tmin <- min(Table$Date)
  # Regularize it
  t_reg <- regul(x=Table$Date, y=Table$CHLOROA, xmin=tmin, n=nrow(Table), deltat=interval)
  #windows(title = "regul step 1")
  #plot(t_reg)
  # Line to test different parameter to make the regularization with the less interpolation
  data_regul_ini <- regul.screen(x=as.numeric(Table$Date),
                                 xmin=seq(from=as.numeric(tmin-10), to=as.numeric(tmin+10), by=1),
                                 deltat=seq(from=10, to=20, by=1),
                                 tol = 3
  )
  data_regul <- as.data.frame(data_regul_ini$nbr.match)
  max_index <- which(data_regul == max(data_regul, na.rm = TRUE), arr.ind = TRUE)
  # New starting date and new interval
  newtmin <- as.Date(as.numeric(sub(".*=", "", rownames(data_regul)[max_index[1, 1]])))
  newinterval <- as.numeric(sub(".*=", "", colnames(data_regul)[max_index[1, 2]]))
  # Regularize with this
  data_regul_n <- as.data.frame(data_regul_ini$n)
  new_n <- data_regul_n[max_index[1, 1],max_index[1, 2]]
  t_reg <- regul(x=Table$Date, y=Table$CHLOROA, xmin=newtmin, n=new_n, deltat=newinterval,method="linear")
  #windows(title = "regul step 2")
  #plot(t_reg)
  Table_regul <- as.data.frame(cbind(t_reg$x,t_reg$y))
  colnames(Table_regul) <- c("Date","CHLOROA")
  #windows(title = "Regulated ts")
  #plot(Table_regul,type="o")
  CHLOROA_ts <- extract(t_reg)
  #plot(CHLOROA_ts, type="o")
  # Compute the autocorrelation test
  acf_CHLOROA <- acf(CHLOROA_ts,na.action = na.omit,plot = F)
  ACF_data <- as.data.frame(cbind(acf_CHLOROA$acf,acf_CHLOROA$lag))
  lag_data <- filter(ACF_data,V1<0.20)[1,2]
  Table_regul$year <- format(Table_regul$Date, "%Y")
  sampling <- count(Table_regul, year)
  np <- median(sampling[,2])
  # Deseasonalization with STL
  dec <- stlplus(Table_regul$CHLOROA,t=Table_regul$Date, n.p=np, s.window="periodic", t.window=NULL)
  plot(dec, scales=list(y="free"))
  Table_regul$CHLOROA_noseason <- dec$data$raw - dec$data$seasonal
  
  residual <- dec$data$remainder
  plot(acf(residual,na.action = na.pass,plot=F),main = paste0("Residuals structure ",station))
  
  outliers_CHLOROA<- outliers(Table_regul$CHLOROA_noseason)
  #windows(title = "Outliers")
  plot(y=Table_regul$CHLOROA_noseason,x=Table_regul$Date, type="l",main = paste0("Outliers",station))
  points(
    # the x coordinate is the index
    x=Table_regul[outliers_CHLOROA,]$Date,
    # the y coordinate is the value
    y=Table_regul$CHLOROA_noseason[outliers_CHLOROA],
    # make points solid and red
    pch=16, col="red")
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  # Indicate if this date is an outlier or not
  Table_regul$Code_point_Libelle <- station
  Table_regul$Outlier <- "NON"
  Table_regul[outliers_CHLOROA,]$Outlier <- "OUI"
  # Save it
  write.csv2(Table_regul,file=paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), row.names = FALSE,dec = ".")
  
}


# Find the date that corresponds to the outliers: it is the date logically closest to the point
# Loading data
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

# Loading the result of Script_regul_deseason_trend_change.R for all the stations
{i=1
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_Toulon <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                            delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=2
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_Ansecarteau <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                                 delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=3
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_Antifer <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                             delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=4
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_Atso <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                          delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=5
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_Auger <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                           delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=6
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_Barcares <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                              delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=7
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_Boisdelachaise <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                                    delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=8
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_Bouzigues <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                               delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=9
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_Cabourg <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                             delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=10
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_Calvi <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                           delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=11
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_Dianacentre <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                                 delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  i=12
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_Géfosse <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                             delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=13
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_Cornard <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                             delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=14
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_Hebihens <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                              delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=15
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_Loguivy <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                             delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=16
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_MenerRoue <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                               delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=17
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_OuestLoscolo <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                                  delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=18
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_ParcLeucate <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                                 delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=19
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_Boulogne <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                              delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=20
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_Setemer <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                             delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  i=21
  station <- levels(as.factor(data$Code_point_Libelle))[i]
  print(station)
  nom_fichier <- paste0("Outliers_regularise_desaisonalise_CHLOROA",station)
  nom_fichier <- paste0(nom_fichier)
  
  data_Teychan <- read_delim(paste0("output/tableaux/CHLAprocessing/",nom_fichier,".csv"), 
                             delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = ""), trim_ws = TRUE)
  
  data_outliers_change <- rbind( data_Ansecarteau  ,    data_Antifer        ,  data_Atso       ,      data_Auger     ,     
                                 data_Barcares     ,    data_Boisdelachaise ,  data_Boulogne   ,      data_Cabourg   ,       data_Calvi,          
                                 data_Cornard     ,     data_Dianacentre   ,   data_Géfosse   ,       data_Hebihens ,        data_Loguivy  ,      
                                 data_MenerRoue   ,     data_OuestLoscolo  ,   data_Setemer   ,       data_Teychan  ,data_Toulon, data_ParcLeucate, data_Bouzigues)
  
}


data_outliers <- filter(data_outliers_change, Outlier =="OUI")

# Find the closest date to the interpolate value which correspond the best to the outlier date
for (k in 1:21){
  data_Date_find <- c("","")
  data_Date_find <- as.data.frame(data_Date_find)
  
  data_Date_tosearch <- c("","")
  data_Date_tosearch <- as.data.frame(data_Date_tosearch)
  j=3
  station <- levels(as.factor(data_outliers$Code_point_Libelle))[k]
  Tabletosearch <- filter(data_outliers, Code_point_Libelle == station)
  Tabletofind <- filter(data, Code_point_Libelle == station)
  for (j in 1:nrow(Tabletosearch)){
    datetosearch <- as.Date(levels(as.factor(Tabletosearch$Date))[j])
    datefind <- Tabletofind$Date[which.min(abs(Tabletofind$Date - datetosearch))]
    
    data_Date_tosearch[j,1] <- as.numeric(datetosearch)
    data_Date_find[j,1] <- as.numeric(datefind)
    data_Date_tosearch[j,2] <- station
    data_Date_find[j,2] <- station
    
  }
  data_Date_ok <- cbind(data_Date_tosearch,data_Date_find)
  colnames(data_Date_ok) <- c("Datetosearch","Code_point_Libelle","Datefind","Station")
  data_Date_ok$Datetosearch <- as.Date(as.numeric(data_Date_ok$Datetosearch))
  data_Date_ok$Datefind <- as.Date(as.numeric(data_Date_ok$Datefind))
  
  nom_fichier <- paste0("Outliers_matchingdate",station)
  nom_fichier <- paste0(nom_fichier)
  write.csv2(data_Date_ok,file=paste0("output/tableaux/CHLAprocessing/Matchingdate_outliers/",nom_fichier,".csv"), row.names = FALSE,dec = ".")
  
}


# Indicate on the initial data where the outliers are
# Import data

files <- list.files("output/tableaux/CHLAprocessing/Matchingdate_outliers/", full.names=TRUE)
# count how many we have
length(files)
realdate <- tibble()
# Import all the data on one dataframe
for ( i in 1:length(files) ) {
  message("I am reading file number ", i)
  # get the file path
  file <- files[i]
  # read the file
  temp <- read_delim(file, 
                     delim = ";", escape_double = FALSE, trim_ws = TRUE)
  # combine it with the previously read data
  realdate <- bind_rows(realdate, temp)
}
# At this date there is an outliers
realdate$Outlier <- "OUI"
realdate_forfusion <- dplyr::select(realdate,Code_point_Libelle,Datefind,Outlier)
colnames(realdate_forfusion) <- c("Code_point_Libelle","Date","Outlier")

# Join the data and the info of there is a bloom !
# Load data
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)
data_withoutliers <- left_join(data,realdate_forfusion)
# The fact of having made the correspondence between "false dates" and real ones means that several false
# can correspond to only 1 "real date" which is therefore duplicated

# Delete the duplicates
doublons_final <- data_withoutliers[duplicated(data_withoutliers$ID.interne.passage) |
                                      duplicated(data_withoutliers$ID.interne.passage, fromLast = TRUE), ]

resultat_filtre_final <- doublons_final %>%
  filter(duplicated(ID.interne.passage) | n()==1)

data_withoutliers_unique <- subset(data_withoutliers, !(ID.interne.passage %in% unique(doublons_final$ID.interne.passage)))
data_withoutliers_ok <- bind_rows(data_withoutliers_unique,resultat_filtre_final)

# Still have duplicates
doublons_final <- data_withoutliers_ok[duplicated(data_withoutliers_ok$ID.interne.passage) |
                                         duplicated(data_withoutliers_ok$ID.interne.passage, fromLast = TRUE), ]

resultat_filtre_final <- doublons_final %>%
  filter(duplicated(ID.interne.passage) | n()==1)

data_withoutliers_unique <- subset(data_withoutliers_ok, !(ID.interne.passage %in% unique(doublons_final$ID.interne.passage)))
data_withoutliers_ok <- bind_rows(data_withoutliers_unique,resultat_filtre_final)
data_withoutliers_ok <- data_withoutliers_ok |>
  group_by(Code.Region, Code_point_Libelle, lon, lat, Year, Month, Date, ID.interne.passage, Prelevement.niveau)

# Still have duplicates....
doublons_final <- data_withoutliers_ok[duplicated(data_withoutliers_ok$ID.interne.passage) |
                                         duplicated(data_withoutliers_ok$ID.interne.passage, fromLast = TRUE), ]

resultat_filtre_final <- doublons_final %>%
  filter(duplicated(ID.interne.passage) | n()==1)

data_withoutliers_unique <- subset(data_withoutliers_ok, !(ID.interne.passage %in% unique(doublons_final$ID.interne.passage)))

data_withoutliers_ok <- bind_rows(data_withoutliers_unique,resultat_filtre_final)

data_withoutliers_ok <- data_withoutliers_ok |>
  group_by(Code.Region, Code_point_Libelle, lon, lat, Year, Month, Date, ID.interne.passage, Prelevement.niveau)

# That's ok now, save it
write.csv2(data_withoutliers_ok,file="output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_blooms_final.csv", row.names = FALSE,dec = ".")


### Associate blooming genus with bloom information ###
data <- read_delim("output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_blooms_final.csv", 
                   delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                       grouping_mark = ""), trim_ws = TRUE)

data <- dplyr::select(data,-Outlier)

##########################################
# Table_bloom_R_v3c.csv was made by hand #
##########################################
# Loading information about the different blooms
Table_bloom_R <- read_delim("output/tableaux/Table_bloom_R_v3c.csv", 
                            delim = ";", escape_double = FALSE, col_types = cols(Date = col_date(format = "%d/%m/%Y")), 
                            trim_ws = TRUE)

Table_bloom_R <- Table_bloom_R[complete.cases(Table_bloom_R$Code_point_Libelle),] 
Table_bloom_R <- dplyr::select(Table_bloom_R, Code_point_Libelle, Date,Abdtot:Bloom_Genre)


data_ok <- left_join(data,Table_bloom_R, join_by(Code_point_Libelle, Date))
# Save it
write.csv2(data_ok,file="output/data_modif/Table_FLORTOT_Surf_0722_COM_period_Stselect_hydro_phyto_chloro_phylum_period15_chlafilter_cluster5_blooms_caracterised_final.csv", row.names = FALSE,dec = ".")

