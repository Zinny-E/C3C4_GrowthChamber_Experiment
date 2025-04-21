###############################################################################
## Libraries
###############################################################################
library(dplyr)
library(tidyverse)
library(LeafArea)
library(plantecophys)
library(car)
library(emmeans)
library(readxl)

###############################################################################
## Import files
###############################################################################
#biomass_area <- read.csv("../data_sheets/NxCO2_tla_biomass_data.csv")
#id <- read.csv("../data_sheets/NxCO2_id_datasheet.csv")
isotopes <- read.csv("../CO2xTempleafIsotopes.csv")
totalN <- read.csv("../totalN.csv")
d13c.air <- read.csv("../NxCO2xI_d13c_air.csv")


###############################################################################
## Load propN functions
###############################################################################
source("../r_functions/calc_chi.R")

###############################################################################
## Focal leaf area
###############################################################################
## Calculate leaf area
ij.path <- "/Users/zinny/Desktop/ImageJ.app/"

## Calculate leaf disk area
imagepath.focal <- "/Users/zinny/git/C3C4_GrowthChamber_Experiment/leaf_area/"

focal.area <- run.ij(path.imagej = ij.path,
                     set.directory = imagepath.focal,
                     distance.pixel = 120.15, # done from set scale
                     known.distance = 1, low.size = 0.05,
                     set.memory = 30)

names(focal.area)[1:2] <- c("id", "focal.area")
focal.area$id <- gsub("_focal", "", focal.area$id)


#change ID tolowercase
focal.area <- focal.area%>%
  mutate(id = tolower(id)) %>%
  rename(ID = id)


#merge to isotope data
isotopes <-  isotopes %>%
  left_join(focal.area, by = "ID")


###############################################################################
## calculate nmass leaf N content
###############################################################################
##rename columns
totalN <- totalN %>% rename(ID = 'X.ID')

totalN$ID <- gsub("-", "_", totalN$ID)

isotopes <- isotopes %>%
  left_join(totalN %>% select(ID, total_N, sample_weight), by = "ID")


#total N = microgram, sample weight = mg
isotopes <- isotopes %>%
  mutate(nmass = as.numeric(total_N/sample_weight)/1000)



###############################################################################
## Calculate d13C in air between the two CO2 treatments
###############################################################################
d13c_air <- d13c.air %>%
  separate(id, into = c("name", "chamber", "co2", "rep")) %>%
  group_by(co2, chamber) %>%
  summarize(d13c.air = mean(d13c_air),
            co2.ppm = mean(co2_ppm)) %>%
  mutate(co2 = ifelse(co2.ppm < 800, 420, 1000)) %>%
  select(co2_cat = co2, everything(), -co2.ppm) %>%
  data.frame()

head(isotopes)

d13c_air$co2_cat <- as.double(d13c_air$co2_cat)


###

###############################################################################
## Calculate chi
###############################################################################
delta <- isotopes %>%
  mutate(chamber = str_extract(ID, "[^_]+$")) %>%
  dplyr::select(ID, chamber, leaf.d13c, co2, ps_pathway) %>%
  mutate(co2_cat = ifelse(co2 == "EC", 1000, 420)) %>%
  left_join(d13c_air, by = c("chamber", "co2_cat")) %>%
  mutate(delta = ifelse(ps_pathway == "C3", 
                  calc_chi_c3(leaf.d13c = as.numeric(leaf.d13c), 
                              air = as.numeric(d13c.air))[[1]],
                  calc_chi_c4(leaf.d13c = as.numeric(leaf.d13c), 
                              air = as.numeric(d13c.air))[[1]]))%>%
  dplyr::select(ID, delta)
  

isotopes <- isotopes %>% 
  left_join(delta,  by = "ID")



  
  ###############################################################################
  ## Calculate narea and marea
  ###############################################################################
  
  ## Leaf N content
isotopes <- isotopes %>% 
  mutate(marea = as.numeric(leaf.biomass)/ (focal.area / 10000)) %>% ###convert cm2 to m2
   mutate (narea = nmass * marea)

  
  
##merge file to bigdata
complieddata <- read.csv("../TXCO2.clean_complied_data.csv")

complieddata <- complieddata %>%
  left_join(isotopes %>% select(ID, ps_pathway, biomass.g., total_N, delta, 
                                nmass, marea, narea), by = c("ID", "ps_pathway"))

  

  
##write csv 
write.csv(complieddata, "TxCO2_fulldataset.csv")
  
  

