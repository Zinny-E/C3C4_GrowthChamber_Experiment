###############################################################################
## Libraries
###############################################################################
library(dplyr)
library(tidyverse)
library(plantecophys)

###############################################################################
## Import files
###############################################################################
chamber2 <- read.csv("~/git/C3C4_GrowthChamber_Experiment/licor_cleaned/TxCO2_combined_datasheets/TXCO2_co2_resp_chamber_2.csv")
##rd.chamber5 <- read.csv("../TXCO2_datasheets/TxCO2_rd_chamber5.csv")

aci_prep_chamber2 <- chamber2 %>%
  group_by(id) %>%
  dplyr::select(id, machine, A, Ci, Ca, gsw, 
                CO2_s,	CO2_r,	H2O_s,	H2O_r,
                Qin, VPDleaf, Flow,	Tair, Tleaf) %>%
  mutate(keep.row = "yes" ) %>%
  mutate(id = case_when(
    Tleaf >= 27.4 ~ paste0(id, "_27.5"),
    TRUE ~ as.character(id)
  )) %>%
  mutate(T_leaf =  mean(Tleaf, na.rm = TRUE))
data.frame()


## setting setpoint for temperature
## muatate(temp.setpoint = ifelse(tleaf>19 & tleaf <21, 20))
write_csv(aci_prep_chamber2, "../TxCO2_datasheets/TxCO2_co2_resp_chamber2.csv")

###############################################################################
## Load custom fxns
###############################################################################
source("/Users/zinny/git/r_functions/stomatal_limitation.R")

#####################################################################
# A/Ci curves
#####################################################################

#######################################
# Ambient CO2 low Temperature (C3)
#######################################
ely_can5_t3_ch2 <- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "ely_can5_t3_ch2" &
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can5_t3_ch2)
aci.coefs[48,] <- data.frame(id = "ely_can5_t3_ch2", pathway ="C3", t(coef(ely_can5_t3_ch2)))


ely_can5_t3_ch2_27.5 <- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "ely_can5_t3_ch2_27.5" &
                                                  Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can5_t3_ch2_27.5)
aci.coefs[49,] <- data.frame(id = "ely_can5_t3_ch2_27.5", pathway ="C3", t(coef(ely_can5_t3_ch2_27.5)))


pas_smi7_t4_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "pas_smi7_t4_ch2" &
                                                  Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi7_t4_ch2)
aci.coefs[50,] <- data.frame(id = "pas_smi7_t4_ch2", pathway ="C3", t(coef(pas_smi7_t4_ch2)))


pas_smi7_t4_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "pas_smi7_t4_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi7_t4_ch2_27.5)
aci.coefs[51,] <- data.frame(id = "pas_smi7_t4_ch2_27.5", pathway ="C3", t(coef(pas_smi7_t4_ch2_27.5)))


pas_smi1_t1_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "pas_smi1_t1_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi1_t1_ch2)
aci.coefs[52,] <- data.frame(id = "pas_smi1_t1_ch2", pathway ="C3", t(coef(pas_smi1_t1_ch2)))


pas_smi1_t1_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "pas_smi1_t1_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi1_t1_ch2_27.5)
aci.coefs[53,] <- data.frame(id = "pas_smi1_t1_ch2_27.5", pathway ="C3", t(coef(pas_smi1_t1_ch2_27.5)))


pas_smi4_t2_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "pas_smi4_t2_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi4_t2_ch2)
aci.coefs[54,] <- data.frame(id = "pas_smi4_t2_ch2", pathway ="C3", t(coef(pas_smi4_t2_ch2)))


pas_smi4_t2_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "pas_smi4_t2_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi4_t2_ch2_27.5)
aci.coefs[55,] <- data.frame(id = "pas_smi4_t2_ch2_27.5", pathway ="C3", t(coef(pas_smi4_t2_ch2_27.5)))


pas_smi6_t3_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "pas_smi6_t3_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi6_t3_ch2)
aci.coefs[56,] <- data.frame(id = "pas_smi6_t3_ch2", pathway ="C3", t(coef(pas_smi6_t3_ch2)))



pas_smi6_t3_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "pas_smi6_t3_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi6_t3_ch2_27.5)
aci.coefs[57,] <- data.frame(id = "pas_smi6_t3_ch2_27.5", pathway ="C3", t(coef(pas_smi6_t3_ch2_27.5)))


pas_smi8_t4_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "pas_smi8_t4_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi8_t4_ch2)
aci.coefs[58,] <- data.frame(id = "pas_smi8_t4_ch2", pathway ="C3", t(coef(pas_smi8_t4_ch2)))


pas_smi8_t4_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "pas_smi8_t4_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi8_t4_ch2_27.5)
aci.coefs[59,] <- data.frame(id = "pas_smi8_t4_ch2_27.5", pathway ="C3", t(coef(pas_smi8_t4_ch2_27.5)))


ely_can1_t1_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "ely_can1_t1_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can1_t1_ch2)
aci.coefs[60,] <- data.frame(id = "ely_can1_t1_ch2", pathway ="C3", t(coef(ely_can1_t1_ch2)))



ely_can1_t1_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "ely_can1_t1_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can1_t1_ch2_27.5)
aci.coefs[61,] <- data.frame(id = "ely_can1_t1_ch2_27.5", pathway ="C3", t(coef(ely_can1_t1_ch2_27.5)))



ely_can2_t1_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "ely_can2_t1_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can2_t1_ch2)
aci.coefs[62,] <- data.frame(id = "ely_can2_t1_ch2", pathway ="C3", t(coef(ely_can2_t1_ch2)))


ely_can2_t1_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "ely_can2_t1_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can2_t1_ch2_27.5)
aci.coefs[63,] <- data.frame(id = "ely_can2_t1_ch2_27.5", pathway ="C3", t(coef(ely_can2_t1_ch2_27.5)))


ely_can3_t2_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "ely_can3_t2_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can3_t2_ch2)
aci.coefs[64,] <- data.frame(id = "ely_can3_t2_ch2", pathway ="C3", t(coef(ely_can3_t2_ch2)))


ely_can3_t2_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "ely_can3_t2_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can3_t2_ch2_27.5)
aci.coefs[65,] <- data.frame(id = "ely_can3_t2_ch2_27.5", pathway ="C3", t(coef(ely_can3_t2_ch2_27.5)))


ely_can4_t2_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "ely_can4_t2_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can4_t2_ch2)
aci.coefs[66,] <- data.frame(id = "ely_can4_t2_ch2", pathway ="C3", t(coef(ely_can4_t2_ch2)))




ely_can4_t2_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "ely_can4_t2_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can4_t2_ch2_27.5)
aci.coefs[67,] <- data.frame(id = "ely_can4_t2_ch2_27.5", pathway ="C3", t(coef(ely_can4_t2_ch2_27.5)))


ely_can6_t3_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "ely_can6_t3_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can6_t3_ch2)
aci.coefs[68,] <- data.frame(id = "ely_can6_t3_ch2", pathway ="C3", t(coef(ely_can3_t2_ch2)))


ely_can6_t3_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "ely_can6_t3_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can6_t3_ch2_27.5)
aci.coefs[69,] <- data.frame(id = "ely_can6_t3_ch2_27.5", pathway ="C3", t(coef(ely_can3_t2_ch2_27.5)))



ely_can7_t4_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "ely_can7_t4_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can7_t4_ch2)
aci.coefs[70,] <- data.frame(id = "ely_can7_t4_ch2", pathway ="C3", t(coef(ely_can7_t4_ch2)))


ely_can7_t4_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "ely_can7_t4_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can7_t4_ch2_27.5)
aci.coefs[71,] <- data.frame(id = "ely_can7_t4_ch2_27.5", pathway ="C3", t(coef(ely_can7_t4_ch2_27.5)))



ely_can8_t4_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "ely_can8_t4_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can8_t4_ch2)
aci.coefs[72,] <- data.frame(id = "ely_can8_t4_ch2", pathway ="C3", t(coef(ely_can8_t4_ch2)))


ely_can8_t4_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "ely_can8_t4_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can8_t4_ch2_27.5)
aci.coefs[73,] <- data.frame(id = "ely_can8_t4_ch2_27.5", pathway ="C3", t(coef(ely_can8_t4_ch2_27.5)))


poa_pra1_t2_ch5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "poa_pra1_t2_ch5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra1_t2_ch5)
aci.coefs[74,] <- data.frame(id = "poa_pra1_t2_ch5", pathway ="C3", t(coef(poa_pra1_t2_ch5)))


poa_pra1_t2_ch5_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "poa_pra1_t2_ch5_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra1_t2_ch5_27.5)
aci.coefs[75,] <- data.frame(id = "poa_pra1_t2_ch5_27.5", pathway ="C3", t(coef(poa_pra1_t2_ch5_27.5)))


poa_pra2_t1_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "poa_pra2_t1_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra2_t1_ch2)
aci.coefs[76,] <- data.frame(id = "poa_pra2_t1_ch2", pathway ="C3", t(coef(poa_pra2_t1_ch2)))



poa_pra2_t1_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "poa_pra2_t1_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra2_t1_ch2_27.5)
aci.coefs[77,] <- data.frame(id = "poa_pra2_t1_ch2_27.5", pathway ="C3", t(coef(poa_pra2_t1_ch2_27.5)))


poa_pra3_t2_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "poa_pra3_t2_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra3_t2_ch2)
aci.coefs[78,] <- data.frame(id = "poa_pra3_t2_ch2", pathway ="C3", t(coef(poa_pra3_t2_ch2)))


poa_pra3_t2_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "poa_pra3_t2_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra3_t2_ch2_27.5)
aci.coefs[79,] <- data.frame(id = "poa_pra3_t2_ch2_27.5", pathway ="C3", t(coef(poa_pra3_t2_ch2_27.5)))


poa_pra4_t2_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "poa_pra4_t2_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra4_t2_ch2)
aci.coefs[80,] <- data.frame(id = "poa_pra4_t2_ch2", pathway ="C3", t(coef(poa_pra4_t2_ch2)))


poa_pra4_t2_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "poa_pra4_t2_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra4_t2_ch2_27.5)
aci.coefs[81,] <- data.frame(id = "poa_pra4_t2_ch2_27.5", pathway ="C3", t(coef(poa_pra4_t2_ch2_27.5)))


poa_pra5_t3_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "poa_pra5_t3_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra5_t3_ch2)
aci.coefs[82,] <- data.frame(id = "poa_pra5_t3_ch2", pathway ="C3", t(coef(poa_pra5_t3_ch2)))



poa_pra5_t3_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "poa_pra5_t3_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra5_t3_ch2_27.5)
aci.coefs[83,] <- data.frame(id = "poa_pra5_t3_ch2_27.5", pathway ="C3", t(coef(poa_pra5_t3_ch2_27.5)))


poa_pra6_t3_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "poa_pra6_t3_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra6_t3_ch2)
aci.coefs[84,] <- data.frame(id = "poa_pra6_t3_ch2", pathway ="C3", t(coef(poa_pra6_t3_ch2)))



poa_pra6_t3_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "poa_pra6_t3_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra6_t3_ch2_27.5)
aci.coefs[85,] <- data.frame(id = "poa_pra6_t3_ch2_27.5", pathway ="C3", t(coef(poa_pra6_t3_ch2_27.5)))


poa_pra8_t4_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "poa_pra8_t4_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra8_t4_ch2)
aci.coefs[86,] <- data.frame(id = "poa_pra8_t4_ch2", pathway ="C3", t(coef(poa_pra8_t4_ch2)))



poa_pra8_t4_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "poa_pra8_t4_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra8_t4_ch2_27.5)
aci.coefs[87,] <- data.frame(id = "poa_pra8_t4_ch2_27.5", pathway ="C3", t(coef(poa_pra8_t4_ch2_27.5)))



pas_smi5_t3_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "pas_smi5_t3_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi5_t3_ch2)
aci.coefs[88,] <- data.frame(id = "pas_smi5_t3_ch2", pathway ="C3", t(coef(pas_smi5_t3_ch2)))


pas_smi5_t3_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "pas_smi5_t3_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi5_t3_ch2_27.5)
aci.coefs[89,] <- data.frame(id = "pas_smi5_t3_ch2_27.5", pathway ="C3", t(coef(pas_smi5_t3_ch2_27.5)))


poa_pra7_t4_ch2<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "poa_pra7_t4_ch2" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra7_t4_ch2)
aci.coefs[90,] <- data.frame(id = "poa_pra7_t4_ch2", pathway ="C3", t(coef(poa_pra7_t4_ch2)))



poa_pra7_t4_ch2_27.5<- aci_prep_chamber2 %>% filter(keep.row == "yes" & id == "poa_pra7_t4_ch2_27.5" &
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_leaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd25"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra7_t4_ch2_27.5)
aci.coefs[91,] <- data.frame(id = "poa_pra7_t4_ch2_27.5", pathway ="C3", t(coef(poa_pra7_t4_ch2_27.5)))




