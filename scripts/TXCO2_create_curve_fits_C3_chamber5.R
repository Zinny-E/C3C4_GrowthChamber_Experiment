###############################################################################
## Libraries
###############################################################################
library(dplyr)
library(tidyverse)
library(plantecophys)

###############################################################################
## Import files
###############################################################################
chamber5 <- read.csv("../TXCO2_datasheets/TxCO2_resp_chamber5.csv")
rd.chamber5 <- read.csv("../TXCO2_datasheets/TxCO2_rd_chamber5.csv")

aci_prep <- chamber5 %>%
  group_by(id) %>%
  dplyr::select(id, machine, A, Ci, Ca, gsw, 
                CO2_s,	CO2_r,	H2O_s,	H2O_r,
                Qin, VPDleaf, Flow,	Tair, Tleaf,	tleaf) %>%
  mutate(keep.row = "yes") %>%
    mutate(id = case_when(
      tleaf >= 27.5 ~ paste0(id, "_27.5"),
      TRUE ~ as.character(id)
    ))
  data.frame()

#write_csv(aci_prep, "../TXCO2_datasheets/chamber5.csv")

###############################################################################
## Load custom fxns
###############################################################################
source("/Users/zinny/git/r_functions/stomatal_limitation.R")
  
#####################################################################
# A/Ci curves
#####################################################################
  
#######################################
# Elevated CO2 Low Tempertaure (C3)
#######################################
pas_smi25_t1_ch5 <- aci_prep %>% filter(keep.row == "yes" & id == "pas_smi25_t1_ch5" &
                                          Ci < 1200) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = 19.9),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi25_t1_ch5)

aci.coefs <- data.frame(id = "pas_smi25_t1_ch5", pathway ="C3", t(coef(pas_smi25_t1_ch5)))



pas_smi25_t1_ch5_27.5 <- aci_prep %>% filter(keep.row == "yes" & id == "pas_smi25_t1_ch5_27.5" &
                                           Ci < 1200) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = 27.5),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi25_t1_ch5_27.5)
aci.coefs[2,] <- data.frame(id = "pas_smi25_t1_ch5_27.5", pathway ="C3", t(coef(pas_smi25_t1_ch5_27.5)))


pas_smi26_t1_ch5<- aci_prep %>% filter(keep.row == "yes" & id == "pas_smi26_t1_ch5" &
                                                Ci < 1200) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = 19.9),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi26_t1_ch5)
aci.coefs[3,] <- data.frame(id = "pas_smi26_t1_ch5", pathway ="C3", t(coef(pas_smi26_t1_ch5)))



pas_smi28_t2_ch5_redo<- aci_prep %>% filter(keep.row == "yes" & id == "pas_smi28_t2_ch5_redo" &
                                         Ci < 1200) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = 19.9),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi28_t2_ch5_redo)
aci.coefs[11,] <- data.frame(id = "pas_smi28_t2_ch5_redo", pathway ="C3", t(coef(pas_smi28_t2_ch5_redo)))

ely_can25_t1_ch5 <- aci_prep %>% filter(keep.row == "yes" & id == "ely_can25_t1_ch5" &
                                                Ci < 1200) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = 20),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can25_t1_ch5)
aci.coefs[4,] <- data.frame(id = "ely_can25_t1_ch5", pathway ="C3", t(coef(ely_can25_t1_ch5)))


ely_can25_t1_ch5_27.5 <- aci_prep %>% filter(keep.row == "yes" & id == "ely_can25_t1_ch5_27.5" &
                                           Ci < 850) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = 27.5),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can25_t1_ch5_27.5)
aci.coefs[5,] <- data.frame(id = "ely_can25_t1_ch5_27.5", pathway ="C3", t(coef(ely_can25_t1_ch5_27.5)))



poa_pra25_t1_ch5 <- aci_prep %>% filter(keep.row == "yes" & id == "poa_pra25_t1_ch5" &
                                           Ci < 850) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = 20),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra25_t1_ch5)
aci.coefs[6,] <- data.frame(id = "poa_pra25_t1_ch5", pathway ="C3", t(coef(poa_pra25_t1_ch5)))



poa_pra25_t1_ch5_27.5 <- aci_prep %>% filter(keep.row == "yes" & id == "poa_pra25_t1_ch5_27.5" &
                                           Ci < 1200) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = 27.5),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra25_t1_ch5_27.5)
aci.coefs[7,] <- data.frame(id = "poa_pra25_t1_ch5_27.5", pathway ="C3", t(coef(poa_pra25_t1_ch5_27.5)))


poa_pra25_t1_ch5 <- aci_prep %>% filter(keep.row == "yes" & id == "poa_pra25_t1_ch5" &
                                           Ci < 1200) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = 20),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra25_t1_ch5)
aci.coefs[6,] <- data.frame(id = "poa_pra25_t1_ch5", pathway ="C3", t(coef(poa_pra25_t1_ch5)))


poa_pra30_t3_ch5 <- aci_prep %>% filter(keep.row == "yes" & id == "poa_pra30_t3_ch5" &
                                           Ci < 1200) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = 20),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra30_t3_ch5)
aci.coefs[7,] <- data.frame(id = "poa_pra30_t3_ch5", pathway ="C3", t(coef(poa_pra30_t3_ch5)))



poa_pra30_t3_ch5_27.5 <- aci_prep %>% filter(keep.row == "yes" & id == "poa_pra30_t3_ch5_27.5" &
                                           Ci < 1200) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = 20),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra30_t3_ch5_27.5)
aci.coefs[8,] <- data.frame(id = "poa_pra30_t3_ch5_27.5", pathway ="C3", t(coef(poa_pra30_t3_ch5_27.5)))


poa_pra31_t4_ch5 <- aci_prep %>% filter(keep.row == "yes" & id == "poa_pra31_t4_ch5" &
                                           Ci < 1200) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = 20),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra31_t4_ch5)
aci.coefs[9,] <- data.frame(id = "poa_pra31_t4_ch5", pathway ="C3", t(coef(poa_pra31_t4_ch5)))


poa_pra31_t4_ch5_27.5 <- aci_prep%>% filter(keep.row == "yes" & id == "poa_pra31_t4_ch5_27.5" &
                                           Ci < 1200) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = 27.5),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra31_t4_ch5_27.5)
aci.coefs[10,] <- data.frame(id = "poa_pra31_t4_ch5_27.5", pathway ="C3", t(coef(poa_pra31_t4_ch5_27.5)))


