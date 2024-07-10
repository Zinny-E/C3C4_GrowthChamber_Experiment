###############################################################################
## Libraries
###############################################################################
library(dplyr)
library(tidyverse)
library(plantecophys)

###############################################################################
## Import files
###############################################################################
chamber4 <- read.csv("~/git/C3C4_GrowthChamber_Experiment/TxCO2_licorcleaned/TxCO2combinedataset/TXCO2_co2_resp_chamber_4.csv") %>%
  mutate(temp.setpoint = ifelse(Tleaf > 19 & Tleaf < 21, 
                                20,
                                ifelse(Tleaf > 24 & Tleaf < 29,
                                       27.5,
                                       ifelse(Tleaf > 30 & Tleaf < 36,
                                              35,
                                              NA))))



rd_chamber4 <- read.csv("~/git/C3C4_GrowthChamber_Experiment/TxCO2_licorcleaned/TxCO2combinedataset/TxCO2_rd_chamber_4.csv") %>%
  mutate(temp.setpoint = ifelse(Tleaf > 19 & Tleaf < 21, 
                                20,
                                ifelse(Tleaf > 24 & Tleaf < 29,
                                       27.5,
                                       ifelse(Tleaf > 30 & Tleaf < 36,
                                              35,
                                              NA))))


##grouping respiration by id and temp setpoint and then finding the mean
rd_chamber4_group_by <- group_by(rd_chamber4, id, temp.setpoint, chamber, machine)
rd_chamber4_mean <- summarise(rd_chamber4_group_by, rd_mean = mean(rd, na.rm = T))
head(rd_chamber4_mean)


###############################################################################
## Load custom fxns
###############################################################################
source("/Users/zinny/git/r_functions/stomatal_limitation.R")

###############################################################################
## Prep data frame to fit A/Ci curves
###############################################################################
aci.prep.chamber4 <- chamber4  %>%
  group_by(id) %>%
  dplyr::select(id, machine, A, Ci, Ca, gsw, 
                CO2_s,	CO2_r,	H2O_s,	H2O_r,
                Qin, VPDleaf, Flow,	Tair, Tleaf, temp.setpoint) %>%
  arrange(id) %>%
  left_join(rd_chamber4_mean, by = c("id", "temp.setpoint")) %>%
  filter(Qin > 1000) %>% # remove strange Qin values - seems to be Rd values merged into compiled sheet
  dplyr::select (-Tleaf) %>%
  group_by(id) %>%
  mutate(keep.row = "yes") %>%
  mutate(id = case_when(
    temp.setpoint <= 27.5 ~ paste0(id, "_27.5"), #adding _27.5 to the id's that were measured at 27.5
    TRUE ~ as.character(id)
  )) %>%
  data.frame()


#aci.prep.chamber4<- rename(aci.prep.chamber4, Tleaf = temp.setpoint, rd = rd_mean)

write.csv(aci.prep.chamber4, "../C3C4_GrowthChamber_Experiment/TxCO2_licorcleaned/TxCO2combinedataset/aci.prep.chamber4.csv", 
          row.names = FALSE)


#######################################
# Chamber 4 high CO2 High Temperature
#######################################

ely_can19_t2_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "ely_can19_t2_ch4_27.5" &
                                                       Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can19_t2_ch4_27.5)
#summary(ely_can19_t2_ch4_27.5)
chamber4_aci.coefs <- data.frame(id = "ely_can19_t2_ch4_27.5", pathway = "C3", t(coef(ely_can19_t2_ch4_27.5)))


ely_can19_t2_ch4 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "ely_can19_t2_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can19_t2_ch4)
#summary(ely_can19_t2_ch4)
chamber4_aci.coefs[2,] <- c(id = "ely_can19_t2_ch4", pathway = "C3", t(coef(ely_can19_t2_ch4)))

ely_can20_t2_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "ely_can20_t2_ch4_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can20_t2_ch4_27.5)
#summary(ely_can20_t2_ch4_27.5)
chamber4_aci.coefs[3,] <- c(id = "ely_can20_t2_ch4_27.5", pathway = "C3", t(coef(ely_can20_t2_ch4_27.5)))

ely_can20_t2_ch4 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "ely_can20_t2_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

aci.prep.chamber4$keep.row[c(2231,2232)] <- "no" # removes points that are outliers

plot(ely_can20_t2_ch4)
#summary(ely_can20_t2_ch4)
chamber4_aci.coefs[4,] <- c(id = "ely_can20_t2_ch4", pathway = "C3", t(coef(ely_can20_t2_ch4)))


ely_can21_t3_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "ely_can21_t3_ch4_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can21_t3_ch4_27.5)
#summary(ely_can21_t3_ch4_27.5)
chamber4_aci.coefs[5,] <- c(id = "ely_can21_t3_ch4_27.5", pathway = "C3", t(coef(ely_can21_t3_ch4_27.5)))


ely_can21_t3_ch4<- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "ely_can21_t3_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can21_t3_ch4)
#summary(ely_can20_t2_ch4_27.5)
chamber4_aci.coefs[6,] <- c(id = "ely_can21_t3_ch4", pathway = "C3", t(coef(ely_can21_t3_ch4)))


ely_can22_t3_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "ely_can22_t3_ch4_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


##aci.prep.chamber4$keep.row[c(2228,2229,2230,2231,2232)] <- "no" # removes points that are outliers

plot(ely_can22_t3_ch4_27.5)
#summary(ely_can22_t3_ch4_27.5)
chamber4_aci.coefs[7,] <- c(id = "ely_can22_t3_ch4_27.5", pathway = "C3", t(coef(ely_can22_t3_ch4_27.5)))

ely_can22_t3_ch4 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "ely_can22_t3_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)




plot(ely_can22_t3_ch4)
#summary(ely_can22_t3_ch4)
chamber4_aci.coefs[8,] <- c(id = "ely_can22_t3_ch4", pathway = "C3", t(coef(ely_can22_t3_ch4)))


ely_can23_t4_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "ely_can23_t4_ch4_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can23_t4_ch4_27.5)
#summary(ely_can23_t4_ch4_27.5)
chamber4_aci.coefs[9,] <- c(id = "ely_can23_t4_ch4_27.5", pathway = "C3", t(coef(ely_can23_t4_ch4_27.5)))


ely_can23_t4_ch4 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "ely_can23_t4_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

##aci.prep.chamber4$keep.row[c(2805,2806,2807,2808)] <- "no" # removes points that are outliers

plot(ely_can23_t4_ch4)
#summary(ely_can23_t4_ch4)
chamber4_aci.coefs[10,] <- c(id = "ely_can23_t4_ch4", pathway = "C3", t(coef(ely_can23_t4_ch4)))


ely_can24_t4_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "ely_can24_t4_ch4_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

##aci.prep.chamber4$keep.row[c(2894,2895,2896,2897,2898,2899,2900,2901,2902,2903,
 #                             2904)] <- "no" # removes points that are outliers

plot(ely_can24_t4_ch4_27.5)
#summary(ely_can24_t4_ch4_27.5)
chamber4_aci.coefs[11,] <- c(id = "ely_can24_t4_ch4_27.5", pathway = "C3", t(coef(ely_can24_t4_ch4_27.5)))


ely_can24_t4_ch4<- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "ely_can24_t4_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can24_t4_ch4)
#summary(ely_can24_t4_ch4)
chamber4_aci.coefs[12,] <- c(id = "ely_can24_t4_ch4", pathway = "C3", t(coef(ely_can24_t4_ch4)))



pas_smi18_t1_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "Pas_smi18_t1_ch4_27.5" &
                                                  Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

aci.prep.chamber4$keep.row[c(5,6,7,8,9,10,11,12,13,14,
                             15,16,17,18,19,20,21,22,23,24,
                             25,26,27,28,29,30,31,32,33,34,
                             35,36,37,38,39,40,41,42,43,44,
                             45,46,47,48,49,50,51,52,53,54,
                             55,56,57,58,59,60,61,62,63,64,
                             65,66,67,68,69,70,71,72,73,74,
                             75,76,77,78,79,80,81,82,83,84,
                             85,86,87,88,89,90,91,92,93,94,
                             95,96)] <- "no" # removes points that are outliers

plot(pas_smi18_t1_ch4_27.5)
chamber4_aci.coefs[13,] <- c(id = "pas_smi18_t1_ch4_27.5", pathway = "C3", t(coef(pas_smi18_t1_ch4_27.5)))


pas_smi18_t1_ch4<- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "Pas_smi18_t1_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)



plot(pas_smi18_t1_ch4_27.5)
chamber4_aci.coefs[14,] <- c(id = "pas_smi18_t1_ch4_27.5", pathway = "C3", t(coef(pas_smi18_t1_ch4_27.5)))


pas_smi20_t2_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "pas_smi20_t2_ch4_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi20_t2_ch4_27.5)
chamber4_aci.coefs[15,] <- c(id = "pas_smi20_t2_ch4_27.5", pathway = "C3", t(coef(pas_smi20_t2_ch4_27.5)))


pas_smi20_t2_ch4 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "pas_smi20_t2_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi20_t2_ch4)
chamber4_aci.coefs[16,] <- c(id = "pas_smi20_t2_ch4", pathway = "C3", t(coef(pas_smi20_t2_ch4)))


pas_smi21_t3_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "pas_smi21_t3_ch4_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi21_t3_ch4_27.5)
chamber4_aci.coefs[17,] <- c(id = "pas_smi21_t3_ch4_27.5", pathway = "C3", t(coef(pas_smi21_t3_ch4_27.5)))


pas_smi21_t3_ch4 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "pas_smi21_t3_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi21_t3_ch4)
chamber4_aci.coefs[18,] <- c(id = "pas_smi21_t3_ch4", pathway = "C3", t(coef(pas_smi21_t3_ch4)))



pas_smi22_t3_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "pas_smi22_t3_ch4_27.5" &
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi22_t3_ch4_27.5)
chamber4_aci.coefs[19,] <- c(id = "pas_smi22_t3_ch4_27.5", pathway = "C3", t(coef(pas_smi22_t3_ch4_27.5)))



pas_smi22_t3_ch4 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "pas_smi22_t3_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi22_t3_ch4)
chamber4_aci.coefs[20,] <- c(id = "pas_smi22_t3_ch4", pathway = "C3", t(coef(pas_smi22_t3_ch4)))


pas_smi23_t4_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "pas_smi23_t4_ch4_27.5" &
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi23_t4_ch4_27.5)
chamber4_aci.coefs[21,] <- c(id = "pas_smi23_t4_ch4_27.5", pathway = "C3", t(coef(pas_smi23_t4_ch4_27.5)))


pas_smi23_t4_ch4 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "pas_smi23_t4_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi23_t4_ch4)
chamber4_aci.coefs[22,] <- c(id = "pas_smi23_t4_ch4", pathway = "C3", t(coef(pas_smi23_t4_ch4)))



pas_smi24_t4_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "pas_smi24_t4_ch4_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi24_t4_ch4_27.5)
chamber4_aci.coefs[23,] <- c(id = "pas_smi24_t4_ch4_27.5", pathway = "C3", t(coef(pas_smi24_t4_ch4_27.5)))



pas_smi24_t4_ch4 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "pas_smi24_t4_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi24_t4_ch4)
chamber4_aci.coefs[24,] <- c(id = "pas_smi24_t4_ch4", pathway = "C3", t(coef(pas_smi24_t4_ch4)))

poa_pra17_t1_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "poa_pra17_t1_ch4_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra17_t1_ch4_27.5)
chamber4_aci.coefs[25,] <- c(id = "poa_pra17_t1_ch4_27.5", pathway = "C3", t(coef(poa_pra17_t1_ch4_27.5)))

poa_pra17_t1_ch4 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "poa_pra17_t1_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra17_t1_ch4)
chamber4_aci.coefs[26,] <- c(id = "poa_pra17_t1_ch4", pathway = "C3", t(coef(poa_pra17_t1_ch4)))


poa_pra18_t1_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "poa_pra18_t1_ch4_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra18_t1_ch4_27.5)
chamber4_aci.coefs[27,] <- c(id = "poa_pra18_t1_ch4_27.5", pathway = "C3", t(coef(poa_pra18_t1_ch4_27.5)))




poa_pra18_t1_ch4 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "poa_pra18_t1_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra18_t1_ch4)
chamber4_aci.coefs[28,] <- c(id = "poa_pra118_t1_ch4", pathway = "C3", t(coef(poa_pra17_t1_ch4)))


poa_pra19_t2_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "poa_pra19_t2_ch4_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra19_t2_ch4_27.5)
chamber4_aci.coefs[29,] <- c(id = "poa_pra19_t2_ch4_27.5", pathway = "C3", t(coef(poa_pra19_t2_ch4_27.5)))


poa_pra19_t2_ch4<- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "poa_pra19_t2_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra19_t2_ch4)
chamber4_aci.coefs[30,] <- c(id = "poa_pra19_t2_ch4", pathway = "C3", t(coef(poa_pra19_t2_ch4)))


poa_pra20_t2_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "poa_pra20_t2_ch4_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra20_t2_ch4_27.5)
chamber4_aci.coefs[31,] <- c(id = "poa_pra20_t2_ch4_27.5", pathway = "C3", t(coef(poa_pra20_t2_ch4_27.5)))



poa_pra20_t2_ch4 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "poa_pra20_t2_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra20_t2_ch4)
chamber4_aci.coefs[32,] <- c(id = "poa_pra20_t2_ch4", pathway = "C3", t(coef(poa_pra20_t2_ch4)))

poa_pra21_t3_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "poa_pra21_t3_ch4_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra21_t3_ch4_27.5)
chamber4_aci.coefs[33,] <- c(id = "poa_pra21_t3_ch4_27.5", pathway = "C3", t(coef(poa_pra21_t3_ch4_27.5)))


poa_pra21_t3_ch4<- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "poa_pra21_t3_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra21_t3_ch4)
chamber4_aci.coefs[34,] <- c(id = "poa_pra21_t3_ch4", pathway = "C3", t(coef(poa_pra21_t3_ch4)))


poa_pra23_t4_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "poa_pra23_t4_ch4_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra23_t4_ch4_27.5)
chamber4_aci.coefs[35,] <- c(id = "poa_pra23_t4_ch4_27.5", pathway = "C3", t(coef(poa_pra23_t4_ch4_27.5)))



poa_pra23_t4_ch4<- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "poa_pra23_t4_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra23_t4_ch4)
chamber4_aci.coefs[36,] <- c(id = "poa_pra23_t4_ch4", pathway = "C3", t(coef(poa_pra23_t4_ch4)))


poa_pra24_t4_ch4_27.5 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "poa_pra24_t4_ch4_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra24_t4_ch4_27.5)
chamber4_aci.coefs[37,] <- c(id = "poa_pra24_t4_ch4_27.5", pathway = "C3", t(coef(poa_pra24_t4_ch4_27.5)))


poa_pra24_t4_ch4 <- aci.prep.chamber4 %>% filter(keep.row == "yes" & id == "poa_pra24_t4_ch4" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra24_t4_ch4)
chamber4_aci.coefs[38,] <- c(id = "poa_pra24_t4_ch4", pathway = "C3", t(coef(poa_pra24_t4_ch4)))


##write aci_coefs data

write.csv(chamber4_aci.coefs, "~/git/C3C4_GrowthChamber_Experiment/TxCO2_licorcleaned/TxCO2combinedataset/aci.coefs_chamber4.csv", 
          row.names = FALSE)


