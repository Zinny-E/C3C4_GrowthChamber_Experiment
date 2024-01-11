###############################################################################
## Libraries
###############################################################################
library(dplyr)
library(tidyverse)
library(plantecophys)

###############################################################################
## Import files
###############################################################################
chamber2 <- read.csv("~/git/C3C4_GrowthChamber_Experiment/licor_cleaned/TxCO2_combined_datasheets/TXCO2_co2_resp_chamber_2.csv") %>%
    mutate(temp.setpoint = ifelse(Tleaf > 19 & Tleaf < 21, 
                                  20,
                                  ifelse(Tleaf > 27 & Tleaf < 28,
                                         27.5,
                                         ifelse(Tleaf > 34 & Tleaf < 36,
                                                35,
                                                NA))))



rd_chamber2 <- read.csv("~/git/C3C4_GrowthChamber_Experiment/licor_cleaned/TxCO2_combined_datasheets/TxCO2_rd_chamber_2.csv") %>%
    mutate(temp.setpoint = ifelse(Tleaf > 19 & Tleaf < 21, 
                                  20,
                                  ifelse(Tleaf > 27 & Tleaf < 28,
                                         27.5,
                                         ifelse(Tleaf > 34 & Tleaf < 36,
                                                35,
                                                NA))))
##grouping respiration by id and temp setpoint and then finding the mean
rd_chamber2_group_by <- group_by(rd_chamber2, id, temp.setpoint, chamber, machine)
rd_chamber2_mean <- summarise(rd_chamber2_group_by, rd_mean = mean(rd, na.rm = T))
head(rd_chamber2_mean)



###############################################################################
## Load custom fxns
###############################################################################
source("/Users/zinny/git/r_functions/stomatal_limitation.R")

###############################################################################
## Prep data frame to fit A/Ci curves
###############################################################################
aci.prep <- chamber2  %>%
  group_by(id) %>%
  dplyr::select(id, machine, A, Ci, Ca, gsw, 
                CO2_s,	CO2_r,	H2O_s,	H2O_r,
                Qin, VPDleaf, Flow,	Tair, Tleaf, temp.setpoint) %>%
  arrange(id) %>%
  left_join(rd_chamber2_mean, by = c("id", "temp.setpoint")) %>%
  filter(Qin > 1000) %>% # remove strange Qin values - seems to be Rd values merged into compiled sheet
  dplyr::select (-Tleaf) %>%
  group_by(id) %>%
  mutate(keep.row = "yes") %>%
  mutate(id = case_when(
    temp.setpoint >= 27.5 ~ paste0(id, "_27.5"), #addinf _27.5 to the id's that were measured at 27.5
    TRUE ~ as.character(id)
  )) %>%
  data.frame()


aci.prep<- rename(aci.prep, Tleaf = temp.setpoint, rd = rd_mean)

##
b <- chamber2_aci.coefs
a <- ifelse(grepl("_27.5", chamber2_aci.coefs$id), 27.5, 20)
b$temperature <- a
b <- na.omit(b)
d <- b %>% group_by(temperature) %>%
  summarise(meanVcmax=mean(as.numeric(Vcmax)))

class(b$Vcmax)
ggplot(d, aes(temperature, Vcmax)) +
  geom_point()
#aci.prep <- left_join(chamber2, rd_chamber2_mean, by = c("id", "temp.setpoint"))
#head(aci.prep)

#aci.prep$keep.row[c(,)] <- "no" # removes points that are outliers

##write.csv(aci.prep, "../licor_cleaned/TxCO2_combined_datasheets/aci.prep.chamber2.csv", row.names = FALSE)
#####################################################################
# A/Ci curves
#####################################################################

#######################################
# Chamber 2 Ambient CO2 Low Tempertaure
#######################################

ely_can1_t1_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "ely_can1_t1_ch2" &
                                                  Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
#summary(ely_can1_t1_ch2)
plot(ely_can1_t1_ch2)
chamber2_aci.coefs <- data.frame(id = "ely_can1_t1_ch2", pathway = "C3", t(coef(ely_can1_t1_ch2)))


ely_can1_t1_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "ely_can1_t1_ch2_27.5" &
                                              Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can1_t1_ch2_27.5)
chamber2_aci.coefs[2,] <- c(id = "ely_can1_t1_ch2_27.5", pathway = "C3", t(coef(ely_can1_t1_ch2_27.5)))

#no curve made, make curve..
ely_can2_t1_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "ely_can2_t1_ch2" &
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
summary(ely_can2_t1_ch2)
#aci.prep$keep.row[c(1438,1440,1479,1489,1487,1502)] <- "no" # removes points that are outliers
plot(ely_can2_t1_ch2)
chamber2_aci.coefs[3,] <- data.frame(id = "ely_can2_t1_ch2", pathway = "C3",t(coef(ely_can2_t1_ch2)))


ely_can2_t1_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "ely_can2_t1_ch2_27.5" &
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
summary(ely_can2_t1_ch2_27.5)
aci.prep$keep.row[c(1580,1581,1582, 1583, 1584,1585,1586,1587,1588,1589,1590,1591,1576,1577,1536,1579,1578, 1567,1535,1562)] <- "no" # removes points that are outliers
plot(ely_can2_t1_ch2_27.5)
chamber2_aci.coefs[4,] <- data.frame(id = "ely_can2_t1_ch2_27.5",pathway = "C3", t(coef(ely_can2_t1_ch2_27.5)))


ely_can3_t2_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "ely_can3_t2_ch2" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can3_t2_ch2)
chamber2_aci.coefs[5,] <- c(id = "ely_can3_t2_ch2", pathway = "C3",t(coef(ely_can3_t2_ch2)))


ely_can3_t2_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "ely_can3_t2_ch2_27.5" & 
                                              Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can3_t2_ch2_27.5)
chamber2_aci.coefs[6,] <- c(id = "ely_can3_t2_ch2_27.5",pathway = "C3", t(coef(ely_can3_t2_ch2_27.5)))


ely_can4_t2_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "ely_can4_t2_ch2"& 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can4_t2_ch2)
chamber2_aci.coefs[7,] <- c(id = "ely_can4_t2_ch2",pathway = "C3", t(coef(ely_can4_t2_ch2)))

#redo
ely_can4_t2_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "ely_can4_t2_ch2_27.5"& 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
#aci.prep$keep.row[c(1958)] <- "no" # removes points that are outliers
plot(ely_can4_t2_ch2_27.5)
chamber2_aci.coefs[8,] <- c(id = "ely_can4_t2_ch2_27.5", t(coef(ely_can4_t2_ch2_27.5)))

ely_can5_t3_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "ely_can5_t3_ch2"& 
                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
aci.prep$keep.row[c(2016)] <- "no" # removes points that are outliers
plot(ely_can5_t3_ch2)
chamber2_aci.coefs[9,] <- c(id = "ely_can5_t3_ch2", t(coef(ely_can5_t3_ch2)))



ely_can5_t3_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "ely_can5_t3_ch2_27.5"& 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
#aci.prep$keep.row[c(2110,2111,2112,2113.2114)] <- "no" # removes points that are outliers
plot(ely_can5_t3_ch2_27.5)
chamber2_aci.coefs[10,] <- c(id = "ely_can5_t3_ch2_27.5", pathway = "C3",t(coef(ely_can5_t3_ch2_27.5)))


ely_can6_t3_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "ely_can6_t3_ch2" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can6_t3_ch2)
chamber2_aci.coefs[11,] <- c(id = "ely_can6_t3_ch2", pathway = "C3",t(coef(ely_can6_t3_ch2)))


ely_can6_t3_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "ely_can6_t3_ch2_27.5"& 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can6_t3_ch2_27.5)
chamber2_aci.coefs[12,] <- c(id = "ely_can6_t3_ch2_27.5", pathway = "C3",t(coef(ely_can6_t3_ch2_27.5)))


ely_can7_t4_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "ely_can7_t4_ch2"& 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
#aci.prep$keep.row[c(2464,2447, 2467,2470, 2472,2450, 2452,2468,2466,2463,2465,2471, 2446,2448,2462,2448,2449,2446,2461,2460)] <- "no" # removes points that are outliers
plot(ely_can7_t4_ch2)
chamber2_aci.coefs[13,] <- c(id = "ely_can7_t4_ch2", pathway = "C3",t(coef(ely_can7_t4_ch2)))


ely_can7_t4_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "ely_can7_t4_ch2_27.5"& 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can7_t4_ch2_27.5)
chamber2_aci.coefs[14,] <- c(id = "ely_can7_t4_ch2_27.5",pathway = "C3", t(coef(ely_can7_t4_ch2_27.5)))


ely_can8_t4_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "ely_can8_t4_ch2"& 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can8_t4_ch2)
chamber2_aci.coefs[15,] <- c(id = "ely_can8_t4_ch2", pathway = "C3",t(coef(ely_can8_t4_ch2)))



ely_can8_t4_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "ely_can8_t4_ch2_27.5"& 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can8_t4_ch2_27.5)
chamber2_aci.coefs[16,] <- c(id = "ely_can8_t4_ch2_27.5", pathway = "C3",t(coef(ely_can8_t4_ch2_27.5)))



pas_smi1_t1_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "pas_smi1_t1_ch2"& 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi1_t1_ch2)
chamber2_aci.coefs[17,] <- c(id = "pas_smi1_t1_ch2", pathway = "C3",t(coef(pas_smi1_t1_ch2)))



pas_smi1_t1_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "pas_smi1_t1_ch2_27.5"& 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi1_t1_ch2_27.5)
chamber2_aci.coefs[18,] <- c(id = "pas_smi1_t1_ch2_27.5", pathway = "C3",t(coef(pas_smi1_t1_ch2_27.5)))


pas_smi4_t2_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "pas_smi4_t2_ch2" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi4_t2_ch2)
chamber2_aci.coefs[19,] <- c(id = "pas_smi4_t2_ch2", pathway = "C3",t(coef(pas_smi4_t2_ch2)))


pas_smi4_t2_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "pas_smi4_t2_ch2_27.5" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi4_t2_ch2_27.5)
chamber2_aci.coefs[20,] <- c(id = "pas_smi4_t2_ch2_27.5", pathway = "C3",t(coef(pas_smi4_t2_ch2_27.5)))


pas_smi5_t3_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "pas_smi5_t3_ch2" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi5_t3_ch2)
chamber2_aci.coefs[21,] <- c(id = "pas_smi5_t3_ch2",pathway = "C3", t(coef(pas_smi5_t3_ch2)))


pas_smi5_t3_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "pas_smi5_t3_ch2_27.5" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi5_t3_ch2_27.5)
chamber2_aci.coefs[22,] <- c(id = "pas_smi5_t3_ch2_27.5", pathway = "C3",t(coef(pas_smi5_t3_ch2_27.5)))


pas_smi6_t3_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "pas_smi6_t3_ch2" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi6_t3_ch2)
chamber2_aci.coefs[23,] <- c(id = "pas_smi6_t3_ch2", pathway = "C3",t(coef(pas_smi6_t3_ch2)))


pas_smi6_t3_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "pas_smi6_t3_ch2_27.5" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi6_t3_ch2_27.5)
chamber2_aci.coefs[24,] <- c(id = "pas_smi6_t3_ch2_27.5",pathway = "C3", t(coef(pas_smi6_t3_ch2_27.5)))


pas_smi7_t4_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "pas_smi7_t4_ch2" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi7_t4_ch2)
chamber2_aci.coefs[25,] <- c(id = "pas_smi7_t4_ch2", pathway = "C3",t(coef(pas_smi7_t4_ch2)))


pas_smi7_t4_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "pas_smi7_t4_ch2_27.5" & 
                                         Ci < 850) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
#aci.prep$keep.row[c(3649,3648,3647,3646)] <- "no" # removes points that are outliers

plot(pas_smi7_t4_ch2_27.5)
chamber2_aci.coefs[26,] <- c(id = "pas_smi7_t4_ch2_27.5",pathway = "C3", t(coef(pas_smi7_t4_ch2_27.5)))


pas_smi8_t4_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "pas_smi8_t4_ch2" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi8_t4_ch2)
chamber2_aci.coefs[27,] <- c(id = "pas_smi8_t4_ch2", pathway = "C3",t(coef(pas_smi8_t4_ch2)))


pas_smi8_t4_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "pas_smi8_t4_ch2_27.5" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi8_t4_ch2_27.5)
chamber2_aci.coefs[28,] <- c(id = "pas_smi8_t4_ch2_27.5",pathway = "C3", t(coef(pas_smi8_t4_ch2_27.5)))

poa_pra1_t1_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "poa_pra1_t1_ch2" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

#aci.prep$keep.row[c(3937, 3935,3936)] <- "no" # removes points that are outliers
plot(poa_pra1_t1_ch2)
chamber2_aci.coefs[29,] <- c(id = "poa_pra1_t1_ch2",pathway = "C3", t(coef(poa_pra1_t1_ch2)))

poa_pra1_t1_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "poa_pra1_t1_ch2_27.5" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

#aci.prep$keep.row[c(4033, 4032, 4031)] <- "no" # removes points that are outliers
plot(poa_pra1_t1_ch2_27.5)
chamber2_aci.coefs[30,] <- c(id = "poa_pra1_t1_ch2_27.5", pathway = "C3",t(coef(poa_pra1_t1_ch2_27.5)))


poa_pra2_t1_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "poa_pra2_t1_ch2" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

#aci.prep$keep.row[c(4128,4126,4129,4127,4125)] <- "no" # removes points that are outliers
plot(poa_pra2_t1_ch2)
chamber2_aci.coefs[31,] <- c(id = "poa_pra2_t1_ch2", pathway = "C3",t(coef(poa_pra2_t1_ch2)))


poa_pra2_t1_ch2_27.5<- aci.prep %>% filter(keep.row == "yes" & id == "poa_pra2_t1_ch2_27.5" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

#aci.prep$keep.row[c(4225,4224)] <- "no" # removes points that are outliers
plot(poa_pra2_t1_ch2_27.5)
chamber2_aci.coefs[32,] <- c(id = "poa_pra2_t1_ch2_27.5", pathway = "C3",t(coef(poa_pra2_t1_ch2_27.5)))

poa_pra3_t2_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "poa_pra3_t2_ch2" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra3_t2_ch2)
chamber2_aci.coefs[33,] <- c(id = "poa_pra3_t2_ch2", pathway = "C3",t(coef(poa_pra3_t2_ch2)))

poa_pra3_t2_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "poa_pra3_t2_ch2_27.5" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra3_t2_ch2_27.5)
chamber2_aci.coefs[34,] <- c(id = "poa_pra3_t2_ch2_27.5", pathway = "C3",t(coef(poa_pra3_t2_ch2_27.5)))

poa_pra4_t2_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "poa_pra4_t2_ch2" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra4_t2_ch2)
chamber2_aci.coefs[35,] <- c(id = "poa_pra4_t2_ch2", pathway = "C3",t(coef(poa_pra4_t2_ch2)))

poa_pra4_t2_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "poa_pra4_t2_ch2_27.5" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra4_t2_ch2_27.5)
chamber2_aci.coefs[36,] <- c(id = "poa_pra4_t2_ch2_27.5", pathway = "C3",t(coef(poa_pra4_t2_ch2_27.5)))


poa_pra5_t3_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "poa_pra5_t3_ch2" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra5_t3_ch2)
chamber2_aci.coefs[37,] <- c(id = "poa_pra5_t3_ch2",pathway = "C3", t(coef(poa_pra5_t3_ch2)))


poa_pra5_t3_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "poa_pra5_t3_ch2_27.5" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra5_t3_ch2_27.5)
chamber2_aci.coefs[38,] <- c(id = "poa_pra5_t3_ch2_27.5", pathway = "C3",t(coef(poa_pra5_t3_ch2_27.5)))

poa_pra6_t3_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "poa_pra6_t3_ch2" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra6_t3_ch2)
chamber2_aci.coefs[39,] <- c(id = "poa_pra6_t3_ch2", pathway = "C3",t(coef(poa_pra6_t3_ch2)))

poa_pra6_t3_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "poa_pra6_t3_ch2_27.5" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra6_t3_ch2_27.5)
chamber2_aci.coefs[40,] <- c(id = "poa_pra6_t3_ch2_27.5", pathway = "C3",t(coef(poa_pra6_t3_ch2_27.5)))


poa_pra7_t4_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "poa_pra7_t4_ch2" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra7_t4_ch2)
chamber2_aci.coefs[41,] <- c(id = "poa_pra7_t4_ch2", pathway = "C3",t(coef(poa_pra7_t4_ch2)))

poa_pra7_t4_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "poa_pra7_t4_ch2_27.5" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra7_t4_ch2_27.5)
chamber2_aci.coefs[42,] <- c(id = "poa_pra7_t4_ch2_27.5",pathway = "C3", t(coef(poa_pra7_t4_ch2_27.5)))

poa_pra8_t4_ch2 <- aci.prep %>% filter(keep.row == "yes" & id == "poa_pra8_t4_ch2" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra8_t4_ch2)
chamber2_aci.coefs[43,] <- c(id = "poa_pra8_t4_ch2", pathway = "C3",t(coef(poa_pra8_t4_ch2)))

poa_pra8_t4_ch2_27.5 <- aci.prep %>% filter(keep.row == "yes" & id == "poa_pra8_t4_ch2_27.5") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra8_t4_ch2_27.5)
chamber2_aci.coefs[44,] <- c(id = "poa_pra8_t4_ch2_27.5", pathway = "C3",t(coef(poa_pra8_t4_ch2_27.5)))


##write aci_coefs data

write.csv(chamber2_aci.coefs, "~/git/C3C4_GrowthChamber_Experiment/TxCO2_datasheets/aci.coefs_chamber2.csv", row.names = FALSE)


#####################################################################
# A/Ci curve temp standardization
#####################################################################
aci.coefs$Vcmax <- as.numeric(aci.coefs$Vcmax)
aci.coefs$Jmax <- as.numeric(aci.coefs$Jmax)
aci.coefs$Rd <- as.numeric(aci.coefs$Rd)
aci.coefs$TPU <- as.numeric(aci.coefs$TPU)

aci.coefs[, c(2:5)] <- round(aci.coefs[, c(2:5)], digits = 3)

aci.fits <- aci.coefs %>% left_join(aci.temps) %>%
  mutate(vcmax25 = temp_standardize(Vcmax, "Vcmax", standard.to = 25,
                                    tLeaf = Tleaf, tGrow = 22.5),
         jmax25 = temp_standardize(Jmax, "Jmax", standard.to = 25,
                                   tLeaf = Tleaf, tGrow = 22.5),
         jmax25.vcmax25 = jmax25 / vcmax25) %>%
  dplyr::select(id, tleaf = Tleaf, vcmax25, jmax25, jmax25.vcmax25, rd25 = Rd, tpu = TPU) %>%
  mutate_if(is.numeric, round, 3)

#####################################################################
# Extract snapshot measurements at 420 ppm CO2
#####################################################################
anet <- aci.prep %>%
  filter(CO2_r > 419.5 & CO2_r < 420.5) %>%
  group_by(id) %>%
  summarize(anet = mean(A),
            gsw = mean(gsw),
            ci.ca = mean(Ci) / mean(Ca)) %>%
  filter(id != "a_n_630_141" & id != "a_y_280_100" & id != "a_y_350_101") %>%
  mutate(id = ifelse(id == "a_n_630_141_b", "a_n_630_141", id),
         id = ifelse(id == "a_y_280_100_b", "a_y_280_100", id),
         id = ifelse(id == "a_y_350_101_b", "a_y_350_101", id))

#####################################################################
# Extract snapshot measurements at growth CO2 concentration
#####################################################################
agrowth_prep <- aci.prep %>%
  filter(id != "a_n_630_141" & id != "a_y_280_100" & id != "a_y_350_101") %>%
  mutate(id = ifelse(id == "a_n_630_141_b", "a_n_630_141", id),
         id = ifelse(id == "a_y_280_100_b", "a_y_280_100", id),
         id = ifelse(id == "a_y_350_101_b", "a_y_350_101", id)) %>%
  separate(col = id, into = c("co2", "inoc", "n.ppm", "rep"), remove = FALSE)

agrowth_amb <- agrowth_prep %>% filter(co2 == "a") %>%
  group_by(id) %>%
  filter(CO2_r > 419.5 & CO2_r < 420.5) %>%
  summarize(anet.growth = mean(A),
            gsw.growth = mean(gsw))

agrowth_elv <- agrowth_prep %>% filter(co2 == "e") %>%
  group_by(id) %>%
  filter(CO2_r > 990 & CO2_r < 1010) %>%
  summarize(anet.growth = mean(A),
            gsw.growth = mean(gsw))
agrowth <- agrowth_amb %>% full_join(agrowth_elv)

#####################################################################
# Merge snapshot measurements with A/Ci rate estimates
#####################################################################
photo.data <- anet %>% full_join(agrowth) %>% full_join(aci.fits) %>% 
  mutate_if(is.numeric, round, 3) %>%
  dplyr::select(id, tleaf, anet, anet.growth, gsw, gsw.growth, ci.ca, vcmax25, jmax25, jmax25.vcmax25, rd25, tpu)

#####################################################################
# Write file
#####################################################################
write.csv(photo.data, "../data_sheets/NxCO2_photo_data.csv", row.names = FALSE)



