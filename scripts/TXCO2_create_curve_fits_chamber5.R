###############################################################################
## Libraries
###############################################################################
library(dplyr)
library(tidyverse)
library(plantecophys)

###############################################################################
## Import files
###############################################################################
chamber5<- read.csv("~/git/C3C4_GrowthChamber_Experiment/licor_cleaned/TxCO2_combined_datasheets/TXCO2_co2_resp_chamber_5.csv") %>%
  mutate(temp.setpoint = ifelse(Tleaf > 19 & Tleaf < 21, 
                                20,
                                ifelse(Tleaf > 27 & Tleaf < 28,
                                       27.5,
                                       ifelse(Tleaf > 34 & Tleaf < 36,
                                              35,
                                              NA))))



rd_chamber5 <- read.csv("~/git/C3C4_GrowthChamber_Experiment/licor_cleaned/TxCO2_combined_datasheets/TxCO2_rd_chamber_5.csv") %>%
  mutate(temp.setpoint = ifelse(Tleaf > 19 & Tleaf < 21, 
                                20,
                                ifelse(Tleaf > 27 & Tleaf < 28,
                                       27.5,
                                       ifelse(Tleaf > 34 & Tleaf < 36,
                                              35,
                                              NA))))

rd_chamber5_group_by <- group_by(rd_chamber5, id, temp.setpoint, chamber, machine)
rd_chamber5_mean <- summarise(rd_chamber5_group_by, rd_mean = mean(rd, na.rm = T))
head(rd_chamber5_mean)



###############################################################################
## Load custom fxns
###############################################################################
source("/Users/zinny/git/r_functions/stomatal_limitation.R")

###############################################################################
## Prep data frame to fit A/Ci curves
###############################################################################
aci_prep_chamber5 <- chamber5  %>%
  group_by(id) %>%
  dplyr::select(id, machine, A, Ci, Ca, gsw, 
                CO2_s,	CO2_r,	H2O_s,	H2O_r,
                Qin, VPDleaf, Flow,	Tair, Tleaf, temp.setpoint) %>%
  arrange(id) %>%
  left_join(rd_chamber5_mean, by = c("id", "temp.setpoint")) %>%
  filter(Qin > 1000) %>% # remove strange Qin values - seems to be Rd values merged into compiled sheet
  dplyr::select (-Tleaf) %>%
  group_by(id) %>%
  mutate(keep.row = "yes") %>%
  mutate(id = case_when(
    temp.setpoint >= 27.5 ~ paste0(id, "_27.5"), #addinf _27.5 to the id's that were measured at 27.5
    TRUE ~ as.character(id)
  )) %>%
  data.frame()

aci_prep_chamber5 <- rename(aci_prep_chamber5, Tleaf = temp.setpoint, rd = rd_mean)



#aci.prep$keep.row[c(,)] <- "no" # removes points that are outliers
  
#####################################################################
# A/Ci curves
#####################################################################
  
#######################################
# Chamber 5 Elevated CO2 Low Tempertaure
#######################################
ely_can25_t1_ch5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "ely_can25_t1_ch5" &
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can25_t1_ch5)
chamber5_aci.coefs <- data.frame(id = "ely_can25_t1_ch5", pathway = "C3", t(coef(ely_can25_t1_ch5)))


ely_can25_t1_ch5_27.5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "ely_can25_t1_ch5_27.5") %>%
  
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can25_t1_ch5_27.5)
chamber5_aci.coefs[2,] <- c(id = "ely_can25_t1_ch5_27.5", pathway ="C3", t(coef(ely_can25_t1_ch5_27.5)))


ely_can26_t1_ch5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "ely_can26_t1_ch5") %>%
                                                   
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can26_t1_ch5)
chamber5_aci.coefs[3,] <- c(id = "ely_can26_t1_ch5", pathway ="C3", t(coef(ely_can26_t1_ch5)))

ely_can26_t1_ch5_27.5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "ely_can26_t1_ch5_27.5") %>%
                                                      
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can26_t1_ch5_27.5)
chamber5_aci.coefs[4,] <- data.frame(id = "ely_can26_t1_ch5_27.5", pathway ="C3", t(coef(ely_can26_t1_ch5_27.5)))


ely_can27_t2_ch5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "ely_can27_t2_ch5") %>%
  
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can27_t2_ch5)
chamber5_aci.coefs[5,] <- data.frame(id = "ely_can27_t2_ch5", pathway ="C3", t(coef(ely_can27_t2_ch5)))



ely_can27_t2_ch5_27.5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "ely_can27_t2_ch5_27.5") %>%
  
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can27_t2_ch5_27.5)
chamber5_aci.coefs[6,] <- c(id = "ely_can27_t2_ch5_27.5", pathway ="C3", t(coef(ely_can27_t2_ch5_27.5)))


ely_can28_t2_ch5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "ely_can28_t2_ch5" &
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(ely_can28_t2_ch5)
chamber5_aci.coefs[7,] <- data.frame(id = "ely_can28_t2_ch5", pathway ="C3", t(coef(ely_can28_t2_ch5)))


ely_can28_t2_ch5_27.5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "ely_can28_t2_ch5_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can28_t2_ch5_27.5)
chamber5_aci.coefs[8,] <- data.frame(id = "ely_can28_t2_ch5_27.5", pathway ="C3", t(coef(ely_can28_t2_ch5_27.5)))


ely_can29_t3_ch5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "ely_can29_t3_ch5" &
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can29_t3_ch5)
chamber5_aci.coefs[9,] <- data.frame(id = "ely_can29_t3_ch5", pathway ="C3", t(coef(ely_can29_t3_ch5)))


ely_can29_t3_ch5_27.5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "ely_can29_t3_ch5_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can29_t3_ch5_27.5)
chamber5_aci.coefs[10,] <- data.frame(id = "ely_can29_t3_ch5_27.5", pathway ="C3", t(coef(ely_can29_t3_ch5_27.5)))


ely_can30_t3_ch5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "ely_can30_t3_ch5" &
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can30_t3_ch5)
aci.coefs[11,] <- data.frame(id = "ely_can30_t3_ch5", pathway ="C3", t(coef(ely_can30_t3_ch5)))


ely_can30_t3_ch5_27.5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "ely_can30_t3_ch5_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can30_t3_ch5_27.5)
chamber5_aci.coefs[12,] <- data.frame(id = "ely_can30_t3_ch5_27.5", pathway ="C3", t(coef(ely_can30_t3_ch5_27.5)))


ely_can31_t4_ch5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "ely_can31_t4_ch5" &
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can31_t4_ch5)
chamber5_aci.coefs[13,] <- data.frame(id = "ely_can31_t4_ch5", pathway ="C3", t(coef(ely_can31_t4_ch5)))


ely_can31_t4_ch5_27.5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "ely_can31_t4_ch5_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can31_t4_ch5_27.5)
chamber5_aci.coefs[14,] <- data.frame(id = "ely_can31_t4_ch5_27.5", pathway ="C3", t(coef(ely_can31_t4_ch5_27.5)))


ely_can32_t4_ch5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "ely_can32_t4_ch5" &
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can32_t4_ch5)
chamber5_aci.coefs[15,] <- data.frame(id = "ely_can32_t4_ch5", pathway ="C3", t(coef(ely_can32_t4_ch5)))


ely_can32_t4_ch5_27.5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "ely_can32_t4_ch5_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can32_t4_ch5_27.5)
chamber5_aci.coefs[16,] <- data.frame(id = "ely_can32_t4_ch5_27.5", pathway ="C3", t(coef(ely_can32_t4_ch5_27.5)))




pas_smi25_t1_ch5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "pas_smi25_t1_ch5" &
                                          Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi25_t1_ch5)
chamber5_aci.coefs[17,] <- data.frame(id = "pas_smi25_t1_ch5", pathway ="C3", t(coef(pas_smi25_t1_ch5)))



pas_smi25_t1_ch5_27.5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "pas_smi25_t1_ch5_27.5" &
                                           Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi25_t1_ch5_27.5)
chamber5_aci.coefs[18,] <- data.frame(id = "pas_smi25_t1_ch5_27.5", pathway ="C3", t(coef(pas_smi25_t1_ch5_27.5)))


pas_smi26_t1_ch5<- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "pas_smi26_t1_ch5" &
                                                Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
#aci_prep_chamber5$keep.row[c(2786,2787,2788,2789,2790,
#                             2791,2792,2733,2794,2795,2796,2797,
 #                           2798,2799,2880,2801,2802,2803,2804,2805,2806,
  #                          2807,2808,2809,2810,2811,2812,2813,2814,
   #                         2815,2817,2816,2818,2819,2874,2875, 2876, 
    #                        2862,2863,2864, 2865,2866, 2867,2868, 2869,2870,2871,2872,2873,
     #                       2877, 2878, 2879,2880,2847,2848,2849,2850,2851,2852,
      #                      2853,2854,2855, 2856,2857, 2858, 2859,2860, 2861)] <- "no" # removes points that are outliers
                             
plot(pas_smi26_t1_ch5)
chamber5_aci.coefs[19,] <- data.frame(id = "pas_smi26_t1_ch5", pathway ="C3", t(coef(pas_smi26_t1_ch5)))


pas_smi26_t1_ch5_27.5<- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "pas_smi26_t1_ch5_27.5") %>%
                                                
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi26_t1_ch5_27.5)
chamber5_aci.coefs[20,] <- data.frame(id = "pas_smi26_t1_ch5_27.5", pathway ="C3", t(coef(pas_smi26_t1_ch5_27.5)))

pas_smi27_t2_ch5<- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "pas_smi27_t2_ch5" &
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi27_t2_ch5)
chamber5_aci.coefs[21,] <- data.frame(id = "pas_smi27_t2_ch5", pathway ="C3", t(coef(pas_smi27_t2_ch5)))

pas_smi27_t2_ch5_27.5<- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "pas_smi27_t2_ch5_27.5" &
                                                  Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi27_t2_ch5_27.5)
chamber5_aci.coefs[22,] <- data.frame(id = "pas_smi27_t2_ch5_27.5", pathway ="C3", t(coef(pas_smi27_t2_ch5_27.5)))


pas_smi28_t2_ch5<- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "pas_smi28_t2_ch5" &
                                                  Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi28_t2_ch5)
chamber5_aci.coefs[23,] <- data.frame(id = "pas_smi28_t2_ch5", pathway ="C3", t(coef(pas_smi28_t2_ch5)))

pas_smi28_t2_ch5_27.5<- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "pas_smi28_t2_ch5_27.5" &
                                                  Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi28_t2_ch5_27.5)
chamber5_aci.coefs[24,] <- data.frame(id = "pas_smi28_t2_ch5_27.5", pathway ="C3", t(coef(pas_smi28_t2_ch5_27.5)))


pas_smi29_t3_ch5<- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "pas_smi29_t3_ch5" &
                                                  Ci < 850) %>% # put at 850 cos thats the curve flattens 
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi29_t3_ch5)
chamber5_aci.coefs[25,] <- data.frame(id = "pas_smi29_t3_ch5", pathway ="C3", t(coef(pas_smi29_t3_ch5)))


pas_smi29_t3_ch5_27.5<- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "pas_smi29_t3_ch5_27.5" &
                                                  Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
#aci_prep_chamber5$keep.row[c(3531,3537,3541,3542,3545,
    #                         3546,3549,3539,3538)] <- "no"
plot(pas_smi29_t3_ch5_27.5)
chamber5_aci.coefs[26,] <- data.frame(id = "pas_smi29_t3_ch5_27.5", pathway ="C3", t(coef(pas_smi29_t3_ch5_27.5)))



pas_smi30_t3_ch5<- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "pas_smi30_t3_ch5" &
                                                       Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

#aci_prep_chamber5$keep.row[c(3576,3577, 3639,3640,3643,3644,
 #                            3623,3626,3624)] <- "no"


plot(pas_smi30_t3_ch5)
chamber5_aci.coefs[27,] <- data.frame(id = "pas_smi30_t3_ch5", pathway ="C3", t(coef(pas_smi30_t3_ch5)))


pas_smi30_t3_ch5_27.5<- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "pas_smi30_t3_ch5_27.5" &
                                                  Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi30_t3_ch5_27.5)
chamber5_aci.coefs[28,] <- data.frame(id = "pas_smi30_t3_ch5_27.5", pathway ="C3", t(coef(pas_smi30_t3_ch5_27.5)))


pas_smi31_t4_ch5<- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "pas_smi31_t4_ch5" &
                                                  Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
aci_prep_chamber5$keep.row[c(3768,3769,3767,3766,3764,
                             3960,3961,3957,3958)] <- "no"


plot(pas_smi31_t4_ch5)
chamber5_aci.coefs[29,]<- data.frame(id = "pas_smi31_t4_ch5", pathway ="C3", t(coef(pas_smi31_t4_ch5)))


pas_smi31_t4_ch5_27.5<- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "pas_smi31_t4_ch5_27.5" &
                                                  Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi31_t4_ch5_27.5)
chamber5_aci.coefs[30,] <- data.frame(id = "pas_smi31_t4_ch5_27.5", pathway ="C3", t(coef(pas_smi31_t4_ch5_27.5)))


pas_smi32_t4_ch5<- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "pas_smi32_t4_ch5" &
                                                  Ci < 1200) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi32_t4_ch5)
chamber5_aci.coefs[31,] <- data.frame(id = "pas_smi32_t4_ch5", pathway ="C3", t(coef(pas_smi32_t4_ch5)))


pas_smi32_t4_ch5_27.5<- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "pas_Smi32_t4_Ch5_27.5" &
                                                  Ci < 1200) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "T_eaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(pas_smi32_t4_ch5_27.5)
chamber5_aci.coefs[32,] <- data.frame(id = "pas_smi32_t4_ch5_27.5", pathway ="C3", t(coef(pas_smi32_t4_ch5_27.5)))




poa_pra25_t1_ch5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "poa_pra25_t1_ch5" &
                                           Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra25_t1_ch5)
chamber5_aci.coefs[33,] <- data.frame(id = "poa_pra25_t1_ch5", pathway ="C3", t(coef(poa_pra25_t1_ch5)))



poa_pra25_t1_ch5_27.5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "poa_pra25_t1_ch5_27.5" &
                                           Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

#aci_prep_chamber5$keep.row[c(4300,4321,4316,4320,4317,4315,4318,4311,4307,4308,
   #                          4321,4319,4316,4314,4313,4312,4311)] <- "no"

plot(poa_pra25_t1_ch5_27.5)
chamber5_aci.coefs[34,] <- data.frame(id = "poa_pra25_t1_ch5_27.5", pathway ="C3", t(coef(poa_pra25_t1_ch5_27.5)))


poa_pra26_t1_ch5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "poa_pra26_t1_ch5" &
                                           Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra26_t1_ch5)
chamber5_aci.coefs[35,] <- data.frame(id = "poa_pra26_t1_ch5", pathway ="C3", t(coef(poa_pra26_t1_ch5)))



poa_pra26_t1_ch5_27.5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "poa_pra26_t1_ch5_27.5" &
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


#aci_prep_chamber5$keep.row[c(4443,4444,4445,4446,4447,4448,4449,4450)] <- "no"

plot(poa_pra26_t1_ch5_27.5)
chamber5_aci.coefs[36,] <- data.frame(id = "poa_pra26_t1_ch5_27.5", pathway ="C3", t(coef(poa_pra26_t1_ch5_27.5)))


poa_pra27_t2_ch5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "poa_pra27_t2_ch5" &
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra27_t2_ch5)
chamber5_aci.coefs[37,] <- data.frame(id = "poa_pra27_t2_ch5", pathway ="C3", t(coef(poa_pra27_t2_ch5)))



poa_pra27_t2_ch5_27.5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "poa_pra27_t2_ch5_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

#aci_prep_chamber5$keep.row[c(4685,4867,4687,4688,4689,4690,4691,4692,4693,4693,4694,4695,4696,
#                          4697,4698,4697,4698,4699,4700,4701,4702,4703,4704,4704,4705,
  #                         4705,4706)] <- "no"

plot(poa_pra27_t2_ch5_27.5)
chamber5_aci.coefs[38,] <- data.frame(id = "poa_pra27_t2_ch5_27.5", pathway ="C3", t(coef(poa_pra27_t2_ch5_27.5)))


poa_pra28_t2_ch5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "poa_pra28_t2_ch5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra28_t2_ch5)
chamber5_aci.coefs[39,] <- data.frame(id = "poa_pra28_t2_ch5", pathway ="C3", t(coef(poa_pra28_t2_ch5)))


poa_pra28_t2_ch5_27.5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "poa_pra28_t2_ch5_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

 plot(poa_pra28_t2_ch5_27.5)
chamber5_aci.coefs[40,] <- data.frame(id = "poa_pra28_t2_ch5_27.5", pathway ="C3", t(coef(poa_pra28_t2_ch5_27.5)))


poa_pra29_t3_ch5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "poa_pra29_t3_ch5" &
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra29_t3_ch5)
chamber5_aci.coefs[41,] <- data.frame(id = "poa_pra29_t3_ch5", pathway ="C3", t(coef(poa_pra29_t3_ch5)))



poa_pra29_t3_ch5_27.5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "poa_pra29_t3_ch5_27.5" &
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

#aci_prep_chamber5$keep.row[c(5012,5013,5014,5014,5015,5016,5017,5018)] <- "no"

plot(poa_pra29_t3_ch5_27.5)
chamber5_aci.coefs[42,] <- data.frame(id = "poa_pra29_t3_ch5_27.5", pathway ="C3", t(coef(poa_pra29_t3_ch5_27.5)))

poa_pra30_t3_ch5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "poa_pra30_t3_ch5" &
                                           Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)
plot(poa_pra30_t3_ch5)
chamber5_aci.coefs[43,] <- data.frame(id = "poa_pra30_t3_ch5", pathway ="C3", t(coef(poa_pra30_t3_ch5)))

 
poa_pra30_t3_ch5_27.5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "poa_pra30_t3_ch5_27.5" &
                                           Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra30_t3_ch5_27.5)
chamber5_aci.coefs[44,] <- data.frame(id = "poa_pra30_t3_ch5_27.5", pathway ="C3", t(coef(poa_pra30_t3_ch5_27.5)))


poa_pra31_t4_ch5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "poa_pra31_t4_ch5" &
                                           Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

#aci_prep_chamber5$keep.row[c(5306)] <- "no"
plot(poa_pra31_t4_ch5)
chamber5_aci.coefs[45,] <- data.frame(id = "poa_pra31_t4_ch5", pathway ="C3", t(coef(poa_pra31_t4_ch5)))


poa_pra31_t4_ch5_27.5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "poa_pra31_t4_ch5_27.5" &
                                           Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd ="rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra31_t4_ch5_27.5)
chamber5_aci.coefs[46,] <- data.frame(id = "poa_pra31_t4_ch5_27.5", pathway ="C3", t(coef(poa_pra31_t4_ch5_27.5)))


poa_pra32_t4_ch5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "poa_pra32_t4_ch5" &
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra32_t4_ch5)
chamber5_aci.coefs[43,] <- data.frame(id = "poa_pra32_t4_ch5", pathway ="C3", t(coef(poa_pra32_t4_ch5)))


poa_pra32_t4_ch5_27.5 <- aci_prep_chamber5 %>% filter(keep.row == "yes" & id == "poa_pra32_t4_ch5_27.5" &
                                                      Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(poa_pra32_t4_ch5_27.5)
chamber5_aci.coefs[44,] <- data.frame(id = "poa_pra32_t4_ch5_27.5", pathway ="C3", t(coef(poa_pra32_t4_ch5_27.5)))

write.csv(chamber5_aci.coefs, "~/git/C3C4_GrowthChamber_Experiment/TxCO2_datasheets/aci.coefs_chamber5.csv", row.names = FALSE)


