###############################################################################
## Libraries
###############################################################################
library(dplyr)
library(tidyverse)
library(plantecophys)

###############################################################################
## Import files
###############################################################################
chamber3 <- read.csv("~/git/C3C4_GrowthChamber_Experiment/TxCO2_licorcleaned/TxCO2combinedataset/TXCO2_co2_resp_chamber_3.csv") %>%
  mutate(temp.setpoint = ifelse(Tleaf > 19 & Tleaf < 21, 
                                20,
                                ifelse(Tleaf > 24 & Tleaf < 29,
                                       27.5,
                                       ifelse(Tleaf > 30 & Tleaf < 36,
                                              35,
                                              NA))))



rd_chamber3 <- read.csv("~/git/C3C4_GrowthChamber_Experiment/TxCO2_licorcleaned/TxCO2combinedataset/TxCO2_rd_chamber_3.csv") %>%
  mutate(temp.setpoint = ifelse(Tleaf > 19 & Tleaf < 21, 
                                20,
                                ifelse(Tleaf > 24 & Tleaf < 29,
                                       27.5,
                                       ifelse(Tleaf > 30 & Tleaf < 36,
                                              35,
                                              NA))))


##grouping respiration by id and temp setpoint and then finding the mean
rd_chamber3_group_by <- group_by(rd_chamber3, id, temp.setpoint, chamber, machine)
rd_chamber3_mean <- summarise(rd_chamber3_group_by, rd_mean = mean(rd, na.rm = T))
head(rd_chamber3_mean)

###############################################################################
## Load custom fxns
###############################################################################
source("/Users/zinny/git/r_functions/stomatal_limitation.R")

###############################################################################
## Prep data frame to fit A/Ci curves
###############################################################################
aci.prep.chamber3 <- chamber3  %>%
  group_by(id) %>%
  dplyr::select(id, machine, A, Ci, Ca, gsw, 
                CO2_s,	CO2_r,	H2O_s,	H2O_r,
                Qin, VPDleaf, Flow,	Tair, Tleaf, temp.setpoint) %>%
  arrange(id) %>%
  left_join(rd_chamber3_mean, by = c("id", "temp.setpoint")) %>%
  filter(Qin > 1000) %>% # remove strange Qin values - seems to be Rd values merged into compiled sheet
  dplyr::select (-Tleaf) %>%
  group_by(id) %>%
  mutate(keep.row = "yes") %>%
  mutate(id = case_when(
    temp.setpoint <= 27.5 ~ paste0(id, "_27.5"), #adding _27.5 to the id's that were measured at 27.5
    TRUE ~ as.character(id)
  )) %>%
  data.frame()


 #aci.prep.chamber3<- rename(aci.prep.chamber3, Tleaf = temp.setpoint, rd = rd_mean)

#aci.prep$keep.row[c(,)] <- "no" # removes points that are outliers

##write.csv(aci.prep.chamber3, "../TxCO2_licorcleaned/TxCO2combinedataset/aci.prep.chamber3.csv", row.names = FALSE)
#####################################################################
# A/Ci curves
#####################################################################

#######################################
# Chamber 3 Ambient CO2 High Temperature
#######################################

ely_can9_t1_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "ely_can9_t1_ch3_27.5" &
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can9_t1_ch3_27.5)
#summary(ely_can9_t1_ch3_27.5)
chamber3_aci.coefs <- data.frame(id = "ely_can9_t1_ch3_27.5", pathway = "C3", t(coef(ely_can9_t1_ch3_27.5)))


ely_can9_t1_ch3 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "ely_can9_t1_ch3" &
                                                       Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can9_t1_ch3)
#summary(ely_can9_t1_ch3)
chamber3_aci.coefs[2,] <- c(id = "ely_can9_t1_ch3", pathway = "C3", t(coef(ely_can9_t1_ch3)))


ely_can10_t1_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "ely_can10_t1_ch3_27.5" &
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

#summary(ely_can10_t1_ch3_27.5)
plot(ely_can10_t1_ch3_27.5)
chamber3_aci.coefs[3,] <- c(id = "ely_can10_t1_ch3_27.5", pathway = "C3", t(coef(ely_can10_t1_ch3_27.5)))


ely_can10_t1_ch3 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "ely_can10_t1_ch3" &
                                              Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can10_t1_ch3)
#summary(ely_can10_t1_ch3)
chamber3_aci.coefs[4,] <- c(id = "ely_can10_t1_ch3", pathway = "C3", t(coef(ely_can10_t1_ch3)))



ely_can11_t2_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "ely_can11_t2_ch3_27.5" &
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can11_t2_ch3_27.5)
#summary(ely_can10_t1_ch3)
chamber3_aci.coefs[5,] <- c(id = "ely_can11_t2_ch3_27.5", pathway = "C3", t(coef(ely_can11_t2_ch3_27.5)))


ely_can11_t2_ch3 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "ely_can11_t2_ch3" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can11_t2_ch3)
#summary(ely_can11_t2_ch3)
chamber3_aci.coefs[6,] <- c(id = "ely_can11_t2_ch3", pathway = "C3", t(coef(ely_can11_t2_ch3)))



ely_can12_t2_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "ely_can12_t2_ch3_27.5" &
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


aci.prep.chamber3$keep.row[c(2498,2499,2500,2501,2502,2503,2504,2505,2506,2507,
                             2508,2509,2510,2511,2512,2513,2514,2515,2516,2517,2518,
                             2519,2520,2521,2522,2523,2524,2525,2526,2527,2528,2529,
                             2530,2531,2532,2533,2534,2535,2536,2537,2538,2539,2540,
                             2541,2542,2543,2544,2545,2546,2547,2548,2549,2550,2551,
                             2552,2553,2554,2555,2556,2557,2558,2559,2560,2561,
                             2562,2563,2564,2565,2566,2567,2568,2569,2570,2571,2572,
                             2573,2574,2575,2576,2577,2578,2579,2580,2581,2582,
                             2583,2584,2585,2586,2587,2588,2589,2590,2591,2592,2593)] <- "no" # removes points that are outliers
plot(ely_can12_t2_ch3_27.5)
#summary(ely_can12_t2_ch3_27.5)
chamber3_aci.coefs[7,] <- c(id = "ely_can12_t2_ch3_27.5", pathway = "C3", t(coef(ely_can12_t2_ch3_27.5)))


ely_can12_t2_ch3 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "ely_can12_t2_ch3" &
                                              Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

aci.prep.chamber3$keep.row[c(2617,2402,2403,2404,2405,2406,2407,2408,2409,2410,
                             2411,2412,2413,2414,2415,2416,2417,2418,2418,2419,
                             2420,2421,2422,2423,2424,2425,2426,2427,2428,
                             2429,2430,2431,2432,2433,2434,2435,2436,2437,2438,
                             2439,2440,2441,2442,2443,2444,2445,2446,2447,2448,
                             2449,2450,2451,2452,2453,2454,2455,2456,2456,2457,
                             2458,2459,2460,2461,2462,2463,2464,2465,2466,2467,
                             2468,2469,2470,2471,2472,2473,2474,2475,2476,2477,2478)] <- "no" # removes points that are outliers

summary(ely_can12_t2_ch3)
plot(ely_can12_t2_ch3)
chamber3_aci.coefs[8,] <- c(id = "ely_can12_t2_ch3",pathway = "C3", t(coef(ely_can12_t2_ch3)))



ely_can13_t3_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "ely_can13_t3_ch3_27.5" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can13_t3_ch3_27.5)
chamber3_aci.coefs[9,] <- c(id = "ely_can13_t3_ch3_27.5", pathway = "C3",t(coef(ely_can13_t3_ch3_27.5)))


ely_can13_t3_ch3 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "ely_can13_t3_ch3" & 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can13_t3_ch3)
chamber3_aci.coefs[10,] <- c(id = "ely_can13_t3_ch3", pathway = "C3",t(coef(ely_can13_t3_ch3)))


ely_can14_t3_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "ely_can14_t3_ch3_27.5" & 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can14_t3_ch3_27.5)
chamber3_aci.coefs[11,] <- c(id = "ely_can14_t3_ch3_27.5", pathway = "C3",t(coef(ely_can14_t3_ch3_27.5)))

ely_can14_t3_ch3 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "ely_can14_t3_ch3" & 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

#aci.prep.chamber3$keep.row[c(2997,2998,2999,3000,3001)] <- "no" # removes points that are outliers
summary(ely_can14_t3_ch3)
plot(ely_can14_t3_ch3)
chamber3_aci.coefs[12,] <- c(id = "ely_can14_t3_ch3", pathway = "C3",t(coef(ely_can14_t3_ch3)))


ely_can15_t4_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "ely_can15_t4_ch3_27.5" & 
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


aci.prep.chamber3$keep.row[c(2835,2834,2833,2866,2867,2868,
                             2869,2870,2871,2872,2873,2874,2875,2875,2876,2877,
                             2878,2879,2880,2881,2882,2883,2884)] <- "no" # removes points that are outliers
summary(ely_can15_t4_ch3_27.5)
plot(ely_can15_t4_ch3_27.5)
chamber3_aci.coefs[13,] <- c(id = "ely_can15_t4_ch3_27.5", pathway = "C3",t(coef(ely_can15_t4_ch3_27.5)))


ely_can15_t4_ch3 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "ely_can15_t4_ch3" & 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can15_t4_ch3)
chamber3_aci.coefs[14,] <- c(id = "ely_can15_t4_ch3", pathway = "C3",t(coef(ely_can15_t4_ch3)))


ely_can16_t4_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "ely_can16_t4_ch3_27.5" & 
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can16_t4_ch3_27.5)
chamber3_aci.coefs[15,] <- c(id = "ely_can16_t4_ch3_27.5", pathway = "C3",t(coef(ely_can16_t4_ch3_27.5)))



ely_can16_t4_ch3 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "ely_can16_t4_ch3" & 
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(ely_can16_t4_ch3)
chamber3_aci.coefs[16,] <- c(id = "ely_can16_t4_ch3", pathway = "C3",t(coef(ely_can16_t4_ch3)))



pas_smi9_t1_ch3_27.5 <- aci.prep.chamber3%>% filter(keep.row == "yes" & id == "pas_smi9_t1_ch3_27.5"& 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)



aci.prep.chamber3$keep.row[c(4910,4911,4944,4945,4946,
            4913,4914,4915,4916,4917,4918,4919,4920,4921,4922)] <- "no" # removes points that are outliers

plot(pas_smi9_t1_ch3_27.5)
chamber3_aci.coefs[17,] <- c(id = "pas_smi9_t1_ch3_27.5", pathway = "C3",t(coef(pas_smi9_t1_ch3_27.5)))



pas_smi9_t1_ch3 <- aci.prep.chamber3%>% filter(keep.row == "yes" & id == "pas_smi9_t1_ch3"& 
                                                      Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(pas_smi9_t1_ch3)
chamber3_aci.coefs[18,] <- c(id = "pas_smi9_t1_ch3", pathway = "C3",t(coef(pas_smi9_t1_ch3)))


pas_smi10_t1_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "pas_smi10_t1_ch3_27.5"& 
                                                 Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)



#aci.prep.chamber3$keep.row[c()] <- "no" # removes points that are outliers

plot(pas_smi10_t1_ch3_27.5)
chamber3_aci.coefs[19,] <- c(id = "pas_smi10_t1_ch3_27.5", pathway = "C3",t(coef(pas_smi10_t1_ch3_27.5)))



pas_smi10_t1_ch3 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "pas_smi10_t1_ch3"& 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi10_t1_ch3)
chamber3_aci.coefs[20,] <- c(id = "pas_smi10_t1_ch3", pathway = "C3",t(coef(pas_smi10_t1_ch3)))


pas_smi11_t2_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "pas_smi11_t2_ch3_27.5"& 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(pas_smi11_t2_ch3_27.5)
chamber3_aci.coefs[21,] <- c(id = "pas_smi11_t2_ch3_27.5", pathway = "C3",t(coef(pas_smi11_t2_ch3_27.5)))


pas_smi11_t2_ch3 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "pas_smi11_t2_ch3"& 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(pas_smi11_t2_ch3)
chamber3_aci.coefs[22,] <- c(id = "pas_smi11_t2_ch3", pathway = "C3",t(coef(pas_smi11_t2_ch3)))



pas_smi12_t2_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "pas_smi12_t2_ch3_27.5"& 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


aci.prep.chamber3$keep.row[c()] <- "no" # removes points that are outliers

plot(pas_smi12_t2_ch3_27.5)
chamber3_aci.coefs[23,] <- c(id = "pas_smi12_t2_ch3_27.5", pathway = "C3",t(coef(pas_smi12_t2_ch3_27.5)))




pas_smi12_t2_ch3<- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "pas_smi12_t2_ch3"& 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi12_t2_ch3)
chamber3_aci.coefs[24,] <- c(id = "pas_smi12_t2_ch3", pathway = "C3",t(coef(pas_smi12_t2_ch3)))



pas_smi13_t3_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "pas_smi13_t3_ch3_27.5"& 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(pas_smi13_t3_ch3_27.5)
chamber3_aci.coefs[25,] <- c(id = "pas_smi13_t3_ch3_27.5", pathway = "C3",t(coef(pas_smi13_t3_ch3_27.5)))


pas_smi13_t3_ch3<- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "pas_smi13_t3_ch3"& 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

plot(pas_smi13_t3_ch3)
chamber3_aci.coefs[26,] <- c(id = "pas_smi13_t3_ch3", pathway = "C3",t(coef(pas_smi13_t3_ch3)))


pas_smi14_t3_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "pas_smi14_t3_ch3_27.5"& 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(pas_smi14_t3_ch3_27.5)
chamber3_aci.coefs[27,] <- c(id = "pas_smi14_t3_ch3_27.5", pathway = "C3",t(coef(pas_smi14_t3_ch3_27.5)))



pas_smi14_t3_ch3 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "pas_smi14_t3_ch3" & 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(pas_smi14_t3_ch3)
chamber3_aci.coefs[28,] <- c(id = "pas_smi14_t3_ch3", pathway = "C3",t(coef(pas_smi14_t3_ch3)))


pas_smi15_t4_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "pas_smi15_t4_ch3_27.5"& 
                                                  Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(pas_smi15_t4_ch3_27.5)
chamber3_aci.coefs[29,] <- c(id = "pas_smi15_t4_ch3_27.5", pathway = "C3",t(coef(pas_smi15_t4_ch3_27.5)))




pas_smi15_t4_ch3 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "pas_smi15_t4_ch3"& 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(pas_smi15_t4_ch3)
chamber3_aci.coefs[30,] <- c(id = "pas_smi15_t4_ch3", pathway = "C3",t(coef(pas_smi15_t4_ch3)))



pas_smi16_t4_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "pas_smi16_t4_ch3_27.5"& 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(pas_smi16_t4_ch3_27.5)
chamber3_aci.coefs[31,] <- c(id = "pas_smi16_t4_ch3_27.5", pathway = "C3",t(coef(pas_smi16_t4_ch3_27.5)))


pas_smi16_t4_ch3 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "pas_smi16_t4_ch3"& 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(pas_smi16_t4_ch3)
chamber3_aci.coefs[32,] <- c(id = "pas_smi16_t4_ch3", pathway = "C3",t(coef(pas_smi16_t4_ch3)))


poa_pra9_t1_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "poa_pra9_t1_ch3_27.5" & 
                                         Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

aci.prep.chamber3$keep.row[c(6748,6749,6750,6751,6752,6753,6754,6755,6756,6757,
                          6758,6759,6760,6761,6762,6763,6764,6765,6766,6767,
                          6768,6769,6770,6771,6772,6773,6774,6775,6776,6777,
                              6778,6779,6780,6781,6782,6783,6784,6785,6786,6787,
                              6788,6789,6790,6791,6792,6793,6794,6795,6796,6797,
                                6798,6799,6800,6801,6802,6803,6804,6805,6806,6807,
                              6808,6809,6810,6811,6812,6813,6814,6815,6816,6817,
                             6818,6819,6820,6821,6822,6823,6824,6825,6826,6827,
                              6828,6829,6830,6831,6832,6833,6834,6835,6836,6837,
                                              6838,6839,6840,6841,6842,6843)] <- "no" # removes points that are outliers

plot(poa_pra9_t1_ch3_27.5)
chamber3_aci.coefs[33,] <- c(id = "poa_pra9_t1_ch3_27.5",pathway = "C3", t(coef(poa_pra9_t1_ch3_27.5)))



poa_pra9_t1_ch3 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "poa_pra9_t1_ch3"& 
                                                       Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra9_t1_ch3)
chamber3_aci.coefs[34,] <- c(id = "poa_pra9_t1_ch3",pathway = "C3", t(coef(poa_pra9_t1_ch3)))


poa_pra10_t1_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "poa_pra10_t1_ch3_27.5"& 
                                                  Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

aci.prep.chamber3$keep.row[c(5115,5116,5117,5118,5119,5120,5121,5122,5123,5124,
                             5125,5126,5127,5128,5129,5130,5131,5132,5133,5134,
                             5135,5136,5137,5138,5139,5140,5141,5142,5143,5144,
                             5145,5146,5147,5148,5149,5150,5151,5152,5153,5154,
                             5155,5156,5157,5158,5159,5160,5161,5162,5163,5164,
                             5165,5166,5167,5168,5169,5170,5171,5172,5173,5174,
                             5175,5176,5177,5178,5179,5180,5181,5182,5183,5184,
                             5185,5186,5187,5188,5189,5190,5191,5192,5193,5194,
                             5195,5196,5197,5198,5199,5200,5201,5202,5203,5204,
                             5205,5206,5207,5208,5209,5210)] <- "no" # removes points that are outliers


plot(poa_pra10_t1_ch3_27.5)
chamber3_aci.coefs[35,] <- c(id = "poa_pra10_t1_ch3_27.5",pathway = "C3", t(coef(poa_pra10_t1_ch3_27.5)))


poa_pra10_t1_ch3 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "poa_pra10_t1_ch3"& 
                                                  Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

aci.prep.chamber3$keep.row[c(5211,5212,5213,5214,5215,5216,5217,5218,5219,5220,
                             5221,5222,5223,5224,5225,5226,5227,5228,5229,5230,
                             5231,5232,5233,5234,5235,5236,5237,5238,5239,5240,
                             5241,5242,5243,5244,5245,5246,5247,5248,5249,5250,
                             5251,5252,5253,5254,5255,5256,5257,5258,5259,
                             5260,5261,5262,5263,5264,5265,5266,5267,5268,5269,
                             5270,5271,5272,5273,5274,5275,5276,5277,5278,5279,
                             5280,5281,5282,5283,5284,5285,5286,5287,5288,5289,
                             5290,5291,5292,5293,5294,5295,5296,5297,5298,5299,
                             5300,5301,5302,5303,5304,5305,5306)] <- "no" # removes points that are outliers

plot(poa_pra10_t1_ch3)
chamber3_aci.coefs[36,] <- c(id = "poa_pra10_t1_ch3",pathway = "C3", t(coef(poa_pra10_t1_ch3)))


poa_pra12_t2_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "poa_pra12_t2_ch3_27.5"& 
                                                  Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

aci.prep.chamber3$keep.row[c(5691,5692,5746,5748, 5521,5522,
                             5693,5694,5695,5696,5697,5698,5699,5700,5701,5702,
                             5703,5704,5705,5706,5707,5708,5709,5710,5711,5712,
                             5713,5714,5715,5716,5717,5718,5719,5720,5721,5722,
                             5723,5724,5725,5726,5727,5728,5729,5730,5731,5732,
                             5733,5734,5735,5736,5737,5738,5739,5740,5741,5742,
                             5743,5744,5745,5746,5747,5748,5749,5750,5751,5752,
                             5753,5754,5755,5756,5757,5758,5759,5760,5761,5762,5763,5764,
                             5765,5766,5767,5768,5769,5770,5771,5772,5773,5774,
                             5775,5776,5777,5778,5779,5780,5781,5782,5783,5784,
                             5785,5786)] <- "no" # removes points that are outliers


plot(poa_pra12_t2_ch3_27.5)
chamber3_aci.coefs[37,] <- c(id = "poa_pra12_t2_ch3_27.5",pathway = "C3", t(coef(poa_pra12_t2_ch3_27.5)))


poa_pra12_t2_ch3<- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "poa_pra12_t2_ch3"& 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

aci.prep.chamber3$keep.row[c(5787,5788,5789,5790,5791,5792,5793,5794,5795,5796,
                             5797,5798,5799,5800,5801,5802,5803,5804,5805,5806,
                             5807,5808,5809,5810,5811,5812,5813,5814,5815,5816,
                             5817,5818,5819,5820,5821,5822,5823,5824,5825,5826,
                             5827,5828,5829,5830,5831,5832,5833,5834,5835,5836,
                             5837,5838,5839,5840,5841,5842,5843,5844,5845,5846,
                             5847,5848,5849,5850,5851,5852,5853,5854,5855,5856,
                             5857,5858,5859,5860,5861,5862,5863,5864,5865,5866,
                             5867,5868,5869,5870,5871,5872,5873,5874,5875,5876,
                             5877,5878,5879,5880,5881,5882,5883)] <- "no" # removes points that are outliers

plot(poa_pra12_t2_ch3)
chamber3_aci.coefs[38,] <- c(id = "poa_pra12_t2_ch3",pathway = "C3", t(coef(poa_pra12_t2_ch3)))


poa_pra13_t3_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "poa_pra13_t3_ch3_27.5"& 
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra13_t3_ch3_27.5)
chamber3_aci.coefs[39,] <- c(id = "poa_pra13_t3_ch3_27.5",pathway = "C3", t(coef(poa_pra13_t3_ch3_27.5)))


poa_pra13_t3_ch3 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "poa_pra13_t3_ch3"& 
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)

aci.prep.chamber3$keep.row[c(6000,6001,6002,6005,6003,6004,6030,6031,6029,
                             6014,6007,6008,6010,6013,6011,6012, 6008,6009, 6006,
                             6007,5983,5984,6028,6027,6033)] <- "no" # removes points that are outliers

plot(poa_pra13_t3_ch3)
chamber3_aci.coefs[40,] <- c(id = "poa_pra13_t3_ch3",pathway = "C3", t(coef(poa_pra13_t3_ch3)))

  
poa_pra14_t3_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "poa_pra14_t3_ch3_27.5"& 
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra14_t3_ch3_27.5)
chamber3_aci.coefs[41,] <- c(id = "poa_pra14_t3_ch3_27.5",pathway = "C3", t(coef(poa_pra14_t3_ch3_27.5)))



poa_pra14_t3_ch3 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "poa_pra14_t3_ch3"& 
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra14_t3_ch3)
chamber3_aci.coefs[42,] <- c(id = "poa_pra14_t3_ch3",pathway = "C3", t(coef(poa_pra14_t3_ch3)))


poa_pra15_t4_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "poa_pra15_t4_ch3_27.5"& 
                                                   Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra15_t4_ch3_27.5)
chamber3_aci.coefs[43,] <- c(id = "poa_pra15_t4_ch3_27.5",pathway = "C3", t(coef(poa_pra15_t4_ch3_27.5)))

poa_pra15_t4_ch3 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "poa_pra15_t4_ch3"& 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra15_t4_ch3)
chamber3_aci.coefs[44,] <- c(id = "poa_pra15_t4_ch3",pathway = "C3", t(coef(poa_pra15_t4_ch3)))


poa_pra16_t4_ch3_27.5 <- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "poa_pra16_t4_ch3_27.5"& 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra16_t4_ch3_27.5)
chamber3_aci.coefs[45,] <- c(id = "poa_pra16_t4_ch3_27.5",pathway = "C3", t(coef(poa_pra16_t4_ch3_27.5)))


poa_pra16_t4_ch3<- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "poa_pra16_t4_ch3"& 
                                                        Ci < 2000) %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "Tleaf",
                         Ci = "Ci",
                         PPFD = "Qin", 
                         Rd = "rd"),
         fitTPU = TRUE, Tcorrect = FALSE, useRd = TRUE)


plot(poa_pra16_t4_ch3)
chamber3_aci.coefs[46,] <- c(id = "poa_pra16_t4_ch3",pathway = "C3", t(coef(poa_pra16_t4_ch3)))


##write aci_coefs data

write.csv(chamber3_aci.coefs, "~/git/C3C4_GrowthChamber_Experiment/TxCO2_licorcleaned/TxCO2combinedataset/aci.coefs_chamber3.csv", row.names = FALSE)





