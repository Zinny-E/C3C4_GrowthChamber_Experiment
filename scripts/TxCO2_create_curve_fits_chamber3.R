###############################################################################
## Libraries
###############################################################################
library(dplyr)
library(tidyverse)
library(plantecophys)

###############################################################################
## Import files
###############################################################################
chamber3 <- read.csv("~/git/C3C4_GrowthChamber_Experiment/licor_cleaned/TxCO2_combined_datasheets/TXCO2_co2_resp_chamber_3.csv") %>%
  mutate(temp.setpoint = ifelse(Tleaf > 19 & Tleaf < 21, 
                                20,
                                ifelse(Tleaf > 24 & Tleaf < 29,
                                       27.5,
                                       ifelse(Tleaf > 30 & Tleaf < 36,
                                              35,
                                              NA))))



rd_chamber3 <- read.csv("~/git/C3C4_GrowthChamber_Experiment/licor_cleaned/TxCO2_combined_datasheets/TxCO2_rd_chamber_3.csv") %>%
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

##write.csv(aci.prep.chamber3, "../licor_cleaned/TxCO2_combined_datasheets/aci.prep.chamber3.csv", row.names = FALSE)
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


#aci.prep.chamber3$keep.row[c(2211,2212,2213,2214,2215,2216,2217,2218,2219,2220,2221,2222,
#2223,2224,2225,2226,2227,2228,2229,2230,2231,2232,2233,2234,2235,2236,2237,2238,
#2239,2240,224,2242,2243,2244,2245,2246,2247,2248,2249,2250,2251,2252,2253,2254,
#2255,2256,2257,2258,2259,2260,2261,2262,2263,2264,2265,2266,2267,2268,2269,2270,
#2271,2272,2273,2274,2275,2276,2277,2278,2279,2280,2281,2282,2283,2284,2285,2286,2287,2288,
#2289,2290,2291,2292,2293,2294,2295,2296,2297,2298,2299,2300,2301,2302,2303,2304,2305,2306)] <- "no" # removes points that are outliers
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

aci.prep.chamber3$keep.row[c(2115,2116,2117,2118,2119,2120,2121,2122,2123,2124,2125,2126,
                    2127,2128,2129,2130,2131,2132,2133,2134,2135,2136,2137,2138,2139,2140,
                    2141,2142,2143,2144,2145,2146,2147,2148,2149,2150,2151,2152,2153,2154,
                    2155,2156,2157,2158,2159,2160,2161,2162,2163,2164,2165,2166,2167,2168,
                    2169,2170,2171,2172,2173,2174,2175,2176,2177,2178,2179,2180,2181,2182,2183,
                    2184,2185,2186,2187,2188,2189,2190,2191,2192,2193,2194,
                    2195,2196,2197,219,2199,2200,2201,2202,2203,2204,2205,2206,
                    2207,2208,2209,2210,2330)] <- "no" # removes points that are outliers

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

#aci.prep.chamber3$keep.row[c(2714,2713,2712,2711,2758,2757,2710)] <- "no" # removes points that are outliers
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



#aci.prep.chamber3$keep.row[c(4648,4647,4662,4661)] <- "no" # removes points that are outliers

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



aci.prep.chamber3$keep.row[c(4758,4757,4756)] <- "no" # removes points that are outliers

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



aci.prep.chamber3$keep.row[c(3417,3418,3419,3420,3421,3422,3423,3424,3425,3426,
                             3427,3428,3429,3430,3431,3432,3433,3434,3435,3436,3437,
                             3438,3439,3440,3441,3442,3443,3444,3445,3446,3447,3448,
                             3449,3450,3451,3452,3453,3454,3455,3456,3457,3458,3459,
                             3460,3461,3462,3463,3464,3465,3466,3467,3468,
                             3469,3470,3471,3472,3473,3474,3475,3476,3477,3478,3479,
                             3480,3481,3482,3483,3484,3485,3486,3415,3416,3417)] <- "no" # removes points that are outliers

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


aci.prep.chamber3$keep.row[c(3929,3930,3931,3932,3933,3934,3935,3936,3937,3938,
                             3939,3940,3941,3942,3943,3944)] <- "no" # removes points that are outliers

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


aci.prep.chamber3$keep.row[c(3929,3930,3931,3932,3933,3934,3935,3936,3937,3938,
                             3939,3940,3941,3942,3943,3944)] <- "no" # removes points that are outliers

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



pas_smi14_t3_ch3<- aci.prep.chamber3 %>% filter(keep.row == "yes" & id == "pas_smi14_t3_ch3"& 
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

#aci.prep.chamber3$keep.row[c(6272,6273,6274,6275,6276,6277,6278,6279,6280,6281,6282,6283,6284,6285,6286,6287,6288,
#6289,6290,6291,6292,6293,6294,6295,6296,6297,6298,6299,6300,6301,6302,6303,6304,
#6305,6306,6307,6308,6309,6310,6311,6312,6313,6314,6315,6316,6317,6318,6319,6320,
#6321,6322,6323,6324,6325,6326,6327,6328,6329,6330,6331,6332,6333,6334,6335,6336,
#6337,6338,6339,6340,6341,6342,6343,6344,6345,6346,6347,6348,6349,6350,6351,6352,6353,
#6354,6355,6356,6357,6358,6359,6360,6361,6362,6363,6364,6365,6366,6367,6368,6369,6370, 
#6371,6372,6373,6374,6375,6376,6377,6378,6379,6380,6381,6382,6383,6384,6385,6386,6387,
#6388,6389,6390,6391,6392,6393,6394,6395,6396,6397,6398,6399,6400,6401,6402,6403,6404,
#6405,6406,6407,6408,6409,6410,6411,6412,6413,6414,6415,6416,6417,6418,6419,6420,6421,
#6422,6423,6424,6425,6426,6427,6428,6429,6430,6431,6432,6433,6434,6435,6436,6437,6438,
#6439,6440,6441,6442,6443,6444,6445,6446,6447,6448,6449,6450,6451,6452,6453,6454,6455,
#6456,6457,6458,6459,6460,6461,6462,6463)] <- "no" # removes points that are outliers

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

aci.prep.chamber3$keep.row[c(5215,5216,5217,5218,5219,5220,5221,5222,5223,5224,5225,5226,5227,5228,5229,5230,5231,5232,5233,
                             5234,5235,5236,5237,5238,5239,5240,5241,5242,5243,5244,5245,5246,5247,5248,5249,5250,5251,5252,5253,	
                             5254,5255,5256,5257,5258,5259,5260,5261,5262,5263,5264,5265,5266,5267,5268,5269,5270,5271,5272,5273,
                             5274,5275,5276,5277,5278,5279,5280,5281,5282,5283,5284,5285,5286,5287,5288,5289,5290,5291,5292,5293,
                             5294,5295,5296,5297,5298,5299,5300,5301,5302,5303,5304,5305,5306,5307,5308,5309,5310)] <- "no" # removes points that are outliers


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

aci.prep.chamber3$keep.row[c(5311,5312,5313,5314,5315,5316,5317,5318,5319,5320,5321,5322,5323,5324,
                             5325,5326,5327,5328,5329,5330,5331,5332,5333,5334,5335,5336,5337,5338,5339,
                             5340,5341,5342,5343,5344,5345,5346,5347,5348,5349,5350,5351,5352,5353,5354,5355,
                             5356,5357,5358,5359,5360,5361,5362,5363,5364,5365,5366,5367,5368,5369,5370,5371,
                             5372,5373,5374,5375,5376,5377,5378,5379,5380,5381,5382,5383,5384,5385,5386,5387,
                             5388,5389,5390,5391,5392,5393,5394,5395,5396,5397,5398,5399,5400,5401,5402,5403,5404,5405,5406,5407)] <- "no" # removes points that are outliers

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

write.csv(chamber3_aci.coefs, "~/git/C3C4_GrowthChamber_Experiment/TxCO2_datasheets/aci.coefs_chamber3.csv", row.names = FALSE)


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



