###############################################################################
## Libraries
###############################################################################
library(dplyr)
library(tidyverse)



###change id to lowercase
names(df)[names(df) == "ID"] <- "id"

####################################################
#extract anet_420, gsw and ci.ca
####################################################
aci.prep_ch2 <- read.csv("../TxCO2_licorcleaned/TxCO2combinedataset/aci.prep_chamber2.csv")
aci.prep_ch3 <- read.csv("../TxCO2_licorcleaned/TxCO2combinedataset/aci.prep.chamber3.csv")

anet_ch2 <- aci.prep_ch2 %>%
  filter(CO2_r > 419.5 & CO2_r < 420.5) %>%
  group_by(id) %>%
  filter(!grepl("_27\\.5$", id)) %>%
  summarize(anet_growth = mean(A),
            gsw = mean(gsw),
            ci.ca = mean(Ci) / mean(Ca))



anet_ch3 <- aci.prep_ch3 %>%
  filter(CO2_r > 419.5 & CO2_r < 420.5) %>%
  group_by(id) %>%
  filter(!grepl("_27\\.5$", id)) %>%
  summarize(anet_growth = mean(A),
            gsw = mean(gsw),
            ci.ca = mean(Ci) / mean(Ca))


####################################################
#extract anet_1000, gsw and ci.ca
####################################################
aci.prep_ch4 <- read.csv("../TxCO2_licorcleaned/TxCO2combinedataset/aci.prep.chamber4.csv")
aci.prep_ch5 <- read.csv("../TxCO2_licorcleaned/TxCO2combinedataset/aci.prep.chamber5.csv")


anet_ch4 <- aci.prep_ch4 %>%
  filter(CO2_r > 990 & CO2_r < 1010) %>%
  group_by(id) %>%
  filter(!grepl("_27\\.5$", id)) %>%
  summarize(anet_growth = mean(A),
            gsw = mean(gsw),
            ci.ca = mean(Ci) / mean(Ca))



anet_ch5 <- aci.prep_ch5 %>%
  filter(CO2_r > 990 & CO2_r < 1010) %>%
  group_by(id) %>%
  filter(!grepl("_27\\.5$", id)) %>%
  summarize(anet_growth = mean(A),
            gsw = mean(gsw),
            ci.ca = mean(Ci) / mean(Ca))



####Join dataset
df <- df %>%
  left_join(anet_ch2, by = "id")

df <- df %>%
  left_join(anet_ch3, by = "id") %>%
  mutate(
    anet_growth = coalesce(anet_growth.x, anet_growth.y),
    gsw = coalesce(gsw.x, gsw.y),
    ci.ca = coalesce(ci.ca.x, ci.ca.y)
  ) %>%
  dplyr::select(-matches("\\.x$"), -matches("\\.y$")) %>%
  left_join(anet_ch4, by = "id") %>%
  mutate(
    anet_growth = coalesce(anet_growth.x, anet_growth.y),
    gsw = coalesce(gsw.x, gsw.y),
    ci.ca = coalesce(ci.ca.x, ci.ca.y)
  ) %>%
  dplyr::select(-matches("\\.x$"), -matches("\\.y$")) %>%
  left_join(anet_ch5, by = "id") %>%
  mutate(
    anet_growth = coalesce(anet_growth.x, anet_growth.y),
    gsw = coalesce(gsw.x, gsw.y),
    ci.ca = coalesce(ci.ca.x, ci.ca.y)
  ) %>%
  dplyr::select(-matches("\\.x$"), -matches("\\.y$"))









