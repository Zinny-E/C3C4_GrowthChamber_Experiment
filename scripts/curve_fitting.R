library(dplyr)
library(plantecophys)
library(tidyr)
library(ggplot2)
library(tibble)
library(minpack.lm)
library(nlme)
library(car)
library(emmeans)
library(ggpubr)

#########################################################
# reading fit_curve and pred_nls  function

source("~/Desktop/licor data/aci_curve_fit.R")
source("~/Desktop/licor data/pred_nls.R")

###########################################################
# equations and constants for temperature response curves

medlyn_Vcmax <- formula(Vcmax ~ kopt * (Hd * exp((Ha * (Tleaf - Topt))/(Tleaf * R * Topt)))/(Hd - (Ha * (1 - exp((Hd *(Tleaf - Topt))/(Tleaf * R * Topt))))))
medlyn_Jmax <- formula(Jmax ~ kopt * (Hd * exp((Ha * (Tleaf - Topt))/(Tleaf * R * Topt)))/(Hd - (Ha * (1 - exp((Hd *(Tleaf - Topt))/(Tleaf * R * Topt))))))

R <- 8.314
Hd = 200000

##################################################################

# read  and clean merged licor data
gas_exchange <- read.csv("~/Desktop/licor data/gas_echange.csv")

gas_exchange$plants  = gsub("^(\\d+)_.*", "\\1", gas_exchange$id)

resp_avg <- read.csv("~/Desktop/licor data/respiration_average.csv") %>%
  dplyr::rename(obs = X, plants = id)
resp_avg$plants <- as.character(resp_avg$plants)

# merge aci an respiration data
final_aci_data <- left_join(gas_exchange, resp_avg, by = "plants")


# create a column for temp setpoint
final_aci_data$temp_setpoint <- round(final_aci_data$Tleaf, digits = 0)
unique(final_aci_data$temp_setpoint)
# it has 17 and 45, which was not actual set point
# adjusting temp setpoint for values differnt from the actual setpoint

final_aci_data$temp_setpoint <- 
  ifelse(final_aci_data$temp_setpoint == 17, 16, final_aci_data$temp_setpoint)
final_aci_data$temp_setpoint <- 
  ifelse(final_aci_data$temp_setpoint == 45, 46, final_aci_data$temp_setpoint)
unique(final_aci_data$temp_setpoint)

final_aci_data <- final_aci_data %>%
  select(id, machine, A, Ci, Ca, gsw, 
         CO2_s,	CO2_r,	H2O_s,	H2O_r,
         Qin, VPDleaf, Flow,	rd, Tair, temp_setpoint) %>%
  dplyr::rename(Tleaf = temp_setpoint) %>%
  arrange(id)

########################################################
# fit curves one by one

# Plant 1
# Aci curves
# 16 degrees
p1_16 <- fitaci(filter(final_aci_data, id == "1_16"),
       varnames = list(ALEAF = "A",
                       Tleaf = "Tleaf",
                       Ci = "Ci",
                       PPFD = "Qin"),
       fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)

plot(p1_16)

coef_1_16 <- data.frame(id = "1_16", t(coef(p1_16)))

# 22 degrees
p1_22 <- fitaci(filter(final_aci_data, id == "1_22"),
                varnames = list(ALEAF = "A",
                                Tleaf = "Tleaf",
                                Ci = "Ci",
                                PPFD = "Qin"),
                fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)

plot(p1_22)

coef_1_22 <- data.frame(id = "1_22", t(coef(p1_22)))

# 28 degrees
p1_28 <- fitaci(filter(final_aci_data, id == "1_28"),
                varnames = list(ALEAF = "A",
                                Tleaf = "Tleaf",
                                Ci = "Ci",
                                PPFD = "Qin"),
                fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)

plot(p1_28)

coef_1_28 <- data.frame(id = "1_28", t(coef(p1_28)))

# 34 degrees
p1_34 <- fitaci(filter(final_aci_data, id == "1_34"),
                varnames = list(ALEAF = "A",
                                Tleaf = "Tleaf",
                                Ci = "Ci",
                                PPFD = "Qin"),
                fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)

plot(p1_34)

coef_1_34 <- data.frame(id = "1_34", t(coef(p1_34)))

# 40 degrees
p1_40 <- fitaci(filter(final_aci_data, id == "1_40"),
                varnames = list(ALEAF = "A",
                                Tleaf = "Tleaf",
                                Ci = "Ci",
                                PPFD = "Qin"),
                fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)

plot(p1_40)

coef_1_40 <- data.frame(id = "1_40", t(coef(p1_40)))

# 46 degress
# 1_46 had data for two curves. So fitting one

data_1_46 <- filter(final_aci_data, id == "1_46")[1:96, ]
# first 96 rows gives better curves
data_1_46 <- filter(data_1_46, !(Ci >= 350 & Ci <= 600)) 
# removing points that represents the depression in the curve

p1_46 <- fitaci(data_1_46,varnames = list(ALEAF = "A",
                                    Tleaf = "Tleaf",
                                    Ci = "Ci",
                                    PPFD = "Qin"),
                  fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(p1_46)

coef_1_46 <- data.frame(id = "1_46", t(coef(p1_46)))

coef_1 <- rbind(coef_1_16, coef_1_22, coef_1_28, coef_1_34, coef_1_40, coef_1_46)
coef_1$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_1$id))
coef_1$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_1$id)
coef_1 <- coef_1 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_1 <- nlsLM(medlyn_Vcmax, data = coef_1,
                   start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_1)
# note: T opt is within the measured temperature range
# extracting parameters
t_vcmax_para_1 <- pred_nls(t_model_1, coef_1)

# temperature response curve for Jmax

t_model_Jmax_1 <- nlsLM(medlyn_Jmax, data = coef_1,
                   start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_1)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_1 <- pred_nls(t_model_Jmax_1, coef_1)



## Plant 2
p2_16 <- fit_curve(id_name =  "2_16")
plot(p2_16[[1]])
p2_22 <- fit_curve(id_name =  "2_22")
plot(p2_22[[1]])
p2_28 <- fit_curve(id_name =  "2_28")
plot(p2_28[[1]])
p2_34 <- fit_curve(id_name =  "2_34")
plot(p2_34[[1]])
p2_40 <- fit_curve(id_name =  "2_40")
plot(p2_40[[1]])
# fit_curve(id_name =  "2_46") have no data

coef_2 <- rbind(p2_16[[2]], p2_22[[2]], p2_28[[2]], p2_34[[2]], p2_40[[2]])
# note: has higher vcmax values

coef_2$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_2$id))
coef_2$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_2$id)
coef_2 <- coef_2 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_2 <- nlsLM(medlyn_Vcmax, data = coef_2,
                   start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_2)
# note: T opt is within the measured temperature range
t_vcmax_para_2 <- pred_nls(t_model_2, coef_2)

# temperature response curve for Jmax

t_model_Jmax_2 <- nlsLM(medlyn_Jmax, data = coef_2,
                        start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_2)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_2 <- pred_nls(t_model_Jmax_2, coef_2)

# Plant 3
p3_16 <- fit_curve(id_name =  "3_16")
plot(p3_16[[1]])
p3_22 <- fit_curve(id_name =  "3_22")
plot(p3_22[[1]])
p3_28 <- fit_curve(id_name =  "3_28")
plot(p3_28[[1]])
p3_34 <- fit_curve(id_name =  "3_34")
plot(p3_34[[1]])
p3_40 <- fit_curve(id_name =  "3_40")
plot(p3_40[[1]])
p3_46 <- fit_curve(id_name =  "3_46")
plot(p3_46[[1]]) # refitting the curves with some adjustments

data_3_46 <- filter(final_aci_data, id == "3_46")
data_3_46 <- filter(data_3_46, !(Ci >= 350 & Ci <= 600)) 
# removing points that represents the depression in the curve

curve_p3_46 <- fitaci(data_3_46,varnames = list(ALEAF = "A",
                                          Tleaf = "Tleaf",
                                          Ci = "Ci",
                                          PPFD = "Qin"),
                fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(curve_p3_46)

coef_3_46 <- data.frame(id = "3_46", t(coef(curve_p3_46)))

p3_46 <- list(curve_p3_46, coef_3_46)

coef_3 <- rbind(p3_16[[2]], p3_22[[2]], p3_28[[2]], p3_34[[2]], p3_40[[2]], p3_46[[2]])
coef_3$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_3$id))
coef_3$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_3$id)
coef_3 <- coef_3 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_3 <- nlsLM(medlyn_Vcmax, data = coef_3,
                   start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_3)
# note: T opt is within the measured temperature range
t_vcmax_para_3 <- pred_nls(t_model_3, coef_3)

# temperature response curve for Jmax

t_model_Jmax_3 <- nlsLM(medlyn_Jmax, data = coef_3,
                        start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_3)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_3 <- pred_nls(t_model_Jmax_3, coef_3)

# Plant 4
p4_16 <- fit_curve(id_name =  "4_16_re")
plot(p4_16[[1]]) # have no Jmax limitation

data_4_16 <- filter(final_aci_data, id == "4_16_re")

curve_p4_16 <- fitaci(data_4_16, varnames = list(ALEAF = "A",
                                                Tleaf = "Tleaf",
                                                Ci = "Ci",
                                                PPFD = "Qin"),
                      citransition = 300,
                      fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(curve_p4_16)

coef_4_16 <- data.frame(id = "4_16", t(coef(curve_p4_16)))

p4_16 <- list(curve_p4_16, coef_4_16)

p4_22 <- fit_curve(id_name =  "4_22_re")
plot(p4_22[[1]])
p4_28 <- fit_curve(id_name =  "4_28")
plot(p4_28[[1]])
p4_34 <- fit_curve(id_name =  "4_34")
plot(p4_34[[1]])
p4_40 <- fit_curve(id_name =  "4_40")
plot(p4_40[[1]])
p4_46 <- fit_curve(id_name =  "4_46")
plot(p4_46[[1]])

coef_4 <- rbind(p4_16[[2]], p4_22[[2]], p4_28[[2]], p4_34[[2]], p4_40[[2]], p4_46[[2]])
coef_4[2,1] <- "4_22"

coef_4$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_4$id))
coef_4$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_4$id)
coef_4 <- coef_4 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_4 <- nlsLM(medlyn_Vcmax, data = coef_4,
                   start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_4)
# note: T opt is within the measured temperature range
t_vcmax_para_4 <- pred_nls(t_model_4, coef_4)

# temperature response curve for Jmax

t_model_Jmax_4 <- nlsLM(medlyn_Jmax, data = coef_4,
                        start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_4)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_4 <- pred_nls(t_model_Jmax_4, coef_4)

# Plant 9
# fit_curve(id_name =  "9_16")
# fit_curve(id_name =  "9_22")
# fit_curve(id_name =  "9_28")
# fit_curve(id_name =  "9_34")
# fit_curve(id_name =  "9_40")
# fit_curve(id_name =  "9_46")

# Plant 10
p10_16 <- fit_curve(id_name =  "10_16")
plot(p10_16[[1]])
p10_22 <- fit_curve(id_name =  "10_22")
plot(p10_22[[1]])
p10_28 <- fit_curve(id_name =  "10_28")
plot(p10_28[[1]])
p10_34 <- fit_curve(id_name =  "10_34")
plot(p10_34[[1]])
p10_40 <- fit_curve(id_name =  "10_40")
plot(p10_40[[1]])
p10_46 <- fit_curve(id_name =  "10_46")
plot(p10_46[[1]])

coef_10 <- rbind(p10_16[[2]], p10_22[[2]], p10_28[[2]], p10_34[[2]], p10_40[[2]], p10_46[[2]])

# comments: Anet is low 
# vcmax values are also low

coef_10$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_10$id))
coef_10$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_10$id)
coef_10 <- coef_10 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_10 <- nlsLM(medlyn_Vcmax, data = coef_10,
                   start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_10)
# note: T opt is within the measured temperature range
t_vcmax_para_10 <- pred_nls(t_model_10, coef_10)

# temperature response curve for Jmax

t_model_Jmax_10 <- nlsLM(medlyn_Jmax, data = coef_10,
                        start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_10)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_10 <- pred_nls(t_model_Jmax_10, coef_10)

# Plant 11
p11_16 <- fit_curve(id_name =  "11_16")
plot(p11_16[[1]])
p11_22 <- fit_curve(id_name =  "11_22")
plot(p11_22[[1]])
p11_28 <- fit_curve(id_name =  "11_28")
plot(p11_28[[1]])
p11_34 <- fit_curve(id_name =  "11_34")
plot(p11_34[[1]])
p11_40 <- fit_curve(id_name =  "11_40")
plot(p11_40[[1]])
# p11_46 <- fit_curve(id_name =  "11_46")
# plot(p11_46[[1]]) # refitting this curve with some modifications


data_11_46 <- filter(final_aci_data, id == "11_46")
data_11_46 <- filter(data_11_46, !(Ci >= 320 & Ci <= 500)) 
# removing points that represents the depression in the curve

curve_p11_46 <- fitaci(data_11_46,varnames = list(ALEAF = "A",
                                                Tleaf = "Tleaf",
                                                Ci = "Ci",
                                                PPFD = "Qin"),
                       citransition = 320,
                      fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(curve_p11_46)

coef_11_46 <- data.frame(id = "11_46", t(coef(curve_p11_46)))

p11_46 <- list(curve_p11_46, coef_11_46)

coef_11 <- rbind(p11_16[[2]], p11_22[[2]], p11_28[[2]], p11_34[[2]], p11_40[[2]], p11_46[[2]])

# Note: Vcmax low and still going up, probably Topt will be very high

coef_11$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_11$id))
coef_11$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_11$id)
coef_11 <- coef_11 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_11 <- nlsLM(medlyn_Vcmax, data = coef_11,
                   start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_11)
# note: T opt is within the measured temperature range
t_vcmax_para_11 <- pred_nls(t_model_11, coef_11)

# temperature response curve for Jmax

t_model_Jmax_11 <- nlsLM(medlyn_Jmax, data = coef_11,
                        start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_11)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_11 <- pred_nls(t_model_Jmax_11, coef_11)

# Plant 12

p12_16 <- fit_curve(id_name =  "12_16", citrans = 450)
plot(p12_16[[1]])
p12_22 <- fit_curve(id_name =  "12_22")
plot(p12_22[[1]])
p12_28 <- fit_curve(id_name =  "12_28")
plot(p12_28[[1]])
p12_34 <- fit_curve(id_name =  "12_34")
plot(p12_34[[1]])
p12_40 <- fit_curve(id_name =  "12_40")
plot(p12_40[[1]])
# fit_curve(id_name =  "12_46")

data_12_46 <- filter(final_aci_data, id == "12_46")
data_12_46 <- data_12_46[1:96,]

curve_p12_46 <- fitaci(data_12_46,varnames = list(ALEAF = "A",
                                    Tleaf = "Tleaf",
                                    Ci = "Ci",
                                    PPFD = "Qin"),
                      citransition = 500,
                      fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)

plot(curve_p12_46)
coef_12_46 <- data.frame(id = "12_46", t(coef(curve_p12_46)))

p12_46 <- list(curve_p12_46, coef_12_46)

coef_12 <- rbind(p12_16[[2]], p12_22[[2]], p12_28[[2]], p12_34[[2]], p12_40[[2]], p12_46[[2]])

coef_12$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_12$id))
coef_12$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_12$id)
coef_12 <- coef_12 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_12 <- nlsLM(medlyn_Vcmax, data = coef_12,
                   start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_12)
# note: T opt is within the measured temperature range
t_vcmax_para_12 <- pred_nls(t_model_12, coef_12)

# temperature response curve for Jmax

t_model_Jmax_12 <- nlsLM(medlyn_Jmax, data = coef_12,
                        start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_12)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_12 <- pred_nls(t_model_Jmax_12, coef_12)

# Plant 17
p17_16 <- fit_curve(id_name =  "17_16", citrans = 500)
plot(p17_16[[1]])
p17_22 <- fit_curve(id_name =  "17_22")
plot(p17_22[[1]])
p17_28 <- fit_curve(id_name =  "17_28")
plot(p17_28[[1]])
p17_34 <- fit_curve(id_name =  "17_34")
plot(p17_34[[1]])
p17_40 <- fit_curve(id_name =  "17_40")
plot(p17_40[[1]])
# p17_46 <- fit_curve(id_name =  "17_46")
# plot(p17_46[[1]])
# don't have a good curve

coef_17 <- rbind(p17_16[[2]], p17_22[[2]], p17_28[[2]], p17_34[[2]], p17_40[[2]])

# Vcmax, Jmax too low. Probably will exclude for analysis.

coef_17$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_17$id))
coef_17$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_17$id)
coef_17 <- coef_17 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_17 <- nlsLM(medlyn_Vcmax, data = coef_17,
                   start = list(kopt = 50, Ha = 115000, Topt = 360))

# Print the summary of the model
summary(t_model_17)
# note: didn't fit well



# Plant 18
p18_16 <- fit_curve(id_name =  "18_16", citrans = 450)
plot(p18_16[[1]])

p18_22 <- fit_curve(id_name =  "18_22")
plot(p18_22[[1]])

p18_28 <- fit_curve(id_name =  "18_28")
plot(p18_28[[1]])

p18_34 <- fit_curve(id_name =  "18_34")
plot(p18_34[[1]])

p18_40 <- fit_curve(id_name =  "18_40", citrans = 400)
plot(p18_40[[1]])

# p18_46 <- fit_curve(id_name =  "18_46re")
# plot(p18_46[[1]])
# refitting curve with some modification in data

data_18_46 <- filter(final_aci_data, id == "18_46re")
data_18_46 <- filter(data_18_46, !(Ci >= 360 & Ci <= 600)) 
# removing points that represents the depression in the curve

curve_p18_46 <- fitaci(data_18_46,varnames = list(ALEAF = "A",
                                                  Tleaf = "Tleaf",
                                                  Ci = "Ci",
                                                  PPFD = "Qin"),
                       fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(curve_p18_46)

coef_18_46 <- data.frame(id = "18_46", t(coef(curve_p18_46)))

p18_46 <- list(curve_p18_46, coef_18_46)


coef_18 <- rbind(p18_16[[2]], p18_22[[2]], p18_28[[2]], p18_34[[2]], p18_40[[2]], p18_46[[2]])

coef_18$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_18$id))
coef_18$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_18$id)
coef_18 <- coef_18 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_18 <- nlsLM(medlyn_Vcmax, data = coef_18,
                   start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_18)
# note: T opt is within the measured temperature range
t_vcmax_para_18 <- pred_nls(t_model_18, coef_18)

# temperature response curve for Jmax

t_model_Jmax_18 <- nlsLM(medlyn_Jmax, data = coef_18,
                        start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_18)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_18 <- pred_nls(t_model_Jmax_18, coef_18)

# Plant 19
p19_16 <- fit_curve(id_name =  "19_16")
plot(p19_16[[1]])

p19_22 <- fit_curve(id_name =  "19_22")
plot(p19_22[[1]])

p19_28 <- fit_curve(id_name =  "19_28")
plot(p19_28[[1]])

p19_34 <- fit_curve(id_name =  "19_34")
plot(p19_34[[1]])

p19_40 <- fit_curve(id_name =  "19_40")
plot(p19_40[[1]])

# p19_46 <- fit_curve(id_name =  "19_46")
# plot(p19_46[[1]])

coef_19 <- rbind(p19_16[[2]], p19_22[[2]], p19_28[[2]], p19_34[[2]], p19_40[[2]])

# Note: very low Vcmax and Jmax, probably exclude from further analysis

coef_19$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_19$id))
coef_19$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_19$id)
coef_19 <- coef_19 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_19 <- nlsLM(medlyn_Vcmax, data = coef_19,
                   start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_19)
# note: T opt is within the measured temperature range
# did not fit well

# temperature response curve for Jmax

t_model_Jmax_19 <- nlsLM(medlyn_Jmax, data = coef_19,
                         start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_19)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_19 <- pred_nls(t_model_Jmax_19, coef_19)

# Plant 20
p20_16 <- fit_curve(id_name =  "20_16re", citrans = 500)
plot(p20_16[[1]])

p20_22 <- fit_curve(id_name =  "20_22", citrans = 400)
plot(p20_22[[1]])

# p20_28 <- fit_curve(id_name =  "20_28")

data_20_28 <- filter(final_aci_data, id == "20_28")
data_20_28 <- data_20_28[97:192,]

curve_p20_28 <- fitaci(data_20_28,varnames = list(ALEAF = "A",
                                    Tleaf = "Tleaf",
                                    Ci = "Ci",
                                    PPFD = "Qin"),
                  fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)

plot(curve_p20_28)

coef_20_28 <- data.frame(id = "20_28", t(coef(curve_p20_28)))

p20_28 <- list(curve_p20_28, coef_20_28)

p20_34 <- fit_curve(id_name =  "20_34", citrans = 350)
plot(p20_34[[1]])

p20_40 <- fit_curve(id_name =  "20_40")
plot(p20_40[[1]])

# p20_46 <- fit_curve(id_name =  "20_46")
# plot(p20_46[[1]])
# refitting curve with some adjustments

data_20_46 <- filter(final_aci_data, id == "20_46")
data_20_46 <- filter(data_20_46, !(Ci >= 340 & Ci <= 500)) 
# removing points that represents the depression in the curve

curve_p20_46 <- fitaci(data_20_46,varnames = list(ALEAF = "A",
                                                  Tleaf = "Tleaf",
                                                  Ci = "Ci",
                                                  PPFD = "Qin"),
                       fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(curve_p20_46)

coef_20_46 <- data.frame(id = "20_46", t(coef(curve_p20_46)))

p20_46 <- list(curve_p20_46, coef_20_46)


coef_20 <- rbind(p20_16[[2]], p20_22[[2]], p20_28[[2]], p20_34[[2]], p20_40[[2]], p20_46[[2]])
coef_20[1,1] <- "20_16"
coef_20$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_20$id))
coef_20$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_20$id)
coef_20 <- coef_20 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_20 <- nlsLM(medlyn_Vcmax, data = coef_20,
                   start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_20)
# note: T opt is within the measured temperature range

t_vcmax_para_20 <- pred_nls(t_model_20, coef_20)

# temperature response curve for Jmax

t_model_Jmax_20 <- nlsLM(medlyn_Jmax, data = coef_20,
                        start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_20)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_20 <- pred_nls(t_model_Jmax_20, coef_20)

# Plant 25
p25_16 <- fit_curve(id_name =  "25_16")
plot(p25_16[[1]])

p25_22 <- fit_curve(id_name =  "25_22")
plot(p25_22[[1]])

p25_28 <- fit_curve(id_name =  "25_28")
plot(p25_28[[1]])

p25_34 <- fit_curve(id_name =  "25_34")
plot(p25_34[[1]])

p25_40 <- fit_curve(id_name =  "25_40")
plot(p25_40[[1]])

# p25_46 <- fit_curve(id_name =  "25_46")
# plot(p25_46[[1]])

coef_25 <- rbind(p25_16[[2]], p25_22[[2]], p25_28[[2]], p25_34[[2]], p25_40[[2]])
coef_25$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_25$id))
coef_25$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_25$id)
coef_25 <- coef_25 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_25 <- nlsLM(medlyn_Vcmax, data = coef_25,
                   start = list(kopt = 200, Ha = 115000, Topt = 360))

# Print the summary of the model
summary(t_model_25)
# note: T opt is outside the the measured temperature range

# temperature response curve for Jmax

t_model_Jmax_25 <- nlsLM(medlyn_Jmax, data = coef_25,
                        start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_25)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_25 <- pred_nls(t_model_Jmax_25, coef_25)

# Plant 26
p26_16 <- fit_curve(id_name =  "26_16", citrans = 400)
plot(p26_16[[1]])

p26_22 <- fit_curve(id_name =  "26_22")
plot(p26_22[[1]])

p26_28 <- fit_curve(id_name =  "26_28")
plot(p26_28[[1]])

p26_34 <- fit_curve(id_name =  "26_34")
plot(p26_34[[1]])

p26_40 <- fit_curve(id_name =  "26_40", citrans = 300)
plot(p26_40[[1]])

# p26_46 <- fit_curve(id_name =  "26_46")
# plot(p26_46[[1]])

coef_26 <- rbind(p26_16[[2]], p26_22[[2]], p26_28[[2]], p26_34[[2]], p26_40[[2]])

coef_26$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_26$id))
coef_26$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_26$id)
coef_26 <- coef_26 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_26 <- nlsLM(medlyn_Vcmax, data = coef_26,
                   start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_26)
# note: T opt is outside the measured temperature range

# temperature response curve for Jmax

t_model_Jmax_26 <- nlsLM(medlyn_Jmax, data = coef_26,
                        start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_26)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_26 <- pred_nls(t_model_Jmax_26, coef_26)

# Plant 27
p27_16 <- fit_curve(id_name =  "27_16")
plot(p27_16[[1]])

p27_22 <- fit_curve(id_name =  "27_22")
plot(p27_22[[1]])

p27_28 <- fit_curve(id_name =  "27_28")
plot(p27_28[[1]])

p27_34 <- fit_curve(id_name =  "27_34")
plot(p27_34[[1]])

p27_40 <- fit_curve(id_name =  "27_40", citrans = 250)
plot(p27_40[[1]])

# p27_46 <- fit_curve(id_name =  "27_46")
# plot(p27_46[[1]])
# refiitng the curve with some modifications

data_27_46 <- filter(final_aci_data, id == "27_46")
data_27_46 <- filter(data_27_46, !(Ci >= 340 & Ci <= 600)) 
# removing points that represents the depression in the curve

curve_p27_46 <- fitaci(data_27_46,varnames = list(ALEAF = "A",
                                                  Tleaf = "Tleaf",
                                                  Ci = "Ci",
                                                  PPFD = "Qin"),
                       citransition = 400,
                       fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(curve_p27_46)

coef_27_46 <- data.frame(id = "27_46", t(coef(curve_p27_46)))

p27_46 <- list(curve_p27_46, coef_27_46)

coef_27 <- rbind(p27_16[[2]], p27_22[[2]], p27_28[[2]], p27_34[[2]], p27_40[[2]], p27_46[[2]])

coef_27$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_27$id))
coef_27$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_27$id)
coef_27 <- coef_27 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_27 <- nlsLM(medlyn_Vcmax, data = coef_27,
                   start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_27)
# note: T opt is within the measured temperature range
t_vcmax_para_27 <- pred_nls(t_model_27, coef_27)

# temperature response curve for Jmax

t_model_Jmax_27 <- nlsLM(medlyn_Jmax, data = coef_27,
                        start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_27)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_27 <- pred_nls(t_model_Jmax_27, coef_27)


# Plant 28
p28_16 <- fit_curve(id_name =  "28_16")
plot(p28_16[[1]])

p28_22 <- fit_curve(id_name =  "28_22")
plot(p28_22[[1]])

p28_28 <- fit_curve(id_name =  "28_28")
plot(p28_28[[1]])

p28_34 <- fit_curve(id_name =  "28_34")
plot(p28_34[[1]])

p28_40 <- fit_curve(id_name =  "28_40")
plot(p28_40[[1]])

# p28_46 <- fit_curve(id_name =  "28_46")
# plot(p28_46[[1]])

data_28_46 <- filter(final_aci_data, id == "28_46")
data_28_46 <- filter(data_28_46, !(Ci >= 340 & Ci <= 550)) 
# removing points that represents the depression in the curve

curve_p28_46 <- fitaci(data_28_46,varnames = list(ALEAF = "A",
                                                  Tleaf = "Tleaf",
                                                  Ci = "Ci",
                                                  PPFD = "Qin"),
                       citransition = 450,
                       fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(curve_p28_46)

coef_28_46 <- data.frame(id = "28_46", t(coef(curve_p28_46)))

p28_46 <- list(curve_p28_46, coef_28_46)

coef_28 <- rbind(p28_16[[2]], p28_22[[2]], p28_28[[2]], p28_34[[2]], p28_40[[2]], p28_46[[2]])
coef_28$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_28$id))
coef_28$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_28$id)
coef_28 <- coef_28 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_28 <- nlsLM(medlyn_Vcmax, data = coef_28,
                   start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_28)
# note: T opt is within the measured temperature range
t_vcmax_para_28 <- pred_nls(t_model_28, coef_28)

# temperature response curve for Jmax

t_model_Jmax_28 <- nlsLM(medlyn_Jmax, data = coef_28,
                        start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_28)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_28 <- pred_nls(t_model_Jmax_28, coef_28)


# plant 33 to 38
# read data

final_33_38_data <- read.csv("~/Desktop/licor data/final_33_38_data.csv")

# create a column for temp setpoint
final_33_38_data$temp_setpoint <- round(final_33_38_data$Tleaf, digits = 0)
unique(final_33_38_data$temp_setpoint)
# it has 45, which was not actual set point
# adjusting temp setpoint for values differnt from the actual setpoint

final_33_38_data$temp_setpoint <- 
  ifelse(final_33_38_data$temp_setpoint == 45, 46, final_33_38_data$temp_setpoint)
unique(final_33_38_data$temp_setpoint)

final_33_38_data <- final_33_38_data %>%
  select(-Tleaf) %>%
  dplyr::rename(Tleaf = temp_setpoint) %>%
  arrange(id)

final_aci_data <- final_aci_data %>%
  select(-rd)
final_aci_data <- rbind(final_aci_data, final_33_38_data)


# now fit the aci curves

# plant 33

p33_16 <- fit_curve(id_name =  "33_16", citrans = 400)
plot(p33_16[[1]])

p33_22 <- fit_curve(id_name =  "33_22")
plot(p33_22[[1]])

p33_28 <- fit_curve(id_name =  "33_28")
plot(p33_28[[1]])

p33_34 <- fit_curve(id_name =  "33_34")
plot(p33_34[[1]])

p33_40 <- fit_curve(id_name =  "33_40")
plot(p33_40[[1]])

p33_46 <- fit_curve(id_name =  "33_46")
plot(p33_46[[1]])

coef_33 <- rbind(p33_16[[2]], p33_22[[2]], p33_28[[2]], p33_34[[2]], p33_40[[2]], p33_46[[2]])
coef_33$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_33$id))
coef_33$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_33$id)
coef_33 <- coef_33 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_33 <- nlsLM(medlyn_Vcmax, data = coef_33,
                    start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_33)
# note: T opt is within the measured temperature range
t_vcmax_para_33 <- pred_nls(t_model_33, coef_33)

# temperature response curve for Jmax

t_model_Jmax_33 <- nlsLM(medlyn_Jmax, data = coef_33,
                        start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_33)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_33 <- pred_nls(t_model_Jmax_33, coef_33)

# plant 34

p34_16 <- fit_curve(id_name =  "34_16")
plot(p34_16[[1]])

p34_22 <- fit_curve(id_name =  "34_22")
plot(p34_22[[1]])

p34_28 <- fit_curve(id_name =  "34_28")
plot(p34_28[[1]])

p34_34 <- fit_curve(id_name =  "34_34")
plot(p34_34[[1]])

p34_40 <- fit_curve(id_name =  "34_40")
plot(p34_40[[1]])

p34_46 <- fit_curve(id_name =  "34_46")
plot(p34_46[[1]])

coef_34 <- rbind(p34_16[[2]], p34_22[[2]], p34_28[[2]], p34_34[[2]], p34_40[[2]], p34_46[[2]])
coef_34$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_34$id))
coef_34$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_34$id)
coef_34 <- coef_34 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_34 <- nlsLM(medlyn_Vcmax, data = coef_34,
                    start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_34)
# note: T opt is within the measured temperature range
t_vcmax_para_34 <- pred_nls(t_model_34, coef_34)

# temperature response curve for Jmax

# encounter an error
# t_model_Jmax_34 <- nlsLM(medlyn_Jmax, data = coef_34,
#                         start = list(kopt = 10, Ha = 115000, Topt = 300))
# 
# # Print the summary of the model
# summary(t_model_Jmax_1)
# # note: T opt is within the measured temperature range
# # extracting parameters
# t_Jmax_para_1 <- pred_nls(t_model_Jmax_1, coef_1)

# plant 35

p35_16 <- fit_curve(id_name =  "35_16", citrans = 400)
plot(p35_16[[1]])

p35_22 <- fit_curve(id_name =  "35_22")
plot(p35_22[[1]])

p35_28 <- fit_curve(id_name =  "35_28")
plot(p33_28[[1]])

p35_34 <- fit_curve(id_name =  "35_34")
plot(p33_34[[1]])

# p35_40 <- fit_curve(id_name =  "35_40")
# plot(p35_40[[1]])
# refitting with some adjustments

data_35_40 <- filter(final_aci_data, id == "35_40")
data_35_40 <- filter(data_35_40, !(Ci >= 350 & Ci <= 500)) 
# removing points that represents the depression in the curve

curve_p35_40 <- fitaci(data_35_40,varnames = list(ALEAF = "A",
                                                  Tleaf = "Tleaf",
                                                  Ci = "Ci",
                                                  PPFD = "Qin"),
                       citransition = 200,
                       fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(curve_p35_40)

coef_35_40 <- data.frame(id = "35_46", t(coef(curve_p35_40)))

p35_40 <- list(curve_p35_40, coef_35_40)


# p35_46 <- fit_curve(id_name =  "35_46")
# plot(p35_46[[1]])
# refitting with modifications


data_35_46 <- filter(final_aci_data, id == "35_46")
data_35_46 <- filter(data_35_46, !(Ci >= 400 & Ci <= 500)) 
# removing points that represents the depression in the curve

curve_p35_46 <- fitaci(data_35_46,varnames = list(ALEAF = "A",
                                                  Tleaf = "Tleaf",
                                                  Ci = "Ci",
                                                  PPFD = "Qin"),
                       fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(curve_p35_46)

coef_35_46 <- data.frame(id = "35_40", t(coef(curve_p35_46)))

p35_46 <- list(curve_p35_46, coef_35_46)

coef_35 <- rbind(p35_16[[2]], p35_22[[2]], p35_28[[2]], p35_34[[2]], p35_46[[2]], p35_40[[2]])
coef_35$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_35$id))
coef_35$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_35$id)
coef_35 <- coef_35 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_35 <- nlsLM(medlyn_Vcmax, data = coef_35,
                    start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_35)
# note: T opt is within the measured temperature range
t_vcmax_para_35 <- pred_nls(t_model_35, coef_35)

# temperature response curve for Jmax

t_model_Jmax_35 <- nlsLM(medlyn_Jmax, data = coef_35,
                         start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_35)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_35 <- pred_nls(t_model_Jmax_35, coef_35)

# plant 36

p36_16 <- fit_curve(id_name =  "36_16", citrans = 400)
plot(p36_16[[1]])

p36_22 <- fit_curve(id_name =  "36_22")
plot(p36_22[[1]])

p36_28 <- fit_curve(id_name =  "36_28")
plot(p36_28[[1]])

p36_34 <- fit_curve(id_name =  "36_34")
plot(p36_34[[1]])

p36_40 <- fit_curve(id_name =  "36_40")
plot(p36_40[[1]])

p36_46 <- fit_curve(id_name =  "36_46")
plot(p36_46[[1]])

coef_36 <- rbind(p36_16[[2]], p36_22[[2]], p36_28[[2]], p36_34[[2]], p36_40[[2]], p36_46[[2]])
coef_36$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_36$id))
coef_36$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_36$id)
coef_36 <- coef_36 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_36 <- nlsLM(medlyn_Vcmax, data = coef_36,
                    start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_36)
# note: T opt is within the measured temperature range
t_vcmax_para_36 <- pred_nls(t_model_36, coef_36)

# temperature response curve for Jmax

t_model_Jmax_36 <- nlsLM(medlyn_Jmax, data = coef_36,
                         start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_36)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_36 <- pred_nls(t_model_Jmax_36, coef_36)

# plant 37

p37_16 <- fit_curve(id_name =  "37_16", citrans = 400)
plot(p37_16[[1]])

p37_22 <- fit_curve(id_name =  "37_22")
plot(p37_22[[1]])

p37_28 <- fit_curve(id_name =  "37_28")
plot(p37_28[[1]])

p37_34 <- fit_curve(id_name =  "37_34")
plot(p37_34[[1]])

p37_40 <- fit_curve(id_name =  "37_40")
plot(p37_40[[1]])

# p37_46 <- fit_curve(id_name =  "37_46")
# plot(p37_46[[1]])
# refitting with some adjustments

data_37_46 <- filter(final_aci_data, id == "37_46")
data_37_46 <- filter(data_37_46, !(Ci >= 350 & Ci <= 550)) 
# removing points that represents the depression in the curve

curve_p37_46 <- fitaci(data_37_46,varnames = list(ALEAF = "A",
                                                  Tleaf = "Tleaf",
                                                  Ci = "Ci",
                                                  PPFD = "Qin"),
                       fitTPU = TRUE, Tcorrect = FALSE, useRd = FALSE)
plot(curve_p37_46)

coef_37_46 <- data.frame(id = "37_46", t(coef(curve_p37_46)))

p37_46 <- list(curve_p37_46, coef_37_46)


coef_37 <- rbind(p37_16[[2]], p37_22[[2]], p37_28[[2]], p37_34[[2]], p37_40[[2]], p37_46[[2]])
coef_37$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_37$id))
coef_37$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_37$id)
coef_37 <- coef_37 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM

t_model_37 <- nlsLM(medlyn_Vcmax, data = coef_37,
                    start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_37)
# note: T opt is within the measured temperature range
t_vcmax_para_37 <- pred_nls(t_model_37, coef_37)

# temperature response curve for Jmax

t_model_Jmax_37 <- nlsLM(medlyn_Jmax, data = coef_37,
                         start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_37)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_37 <- pred_nls(t_model_Jmax_37, coef_37)

# plant 38

p38_16 <- fit_curve(id_name =  "38_16", citrans = 400)
plot(p38_16[[1]])

p38_22 <- fit_curve(id_name =  "38_22")
plot(p38_22[[1]])

p38_28 <- fit_curve(id_name =  "38_28")
plot(p38_28[[1]])

p38_34 <- fit_curve(id_name =  "38_34")
plot(p38_34[[1]])

p38_40 <- fit_curve(id_name =  "38_40", citrans = 300)
plot(p38_40[[1]])

p38_46 <- fit_curve(id_name =  "38_46", citrans = 600)
plot(p38_46[[1]])

coef_38 <- rbind(p38_16[[2]], p38_22[[2]], p38_28[[2]], p38_34[[2]], p38_40[[2]], p38_46[[2]])
coef_38$Tleaf <- as.integer(gsub(".*_(\\d+)$", "\\1", coef_38$id))
coef_38$plant_id <- gsub("^(\\d+)_.*", "\\1", coef_38$id)
coef_38 <- coef_38 %>%
  mutate(Tleaf = Tleaf + 273.15)

# temperature response curve for Vcmax
# Fit the model using nlsLM
# encounter an error 

# t_model_38 <- nlsLM(medlyn_Vcmax, data = coef_38,
#                     start = list(kopt = 50, Ha = 115000, Topt = 310))
# 
# # Print the summary of the model
# summary(t_model_38)
# # note: T opt is beyond the measured temperature range
# t_vcmax_para_38 <- pred_nls(t_model_38, coef_38)

# temperature response curve for Jmax

t_model_Jmax_33 <- nlsLM(medlyn_Jmax, data = coef_33,
                         start = list(kopt = 50, Ha = 115000, Topt = 310))

# Print the summary of the model
summary(t_model_Jmax_33)
# note: T opt is within the measured temperature range
# extracting parameters
t_Jmax_para_33 <- pred_nls(t_model_Jmax_33, coef_33)

###########################################################

# compile all the temperature response paramters for Vcmax

t_vcmax_para <- rbind(t_vcmax_para_1[[1]], t_vcmax_para_2[[1]], t_vcmax_para_3[[1]],
                           t_vcmax_para_4[[1]], t_vcmax_para_10[[1]], t_vcmax_para_11[[1]], 
                           t_vcmax_para_12[[1]], t_vcmax_para_18[[1]], t_vcmax_para_20[[1]],
                           t_vcmax_para_27[[1]], t_vcmax_para_28[[1]], t_vcmax_para_33[[1]],
                           t_vcmax_para_34[[1]], t_vcmax_para_35[[1]], t_vcmax_para_36[[1]], 
                           t_vcmax_para_37[[1]])

t_vcmax_predicted <- rbind(t_vcmax_para_1[[2]], t_vcmax_para_2[[2]], t_vcmax_para_3[[2]],
                      t_vcmax_para_4[[2]], t_vcmax_para_10[[2]], t_vcmax_para_11[[2]], 
                      t_vcmax_para_12[[2]], t_vcmax_para_18[[2]], t_vcmax_para_20[[2]],
                      t_vcmax_para_27[[2]], t_vcmax_para_28[[2]], t_vcmax_para_33[[2]],
                      t_vcmax_para_34[[2]], t_vcmax_para_35[[2]], t_vcmax_para_36[[2]], 
                      t_vcmax_para_37[[2]])

t_vcmax_predicted$Tleaf <- gsub(".*_(\\d+)$", "\\1", t_vcmax_predicted$id)
t_vcmax_predicted$plant_id <- gsub("^(\\d+)_.*", "\\1", t_vcmax_predicted$id)

# reading meta data
meta <- read.csv("~/Desktop/licor data/metadata.csv")

# combine with t_vcmax_predicted
t_vcmax_predicted <- merge(t_vcmax_predicted, meta, by = "plant_id")

# combine with t_vcmax_para

t_vcmax_para <- merge(t_vcmax_para, meta, by = "plant_id")

##############################
# Kopt model
Kopt_vcmax <- t_vcmax_para %>%
  filter(parameters == "kopt")

Vcmax_kopt_model <- lm(Estimate ~ growth_temp * irrigation_trt, data = Kopt_vcmax)

Anova(Vcmax_kopt_model, type = 3)
summary(Vcmax_kopt_model)
plot(Vcmax_kopt_model)
hist(residuals(Vcmax_kopt_model))
shapiro.test(residuals(Vcmax_kopt_model)) # normally distributed

pairs(emmeans(Vcmax_kopt_model, ~ growth_temp:irrigation_trt))

mean_kopt <- as.data.frame(emmeans(Vcmax_kopt_model, ~ growth_temp:irrigation_trt))

ggplot(mean_kopt, aes(growth_temp, emmean, color = irrigation_trt)) +
  geom_point(size = 5, position = position_dodge(width = 1), alpha = 1) +
  geom_point(data = Kopt_vcmax, aes(x = growth_temp, y = Estimate, 
                                    color = irrigation_trt), size = 2.5,
             alpha = 0.5,
             position = position_jitterdodge()) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE, color = irrigation_trt),
                width = 0.2, linewidth = 1, 
                position = position_dodge(width = 1)) +
  scale_color_manual(values = c("brown", "blue")) +
  labs(x = NULL, y = " Kopt") +
  
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank(),
        legend.position = c(0.15, 0.9))

#############################

# Topt model
Topt_vcmax <- t_vcmax_para %>%
  filter(parameters == "Topt")

Vcmax_topt_model <- lm(Estimate ~ growth_temp * irrigation_trt, data = Topt_vcmax)

Anova(Vcmax_topt_model, type = 3)
summary(Vcmax_topt_model)
hist(residuals(Vcmax_topt_model))
shapiro.test(residuals(Vcmax_topt_model)) # normal distribution

pairs(emmeans(Vcmax_topt_model, ~ growth_temp:irrigation_trt))

mean_topt <- as.data.frame(emmeans(Vcmax_topt_model, ~ growth_temp:irrigation_trt))

ggplot(mean_topt, aes(growth_temp, emmean, color = irrigation_trt)) +
  geom_point(size = 5, position = position_dodge(width = 1), alpha = 1) +
  geom_point(data = Topt_vcmax, aes(x = growth_temp, y = Estimate, 
                                     color = irrigation_trt), size = 2.5,
              alpha = 0.5,
              position = position_jitterdodge()) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE, color = irrigation_trt),
                width = 0.2, linewidth = 1, 
                position = position_dodge(width = 1)) +
  scale_color_manual(values = c("brown", "blue")) +
  labs(x = NULL, y = expression(Topt ~ of ~ Vcmax ~ (degree*C))) +
  
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.9)) +
  scale_y_continuous(labels = function (x) x - 273)



# Ha model
Ha_vcmax <- t_vcmax_para %>%
  filter(parameters == "Ha")

Vcmax_Ha_model <- lm(Estimate ~ growth_temp * irrigation_trt, data = Ha_vcmax)

Anova(Vcmax_Ha_model, type = 3)
summary(Vcmax_Ha_model)
hist(residuals(Vcmax_Ha_model))
shapiro.test(residuals(Vcmax_Ha_model)) # non normal distribution
# refit the model with glm

Vcmax_Ha_model <- glm(Estimate ~ growth_temp * irrigation_trt, data = Ha_vcmax,
                      family = inverse.gaussian(link = "1/mu^2"))

Anova(Vcmax_Ha_model, type = 3)
summary(Vcmax_Ha_model)
hist(residuals(Vcmax_Ha_model))
shapiro.test(residuals(Vcmax_Ha_model)) 

pairs(emmeans(Vcmax_Ha_model, ~ growth_temp:irrigation_trt))

mean_Ha <- as.data.frame(emmeans(Vcmax_Ha_model, ~ growth_temp:irrigation_trt), 
                         type ="response")

ggplot(mean_Ha, aes(growth_temp, response, color = irrigation_trt)) +
  geom_point(size = 5, position = position_dodge(width = 1), alpha = 1) +
  geom_point(data = Ha_vcmax, aes(x = growth_temp, y = Estimate, 
                                    color = irrigation_trt), size = 2.5,
             alpha = 0.5,
             position = position_jitterdodge()) +
  geom_errorbar(aes(ymin = response - SE, ymax = response + SE, color = irrigation_trt),
                width = 0.2, linewidth = 1, 
                position = position_dodge(width = 1)) +
  scale_color_manual(values = c("brown", "blue")) +
  labs(x = NULL, y = " Ha") +
  
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank(),
        legend.position = c(0.2, 0.8))


# create a data frame from mean Kopt, Topt, and Ha

mean_parameters <- as.data.frame(list(Kopt = mean_kopt$emmean,
                                 Ha = mean_Ha$emmean,
                                 Topt = mean_topt$emmean,
                                 treatments = c("LT+D", "HT+D",
                                                "LT+WW", "HT+WW")))

# Plotting
tleaf_seq <- seq(285, 320, 1)

t_parameter_ext_LTD <- mean_parameters[1,1] * (200000 * exp((mean_parameters[1,2] * (tleaf_seq - mean_parameters[1,3]))/(tleaf_seq * R * mean_parameters[1,3])))/(200000 - (mean_parameters[1,2] * (1 - exp((200000 *(tleaf_seq - mean_parameters[1,3]))/(tleaf_seq * R * mean_parameters[1,3])))))
t_parameter_ext_HTD <- mean_parameters[2,1] * (200000 * exp((mean_parameters[2,2] * (tleaf_seq - mean_parameters[2,3]))/(tleaf_seq * R * mean_parameters[2,3])))/(200000 - (mean_parameters[2,2] * (1 - exp((200000 *(tleaf_seq - mean_parameters[2,3]))/(tleaf_seq * R * mean_parameters[2,3])))))
t_parameter_ext_LTWW <- mean_parameters[3,1] * (200000 * exp((mean_parameters[3,2] * (tleaf_seq - mean_parameters[3,3]))/(tleaf_seq * R * mean_parameters[3,3])))/(200000 - (mean_parameters[3,2] * (1 - exp((200000 *(tleaf_seq - mean_parameters[3,3]))/(tleaf_seq * R * mean_parameters[3,3])))))
t_parameter_ext_HTWW <- mean_parameters[4,1] * (200000 * exp((mean_parameters[4,2] * (tleaf_seq - mean_parameters[4,3]))/(tleaf_seq * R * mean_parameters[4,3])))/(200000 - (mean_parameters[4,2] * (1 - exp((200000 *(tleaf_seq - mean_parameters[4,3]))/(tleaf_seq * R * mean_parameters[4,3])))))

t_parameter_ext_df <- data.frame("Tleaf" = tleaf_seq,
                                 "LT_D" = t_parameter_ext_LTD,
                                 "HT_D" = t_parameter_ext_HTD,
                                 "LT_WW" = t_parameter_ext_LTWW,
                                 "HT_WW" = t_parameter_ext_HTWW)

t_parameter_ext_df <- t_parameter_ext_df %>%
  pivot_longer(!Tleaf, names_to = "treatments", values_to = "Vcmax")

t_parameter_ext_df <- t_parameter_ext_df %>%
  mutate(Tleaf = Tleaf - 273.15)


# vcmax model
vcmax_model <- lm(average ~ Tleaf * treatment ,
                   data = t_vcmax_predicted)

anova(vcmax_model)

mean_vcmax <- as.data.frame(emmeans(vcmax_model, ~Tleaf))
mean_vcmax$tleaf <- as.numeric(as.character(mean_vcmax$Tleaf))


(vcmax_graph <- ggplot(mean_vcmax, aes(tleaf, emmean)) +
  geom_point(size = 4) +
  geom_line(data = t_parameter_ext_df,
            aes(x = Tleaf, y = Vcmax, color = treatments), linewidth = 1) +
  theme_bw(base_size = 16) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.5, linewidth = 1) +
  scale_x_continuous(limits = c(10,50),
                     breaks = seq(0,50,5)) +
  labs(x = expression(Leaf ~ Temperature ~ (degree*C)),
       y = expression(Vcmax ~ (mu*mol ~ m^{-2} ~ s^{-1}))) +
  
  scale_color_manual(values = c("turquoise", "orange","magenta", "brown" ),
                     label = c("High Temperature, Drought", "High Temperature, Irrigated",
                                       "Ambient Temperature, Drought", "Ambient Temperature, Irrigated")) +
  theme(legend.title = element_blank(),
        legend.position = c(0.15, 0.84),
        panel.grid = element_blank()))


##############################################################################

# compile all the temperature response parameters for Jmax

t_Jmax_para <- rbind(t_Jmax_para_1[[1]], t_Jmax_para_2[[1]], t_Jmax_para_3[[1]],
                      t_Jmax_para_4[[1]], t_Jmax_para_10[[1]], t_Jmax_para_11[[1]], 
                      t_Jmax_para_12[[1]], t_Jmax_para_18[[1]], t_Jmax_para_20[[1]],
                      t_Jmax_para_27[[1]], t_Jmax_para_28[[1]], t_Jmax_para_33[[1]],
                      t_Jmax_para_19[[1]], t_Jmax_para_35[[1]], t_Jmax_para_36[[1]], 
                      t_Jmax_para_37[[1]])

t_Jmax_predicted <- rbind(t_Jmax_para_1[[2]], t_Jmax_para_2[[2]], t_Jmax_para_3[[2]],
                           t_Jmax_para_4[[2]], t_Jmax_para_10[[2]], t_Jmax_para_11[[2]], 
                           t_Jmax_para_12[[2]], t_Jmax_para_18[[2]], t_Jmax_para_20[[2]],
                           t_Jmax_para_27[[2]], t_Jmax_para_28[[2]], t_Jmax_para_33[[2]],
                           t_Jmax_para_19[[2]], t_Jmax_para_35[[2]], t_Jmax_para_36[[2]], 
                           t_Jmax_para_37[[2]])

t_Jmax_predicted$Tleaf <- gsub(".*_(\\d+)$", "\\1", t_Jmax_predicted$id)
t_Jmax_predicted$plant_id <- gsub("^(\\d+)_.*", "\\1", t_Jmax_predicted$id)

# combine with t_Jmax_predicted
t_Jmax_predicted <- merge(t_Jmax_predicted, meta, by = "plant_id")

# combine with t_vcmax_para

t_Jmax_para <- merge(t_Jmax_para, meta, by = "plant_id")

##############################
# Kopt model
Kopt_Jmax <- t_Jmax_para %>%
  filter(parameters == "kopt")

Jmax_kopt_model <- glm(Estimate ~ growth_temp * irrigation_trt, 
                       data = Kopt_Jmax, family = "inverse.gaussian") 

AIC(Jmax_kopt_model)


Anova(Jmax_kopt_model, type = 3)
summary(Jmax_kopt_model)
hist(residuals(Jmax_topt_model))

mean_kopt_jmax <- as.data.frame(emmeans(Jmax_kopt_model, ~ growth_temp:irrigation_trt))

ggplot(mean_kopt_jmax, aes(growth_temp, emmean, color = irrigation_trt)) +
  geom_point(size = 5, position = position_dodge(width = 1), alpha = 1) +
  geom_point(data = Kopt_Jmax, aes(x = growth_temp, y = Estimate, 
                                    color = irrigation_trt), size = 2.5,
             alpha = 0.5,
             position = position_jitterdodge()) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = irrigation_trt),
                width = 0.2, linewidth = 1, 
                position = position_dodge(width = 1)) +
  scale_color_manual(values = c("brown", "blue")) +
  labs(x = NULL, y = " Kopt") +
  
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank(),
        legend.position = c(0.15, 0.9))

#############################

# Topt model for Jmax
Topt_Jmax <- t_Jmax_para %>%
  filter(parameters == "Topt")

Jmax_topt_model <- glm(Estimate ~ growth_temp * irrigation_trt, data = Topt_Jmax,
                       family = "inverse.gaussian")

Anova(Jmax_topt_model, type = 3)
summary(Jmax_topt_model)

mean_topt_jmax <- as.data.frame(emmeans(Jmax_topt_model, ~ growth_temp:irrigation_trt, 
                                        type = "response"))

ggplot(mean_topt_jmax, aes(growth_temp, response, color = irrigation_trt)) +
  geom_point(size = 5, position = position_dodge(width = 1), alpha = 1) +
  geom_point(data = Topt_Jmax, aes(x = growth_temp, y = Estimate, 
                                    color = irrigation_trt), size = 2.5,
             alpha = 0.5,
             position = position_jitterdodge()) +
  geom_errorbar(aes(ymin = response - SE, ymax = response + SE, color = irrigation_trt),
                width = 0.2, linewidth = 1, 
                position = position_dodge(width = 1)) +
  scale_color_manual(values = c("brown", "blue")) +
  labs(x = NULL, y = expression(Topt ~ or ~ Jmax ~ (degree*C))) +
  
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.9)) +
  scale_y_continuous(labels = function (x) x - 273)



# Ha model
Ha_Jmax <- t_Jmax_para %>%
  filter(parameters == "Ha")

Jmax_Ha_model <- lm(Estimate ~ growth_temp * irrigation_trt, data = Ha_Jmax)

Anova(Jmax_Ha_model, type = 3)
summary(Jmax_Ha_model)

mean_Ha_jmax <- as.data.frame(emmeans(Jmax_Ha_model, ~ growth_temp:irrigation_trt))

ggplot(mean_Ha_jmax, aes(growth_temp, emmean, color = irrigation_trt)) +
  geom_point(size = 5, position = position_dodge(width = 1), alpha = 1) +
  geom_point(data = Ha_Jmax, aes(x = growth_temp, y = Estimate, 
                                  color = irrigation_trt), size = 2.5,
             alpha = 0.5,
             position = position_jitterdodge()) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = irrigation_trt),
                width = 0.2, linewidth = 1, 
                position = position_dodge(width = 1)) +
  scale_color_manual(values = c("brown", "blue")) +
  labs(x = NULL, y = " Ha") +
  
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.9))


# create a data frame from mean Kopt, Topt, and Ha

mean_parameters_jmax <- as.data.frame(list(Kopt = mean_kopt_jmax$emmean,
                                      Ha = mean_Ha_jmax$emmean,
                                      Topt = mean_topt_jmax$emmean,
                                      treatments = c("LT+D", "HT+D",
                                                     "LT+WW", "HT+WW")))

# Plotting
tleaf_seq <- seq(285, 320, 1)

t_parameter_Jmax_ext_LTD <- mean_parameters_jmax[1,1] * (200000 * exp((mean_parameters_jmax[1,2] * (tleaf_seq - mean_parameters_jmax[1,3]))/(tleaf_seq * R * mean_parameters_jmax[1,3])))/(200000 - (mean_parameters_jmax[1,2] * (1 - exp((200000 *(tleaf_seq - mean_parameters_jmax[1,3]))/(tleaf_seq * R * mean_parameters_jmax[1,3])))))
t_parameter_Jmax_ext_HTD <- mean_parameters_jmax[2,1] * (200000 * exp((mean_parameters_jmax[2,2] * (tleaf_seq - mean_parameters_jmax[2,3]))/(tleaf_seq * R * mean_parameters_jmax[2,3])))/(200000 - (mean_parameters_jmax[2,2] * (1 - exp((200000 *(tleaf_seq - mean_parameters_jmax[2,3]))/(tleaf_seq * R * mean_parameters_jmax[2,3])))))
t_parameter_Jmax_ext_LTWW <- mean_parameters_jmax[3,1] * (200000 * exp((mean_parameters_jmax[3,2] * (tleaf_seq - mean_parameters_jmax[3,3]))/(tleaf_seq * R * mean_parameters_jmax[3,3])))/(200000 - (mean_parameters_jmax[3,2] * (1 - exp((200000 *(tleaf_seq - mean_parameters_jmax[3,3]))/(tleaf_seq * R * mean_parameters_jmax[3,3])))))
t_parameter_Jmax_ext_HTWW <- mean_parameters_jmax[4,1] * (200000 * exp((mean_parameters_jmax[4,2] * (tleaf_seq - mean_parameters_jmax[4,3]))/(tleaf_seq * R * mean_parameters_jmax[4,3])))/(200000 - (mean_parameters_jmax[4,2] * (1 - exp((200000 *(tleaf_seq - mean_parameters_jmax[4,3]))/(tleaf_seq * R * mean_parameters_jmax[4,3])))))

t_parameter_Jmax_ext_df <- data.frame("Tleaf" = tleaf_seq,
                                 "LT_D" = t_parameter_Jmax_ext_LTD,
                                 "HT_D" = t_parameter_Jmax_ext_HTD,
                                 "LT_WW" = t_parameter_Jmax_ext_LTWW,
                                 "HT_WW" = t_parameter_Jmax_ext_HTWW)

t_parameter_Jmax_ext_df <- t_parameter_Jmax_ext_df %>%
  pivot_longer(!Tleaf, names_to = "treatments", values_to = "Jmax")

t_parameter_Jmax_ext_df <- t_parameter_Jmax_ext_df %>%
  mutate(Tleaf = Tleaf - 273.15)


# Jmax model
Jmax_model <- lm(average ~ Tleaf * treatment ,
                  data = t_Jmax_predicted)

Anova(Jmax_model, type = 3)

mean_Jmax <- as.data.frame(emmeans(Jmax_model, ~Tleaf))
mean_Jmax$tleaf <- as.numeric(as.character(mean_Jmax$Tleaf))


(Jmax_graph <- ggplot(mean_Jmax, aes(tleaf, emmean)) +
  geom_point(size = 4) +
  geom_line(data = t_parameter_Jmax_ext_df,
            aes(x = Tleaf, y = Jmax, color = treatments), linewidth = 1) +
  theme_bw(base_size = 16) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.5, linewidth = 1) +
  scale_x_continuous(limits = c(10,50),
                     breaks = seq(0,50,5)) +
  labs(x = expression(Leaf ~ Temperature ~ (degree*C)),
       y = expression(Jmax ~ (mu*mol ~ m^{-2} ~ s^{-1}))) +
 
   scale_color_manual(values = c("turquoise", "orange","magenta", "brown" ),
                     label = c("High Temperature, Drought", "High Temperature, Irrigated",
                               "Ambient Temperature, Drought", "Ambient Temperature, Irrigated")) +
  theme(legend.title = element_blank(),
        legend.position = c(0.15, 0.84),
        panel.grid = element_blank()))


ggarrange(vcmax_graph,
          Jmax_graph,
          common.legend = TRUE)


