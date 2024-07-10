# example script for fitting ACi curves for C3 and C4 plants as implemented in the SLCE dataset

############################################
############################################
############################################
## C3 example
############################################
############################################
############################################

############################################
# load plantecophys package to fit C3 curves
############################################
library(plantecophys)

############################################
# set parameters
############################################
id_c3 <- 'Auburn_Iopa_1' # identifier for C3 curve of interest
b_tresp_R <- 0.1781075 # value for b parameter in temperature response curve for respiration (From Smith and Dukes (2017), doi = 10.1111/gcb.13735)
c_tresp_R <- -0.00179152 # value for c parameter in temperature response curve for respiration (From Smith and Dukes (2017), doi = 10.1111/gcb.13735)

############################################
# read in ACi curve
############################################
aci_c3 <- read.csv('LCE_ACi_curves/Auburn_Iopa_1.csv') # set path within local environment

############################################
# read in SLCE_data
############################################
slce <- read.csv('LCE_data.csv') # set path to local environment
slce_c3 = subset(slce, aci_id == id_c3) # slce data for individual of interest

############################################
# extract dark respiration value from SLCE_data and add to Aci dataset
############################################
rd_c3 <- slce_c3$Rd # find Rd value for curve of interest
rd_c3_tphoto <- rd_c3 / exp((b_tresp_R * (slce_c3$Tleaf_R - slce_c3$Tleaf_photo)) + (c_tresp_R * (slce_c3$Tleaf_R^2 - slce_c3$Tleaf_photo^2))) # adjust Rd to similar temperature as the A/Ci data
aci_c3$Rd <- rd_c3_tphoto # add dark respiration to ACi data

############################################
# fit Aci curve
############################################
fit_c3 <- fitaci(aci_c3, varnames = list(ALEAF = "Photo", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Pari", Rd = "Rd"), 
                 Tcorrect=FALSE, useRd=TRUE) # fit curve

############################################
############################################
############################################
## C4 example
############################################
############################################
############################################

############################################
# load minpack.lm package to fit C4 curves using nlsLM
############################################
library(minpack.lm)
library(dplyr)
############################################
# load implementation of C4 ACi curve fitter (based on von Caemmerer 2000 and function AciC4 from plantecophys)
############################################
#general inputs
O2=210
FRM=0.5
alpha=0
Q10=2
x=0.4
THETA=0.7
low_gammastar <- 1.93e-4 # Half the reciprocal for Rubisco specificity (NOT CO2 compensation point)
Vpr=80
gbs= 3e-3

# enzyme-limited photosynthesis function
A.enzyme.func=Photo ~ (-(-(((pmin(Ci * Vpmax/(Ci + Kp), Vpr)) - Rm + gbs * Ci) + (Vcmax - Rd) + gbs * 
                             K +((alpha/0.047) * (low_gammastar * Vcmax + Rd * Kc/Ko)))) - 
                         sqrt((-(((pmin(Ci * Vpmax/(Ci + Kp), Vpr)) - Rm + gbs * Ci) + (Vcmax - Rd) + gbs * 
                                   K +((alpha/0.047) * (low_gammastar * Vcmax + Rd * Kc/Ko))))^2 - 4 * a.c * 
                                ((Vcmax - Rd) * ((pmin(Ci * Vpmax/(Ci + Kp), Vpr)) - Rm + gbs * Ci) - 
                                   (Vcmax * gbs * low_gammastar * O2 + Rd * gbs * K))))/(2 * a.c)-Rd

# light-limited photosynthesis function
A.light.func=Photo ~ (-(-((x * ((1/(2 * THETA)) * 
                                  (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 4 * THETA * Qp2 * Jmax)))/2 - 
                             Rm + gbs * Ci) + ((1 - x) * ((1/(2 * THETA)) * 
                                                            (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 4 * THETA * Qp2 * Jmax)))/3 - Rd) +
                            gbs * (7 * low_gammastar * O2/3) + alpha * low_gammastar/0.047 *((1 - x) * 
                                                                                               ((1/(2 * THETA)) * 
                                                                                                  (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 4 * THETA * Qp2 * Jmax)))/3 + 7*Rd/3))) - 
                        sqrt((-((x * ((1/(2 * THETA)) * (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 4 * THETA * Qp2 * Jmax)))/2 - 
                                   Rm + gbs * Ci) + ((1 - x) * ((1/(2 * THETA)) * 
                                                                  (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 4 * THETA * Qp2 * Jmax)))/3 - Rd) +
                                  gbs * (7 * low_gammastar * O2/3) + alpha * low_gammastar/0.047 *((1 - x) * 
                                                                                                     ((1/(2 * THETA)) * 
                                                                                                        (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 4 * THETA * Qp2 * Jmax)))/3 + 7*Rd/3)))^2 - 
                               4 * a.j * (((x * ((1/(2 * THETA)) * (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 4 * THETA * Qp2 * Jmax)))/2 - 
                                              Rm + gbs * Ci) * ((1 - x) * ((1/(2 * THETA)) * (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 4 * 
                                                                                                                  THETA * Qp2 * Jmax)))/3 - Rd)) -gbs * 
                                            low_gammastar * O2 * ((1 - x) * ((1/(2 * THETA)) * 
                                                                               (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 4 * THETA * Qp2 * Jmax)))/3 + 7 * Rd/3))))/(2 * a.j)-Rd

############################################
# set parameters
############################################
id_c4 <- 'Auburn_Iopa_1' # identifier for C4 curve of interest
b_tresp_R <- 0.1781075 # value for b parameter in temperature response curve for respiration (From Smith and Dukes (2017), doi = 10.1111/gcb.13735)
c_tresp_R <- -0.00179152 # value for c parameter in temperature response curve for respiration (From Smith and Dukes (2017), doi = 10.1111/gcb.13735)

############################################
# read in ACi curve
############################################
#aci_c4 <- read.csv('LCE_ACi_curves/UMBS_Arub_3.csv') # set path within local environment
aci_c4 <- read.csv('Chamber.5/bou_cur31_t4_ch5.csv')
aci_c4 <- aci_c4 %>% 
              rename(Photo = A)
  ############################################
# read in SLCE_data
############################################
#slce <- read.csv('LCE_data.csv') # set path to local environment
#slce_c4 = subset(slce, aci_id == id_c4) # slce data for individual of interest

############################################
# extract dark respiration value from SLCE_data and add to Aci dataset
############################################
#rd_c4 <- slce_c4$Rd # find Rd value for curve of interest
#rd_c4_tphoto <- rd_c4 / exp((b_tresp_R * (slce_c4$Tleaf_R - slce_c4$Tleaf_photo)) + 
  #                            (c_tresp_R * (slce_c4$Tleaf_R^2 - slce_c4$Tleaf_photo^2))) # adjust Rd to similar temperature as the A/Ci data
#aci_c4$Rd <- rd_c4_tphoto # add dark respiration to ACi data

############################################
# set fit parameters (see documentation for function AciC4 from plantecophys)
############################################
#Tleaf= mean(aci_c4$Tleaf)
#Rd= rd_c4_tphoto
#PPFD=mean(aci_c4$Pari)

Tleaf= aci_c4$Tleaf[1]
Rd =aci_c4$Rd[1]
PPFD=mean(aci_c4$Qin)
# Michaelis-Menten coefficients for CO2 (Kc, mu mol mol-1) and O (Ko, mmol mol-1) and combined (K)
Kc <- 650*Q10^((Tleaf-25)/10)
Kp <- 80*Q10^((Tleaf-25)/10)
Ko <- 450*Q10^((Tleaf-25)/10)
K <- Kc*(1+O2/Ko)
Rm <-  FRM*Rd # Day leaf respiration, umol m-2 s-1
Qp2 <- PPFD*0.85*(1-0.15)/2 # Non-rectangular hyperbola describing light effect on electron transport rate (J)
a.c <- 1 - (alpha*Kc)/(0.047*Ko)
a.j <- 1 - 7 * low_gammastar * alpha/(3 * 0.047)

############################################
# plot data to estimate transition point & find irregularities
############################################
plot(aci_c4$Photo~aci_c4$Ci)
ci_trans <- 750   # estimated Ci transition point
# note the need to remove Ci point below 0

############################################
# fit the enzyme and light limited portions of the curve separately
############################################
fit_enzyme <- nlsLM(A.enzyme.func, data = subset(aci_c4, Ci < ci_trans & Ci > 0), 
                    start=list(Vcmax=30,Vpmax=100), control=nls.control(maxiter=500, minFactor=1/10000)) # fit the enzyme limited portion of the curve (Vcmax and Vpmax)
fit_light = nlsLM(A.light.func, data= subset(aci_c4, Ci >= ci_trans), 
                  start = list(Jmax = 130), control = nls.control(maxiter = 500, minFactor = 1/10000)) # fit the light-limted portion of the curve (Jmax)


coef(fit_enzyme)

coef(fit_light)

Vcmax <- coef(fit_enzyme)["Vcmax"]
Vpmax <- coef(fit_enzyme)["Vpmax"]
Jmax <- coef(fit_light)["Jmax"]


id_name <- 'sor_nut25_t1_ch5_27.5'



# Creating the dataframe
chamber5.coeff<- data.frame("ID" = id_name,
                             pathway = "C4",
                              "Vcmax" = Vcmax,
                              "Vpmax" = Vpmax,
                              "Jmax" = Jmax,
                            "ci_trans" = ci_trans)


# Creating the dataframe
chamber5.coeff[17,] <- data.frame("ID" = id_name,
                             pathway = "C4",
                             "Vcmax" = Vcmax,
                             "Vpmax" = Vpmax,
                             "Jmax" = Jmax,
                             "ci_trans" = ci_trans)



write.csv(chamber2.coeff, "~/git_repo/LCE/aci_coeff/aci.coeff_chamber2.csv", row.names = FALSE)







data.c4 <- read.csv("~/git_repo/LCE/aci_coeff/aci.coeff_chamber2.csv")



# Plotting the entire curve
plot(aci_c4$Ci, aci_c4$Photo, 
     xlab = "Ci", 
     ylab = "A", 
     pch = 19, 
     cex = 0.5)

# Add fitted curves
points(aci_c4$Ci[aci_c4$Ci < ci_trans & aci_c4$Ci > 0], 
       predict(fit_enzyme, newdata = aci_c4[aci_c4$Ci < ci_trans & aci_c4$Ci > 0,]), 
       col = "blue", lwd = 2)

points(aci_c4$Ci[aci_c4$Ci < ci_trans & aci_c4$Ci > 0], 
       predict(fit_enzyme), 
       col = "blue", lwd = 2)

points(aci_c4$Ci[aci_c4$Ci >= ci_trans], 
       predict(fit_light, newdata = aci_c4[aci_c4$Ci >= ci_trans,]), 
       col = "red", lwd = 2)
# Add a vertical line for ci_trans
abline(v = ci_trans, col = "black", lwd = 2, lty = 2)

#output
Vcmax <- coef(fit_enzyme)["Vcmax"]
Vpmax <- coef(fit_enzyme)["Vpmax"]
Jmax <- coef(fit_light)["Jmax"]






library(ggplot2)
library(dplyr)
ggplot(data.c4, aes(x = temp, y = Vcmax, color = as.factor(ID))) +
  geom_boxplot() + # Use points to represent data
  labs(x = "Temperature", y = "Vcmax", color = "ID") +
  theme_minimal() # Optional: Use a minimal theme for the plot
  

