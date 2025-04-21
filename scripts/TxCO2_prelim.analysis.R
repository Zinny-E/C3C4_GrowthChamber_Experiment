#############################
#load libraries
##############################
##libraries
library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(multcomp)
library(MuMIn)
library(multcompView)

# Turn off digit rounding in emmean args
emm_options(opt.digits = FALSE)

############Load complied datasheet
df <- read.csv("C:/Users/NickAdmin/OneDrive - Texas Tech University/git/C3C4_GrowthChamber_Experiment/TxCO2_fulldataset.csv", 
               na.strings = "NA")


# add in columns for analysis/plotting for growth temperature
## co2 and temperature treatments
#################################
##Vpmax and amax
#################################
df$co2_trt[df$trt == 'AC+LT' | df$trt == 'AC+HT'] <- 'AC'
df$co2_trt[df$trt == 'EC+LT' | df$trt == 'EC+HT'] <- 'EC'
df$temp_trt[df$trt == 'AC+LT' | df$trt == 'EC+LT'] <- 'LT'
df$temp_trt[df$trt == 'AC+HT' | df$trt == 'EC+HT'] <- 'HT'
df$vpmaxLAGT <- NA
df$vpmaxLAGT[df$temp_trt == 'HT'] <- df$vpmaxLA35[df$temp_trt == 'HT']
df$vpmaxLAGT[df$temp_trt == 'LT'] <- df$vpmaxLA20[df$temp_trt == 'LT']
df$temp_trt[df$trt == 'AC+HT' | df$trt == 'EC+HT'] <- 'HT'
df$amaxLAGT <- NA
df$amaxLAGT[df$temp_trt == 'HT'] <- df$amaxLA35[df$temp_trt == 'HT']
df$amaxLAGT[df$temp_trt == 'LT'] <- df$amaxLA20[df$temp_trt == 'LT']

###########################################
##vcmax and jmax
#############################################
df$vcmaxGT <- NA
df$vcmaxGT[df$temp_trt == 'HT'] <- df$vcmax35[df$temp_trt == 'HT']
df$vcmaxGT[df$temp_trt == 'LT'] <- df$vcmax20[df$temp_trt == 'LT']
df$jmaxGT <- NA
df$jmaxGT[df$temp_trt == 'HT'] <- df$jmax35[df$temp_trt == 'HT']
df$jmaxGT[df$temp_trt == 'LT'] <- df$jmax20[df$temp_trt == 'LT']

###########################################
##anet420
#############################################
df$anet420 <- NA
df$anet1000 <- NA
df$anet420[df$co2_trt == 'AC'] <- df$anet_growth[df$co2_trt == 'AC']
df$anet1000[df$co2_trt == 'EC'] <- df$anet_growth[df$co2_trt == 'EC']


#################################jmax27.5/vcmax27.5
df <- df %>%
  mutate(jmax_vcmax = as.numeric(jmax27.5) / as.numeric(vcmax27.5))



#################################jmax27.5/vcmax27.5
df <- df %>%
  mutate(amax_vpmax = as.numeric(amaxLA27.5) / as.numeric(vpmaxLA27.5))



#df.long <- pivot_longer(df,cols = c("vcmax20", "vcmax27.5", "vcmax35"),
#                       names_to = "growth.temp", values_to = "vcmax")


###########################statistical analysis##################
############################################
#vcmax27.5 #########vcmaxCT
############################################
##in chr changing to numeric 
df$vcmax27.5 <- as.numeric(df$vcmax27.5)

vcmaxCT_lm <- lm(vcmax27.5 ~ co2_trt*temp_trt*species,
                 data = subset(df, ps_pathway == "C3"),
                 na.action = na.exclude)



# Check model assumptions
plot(resid(vcmaxCT_lm) ~ fitted(vcmaxCT_lm))
qqnorm(residuals(vcmaxCT_lm))
qqline(residuals(vcmaxCT_lm))
hist(residuals(vcmaxCT_lm))
shapiro.test(residuals(vcmaxCT_lm))
outlierTest(vcmaxCT_lm)

######################model results################
Anova(vcmaxCT_lm)
summary(vcmaxCT_lm)
r.squaredGLMM(vcmaxCT_lm)

###########Individual effects#########################
emmeans(vcmaxCT_lm, ~species)
emmeans(vcmaxCT_lm, ~co2_trt)
emmeans(vcmaxCT_lm, ~temp_trt*species)
emmeans(vcmaxCT_lm, ~temp_trt)
emmeans(vcmaxCT_lm, ~co2_trt*temp_trt)

#test(emtrends(vpmaxLAGT_lm, ~1, "growth_temp"))


############Post-hoc tests
cld(emmeans(vcmaxCT_lm, pairwise~species))
cld(emmeans(vcmaxCT_lm, pairwise~temp_trt))
cld(emmeans(vcmaxCT_lm, pairwise~co2_trt))
cld(emmeans(vcmaxCT_lm, pairwise~temp_trt*co2_trt*species))



.############################################################
################jmax27.5(jmaxCT)##################
###########################################################
df$jmax27.5 <- as.numeric(df$jmax27.5)

jmaxCT_lm <- lm(log(jmax27.5) ~ co2_trt*temp_trt*species, 
                 data =  subset(df, ps_pathway == "C3"), 
                na.action = na.exclude)


# Check model assumptions
plot(resid(jmaxCT_lm) ~ fitted(jmaxCT_lm))
qqnorm(residuals(jmaxCT_lm))
qqline(residuals(jmaxCT_lm))
hist(residuals(jmaxCT_lm))
shapiro.test(residuals(jmaxCT_lm))
outlierTest(jmaxCT_lm)


######################model results################
Anova(jmaxCT_lm)
summary(jmaxCT_lm)
r.squaredGLMM(jmaxCT_lm)


coef(jmaxCT_lm)
###########Individual effects#########################
emmeans(jmaxCT_lm, ~species)
emmeans(jmaxCT_lm, pairwise~co2_trt, type = "response")
emmeans(jmaxCT_lm, pairwise~temp_trt, type = "response")
emmeans(jmaxCT_lm, pairwise~temp_trt*co2_trt*species, type = "response")
emmeans(jmaxCT_lm, ~temp_trt*species, type = "response")
#test(emtrends(vpmaxLAGT_lm, ~1, "growth_temp"))

############Post-hoc tests
cld(emmeans(jmaxCT_lm, pairwise~co2_trt, type = "response"))
cld(emmeans(jmaxCT_lm, pairwise~temp_trt, type = "response"))
cld(emmeans(jmaxCT_lm, pairwise~temp_trt*species, type = "response"))
cld(emmeans(jmaxCT_lm, pairwise~co2_trt*temp_trt*species, type = "response"))



############################################################
###################vcmaxGT##############
vcmaxGT_lm <- lm(log(vcmaxGT) ~ co2_trt*temp_trt*species, 
                 data =  subset(df, ps_pathway == "C3"), na.action = na.exclude)



####check assumptions
plot(resid(vcmaxGT_lm) ~ fitted(vcmaxGT_lm))
qqnorm(residuals(vcmaxGT_lm))
qqline(residuals(vcmaxGT_lm))
hist(residuals(vcmaxGT_lm))
shapiro.test(residuals(vcmaxGT_lm))
outlierTest(vcmaxGT_lm)



################model results###############
Anova(vcmaxGT_lm)
summary(vcmaxGT_lm)
r.squaredGLMM(vcmaxGT_lm)


###########Individual effects#########################
emmeans(vcmaxGT_lm, ~species)
emmeans(vcmaxGT_lm, ~temp_trt)
emmeans(vcmaxGT_lm, ~temp_trt*species)
emmeans(vcmaxGT_lm, ~co2_trt)
emmeans(vcmaxGT_lm, ~temp_trt*co2_trt)


############Post-hoc tests
cld(emmeans(vcmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"))
cld(emmeans(vcmaxGT_lm, pairwise~temp_trt, type = "response"))
cld(emmeans(vcmaxGT_lm, pairwise~co2_trt, type = "response"))
cld(emmeans(vcmaxGT_lm, pairwise~co2_trt*temp_trt*species, type = "response"))



############################################################
##########jmaxGT########
jmaxGT_lm <- lm(log(jmaxGT) ~ co2_trt*temp_trt*species, 
                data =  subset(df, ps_pathway == "C3"), na.action = na.exclude)


####check assumptions
plot(resid(jmaxGT_lm) ~ fitted(jmaxGT_lm))
qqnorm(residuals(jmaxGT_lm))
qqline(residuals(jmaxGT_lm))
hist(residuals(jmaxGT_lm))
shapiro.test(residuals(jmaxGT_lm))
outlierTest(jmaxGT_lm)
            
            
            
###########################model results
Anova(jmaxGT_lm)
summary(jmaxGT_lm)
r.squaredGLMM(jmaxGT_lm)


###################individual effects##############
emmeans(jmaxGT_lm, ~species)
emmeans(jmaxGT_lm, ~co2_trt)
emmeans(jmaxGT_lm, ~temp_trt)
emmeans(jmaxGT_lm, ~temp_trt*co2_trt)


############Post-hoc tests
cld(emmeans(jmaxGT_lm, pairwise~temp_trt, type = "response"))
cld(emmeans(jmaxGT_lm, pairwise~temp_trt*species, type = "response"))
cld(emmeans(jmaxGT_lm, pairwise~species, type = "response"))





########################################
####vpmax27.5####################

vpmaxLACT_lm <- lm(log(vpmaxLA27.5) ~ co2_trt*temp_trt*species, 
                   data = subset(df, ps_pathway  ==  "C4"),
                   na.action = na.exclude)

# Check model assumptions
plot(resid(vpmaxLACT_lm) ~ fitted(vpmaxLACT_lm))
qqnorm(residuals(vpmaxLACT_lm))
qqline(residuals(vpmaxLACT_lm))
hist(residuals(vpmaxLACT_lm))
shapiro.test(residuals(vpmaxLACT_lm))
outlierTest(vpmaxLACT_lm)

######################model results
Anova(vpmaxLACT_lm)
summary(vpmaxLACT_lm)
r.squaredGLMM(vpmaxLACT_lm)


###########Individual effects#########################
emmeans(vpmaxLACT_lm, ~species)
emmeans(vpmaxLACT_lm, pairwise~temp_trt, type = "response")
emmeans(vpmaxLACT_lm, ~co2_trt*temp_trt, type = "response")
#test(emtrends(vpmaxLAGT_lm, ~1, "growth_temp"))


# Post-hoc tests
cld(emmeans(vpmaxLACT_lm, pairwise~co2_trt*temp_trt*species, type = "response"))
cld(emmeans(vpmaxLACT_lm, pairwise~temp_trt, type = "response"))




#################################################################
##########################amaxLACT################
amaxLACT_lm <- lm(log(amaxLA27.5) ~ co2_trt*temp_trt*species, 
                  subset(df, ps_pathway  ==  "C4"),
                  na.action = na.exclude)



# Check model assumptions
plot(resid(amaxLACT_lm) ~ fitted(amaxLACT_lm))
qqnorm(residuals(amaxLACT_lm))
qqline(residuals(amaxLACT_lm))
hist(residuals(amaxLACT_lm))
shapiro.test(residuals(amaxLACT_lm))
outlierTest(amaxLACT_lm)

######################model results
Anova(amaxLACT_lm)
summary(amaxLACT_lm)
r.squaredGLMM(amaxLACT_lm)


#####################individual effect
emmeans(amaxLACT_lm, ~species)
emmeans(amaxLACT_lm, ~temp_trt*co2_trt)
emmeans(amaxLACT_lm, ~temp_trt)


# Post-hoc tests
cld(emmeans(amaxLACT_lm, pairwise~species, type = "response"))
cld(emmeans(amaxLACT_lm, pairwise~temp_trt, type = "response"))
cld(emmeans(amaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))






#####################################################
##########vpmaxGT#######################
vpmaxLAGT_lm <- lm(log(vpmaxLAGT) ~ co2_trt*temp_trt*species, 
                   data = subset(df, ps_pathway  ==  "C4"),
                   na.action = na.exclude)



# Check model assumptions
plot(resid(vpmaxLAGT_lm) ~ fitted(vpmaxLAGT_lm))
qqnorm(residuals(vpmaxLAGT_lm))
qqline(residuals(vpmaxLAGT_lm))
hist(residuals(vpmaxLAGT_lm))
shapiro.test(residuals(vpmaxLAGT_lm))
outlierTest(vpmaxLAGT_lm)

######################model results
Anova(vpmaxLAGT_lm)
summary(vpmaxLAGT_lm)
r.squaredGLMM(vpmaxLAGT_lm)


###########Individual effects#########################
emmeans(vpmaxLAGT_lm, ~species)
emmeans(vpmaxLAGT_lm, pairwise~temp_trt, type = "response")
emmeans(vpmaxLAGT_lm, ~temp_trt*co2_trt)
emmeans(vpmaxLAGT_lm, ~co2_trt*species)
#test(emtrends(vpmaxLAGT_lm, ~1, "growth_temp"))


# Post-hoc tests
cld(emmeans(vpmaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"))
cld(emmeans(vpmaxLAGT_lm, pairwise~co2_trt*species, type = "response"))
cld(emmeans(vpmaxLAGT_lm, pairwise~co2_trt, type = "response"))





###############################################################
#amaxLAGT 
################################################################

amaxLAGT_lm <- lm(log(amaxLAGT) ~ co2_trt*temp_trt*species, 
                  data = subset(df, ps_pathway  ==  "C4"), 
                  na.action = na.exclude)



# Check model assumptions
plot(resid(amaxLAGT_lm) ~ fitted(amaxLAGT_lm))
qqnorm(residuals(amaxLAGT_lm))
qqline(residuals(amaxLAGT_lm))
hist(residuals(amaxLAGT_lm))
shapiro.test(residuals(amaxLAGT_lm))
outlierTest(amaxLAGT_lm)



############model results
Anova(amaxLAGT_lm)

###########Individual effects#########################
emmeans(amaxLAGT_lm, ~species)
emmeans(amaxLAGT_lm, ~temp_trt)
emmeans(amaxLAGT_lm, pairwise~co2_trt, type = "response")
emmeans(amaxLAGT_lm, pairwise~temp_trt, type = "response")
emmeans(amaxLAGT_lm, ~temp_trt*co2_trt)


# Post-hoc tests
cld(emmeans(amaxLAGT_lm, pairwise~temp_trt, type = "response"))
cld(emmeans(amaxLAGT_lm, pairwise~co2_trt*species, type = "response"))


 #drop the neg

df <- subset(df, !(id == "bou_cur15_t4_ch3" & anet_growth == -0.7564690))


#####################################################
##########Anetgrowth#######################
anetgrowth_lm <- lm(anet_growth ~ co2_trt*temp_trt*species, 
                     data = df,
                     na.action = na.exclude)



# Check model assumptions
plot(resid(anetgrowth_lm) ~ fitted(anetgrowth_lm))
qqnorm(residuals(anetgrowth_lm))
qqline(residuals(anetgrowth_lm))
hist(residuals(anetgrowth_lm))
shapiro.test(residuals(anetgrowth_lm))
outlierTest(anetgrowth_lm)

######################model results
Anova(anetgrowth_lm)
summary(anetgrowth_lm)

###########Individual effects#########################
emmeans(anetgrowth_lm, ~species)
emmeans(anetgrowth_lm, ~temp_trt)
emmeans(anetgrowth_lm, pairwise~co2_trt)
emmeans(anetgrowth_lm, pairwise~temp_trt)
emmeans(anetgrowth_lm, ~temp_trt*co2_trt)


# Post-hoc tests
cld(emmeans(anetgrowth_lm, pairwise~species))
cld(emmeans(anetgrowth_lm, pairwise~temp_trt))
cld(emmeans(anetgrowth_lm, pairwise~co2_trt))
cld(emmeans(anetgrowth_lm, pairwise~temp_trt*co2_trt))
cld(emmeans(anetgrowth_lm, pairwise~temp_trt*species))
cld(emmeans(anetgrowth_lm, pairwise~co2_trt*species))
cld(emmeans(anetgrowth_lm, pairwise~co2_trt*temp_trt*species))

hist(df$marea)

#####################################################
##########gsw#######################
gsw_lm <- lm(log(gsw) ~ co2_trt*temp_trt*species, 
                    data = df,
                    na.action = na.exclude)



# Check model assumptions
plot(resid(gsw_lm) ~ fitted(gsw_lm))
qqnorm(residuals(gsw_lm))
qqline(residuals(gsw_lm))
hist(residuals(gsw_lm))
shapiro.test(residuals(gsw_lm))
outlierTest(gsw_lm)

######################model results
Anova(gsw_lm)
summary(gsw_lm)

###########Individual effects#########################
emmeans(gsw_lm, ~species)
emmeans(gsw_lm, ~temp_trt)
emmeans(gsw_lm, pairwise~co2_trt)
emmeans(gsw_lm, pairwise~temp_trt)
emmeans(gsw_lm, ~temp_trt*co2_trt)


# Post-hoc tests
cld(emmeans(gsw_lm, pairwise~temp_trt, type = "response"))
cld(emmeans(gsw_lm, pairwise~co2_trt, type = "response"))
cld(emmeans(gsw_lm, pairwise~temp_trt*co2_trt, type = "response"))
cld(emmeans(gsw_lm, pairwise~temp_trt*species, type = "response"))
cld(emmeans(gsw_lm, pairwise~co2_trt*temp_trt*species, type = "response"))







############################################################
#####chi











############################################################
#####Nmass
nmass_lm <- lm(log(nmass) ~ co2_trt*temp_trt*species, 
             data = df,
             na.action = na.exclude)



# Check model assumptions
plot(resid(nmass_lm) ~ fitted(nmass_lm))
qqnorm(residuals(nmass_lm))
qqline(residuals(nmass_lm))
hist(residuals(nmass_lm))
shapiro.test(residuals(nmass_lm))
outlierTest(nmass_lm)

######################model results
Anova(nmass_lm)
summary(nmass_lm)

###########Individual effects#########################
emmeans(nmass_lm, ~species)
emmeans(nmass_lm, ~temp_trt)
emmeans(nmass_lm, pairwise~co2_trt)
emmeans(nmass_lm, pairwise~temp_trt)
emmeans(nmass_lm, ~temp_trt*co2_trt)


# Post-hoc tests
cld(emmeans(nmass_lm, pairwise~species, type = "response"))
cld(emmeans(nmass_lm, pairwise~co2_trt, type = "response"))
cld(emmeans(nmass_lm, pairwise~temp_trt, type = "response"))
cld(emmeans(nmass_lm, pairwise~temp_trt*co2_trt, type = "response"))
cld(emmeans(nmass_lm, pairwise~temp_trt*species, type = "response"))
cld(emmeans(nmass_lm, pairwise~co2_trt*species, type = "response"))



############################################################
#####biomass
biomass_lm <- lm(log(biomass.g.) ~ co2_trt*temp_trt*species, 
               data = df,
               na.action = na.exclude)



# Check model assumptions
plot(resid(biomass_lm) ~ fitted(biomass_lm))
qqnorm(residuals(biomass_lm))
qqline(residuals(biomass_lm))
hist(residuals(biomass_lm))
shapiro.test(residuals(biomass_lm))
outlierTest(biomass_lm)

######################model results
Anova(biomass_lm)
summary(nmass_lm)

###########Individual effects#########################
emmeans(nmass_lm, ~species)
emmeans(nmass_lm, ~temp_trt)
emmeans(nmass_lm, pairwise~co2_trt)
emmeans(nmass_lm, pairwise~temp_trt)
emmeans(nmass_lm, ~temp_trt*co2_trt)


# Post-hoc tests
cld(emmeans(biomass_lm, pairwise~temp_trt, type = "response"))
cld(emmeans(biomass_lm, pairwise~co2_trt, type = "response"))
cld(emmeans(biomass_lm, pairwise~temp_trt*co2_trt, type = "response"))
cld(emmeans(biomass_lm, pairwise~temp_trt*species, type = "response"))
cld(emmeans(biomass_lm, pairwise~co2_trt*species, type = "response"))





############################################################
#####Narea
narea_lm <- lm(log(narea) ~ co2_trt*temp_trt*species, 
               data = df,
               na.action = na.exclude)



# Check model assumptions
plot(resid(narea_lm) ~ fitted(narea_lm))
qqnorm(residuals(narea_lm))
qqline(residuals(narea_lm))
hist(residuals(narea_lm))
shapiro.test(residuals(narea_lm))
outlierTest(narea_lm)

######################model results
Anova(narea_lm)
summary(narea_lm)

###########Individual effects#########################
emmeans(narea_lm, ~species)
emmeans(narea_lm, ~temp_trt)
emmeans(narea_lm, pairwise~co2_trt)
emmeans(narea_lm, pairwise~temp_trt)
emmeans(narea_lm, ~temp_trt*co2_trt)


# Post-hoc tests
cld(emmeans(narea_lm, pairwise~temp_trt, type = "response"))
cld(emmeans(narea_lm, pairwise~co2_trt, type = "response"))
cld(emmeans(narea_lm, pairwise~co2_trt*temp_trt, type = "response"))
cld(emmeans(narea_lm, pairwise~temp_trt*species, type = "response"))
cld(emmeans(narea_lm, pairwise~co2_trt*species, type = "response"))







############################################################
#####marea
marea_lm <- lm(log(marea) ~ co2_trt*temp_trt*species, 
               data = df,
               na.action = na.exclude)



# Check model assumptions
plot(resid(marea_lm) ~ fitted(marea_lm))
qqnorm(residuals(marea_lm))
qqline(residuals(marea_lm))
hist(residuals(marea_lm))
shapiro.test(residuals(marea_lm))
outlierTest(marea_lm)

######################model results
Anova(marea_lm)
summary(marea_lm)

###########Individual effects#########################
emmeans(marea_lm, ~species)
emmeans(marea_lm, ~temp_trt)
emmeans(marea_lm, pairwise~co2_trt)
emmeans(marea_lm, pairwise~temp_trt)
emmeans(marea_lm, ~temp_trt*co2_trt)


# Post-hoc tests
cld(emmeans(marea_lm, pairwise~temp_trt, type = "response"))
cld(emmeans(marea_lm, pairwise~co2_trt, type = "response"))
cld(emmeans(marea_lm, pairwise~temp_trt*co2_trt, type = "response"))
cld(emmeans(marea_lm, pairwise~temp_trt*species, type = "response"))
cld(emmeans(marea_lm, pairwise~co2_trt*species, type = "response"))




# Replace max value(s) column with NA
max_value <- max(df$jmax_vcmax, na.rm = TRUE)  
df$jmax_vcmax[df$jmax_vcmax == max_value] <- NA   
############################################################
#####jmaxvcmax_ratio
jmaxvcmax_lm <- lm(log(jmax_vcmax) ~ co2_trt*temp_trt*species, 
                   data = subset(df, ps_pathway  ==  "C3"),
               na.action = na.exclude)



# Check model assumptions
plot(resid(jmaxvcmax_lm) ~ fitted(jmaxvcmax_lm))
qqnorm(residuals(jmaxvcmax_lm))
qqline(residuals(jmaxvcmax_lm))
hist(residuals(jmaxvcmax_lm))
shapiro.test(residuals(jmaxvcmax_lm))
outlierTest(jmaxvcmax_lm)

######################model results
Anova(jmaxvcmax_lm)
summary(jmaxvcmax_lm)

###########Individual effects#########################
emmeans(jmaxvcmax_lm, ~species)
emmeans(jmaxvcmax_lm, ~temp_trt)
emmeans(jmaxvcmax_lm, pairwise~co2_trt)
emmeans(jmaxvcmax_lm, pairwise~temp_trt)
emmeans(jmaxvcmax_lm, ~temp_trt*co2_trt)


# Post-hoc tests
cld(emmeans(jmaxvcmax_lm, pairwise~temp_trt, type = "response"))
cld(emmeans(jmaxvcmax_lm, pairwise~co2_trt, type = "response"))
cld(emmeans(jmaxvcmax_lm, pairwise~temp_trt*co2_trt, type = "response"))
cld(emmeans(jmaxvcmax_lm, pairwise~temp_trt*species, type = "response"))
cld(emmeans(jmaxvcmax_lm, pairwise~co2_trt*temp_trt*species, type = "response"))






# Replace max value(s) column with NA
max_value <- max(df$amax_vpmax, na.rm = TRUE)  
df$amax_vpmax[df$amax_vpmax == max_value] <- NA
###################################################################
###amaxvpmax ratio
amaxvpmax_lm <- lm(log(amax_vpmax) ~ co2_trt*temp_trt*species, 
               data = df,
               na.action = na.exclude)




plot(resid(amaxvpmax_lm) ~ fitted(amaxvpmax_lm))
qqnorm(residuals(jmaxvcmax_lm))
qqline(residuals(jmaxvcmax_lm))
hist(residuals(jmaxvcmax_lm))
shapiro.test(residuals(jmaxvcmax_lm))
outlierTest(jmaxvcmax_lm)

######################model results
Anova(amaxvpmax_lm)
summary(jmaxvcmax_lm)

###########Individual effects#########################
emmeans(jmaxvcmax_lm, ~species)
emmeans(jmaxvcmax_lm, ~temp_trt)
emmeans(jmaxvcmax_lm, pairwise~co2_trt)
emmeans(jmaxvcmax_lm, pairwise~temp_trt)
emmeans(jmaxvcmax_lm, ~temp_trt*co2_trt)


# Post-hoc tests
cld(emmeans(amaxvpmax_lm, pairwise~temp_trt, type = "response"))
cld(emmeans(amaxvpmax_lm, pairwise~co2_trt, type = "response"))
cld(emmeans(jmaxvcmax_lm, pairwise~temp_trt*co2_trt, type = "response"))
cld(emmeans(amaxvpmax_lm, pairwise~temp_trt*species, type = "response"))
cld(emmeans(jmaxvcmax_lm, pairwise~co2_trt*temp_trt*species, type = "response"))








############################################################
#####delta
delta_lm <- lm(delta ~ co2_trt*temp_trt*species, 
                   data = df,
                   na.action = na.exclude)

# Check model assumptions
plot(resid(delta_lm) ~ fitted(delta_lm))
qqnorm(residuals(delta_lm))
qqline(residuals(delta_lm))
hist(residuals(delta_lm))
shapiro.test(residuals(delta_lm))
outlierTest(delta_lm)

######################model results
Anova(delta_lm)
summary(delta_lm)

###########Individual effects#########################
emmeans(delta_lm, ~species)
emmeans(delta_lm, ~temp_trt)
emmeans(delta_lm, pairwise~co2_trt)
emmeans(delta_lm, pairwise~temp_trt)
emmeans(delta_lm, ~temp_trt*co2_trt)


# Post-hoc tests
cld(emmeans(delta_lm, pairwise~temp_trt, type = "response"))
cld(emmeans(delta_lm, pairwise~co2_trt, type = "response"))
cld(emmeans(delta_lm, pairwise~temp_trt*co2_trt, type = "response"))
cld(emmeans(delta_lm, pairwise~temp_trt*species, type = "response"))
cld(emmeans(delta_lm, pairwise~co2_trt*temp_trt*species, type = "response"))













ggplot(df, aes(x = marea * nmass, y = narea)) +
  geom_point(alpha = 0.6, color = "steelblue", size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +
  labs(x = expression(M[a] ~ "×" ~ N[m]),
       y = expression(N[a]),
       title = "1:1 Relationship") +
  theme_minimal(base_size = 14)




































































































































































##############models with species as random effect######  
#################
############################################
##vcmaxCT(vcmax27.5 )##############
############################################

vcmaxCT_lmer <-  lmer(vcmax27.5 ~ co2_trt*temp_trt + (1|species),
                        data = subset(df, ps_pathway  ==  "C3"), 
                        na.action = na.exclude)



# Check model assumptions
plot(resid(vcmaxCT_lmer) ~ fitted(vcmaxCT_lmer))
qqnorm(residuals(vcmaxCT_lmer))
qqline(residuals(vcmaxCT_lmer))
hist(residuals(vcmaxCT_lmer))
shapiro.test(residuals(vcmaxCT_lmer))
outlierTest(vcmaxCT_lmer)
densityPlot(residuals(vcmaxCT_lmer))


###model results################
Anova(vcmaxCT_lmer)
summary(vcmaxCT_lmer)
r.squaredGLMM(vcmaxCT_lmer)

###########Individual effects#########################
emmeans(vcmaxCT_lmer, pairwise~co2_trt, type = "response")
emmeans(vcmaxCT_lmer, pairwise~temp_trt, type = "response")


############Post-hoc tests
cld(emmeans(vcmaxCT_lmer, pairwise~temp_trt, type = "response"))
cld(emmeans(vcmaxCT_lmer, pairwise~co2_trt, type = "response"))





############################################################
################jmax27.5(jmaxCT)##################
###########################################################
df$jmax27.5 <- as.numeric(df$jmax27.5)

jmaxCT_lmer <- lmer(log(jmax27.5) ~ co2_trt*temp_trt + (1|species),
                    data = subset(df, ps_pathway  ==  "C3"), 
                    na.action = na.exclude)


# Check model assumptions
plot(resid(jmaxCT_lmer) ~ fitted(jmaxCT_lmer))
qqnorm(residuals(jmaxCT_lmer))
qqline(residuals(jmaxCT_lmer))
hist(residuals(jmaxCT_lmer))
shapiro.test(residuals(jmaxCT_lmer))
outlierTest(jmaxCT_lmer)
densityPlot(residuals(jmaxCT_lmer))

######################model results################
Anova(jmaxCT_lmer)
summary(jmaxCT_lmer)
r.squaredGLMM(jmaxCT_lmer)


###########Individual effects#########################
emmeans(jmaxCT_lmer, pairwise~co2_trt, type = "response")
emmeans(jmaxCT_lmer, pairwise~temp_trt, type = "response")


############Post-hoc tests
cld(emmeans(jmaxCT_lmer, pairwise~co2_trt, type = "response"))
cld(emmeans(jmaxCT_lmer, pairwise~temp_trt,type = "response"))



############################################################
#############vcmaxGT#######
#######################################################
vcmaxGT_lmer <- lmer(log(vcmaxGT) ~ co2_trt*temp_trt + (1|species),
                     data = subset(df, ps_pathway  ==  "C3"), 
                     na.action = na.exclude)



####check assumptions
plot(resid(vcmaxGT_lmer) ~ fitted(vcmaxGT_lmer))
qqnorm(residuals(vcmaxGT_lmer))
qqline(residuals(vcmaxGT_lmer))
hist(residuals(vcmaxGT_lmer))
shapiro.test(residuals(vcmaxGT_lmer))
outlierTest(vcmaxGT_lmer)
densityPlot(residuals(vcmaxGT_lmer))


################model results###############
Anova(vcmaxGT_lmer)
summary(vcmaxGT_lmer)
r.squaredGLMM(vcmaxGT_lmer)


###########Individual effects#########################
emmeans(vcmaxGT_lmer, pairwise~temp_trt)


############Post-hoc tests
cld(emmeans(vcmaxGT_lmer, pairwise~temp_trt, type = "response"))




############################################################
##########jmaxGT########
jmaxGT_lmer <- lmer(log(jmaxGT) ~ co2_trt*temp_trt + (1|species),
                  data = subset(df, ps_pathway  ==  "C3"), 
                  na.action = na.exclude)


####check assumptions
plot(resid(jmaxGT_lmer) ~ fitted(jmaxGT_lmer))
qqnorm(residuals(jmaxGT_lmer))
qqline(residuals(jmaxGT_lmer))
hist(residuals(jmaxGT_lmer))
shapiro.test(residuals(jmaxGT_lmer))
outlierTest(jmaxGT_lmer)



###########################model results
Anova(jmaxGT_lmer)
summary(jmaxGT_lmer)
r.squaredGLMM(jmaxGT_lmer)


###################individual effects##############
emmeans(jmaxGT_lmer, ~temp_trt)


############Post-hoc tests
cld(emmeans(jmaxGT_lmer, pairwise~temp_trt, type = "response"))



####C4
########################################
####vpmax27.5####################
vpmaxLACT_lmer <-lmer(vpmaxLA27.5 ~ co2_trt*temp_trt + (1|species),
                    data = subset(df, ps_pathway  ==  "C4"), 
                    na.action = na.exclude)

# Check model assumptions
plot(resid(vpmaxLACT_lmer) ~ fitted(vpmaxLACT_lmer))
qqnorm(residuals(vpmaxLACT_lmer))
qqline(residuals(vpmaxLACT_lmer))
hist(residuals(vpmaxLACT_lmer))
shapiro.test(residuals(vpmaxLACT_lmer))
outlierTest(vpmaxLACT_lmer)

######################model results
Anova(vpmaxLACT_lmer)
summary(vpmaxLACT_lmer)
r.squaredGLMM(vpmaxLACT_lmer)


###########Individual effects#########################
emmeans(vpmaxLACT_lmer, pairwise~temp_trt)
emmeans(vpmaxLACT_lmer, pairwise~co2_trt*temp_trt)



# Post-hoc tests
cld(emmeans(vpmaxLACT_lmer, pairwise~co2_trt*temp_trt))
cld(emmeans(vpmaxLACT_lmer, pairwise~temp_trt))




#################################################################
##########################amaxLACT################
amaxLACT_lmer <- lmer(log(amaxLA27.5) ~ co2_trt*temp_trt + (1|species),
                     data = subset(df, ps_pathway  ==  "C4"), 
                     na.action = na.exclude)



# Check model assumptions
plot(resid(amaxLACT_lmer) ~ fitted(amaxLACT_lmer))
qqnorm(residuals(amaxLACT_lmer))
qqline(residuals(amaxLACT_lmer))
hist(residuals(amaxLACT_lmer))
shapiro.test(residuals(amaxLACT_lmer))
outlierTest(amaxLACT_lmer)

######################model results
Anova(amaxLACT_lmer)
summary(amaxLACT_lmer)
r.squaredGLMM(amaxLACT_lmer)


#####################individual effect
emmeans(amaxLACT_lmer, pairwise~temp_trt)


# Post-hoc tests
cld(emmeans(amaxLACT_lmer, pairwise~temp_trt, type = "response"))



#####################################################
##########vpmaxGT###############
################################################
vpmaxLAGT_lmer <-lmer(vpmaxLAGT ~ co2_trt*temp_trt + (1|species),
                    data = subset(df, ps_pathway  ==  "C4"), 
                    na.action = na.exclude)



# Check model assumptions
plot(resid(vpmaxLAGT_lmer) ~ fitted(vpmaxLAGT_lmer))
qqnorm(residuals(vpmaxLAGT_lmer))
qqline(residuals(vpmaxLAGT_lmer))
hist(residuals(vpmaxLAGT_lmer))
shapiro.test(residuals(vpmaxLAGT_lmer))
outlierTest(vpmaxLAGT_lmer)

######################model results
Anova(vpmaxLAGT_lmer)
summary(vpmaxLAGT_lmer)
r.squaredGLMM(vpmaxLAGT_lmer)


###########Individual effects#########################

emmeans(vpmaxLAGT_lm, pairwise~temp_trt)
emmeans(vpmaxLAGT_lmer, ~temp_trt*co2_trt)



# Post-hoc tests
cld(emmeans(vpmaxLAGT_lmer, pairwise~co2_trt*temp_trt))
cld(emmeans(vpmaxLAGT_lmer, pairwise~co2_trt))






###############################################################
#amaxLAGT 
################################################################

amaxLAGT_lmer <- lmer(log(amaxLAGT) ~ co2_trt*temp_trt + (1|species),
                    data = subset(df, ps_pathway  ==  "C4"), 
                    na.action = na.exclude)



# Check model assumptions
plot(resid(amaxLAGT_lmer) ~ fitted(amaxLAGT_lmer))
qqnorm(residuals(amaxLAGT_lmer))
qqline(residuals(amaxLAGT_lmer))
hist(residuals(amaxLAGT_lmer))
shapiro.test(residuals(amaxLAGT_lmer))
outlierTest(amaxLAGT_lmer)



############model results
Anova(amaxLAGT_lmer)
r.squaredGLMM(amaxLAGT_lmer)
###########Individual effects#########################
emmeans(amaxLAGT_lmer, ~temp_trt)
emmeans(amaxLAGT_lmer, pairwise~co2_trt, type = "response")
emmeans(amaxLAGT_lmer, pairwise~temp_trt, type = "response")



# Post-hoc tests
cld(emmeans(amaxLAGT_lmer, pairwise~co2_trt*temp_trt, type = "response"))
cld(emmeans(amaxLAGT_lmer, pairwise~temp_trt, type = "response"))


#####################Anetgrowth
anetgrowth_lmer <- lmer(anet_growth ~ co2_trt*temp_trt + (1|species),
                      data = subset(df, ps_pathway  ==  "C4"), 
                      na.action = na.exclude)



# Check model assumptions
plot(resid(anetgrowth_lmer) ~ fitted(anetgrowth_lmer))
qqnorm(residuals(anetgrowth_lmer))
qqline(residuals(anetgrowth_lmer))
hist(residuals(anetgrowth_lmer))
shapiro.test(residuals(anetgrowth_lmer))
outlierTest(anetgrowth_lmer)



############model results
Anova(anetgrowth_lmer)
r.squaredGLMM(amaxLAGT_lmer)










