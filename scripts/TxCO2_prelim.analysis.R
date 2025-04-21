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
df <- read.csv("~/git/C3C4_GrowthChamber_Experiment/TXCO2.clean_complied_data.csv", 
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
cld(emmeans(vcmaxCT_lm, pairwise~temp_trt, type = "response"))
cld(emmeans(vcmaxCT_lm, pairwise~co2_trt, type = "response"))
cld(emmeans(vcmaxCT_lm, pairwise~co2_trt*temp_trt*species, type = "response"))
cld(emmeans(vcmaxCT_lm, pairwise~temp_trt*species, type = "response"))




############################################################
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
emmeans(jmaxCT_lm, pairwise~temp_trt*co2_trt, type = "response")
emmeans(jmaxCT_lm, ~temp_trt*co2_trt)
#test(emtrends(vpmaxLAGT_lm, ~1, "growth_temp"))

############Post-hoc tests
cld(emmeans(jmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))
cld(emmeans(jmaxCT_lm, pairwise~temp_trt*species, type = "response"))
cld(emmeans(jmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))
cld(emmeans(jmaxCT_lm, pairwise~temp_trt*species, type = "response"))


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
cld(emmeans(vcmaxGT_lm, pairwise~temp_trt*species, type = "response"))



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
cld(emmeans(jmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"))
cld(emmeans(jmaxGT_lm, pairwise~temp_trt, type = "response"))
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
cld(emmeans(vpmaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))
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
cld(emmeans(amaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))
cld(emmeans(amaxLACT_lm, pairwise~temp_trt, type = "response"))







#####################################################
##########vpmaxGT#######################
vpmaxLAGT_lm <- lm(vpmaxLAGT ~ co2_trt*temp_trt*species, 
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
cld(emmeans(amaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"))
cld(emmeans(amaxLAGT_lm, pairwise~co2_trt*species, type = "response"))

contrast(emmeans(amaxLAGT_lm, ~co2_trt*temp_trt, type = "response"), 
         simple =  "each", combine = TRUE )





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










