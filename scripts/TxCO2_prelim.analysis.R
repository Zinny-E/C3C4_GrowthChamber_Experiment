##plots
#############################
#load libraries
##############################
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(gridExtra)
library(lme4)
library(car)
library(emmeans)
library(multcomp)
library(MuMIn)

# Load compiled datasheet
df <- read.csv("~/git/C3C4_GrowthChamber_Experiment/TXCO2.complied_data.csv", 
               na.strings = "NA")

df <- read.csv("~/Desktop/testdata.csv", 
               na.strings = "NA")

# add in columns for analysis/plotting for growth temperature
## co2 and temperature treatments
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
vcmaxCT_lm <- lm(vcmax27.5 ~ co2_trt*temp_trt*species,
                 data = subset(df, ps_pathway == "C3"),
                 na.action = na.exclude)


df$vcmax27.5 <- as.numeric(df$vcmax27.5)

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
emmeans(vcmaxCT_lm, ~co2_trt, type = "response")
emmeans(vcmaxCT_lm, ~temp_trt*species, type = "response")
emmeans(vcmaxCT_lm, ~temp_trt, type = "response")
emmeans(vcmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response")
emmeans(vcmaxCT_lm, ~temp_trt*co2_trt*species)
#test(emtrends(vpmaxLAGT_lm, ~1, "growth_temp"))


############Post-hoc tests
cld(emmeans(vcmaxCT_lm, pairwise~temp_trt, type = "response"))
cld(emmeans(vcmaxCT_lm, pairwise~co2_trt, type = "response"))
cld(emmeans(vcmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))
cld(emmeans(vcmaxCT_lm, pairwise~temp_trt*species, type = "response"))


##creating letter for grouping used for making plots
cld(emmeans(vcmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))[1, 8]
cld(emmeans(vcmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))[2, 8]
cld(emmeans(vcmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))[3, 8]
cld(emmeans(vcmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))[4, 8]


vcmaxCT_letter <- data.frame(x= c(0.8, 1.2, 1.8, 2.2),
                             temp_trt = c('LT', 'LT', 'HT',  'HT'),
                             co2_trt = c('AC', 'EC', 'AC', 'EC'),
                             y =c(100, 125, 175, 120),
                             group = c(cld(emmeans(vcmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))[1, 8],
                                       cld(emmeans(vcmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))[2, 8],
                                       cld(emmeans(vcmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))[3, 8],
                                       cld(emmeans(vcmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))[4, 8]))



vcmaxCT_letter$letter[vcmaxCT_letter$group == " 1 "] <- "a"
vcmaxCT_letter$letter[vcmaxCT_letter$group == " 12"] <- "ab"
vcmaxCT_letter$letter[vcmaxCT_letter$group == "  2"] <- "b"


#vcmaxCT<- as.data.frame(emmeans(vcmaxCT_lm, ~co2_trt * temp_trt))

############################################################
################jmax27.5(jmaxCT)##################
###########################################################
jmaxCT_lm <- lm(log(jmax27.5) ~ co2_trt*temp_trt*species, 
                 data =  subset(df, ps_pathway == "C3"), 
                na.action = na.exclude)

jmaxCT_lmg <- glm(jmax27.5 ~ co2_trt*temp_trt*species, 
                 data =  subset(df, ps_pathway == "C3"), family = "gaussian",
                 na.action = na.exclude)


df$jmax27.5 <- as.numeric(df$jmax27.5)

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

AIC(jmaxCT_lm, jmaxCT_lmg)
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


##creating letter for grouping
cld(emmeans(jmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))[1, 8]
cld(emmeans(jmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))[2, 8]
cld(emmeans(jmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))[3, 8]
cld(emmeans(jmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))[4, 8]


jmaxCT_letter <- data.frame(x= c(0.8, 1.2, 1.8, 2.2),
                             temp_trt = c('LT', 'LT', 'HT',  'HT'),
                             co2_trt = c('AC', 'EC', 'AC', 'EC'),
                             y =c(175, 205, 300, 200),
                             group = c(cld(emmeans(jmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))[1, 8],
                                       cld(emmeans(jmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))[2, 8],
                                       cld(emmeans(jmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))[3, 8],
                                       cld(emmeans(jmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"))[4, 8]))



jmaxCT_letter$letter[jmaxCT_letter$group == " 1 "] <- "a"
jmaxCT_letter$letter[jmaxCT_letter$group == " 12"] <- "ab"
jmaxCT_letter$letter[jmaxCT_letter$group == "  2"] <- "b"



#jmaxCT<- as.data.frame(emmeans(jmaxCT_lm, ~temp_trt*co2_trt))

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
summary(vcmaxGT_lm)
Anova(vcmaxGT_lm)
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



##creating letter for grouping
cld(emmeans(vcmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[1, 8]
cld(emmeans(vcmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[2, 8]
cld(emmeans(vcmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[3, 8]
cld(emmeans(vcmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[4, 8]


vcmaxGT_letter <- data.frame(x= c(0.8, 1.2, 1.8, 2.2),
                            temp_trt = c('LT', 'LT', 'HT',  'HT'),
                            co2_trt = c('AC', 'EC', 'AC', 'EC'),
                            y =c(75, 80, 150, 210),
                            group = c(cld(emmeans(vcmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[1, 8],
                                      cld(emmeans(vcmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[2, 8],
                                      cld(emmeans(vcmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[3, 8],
                                      cld(emmeans(vcmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[4, 8]))



vcmaxGT_letter$letter[vcmaxGT_letter$group == " 1 "] <- "a"
vcmaxGT_letter$letter[vcmaxGT_letter$group == " 12"] <- "ab"
vcmaxGT_letter$letter[vcmaxGT_letter$group == "  2"] <- "b"

#contrast(emmeans(vcmaxGT_lm, ~temp_trt, type = "response"), 
 #        simple =  "each", combine = TRUE )


#vcmaxGT<- as.data.frame(emmeans(vcmaxGT_lm, ~temp_trt*co2_trt))




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
summary(jmaxGT_lm)
Anova(jmaxGT_lm)
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


##creating letter for grouping
cld(emmeans(jmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[1, 8]
cld(emmeans(jmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[2, 8]
cld(emmeans(jmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[3, 8]
cld(emmeans(jmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[4, 8]


jmaxGT_letter <- data.frame(x= c(0.8, 1.2, 1.8, 2.2),
                             temp_trt = c('LT', 'LT', 'HT',  'HT'),
                             co2_trt = c('AC', 'EC', 'AC', 'EC'),
                             y =c(125, 150, 210, 240),
                             group = c(cld(emmeans(jmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[1, 8],
                                       cld(emmeans(jmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[2, 8],
                                       cld(emmeans(jmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[3, 8],
                                       cld(emmeans(jmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[4, 8]))



jmaxGT_letter$letter[vcmaxGT_letter$group == " 1 "] <- "a"
jmaxGT_letter$letter[vcmaxGT_letter$group == " 12"] <- "ab"
jmaxGT_letter$letter[vcmaxGT_letter$group == "  2"] <- "b"



contrast(emmeans(jmaxCT_lm, ~co2_trt*temp_trt, type = "response"), 
        simple =  "each", combine = TRUE )


##jmaxGT
#jmaxGT<- as.data.frame(emmeans(jmaxCT_lm, ~temp_trt*co2_trt))



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


cld(emmeans(vpmaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))[1, 8]
cld(emmeans(vpmaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))[2, 8]
cld(emmeans(vpmaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))[3, 8]
cld(emmeans(vpmaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))[4, 8]


vpmaxCT_letter <- data.frame(x= c(0.8, 1.2, 1.8, 2.2),
                               temp_trt = c('LT', 'LT', 'HT',  'HT'),
                               co2_trt = c('AC', 'EC', 'AC', 'EC'),
                               y =c(100, 145, 75, 80),
                               group = c(cld(emmeans(vpmaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))[1, 8],
                                         cld(emmeans(vpmaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))[2, 8],
                                         cld(emmeans(vpmaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))[3, 8],
                                         cld(emmeans(vpmaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))[4, 8]))



vpmaxCT_letter$letter[vpmaxCT_letter$group == " 1 "] <- "a"
vpmaxCT_letter$letter[vpmaxCT_letter$group == " 12"] <- "ab"
vpmaxCT_letter$letter[vpmaxCT_letter$group == "  2"] <- "b"


#vpmaxCT<- as.data.frame(emmeans(vpmaxLACT_lm, ~temp_trt*co2_trt))





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
summary(amaxLACT_lm)
Anova(amaxLACT_lm)
r.squaredGLMM(amaxLACT_lm)


#####################individual effect
emmeans(amaxLACT_lm, ~species)
emmeans(amaxLACT_lm, ~temp_trt*co2_trt)
emmeans(amaxLACT_lm, ~temp_trt)


# Post-hoc tests
cld(emmeans(amaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))
cld(emmeans(amaxLACT_lm, pairwise~temp_trt, type = "response"))


###creating dataframe for grouping plot
cld(emmeans(amaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))[1, 8]
cld(emmeans(amaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))[2, 8]
cld(emmeans(amaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))[3, 8]
cld(emmeans(amaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))[4, 8]


amaxCT_letter <- data.frame(x= c(0.8, 1.2, 1.8, 2.2),
                             temp_trt = c('LT', 'LT', 'HT',  'HT'),
                             co2_trt = c('AC', 'EC', 'AC', 'EC'),
                             y =c(100, 185, 75, 80),
                             group = c(cld(emmeans(amaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))[1, 8],
                                       cld(emmeans(amaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))[2, 8],
                                       cld(emmeans(amaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))[3, 8],
                                       cld(emmeans(amaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"))[4, 8]))



amaxCT_letter$letter[amaxCT_letter$group == " 1"] <- "a"
#amaxCT_letter$letter[amaxCT_letter$group == " 2"] <- "b"

##amax27.5
#amaxCT<- as.data.frame(emmeans(amaxLACT_lm, ~temp_trt*co2_trt))




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



cld(emmeans(vpmaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[1, 8]
cld(emmeans(vpmaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[2, 8]
cld(emmeans(vpmaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[3, 8]
cld(emmeans(vpmaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[4, 8]


vpmaxGT_letter <- data.frame(x= c(0.8, 1.2, 1.8, 2.2),
                             temp_trt = c('LT', 'LT', 'HT',  'HT'),
                             co2_trt = c('AC', 'EC', 'AC', 'EC'),
                             y =c(100, 100, 100, 100),
                             group = c(cld(emmeans(vpmaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[1, 8],
                                       cld(emmeans(vpmaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[2, 8],
                                       cld(emmeans(vpmaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[3, 8],
                                       cld(emmeans(vpmaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[4, 8]))



vpmaxGT_letter$letter[vpmaxGT_letter$group == " 1"] <- "a"
#vpmaxGT_letter$letter[vpmaxGT_letter$group == " 12"] <- "ab"
#vpmaxGT_letter$letter[vpmaxGT_letter$group == "  2"] <- "b"




#vpmaxGT<- as.data.frame(emmeans(vpmaxLAGT_lm, ~temp_trt*co2_trt))



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




cld(emmeans(amaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[1, 8]
cld(emmeans(amaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[2, 8]
cld(emmeans(amaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[3, 8]
cld(emmeans(amaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[4, 8]


amaxGT_letter <- data.frame(x= c(0.8, 1.2, 1.8, 2.2),
                             temp_trt = c('LT', 'LT', 'HT',  'HT'),
                             co2_trt = c('AC', 'EC', 'AC', 'EC'),
                             y =c(100, 100, 100, 100),
                             group = c(cld(emmeans(amaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[1, 8],
                                       cld(emmeans(amaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[2, 8],
                                       cld(emmeans(amaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[3, 8],
                                       cld(emmeans(amaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"))[4, 8]))



amaxGT_letter $letter[amaxGT_letter $group == " 1"] <- "a"
#vpmaxGT_letter$letter[vpmaxGT_letter$group == " 12"] <- "ab"
#vpmaxGT_letter$letter[vpmaxGT_letter$group == "  2"] <- "b"




#amaxGT<- as.data.frame(emmeans(amaxLAGT_lm, ~temp_trt*co2_trt))





amaxLAGT_lmer <- lmer(log(amaxLAGT) ~ co2_trt*temp_trt + (1|species),
                  data = subset(df, ps_pathway  ==  "C4"), 
                  na.action = na.exclude)

plot(resid(amaxLAGT_lmer) ~ fitted(amaxLAGT_lmer))
qqnorm(residuals(amaxLAGT_lmer))
qqline(residuals(amaxLAGT_lmer))
hist(residuals(amaxLAGT_lmer))
shapiro.test(residuals(amaxLAGT_lm))
outlierTest(amaxLAGT_lm)

Anova(amaxLAGT_lmer)

cld(emmeans(amaxLAGT_lmer, pairwise~co2_trt*temp_trt, type = "response"))

###############################################################
##plots###
## Add colorblind friendly palette
co2.cols <- c("#2166ac", "#b2182b")
co2.cols <- c("#2166ac", "#b2182b")
full.cols <- c("#b2182b", "#f4a582", "#2166ac", "#92c5d3")

## Create blank plot as spacer plot
blank.plot <- ggplot() + 
  theme_bw() +
  theme(panel.background = element_rect(color = "white",
                                        fill = "white"),
        panel.border = element_rect(color = "white"))







# Adjusting temp_trt levels to have LT on the left and HT on the right
df$temp_trt <- factor(df$temp_trt, levels = c("LT", "HT"))

plot <- ggplot(subset(df, ps_pathway == "C3", !is.na(jmaxGT)), aes(x = temp_trt, y = jmaxGT)) +
  geom_point(shape = 21, size = 3, alpha = 0.9, aes(fill = co2_trt)) +  # Add data points with color representing different co2_trt levels
  geom_smooth(method = "lm", se = FALSE, aes(group = co2_trt)) +  # Add regression lines for each co2_trt level
  scale_fill_discrete(labels = c("ambient", "elevated")) +
  scale_color_manual(values=c("ambient"= "#2166ac", "elevated"="#b2182b")) +
  labs(x = expression(bold("Temperature")), y = expression(bold(italic("J")["maxGT"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression("CO"["2"]*" treatment")) +  # Label the axes and legend
  scale_y_continuous(limits = c(0, 200)) +
  theme_bw(base_size = 22)




#plot for TTABSS
# Your data
# Adjusting temp_trt levels to have LT on the left and HT on the right
vcmaxCT$temp_trt <- factor(vcmaxCT$temp_trt, levels = c("LT", "HT"))

# Plot
ggplot(vcmaxCT, aes(x=temp_trt, y=emmean, color=co2_trt)) +
  geom_point(position = position_dodge(width= 0.75), size=5) +
  geom_errorbar(position = position_dodge(width= 0.75),
                aes(ymin=lower.CL, ymax=upper.CL), width=.2, linewidth = 1.5) +
  scale_color_manual(values=c("AC"="blue", "EC"="red")) +
  labs(x = expression(bold("Temperature")),
       y = expression(bold(italic("V")["cmaxCT"])),
       color = expression(bold("CO"["2"]))) +
  theme_minimal() +
  theme(text = element_text(size=12))


df$temp_trt <- factor(df$temp_trt, levels = c("LT", "HT"))




##boxplot
ggplot(subset(df, ps_pathway == "C3", !is.na(vcmax27.5)), aes(x=temp_trt, color =co2_trt)) +
  geom_boxplot(aes(y = vcmax27.5, color = co2_trt), alpha = 0.5) +
  geom_point(aes(y = vcmaxCT, color = co2_trt), alpha = 0.5,
             position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.3)) +
  geom_point(data = vcmaxCT, aes(x = temp_trt, y = emmean, group = co2_trt),
             size = 3, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = vcmaxCT,
                aes(x = temp_trt, ymin=lower.CL, ymax=upper.CL, group = co2_trt),
                width = 0.2, stat = "identity", linewidth = 1, 
                position = position_dodge(width = 0.75)) + 
  scale_color_manual(labels = c("ambient", "elevated"),
                     values = c("blue", "red")) +
  theme_minimal()






#write.csv(df, "~/git/C3C4_GrowthChamber_Experiment/TxCO2_datasheets/combine.csv", row.names = FALSE)







###############plots################
#vpmax27.5
##################################
ggplot(data = df, aes(x= temp_trt, y = vpmaxLAGT, fill = co2_trt)) +
  geom_boxplot() 

ggplot(data = df, aes(x= temp_trt, y = vpmaxLAGT, fill = co2_trt)) +
  geom_boxplot(outlier.shape = NA) 

ggplot(data = df, aes(x= temp_trt, y = amaxLA27.5, fill = co2_trt)) +
  geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits = c(0, 200))



# Generate a sequence of Ci values
ci_sequence <- seq(from = 5, to = 50, by = 5)  # Adjust range and step size as needed

min_amax <- min(df$vcmaxGT, na.rm = TRUE)
max_amax <- max(df$vcmaxGT, na.rm = TRUE)

prediction_df <- data.frame(
  ci = ci_sequence,
  amax = runif(length(ci_sequence), min = min_amax, max = max_amax)
)



ggplot(data =  subset(df, ps_pathway == "C3"),  aes(x= temp_trt, y = jmax27.5, fill = co2_trt)) +
  geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits = c(0, 200))

plot(df$amaxLA27.5~df$vpmax27.5)
