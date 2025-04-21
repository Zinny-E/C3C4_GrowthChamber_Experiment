#loading packages
library(ggplot2)
library(ggpubr)

#####plotting by trt by species
spp.vcmaxCT <-ggplot(subset(df, ps_pathway == "C3", !is.na(vcmax27.5)), 
       aes(x=temp_trt, y=vcmax27.5, fill = co2_trt)) +
  geom_boxplot() + 
  facet_wrap(vars(species))


spp.jmaxCT <-ggplot(subset(df, ps_pathway == "C3", !is.na(jmax27.5)), 
       aes(x=temp_trt, y=jmax27.5, fill=co2_trt)) +
  geom_boxplot() +
  facet_wrap(vars(species))

spp.vcmaxGT <-ggplot(subset(df, ps_pathway == "C3", !is.na(vcmaxGT)), 
       aes(x=temp_trt, y=vcmaxGT, fill=co2_trt)) +
  geom_boxplot() +
  facet_wrap(vars(species))

spp.jmaxGT <-ggplot(subset(df, ps_pathway == "C3", !is.na(jmaxGT)), 
       aes(x=temp_trt, y=jmaxGT, fill=co2_trt)) +
  geom_boxplot() +
  facet_wrap(vars(species)) +
  scale_y_continuous(limits = c(0, 450)) # trying to take off outlier, not the best might have to consider another way to deal with it

###arrange plots
ggarrange(spp.vcmaxCT, spp.jmaxCT,  spp.vcmaxGT, spp.jmaxGT, ncol = 2, nrow = 2)

####C4 species
spp.vpmaxCT <- ggplot(subset(df, ps_pathway == "C4", !is.na(vpmaxLA27.5)), 
       aes(x=temp_trt, y=vpmaxLA27.5, fill=co2_trt)) +
  geom_boxplot() +
  facet_wrap(vars(species))

spp.amaxCT <- ggplot(subset(df, ps_pathway == "C4", !is.na(amaxLA27.5)), 
       aes(x=temp_trt, y=amaxLA27.5, fill=co2_trt)) +
  geom_boxplot() +
  facet_wrap(vars(species))

spp.vpmaxGT <- ggplot(subset(df, ps_pathway == "C4", !is.na(vpmaxLAGT)), 
       aes(x=temp_trt, y=vpmaxLAGT, fill=co2_trt)) +
  geom_boxplot() +
  facet_wrap(vars(species))


spp.amaxGT <-ggplot(subset(df, ps_pathway == "C4", !is.na(amaxLAGT)),
       aes(x=temp_trt, y=amaxLAGT, fill=co2_trt)) +
  geom_boxplot() +
  facet_wrap(vars(species))


###arrange plots
ggarrange(spp.vpmaxCT, spp.amaxCT,  spp.vpmaxGT, spp.amaxGT, ncol = 2, nrow = 2)

########################vpmaxCT###############
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





df$temp_trt <- factor(df$temp_trt, levels = c("LT", "HT"))

#####################################################
#########plotting the model result####################
##############C3######################################
ggplot(subset(df, ps_pathway == "C3", !is.na(vcmax27.5)), 
       aes(x=temp_trt, y=vcmax27.5, fill = co2_trt)) +
  geom_boxplot() +
  geom_text(data = vcmaxCT_letter, aes(x= x, y= y, label = letter), size = 6) +
  scale_fill_manual(values=c("#2166ac", "#b2182b"), 
                    labels = c("ambient", "elevated"))


ggplot(subset(df, ps_pathway == "C3", !is.na(jmax27.5)), 
       aes(x=temp_trt, y=jmax27.5, fill=co2_trt)) +
  geom_boxplot() +
  geom_text(data = jmaxCT_letter, aes(x= x, y= y, label = letter), size = 6) +
  scale_fill_manual(values=c("#2166ac", "#b2182b"), 
                    labels = c("ambient", "elevated"))

ggplot(subset(df, ps_pathway == "C3", !is.na(vcmaxGT)), 
       aes(x=temp_trt, y=vcmaxGT, fill=co2_trt)) +
  geom_boxplot() +
  geom_text(data = vcmaxGT_letter, aes(x= x, y= y, label = letter), size = 6) +
  scale_fill_manual(values=c("#2166ac", "#b2182b"), 
                    labels = c("ambient", "elevated"))

ggplot(subset(df, ps_pathway == "C3", !is.na(jmaxGT)), 
       aes(x=temp_trt, y=jmaxGT, fill=co2_trt)) +
  geom_boxplot() +
  geom_text(data = jmaxGT_letter, aes(x= x, y= y, label = letter), size = 6) +
  scale_fill_manual(values=c("#2166ac", "#b2182b"), 
                    labels = c("ambient", "elevated"))






############################################################
####C4########################
ggplot(subset(df, ps_pathway == "C4", !is.na(vpmaxLA27.5)), 
       aes(x=temp_trt, y=vpmaxLA27.5, fill=co2_trt)) +
  geom_boxplot() +
  geom_text(data = vpmaxCT_letter, aes(x= x, y= y, label = letter), size = 6) +
  scale_fill_manual(values=c("#2166ac", "#b2182b"), 
                    labels = c("ambient", "elevated"))



ggplot(subset(df, ps_pathway == "C4", !is.na(amaxLA27.5)), 
       aes(x=temp_trt, y=amaxLA27.5, fill=co2_trt)) +
  geom_boxplot() +
  geom_text(data = amaxCT_letter, aes(x= x, y= y, label = letter), size = 6) +
  scale_fill_manual(values=c("#2166ac", "#b2182b"), 
                    labels = c("ambient", "elevated"))
  

ggplot(subset(df, ps_pathway == "C4", !is.na(vpmaxLAGT)),
       aes(x=temp_trt, y=vpmaxLAGT, fill=co2_trt)) +
  geom_boxplot() +
  geom_text(data = vpmaxGT_letter, aes(x= x, y= y, label = letter), size = 6) +
  scale_fill_manual(values=c("#2166ac", "#b2182b"), 
                    labels = c("ambient", "elevated"))


ggplot(subset(df, ps_pathway == "C4", !is.na(amaxLAGT)), 
       aes(x=temp_trt, y=amaxLAGT, fill=co2_trt)) +
  geom_boxplot() +
  geom_text(data = amaxGT_letter, aes(x= x, y= y, label = letter), size = 6) +
  scale_fill_manual(values=c("#2166ac", "#b2182b"), 
                    labels = c("ambient", "elevated"))
