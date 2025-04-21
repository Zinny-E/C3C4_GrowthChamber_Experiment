library(ggplot2)
library(tidyverse)
library(ggpubr)
library(gridExtra)

##creating letter for grouping used for making plots
#####################vcmaxCT######################

vcmaxCT_letter <- cld(emmeans(vcmaxCT_lm, pairwise~co2_trt*temp_trt), Letters = letters) %>%
  mutate(.group = trimws(.group, "both")) %>% data.frame()


############################vcmaxGT########################
vcmaxGT_letter <- cld(emmeans(vcmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"),
                      Letters = letters) %>%
  mutate(.group = trimws(.group, "both")) %>% data.frame()


##########################jmaxCT####################
##creating letter for grouping
jmaxCT_letter <- cld(emmeans(jmaxCT_lm, pairwise~co2_trt*temp_trt, type = "response"), 
                     Letters = letters) %>%
  mutate(.group = trimws(.group, "both")) %>% data.frame()


###################jmaxGT######################
jmaxGT_letter <- cld(emmeans(jmaxGT_lm, pairwise~co2_trt*temp_trt, type = "response"), 
                     Letters = letters) %>%
  mutate(.group = trimws(.group, "both")) %>% data.frame()


########################vpmaxCT###############
vpmaxCT_letter <- cld(emmeans(vpmaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"),
                      Letters = letters) %>%
  mutate(.group = trimws(.group, "both")) %>% data.frame()


####################vpmaxGT######################
vpmaxGT_letter <- cld(emmeans(vpmaxLAGT_lm, pairwise~co2_trt*temp_trt),
                      Letters = letters) %>%
  mutate(.group = trimws(.group, "both")) %>% data.frame()


######################amaxCT######################
amaxCT_letter <- cld(emmeans(amaxLACT_lm, pairwise~co2_trt*temp_trt, type = "response"),
                      Letters = letters) %>%
  mutate(.group = trimws(.group, "both")) %>% data.frame()


####################amaxGT#################
amaxGT_letter <- cld(emmeans(amaxLAGT_lm, pairwise~co2_trt*temp_trt, type = "response"),
                      Letters = letters) %>%
  mutate(.group = trimws(.group, "both")) %>% data.frame()


########changing temp_trt level##
# Adjusting temp_trt levels to have LT on the left and HT on the right
df$temp_trt <- factor(df$temp_trt, levels = c("LT", "HT"))

##########################plots###
figtheme <- theme_minimal(base_size = 18) +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(size = 1.5, fill = NA),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background=element_blank(),
        legend.title = element_text(face = "bold"),
        axis.ticks.length = unit(0.25, "cm"),
        panel.grid.minor.y = element_blank(),
        legend.text.align = 0)

############### C3 plots################
##################################################
vcmaxCT.plot <- ggplot(subset(df, ps_pathway == "C3", !is.na(vcmax27.5)), 
                       aes(x=temp_trt, color =co2_trt)) +
  geom_jitter(aes(y = vcmax27.5, color = co2_trt), alpha = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  geom_point(data = vcmaxCT_letter, aes(x = temp_trt, y = emmean, group = co2_trt),
             size = 5, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = vcmaxCT_letter,
                aes(x = temp_trt, ymin=lower.CL, ymax=upper.CL, group = co2_trt),
                width = 0.23, stat = "identity", linewidth = 1, 
                position = position_dodge(width = 0.75)) + 
  geom_text(data = vcmaxCT_letter, aes(x= temp_trt, y= 200, label = .group), 
            position = position_dodge(width = 0.75))  +
  scale_color_manual(labels = c("ambient", "elevated"),
                     values = c("#2c46a5", "#FF9200")) +
  labs(x = expression(bold("Temperature")), 
       y = expression(bold(italic("V")["cmax27.5"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression("CO"["2"]*" treatment")) +  # Label the axes and legend
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))

vcmaxCT.plot


vcmaxGT.plot <- ggplot(subset(df, ps_pathway == "C3", !is.na(vcmaxGT)), 
                       aes(x=temp_trt, color =co2_trt)) +
  geom_jitter(aes(y = vcmaxGT, color = co2_trt), alpha = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  geom_point(data = vcmaxGT_letter, aes(x = temp_trt, y = response, group = co2_trt),
             size = 5, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = vcmaxGT_letter,
                aes(x = temp_trt, ymin=lower.CL, ymax=upper.CL, group = co2_trt),
                width = 0.23, stat = "identity", linewidth = 1, 
                position = position_dodge(width = 0.75)) + 
  geom_text(data = vcmaxGT_letter, aes(x= temp_trt, y= 350, label = .group), 
            position = position_dodge(width = 0.75) ) +
  scale_color_manual(labels = c("ambient", "elevated"),
                     values = c("#2c46a5", "#FF9200")) +
  labs(x = expression(bold("Temperature")), 
       y = expression(bold(italic("V")["cmaxGT"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression("CO"["2"]*" treatment")) +  # Label the axes and legend
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))

vcmaxGT.plot


###combining the plots using ggarrange#####################
# png("[insert path here]",
#     width = 13, height = 10, units = 'in', res = 600)
ggarrange(vcmaxCT.plot, vcmaxGT.plot,
          align = "h", common.legend = T, legend = "right",
          labels = "AUTO", vjust = 1.5,
          font.label = list(size = 18, face = "bold"))

# dev.off()



jmaxCT.plot <- ggplot(subset(df, ps_pathway == "C3", !is.na(jmax27.5)), 
                      aes(x=temp_trt, color =co2_trt)) +
  geom_jitter(aes(y = jmax27.5, color = co2_trt), alpha = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  geom_point(data = jmaxCT_letter, aes(x = temp_trt, y = response, group = co2_trt),
             size = 5, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = jmaxCT_letter,
                aes(x = temp_trt, ymin=lower.CL, ymax=upper.CL, group = co2_trt),
                width = 0.2, stat = "identity", linewidth = 1, 
                position = position_dodge(width = 0.75)) + 
  geom_text(data = jmaxCT_letter, aes(x = temp_trt, y = 350, label = .group),
            position = position_dodge(width = 0.75)) +
  scale_color_manual(labels = c("ambient", "elevated"),
                     values = c("#2c46a5", "#FF9200")) +
  labs(x = expression(bold("Temperature")), 
       y = expression(bold(italic("J")["max27.5"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression("CO"["2"]*" treatment")) +  # Label the axes and legend
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))

jmaxCT.plot


jmaxGT.plot <- ggplot(subset(df, ps_pathway == "C3", !is.na(jmaxGT)), 
                      aes(x=temp_trt, color =co2_trt)) +
  geom_jitter(aes(y = jmaxGT, color = co2_trt), alpha = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  geom_point(data = jmaxGT_letter, aes(x = temp_trt, y = response, group = co2_trt),
             size = 5, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = jmaxGT_letter,
                aes(x = temp_trt, ymin=lower.CL, ymax=upper.CL, group = co2_trt),
                width = 0.2, stat = "identity", linewidth = 1, 
                position = position_dodge(width = 0.75)) + 
  geom_text(data = jmaxGT_letter, aes(x= temp_trt, y= 400, label = .group), 
            position = position_dodge(width = 0.75)) +
  scale_color_manual(labels = c("ambient", "elevated"),
                     values = c("#2c46a5", "#FF9200")) +
  scale_y_continuous(limits = c(0, 450)) +
  labs(x = expression(bold("Temperature")), 
       y = expression(bold(italic("J")["maxGT"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression("CO"["2"]*" treatment")) +  # Label the axes and legend
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))

jmaxGT.plot


###combining the plots using ggarrange#####################
# png("[insert path here]",
#     width = 13, height = 10, units = 'in', res = 600)
ggarrange(jmaxCT.plot, jmaxGT.plot,
          align = "h", common.legend = T, legend = "right",
          labels = "AUTO", vjust = 1.5,
          font.label = list(size = 18, face = "bold"))

# dev.off()



############### C4 plots################
##################################################
vpmaxCT.plot <- ggplot(subset(df, ps_pathway == "C4", !is.na(vpmaxLA27.5)), 
                       aes(x=temp_trt, color =co2_trt)) +
  geom_jitter(aes(y = vpmaxLA27.5, color = co2_trt), alpha = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  geom_point(data = vpmaxCT_letter, aes(x = temp_trt, y = response, group = co2_trt),
             size = 5, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = vpmaxCT_letter,
                aes(x = temp_trt, ymin=lower.CL, ymax=upper.CL, group = co2_trt),
                width = 0.2, stat = "identity", linewidth = 1, 
                position = position_dodge(width = 0.75)) + 
  geom_text(data = vpmaxCT_letter, aes(x= temp_trt, y= 200, label = .group), 
            position = position_dodge(width = 0.75))  +
  scale_color_manual(labels = c("ambient", "elevated"),
                     values = c("#2c46a5", "#FF9200")) +
  labs(x = expression(bold("Temperature")), 
       y = expression(bold(italic("V")["pmax27.5"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression("CO"["2"]*" treatment")) +  # Label the axes and legend
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))

vpmaxCT.plot



vpmaxGT.plot <- ggplot(subset(df, ps_pathway == "C4", !is.na(vpmaxLAGT)), 
                       aes(x=temp_trt, color =co2_trt)) +
  geom_jitter(aes(y = vpmaxLAGT, color = co2_trt), alpha = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  geom_point(data = vpmaxGT_letter, aes(x = temp_trt, y = emmean, group = co2_trt),
             size = 5, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = vpmaxGT_letter,
                aes(x = temp_trt, ymin=lower.CL, ymax=upper.CL, group = co2_trt),
                width = 0.2, stat = "identity", linewidth = 1, 
                position = position_dodge(width = 0.75)) + 
  geom_text(data = vpmaxGT_letter, aes(x= temp_trt, y= 200, label = .group), 
            position = position_dodge(width = 0.75))  +
  scale_color_manual(labels = c("ambient", "elevated"),
                     values = c("#2c46a5", "#FF9200")) +
  labs(x = expression(bold("Temperature")), 
       y = expression(bold(italic("V")["pmaxGT"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression("CO"["2"]*" treatment")) +  # Label the axes and legend
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))

vpmaxGT.plot


###combining the plots using ggarrange#####################
# png("[insert path here]",
#     width = 13, height = 10, units = 'in', res = 600)
ggarrange(vpmaxCT.plot, vpmaxGT.plot,
          align = "h", common.legend = T, legend = "right",
          labels = "AUTO", vjust = 1.5,
          font.label = list(size = 18, face = "bold"))

# dev.off()




amaxCT.plot <- ggplot(subset(df, ps_pathway == "C4", !is.na(amaxLA27.5)), 
                       aes(x=temp_trt, color =co2_trt)) +
  geom_jitter(aes(y = amaxLA27.5, color = co2_trt), alpha = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  geom_point(data = amaxCT_letter, aes(x = temp_trt, y = response, group = co2_trt),
             size = 3.5, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = amaxCT_letter,
                aes(x = temp_trt, ymin=lower.CL, ymax=upper.CL, group = co2_trt),
                width = 0.2, stat = "identity", linewidth = 1, 
                position = position_dodge(width = 0.75)) + 
  geom_text(data = amaxCT_letter, aes(x= temp_trt, y= 200, label = .group), 
            position = position_dodge(width = 0.75))  +
  scale_color_manual(labels = c("ambient", "elevated"),
                     values = c("#2c46a5", "#FF9200")) +
  labs(x = expression(bold("Temperature")), 
       y = expression(bold(italic("A")["max27.5"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression("CO"["2"]*" treatment")) +  # Label the axes and legend
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))

amaxCT.plot



amaxGT.plot <- ggplot(subset(df, ps_pathway == "C4", !is.na(amaxLAGT)), 
                      aes(x=temp_trt, color =co2_trt)) +
  geom_jitter(aes(y = amaxLAGT, color = co2_trt), alpha = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  geom_point(data = amaxGT_letter, aes(x = temp_trt, y = response, group = co2_trt),
             size = 3.5, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = amaxGT_letter,
                aes(x = temp_trt, ymin=lower.CL, ymax=upper.CL, group = co2_trt),
                width = 0.2, stat = "identity", linewidth = 1, 
                position = position_dodge(width = 0.75)) + 
  geom_text(data = amaxGT_letter, aes(x= temp_trt, y= 200, label = .group), 
            position = position_dodge(width = 0.75))  +
  scale_color_manual(labels = c("ambient", "elevated"),
                     values = c("#2c46a5", "#FF9200")) +
  labs(x = expression(bold("Temperature")), 
       y = expression(bold(italic("A")["maxGT"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression("CO"["2"]*" treatment")) +  # Label the axes and legend
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))

amaxGT.plot


###combining the plots using ggarrange#####################
# png("[insert path here]",
#     width = 13, height = 10, units = 'in', res = 600)
ggarrange(amaxCT.plot, amaxGT.plot,
          align = "h", common.legend = T, legend = "right",
          labels = "AUTO", vjust = 1.5,
          font.label = list(size = 18, face = "bold"))

# dev.off()




















###vcmac27.5 boxplot
ggplot(subset(df, ps_pathway == "C3", !is.na(vcmax27.5)), 
       aes(x=temp_trt, y=vcmax27.5, fill = co2_trt)) +
  geom_boxplot() +
  geom_text(data = vcmaxCT_letter, aes(x= x, y= y, label = letter), size = 6) +
  scale_fill_manual(values=c("#2166ac", "#b2182b"), 
                    labels = c("ambient", "elevated"))