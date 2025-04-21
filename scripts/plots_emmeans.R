library(ggplot2)
library(tidyverse)
library(ggpubr)
library(gridExtra)



########changing temp_trt level##
# Adjusting temp_trt levels to have LT on the left and HT on the right
df$temp_trt <- factor(df$temp_trt, levels = c("LT", "HT"))
df$co2_trt <- factor(df$co2_trt, levels = c("AC", "EC"))



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




vcmaxCTSxCxTletter <- cld(emmeans(vcmaxCT_lm, pairwise~co2_trt*temp_trt*species), 
                         Letters = letters) %>%
  mutate(.group = trimws(.group, "both")) %>% data.frame()



ggplot(subset(df, ps_pathway == "C3", !is.na(vcmax27.5)), 
       aes(x = species, y = vcmax27.5, color = temp_trt, shape = co2_trt)) +  # Combine temp & CO2
  geom_jitter(aes(y = vcmax27.5, color = species), alpha = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  geom_point(data = vcmaxCTSxCxTletter, aes(x = interaction(temp_trt, co2_trt), y = emmean, group = species),
             size = 5, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = vcmaxCTSxCxTletter,
                aes(x = interaction(temp_trt, co2_trt), ymin = lower.CL, ymax = upper.CL, group = species),
                width = 0.23, stat = "identity", linewidth = 1, 
                position = position_dodge(width = 0.75)) + 
  geom_text(data = vcmaxCTSxCxTletter, aes(x = interaction(temp_trt, co2_trt), y = 350, label = .group), 
            position = position_dodge(width = 0.75), size = 5) +  # Adjust text size
  scale_color_manual(labels = c("Elymus", "Pascopyrum", "Poa"),
                     values = c("#018571", "#E66101", "#5E3C99")) +
  labs(x = expression(bold("Temperature & CO₂ Treatment")), 
       y = expression(bold(italic("Vc")["max27.5"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression("Species")) +  # Label the axes and legend
  figtheme +
  theme_minimal() +  # Ensure a clean layout
  theme(
    text = element_text(size=12),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(hjust = 0.5)
  )
##creating letter for grouping used for making plots
#####################anetgrowth######################

anetgrowth_letter <- cld(emmeans(anetgrowth_lm, pairwise~co2_trt*temp_trt),
                         Letters = letters) %>%
  mutate(.group = trimws(.group, "both")) %>% data.frame()


anetgrowth.plot <- ggplot(subset(df, !is.na(anet_growth)), 
                          aes(x=temp_trt, color =co2_trt)) +
  geom_jitter(aes(y = anet_growth, color = co2_trt), alpha = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  geom_point(data = anetgrowth_letter, aes(x = temp_trt, y = emmean, group = co2_trt),
             size = 5, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = anetgrowth_letter,
                aes(x = temp_trt, ymin=lower.CL, ymax=upper.CL, group = co2_trt),
                width = 0.23, stat = "identity", linewidth = 1, 
                position = position_dodge(width = 0.75)) + 
  geom_text(data = anetgrowth_letter, aes(x= temp_trt, y= 75, label = .group), 
            position = position_dodge(width = 0.75))  +
  scale_color_manual(labels = c("ambient", "elevated"),
                     values = c("#2c46a5", "#FF9200")) +
  labs(x = expression(bold("Temperature")), 
       y = expression(bold(italic("A")["netgrowth"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression("CO"["2"]*" treatment")) +  # Label the axes and legend
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))

anetgrowth.plot


######interaction###
##creating letter for grouping
anetgrowthCxS_letter <- cld(emmeans(anetgrowth_lm, pairwise~temp_trt*species),
                            Letters = letters) %>%
  mutate(.group = trimws(.group, "both")) %>% data.frame()


anetgrowthCxS.plot <- ggplot(subset(df, !is.na(anet_growth)), 
                             aes(x = temp_trt, color = species)) +
  geom_jitter(aes(y = anet_growth), alpha = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  geom_point(data = anetgrowthCxS_letter, aes(x = temp_trt, y = emmean, group = species),
             size = 5, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = anetgrowthCxS_letter,
                aes(x = temp_trt, ymin = lower.CL, ymax = upper.CL, group = species),
                width = 0.23, linewidth = 1, 
                position = position_dodge(width = 0.75)) + 
  geom_text(data = anetgrowthCxS_letter, aes(x= temp_trt, y= 75, label = .group), 
            position = position_dodge(width = 0.75))  +
  scale_color_manual(
    values = c(
      "Bouteloua Curtipendula" = "#018571",
      "Pascopyrum Smithii" = "#E66101",
      "Schizachyrium Scoparium" = "#5E3C99",
      "Elymus Canadensis" = "#8B8B00",
      "Poa Pratensis" = "#8B0000",
      "Sorghastrum nutans" = "#556B2F"
    )
  ) +
  labs(x = expression(bold("Temperature")), 
       y = expression(bold(italic("A")["netgrowth"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression("species")) +  # Label the axes and legend
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))

anetgrowthCxS.plot




###################################################################
###gsw#########
gsw_letter <- cld(emmeans(gsw_lm, pairwise~co2_trt*temp_trt, type = "response"),
                  Letters = letters) %>%
  mutate(.group = trimws(.group, "both")) %>% data.frame()


gsw.plot <- ggplot(subset(df, !is.na(gsw)), 
                   aes(x=temp_trt, color =co2_trt)) +
  geom_jitter(aes(y = gsw, color = co2_trt), alpha = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  geom_point(data = gsw_letter , aes(x = temp_trt, y = response, group = co2_trt),
             size = 5, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = gsw_letter,
                aes(x = temp_trt, ymin=lower.CL, ymax=upper.CL, group = co2_trt),
                width = 0.23,linewidth = 1, 
                position = position_dodge(width = 0.75)) + 
  geom_text(data = gsw_letter, aes(x= temp_trt, y= 0.75, label = .group), 
            position = position_dodge(width = 0.75))  +
  scale_color_manual(labels = c("ambient", "elevated"),
                     values = c("#2c46a5", "#FF9200")) +
  labs(x = expression(bold("Temperature")), 
       y = expression(bold(italic("G")["sw"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression("CO"["2"]*" treatment")) +  # Label the axes and legend
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))

gsw.plot




####################################################Nmass

nmass_letter <- cld(emmeans(nmass_lm, pairwise~co2_trt, type = "response"),
                    Letters = letters) %>%
  mutate(.group = trimws(.group, "both")) %>% data.frame()


nmass_plot <- ggplot(subset(df, !is.na(nmass)), 
                     aes(x=co2_trt, color =co2_trt)) +
  geom_jitter(aes(y = nmass, color = co2_trt), alpha = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  geom_point(data = nmass_letter , aes(x = co2_trt, y = response, group = co2_trt),
             size = 5, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = nmass_letter,
                aes(x = co2_trt, ymin=lower.CL, ymax=upper.CL, group = co2_trt),
                width = 0.23,linewidth = 1, 
                position = position_dodge(width = 0.75)) + 
  geom_text(data = nmass_letter, aes(x= co2_trt, y = 0.05, label = .group), 
            position = position_dodge(width = 0.75))  +
  scale_color_manual(labels = c("ambient", "elevated"),
                     values = c("#2c46a5", "#FF9200")) +
  labs(x = expression(bold("Temperature")), 
       y = expression(bold(italic("N")["mass"]*" ("*g*"N"^"-2"*" g"^"-1"*")")),
       fill = expression("CO"["2"]*" treatment")) +  # Label the axes and legend
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))

nmass_plot





nmassCxS_letter <- cld(emmeans(nmass_lm, pairwise~co2_trt*species, type = "response"),
                       Letters = letters) %>%
  mutate(.group = trimws(.group, "both")) %>% data.frame()




nmassCxS.plot <- ggplot(subset(df, !is.na(nmass)), 
                        aes(x = co2_trt, color = species)) +
  geom_jitter(aes(y = nmass), alpha = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  geom_point(data = nmass_letter, aes(x = co2_trt, y = response, group = species),
             size = 5, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = nmass_letter,
                aes(x = co2_trt, ymin = lower.CL, ymax = upper.CL, group = species),
                width = 0.23, linewidth = 1, 
                position = position_dodge(width = 0.75)) +  # Ensure error bars align
  geom_text(data = nmass_letter, 
            aes(x = co2_trt, y = max(nmass_letter$upper.CL) + 0.01, label = .group),  # Dynamically place letters on top
            position = position_dodge(width = 0.75)) +  # Make text larger for readability
  scale_color_manual(
    values = c(
      "Bouteloua Curtipendula" = "#018571",
      "Pascopyrum Smithii" = "#E66101",
      "Schizachyrium Scoparium" = "#5E3C99",
      "Elymus Canadensis" = "#8B8B00",
      "Poa Pratensis" = "#8B0000",
      "Sorghastrum nutans" = "#556B2F"
    )
  ) +
  labs(x = expression(bold("CO2 treatment")), 
       y = expression(bold(italic("N")["mass"]*" ("*g*"N"^"-2"*" g"^"-1"*")")),
       fill = expression("species")) +  # Label the axes and legend
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))

nmassCxS.plot

#########################################################################
##############narea
narea_letter <- cld(emmeans(narea_lm, pairwise~co2_trt*temp_trt, type = "response"),
                    Letters = letters) %>%
  mutate(.group = trimws(.group, "both")) %>% data.frame()


narea.plot <- ggplot(subset(df, !is.na(narea)), 
                     aes(x=temp_trt, color =co2_trt)) +
  geom_jitter(aes(y = narea, color = co2_trt), alpha = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  geom_point(data =narea_letter , aes(x = temp_trt, y = response, group = co2_trt),
             size = 5, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = narea_letter,
                aes(x = temp_trt, ymin=lower.CL, ymax=upper.CL, group = co2_trt),
                width = 0.23,linewidth = 1, 
                position = position_dodge(width = 0.75)) + 
  geom_text(data = narea_letter, aes(x= temp_trt, y= 1.5, label = .group), 
            position = position_dodge(width = 0.75))  +
  scale_color_manual(labels = c("ambient", "elevated"),
                     values = c("#2c46a5", "#FF9200")) +
  labs(x = expression(bold("Temperature")), 
       y = expression(bold(italic("N")["area"]*" ("*g*"m"^"-2"*")")),
       fill = expression("CO"["2"]*" treatment")) +  # Label the axes and legend
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))

narea.plot


#########################################################################################
###Marea

marea_letter <- cld(emmeans(marea_lm, pairwise~co2_trt*temp_trt, type = "response"),
                    Letters = letters) %>%
  mutate(.group = trimws(.group, "both")) %>% data.frame()


marea.plot <- ggplot(subset(df, !is.na(marea)), 
                     aes(x=temp_trt, color =co2_trt)) +
  geom_jitter(aes(y = marea, color = co2_trt), alpha = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  geom_point(data =marea_letter , aes(x = temp_trt, y = response, group = co2_trt),
             size = 5, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = marea_letter,
                aes(x = temp_trt, ymin=lower.CL, ymax=upper.CL, group = co2_trt),
                width = 0.23,linewidth = 1, 
                position = position_dodge(width = 0.75)) + 
  geom_text(data = marea_letter, aes(x= temp_trt, y= 100, label = .group), 
            position = position_dodge(width = 0.75))  +
  scale_color_manual(labels = c("ambient", "elevated"),
                     values = c("#2c46a5", "#FF9200")) +
  labs(x = expression(bold("Temperature")), 
       y = expression(bold(italic("M")["area"]*" ("*g*"m"^"-2"*")")),
       fill = expression("CO"["2"]*" treatment")) +  # Label the axes and legend
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))

marea.plot



################################################################
###biomass
biomass_letter <- cld(emmeans(biomass_lm, pairwise~co2_trt*temp_trt, type = "response"),
                      Letters = letters) %>%
  mutate(.group = trimws(.group, "both")) %>% data.frame()


biomass_plot <- ggplot(subset(df, !is.na(biomass.g.)), 
                       aes(x=temp_trt, color =co2_trt)) +
  geom_jitter(aes(y = biomass.g., color = co2_trt), alpha = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  geom_point(data = biomass_letter, aes(x = temp_trt, y = response, group = co2_trt),
             size = 5, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = biomass_letter,
                aes(x = temp_trt, ymin=lower.CL, ymax=upper.CL, group = co2_trt),
                width = 0.23, stat = "identity", linewidth = 1, 
                position = position_dodge(width = 0.75)) + 
  geom_text(data = biomass_letter, aes(x= temp_trt, y= 15, label = .group), 
            position = position_dodge(width = 0.75))  +
  scale_color_manual(labels = c("ambient", "elevated"),
                     values = c("#2c46a5", "#FF9200")) +
  labs(x = expression(bold("Temperature")), 
       y = expression(bold("Biomass (g)")),
       fill = expression("CO"["2"]*" treatment")) +  # Label the axes and legend
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))

biomass_plot





biomassCxS_letter <- cld(emmeans(biomass_lm, pairwise~temp_trt*species, type = "response"),
                         Letters = letters) %>%
  mutate(.group = trimws(.group, "both")) %>% data.frame()


biomassCxS.plot <- ggplot(subset(df, !is.na(biomass.g.)), 
                          aes(x = temp_trt, color = species)) +
  geom_jitter(aes(y = biomass.g.), alpha = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  geom_point(data = biomassCxS_letter, aes(x = temp_trt, y = response, group = species),
             size = 5, shape = 15, position = position_dodge(width = 0.75)) +
  geom_errorbar(data = biomassCxS_letter,
                aes(x = temp_trt, ymin = lower.CL, ymax = upper.CL, group = species),
                width = 0.23, linewidth = 1, 
                position = position_dodge(width = 0.75)) + 
  geom_text(data = biomassCxS_letter, aes(x= temp_trt, y= 15, label = .group), 
            position = position_dodge(width = 0.75))  +
  scale_color_manual(
    values = c(
      "Bouteloua Curtipendula" = "#018571",
      "Pascopyrum Smithii" = "#E66101",
      "Schizachyrium Scoparium" = "#5E3C99",
      "Elymus Canadensis" = "#8B8B00",
      "Poa Pratensis" = "#8B0000",
      "Sorghastrum nutans" = "#556B2F"
    )
  ) +
  labs(x = expression(bold("Temperature")), 
       y = expression(bold(italic("Biomass")["g"])),
       fill = expression("species")) +  # Label the axes and legend
  figtheme +
  theme(text = element_text(size=12)) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0.5))

biomassCxS.plot

