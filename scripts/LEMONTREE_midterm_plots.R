# Load libraries
library(ggplot2)
library(dplyr)
library(ggpubr)

##setwd(#"/Users/eaperkowski/git/optimal_vcmax_R/")

setwd("~/git_repo/optimal_vcmax_R")

# Load optimal vcmax R script and functions therein
source("calc_optimal_vcmax.R")
sourceDirectory("functions")

# Do optimality model simulations for C4 species across range in temps and 
# binned into two CO2 treatments
df.c4.420 <- calc_optimal_vcmax(pathway = "C4", tg_c = seq(20, 35, 1), cao = 420)
df.c4.1000 <- calc_optimal_vcmax(pathway = "C4", tg_c = seq(20, 35, 1), cao = 1000)

# Do optimality model simulations for C3 species across range in temps and 
# binned into two CO2 treatments
df.c3.420 <- calc_optimal_vcmax(pathway = "C3", tg_c = seq(20, 35, 1), cao = 420)
df.c3.1000 <- calc_optimal_vcmax(pathway = "C3", tg_c = seq(20, 35, 1), cao = 1000)

# Merge all simulations into central data frame
df.total <- df.c3.420 %>% full_join(df.c3.1000) %>% full_join(df.c4.420) %>%
  full_join(df.c4.1000)

# Net photosynthesis plot
ac.plot <- ggplot(data = df.total, 
                  aes(x = tg_c, y = Ac, color = pathway, linetype = factor(cao))) +
  geom_line(size = 2) +
  scale_color_manual(values = c("blue", "red")) +  # Specify colors
  scale_linetype_manual(values = c("dotted", "solid")) +  # Specify linetypes
  labs(x = expression("Growth temperature ("*degree*"C)"),
       y = "Photosynthesis",
       color = "Photo. pathway",
       linetype = expression("CO"[2])) +
  theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))


# Vcmax plot
vcmax.plot <- ggplot(data = df.total, 
       aes(x = tg_c, y = vcmax25, color = pathway, linetype = factor(cao))) +
  geom_line(size = 2) +
  scale_color_manual(values = c("blue", "red")) +  # Specify colors
  scale_linetype_manual(values = c("dotted", "solid")) +  # Specify linetypes
  labs(x = expression("Growth temperature ("*degree*"C)"),
       y = expression(bold(italic("V")["cmax25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       color = "Photo. pathway",
       linetype = expression("CO"[2])) +
  theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  theme(legend.position = "top")

# Jmax plot
jmax.plot <- ggplot(data = df.total, aes(x = tg_c, y = jmax25, color = pathway, linetype = factor(cao))) +
  geom_line(size = 2) +
  scale_color_manual(values = c("blue", "red")) +  # Specify colors
  scale_linetype_manual(values = c("dotted", "solid")) +  # Specify linetypes
  labs(x = expression("Growth temperature ("*degree*"C)"),
       y = expression(bold(italic("J")["max25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       color = "Photo. pathway",
       linetype = expression("CO"[2])) +
  theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  theme(legend.position = "top")

# Chi plot
chi.plot <- ggplot(data = df.total, aes(x = tg_c, y = chi, color = pathway, 
                            linetype = factor(cao))) +
  geom_line(size = 2) +
  scale_color_manual(values = c("blue", "red")) +  # Specify colors
  scale_linetype_manual(values = c("dotted", "solid")) +  # Specify linetypes
  labs(x = expression("Growth temperature ("*degree*"C)"),
       y = expression(chi),
       color = "Photo. pathway",
       linetype = expression("CO"[2])) +
  theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  theme(legend.position = "top") 

# Vpmax plot
vpmax.plot <- ggplot(data = subset(df.total,pathway=="C4"), aes(x = tg_c, y = vpmax, color = pathway, 
                            linetype = factor(cao))) +
  geom_line(size = 2) +
  scale_color_manual(values = c("red")) +  # Specify colors
  scale_linetype_manual(values = c("dotted", "solid")) +  # Specify linetypes
  labs(x = expression("Growth temperature ("*degree*"C)"),
       y = expression(bold(italic("V")["pmax25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       color = "Photo. pathway",
       linetype = expression("CO"[2])) +
  theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  theme(legend.position = "top")

# Write photosynthesis figure
png("~/Desktop/Lemontree/zinny_photosynthesis.png", width = 8, height = 5, 
    units = "in", res = 600)
ac.plot
dev.off()

# Write mechanism figure
png("~/Desktop/Lemontree/zinny_mechanism.png", width = 20, height = 7, 
    units = "in", res = 600)
ggarrange(vcmax.plot, jmax.plot, chi.plot, vpmax.plot, 
          nrow = 1, ncol = 4, common.legend = TRUE, legend = "top") 
dev.off()


## calculate difference between c3 and c4 photosynthesis under different conditions

diff_photo_lowT_lowCO2 <- ((df.c4.420$Al[1] - df.c3.420$Al[1])/df.c4.420$Al[1]) * 100 #c4 percent difference at 20C
diff_photo_lowT_highCO2 <- ((df.c4.1000$Al[1] - df.c3.1000$Al[1])/df.c4.1000$Al[1]) * 100 #c4 percent difference at 20C
diff_photo_highT_lowCO2 <- ((df.c4.420$Al[16] - df.c3.420$Al[16])/df.c4.420$Al[16]) * 100 #c4 percent difference at 35C
diff_photo_highT_highCO2 <- ((df.c4.1000$Al[16] - df.c3.1000$Al[16])/df.c4.1000$Al[16]) * 100 #c4 percent difference at 35C




