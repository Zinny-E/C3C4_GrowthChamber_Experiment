# Load libraries
library(ggplot2)
library(dplyr)
library(ggpubr)

##Optimality model hypothesis plots

setwd("~/git_repo/optimal_vcmax_R")

# Load optimal vcmax R script and functions therein
source("calc_optimal_vcmax.R")
sourceDirectory("functions")

# Do optimality model simulations for C4 species across range in temps and 
# binned into two CO2 treatments
df.c4.420 <- calc_optimal_vcmax(pathway = "C4", tg_c = seq(10, 45, 1), cao = 420)
df.c4.1000 <- calc_optimal_vcmax(pathway = "C4", tg_c = seq(10, 45, 1), cao = 1000)

# Do optimality model simulations for C3 species across range in temps and 
# binned into two CO2 treatments
df.c3.420 <- calc_optimal_vcmax(pathway = "C3", tg_c = seq(10, 45, 1), cao = 420)
df.c3.1000 <- calc_optimal_vcmax(pathway = "C3", tg_c = seq(10, 45, 1), cao = 1000)

# Merge all simulations into central data frame
df.total <- df.c3.420 %>% full_join(df.c3.1000) %>% full_join(df.c4.420) %>%
  full_join(df.c4.1000)

#merge just the c4 table
df.total_c4 <- df.c4.420 %>% full_join(df.c4.1000) 

# Net photosynthesis plot (amax)
amax.plot <- ggplot(data = subset(df.total_c4,pathway=="C4"), 
                  aes(x = tg_c, y = amax/amax[16], color = pathway, linetype = factor(cao))) +
  geom_line(size = 2) + coord_cartesian(ylim = c(0.5,2.0)) +
  scale_color_manual(values = c("red", "blue")) +  # Specify colors
  scale_linetype_manual(values = c("dotted", "solid")) +  # Specify linetypes
  scale_x_continuous(breaks = c(10, 15,20, 25, 30, 35, 40)) + # Specify breaks for specific temperature values # Specify linetypes
  labs(x = expression("Growth temperature ("*degree*"C)"),
       y = expression(bold(italic("A")["max25"])),
       color = "Photo. pathway",
       linetype = expression("CO"[2])) +
  theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  theme(legend.position = "top") 


# Vcmax plot
vcmax.plot <- ggplot(data = subset(df.total,pathway=="C3"), 
       aes(x = tg_c, y = vcmax25/vcmax25[16], color = pathway, linetype = factor(cao))) +
  geom_line(size = 2) + coord_cartesian(ylim = c(0.5,2.0)) +
  scale_color_manual(values = c("blue", "red")) +  # Specify colors
  scale_linetype_manual(values = c("dotted", "solid")) +  # Specify linetypes
  scale_x_continuous(breaks = c(10, 15,20, 25, 30, 35, 40)) + # Specify breaks for specific temperature values # Specify linetypes
  labs(x = expression("Growth temperature ("*degree*"C)"),
       y = expression(bold(italic("V")["cmax25"])),
       color = "Photo. pathway",
       linetype = expression("CO"[2])) +
  theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  theme(legend.position = "top")

# Jmax plot
jmax.plot <- ggplot(data = subset(df.total,pathway=="C3"), aes(x = tg_c, y = jmax25/jmax25[16], color = pathway, linetype = factor(cao))) +
  geom_line(size = 2) + coord_cartesian(ylim = c(0.5,2.0)) +
  scale_color_manual(values = c("blue", "red")) +  # Specify colors
  scale_linetype_manual(values = c("dotted", "solid")) + # Specify linetypes
  scale_x_continuous(breaks = c(10, 15,20, 25, 30, 35, 40)) + # Specify breaks for specific temperature values # Specify linetypes
  labs(x = expression("Growth temperature ("*degree*"C)"),
       y =  expression(bold(italic(" Investment in J")["max25"])),
       color = "Photo. pathway",
       linetype = expression("CO"[2])) +
  theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  theme(legend.position = "top")

# Chi plot
chi.plot.c3 <- ggplot(data = subset(df.total,pathway=="C3"), aes(x = tg_c, y = chi/chi[16], color = pathway, 
                            linetype = factor(cao))) +
  geom_line(size = 2) +
  scale_color_manual(values = c("blue", "red")) +  # Specify colors
  scale_linetype_manual(values = c("dotted", "solid")) +  # Specify linetypes
  scale_x_continuous(breaks = c(10, 15,20, 25, 30, 35, 40)) + # Specify breaks for specific temperature values # Specify linetypes
  labs(x = expression("Growth temperature ("*degree*"C)"),
       y = expression(chi),
       color = "Photo. pathway",
       linetype = expression("CO"[2])) +
  theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  theme(legend.position = "top") 

# Chi plot
chi.plot.c4 <- ggplot(data = subset(df.total_c4,pathway=="C4"),  aes(x = tg_c, y = chi/chi[16], color = pathway, 
                                        linetype = factor(cao))) +
  geom_line(size = 2) +
  scale_color_manual(values = c("red", "blue")) +  # Specify colors
  scale_linetype_manual(values = c("dotted", "solid")) +  # Specify linetypes
  scale_x_continuous(breaks = c(10, 15,20, 25, 30, 35, 40)) + # Specify breaks for specific temperature values # Specify linetypes
  labs(x = expression("Growth temperature ("*degree*"C)"),
       y = expression(chi),
       color = "Photo. pathway",
       linetype = expression("CO"[2])) +
  theme_bw(base_size = 20) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  theme(legend.position = "top") 

# Vpmax plot
vpmax.plot <- ggplot(data = subset(df.total_c4,pathway=="C4"), aes(x = tg_c, y = vpmax25/vpmax25[16], color = pathway, 
                            linetype = factor(cao))) +
  geom_line(size = 2) + coord_cartesian(ylim = c(0.5,2.0)) +
  scale_color_manual(values = c("red","blue")) +  # Specify colors
  scale_linetype_manual(values = c("dotted", "solid")) +  # Specify linetypes
  scale_x_continuous(breaks = c(10, 15,20, 25, 30, 35, 40)) + # Specify breaks for specific temperature values # Specify linetypes
  labs(x = expression("Growth temperature ("*degree*"C)"),
       y = expression(bold(italic("V")["pmax25"])),
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
png("~/Desktop/Lemontree_presentations/dissertation_hypothesis.png", width = 20, height = 7, 
    units = "in", res = 600)
ggarrange(vcmax.plot,  jmax.plot, chi.plot.c3,  vpmax.plot, amax.plot, chi.plot.c4,
          nrow = 2, ncol = 3, common.legend = TRUE, legend = "right") 
dev.off()


## calculate difference between c3 and c4 photosynthesis under different conditions

diff_photo_lowT_lowCO2 <- ((df.c4.420$Al[1] - df.c3.420$Al[1])/df.c4.420$Al[1]) * 100 #c4 percent difference at 20C
diff_photo_lowT_highCO2 <- ((df.c4.1000$Al[1] - df.c3.1000$Al[1])/df.c4.1000$Al[1]) * 100 #c4 percent difference at 20C
diff_photo_highT_lowCO2 <- ((df.c4.420$Al[16] - df.c3.420$Al[16])/df.c4.420$Al[16]) * 100 #c4 percent difference at 35C
diff_photo_highT_highCO2 <- ((df.c4.1000$Al[16] - df.c3.1000$Al[16])/df.c4.1000$Al[16]) * 100 #c4 percent difference at 35C





 ggplot(data = subset(df.total,pathway=="C4"), aes(x = vpmax25 , y = amax, color = pathway, 
                                                                linetype = factor(cao))) +
  geom_line(size = 2)


 

 
 
 df.c4.20 <- calc_optimal_vcmax(pathway = "C4", cao = seq(250, 1000, 50), tg_c = 20)
 df.c4.35 <- calc_optimal_vcmax(pathway = "C4", cao = seq(250, 1000, 50), tg_c = 35)
 df.c4.27.5 <- calc_optimal_vcmax(pathway = "C4", cao = seq(250, 1000, 50), tg_c = 27.5)
 
 # Merge all simulations into central data frame
 df.total.temp <- df.c4.20 %>% full_join(df.c4.35) %>% full_join(df.c4.27.5) 
 
 
 ggplot(data = df.total.temp, aes(x = tg_c, y = vpmax25, color = as.factor(tg_c))) +
   geom_boxplot()
 