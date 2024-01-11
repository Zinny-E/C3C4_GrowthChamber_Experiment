# Compiling script for CO2 response curves and dark respiration values
# Note that all paths in the first part of this script assume working directory 
# root is the path to this script in the "C3C4_GrowthChamber_Experiment" git
# repository. Path assigned via Apple operating system, may differ  on Windows 
# operating systems.

###############################################################################
## Load readLicorData package
###############################################################################
library(readLicorData)
library(dplyr)

###############################################################################
## Load co2 response curve data files for chambers. Adding column to indicate whether
## curves were done from what chamber . 
###############################################################################

# Albert
eco2_HT_ch4d1_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-02-1143_chamber4") %>%
  mutate(chamber = 4)
#  write.csv(eco2_HT_ch4d1_alb, "../licor_cleaned/chamber_4/eco2_HT_ch4D1_alb.csv", row.names = FALSE)


eco2_HT_ch4d2_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-03-1017_chamber4_day2") %>%
  mutate(chamber = 4)
# write.csv(eco2_HT_ch4d2_alb, "../licor_cleaned/chamber_4/eco2_HT_ch4D2_alb.csv", row.names = FALSE)

eco2_HT_ch4d3_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-06-1027_chamber4_day3") %>%
  mutate(chamber = 4)
# write.csv(eco2_HT_ch4d3_alb, "../licor_cleaned/chamber_4/eco2_HT_ch4D3_alb.csv", row.names = FALSE)

aco2_HT_ch3d4_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-08-0948_chamber3_day4") %>%
  mutate(chamber = 3)
#  write.csv(aco2_HT_ch3d4_alb, "../licor_cleaned/chamber_3/aco2_HT_ch3D4_alb.csv", row.names = FALSE)
 
aco2_HT_ch3d5_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-10-1009_chamber3_day5") %>%
   mutate(chamber = 3)
# write.csv(aco2_HT_ch3d5_alb, "../licor_cleaned/chamber_3/aco2_HT_ch3D5_alb.csv", row.names = FALSE)
 
aco2_HT_ch3d6_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-14-1040_chamber3_day6") %>%
   mutate(chamber = 3)
# write.csv(aco2_HT_ch3d6_alb, "../licor_cleaned/chamber_3/aco2_HT_ch3D6_alb.csv", row.names = FALSE)
 
eco2_LT_ch5d1_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-21-1207_chamber5_d1") %>%
  mutate(chamber = 5)
# write.csv(eco2_LT_ch5d1_alb, "../licor_cleaned/chamber_5/eco2_LT_ch5D1_alb.csv", row.names = FALSE)
 
 eco2_LT_ch5d2_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-22-1036_chamber5_d2") %>%
   mutate(chamber = 5)
# write.csv(eco2_LT_ch5d2_alb, "../licor_cleaned/chamber_5/eco2_LT_ch5D2_alb.csv", row.names = FALSE)

eco2_LT_ch5d3_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-23-1032_chamber5_d3") %>%
   mutate(chamber = 5)
# write.csv(eco2_LT_ch5d3_alb, "../licor_cleaned/chamber_5/eco2_LT_ch5D3_alb.csv", row.names = FALSE)

aco2_LT_ch2d1_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-24-1026_chamber2_d1") %>%
  mutate(chamber = 2)
#  write.csv(aco2_LT_ch2d1_alb, "../licor_cleaned/chamber_2/aco2_LT_ch2D1_alb.csv", row.names = FALSE)
 
aco2_LT_ch2d2_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-28-1240_chamber2_d2") %>%
  mutate(chamber = 2)
# write.csv(aco2_LT_ch2d2_alb, "../licor_cleaned/chamber_2/aco2_LT_ch2D2_alb.csv", row.names = FALSE)

 aco2_LT_ch2d3_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-29-1216_chamber2_d3") %>%
   mutate(chamber = 2)
# write.csv(aco2_LT_ch2d3_alb, "../licor_cleaned/chamber_2/aco2_LT_ch2D3_alb.csv", row.names = FALSE)

 aco2_LT_ch2d4_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-30-1059_chamber2_d4") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d4_alb, "../licor_cleaned/chamber_2/aco2_LT_ch2D4_alb.csv", row.names = FALSE)
 
 aco2_LT_ch2d5_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-31-1105_chamber2_d5") %>%
   mutate(chamber = 2)
# write.csv(aco2_LT_ch2d5_alb, "../licor_cleaned/chamber_2/aco2_LT_ch2D5_alb.csv", row.names = FALSE)
 
 aco2_LT_ch2d6_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-04-03-1131_chamber2_d6") %>%
   mutate(chamber = 2)
# write.csv(aco2_LT_ch2d6_alb, "../licor_cleaned/chamber_2/aco2_LT_ch2D6_alb.csv", row.names = FALSE)
 
 aco2_LT_ch2d7_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-04-04-1107_chamber2_d7") %>%
   mutate(chamber = 2)
# write.csv(aco2_LT_ch2d7_alb, "../licor_cleaned/chamber_2/aco2_LT_ch2D7_alb.csv", row.names = FALSE)
 
 aco2_LT_ch2d8_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-04-05-1356_chamber2_d8") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d8_alb, "../licor_cleaned/chamber_2/aco2_LT_ch2D8_alb.csv", row.names = FALSE)
 
 
# Stan
 eco2_HT_ch4d1__stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-01-1327_chamber4") %>%
   mutate(chamber = 4)
  # write.csv(eco2_HT_ch4d1__stan, "../licor_cleaned/chamber_4/eco2_HT_ch4D1__stan.csv", row.names = FALSE)

 eco2_HT_ch4d1_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-02-1141_chamber4") %>%
   mutate(chamber = 4)
 # write.csv(eco2_HT_ch4d1_stan, "../licor_cleaned/chamber_4/eco2_HT_ch4D1_stan.csv", row.names = FALSE)
 
 eco2_HT_ch4d2_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-03-1013_chamber4_day2") %>%
   mutate(chamber = 4)
 # write.csv(eco2_HT_ch4d2_stan, "../licor_cleaned/chamber_4/eco2_HT_ch4D2_stan.csv", row.names = FALSE)

 eco2_HT_ch4d3_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-06-1025_chamber4_day3") %>%
   mutate(chamber = 4)
#  write.csv(eco2_HT_ch4d3_stan, "../licor_cleaned/chamber_4/eco2_HT_ch4D3_stan.csv", row.names = FALSE)
 
 aco2_HT_ch3d4_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-08-0949_chamber3_day4") %>%
   mutate(chamber = 3)
 #  write.csv(aco2_HT_ch3d4_stan, "../licor_cleaned/chamber_3/aco2_HT_ch3D4_stan.csv", row.names = FALSE)
 
 aco2_HT_ch3d5_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-10-1007_chamber3_day5") %>%
   mutate(chamber = 3)
 # write.csv(aco2_HT_ch3d5_stan, "../licor_cleaned/chamber_3/aco2_HT_ch3D5_stan.csv", row.names = FALSE)
 
 aco2_HT_ch3d6_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-14-1038_chamber3_day6") %>%
   mutate(chamber = 3)
 # write.csv(aco2_HT_ch3d6_stan, "../licor_cleaned/chamber_3/aco2_HT_ch3D6_stan.csv", row.names = FALSE)
 
 eco2_LT_ch5d1_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-21-1207_chamber5_d1") %>%
   mutate(chamber = 5)
 # write.csv(eco2_LT_ch5d1_stan, "../licor_cleaned/chamber_5/eco2_LT_ch5D1_stan.csv", row.names = FALSE)
 
 eco2_LT_ch5d2_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-22-1036_chamber5_d2") %>%
   mutate(chamber = 5)
 # write.csv(eco2_LT_ch5d2_stan, "../licor_cleaned/chamber_5/eco2_LT_ch5D2_stan.csv", row.names = FALSE)
 
 eco2_LT_ch5d3_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-23-1031_chamber5_d3") %>%
   mutate(chamber = 5)
 # write.csv(eco2_LT_ch5d3_stan, "../licor_cleaned/chamber_5/eco2_LT_ch5D3_stan.csv", row.names = FALSE) 
 
 aco2_LT_ch2d1_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-24-1025_chamber2_d1") %>%
   mutate(chamber = 2)
 #  write.csv(aco2_LT_ch2d1_stan, "../licor_cleaned/chamber_2/aco2_LT_ch2D1_stan.csv", row.names = FALSE)
 
 aco2_LT_ch2d2_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-28-1253_chamber2_d2") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d2_stan, "../licor_cleaned/chamber_2/aco2_LT_ch2D2_stan.csv", row.names = FALSE)
 
 aco2_LT_ch2d3_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-29-1215_chamber2_d3") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d3_stan, "../licor_cleaned/chamber_2/aco2_LT_ch2D3_stan.csv", row.names = FALSE)
 
 aco2_LT_ch2d4_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-30-1058_chamber2_d4") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d4_stan, "../licor_cleaned/chamber_2/aco2_LT_ch2D4_stan.csv", row.names = FALSE)
 
 aco2_LT_ch2d5_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-31-1105_chamber2_d5") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d5_stan, "../licor_cleaned/chamber_2/aco2_LT_ch2D5_stan.csv", row.names = FALSE)
 
 aco2_LT_ch2d6_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-04-03-1126_chamber2_d6") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d6_stan, "../licor_cleaned/chamber_2/aco2_LT_ch2D6_stan.csv", row.names = FALSE)
 
 aco2_LT_ch2d7_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-04-04-1106_chamber2_d7") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d7_stan, "../licor_cleaned/chamber_2/aco2_LT_ch2D7_stan.csv", row.names = FALSE)
 
 aco2_LT_ch2d8_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-04-05-1355_chamber2_d8") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d8_stan, "../licor_cleaned/chamber_2/aco2_LT_ch2D8_stan.csv", row.names = FALSE)
 
 
# Ozzie
 eco2_HT_ch4d1_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-01-1335_chamber4") %>%
   mutate(chamber = 4)
 # write.csv(eco2_HT_ch4d1_ozzie, "../licor_cleaned/chamber_4/eco2_HT_ch4D1_ozzie.csv", row.names = FALSE)
 
 eco2_HT_ch4d1.1_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-02-1454_chamber4") %>%
   mutate(chamber = 4)
 # write.csv(eco2_HT_ch4d1.1_ozzie, "../licor_cleaned/chamber_4/eco2_HT_ch4D1.1_ozzie.csv", row.names = FALSE)
 
 eco2_HT_ch4d2_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-03-1021_chamber4_day2") %>%
   mutate(chamber = 4)
 # write.csv(eco2_HT_ch4d2_ozzie, "../licor_cleaned/chamber_4/eco2_HT_ch4D2_ozzie.csv", row.names = FALSE)

 eco2_HT_ch4d2.1_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-03-1022_chamber4_day2") %>%
   mutate(chamber = 4)
 # write.csv(eco2_HT_ch4d2.1_ozzie, "../licor_cleaned/chamber_4/eco2_HT_ch4D2.1_ozzie.csv", row.names = FALSE)
 
 eco2_HT_ch4d3_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-06-1031_chamber4_day3") %>%
   mutate(chamber = 4)
 # write.csv(eco2_HT_ch4d3_ozzie, "../licor_cleaned/chamber_4/eco2_HT_ch4D3_ozzie.csv", row.names = FALSE)
 
 aco2_HT_ch3d4_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-08-0954_chamber3_day4") %>%
   mutate(chamber = 3)
 # write.csv(aco2_HT_ch3d4_ozzie, "../licor_cleaned/chamber_3/aco2_HT_ch3D4_ozzie.csv", row.names = FALSE)
 
 aco2_HT_ch3d4.1_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-08-1128_chamber3_day4") %>%
   mutate(chamber = 3)
 # write.csv(aco2_HT_ch3d4.1_ozzie, "../licor_cleaned/chamber_3/aco2_HT_ch3D4.1_ozzie.csv", row.names = FALSE)
 
 aco2_HT_ch3d4.2_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-08-1326_chamber3_day4_") %>%
   mutate(chamber = 3)
 # write.csv(aco2_HT_ch3d4.2_ozzie, "../licor_cleaned/chamber_3/aco2_HT_ch3D4.2_ozzie.csv", row.names = FALSE)
 
 aco2_HT_ch3d5_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-10-1008_chamber3_day5") %>%
   mutate(chamber = 3)
 # write.csv(aco2_HT_ch3d5_ozzie, "../licor_cleaned/chamber_3/aco2_HT_ch3D5_ozzie.csv", row.names = FALSE)
 
 aco2_HT_ch3d5.1_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-10-1013_chamber3_day5") %>%
   mutate(chamber = 3)
 # write.csv(aco2_HT_ch3d5.1_ozzie, "../licor_cleaned/chamber_3/aco2_HT_ch3D5.1_ozzie.csv", row.names = FALSE)
 
 aco2_HT_ch3d6_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-14-1043_chamber3_day6") %>%
   mutate(chamber = 3)
 # write.csv(aco2_HT_ch3d6_ozzie, "../licor_cleaned/chamber_3/aco2_HT_ch3D6_ozzie.csv", row.names = FALSE)
 
 aco2_HT_ch3d6.1_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-14-1044_chamber3_day6") %>%
   mutate(chamber = 3)
 # write.csv(aco2_HT_ch3d6.1_ozzie, "../licor_cleaned/chamber_3/aco2_HT_ch3D6.1_ozzie.csv", row.names = FALSE)
 
 eco2_LT_ch5d1_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-21-1212_chamber5_d1") %>%
   mutate(chamber = 5)
 # write.csv(eco2_LT_ch5d1_ozzie, "../licor_cleaned/chamber_5/eco2_LT_ch5D1_ozzie.csv", row.names = FALSE)
 
 eco2_LT_ch5d2_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-22-1041_chamber5_d2") %>%
   mutate(chamber = 5)
 # write.csv(eco2_LT_ch5d2_ozzie, "../licor_cleaned/chamber_5/eco2_LT_ch5D2_ozzie.csv", row.names = FALSE)
 
 eco2_LT_ch5d3_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-23-1034_chamber5_d3") %>%
   mutate(chamber = 5)
 # write.csv(eco2_LT_ch5d3_ozzie, "../licor_cleaned/chamber_5/eco2_LT_ch5D3_ozzie.csv", row.names = FALSE) 
 
 aco2_LT_ch2d1_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-24-1027_chamber2_d1") %>%
   mutate(chamber = 2)
 #  write.csv(aco2_LT_ch2d1_ozzie, "../licor_cleaned/chamber_2/aco2_LT_ch2D1_ozzie.csv", row.names = FALSE)
 
 aco2_LT_ch2d2_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-28-1255_chamber2_d2") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d2_ozzie, "../licor_cleaned/chamber_2/aco2_LT_ch2D2_ozzie.csv", row.names = FALSE)
 
 aco2_LT_ch2d3_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-29-1218_chamber2_d3") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d3_ozzie, "../licor_cleaned/chamber_2/aco2_LT_ch2D3_ozzie.csv", row.names = FALSE)
 
 aco2_LT_ch2d4_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-30-1101_chamber2_d4") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d4_ozzie, "../licor_cleaned/chamber_2/aco2_LT_ch2D4_ozzie.csv", row.names = FALSE)
 
 aco2_LT_ch2d5_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-31-1109_chamber2_d5") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d5_ozzie, "../licor_cleaned/chamber_2/aco2_LT_ch2D5_ozzie.csv", row.names = FALSE)
 
 aco2_LT_ch2d7_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-04-04-1108_chamber2_d7") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d7_ozzie, "../licor_cleaned/chamber_2/aco2_LT_ch2D7_ozzie.csv", row.names = FALSE)
 
 aco2_LT_ch2d8_ozzie <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-04-05-1401_chamber2_d8") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d8_ozzie, "../licor_cleaned/chamber_2/aco2_LT_ch2D8_ozzie.csv", row.names = FALSE)
 
 
 #Gibson
 eco2_HT_ch4d1_gibson <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-02-1454_chamber4") %>%
   mutate(chamber = 4)
 # write.csv(eco2_HT_ch4d1_gibson, "../licor_cleaned/chamber_4/eco2_HT_ch4D1_gibson.csv", row.names = FALSE)
 
 eco2_HT_ch4d2_gibson <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-03-1021_chamber4_day2") %>%
   mutate(chamber = 4)
 # write.csv(eco2_HT_ch4d2_gibson, "../licor_cleaned/chamber_4/eco2_HT_ch4D2_gibson.csv", row.names = FALSE)
 
 eco2_HT_ch4d3_gibson <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-06-1031_chamber4_day3") %>%
   mutate(chamber = 4)
 # write.csv(eco2_HT_ch4d3_gibson, "../licor_cleaned/chamber_4/eco2_HT_ch4D3_gibson.csv", row.names = FALSE)

 aco2_HT_ch3d4_gibson <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-01-1335_chamber4") %>%
   mutate(chamber = 3)
 #  write.csv(aco2_HT_ch3d4_gibson, "../licor_cleaned/chamber_3/aco2_HT_ch3D4_gibson.csv", row.names = FALSE)
 
 aco2_HT_ch3d5_gibson <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-10-1013_chamber3_day5") %>%
   mutate(chamber = 3)
 # write.csv(aco2_HT_ch3d5_gibson, "../licor_cleaned/chamber_3/aco2_HT_ch3D5_gibson.csv", row.names = FALSE)

 aco2_HT_ch3d6_gibson <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-14-1044_chamber3_day6") %>%
   mutate(chamber = 3)
 # write.csv(aco2_HT_ch3d6_gibson, "../licor_cleaned/chamber_3/aco2_HT_ch3D6_gibson.csv", row.names = FALSE) 
 
 eco2_LT_ch5d1_gibson <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-21-1209_chamber5_d1") %>%
   mutate(chamber = 5)
 # write.csv(eco2_LT_ch5d1_gibson, "../licor_cleaned/chamber_5/eco2_LT_ch5D1_gibson.csv", row.names = FALSE)
 
 eco2_LT_ch5d2_gibson <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-22-1039_chamber5_d2") %>%
   mutate(chamber = 5)
 # write.csv(eco2_LT_ch5d2_gibson, "../licor_cleaned/chamber_5/eco2_LT_ch5D2_gibson.csv", row.names = FALSE)
 
 eco2_LT_ch5d3_gibson <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-23-1033_chamber5_d3") %>%
   mutate(chamber = 5)
 # write.csv(eco2_LT_ch5d3_gibson, "../licor_cleaned/chamber_5/eco2_LT_ch5D3_gibson.csv", row.names = FALSE) 
 
 aco2_LT_ch2d1_gibson <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-24-1028_chamber2_d1") %>%
   mutate(chamber = 2)
 #  write.csv(aco2_LT_ch2d1_gibson, "../licor_cleaned/chamber_2/aco2_LT_ch2D1_gibson.csv", row.names = FALSE)
 
 aco2_LT_ch2d4_gibson <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-30-1102_chamber2_d4") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d4_gibson, "../licor_cleaned/chamber_2/aco2_LT_ch2D4_gibson.csv", row.names = FALSE)
 
 aco2_LT_ch2d5_gibson <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-31-1108_chamber2_d5") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d5_gibson, "../licor_cleaned/chamber_2/aco2_LT_ch2D5_gibson.csv", row.names = FALSE)
 
 aco2_LT_ch2d8_gibson <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-04-05-1359_chamber2_d8") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d8_gibson, "../licor_cleaned/chamber_2/aco2_LT_ch2D8_gibson.csv", row.names = FALSE)
 
 
 #New
 eco2_HT_ch4d3_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-06-1106_chamber4_day3") %>%
   mutate(chamber = 4)
 # write.csv(eco2_HT_ch4d3_new, "../licor_cleaned/chamber_4/eco2_HT_ch4D3_new.csv", row.names = FALSE)
 
 aco2_HT_ch3d4_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-08-0952_chamber3_day4") %>%
   mutate(chamber = 3)
 #  write.csv(aco2_HT_ch3d4_new, "../licor_cleaned/chamber_3/aco2_HT_ch3D4_new.csv", row.names = FALSE)
 
 aco2_HT_ch3d5_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-10-1011_chamber3_day5") %>%
   mutate(chamber = 3)
 # write.csv(aco2_HT_ch3d5_new, "../licor_cleaned/chamber_3/aco2_HT_ch3D5_new.csv", row.names = FALSE)
 
 aco2_HT_ch3d6_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-14-1046_chamber3_day6") %>%
   mutate(chamber = 3)
 # write.csv(aco2_HT_ch3d6_new, "../licor_cleaned/chamber_3/aco2_HT_ch3D6_new.csv", row.names = FALSE) 
 
 eco2_LT_ch5d1_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-21-1213_chamber5_d1") %>%
   mutate(chamber = 5)
 # write.csv(eco2_LT_ch5d1_new, "../licor_cleaned/chamber_5/eco2_LT_ch5D1_new.csv", row.names = FALSE)
 
 eco2_LT_ch5d2_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-22-1042_chamber5_d2") %>%
   mutate(chamber = 5)
 # write.csv(eco2_LT_ch5d2_new, "../licor_cleaned/chamber_5/eco2_LT_ch5D2_new.csv", row.names = FALSE)
 
 eco2_LT_ch5d3_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-23-1035_chamber5_d3") %>%
   mutate(chamber = 5)
 # write.csv(eco2_LT_ch5d3_new, "../licor_cleaned/chamber_5/eco2_LT_ch5D3_new.csv", row.names = FALSE) 
 
 aco2_LT_ch2d1_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-24-1028_chamber2_d1") %>%
   mutate(chamber = 2)
 #  write.csv(aco2_LT_ch2d1_new, "../licor_cleaned/chamber_2/aco2_LT_ch2D1_new.csv", row.names = FALSE)
 
 aco2_LT_ch2d2_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-28-1346_chamber2_d2") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d2_new, "../licor_cleaned/chamber_2/aco2_LT_ch2D2_new.csv", row.names = FALSE)
 
 aco2_LT_ch2d3_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-29-1219_chamber2_d3") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d3_new, "../licor_cleaned/chamber_2/aco2_LT_ch2D3_new.csv", row.names = FALSE)
 
 aco2_LT_ch2d4_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-08-0952_chamber3_day4") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d4_new, "../licor_cleaned/chamber_2/aco2_LT_ch2D4_new.csv", row.names = FALSE)
 
 aco2_LT_ch2d6_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-04-03-1129_chamber2_d6") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d6_new, "../licor_cleaned/chamber_2/aco2_LT_ch2D6_new.csv", row.names = FALSE)
 
 aco2_LT_ch2d7_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-04-04-1118_chamber2_d7") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d7_new, "../licor_cleaned/chamber_2/aco2_LT_ch2D7_new.csv", row.names = FALSE)
 
 aco2_LT_ch2d8_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-04-05-1401_chamber2_d8") %>%
   mutate(chamber = 2)
 # write.csv(aco2_LT_ch2d8_new, "../licor_cleaned/chamber_2/aco2_LT_ch2D8_new.csv", row.names = FALSE)

###############################################################################
## Load dark respiration data files for machine. Adding column to indicate 
## whether curves were done on which chamber 
###############################################################################

# Albert
 aco2_rd_ch2d1_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-24-1509_chamber2_d1_dr") %>%
   mutate(chamber = 2)
# write.csv(aco2_rd_ch2d1_alb, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d1_alb.csv", row.names = FALSE)
 
aco2_rd_ch2d3_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-29-1527_chamber2_d3_dr") %>%
   mutate(chamber = 2)
# write.csv(aco2_rd_ch2d3_alb, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d3_alb.csv", row.names = FALSE)
 
aco2_rd_ch2d4_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-30-1609_chamber2_d4_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d4_alb, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d4_alb.csv", row.names = FALSE)

aco2_rd_ch2d5_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-31-1517_chamber2_d5_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d5_alb, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d5_alb.csv", row.names = FALSE)

aco2_rd_ch2d6_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-04-03-1511_chamber2_d6_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d6_alb, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d6_alb.csv", row.names = FALSE)

aco2_rd_ch2d7_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-04-04-1543_chamber2_d7_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d7_alb, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d7_alb.csv", row.names = FALSE)

aco2_rd_ch2d8_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-04-05-1631_chamber2_d8_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d8_alb, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d8_alb.csv", row.names = FALSE)

 eco2_rd_ch4d1_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-02-1453_chamber4_dr") %>%
  mutate(chamber = 4)
# write.csv(eco2_rd_ch4d1_alb, "../licor_cleaned/chamber_4/dark_resp/eco2_rd_ch4d1_alb.csv", row.names = FALSE)

 eco2_rd_ch4d2_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-03-1507_chamber4_dr_day2") %>%
   mutate(chamber = 4)
 # write.csv(eco2_rd_ch4d2_alb, "../licor_cleaned/chamber_4/dark_resp/eco2_rd_ch4d2_alb.csv", row.names = FALSE)

 eco2_rd_ch4d3_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-06-1519_chamber4_dr_day3") %>%
   mutate(chamber = 4)
 # write.csv(eco2_rd_ch4d3_alb, "../licor_cleaned/chamber_4/dark_resp/eco2_rd_ch4d3_alb.csv", row.names = FALSE)

 aco2_rd_ch3d4_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-08-1432_chamber3_dr_day4") %>%
   mutate(chamber = 3)
 # write.csv(aco2_rd_ch3d4_alb, "../licor_cleaned/chamber_3/dark_resp/aco2_rd_ch3d4_alb.csv", row.names = FALSE)

aco2_rd_ch3d5_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-10-1502_chamber3_dr_day5") %>%
   mutate(chamber = 3)
# write.csv(aco2_rd_ch3d5_alb, "../licor_cleaned/chamber_3/dark_resp/aco2_rd_ch3d5_alb.csv", row.names = FALSE) 
 
aco2_rd_ch3d5_2_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-10-1529_chamber3_dr_day5_2") %>%
  mutate(chamber = 3)
# write.csv(aco2_rd_ch3d5_2_alb, "../licor_cleaned/chamber_3/dark_resp/aco2_rd_ch3d5_2_alb.csv", row.names = FALSE) 

aco2_rd_ch3d6_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-14-1521_chamber3_dr_day6") %>%
  mutate(chamber = 3)
# write.csv(aco2_rd_ch3d6_alb, "../licor_cleaned/chamber_3/dark_resp/aco2_rd_ch3d6_alb.csv", row.names = FALSE)

eco2_rd_ch5d1_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-21-2038_chamber5_d1_dr") %>%
  mutate(chamber = 5)
# write.csv(eco2_rd_ch5d1_alb, "../licor_cleaned/chamber_5/dark_resp/eco2_rd_ch5d1_alb.csv", row.names = FALSE)

eco2_rd_ch5d2_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-22-1928_chamber5_d2_dr") %>%
  mutate(chamber = 5)
# write.csv(eco2_rd_ch5d2_alb, "../licor_cleaned/chamber_5/dark_resp/eco2_rd_ch5d2_alb.csv", row.names = FALSE)

eco2_rd_ch5d3_alb <- licorData(location = "../licor_raw/Zinny_GC_Albert/2023-03-23-2000_chamber5_d3_dr") %>%
  mutate(chamber = 5)
# write.csv(eco2_rd_ch5d3_alb, "../licor_cleaned/chamber_5/dark_resp/eco2_rd_ch5d3_alb.csv", row.names = FALSE)

# Stan
eco2_rd_ch4d1_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-02-1606_chamber4_dr") %>%
  mutate(chamber = 4)
# write.csv(eco2_rd_ch4d1_stan, "../licor_cleaned/chamber_4/dark_resp/eco2_rd_ch4d1_stan.csv", row.names = FALSE)

eco2_rd_ch4d1.1_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-02-1728_chamber4_dr2") %>%
  mutate(chamber = 4)
# write.csv(eco2_rd_ch4d1.1_stan, "../licor_cleaned/chamber_4/dark_resp/eco2_rd_ch4d1.1_stan.csv", row.names = FALSE)

eco2_rd_ch4d2_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-03-1507_chamber4_dr_day2") %>%
  mutate(chamber = 4)
# write.csv(eco2_rd_ch4d2_stan, "../licor_cleaned/chamber_4/dark_resp/eco2_rd_ch4d2_stan.csv", row.names = FALSE)

aco2_rd_ch3d5_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-10-1502_chamber3_dr_day5") %>%
  mutate(chamber = 3)
# write.csv(aco2_rd_ch3d5_stan, "../licor_cleaned/chamber_3/dark_resp/aco2_rd_ch3d5_stan.csv", row.names = FALSE) 

aco2_rd_ch3d6_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-14-1511_chamber3_dr_day6") %>%
  mutate(chamber = 3)
# write.csv(aco2_rd_ch3d6_stan, "../licor_cleaned/chamber_3/dark_resp/aco2_rd_ch3d6_stan.csv", row.names = FALSE) 

eco2_rd_ch5d1_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-21-2038_chamber5_d1_dr") %>%
  mutate(chamber = 5)
# write.csv(eco2_rd_ch5d1_stan, "../licor_cleaned/chamber_5/dark_resp/eco2_rd_ch5d1_stan.csv", row.names = FALSE)

eco2_rd_ch5d2_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-22-1926_chamber5_d2_dr") %>%
  mutate(chamber = 5)
# write.csv(eco2_rd_ch5d2_stan, "../licor_cleaned/chamber_5/dark_resp/eco2_rd_ch5d2_stan.csv", row.names = FALSE)

eco2_rd_ch5d3_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-23-2000_chamber5_d3_dr") %>%
  mutate(chamber = 5)
# write.csv(eco2_rd_ch5d3_stan, "../licor_cleaned/chamber_5/dark_resp/eco2_rd_ch5d3_stan.csv", row.names = FALSE)

aco2_rd_ch2d1_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-24-1509_chamber2_d1_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d1_stan, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d1_stan.csv", row.names = FALSE)

aco2_rd_ch2d2_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-28-1642_chamber2_d2_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d2_stan, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d2_stan.csv", row.names = FALSE)

aco2_rd_ch2d3_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-29-1527_chamber2_d3_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d3_stan, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d3_stan.csv", row.names = FALSE)

aco2_rd_ch2d4_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-30-1609_chamber2_d4_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d4_stan, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d4_stan.csv", row.names = FALSE)

aco2_rd_ch2d5_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-03-31-1516_chamber2_d5_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d5_stan, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d5_stan.csv", row.names = FALSE)

aco2_rd_ch2d6_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-04-03-1510_chamber2_d6_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d6_stan, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d6_stan.csv", row.names = FALSE)

aco2_rd_ch2d7_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-04-04-1542_chamber2_d7_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d7_stan, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d7_stan.csv", row.names = FALSE)

aco2_rd_ch2d8_stan <- licorData(location = "../licor_raw/Zinny_GC_Stan/2023-04-05-1629_chamber2_d8_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d8_stan, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d8_stan.csv", row.names = FALSE)


# Ozzie
eco2_rd_ch4d1_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-01-1622_chamber4_dr") %>%
  mutate(chamber = 4)
# write.csv(eco2_rd_ch4d1_ozz, "../licor_cleaned/chamber_4/dark_resp/eco2_rd_ch4d1_ozz.csv", row.names = FALSE)

eco2_rd_ch4d1_1_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-02-1546_chamber4_dr") %>%
  mutate(chamber = 4)
# write.csv(eco2_rd_ch4d1_1_ozz, "../licor_cleaned/chamber_4/dark_resp/eco2_rd_ch4d1_1_ozz.csv", row.names = FALSE)

eco2_rd_ch4d2_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-03-1507_chamber4_dr_day2") %>%
  mutate(chamber = 4)
# write.csv(eco2_rd_ch4d2_ozz, "../licor_cleaned/chamber_4/dark_resp/eco2_rd_ch4d2_ozz.csv", row.names = FALSE)

eco2_rd_ch4d2_2_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-03-1509_chamber4_dr_day2") %>%
  mutate(chamber = 4)
# write.csv(eco2_rd_ch4d2_2_ozz, "../licor_cleaned/chamber_4/dark_resp/eco2_rd_ch4d2_2_ozz.csv", row.names = FALSE)


eco2_rd_ch4d3_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-06-1551_chamber4_dr_day3") %>%
  mutate(chamber = 4)
# write.csv(eco2_rd_ch4d3_ozz, "../licor_cleaned/chamber_4/dark_resp/eco2_rd_ch4d3_ozz.csv", row.names = FALSE)


aco2_rd_ch3d4_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-08-1606_chamber3_dr_day4") %>%
  mutate(chamber = 3)
# write.csv(aco2_rd_ch3d4_ozz, "../licor_cleaned/chamber_3/dark_resp/aco2_rd_ch3d4_ozz.csv", row.names = FALSE) 

aco2_rd_ch3d5_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-10-1504_chamber3_dr_day5") %>%
  mutate(chamber = 3)
# write.csv(aco2_rd_ch3d5_ozz, "../licor_cleaned/chamber_3/dark_resp/aco2_rd_ch3d5_ozz.csv", row.names = FALSE) 

aco2_rd_ch3d5_2_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-10-1507_chamber3_dr_day5") %>%
  mutate(chamber = 3)
# write.csv(aco2_rd_ch3d5_2_ozz, "../licor_cleaned/chamber_3/dark_resp/aco2_rd_ch3d5_2_ozz.csv", row.names = FALSE) 

aco2_rd_ch3d6_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-14-1513_chamber3_dr_day6") %>%
  mutate(chamber = 3)
# write.csv(aco2_rd_ch3d6_ozz, "../licor_cleaned/chamber_3/dark_resp/aco2_rd_ch3d6_ozz.csv", row.names = FALSE)

aco2_rd_ch3d6.1_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-14-1514_chamber3_dr_day6") %>%
  mutate(chamber = 3)
# write.csv(aco2_rd_ch3d6.1_ozz, "../licor_cleaned/chamber_3/dark_resp/aco2_rd_ch3d6.1_ozz.csv", row.names = FALSE)

eco2_rd_ch5d1_1_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-21-2041_chamber5_d1_dr") %>%
  mutate(chamber = 5)
# write.csv(eco2_rd_ch5d1_1_ozz, "../licor_cleaned/chamber_5/dark_resp/eco2_rd_ch5d1_1_ozz.csv", row.names = FALSE)

eco2_rd_ch5d2_2_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-22-1929_chamber5_d2_dr") %>%
  mutate(chamber = 5)
# write.csv(eco2_rd_ch5d2_2_ozz, "../licor_cleaned/chamber_5/dark_resp/eco2_rd_ch5d2_2_ozz.csv", row.names = FALSE)

eco2_rd_ch5d3_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-23-2003_chamber5_d3_dr") %>%
  mutate(chamber = 5)
# write.csv(eco2_rd_ch5d3_ozz, "../licor_cleaned/chamber_5/dark_resp/eco2_rd_ch5d3_ozz.csv", row.names = FALSE)

eco2_rd_ch5d3_3_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-23-2041_chamber5_d3_dr") %>%
  mutate(chamber = 5)
# write.csv(eco2_rd_ch5d3_3_ozz, "../licor_cleaned/chamber_5/dark_resp/eco2_rd_ch5d3_3_ozz.csv", row.names = FALSE)

aco2_rd_ch2d1_1_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-24-1531_chamber2_d1_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d1_1_ozz, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d1_1_ozz.csv", row.names = FALSE)

aco2_rd_ch2d2_2_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-28-1647_chamber2_d2_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d2_2_ozz, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d2_2_ozz.csv", row.names = FALSE)

aco2_rd_ch2d3_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-29-1530_chamber2_d3_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d3_ozz, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d3_ozz.csv", row.names = FALSE)

aco2_rd_ch2d4_4_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-30-1614_chamber2_d4_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d4_4_ozz, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d4_4_ozz.csv", row.names = FALSE)

aco2_rd_ch2d5_5_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-03-31-1520_chamber2_d5_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d5_5_ozz, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d5_5_ozz.csv", row.names = FALSE)

aco2_rd_ch2d7_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-04-04-1545_chamber2_d7_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d7_ozz, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d7_ozz.csv", row.names = FALSE)

aco2_rd_ch2d8_8_ozz <- licorData(location = "../licor_raw/Zinny_GC_OzzieGibson/2023-04-05-1633_chamber2_d8_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d8_8_ozz, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d8_8_ozz.csv", row.names = FALSE)


#New
eco2_rd_ch4d1_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-06-1523_chamber4_dr_day3") %>%
  mutate(chamber = 4)
# write.csv(eco2_rd_ch4d1_new, "../licor_cleaned/chamber_4/dark_resp/eco2_rd_ch4d1_new.csv", row.names = FALSE)

aco2_rd_ch3d4_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-08-1436_chamber3_dr_day4") %>%
  mutate(chamber = 3)
# write.csv(aco2_rd_ch3d4_new, "../licor_cleaned/chamber_3/dark_resp/aco2_rd_ch3d4_new.csv", row.names = FALSE)

aco2_rd_ch3d5_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-10-1509_chamber3_dr_day5") %>%
  mutate(chamber = 3)
# write.csv(aco2_rd_ch3d5_new, "../licor_cleaned/chamber_3/dark_resp/aco2_rd_ch3d5_new.csv", row.names = FALSE)

aco2_rd_ch3d6_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-14-1527_chamber3_dr_day6") %>%
  mutate(chamber = 3)
# write.csv(aco2_rd_ch3d6_new, "../licor_cleaned/chamber_3/dark_resp/aco2_rd_ch3d6_new.csv", row.names = FALSE)

eco2_rd_ch5d1_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-21-2042_chamber5_d1_dr") %>%
  mutate(chamber = 5)
# write.csv(eco2_rd_ch5d1_new, "../licor_cleaned/chamber_5/dark_resp/eco2_rd_ch5d1_new.csv", row.names = FALSE)

eco2_rd_ch5d2_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-22-1930_chamber5_d2_dr") %>%
  mutate(chamber = 5)
# write.csv(eco2_rd_ch5d2_new, "../licor_cleaned/chamber_5/dark_resp/eco2_rd_ch5d2_new.csv", row.names = FALSE)

eco2_rd_ch5d3_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-23-2004_chamber5_d3_dr") %>%
  mutate(chamber = 5)
# write.csv(eco2_rd_ch5d3_new, "../licor_cleaned/chamber_5/dark_resp/eco2_rd_ch5d3_new.csv", row.names = FALSE)

aco2_rd_ch2d1_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-24-1532_chamber2_d1_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d1_new, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d1_new.csv", row.names = FALSE)

aco2_rd_ch2d2_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-28-1657_chamber2_d2_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d2_new, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d2_new.csv", row.names = FALSE)

aco2_rd_ch2d3_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-29-1532_chamber2_d3_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d3_new, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d3_new.csv", row.names = FALSE)

aco2_rd_ch2d4_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-30-1615_chamber2_d4_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d4_new, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d4_new.csv", row.names = FALSE)

aco2_rd_ch2d5_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-03-31-1521_chamber2_d5_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d5_new, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d5_new.csv", row.names = FALSE)

aco2_rd_ch2d6_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-04-03-1516_chamber2_d6_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d6_new, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d6_new.csv", row.names = FALSE)

aco2_rd_ch2d7_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-04-04-1547_chamber2_d7_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d7_new, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d7_new.csv", row.names = FALSE)

aco2_rd_ch2d8_new <- licorData(location = "../licor_raw/Zinny_GC_New/2023-04-05-1635_chamber2_d8_dr") %>%
  mutate(chamber = 2)
# write.csv(aco2_rd_ch2d8_new, "../licor_cleaned/chamber_2/dark_resp/aco2_rd_ch2d8_new.csv", row.names = FALSE)

###############################################################################
## Merge co2 response curves into single file. Useful for 'fitacis' when 
## fitting multiple curves
###############################################################################
# NOTE: Using list.files notation to avoid common merge conflict with 
# readLicorData package. Cols seem to be assigned different classes when
# cleaned through 'licorData', which makes merging files difficult/unnecessarily
# time consuming.for the list.file change chamber name
# Reloading files into list of data frames, then merging through
# reshape::merge_all() seems to do the trick.

# List files
file.list <- list.files("../licor_cleaned/chamber_4",
                        recursive = TRUE,
                        pattern = "\\.csv$",
                        full.names = TRUE)
file.list <- setNames(file.list, stringr::str_extract(basename(file.list), 
                                                      '.*(?=\\.csv)'))

# Merge list of data frames, arrange by machine, measurement type, id, and time elapsed
merged_curves <- lapply(file.list, read.csv) %>%
  reshape::merge_all() %>%
  arrange(machine, chamber, id, elapsed)

co2_resp_chamber_2 <- merged_curves %>%
  filter(chamber == 2)
write.csv(co2_resp_chamber_2, "../licor_cleaned/TxCO2_combined_datasheets/TXCO2_co2_resp_chamber_2.csv", row.names = FALSE)

co2_resp_chamber_3 <- merged_curves %>%
  filter(chamber == 3)
write.csv(co2_resp_chamber_3, "../licor_cleaned/TxCO2_combined_datasheets/TXCO2_co2_resp_chamber_3.csv", row.names = FALSE)

co2_resp_chamber_4 <- merged_curves %>%
  filter(chamber == 4)
write.csv(co2_resp_chamber_4, "../licor_cleaned/TxCO2_combined_datasheets/TXCO2_co2_resp_chamber_4.csv", row.names = FALSE)

co2_resp_chamber_5 <- merged_curves %>%
  filter(chamber == 5)
write.csv(co2_resp_chamber_5, "../licor_cleaned/datasheets/TXCO2_co2_resp_chamber_5.csv", row.names = FALSE)


###############################################################################
## Merge rd files into single file. Useful for 'fitacis' when fitting multiple
## curves
###############################################################################
# NOTE: Using list.files notation to avoid common merge conflict with 
# readLicorData package. Cols seem to be assigned different classes when
# cleaned through 'licorData', which makes merging files difficult/unnecessarily
# time consuming. Reloading files into list of data frames, then merging through
# reshape::merge_all() seems to do the trick.


# List files
file.list.rd <- list.files("../licor_cleaned/chamber_4/dark_resp",
                        recursive = TRUE,
                        pattern = "\\.csv$",
                        full.names = TRUE)
file.list.rd <- setNames(file.list.rd, stringr::str_extract(basename(file.list.rd), 
                                                      '.*(?=\\.csv)'))

# Merge list of data frames, arrange by machine, measurement type, id, and time elapsed
rd <- lapply(file.list.rd, read.csv) %>%
  reshape::merge_all()

rd.chamber_4 <- rd %>%
  filter(chamber == 4) %>%
  group_by(id, chamber) %>%
  mutate(A = ifelse(A > 0, NA, A),
         rd = abs(A)) %>%
  select(id,rd,Tleaf,machine,chamber) %>%
  arrange(id)

write.csv(rd.chamber_4, "../licor_cleaned/TxCO2_combined_datasheets/TxCO2_rd_chamber_4.csv", row.names = FALSE)

rd.chamber_2 <- rd %>%
  filter(chamber == 2) %>%
  group_by(id, chamber) %>%
  mutate(A = ifelse(A > 0, NA, A),
         rd = abs(A)) %>%
  select(id,rd,Tleaf,machine,chamber) %>%
  arrange(id)

write.csv(rd.chamber_2, "../licor_cleaned/datasheets/TxCO2_rd_chamber_2.csv", row.names = FALSE)

## End of data cleaning, ready for curve fitting ##