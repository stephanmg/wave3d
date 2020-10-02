library(ggplot2)
library(gridExtra)

# destination = soma
destination1 <- read.csv("meas_meas4_ca_cyt", sep="\t")
destination2 <- read.csv("meas_meas3_ca_cyt", sep="\t")
# source = neurite tip
source <- read.csv("meas_meas5_ca_cyt", sep="\t")
# path length
length1 <- 27.490043661338
length2 <- 18.491816652675
length1 <- 10
length2 <- 20

# name columns
colnames(source) <- c("Time", "Concentration")
colnames(destination1) <- c("Time", "Concentration")
colnames(destination2) <- c("Time", "Concentration")

# calculate wave speed
maxDestIndex1 <- which.max(destination1$Concentration)
maxDestIndex2 <- which.max(destination2$Concentration)
maxSourceIndex <- which.max(source$Concentration)
cat(paste("        Time passed: ", (source[maxSourceIndex, 1]-destination1[maxDestIndex1, 1]), " [s]\n"))
cat(paste("    Distance (left): " , length1, " [µm]\n"))
cat(paste("   Distance (right): " , length2, " [µm]\n"))
cat(paste("  Wave speed (left): ", length1/(source[maxSourceIndex, 1]-destination1[maxDestIndex1, 1]), " [µm/s]\n"))
cat(paste(" Wave speed (right): ", length2/(source[maxSourceIndex, 1]-destination2[maxDestIndex2, 1]), " [µm/s]\n"))
