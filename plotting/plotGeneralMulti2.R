## Measurement example (compare measurements at specific measurements)

library(ggplot2)
library(gridExtra)

scaleFUN <- function(x) formatC(x, format="e", digits=4)

args = commandArgs(trailingOnly=TRUE)
logScale = F
endTime = -1
outputname = "test.png"
if (length(args) > 4) {
  logScale = args[1]
  endTime = as.double(args[2])
  outputname = args[3]
} else {
  if (length(args)-4 %% 2 != 0) {
    stop("Number of measurements does not match. For each measurement subset we need ER and CYT")
  } 
  stop("Error: Required number of arguments is at least 5: LOG_SCALE END_TIME OUTPUT_NAME CYT_MEA ER_MEA")
}

a = (length(args)-3)/2
cyt <- data.frame(Time=c(), Concentration=c(), Measurement=c())
j=0
for (i in 4:(4+a-1)) {
  data <- read.csv(args[i], sep="\t")
  colnames(data) <- c("Time", "Concentration")
  measurement <- rep(factor(j), length(data))
  data <- data.frame(Time=data$Time, Concentration=data$Concentration, Measurement=measurement)
  cyt <- rbind(cyt, data)
  j <- j + 1
}


er <- data.frame(Time=c(), Concentration=c(), Measurement=c())
j = 0 
for (i in (4+a):(length(args))) {
  data <- read.csv(args[i], sep="\t")
 # print(args[i])
  colnames(data) <- c("Time", "Concentration")
  measurement <- rep(factor(j), length(data))
  data <- data.frame(Time=data$Time, Concentration=data$Concentration, Measurement=measurement)
  er <- rbind(er, data)
  j <- j  + 1
}

## cytosol
cytosol <- ggplot(data=cyt, aes(x=Time, y=Concentration, group=Measurement, colour=Measurement))
cytosol <- cytosol + geom_line()
#cytosol <- cytosol + geom_point()
cytosol <- cytosol + ggtitle(bquote('Calcium wave experiment (' ~ bold('ER') ~ ' vs ' ~ bold('Cyt') ~ ') with single spike before 1st BP'))
cytosol <- cytosol + xlab("Time [s]")
cytosol <- cytosol + ylab("Calcium concentration [µM]")
cytosol <- cytosol + ylab(bquote(Ca^'2+' ~~ "concentration [M]"))
cytosol <- cytosol + labs(subtitle="Measurement subsets in Cytosol (Tip -> Soma)")
cytosol <- cytosol + theme_bw()
cytosol <- cytosol + theme(plot.subtitle=element_text(face="bold"))
cytosol <- cytosol + scale_y_continuous(labels=scaleFUN)
if (logScale) {
  cytosol <- cytosol + scale_x_log10(labels=scaleFUN)
}
cytosol <- cytosol + geom_vline(xintercept = endTime, linetype="solid", color = "green", size=0.5)
maxTime = cyt[which.max(cyt$Concentration),]$Time
cytosol <- cytosol + annotate("text", x = maxTime, y = max(cyt$Concentration), label = scaleFUN(max(cyt$Concentration)), hjust = -0.25)


## ER
endoplasmic_reticulum <- ggplot(data=er, aes(x=Time, y=Concentration, group=Measurement, colour=Measurement))
endoplasmic_reticulum <- endoplasmic_reticulum + geom_line()
#endoplasmic_reticulum <- endoplasmic_reticulum + geom_point()
endoplasmic_reticulum <- endoplasmic_reticulum + ggtitle(bquote('Calcium wave experiment (' ~ bold('ER') ~ ' vs ' ~ bold('Cyt') ~ ') with single spike before 1st BP'))
endoplasmic_reticulum <- endoplasmic_reticulum + xlab("Time [s]")
endoplasmic_reticulum <- endoplasmic_reticulum + ylab("Calcium concentration [µM]")
endoplasmic_reticulum <- endoplasmic_reticulum + ylab(bquote(Ca^'2+' ~~ "concentration [M]"))
endoplasmic_reticulum <- endoplasmic_reticulum + labs(subtitle="Measurement subsets in ER (Tip -> Soma)")
endoplasmic_reticulum <- endoplasmic_reticulum + theme_bw()
endoplasmic_reticulum <- endoplasmic_reticulum + theme(plot.subtitle=element_text(face="bold"))
endoplasmic_reticulum <- endoplasmic_reticulum + scale_y_continuous(labels=scaleFUN)
if (logScale) {
  endoplasmic_reticulum <- endoplasmic_reticulum + scale_x_log10(labels=scaleFUN)
}
endoplasmic_reticulum <- endoplasmic_reticulum + geom_vline(xintercept = endTime, linetype="solid", color = "green", size=0.5)
maxTime = er[which.max(er$Concentration),]$Time
endoplasmic_reticulum <- endoplasmic_reticulum + annotate("text", x = maxTime, y = max(er$Concentration), label = scaleFUN(max(er$Concentration)), hjust = -0.25)


ggsave(filename="test.png", arrangeGrob(cytosol, endoplasmic_reticulum))
