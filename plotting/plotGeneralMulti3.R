## refinement example (compare refinements at specific measurements)

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
    stop("Number of measurements does not match. For each cytosol measurement we need a corresponding endoplasmic reticulum measurement!")
  } 
  stop("Error: Required number of arguments is at least 5: LOG_SCALE END_TIME OUTPUT_NAME CYT_MEA ER_MEA")
}

a = (length(args)-3)/2
cyt <- data.frame(Time=c(), Concentration=c(), Refinement=c())
j=0
for (i in 4:(4+a-1)) {
  data <- read.csv(args[i], sep="\t", nrows=1750)
  colnames(data) <- c("Time", "Concentration")
  refinement <- rep(factor(j), length(data))
  data <- data.frame(Time=data$Time, Concentration=data$Concentration, Refinement=refinement)
  cyt <- rbind(cyt, data)
  j <- j + 1
}

pdist <- function(df) {
  do.call(rbind, combn(unique(df$Refinement), 2, function(x) {
    df1 <- subset(df, Refinement == x[1])
    df2 <- subset(df, Refinement == x[2])
    df3 <- merge(df1, df2, by = 'Time')
    Concentration <- abs(df3$Concentration.x - df3$Concentration.y)
    data.frame(combn = paste(x, collapse = ','), 
             time = df3$Time[which.max(Concentration)],
             newPos = max(df3$Concentration.y),
             newPos2 = max(df3$Concentration.x),
             max_difference = max(Concentration))
  }, simplify = FALSE))
}

pdist_annotate <- function(df, plot, hjust) {
  annotate("text", x = as.double(df[2]), y = as.double(df[3]), label = paste("||", df[1], "||_2 = ", scaleFUN(as.double(df[5])), sep=" "), hjust=hjust, size=2) 
}

pdist_lines <- function(df, plot, hjust) {
  geom_vline(xintercept = as.double(df[2]), linetype="dashed", color = "black", size=0.25)
}

er <- data.frame(Time=c(), Concentration=c(), Refinement=c())
j = 0 
for (i in (4+a):(length(args))) {
  data <- read.csv(args[i], sep="\t", nrows=1750)
 # print(args[i])
  colnames(data) <- c("Time", "Concentration")
  refinement <- rep(factor(j), length(data))
  data <- data.frame(Time=data$Time, Concentration=data$Concentration, Refinement=refinement)
  er <- rbind(er, data)
  j <- j + 1
}

## cytosol
cytosol <- ggplot(data=cyt, aes(x=Time, y=Concentration, group=Refinement, colour=Refinement))
cytosol <- cytosol + geom_line()
#cytosol <- cytosol + geom_point()
cytosol <- cytosol + ggtitle(bquote('Calcium wave experiment (' ~ bold('ER') ~ ' vs ' ~ bold('Cyt') ~ ') with single spike before 1st BP'))
cytosol <- cytosol + xlab("Time [s]")
cytosol <- cytosol + ylab("Calcium concentration [µM]")
cytosol <- cytosol + ylab(bquote(Ca^'2+' ~~ "concentration [M]"))
cytosol <- cytosol + labs(subtitle="Cytosol (At neurite tip)")
cytosol <- cytosol + theme_bw()
cytosol <- cytosol + theme(plot.subtitle=element_text(face="bold"))
cytosol <- cytosol + scale_y_continuous(labels=scaleFUN)
if (logScale) {
  cytosol <- cytosol + scale_x_log10(labels=scaleFUN)
}
cytosol <- cytosol + geom_vline(xintercept = endTime, linetype="solid", color = "green", size=0.5)
maxTime = endTime
if (endTime == -1) {
  maxTime = cyt[which.max(cyt$Concentration),]$Time
}
cytosol <- cytosol + coord_cartesian(xlim = c(0, maxTime))
cytosol <- cytosol + annotate("text", x = maxTime, y = max(cyt$Concentration), label = scaleFUN(max(cyt$Concentration)), hjust = 0.75)
cytosol <- cytosol + scale_color_hue(labels=c("0 ref: 2,304", "1 ref: 18,432", "2 ref: 147,456" , "3 ref: 1,179,648" , "4 ref: 9,437,184", "5 ref: 75,497,472"))

## ER
endoplasmic_reticulum <- ggplot(data=er, aes(x=Time, y=Concentration, group=Refinement, colour=Refinement))
endoplasmic_reticulum <- endoplasmic_reticulum + geom_line()
#endoplasmic_reticulum <- endoplasmic_reticulum + geom_point()
endoplasmic_reticulum <- endoplasmic_reticulum + ggtitle(bquote('Calcium wave experiment (' ~ bold('ER') ~ ' vs ' ~ bold('Cyt') ~ ') with single spike before 1st BP'))
endoplasmic_reticulum <- endoplasmic_reticulum + xlab("Time [s]")
endoplasmic_reticulum <- endoplasmic_reticulum + ylab("Calcium concentration [µM]")
endoplasmic_reticulum <- endoplasmic_reticulum + ylab(bquote(Ca^'2+' ~~ "concentration [M]"))
endoplasmic_reticulum <- endoplasmic_reticulum + labs(subtitle="endoplasmic_reticulum (At neurite tip)")
endoplasmic_reticulum <- endoplasmic_reticulum + theme_bw()
endoplasmic_reticulum <- endoplasmic_reticulum + theme(plot.subtitle=element_text(face="bold"))
endoplasmic_reticulum <- endoplasmic_reticulum + scale_y_continuous(labels=scaleFUN)
if (logScale) {
  endoplasmic_reticulum <- endoplasmic_reticulum + scale_x_log10(labels=scaleFUN)
}
endoplasmic_reticulum <- endoplasmic_reticulum + geom_vline(xintercept = endTime, linetype="solid", color = "green", size=0.5)
maxTime = endTime
if (endTime == -1) {
  maxTime = er[which.max(er$Concentration),]$Time
}
endoplasmic_reticulum <- endoplasmic_reticulum + coord_cartesian(xlim = c(0, maxTime))
endoplasmic_reticulum <- endoplasmic_reticulum + annotate("text", x = maxTime, y = max(er$Concentration), label = scaleFUN(max(er$Concentration)), hjust = 0.75)

temp = pdist(cyt)
print(temp)
# Use this to annotate distances between refienemtns 1,2, 2,3, 3,4 and so forth
#finalTable <- data.frame()
# for (i in 1:numRefs+1) {
#  from = i
#  to = i + 1
#  tmp = temp[temp[grep(paste(from, to, sep=","), temp$combn),]
#  finalTable <- rbind(finalTable, tmp)
#}
annos = apply(temp[grep("0,", temp$combn),], 1, pdist_annotate, cytosol, 0)
mylines = apply(temp[grep("0,", temp$combn),], 1, pdist_lines, cytosol, 0)
for (i in 0:length(annos)) {
  cytosol <- cytosol + annos[i+1]
  cytosol <- cytosol + mylines[i+1]
}

temp = pdist(er)
annos = apply(temp[grep("0,", temp$combn),], 1, pdist_annotate, endoplasmic_reticulum, 0)
mylines = apply(temp[grep("0,", temp$combn),], 1, pdist_lines, endoplasmic_reticulum, 0)
for (i in 0:length(annos)) {
  endoplasmic_reticulum <- endoplasmic_reticulum + annos[i+1]
  endoplasmic_reticulum <- endoplasmic_reticulum + mylines[i+1]
}

endoplasmic_reticulum <- endoplasmic_reticulum + scale_color_hue(labels=c("0 ref: 2,304", "1 ref: 18,432", "2 ref: 147,456" , "3 ref: 1,179,648" , "4 ref: 9,437,184", "5 ref: 75,497,472"))

ggsave(filename=outputname, arrangeGrob(cytosol, endoplasmic_reticulum))
