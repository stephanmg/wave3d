library(ggplot2)
library(gridExtra)

scaleFUN <- function(x) formatC(x, format="e", digits=4)

args = commandArgs(trailingOnly=TRUE)
logScale = F
endTime = -1
if (length(args)==2) {
  logScale = args[1]
  endTime = as.double(args[2])
}

# destination = soma
destination <- read.csv("meas_meas1ER_ca_er", sep="\t")
# source = neurite tip
source <- read.csv("meas_meas1_ca_cyt", sep="\t")

colnames(source) <- c("Time", "Concentration")
colnames(destination) <- c("Time", "Concentration")
src <- ggplot(data=source, aes(x=Time, y=Concentration))
src <- src + geom_line()
#src <- src + geom_point()
src <- src + ggtitle(bquote('Calcium wave experiment (' ~ bold('ER') ~ ' vs ' ~ bold('Cyt') ~ ') with single spike before 1st BP'))
src <- src + xlab("Time [s]")
src <- src + ylab(bquote(Ca^'2+' ~~ "concentration [M]"))
src <- src + labs(subtitle="Cytosol (Before branching point right)")
src <- src + theme_bw()
src <- src + theme(plot.subtitle=element_text(face="bold"))
src <- src + scale_y_continuous(labels=scaleFUN)
if (logScale) {
  src <- src + scale_x_log10(labels=scaleFUN)
}
src <- src + geom_vline(xintercept = endTime, linetype="solid", color = "green", size=0.5)
# src <- src + scale_y_continuous(labels=formatC(seq(min(source$Concentration,destination$Concentration), 
# max(source$Concentration,destination$Concentration), by=0.0000001), format="e", digits=2), 
#breaks=seq(min(source$Concentration), max(source$Concentration), by=0.0000001))
src <- src + geom_hline(yintercept=max(source$Concentration), linetype='dotted', col = 'red')
maxTime = source[which.max(source$Concentration),]$Time
src <- src + annotate("text", x = maxTime, y = max(source$Concentration), label = scaleFUN(max(source$Concentration)), hjust = -0.25)

dest <- ggplot(data=destination, aes(x=Time, y=Concentration))
dest <- dest + geom_line()
#dest <- dest + geom_point()
dest <- dest + xlab("Time [s]")
dest <- dest + ylab(bquote(Ca^'2+' ~~ "concentration [M]"))
dest <- dest + labs(subtitle="ER (Before branching point right)")
dest <- dest + theme_bw()
dest <- dest + theme(plot.subtitle=element_text(face="bold"))
dest <- dest + scale_y_continuous(labels=scaleFUN)
if (logScale) {
  dest <- dest + scale_x_log10(labels=scaleFUN)
}
dest <- dest + geom_vline(xintercept = endTime, linetype="solid", color = "green", size=0.5)
#dest <- dest + scale_y_continuous(labels=formatC(seq(min(source$Concentration,destination$Concentration), 
#max(source$Concentration,destination$Concentration), by=0.0000001), format="e", digits=2), 
#breaks=seq(min(source$Concentration), max(source$Concentration), by=0.0000001))

ggsave(filename="ca_cyt_vs_er_branch_right.png", arrangeGrob(src, dest))
