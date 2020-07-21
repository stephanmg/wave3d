library(ggplot2)
library(gridExtra)

scaleFUN <- function(x) formatC(x, format="e", digits=2)

source <- read.csv("ttttttt_meas2_ca_cyt", sep="\t")
destination <- read.csv("ttttttt_meas2r_ca_cyt", sep="\t")
colnames(source) <- c("Time", "Concentration")
colnames(destination) <- c("Time", "Concentration")
src <- ggplot(data=source, aes(x=Time, y=Concentration))
src <- src + geom_line()
src <- src + geom_point()
src <- src + ggtitle(bquote('Calcium wave experiment (' ~ bold('neurite tip') ~ ' -> ' ~ bold('soma') ~ ')'))
src <- src + xlab("Time [s]")
src <- src + ylab("Calcium concentration [µM]")
src <- src + labs(subtitle="Neurite tip")
src <- src + theme_bw()
src <- src + theme(plot.subtitle=element_text(face="bold"))
src <- src + scale_y_continuous(labels=scaleFUN)
# src <- src + scale_y_continuous(labels=formatC(seq(min(source$Concentration,destination$Concentration), 
# max(source$Concentration,destination$Concentration), by=0.0000001), format="e", digits=2), 
#breaks=seq(min(source$Concentration), max(source$Concentration), by=0.0000001))

dest <- ggplot(data=destination, aes(x=Time, y=Concentration))
dest <- dest + geom_line()
dest <- dest + geom_point()
dest <- dest + xlab("Time [s]")
dest <- dest + ylab("Calcium concentration [µM]")
dest <- dest + labs(subtitle="Soma")
dest <- dest + theme_bw()
dest <- dest + theme(plot.subtitle=element_text(face="bold"))
dest <- dest + scale_y_continuous(labels=scaleFUN)
#dest <- dest + scale_y_continuous(labels=formatC(seq(min(source$Concentration,destination$Concentration), 
#max(source$Concentration,destination$Concentration), by=0.0000001), format="e", digits=2), 
#breaks=seq(min(source$Concentration), max(source$Concentration), by=0.0000001))

ggsave(filename="meas2_ca_cyt.png", arrangeGrob(src, dest))
