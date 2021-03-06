library(ggplot2)
library(gridExtra)

scaleFUN <- function(x) formatC(x, format="e", digits=4)

# destination = soma
destination <- read.csv("meas_meas2r_ca_cyt", sep="\t")
# source = neurite tip
source <- read.csv("meas_meas2_ca_cyt", sep="\t")

colnames(source) <- c("Time", "Concentration")
colnames(destination) <- c("Time", "Concentration")
src <- ggplot(data=source, aes(x=Time, y=Concentration))
src <- src + geom_line()
src <- src + geom_point()
src <- src + ggtitle(bquote('Calcium wave experiment (' ~ bold('neurite tip') ~ ' -> ' ~ bold('soma') ~ ') with single spike'))
src <- src + xlab("Time [s]")
src <- src + ylab("Calcium concentration [M]")
src <- src  + ylab(bquote(Ca^'2+' ~~ "concentration [M]"))
src <- src + labs(subtitle="Neurite tip")
src <- src + theme_bw()
src <- src + theme(plot.subtitle=element_text(face="bold"))
src <- src + scale_y_continuous(labels=scaleFUN)
src <- src + scale_x_log10(labels=scaleFUN)
src <- src + geom_vline(xintercept = 100, linetype="solid", color = "green", size=0.5)
# src <- src + scale_y_continuous(labels=formatC(seq(min(source$Concentration,destination$Concentration), 
# max(source$Concentration,destination$Concentration), by=0.0000001), format="e", digits=2), 
#breaks=seq(min(source$Concentration), max(source$Concentration), by=0.0000001))

dest <- ggplot(data=destination, aes(x=Time, y=Concentration))
dest <- dest + geom_line()
dest <- dest + geom_point()
dest <- dest + xlab("Time [s]")
dest <- dest + ylab(bquote(Ca^'2+' ~~ "concentration [M]"))
dest <- dest + labs(subtitle="Soma")
dest <- dest + theme_bw()
dest <- dest + theme(plot.subtitle=element_text(face="bold"))
dest <- dest + scale_y_continuous(labels=scaleFUN)
dest <- dest + scale_x_log10(labels=scaleFUN)
dest <- dest + geom_vline(xintercept = 100, linetype="solid", color = "green", size=0.5)
#dest <- dest + scale_y_continuous(labels=formatC(seq(min(source$Concentration,destination$Concentration), 
#max(source$Concentration,destination$Concentration), by=0.0000001), format="e", digits=2), 
#breaks=seq(min(source$Concentration), max(source$Concentration), by=0.0000001))

ggsave(filename="ca_cyt_soma_vs_source.png", arrangeGrob(src, dest))
