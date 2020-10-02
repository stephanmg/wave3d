# libraries
library(ggplot2)
library(gridExtra)

# use log scale?
logScale <- T

# test if there is at least one argument: if not, return an error
if (length(args) != 4) {
  stop("Exactly three argument must be supplied: CYT (CSV file), ER (CSV file), LOGSCALE (True or False), OUTNAME (With extension)", call.=FALSE)
}

# input files
d = args[1]
s = args[2]
o = args[3]

# helper functions
scaleFUN <- function(x) formatC(x, format="e", digits=4)

# destination = soma
destination <- read.csv(d, sep="\t")
# source = neurite tip
source <- read.csv(s, sep="\t")

colnames(source) <- c("Time", "Concentration")
colnames(destination) <- c("Time", "Concentration")
src <- ggplot(data=source, aes(x=Time, y=Concentration))
src <- src + geom_line()
#src <- src + geom_point()
src <- src + ggtitle(bquote('Calcium wave experiment (' ~ bold('ER') ~ ' vs ' ~ bold('Cyt') ~ ') with single spike before 1st BP'))
src <- src + xlab("Time [s]")
src <- src + ylab(bquote(Ca^'2+' ~~ "concentration [M]"))
src <- src + labs(subtitle="Cytosol (After branching point: 20 µm)")
src <- src + theme_bw()
src <- src + theme(plot.subtitle=element_text(face="bold"))
src <- src + scale_y_continuous(labels=scaleFUN)
#src <- src + scale_x_log10(labels=scaleFUN)
src <- src + geom_vline(xintercept = 200, linetype="solid", color = "green", size=0.5)
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
dest <- dest + labs(subtitle="ER (After branching point)")
dest <- dest + theme_bw()
dest <- dest + theme(plot.subtitle=element_text(face="bold"))
dest <- dest + scale_y_continuous(labels=scaleFUN)
#dest <- dest + scale_x_log10(labels=scaleFUN)
dest <- dest + geom_vline(xintercept = 200, linetype="solid", color = "green", size=0.5)
#dest <- dest + scale_y_continuous(labels=formatC(seq(min(source$Concentration,destination$Concentration), 
#max(source$Concentration,destination$Concentration), by=0.0000001), format="e", digits=2), 
#breaks=seq(min(source$Concentration), max(source$Concentration), by=0.0000001))

ggsave(filename=o, arrangeGrob(src, dest))
