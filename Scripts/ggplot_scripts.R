library(ggplot2)
QD.plot <- ggplot(data = Hardfilter, aes(x=QD)) + geom_density(alpha=0.2)
QD.plot
#generating MQ  PLOT
MQ.plot <- ggplot(data = Hardfilter, aes(x=MQ)) + geom_density(alpha=0.2)
MQ.plot
#Generating SOR plot
SOR.plot <- ggplot(data = Hardfilter, aes(x=SOR)) + geom_density(alpha=0.2)
SOR.plot
#Generating FS
FS.plot <- ggplot(data = Hardfilter, aes(x=FS)) + geom_density(alpha=0.2)
FS.plot + scale_x_log10()
