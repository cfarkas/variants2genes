library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(dplyr)

Control<-read.table("Control.bam", skip="##", fill = TRUE, row.names=NULL)
Control<-subset(Control, select = c(1, 2))
Control<-plyr::rename(Control, c("row.names"="V1", "chrI"="V2"))
head(Control)
dim(Control)
Case<-read.table("Case.bam", skip="##", fill = TRUE, row.names=NULL)
Case<-subset(Case, select = c(1, 2))
Case<-plyr::rename(Case, c("row.names"="V1", "chrI"="V2"))
head(Case)
dim(Case)
########## subsetting chromosomes ###############

Control_chrI<-subset(Control, V1 == "chrI")
Control_chrII<-subset(Control, V1 == "chrII")
Control_chrIII<-subset(Control, V1 == "chrIII")
Control_chrIV<-subset(Control, V1 == "chrIV")
Control_chrV<-subset(Control, V1 == "chrV")
Control_chrVI<-subset(Control, V1 == "chrVI")
Control_chrVII<-subset(Control, V1 == "chrVII")
Control_chrVIII<-subset(Control, V1 == "chrVIII")
Control_chrIX<-subset(Control, V1 == "chrIX")
Control_chrX<-subset(Control, V1 == "chrX")
Control_chrXI<-subset(Control, V1 == "chrXI")
Control_chrXII<-subset(Control, V1 == "chrXII")
Control_chrXIII<-subset(Control, V1 == "chrXIII")
Control_chrXIV<-subset(Control, V1 == "chrXIV")
Control_chrXV<-subset(Control, V1 == "chrXV")
Control_chrXVI<-subset(Control, V1 == "chrXVI")
Control_chrM<-subset(Control, V1 == "chrM")

Case_chrI<-subset(Case, V1 == "chrI")
Case_chrII<-subset(Case, V1 == "chrII")
Case_chrIII<-subset(Case, V1 == "chrIII")
Case_chrIV<-subset(Case, V1 == "chrIV")
Case_chrV<-subset(Case, V1 == "chrV")
Case_chrVI<-subset(Case, V1 == "chrVI")
Case_chrVII<-subset(Case, V1 == "chrVII")
Case_chrVIII<-subset(Case, V1 == "chrVIII")
Case_chrIX<-subset(Case, V1 == "chrIX")
Case_chrX<-subset(Case, V1 == "chrX")
Case_chrXI<-subset(Case, V1 == "chrXI")
Case_chrXII<-subset(Case, V1 == "chrXII")
Case_chrXIII<-subset(Case, V1 == "chrXIII")
Case_chrXIV<-subset(Case, V1 == "chrXIV")
Case_chrXV<-subset(Case, V1 == "chrXV")
Case_chrXVI<-subset(Case, V1 == "chrXVI")
Case_chrM<-subset(Case, V1 == "chrM")

########## subsetting coordinates ###############

Control_chrI<-subset(Control_chrI, select = c(2))
Control_chrII<-subset(Control_chrII, select = c(2))
Control_chrIII<-subset(Control_chrIII, select = c(2))
Control_chrIV<-subset(Control_chrIV, select = c(2))
Control_chrV<-subset(Control_chrV, select = c(2))
Control_chrVI<-subset(Control_chrVI, select = c(2))
Control_chrVII<-subset(Control_chrVII, select = c(2))
Control_chrVIII<-subset(Control_chrVIII, select = c(2))
Control_chrIX<-subset(Control_chrIX, select = c(2))
Control_chrX<-subset(Control_chrX, select = c(2))
Control_chrXI<-subset(Control_chrXI, select = c(2))
Control_chrXII<-subset(Control_chrXII, select = c(2))
Control_chrXIII<-subset(Control_chrXIII, select = c(2))
Control_chrXIV<-subset(Control_chrXIV, select = c(2))
Control_chrXV<-subset(Control_chrXV, select = c(2))
Control_chrXVI<-subset(Control_chrXVI, select = c(2))
Control_chrM<-subset(Control_chrM, select = c(2))

Case_chrI<-subset(Case_chrI, select = c(2))
Case_chrII<-subset(Case_chrII, select = c(2))
Case_chrIII<-subset(Case_chrIII, select = c(2))
Case_chrIV<-subset(Case_chrIV, select = c(2))
Case_chrV<-subset(Case_chrV, select = c(2))
Case_chrVI<-subset(Case_chrVI, select = c(2))
Case_chrVII<-subset(Case_chrVII, select = c(2))
Case_chrVIII<-subset(Case_chrVIII, select = c(2))
Case_chrIX<-subset(Case_chrIX, select = c(2))
Case_chrX<-subset(Case_chrX, select = c(2))
Case_chrXI<-subset(Case_chrXI, select = c(2))
Case_chrXII<-subset(Case_chrXII, select = c(2))
Case_chrXIII<-subset(Case_chrXIII, select = c(2))
Case_chrXIV<-subset(Case_chrXIV, select = c(2))
Case_chrXV<-subset(Case_chrXV, select = c(2))
Case_chrXVI<-subset(Case_chrXVI, select = c(2))
Case_chrM<-subset(Case_chrM, select = c(2))

par(mfrow=c(4,10))

Control_chrI$Condition<-"Control"
Case_chrI$Condition<-"Case"
ChrI <- rbind(Control_chrI, Case_chrI)
p1 <- ggplot(ChrI, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrI") + ylab("reads") +labs(title="ChrI") + theme(legend.position="none")

Control_chrII$Condition <- "Control"
Case_chrII$Condition <- "Case"
ChrII <- rbind(Control_chrII, Case_chrII)
p2 <- ggplot(ChrII, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrII") + ylab("reads") +labs(title="ChrII") + theme(legend.position="none")

Control_chrIII$Condition <- "Control"
Case_chrIII$Condition <- "Case"
ChrIII <- rbind(Control_chrIII, Case_chrIII)
p3 <- ggplot(ChrIII, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrIII") + ylab("reads") +labs(title="ChrIII") + theme(legend.position="none")

Control_chrIV$Condition<-"Control"
Case_chrIV$Condition<-"Case"
ChrIV <- rbind(Control_chrIV, Case_chrIV)
p4 <- ggplot(ChrIV, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrIV") + ylab("reads") +labs(title="ChrIV") + theme(legend.position="none")

Control_chrV$Condition<-"Control"
Case_chrV$Condition<-"Case"
ChrV <- rbind(Control_chrV, Case_chrV)
p5 <- ggplot(ChrV, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrV") + ylab("reads") +labs(title="ChrV") + theme(legend.position="none")

Control_chrVI$Condition<-"Control"
Case_chrVI$Condition<-"Case"
ChrVI <- rbind(Control_chrVI, Case_chrVI)
p6 <- ggplot(ChrVI, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrVI") + ylab("reads") +labs(title="ChrVI") + theme(legend.position="none")

Control_chrVII$Condition<-"Control"
Case_chrVII$Condition<-"Case"
ChrVII <- rbind(Control_chrVII, Case_chrVII)
p7 <- ggplot(ChrVII, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrVII") + ylab("reads") +labs(title="ChrVII") + theme(legend.position="none")

Control_chrVIII$Condition<-"Control"
Case_chrVIII$Condition<-"Case"
ChrVIII <- rbind(Control_chrVIII, Case_chrVIII)
p8 <- ggplot(ChrVIII, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrVIII") + ylab("reads") +labs(title="ChrVIII") + theme(legend.position="none")

Control_chrIX$Condition<-"Control"
Case_chrIX$Condition<-"Case"
ChrIX <- rbind(Control_chrIX, Case_chrIX)
p9 <- ggplot(ChrIX, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrIX") + ylab("reads") +labs(title="ChrIX") + theme(legend.position="none")

Control_chrX$Condition<-"Control"
Case_chrX$Condition<-"Case"
ChrX <- rbind(Control_chrX, Case_chrX)
p10 <- ggplot(ChrX, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrX") + ylab("reads") +labs(title="ChrX") + theme(legend.position="none") 

Control_chrXI$Condition<-"Control"
Case_chrXI$Condition<-"Case"
ChrXI <- rbind(Control_chrXI, Case_chrXI)
p11 <- ggplot(ChrXI, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrXI") + ylab("reads") +labs(title="ChrXI") + theme(legend.position="none")

Control_chrXII$Condition<-"Control"
Case_chrXII$Condition<-"Case"
ChrXII <- rbind(Control_chrXII, Case_chrXII)
p12 <- ggplot(ChrXII, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrXII") + ylab("reads") +labs(title="ChrXII") + theme(legend.position="none")

Control_chrXIII$Condition<-"Control"
Case_chrXIII$Condition<-"Case"
ChrXIII <- rbind(Control_chrXIII, Case_chrXIII)
p13 <- ggplot(ChrXIII, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrXIII") + ylab("reads") +labs(title="ChrXIII") + theme(legend.position="none")

Control_chrXIV$Condition<-"Control"
Case_chrXIV$Condition<-"Case"
ChrXIV <- rbind(Control_chrXIV, Case_chrXIV)
p14 <- ggplot(ChrXIV, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrXIV") + ylab("reads") +labs(title="ChrXIV") + theme(legend.position="none")

Control_chrXV$Condition<-"Control"
Case_chrXV$Condition<-"Case"
ChrXV <- rbind(Control_chrXV, Case_chrXV)
p15 <- ggplot(ChrXV, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrXV") + ylab("reads") +labs(title="ChrXV") + theme(legend.position="none")

Control_chrXVI$Condition<-"Control"
Case_chrXVI$Condition<-"Case"
ChrXVI <- rbind(Control_chrXVI, Case_chrXVI)
p16 <- ggplot(ChrXVI, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrXVI") + ylab("reads") +labs(title="ChrXVI") + theme(legend.position="none")

Control_chrM$Condition<-"Control"
Case_chrM$Condition<-"Case"
ChrM <- rbind(Control_chrM, Case_chrM)
p17 <- ggplot(ChrM, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrM") + ylab("reads") +labs(title="ChrM") + theme(legend.position="none")

max=max(max(layer_scales(p1)$y$range$range), max(layer_scales(p2)$y$range$range), max(layer_scales(p3)$y$range$range), max(layer_scales(p4)$y$range$range), max(layer_scales(p5)$y$range$range), max(layer_scales(p6)$y$range$range), 
max(layer_scales(p7)$y$range$range), max(layer_scales(p8)$y$range$range), max(layer_scales(p9)$y$range$range), max(layer_scales(p10)$y$range$range), max(layer_scales(p11)$y$range$range), max(layer_scales(p12)$y$range$range), 
max(layer_scales(p13)$y$range$range), max(layer_scales(p14)$y$range$range), max(layer_scales(p15)$y$range$range), max(layer_scales(p16)$y$range$range), max(layer_scales(p17)$y$range$range))

p1 <- ggplot(ChrI, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrI") + ylab("reads") +labs(title="ChrI") + theme(legend.position="none") + ylim(0, max)

p2 <- ggplot(ChrII, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrII") + ylab("reads") +labs(title="ChrII") + theme(legend.position="none") + ylim(0, max)

p3 <- ggplot(ChrIII, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrIII") + ylab("reads") +labs(title="ChrIII") + theme(legend.position="none") + ylim(0, max)

p4 <- ggplot(ChrIV, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrIV") + ylab("reads") +labs(title="ChrIV") + theme(legend.position="none") + ylim(0, max)

p5 <- ggplot(ChrV, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrV") + ylab("reads") +labs(title="ChrV") + theme(legend.position="none") + ylim(0, max)

p6 <- ggplot(ChrVI, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrVI") + ylab("reads") +labs(title="ChrVI") + theme(legend.position="none") + ylim(0, max)

p7 <- ggplot(ChrVII, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrVII") + ylab("reads") +labs(title="ChrVII") + theme(legend.position="none") + ylim(0, max)

p8 <- ggplot(ChrVIII, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrVIII") + ylab("reads") +labs(title="ChrVIII") + theme(legend.position="none") + ylim(0, max)

p9 <- ggplot(ChrIX, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrIX") + ylab("reads") +labs(title="ChrIX") + theme(legend.position="none") + ylim(0, max)

p10 <- ggplot(ChrX, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrX") + ylab("reads") +labs(title="ChrX") + theme(legend.position="none") + ylim(0, max)

p11 <- ggplot(ChrXI, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrXI") + ylab("reads") +labs(title="ChrXI") + theme(legend.position="none") + ylim(0, max)

p12 <- ggplot(ChrXII, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrXII") + ylab("reads") +labs(title="ChrXII") + theme(legend.position="none") + ylim(0, max)

p13 <- ggplot(ChrXIII, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrXIII") + ylab("reads") +labs(title="ChrXIII") + theme(legend.position="none") + ylim(0, max)

p14 <- ggplot(ChrXIV, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrXIV") + ylab("reads") +labs(title="ChrXIV") + theme(legend.position="none") + ylim(0, max)

p15 <- ggplot(ChrXV, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrXV") + ylab("reads") +labs(title="ChrXV") + theme(legend.position="none") + ylim(0, max)

p16 <- ggplot(ChrXVI, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrXVI") + ylab("reads") +labs(title="ChrXVI") + theme(legend.position="none") + ylim(0, max)

p17 <- ggplot(ChrM, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrM") + ylab("reads") +labs(title="ChrM") + theme(legend.position="none") + ylim(0, max)

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lHeight <- sum(legend$Height)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           heights = unit.c(unit(1, "npc") - lHeight, lHeight)))

  grid.newpage()
  grid.draw(combined)

  # return gtable invisibly
  invisible(combined)

}

dsamp <- diamonds[sample(nrow(diamonds), 1000), ]

g <- grid_arrange_shared_legend(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, ncol = 4, nrow=5)
ggsave("sacCer3.pdf", g, width=50, height=52, units="cm")

dev.off()
proc.time()
sessionInfo()
