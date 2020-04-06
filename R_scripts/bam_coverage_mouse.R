library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(dplyr)

Control<-read.table("Control.vcf", skip="##", fill = TRUE, row.names=NULL)
Control<-subset(Control, select = c(1, 2))
Control<-plyr::rename(Control, c("row.names"="V1", "chr1"="V2"))
head(Control)
dim(Control)
Case<-read.table("Case.vcf", skip="##", fill = TRUE, row.names=NULL)
Case<-subset(Case, select = c(1, 2))
Case<-plyr::rename(Case, c("row.names"="V1", "chr1"="V2"))
head(Case)
dim(Case)
########## subsetting chromosomes ###############

Control_chr1<-subset(Control, V1 == "chr1")
Control_chr2<-subset(Control, V1 == "chr2")
Control_chr3<-subset(Control, V1 == "chr3")
Control_chr4<-subset(Control, V1 == "chr4")
Control_chr5<-subset(Control, V1 == "chr5")
Control_chr6<-subset(Control, V1 == "chr6")
Control_chr7<-subset(Control, V1 == "chr7")
Control_chr8<-subset(Control, V1 == "chr8")
Control_chr9<-subset(Control, V1 == "chr9")
Control_chr10<-subset(Control, V1 == "chr10")
Control_chr11<-subset(Control, V1 == "chr11")
Control_chr12<-subset(Control, V1 == "chr12")
Control_chr13<-subset(Control, V1 == "chr13")
Control_chr14<-subset(Control, V1 == "chr14")
Control_chr15<-subset(Control, V1 == "chr15")
Control_chr16<-subset(Control, V1 == "chr16")
Control_chr17<-subset(Control, V1 == "chr17")
Control_chr18<-subset(Control, V1 == "chr18")
Control_chr19<-subset(Control, V1 == "chr19")
Control_chrX<-subset(Control, V1 == "chrX")

Case_chr1<-subset(Case, V1 == "chr1")
Case_chr2<-subset(Case, V1 == "chr2")
Case_chr3<-subset(Case, V1 == "chr3")
Case_chr4<-subset(Case, V1 == "chr4")
Case_chr5<-subset(Case, V1 == "chr5")
Case_chr6<-subset(Case, V1 == "chr6")
Case_chr7<-subset(Case, V1 == "chr7")
Case_chr8<-subset(Case, V1 == "chr8")
Case_chr9<-subset(Case, V1 == "chr9")
Case_chr10<-subset(Case, V1 == "chr10")
Case_chr11<-subset(Case, V1 == "chr11")
Case_chr12<-subset(Case, V1 == "chr12")
Case_chr13<-subset(Case, V1 == "chr13")
Case_chr14<-subset(Case, V1 == "chr14")
Case_chr15<-subset(Case, V1 == "chr15")
Case_chr16<-subset(Case, V1 == "chr16")
Case_chr17<-subset(Case, V1 == "chr17")
Case_chr18<-subset(Case, V1 == "chr18")
Case_chr19<-subset(Case, V1 == "chr19")
Case_chrX<-subset(Case, V1 == "chrX")


########## subsetting coordinates ###############

Control_chr1<-subset(Control_chr1, select = c(2))
Control_chr2<-subset(Control_chr2, select = c(2))
Control_chr3<-subset(Control_chr3, select = c(2))
Control_chr4<-subset(Control_chr4, select = c(2))
Control_chr5<-subset(Control_chr5, select = c(2))
Control_chr6<-subset(Control_chr6, select = c(2))
Control_chr7<-subset(Control_chr7, select = c(2))
Control_chr8<-subset(Control_chr8, select = c(2))
Control_chr9<-subset(Control_chr9, select = c(2))
Control_chr10<-subset(Control_chr10, select = c(2))
Control_chr11<-subset(Control_chr11, select = c(2))
Control_chr12<-subset(Control_chr12, select = c(2))
Control_chr13<-subset(Control_chr13, select = c(2))
Control_chr14<-subset(Control_chr14, select = c(2))
Control_chr15<-subset(Control_chr15, select = c(2))
Control_chr16<-subset(Control_chr16, select = c(2))
Control_chr17<-subset(Control_chr17, select = c(2))
Control_chr18<-subset(Control_chr18, select = c(2))
Control_chr19<-subset(Control_chr19, select = c(2))
Control_chrX<-subset(Control_chrX, select = c(2))

Case_chr1<-subset(Case_chr1, select = c(2))
Case_chr2<-subset(Case_chr2, select = c(2))
Case_chr3<-subset(Case_chr3, select = c(2))
Case_chr4<-subset(Case_chr4, select = c(2))
Case_chr5<-subset(Case_chr5, select = c(2))
Case_chr6<-subset(Case_chr6, select = c(2))
Case_chr7<-subset(Case_chr7, select = c(2))
Case_chr8<-subset(Case_chr8, select = c(2))
Case_chr9<-subset(Case_chr9, select = c(2))
Case_chr10<-subset(Case_chr10, select = c(2))
Case_chr11<-subset(Case_chr11, select = c(2))
Case_chr12<-subset(Case_chr12, select = c(2))
Case_chr13<-subset(Case_chr13, select = c(2))
Case_chr14<-subset(Case_chr14, select = c(2))
Case_chr15<-subset(Case_chr15, select = c(2))
Case_chr16<-subset(Case_chr16, select = c(2))
Case_chr17<-subset(Case_chr17, select = c(2))
Case_chr18<-subset(Case_chr18, select = c(2))
Case_chr19<-subset(Case_chr19, select = c(2))
Case_chrX<-subset(Case_chrX, select = c(2))

par(mfrow=c(2,10))

Control_chr1$Condition<-"Control"
Case_chr1$Condition<-"Case"
Chr1 <- rbind(Control_chr1, Case_chr1)
p1 <- ggplot(Chr1, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr1") + ylab("variants") +labs(title="Chr1") + theme(legend.position="none")

Control_chr2$Condition <- "Control"
Case_chr2$Condition <- "Case"
Chr2 <- rbind(Control_chr2, Case_chr2)
p2 <- ggplot(Chr2, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr2") + ylab("variants") +labs(title="Chr2") + theme(legend.position="none")

Control_chr3$Condition <- "Control"
Case_chr3$Condition <- "Case"
Chr3 <- rbind(Control_chr3, Case_chr3)
p3 <- ggplot(Chr3, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr3") + ylab("variants") +labs(title="Chr3") + theme(legend.position="none")

Control_chr4$Condition<-"Control"
Case_chr4$Condition<-"Case"
Chr4 <- rbind(Control_chr4, Case_chr4)
p4 <- ggplot(Chr4, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr4") + ylab("variants") +labs(title="Chr4") + theme(legend.position="none")

Control_chr5$Condition<-"Control"
Case_chr5$Condition<-"Case"
Chr5 <- rbind(Control_chr5, Case_chr5)
p5 <- ggplot(Chr5, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr5") + ylab("variants") +labs(title="Chr5") + theme(legend.position="none")

Control_chr6$Condition<-"Control"
Case_chr6$Condition<-"Case"
Chr6 <- rbind(Control_chr6, Case_chr6)
p6 <- ggplot(Chr6, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr6") + ylab("variants") +labs(title="Chr6") + theme(legend.position="none")

Control_chr7$Condition<-"Control"
Case_chr7$Condition<-"Case"
Chr7 <- rbind(Control_chr7, Case_chr7)
p7 <- ggplot(Chr7, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr7") + ylab("variants") +labs(title="Chr7") + theme(legend.position="none")

Control_chr8$Condition<-"Control"
Case_chr8$Condition<-"Case"
Chr8 <- rbind(Control_chr8, Case_chr8)
p8 <- ggplot(Chr8, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr8") + ylab("variants") +labs(title="Chr8") + theme(legend.position="none")

Control_chr9$Condition<-"Control"
Case_chr9$Condition<-"Case"
Chr9 <- rbind(Control_chr9, Case_chr9)
p9 <- ggplot(Chr9, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr9") + ylab("variants") +labs(title="Chr9") + theme(legend.position="none")

Control_chr10$Condition<-"Control"
Case_chr10$Condition<-"Case"
Chr10 <- rbind(Control_chr10, Case_chr10)
p10 <- ggplot(Chr10, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr10") + ylab("variants") +labs(title="Chr10") + theme(legend.position="none") 

Control_chr11$Condition<-"Control"
Case_chr11$Condition<-"Case"
Chr11 <- rbind(Control_chr11, Case_chr11)
p11 <- ggplot(Chr11, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr11") + ylab("variants") +labs(title="Chr11") + theme(legend.position="none")

Control_chr12$Condition<-"Control"
Case_chr12$Condition<-"Case"
Chr12 <- rbind(Control_chr12, Case_chr12)
p12 <- ggplot(Chr12, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr12") + ylab("variants") +labs(title="Chr12") + theme(legend.position="none")

Control_chr13$Condition<-"Control"
Case_chr13$Condition<-"Case"
Chr13 <- rbind(Control_chr13, Case_chr13)
p13 <- ggplot(Chr13, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr13") + ylab("variants") +labs(title="Chr13") + theme(legend.position="none")

Control_chr14$Condition<-"Control"
Case_chr14$Condition<-"Case"
Chr14 <- rbind(Control_chr14, Case_chr14)
p14 <- ggplot(Chr14, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr14") + ylab("variants") +labs(title="Chr14") + theme(legend.position="none")

Control_chr15$Condition<-"Control"
Case_chr15$Condition<-"Case"
Chr15 <- rbind(Control_chr15, Case_chr15)
p15 <- ggplot(Chr15, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr15") + ylab("variants") +labs(title="Chr15") + theme(legend.position="none")

Control_chr16$Condition<-"Control"
Case_chr16$Condition<-"Case"
Chr16 <- rbind(Control_chr16, Case_chr16)
p16 <- ggplot(Chr16, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr16") + ylab("variants") +labs(title="Chr16") + theme(legend.position="none")

Control_chr17$Condition<-"Control"
Case_chr17$Condition<-"Case"
Chr17 <- rbind(Control_chr17, Case_chr17)
p17 <- ggplot(Chr17, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr17") + ylab("variants") +labs(title="Chr17") + theme(legend.position="none")

Control_chr18$Condition<-"Control"
Case_chr18$Condition<-"Case"
Chr18 <- rbind(Control_chr18, Case_chr18)
p18 <- ggplot(Chr18, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr18") + ylab("variants") +labs(title="Chr18") + theme(legend.position="none")

Control_chr19$Condition<-"Control"
Case_chr19$Condition<-"Case"
Chr19 <- rbind(Control_chr19, Case_chr19)
p19 <- ggplot(Chr19, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr19") + ylab("variants") +labs(title="Chr19") + theme(legend.position="none")

Control_chrX$Condition<-"Control"
Case_chrX$Condition<-"Case"
ChrX <- rbind(Control_chrX, Case_chrX)
pX <- ggplot(ChrX, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrX") + ylab("variants") +labs(title="ChrX") + theme(legend.position="none")


max=max(max(layer_scales(p1)$y$range$range), max(layer_scales(p2)$y$range$range), max(layer_scales(p3)$y$range$range), max(layer_scales(p4)$y$range$range), max(layer_scales(p5)$y$range$range), max(layer_scales(p6)$y$range$range), 
max(layer_scales(p7)$y$range$range), max(layer_scales(p8)$y$range$range), max(layer_scales(p9)$y$range$range), max(layer_scales(p10)$y$range$range), max(layer_scales(p11)$y$range$range), max(layer_scales(p12)$y$range$range), 
max(layer_scales(p13)$y$range$range), max(layer_scales(p14)$y$range$range),max(layer_scales(p15)$y$range$range), max(layer_scales(p16)$y$range$range), max(layer_scales(p17)$y$range$range), max(layer_scales(p18)$y$range$range), 
max(layer_scales(p19)$y$range$range), max(layer_scales(pX)$y$range$range))

p1 <- ggplot(Chr1, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr1") + ylab("variants") +labs(title="Chr1") + theme(legend.position="none") + ylim(0, max)

p2 <- ggplot(Chr2, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr2") + ylab("variants") +labs(title="Chr2") + theme(legend.position="none") + ylim(0, max)

p3 <- ggplot(Chr3, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr3") + ylab("variants") +labs(title="Chr3") + theme(legend.position="none") + ylim(0, max)

p4 <- ggplot(Chr4, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr4") + ylab("variants") +labs(title="Chr4") + theme(legend.position="none") + ylim(0, max)

p5 <- ggplot(Chr5, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr5") + ylab("variants") +labs(title="Chr5") + theme(legend.position="none") + ylim(0, max)

p6 <- ggplot(Chr6, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr6") + ylab("variants") +labs(title="Chr6") + theme(legend.position="none") + ylim(0, max)

p7 <- ggplot(Chr7, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr7") + ylab("variants") +labs(title="Chr7") + theme(legend.position="none") + ylim(0, max)

p8 <- ggplot(Chr8, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr8") + ylab("variants") +labs(title="Chr8") + theme(legend.position="none") + ylim(0, max)

p9 <- ggplot(Chr9, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr9") + ylab("variants") +labs(title="Chr9") + theme(legend.position="none") + ylim(0, max)

p10 <- ggplot(Chr10, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr10") + ylab("variants") +labs(title="Chr10") + theme(legend.position="none") + ylim(0, max)

p11 <- ggplot(Chr11, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr11") + ylab("variants") +labs(title="Chr11") + theme(legend.position="none") + ylim(0, max)

p12 <- ggplot(Chr12, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr12") + ylab("variants") +labs(title="Chr12") + theme(legend.position="none") + ylim(0, max)

p13 <- ggplot(Chr13, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr13") + ylab("variants") +labs(title="Chr13") + theme(legend.position="none") + ylim(0, max)

p14 <- ggplot(Chr14, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr14") + ylab("variants") +labs(title="Chr14") + theme(legend.position="none") + ylim(0, max)

p15 <- ggplot(Chr15, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr15") + ylab("variants") +labs(title="Chr15") + theme(legend.position="none") + ylim(0, max)

p16 <- ggplot(Chr16, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr16") + ylab("variants") +labs(title="Chr16") + theme(legend.position="none") + ylim(0, max)

p17 <- ggplot(Chr17, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr17") + ylab("variants") +labs(title="Chr17") + theme(legend.position="none") + ylim(0, max)

p18 <- ggplot(Chr18, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr18") + ylab("variants") +labs(title="Chr18") + theme(legend.position="none") + ylim(0, max)

p19 <- ggplot(Chr19, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr19") + ylab("variants") +labs(title="Chr19") + theme(legend.position="none") + ylim(0, max)

pX <- ggplot(ChrX, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrX") + ylab("variants") +labs(title="ChrX") + theme(legend.position="none") + ylim(0, max)


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

g <- grid_arrange_shared_legend(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, pX, ncol = 4, nrow=6)
ggsave("graph.pdf", g, width=50, height=52, units="cm")

dev.off()
proc.time()
sessionInfo()
