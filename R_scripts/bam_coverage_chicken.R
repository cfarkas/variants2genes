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
Control_chr20<-subset(Control, V1 == "chr20")
Control_chr21<-subset(Control, V1 == "chr21")
Control_chr22<-subset(Control, V1 == "chr22")
Control_chr23<-subset(Control, V1 == "chr23")
Control_chr24<-subset(Control, V1 == "chr24")
Control_chr25<-subset(Control, V1 == "chr25")
Control_chr26<-subset(Control, V1 == "chr26")
Control_chr27<-subset(Control, V1 == "chr27")
Control_chr28<-subset(Control, V1 == "chr28")
Control_chr30<-subset(Control, V1 == "chr30")
Control_chr31<-subset(Control, V1 == "chr31")
Control_chr32<-subset(Control, V1 == "chr32")
Control_chr33<-subset(Control, V1 == "chr33")
Control_chrZ<-subset(Control, V1 == "chrZ")

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
Case_chr20<-subset(Case, V1 == "chr20")
Case_chr21<-subset(Case, V1 == "chr21")
Case_chr22<-subset(Case, V1 == "chr22")
Case_chr23<-subset(Case, V1 == "chr23")
Case_chr24<-subset(Case, V1 == "chr24")
Case_chr25<-subset(Case, V1 == "chr25")
Case_chr26<-subset(Case, V1 == "chr26")
Case_chr27<-subset(Case, V1 == "chr27")
Case_chr28<-subset(Case, V1 == "chr28")
Case_chr30<-subset(Case, V1 == "chr30")
Case_chr31<-subset(Case, V1 == "chr31")
Case_chr32<-subset(Case, V1 == "chr32")
Case_chr33<-subset(Case, V1 == "chr33")
Case_chrZ<-subset(Case, V1 == "chrZ")

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
Control_chr20<-subset(Control_chr20, select = c(2))
Control_chr21<-subset(Control_chr21, select = c(2))
Control_chr22<-subset(Control_chr22, select = c(2))
Control_chr23<-subset(Control_chr23, select = c(2))
Control_chr24<-subset(Control_chr24, select = c(2))
Control_chr25<-subset(Control_chr25, select = c(2))
Control_chr26<-subset(Control_chr26, select = c(2))
Control_chr27<-subset(Control_chr27, select = c(2))
Control_chr28<-subset(Control_chr28, select = c(2))
Control_chr30<-subset(Control_chr30, select = c(2))
Control_chr31<-subset(Control_chr31, select = c(2))
Control_chr32<-subset(Control_chr32, select = c(2))
Control_chr33<-subset(Control_chr33, select = c(2))
Control_chrZ<-subset(Control_chrZ, select = c(2))

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
Case_chr20<-subset(Case_chr20, select = c(2))
Case_chr21<-subset(Case_chr21, select = c(2))
Case_chr22<-subset(Case_chr22, select = c(2))
Case_chr23<-subset(Case_chr23, select = c(2))
Case_chr24<-subset(Case_chr24, select = c(2))
Case_chr25<-subset(Case_chr25, select = c(2))
Case_chr26<-subset(Case_chr26, select = c(2))
Case_chr27<-subset(Case_chr27, select = c(2))
Case_chr28<-subset(Case_chr28, select = c(2))
Case_chr30<-subset(Case_chr30, select = c(2))
Case_chr31<-subset(Case_chr31, select = c(2))
Case_chr32<-subset(Case_chr32, select = c(2))
Case_chr33<-subset(Case_chr33, select = c(2))
Case_chrZ<-subset(Case_chrZ, select = c(2))

par(mfrow=c(4,10))

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

Control_chr20$Condition<-"Control"
Case_chr20$Condition<-"Case"
Chr20 <- rbind(Control_chr20, Case_chr20)
p20 <- ggplot(Chr20, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr20") + ylab("variants") +labs(title="Chr20") + theme(legend.position="none")

Control_chr21$Condition<-"Control"
Case_chr21$Condition<-"Case"
Chr21 <- rbind(Control_chr21, Case_chr21)
p21 <- ggplot(Chr21, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr21") + ylab("variants") +labs(title="Chr21") + theme(legend.position="none")

Control_chr22$Condition<-"Control"
Case_chr22$Condition<-"Case"
Chr22 <- rbind(Control_chr22, Case_chr22)
p22 <- ggplot(Chr22, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr22") + ylab("variants") +labs(title="Chr22") + theme(legend.position="none")

Control_chr23$Condition<-"Control"
Case_chr23$Condition<-"Case"
Chr23 <- rbind(Control_chr23, Case_chr23)
p23 <- ggplot(Chr23, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr23") + ylab("variants") +labs(title="Chr23") + theme(legend.position="none")

Control_chr24$Condition<-"Control"
Case_chr24$Condition<-"Case"
Chr24 <- rbind(Control_chr24, Case_chr24)
p24 <- ggplot(Chr24, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr24") + ylab("variants") +labs(title="Chr24") + theme(legend.position="none")

Control_chr25$Condition<-"Control"
Case_chr25$Condition<-"Case"
Chr25 <- rbind(Control_chr25, Case_chr25)
p25 <- ggplot(Chr25, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr25") + ylab("variants") +labs(title="Chr25") + theme(legend.position="none")

Control_chr26$Condition<-"Control"
Case_chr26$Condition<-"Case"
Chr26 <- rbind(Control_chr26, Case_chr26)
p26 <- ggplot(Chr26, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr26") + ylab("variants") +labs(title="Chr26") + theme(legend.position="none")

Control_chr27$Condition<-"Control"
Case_chr27$Condition<-"Case"
Chr27 <- rbind(Control_chr27, Case_chr27)
p27 <- ggplot(Chr27, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr27") + ylab("variants") +labs(title="Chr27") + theme(legend.position="none")

Control_chr28$Condition<-"Control"
Case_chr28$Condition<-"Case"
Chr28 <- rbind(Control_chr28, Case_chr28)
p28 <- ggplot(Chr28, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr28") + ylab("variants") +labs(title="Chr28") + theme(legend.position="none")

Control_chr30$Condition<-"Control"
Case_chr30$Condition<-"Case"
Chr30 <- rbind(Control_chr30, Case_chr30)
p30 <- ggplot(Chr30, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr30") + ylab("variants") +labs(title="Chr30") + theme(legend.position="none")

Control_chr31$Condition<-"Control"
Case_chr31$Condition<-"Case"
Chr31 <- rbind(Control_chr31, Case_chr31)
p31 <- ggplot(Chr31, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr31") + ylab("variants") +labs(title="Chr31") + theme(legend.position="none")

Control_chr32$Condition<-"Control"
Case_chr32$Condition<-"Case"
Chr32 <- rbind(Control_chr32, Case_chr32)
p32 <- ggplot(Chr32, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr32") + ylab("variants") +labs(title="Chr32") + theme(legend.position="none")

Control_chr33$Condition<-"Control"
Case_chr33$Condition<-"Case"
Chr33 <- rbind(Control_chr33, Case_chr33)
p33 <- ggplot(Chr33, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr33") + ylab("variants") +labs(title="Chr33") + theme(legend.position="none")

Control_chrZ$Condition<-"Control"
Case_chrZ$Condition<-"Case"
ChrZ <- rbind(Control_chrZ, Case_chrZ)
pZ <- ggplot(ChrZ, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrZ") + ylab("variants") +labs(title="ChrZ") + theme(legend.position="none")


max=max(max(layer_scales(p1)$y$range$range), max(layer_scales(p2)$y$range$range), max(layer_scales(p3)$y$range$range), max(layer_scales(p4)$y$range$range), max(layer_scales(p5)$y$range$range), max(layer_scales(p6)$y$range$range), 
max(layer_scales(p7)$y$range$range), max(layer_scales(p8)$y$range$range), max(layer_scales(p9)$y$range$range), max(layer_scales(p10)$y$range$range), max(layer_scales(p11)$y$range$range), max(layer_scales(p12)$y$range$range), 
max(layer_scales(p13)$y$range$range), max(layer_scales(p14)$y$range$range), max(layer_scales(p15)$y$range$range), max(layer_scales(p16)$y$range$range), max(layer_scales(p17)$y$range$range), max(layer_scales(p18)$y$range$range), 
max(layer_scales(p19)$y$range$range), max(layer_scales(p20)$y$range$range), max(layer_scales(p21)$y$range$range), max(layer_scales(p22)$y$range$range), max(layer_scales(p23)$y$range$range), max(layer_scales(p24)$y$range$range), max(layer_scales(p25)$y$range$range), max(layer_scales(p26)$y$range$range), max(layer_scales(p27)$y$range$range), max(layer_scales(p28)$y$range$range), max(layer_scales(p30)$y$range$range), max(layer_scales(p31)$y$range$range), max(layer_scales(p32)$y$range$range), max(layer_scales(p33)$y$range$range), max(layer_scales(pZ)$y$range$range))

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

p20 <- ggplot(Chr20, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr20") + ylab("variants") +labs(title="Chr20") + theme(legend.position="none") + ylim(0, max)

p21 <- ggplot(Chr21, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr21") + ylab("variants") +labs(title="Chr21") + theme(legend.position="none") + ylim(0, max)

p22 <- ggplot(Chr22, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr22") + ylab("variants") +labs(title="Chr22") + theme(legend.position="none") + ylim(0, max)

p23 <- ggplot(Chr23, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr23") + ylab("variants") +labs(title="Chr23") + theme(legend.position="none") + ylim(0, max)

p24 <- ggplot(Chr24, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr24") + ylab("variants") +labs(title="Chr24") + theme(legend.position="none") + ylim(0, max)

p25 <- ggplot(Chr25, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr25") + ylab("variants") +labs(title="Chr25") + theme(legend.position="none") + ylim(0, max)

p26 <- ggplot(Chr26, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr26") + ylab("variants") +labs(title="Chr26") + theme(legend.position="none") + ylim(0, max)

p27 <- ggplot(Chr27, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr27") + ylab("variants") +labs(title="Chr27") + theme(legend.position="none") + ylim(0, max)

p28 <- ggplot(Chr28, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr28") + ylab("variants") +labs(title="Chr28") + theme(legend.position="none") + ylim(0, max)

p30 <- ggplot(Chr30, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr30") + ylab("variants") +labs(title="Chr30") + theme(legend.position="none") + ylim(0, max)

p31 <- ggplot(Chr31, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr31") + ylab("variants") +labs(title="Chr31") + theme(legend.position="none") + ylim(0, max)

p32 <- ggplot(Chr32, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr32") + ylab("variants") +labs(title="Chr32") + theme(legend.position="none") + ylim(0, max)

p33 <- ggplot(Chr33, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in Chr33") + ylab("variants") +labs(title="Chr33") + theme(legend.position="none") + ylim(0, max)

pZ <- ggplot(ChrZ, aes(V2, fill = Condition)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity') + xlab("Position in ChrZ") + ylab("variants") +labs(title="ChrZ") + theme(legend.position="none") + ylim(0, max)


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

g <- grid_arrange_shared_legend(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p30, p31, p32, p33, pZ, ncol = 6, nrow=6)
ggsave("graph.pdf", g, width=50, height=52, units="cm")

dev.off()
proc.time()
sessionInfo()
