library(ggplot2)
library(ggsci)
library(gridExtra)

## read data and construct data frames for panel A
A10<-as.matrix(read.csv("./simulation_results/max10_missed1.csv", header = FALSE, sep = ","), byrow = TRUE)
A20<-as.matrix(read.csv("./simulation_results/max20_missed1.csv", header = FALSE, sep = ","), byrow = TRUE)
A50<-as.matrix(read.csv("./simulation_results/max50_missed1.csv", header = FALSE, sep = ","), byrow = TRUE)
A<-rbind(A10, A20, A50)
dfA1 <- data.frame(x = rep(1:20, 60),
                   y = as.vector(t(A)),
                   maxlen = c(rep("10", 400), rep("20", 400), rep("50", 400)))
dfA2 <- data.frame(x = 1:20,
                   y = colMeans(A10),
                   maxlen = rep("10", 20))
dfA3 <- data.frame(x = 1:20,
                   y = colMeans(A20),
                   maxlen = rep("20", 20))
dfA4 <- data.frame(x = 1:20,
                   y = colMeans(A50),
                   maxlen = rep("50", 20))

## read data and construct data frames for panel B
B0<-as.matrix(read.csv("./simulation_results/max50_missed0.csv", header = FALSE, sep = ","), byrow = TRUE)
B1<-as.matrix(read.csv("./simulation_results/max50_missed1.csv", header = FALSE, sep = ","), byrow = TRUE)
B<-rbind(B0, B1)
dfB1 <- data.frame(x = rep(1:20, 40),
                  y = as.vector(t(B)),
                  missed_cleavages = c(rep("0", 400), rep("1", 400)))
dfB2 <- data.frame(x = 1:20,
                  y = colMeans(B[1:20,]),
                  missed_cleavages = rep("0", 20))
dfB3 <- data.frame(x = 1:20,
                   y = colMeans(B[21:33,]),
                   missed_cleavages = rep("1", 20))

## read data and construct data frames for panel C
C4 <- as.matrix(read.csv("./simulation_results/read4_max5_100_missed1.csv", header = FALSE, sep = ","), byrow = TRUE)
C10 <- as.matrix(read.csv("./simulation_results/read10_max5_100_missed1.csv", header = FALSE, sep = ","), byrow = TRUE)
C <- rbind(C4, C10)
dfC1 <- data.frame(x = rep(5*(1:20), 40),
                   y = as.vector(t(C)),
                   n_read = c(rep("4", 400), rep("10", 400)))
dfC2 <- data.frame(x = 5*(1:20),
                   y = colMeans(C4),
                   n_read = rep("4", 20))
dfC3 <- data.frame(x = 5*(1:20),
                   y = colMeans(C10),
                   n_read = rep("10", 20))

## read data and construct data frames for panel D
D4 <- read.csv("./simulation_results/read4_max50_missed1.csv", header = FALSE, sep = ",")
D5 <- read.csv("./simulation_results/read5_max50_missed1.csv", header = FALSE, sep = ",")
D6 <- read.csv("./simulation_results/read6_max50_missed1.csv", header = FALSE, sep = ",")
D <- rbind(D4, D5, D6)
nA<-c(); nR<-c(); nN<-c(); nD<-c(); nC<-c(); nE<-c(); nQ<-c(); nG<-c(); nH<-c(); nI<-c();
nL<-c(); nK<-c(); nM<-c(); nF<-c(); nP<-c(); nS<-c(); nT<-c(); nW<-c(); nY<-c(); nV<-c();
nX<-c()
PLURs <- 20322 - as.numeric(D[,2]) # PLUR = Proteins Lacking Unique Reads
for (i in 1:dim(D)[1]) {
  nA[i] <- sum(charToRaw(D[i,1]) == charToRaw('A'))
  nR[i] <- sum(charToRaw(D[i,1]) == charToRaw('R'))
  nN[i] <- sum(charToRaw(D[i,1]) == charToRaw('N'))
  nD[i] <- sum(charToRaw(D[i,1]) == charToRaw('D'))
  nC[i] <- sum(charToRaw(D[i,1]) == charToRaw('C'))
  nE[i] <- sum(charToRaw(D[i,1]) == charToRaw('E'))
  nQ[i] <- sum(charToRaw(D[i,1]) == charToRaw('Q'))
  nG[i] <- sum(charToRaw(D[i,1]) == charToRaw('G'))
  nH[i] <- sum(charToRaw(D[i,1]) == charToRaw('H'))
  nI[i] <- sum(charToRaw(D[i,1]) == charToRaw('I'))
  nL[i] <- sum(charToRaw(D[i,1]) == charToRaw('L'))
  nK[i] <- sum(charToRaw(D[i,1]) == charToRaw('K'))
  nM[i] <- sum(charToRaw(D[i,1]) == charToRaw('M'))
  nF[i] <- sum(charToRaw(D[i,1]) == charToRaw('F'))
  nP[i] <- sum(charToRaw(D[i,1]) == charToRaw('P'))
  nS[i] <- sum(charToRaw(D[i,1]) == charToRaw('S'))
  nT[i] <- sum(charToRaw(D[i,1]) == charToRaw('T'))
  nW[i] <- sum(charToRaw(D[i,1]) == charToRaw('W'))
  nY[i] <- sum(charToRaw(D[i,1]) == charToRaw('Y'))
  nV[i] <- sum(charToRaw(D[i,1]) == charToRaw('V'))
  nX[i] <- 1
}
min(sum(nA), sum(nR), sum(nN), sum(nD), sum(nC), sum(nE), sum(nQ), sum(nG), sum(nH), sum(nI), sum(nL), sum(nK), sum(nM), sum(nF), sum(nP), sum(nS), sum(nT), sum(nW), sum(nY), sum(nV))
fit <- lm(PLURs ~ nA+nR+nN+nD+nC+nE+nQ+nG+nH+nI+nL+nK+nM+nF+nP+nS+nT+nW+nY+nV)
PLUR_coeffs <- coefficients(fit)[1] + as.vector(coefficients(fit))[2:21] # add intercept
aa_freqs <- read.csv(file = "./simulation_results/amino_acid_frequencies.csv", header = FALSE, sep = ",")
colnames(aa_freqs) <- c("aa", "x")
dfD <- data.frame(x = as.vector(aa_freqs[2]/sum(aa_freqs[2])),
                  y = PLUR_coeffs[1:20],
                  aa <- aa_freqs[1])


## read data and construct data frames for Figure 3
A50<-as.matrix(read.csv("./simulation_results/max50_missed1.csv", header = FALSE, sep = ","), byrow = TRUE)
NI<-as.matrix(read.csv("./simulation_results/non_ideal_max50_missed1_first20.csv", header = FALSE, sep = ","), byrow = TRUE)
NI<-rbind(A50[1,1:20], NI)
dfNI1 <- data.frame(x = rep(1:20, 2),
                    y = as.vector(t(NI)),
                    maxlen = c(rep("ideal", 20), rep("non-ideal",20)))


## plot number of protein with unique reads as function of read amino acids
plotA <- ggplot(NULL, aes(x, y, color = maxlen)) +
  scale_color_nejm() +
  geom_point(
    data = dfA1,
    alpha = 0.4,
    shape = 16,
    size = 3,
    stroke = 0
  ) +
  geom_line(data = dfA2, size = 1) +
  geom_line(data = dfA3, size = 1) +
  geom_line(data = dfA4, size = 1) +
  scale_x_continuous(limits = c(1,20), breaks = c(1:9,10,12,14,16,18,20), minor_breaks = 1:20) +
  scale_y_log10(limits = c(200,21000), breaks = c(300, 1000, 3000, 10000), minor_breaks = rep(1:9, 21)*(10^rep(-10:10, each=9))) +
  xlab('readable amino acids') +
  ylab('proteins lacking unique PPRM') +
  labs(tag = "A", color = "maximum read length") +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(legend.position = c(0.737, 0.812),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = alpha("white", 0)))

## plot number of protein with unique reads as function of read amino acids
plotB <- ggplot(NULL, aes(x, y, color = missed_cleavages)) +
  scale_color_nejm() +
  geom_point(
    data = dfB1,
    alpha = 0.4,
    shape = 16,
    size = 3,
    stroke = 0
  ) +
  geom_line(data = dfB2, size = 1) +
  geom_line(data = dfB3, size = 1) +
  scale_x_continuous(limits = c(1,20), breaks = c(1:9,10,12,14,16,18,20), minor_breaks = 1:20) +
  scale_y_log10(limits = c(200,21000), breaks = c(300, 1000, 3000, 10000), minor_breaks = rep(1:9, 21)*(10^rep(-10:10, each=9))) +
  xlab('readable amino acids') +
  ylab('proteins lacking unique PPRM') +
  labs(tag = "B", color = "missed cleavages") +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(legend.position = c(0.77, 0.85),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = alpha("white", 0)))

## plot number of protein with unique reads as function of read amino acids
plotC <- ggplot(NULL, aes(x, y, color = factor(sort(as.numeric(n_read))))) +
  scale_color_nejm() +
  geom_point(
    data = dfC1,
    alpha = 0.4,
    shape = 16,
    size = 3,
    stroke = 0
  ) +
  geom_line(data = dfC2,
            size = 1) +
  geom_line(data = dfC3,
            size = 1) +
  scale_x_continuous(limits = c(0, 100), breaks = 10*(0:10), minor_breaks = 5*0:20) +
  scale_y_log10(limits = c(200, 21000), breaks = c(300, 1000, 3000, 10000), minor_breaks = rep(1:9, 21)*(10^rep(-10:10, each=9))) +
  xlab("maximum read length") +
  ylab("proteins lacking unique PPRM") +
  labs(tag = "C", color = "readable amino acids") +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(legend.position = c(0.74, 0.85),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = alpha("white", 0)))

## plot number of protein with unique reads as function of read amino acids
plotD <- ggplot(NULL, aes(x, y, label = aa)) +
  scale_color_nejm() +
  geom_point(
    data = dfD,
    alpha = 0.4,
    shape = 16,
    size = 7,
    stroke = 0,
    color = rgb(0,114/256,181/256)
  ) +
  geom_text(data = dfD, aes(fontface = "bold"), nudge_x = 0, nudge_y = 1, color = rgb(0,114/256,181/256)) +

  #ylim(19750, 20100) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = 1:10/100, minor_breaks = NULL) +
  xlab('amino acid frequency') +
  ylab('discriminative power (proteins)') +
  labs(tag = "D") +
  theme_bw() +
  theme(aspect.ratio=1)

## plot ideal vs non-ideal behavior (for read length <=50, 1 missed)
plotNI <- ggplot(NULL, aes(x, y, color = maxlen)) +
  scale_color_nejm() +
  geom_point(
    data = dfNI1,
    alpha = 0.4,
    shape = 16,
    size = 3,
    stroke = 0
  ) +
  #geom_line(data = dfNI2, size = 1) +
  #geom_line(data = dfNI3, size = 1) +
  scale_x_continuous(limits = c(1,20), breaks = c(1:9,10,12,14,16,18,20), minor_breaks = 1:20) +
  scale_y_log10(limits = c(200,21000), breaks = c(300, 1000, 3000, 10000), minor_breaks = rep(1:9, 21)*(10^rep(-10:10, each=9))) +
  xlab('readable amino acids') +
  ylab('proteins lacking unique PPRM') +
  labs(tag = " ", color = "conditions") +
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(legend.position = c(0.83, 0.85),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = alpha("white", 0)))


## plot all four panels A-D
grid.arrange(plotA,plotB,plotC,plotD, ncol=2, padding = NULL)
