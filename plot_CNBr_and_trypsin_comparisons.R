library(ggplot2)
library(ggsci)
library(gridExtra)
setwd('C:/Users/magnus/Desktop/Next Generation Proteomics/For_GitHub')
N_PROTEINS <- 20322 # (maximum) number of proteins that could be identified (for reference)



## read data and construct data frames for panel D
D4 <- read.csv("./simulation_results/CNBr_read4_max50_missed1.csv", header = FALSE, sep = ",")
D5 <- read.csv("./simulation_results/CNBr_read5_max50_missed1.csv", header = FALSE, sep = ",")
D6 <- read.csv("./simulation_results/CNBr_read6_max50_missed1.csv", header = FALSE, sep = ",")
D <- rbind(D4, D5, D6)
nA<-c(); nR<-c(); nN<-c(); nD<-c(); nC<-c(); nE<-c(); nQ<-c(); nG<-c(); nH<-c(); nI<-c();
nL<-c(); nK<-c(); nM<-c(); nF<-c(); nP<-c(); nS<-c(); nT<-c(); nW<-c(); nY<-c(); nV<-c();
nX<-c()
PLURs <- N_PROTEINS - as.numeric(D[,2]) # PLUR = Proteins Lacking Unique Reads
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
dfD1 <- data.frame(x = as.vector(aa_freqs[2]/sum(aa_freqs[2])),
                   y = PLUR_coeffs[1:20],
                   aa <- aa_freqs[1])

D4 <- read.csv("./simulation_results/CNBr_read4_max50_missed0.csv", header = FALSE, sep = ",")
D5 <- read.csv("./simulation_results/CNBr_read5_max50_missed0.csv", header = FALSE, sep = ",")
D6 <- read.csv("./simulation_results/CNBr_read6_max50_missed0.csv", header = FALSE, sep = ",")
D <- rbind(D4, D5, D6)
nA<-c(); nR<-c(); nN<-c(); nD<-c(); nC<-c(); nE<-c(); nQ<-c(); nG<-c(); nH<-c(); nI<-c();
nL<-c(); nK<-c(); nM<-c(); nF<-c(); nP<-c(); nS<-c(); nT<-c(); nW<-c(); nY<-c(); nV<-c();
nX<-c()
PLURs <- N_PROTEINS - as.numeric(D[,2]) # PLUR = Proteins Lacking Unique Reads
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
dfD0 <- data.frame(x = as.vector(aa_freqs[2]/sum(aa_freqs[2])),
                   y = PLUR_coeffs[1:20],
                   aa <- aa_freqs[1])
D4 <- read.csv("./simulation_results/CNBr_read4_max50_missed1.csv", header = FALSE, sep = ",")
D5 <- read.csv("./simulation_results/CNBr_read5_max50_missed1.csv", header = FALSE, sep = ",")
D6 <- read.csv("./simulation_results/CNBr_read6_max50_missed1.csv", header = FALSE, sep = ",")
D <- rbind(D4, D5, D6)
nA<-c(); nR<-c(); nN<-c(); nD<-c(); nC<-c(); nE<-c(); nQ<-c(); nG<-c(); nH<-c(); nI<-c();
nL<-c(); nK<-c(); nM<-c(); nF<-c(); nP<-c(); nS<-c(); nT<-c(); nW<-c(); nY<-c(); nV<-c();
nX<-c()
PLURs <- N_PROTEINS - as.numeric(D[,2]) # PLUR = Proteins Lacking Unique Reads
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
dfD1 <- data.frame(x = as.vector(aa_freqs[2]/sum(aa_freqs[2])),
                   y = PLUR_coeffs[1:20],
                   aa <- aa_freqs[1])
#D0 with trypsin
D4 <- read.csv("./simulation_results/trypsin_read4_max50_missed0.csv", header = FALSE, sep = ",")
D5 <- read.csv("./simulation_results/trypsin_read5_max50_missed0.csv", header = FALSE, sep = ",")
D6 <- read.csv("./simulation_results/trypsin_read6_max50_missed0.csv", header = FALSE, sep = ",")
D <- rbind(D4, D5, D6)
nA<-c(); nR<-c(); nN<-c(); nD<-c(); nC<-c(); nE<-c(); nQ<-c(); nG<-c(); nH<-c(); nI<-c();
nL<-c(); nK<-c(); nM<-c(); nF<-c(); nP<-c(); nS<-c(); nT<-c(); nW<-c(); nY<-c(); nV<-c();
nX<-c()
PLURs <- N_PROTEINS - as.numeric(D[,2]) # PLUR = Proteins Lacking Unique Reads
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
dfD0 <- data.frame(x = as.vector(aa_freqs[2]/sum(aa_freqs[2])),
                   y = PLUR_coeffs[1:20],
                   aa <- aa_freqs[1])



## plot number of protein with unique reads as function of read amino acids for CNBr
plotD0 <- ggplot(NULL, aes(x, y, label = aa)) +
  scale_color_nejm() +
  geom_point(
    data = dfD0,
    alpha = 0.4,
    shape = 16,
    size = 7,
    stroke = 0,
    color = rgb(0,181/256,114/256)
  ) +
  geom_text(data = dfD0, aes(fontface = "bold"), nudge_x = 0, nudge_y = 1, color = rgb(0,181/256,114/256)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = 1:10/100, minor_breaks = NULL) +
  xlab('amino acid frequency') +
  ylab('discriminative power (proteins)') +
  labs(tag = "A") +
  theme_bw() +
  theme(aspect.ratio=1)

plotD1 <- ggplot(NULL, aes(x, y, label = aa)) +
  scale_color_nejm() +
  geom_point(
    data = dfD1,
    alpha = 0.4,
    shape = 16,
    size = 7,
    stroke = 0,
    color = rgb(0,181/256,114/256)
  ) +
  geom_text(data = dfD1, aes(fontface = "bold"), nudge_x = 0, nudge_y = 1, color = rgb(0,181/256,114/256)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = 1:10/100, minor_breaks = NULL) +
  xlab('amino acid frequency') +
  ylab('discriminative power (proteins)') +
  labs(tag = "A") +
  theme_bw() +
  theme(aspect.ratio=1)

## plot all four panels A-D
grid.arrange(plotA, plotD, plotD0, plotD1, ncol = 2, padding = NULL)

