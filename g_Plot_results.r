# analyze and make plots of the simulation results
library(RColorBrewer)

Factor <- "asymmetry: frequency of state 1"
values <- c(0.5, 0.25, 0.1)
nReps <- 20

subsets <- paste("b", 0:4, sep = "")
setNames <- c("rate = 0.25", "rate = 0.5", "rate = 1.0", "rate = 2.0", "rate = 10.0")
Colours <- brewer.pal(11,"Spectral")

#-----------------------------------------------------------------------------------------------------------------------
# first summarize most important ones across sets, then for each set, print details
matList <- list(NULL)
setColours <- c(Colours[1], Colours[4], Colours[6], Colours[9], Colours[11])
for (i in 1:length(subsets)){
    inFile <- list.files(pattern="*.summary")[grep(subsets[i], list.files(pattern="*.summary"))]
    res<- read.table(inFile, header = T)
    res <- res[with(res, order(Fac, Rep)), ]
    matList[[i]] <- matrix(res$RootMedian, nrow = nReps, ncol = length(values), byrow = F, dimnames = list(NULL, values))
}

pdf(file = "Summary_TreeAges.pdf", width = 6.4, height = 5)
nSets <- length(subsets)
nFactors <- length(values)
mat <- matrix(nrow = nReps, ncol = nSets * nFactors)
for (i in 1:nFactors) {
    for (j in 1:nSets){
        mat[ , nSets * (i -1) +j] <- matList[[j]][ ,i]
    }
}
boxplot(mat, xlab = Factor, ylab = "median root age", ylim = c(0,1000), xaxt = "n", frame.plot = F)
shaded <- seq(from = 2, to = nFactors, by = 2)
for (i in shaded ){
  polygon(x = c((i -1) * nSets +0.5, i * nSets +0.5, i * nSets +0.5, (i -1) * nSets +0.5), y = c(-50, -50, 1000, 1000), density = NA, col = "grey80", border = NA)
}
boxplot(mat, xlab = Factor, ylab = "median root age", ylim = c(0,1000), col = setColours, xaxt = "n", frame.plot = F, add = T)
axis(side = 1, at = seq(from = (nSets + 1) / 2, to = nSets*length(values), by = nSets), labels = values)
abline(h = 100, lty = 2, lwd = 2, col = "darkred")

legend(x = 2, y = 1000, legend = setNames, cex = 0.8, fill = setColours, bty = "n")

dev.off()



#-----------------------------------------------------------------------------------------------------------------------
for (i in 1:length(subsets)){
    set <- subsets[i]
    inFile <- list.files(pattern="*.summary")[grep(subsets[i], list.files(pattern="*.summary"))]
    res<- read.table(inFile, header = T)
    res <- res[with(res, order(Fac, Rep)), ]

    #-----------------------------------------------------------------------------------------------------------------------
    # lnL, ASDSF, topo accuracy
    pdf(file = paste(set, "_Summary_lnL_convergence_topoAccuracy.pdf", sep = ""), width = 6.4, height = 12)
    layout(matrix(1:3, ncol=1, nrow = 3) )

    mat <- matrix(res$lnLharmMean, nrow = nReps, ncol = length(values), byrow = F, dimnames = list(NULL, values))
    boxplot(mat, xlab = Factor, ylab = "harmonic mean of lnL", col = Colours[1])

    mat <- matrix(res$ASDSF, nrow = nReps, ncol = length(values), byrow = F, dimnames = list(NULL, values))
    boxplot(mat, xlab = Factor, ylab = "ASDSF topo convergence", col = Colours[2])
    abline(h = 0.01, lty = 3, lwd = 2, col = "darkred")

    mat1 <- matrix(res$propCorr_wFos, nrow = nReps, ncol = length(values), byrow = F, dimnames = list(NULL, values))
    mat2 <- matrix(res$propCorr_woFos, nrow = nReps, ncol = length(values), byrow = F, dimnames = list(NULL, values))
    mat <- matrix(nrow = nReps, ncol = 2*length(values))
    for (i in 1:dim(mat1)[2]) {
        mat[,2*(i-1) +1] <- mat1[, i]
        mat[,2*i] <- mat2[, i]
    }
    boxplot(mat, xlab = Factor, ylab = "proportion of correctly recovered nodes", col = c(Colours[3], Colours[4]), ylim = c(0,1), xaxt = "n")
    axis(side = 1, at = seq(from = 1.5, to = 2*length(values), by = 2), labels = values)

    dev.off()
    #-----------------------------------------------------------------------------------------------------------------------
    # Tree age
    pdf(file = paste(set, "_Summary_TreeAge.pdf", sep = ""), width = 6.4, height = 12)
    layout(matrix(1:3, ncol=1, nrow = 3) )

    mat <- matrix(res$RootMedian, nrow = nReps, ncol = length(values), byrow = F, dimnames = list(NULL, values))
    boxplot(mat, xlab = Factor, ylab = "median root age", col = Colours[5], ylim = c(0,1000))
    abline(h = 100, lty = 3, lwd = 2, col = "darkred")

    mat <- matrix(res$RootUpper, nrow = nReps, ncol = length(values), byrow = F, dimnames = list(NULL, values))
    boxplot(mat, xlab = Factor, ylab = "upper and lower of root age", col = Colours[6], ylim = c(0,1000))
    mat <- matrix(res$RootLower, nrow = nReps, ncol = length(values), byrow = F, dimnames = list(NULL, values))
    boxplot(mat, col = Colours[7], add = T)

    mat <- matrix(res$RootPSRF, nrow = nReps, ncol = length(values), byrow = F, dimnames = list(NULL, values))
    boxplot(mat, xlab = Factor, ylab = "PSRF (convergence) of root age", col = Colours[8])

    dev.off()
    #-----------------------------------------------------------------------------------------------------------------------
    # m1 and m2: rate multiplier for morphology and molecules
    pdf(file = paste(set, "_Summary_RateMultipliers.pdf", sep = ""), width = 6.4, height = 12)
    layout(matrix(1:3, ncol=1, nrow = 3) )

    mat <- matrix(res$m1Median, nrow = nReps, ncol = length(values), byrow = F, dimnames = list(NULL, values))
    boxplot(mat, xlab = Factor, ylab = "median rate multiplier morphology/molecules", col = Colours[9] )
    mat <- matrix(res$m2Median, nrow = nReps, ncol = length(values), byrow = F, dimnames = list(NULL, values))
    boxplot(mat, col = Colours[11], add = T)

    mat <- matrix(res$m1Upper, nrow = nReps, ncol = length(values), byrow = F, dimnames = list(NULL, values))
    boxplot(mat, xlab = Factor, ylab = "upper and lower of morpho rate multiplier", col = Colours[10], ylim = c(0,max(res$m1Upper)))
    mat <- matrix(res$m1Lower, nrow = nReps, ncol = length(values), byrow = F, dimnames = list(NULL, values))
    boxplot(mat, col = Colours[11], add = T)

    mat <- matrix(res$m1PSRF, nrow = nReps, ncol = length(values), byrow = F, dimnames = list(NULL, values))
    boxplot(mat, xlab = Factor, ylab = "PSRF (convergence) of morpho rate multiplier", col = Colours[8])

    dev.off()
}
#-----------------------------------------------------------------------------------------------------------------------