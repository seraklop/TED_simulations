library(phangorn)
library(RColorBrewer)
#-----------------------------------
#settings
nbReps <- 10
treeAge <- 100
nbTaxa <- 12
nbFossils <- nbTaxa/3
nbSites <- 1000
nbMorpho <- 100
propMorpho <- 0.5
molecRate <- 0.25/treeAge
morphoRate <- 0.25/treeAge
TreeFolder <- "../../1_Priors_TreeAge_CR/a_MakeTrees-D-simulateSeqs/"
tipLabelVec <- c(paste("t", 1:nbTaxa, sep = ""), paste("f", 1:nbFossils, sep = ""))

ClockVar <- c(0.00000000001, 0.05, 0.1, 1, 10) # new values
Colours <- brewer.pal(length(ClockVar),"Spectral")
#------------------------------------------------------------------------------- 
#(1) simulate one molecular dataset on each tree
trees <- read.tree(paste(TreeFolder, "Trees_", nbTaxa, "Taxa.tre", sep = ""))

for (i in 1:nbReps){
    treeA <- trees[[i]]
    GeneA <- as.character(simSeq(treeA, l = nbSites, type = "DNA", rate = molecRate) )
    GeneA[grep("f", rownames(GeneA)), ] <- rep("?", nbSites)
    GeneA <- GeneA[tipLabelVec, ]
    GeneA <- as.phyDat(GeneA)
    write.phyDat(GeneA, file=paste("Data_files/Seqs_", nbTaxa, "Taxa_", i, ".phy", sep = ""), nbcol=-1, append=F)
    
}
#------------------------------------------------------------------------------- 
# modify branch lengths for morpho simulations according to lognormal distribution (why not gamma?)
outTreeFile <- paste("NonclockTrees_", nbTaxa, "Taxa.tre")
outStats <- "clockvarStats.txt"
write("ClockVar TreeLength MeanBrlens MedianBrlens MinBrlens MaxBrlens", file = outStats)
pdf(file = "LogNormDistributions.pdf", width = 6.4, height = 5)

for (cl in 1 : length(ClockVar) ){
    # draw branch rate multipliers from a lognormal with mean=1 and var = ClockVar[i]
    mu <- - log(ClockVar[cl] +1) /2 
    sigma <- sqrt(log( ClockVar[cl] +1) )

    lognormFunc <- function (x){
        dlnorm(x, meanlog = mu, sdlog = sigma)
    }
    if (cl == 1){
        curve(lognormFunc, from = 0.0, to = 10, lwd=3, n = 200, col = Colours[cl], ylim = c(0,6))
        abline(v=1, lty=3, col="darkred")
    } else {
        curve(lognormFunc, from = 0.0, to = 10, lwd=3, n = 200, col = Colours[cl], add = T)
    }
    
    for (j in 1:nbReps){
        dataFile <- paste("Data_files/MorphoData_", nbTaxa, "Taxa_", cl, "_", j, ".phy", sep = "")
        ALRVTree <- trees[[j]]
        nbBranches <- length(ALRVTree$edge.length)
        brRateMulti <- rlnorm(nbBranches, meanlog = mu, sdlog = sigma)
        ALRVTree$edge.length <- ALRVTree$edge.length * brRateMulti
        write.tree(ALRVTree, file = outTreeFile, append = T)
        write(paste(ClockVar[i], sum(ALRVTree$edge.length), mean(ALRVTree$edge.length), median(ALRVTree$edge.length), min(ALRVTree$edge.length), max(ALRVTree$edge.length), sep = " "), file = outStats, append = T )

        #sim chars
        nbVarSites <- 0
        morphoDat <- NULL
        while (nbVarSites < nbMorpho ){
            mD <- simSeq(ALRVTree, l = nbMorpho, type = "USER", Q = cbind(c(0,1),c(1,0)), levels = c(0,1), rate = morphoRate)
            mD <- as.data.frame(mD)

            isVariable <- rep(F, dim(mD)[1])
            for ( s in 1:dim(mD)[1] ){
                if (length(unique(as.numeric(mD[s,]))) == 2){
                    isVariable[s] <- T
                }
            }
            if (length(morphoDat) == 0 ){
                morphoDat <- mD[isVariable,]
            } else {
                morphoDat <- rbind(morphoDat, mD[isVariable,])
            }
            nbVarSites <- nbVarSites + sum(isVariable)
        }
        morphoDat <- morphoDat[1:nbMorpho,]
            
        #remove part of the morpho data for the fossils
        for ( fossil in 1: nbFossils ){
            fossilID <- which(names(morphoDat) == paste("f", fossil, sep = "") )  
            morphoDat[[fossilID]][sample(1:nbMorpho, size = nbMorpho*(1-propMorpho), replace = F )] <- "?"
        }
        morphoDat <- morphoDat[, tipLabelVec]
        morphoDat<-as.phyDat(morphoDat, levels=c("1","0"), type="USER")
        write.phyDat(morphoDat, file=dataFile, nbcol=-1, append=T)
        write("", file=dataFile, sep="\n", append=T)
    }
}
    
legend(x = 5, y = 6, legend = ClockVar, col = Colours, lwd = 3, lty = 1)
dev.off()
    
#------------------------------------------------------------------------------- 
#plot some trees
nc_trees <- read.tree(outTreeFile)
pdf("NonClock_trees.pdf", width = 6.4, height = 6.4, onefile = T)
layout(matrix(1:25, nrow = 5, ncol = 5))
for (i in seq(from=1, to=50, by = 2)){
    tree <- nc_trees[i]
    plot(tree[[1]], edge.width = 2, show.tip.label=F, no.margin=T)  
}
dev.off()
