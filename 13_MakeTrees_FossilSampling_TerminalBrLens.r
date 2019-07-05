library(phangorn)
library(TreeSim)

#-----------------------------------
#settings
treeAge <- 100
nbTaxa <- 12
nbReps <- 20
nbFossils <- nbTaxa / 3

Factor <- "length of fossil branches"
terminalBr <- c(0, 5, 10, 20, 30)

youngestFossil <- 5
oldestFossil <- c(65, 40, 15)
treeFilename <- paste("Trees_FossilSampling_", nbTaxa, "_Taxa.tre", sep = "")

#-----------------------------------
# function returning for each branch in the tree the age of its basal and apical node
get.branch.periods <- function(tree, rootAge){
    nbBr <- length(tree$edge.length)
    for (br in 1:nbBr) {
        parent.age <- 0
        parent.node <- tree$edge[br, 1]
        stop <- F
        while(!stop){
            one.up <- which(tree$edge[,2] == parent.node)
            if ( length(one.up) > 0 ){
                parent.age <- parent.age + tree$edge.length[one.up]
                parent.node <- tree$edge[one.up, 1]
            } else {
                stop <- T
            }
        }
        parent.age <- rootAge - parent.age
        child.age <- parent.age - tree$edge.length[br]
        if (br == 1) {
            br.period.list <- list(c(parent.age, child.age))
        } else {
            this.period <- list(c(parent.age, child.age))
            br.period.list <- c(br.period.list, this.period)
        }
    }
    return(br.period.list)
}
#-----------------------------------
# simulate Yule tree under birth rate which maximizes the probability of observing the
# target nb of taxa
# fossils are evenly spaced trough time and attach randomly, on increasingly long 
# terminal branches.

BRate <- log(nbTaxa / 2) / treeAge
treesExtant <- sim.bd.taxa.age(nbTaxa, numbsim = nbReps, BRate, 0, frac=1, age=treeAge, mrca=T)

for (b in 0:(length(oldestFossil) -1)) {
    fossilAges <- seq(from = oldestFossil[b +1], to = youngestFossil, by = -(oldestFossil[b +1] - youngestFossil)/(nbFossils - 1))
    for (t in 1:nbReps) {
        treeExtant <- treesExtant[[t]]
        for (f in 1:length(terminalBr)){
            tree <- treeExtant
            for (k in 1:nbFossils){
                fossil.attach <- fossilAges[k] + terminalBr[f]
                branch.periods <- get.branch.periods(tree, treeAge)
                pot1  <- which(sapply(branch.periods, "[[", 1) > fossil.attach)
                pot2  <- which(sapply(branch.periods, "[[", 2) < fossil.attach)
                br.to.attach.to <- sample(intersect(pot1, pot2), 1)
                node.to.attach.to <- tree$edge[br.to.attach.to, 2]
                age.node <- branch.periods[[br.to.attach.to]][2]
                additBr <- rtree(2)
                fossil.label <- paste("f", k, sep = "")
                additBr$tip.label <- c("toRemove", fossil.label)
                tree <- bind.tree(tree, additBr, where = node.to.attach.to, position = (fossil.attach - age.node) )
                tree <- drop.tip(tree, tip = "toRemove")
                tree$edge.length[which(tree$edge[,2] == which(tree$tip.label == fossil.label))] <- terminalBr[f]
            }
            write.tree(tree, file=paste("b", b, "_", treeFilename, sep = ""), append=T)
        }
    }
}

#------------------------------------------------------------------------------- 
# plotting
#------------------------------------------------------------------------------- 
# 
sets <- paste("b", 0:(length(oldestFossil) -1), "_", sep = "")
pdf("Plots_FossilSamplingTrees.pdf", width = 100, height = 140, onefile = T)
nb <- 7
factor <- length(terminalBr)
layout(matrix(1:(nb*factor), nrow=nb, ncol = factor, byrow = T) )

for (s in sets) {
    all.trees <- read.tree(paste(s, treeFilename, sep = ""))
    for (i in 1:nb) {
        for (j in 1:factor){
            tree <- all.trees[[(i-1) *factor +j]]
            plot(tree, show.tip.label = F)
        }
    }
}
dev.off()

#------------------------------------------------------------------------------- 
