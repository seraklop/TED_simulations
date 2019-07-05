library(phangorn)
library(TreeSim)

#-----------------------------------
#settings
treeAge <- 100
nbTaxa <- c(12, 24, 48, 96)
nbFossils <- nbTaxa / 3
nbReps <- 10 #for now
treeFilename <- "Trees_12-96Taxa.tre"
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
#target nb of taxa
#fossils have fixed ages (here evenly from 15-85 million years) and are all subtended by 
#a 5 My branch.
#attach them to the tree randomly, given the age of the attachment point.

for (i in 1:length(nbTaxa)){
    BRate <- log(nbTaxa[i]/2) / treeAge
    simtrees <- sim.bd.taxa.age(nbTaxa[i], nbReps, BRate, 0, frac=1, age=treeAge, mrca=T)
    
    fossilAges <- seq(from = 85, to = 15, by = -(85-15)/(nbFossils[i] - 1))
    
    for (j in 1:nbReps) {
        tree <- simtrees[[j]]
        for (k in 1:nbFossils[i]){
            fossil.attach <- fossilAges[k] + 5
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
            tree$edge.length[which(tree$edge[,2] == which(tree$tip.label == fossil.label))] <- 5
        }
        write.tree(tree, file=treeFilename, append=T)
    }
}


#------------------------------------------------------------------------------- 
# plotting
#------------------------------------------------------------------------------- 
# 
all.trees <- read.tree(treeFilename)
pdf("Plots_AllTrees_round1.pdf", width = 64, height = 100)
nb <- length(nbTaxa)
toPlot <- 5
layout(matrix(1:(nb*toPlot), ncol=nb, nrow = toPlot) )
for (i in 1:nb) {
    for (j in 1:toPlot){
        tree <- all.trees[[(i-1)*nbReps +j]]
        plot(tree, show.tip.label = F)
    }
}

dev.off()

#------------------------------------------------------------------------------- 
# making input for MB files: fossil calibrations

i=3
fossilAges <- seq(from = 85, to = 15, by = -(85-15)/(nbFossils[i] - 1))
outString <- ""
for (j in 1:length(fossilAges)){
    outString <- paste(outString, " f", j, "=fixed(", fossilAges[j], ")", sep = "")
}
outString
