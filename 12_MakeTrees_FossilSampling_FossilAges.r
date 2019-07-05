library(phangorn)
library(TreeSim)

#-----------------------------------
#settings
treeAge <- 100
nbTaxa <- 12
nbReps <- 20

Factor <- "number of fossils"
Values <- nbTaxa * c(1/12, 1/6, 1/4, 1/3, 1/2)

oldestFossils <- seq(from = 85, to = 15, by = -17.5)
terminalBr <- 5
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
# -------------------------------------------------------------------------------------------
#returns a vector containing the indices of the tips which are descendents of "nodeIndex"
list.descendents.of<-function(tree, nodeIndex){
    if(is.na(nodeIndex)) return(NA)
    if(nodeIndex<1) return(NA)
    
    parent <- tree$edge[,1]
    child <- tree$edge[,2]
    root <- as.integer(parent[!match(parent, child, 0)][1]) #match gives the positions at which the entries of the parent vector occur in the child vector. If =0, then no match and =root node
    tips <- child[is.na(!match(child,parent))]
    
    #if it is a tip
    if(length(which(tips==nodeIndex))>0){
        return(nodeIndex)
    }
    
    is.desc <- c(rep(FALSE,length(tips)))
    
    for (i in 1:length(tips)){
        index <- tips[i]
        while( index != root ){
            index <- parent[child==index]
            if(index == nodeIndex){
                is.desc[i] <- TRUE
                break
            }
        }
    }
    return(tips[is.desc])
}
# -------------------------------------------------------------------------------------------
#finds all branches in a tree that only have taxa starting with certain letters as descendents
get.fossil.branches <- function(tree, pattern = "f"){
    nbBr <- length(tree$edge.length)
    fossil.branches <- NULL
    for (br in 1:nbBr) {
        descs <- tree$tip.label[list.descendents.of(tree, tree$edge[br, 2])]
        if (length(grep(pattern = paste("^", pattern, sep = ""), descs)) == length(descs)){
            fossil.branches <- c(fossil.branches, br)
        }
    }
    return(fossil.branches)
}
#-----------------------------------

#b0: oldest=85 to b4: oldest=15

analyses <- paste("b", 0:4, sep = "")
BRate <- log(nbTaxa / 2) / treeAge
for (b in 0:4) {
    analysis <- analyses[b +1]
    for (t in 1:nbReps) {
        treeExtant <- sim.bd.taxa.age(nbTaxa, 1, BRate, 0, frac=1, age=treeAge, mrca=T)[[1]]
        for (f in 1:length(Values)){
            nbFossils <- Values[f]
            tree <- treeExtant
            if (nbFossils == 1){
                fossilAges <- oldestFossils[b +1]
            } else {
                fossilAges <- seq(from = oldestFossils[b +1], to = oldestFossils[b +1] -5, by = -5/(Values[f] - 1))
            }
        
            for (k in 1:nbFossils){
                fossil.attach <- fossilAges[k] + terminalBr
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
                tree$edge.length[which(tree$edge[,2] == which(tree$tip.label == fossil.label))] <- terminalBr
            }
            write.tree(tree, file=paste(analysis, "_", treeFilename, sep = ""), append=T)
        }
    }
}
#------------------------------------------------------------------------------- 
# plotting
#------------------------------------------------------------------------------- 
# 
sets <- paste("b", 0:4, "_", sep = "")
pdf("Plots_FossilSamplingTrees.pdf", width = 100, height = 140, onefile = T)
nb <- 7
factor <- length(Values)
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
