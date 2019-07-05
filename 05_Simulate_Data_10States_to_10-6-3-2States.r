
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
molecRate <- 0.25/treeAge
morphoRate <- c(0.25, 1, 10) / treeAge
propMorpho <- 0.5
TreeFolder <- "../../1_Priors_TreeAge_CR/a_MakeTrees-D-simulateSeqs/"

nbStates <- c(10, 6, 3, 2) # number of observed states

tipLabelVec <- c(paste("t", 1:nbTaxa, sep = ""), paste("f", 1:nbFossils, sep = ""))
Colours <- brewer.pal(length(nbStates),"Spectral")
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
# We assume that the underlying characters is an ordered, equal-frequency 10-state character. Then we randomly
# unite characters to end up with 'nbStates' states

Qmat_ordered <- matrix(rep(0, 100), ncol = 10, nrow = 10)
Qmat_ordered[1,1] <- -1
Qmat_ordered[1,2] <- 1
Qmat_ordered[10,10] <- -1
Qmat_ordered[10,9] <- 1
for (s in 2:9){
    Qmat_ordered[s,s] <- -2
    Qmat_ordered[s,s-1] <- 1
    Qmat_ordered[s,s+1] <- 1
}
for (b in 0:2){
    for (i in 1 : (length(nbStates)) ){
        for (j in 1:nbReps){
            dataFile <- paste("Data_files/MorphoData_b", b, "_", nbTaxa, "Taxa_",i, "_", j, ".phy", sep = "")
            tree <- trees[[j]]
            nbVarSites <- 0
            morphoDat <- NULL
            while (nbVarSites < nbMorpho ){
                mD <- simSeq(tree, l = nbMorpho, type = "USER", Q = Qmat_ordered, levels = paste(0:9), rate = morphoRate[b +1])
                mD <- as.data.frame(mD)
    
                #combine states into one if needed
                if (nbStates[i] != 10) {
                    currStates <- unique(unlist(mD))
                    tempStates <- letters[1:nbStates[i]]
                    newStates <- (1:nbStates[i]) -1
                    for( st in currStates ){
                        mD[mD == st] <- sample(tempStates, size = 1, replace = T)
                    }
                    for( nst in 1:length(newStates) ){
                        mD[mD == tempStates[nst]] <- newStates[nst]
                    }
                }
    
                #remove part of the morpho data for the fossils
                for ( fossil in 1: nbFossils ){
                    fossilID <- which(names(mD) == paste("f", fossil, sep = "") )  
                    mD[[fossilID]][sample(1:nbMorpho, size = nbMorpho*(1-propMorpho), replace = F )] <- "?"
                }
                
                #only use the variable positions, add these to the matrix
                isVariable <- rep(F, dim(mD)[1])
                for ( s in 1:dim(mD)[1] ){
                    this.states <- unique(as.character(mD[s,]))
                    this.states <- this.states[this.states != "?"]
                    if (length(this.states) >= 2){
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
            
            morphoDat <- morphoDat[, tipLabelVec]
            morphoDat<-as.phyDat(morphoDat, levels=unique(unlist(morphoDat)), type="USER")
            write.phyDat(morphoDat, file=dataFile, nbcol=-1, append=T)
            write("", file=dataFile, sep="\n", append=T)
        }
    }
}
#------------------------------------------------------------------------------- 