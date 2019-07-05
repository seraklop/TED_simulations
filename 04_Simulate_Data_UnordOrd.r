
library(phangorn)
library(RColorBrewer)
#-----------------------------------
#settings
nbReps <- 20
treeAge <- 100
nbTaxa <- 12
nbFossils <- nbTaxa/3
nbSites <- 1000
nbMorpho <- 100
molecRate <- 0.25/treeAge
morphoRate <- c(0.25, 0.5, 1, 2, 10) / treeAge
propMorpho <- 0.5
TreeFolder <- "../../1_Priors_TreeAge_CR/a_MakeTrees-D-simulateSeqs/"

nbStates <- c(2, 3, 6, 10) # number of observed states

tipLabelVec <- c(paste("t", 1:nbTaxa, sep = ""), paste("f", 1:nbFossils, sep = ""))
Colours <- brewer.pal(length(nbStates),"Spectral")
#------------------------------------------------------------------------------- 
#(1) simulate one molecular dataset on each tree
trees <- read.tree(paste(TreeFolder, "Trees_", nbTaxa, "Taxa_20.tre", sep = ""))

for (i in 1:nbReps){
    treeA <- trees[[i]]
    GeneA <- as.character(simSeq(treeA, l = nbSites, type = "DNA", rate = molecRate) )
    GeneA[grep("f", rownames(GeneA)), ] <- rep("?", nbSites)
    GeneA <- GeneA[tipLabelVec, ]
    GeneA <- as.phyDat(GeneA)
    write.phyDat(GeneA, file=paste("Data_files/Seqs_", nbTaxa, "Taxa_", i, ".phy", sep = ""), nbcol=-1, append=F)
}
#------------------------------------------------------------------------------- 
get.ordered.Qmat <- function(nStates){
    if (nStates == 2) {
        Qmat <- rbind(c(-1,1), c(1,-1))
    } else if (nStates > 2){
        Qmat <- matrix(rep(0, nStates^2), ncol = nStates, nrow = nStates)
        Qmat[1,1] <- -1
        Qmat[1,2] <- 1
        Qmat[nStates,nStates] <- -1
        Qmat[nStates, nStates-1] <- 1
        for (s in 2:(nStates -1) ){
            Qmat[s,s] <- -2
            Qmat[s,s-1] <- 1
            Qmat[s,s+1] <- 1
        }
    } else {
        Qmat <- NA
    }
    return(Qmat)
}
#------------------------------------------------------------------------------- 
# We assume that the underlying characters is an ordered, equal-frequency characters

for (b in 0:4){
    for (i in 1 : (length(nbStates)) ){
        for (j in 1:nbReps){
            dataFile <- paste("Data_files/MorphoData_b", b, "_", nbTaxa, "Taxa_",i, "_", j, ".phy", sep = "")
            tree <- trees[[j]]
            nbVarSites <- 0
            morphoDat <- NULL
            Qmat_ordered <- get.ordered.Qmat(nStates = nbStates[i])
            while (nbVarSites < nbMorpho ){
                mD <- simSeq(tree, l = nbMorpho, type = "USER", Q = Qmat_ordered, levels = paste(0:(nbStates[i] -1)), rate = morphoRate[b +1])
                mD <- as.data.frame(mD)
    
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