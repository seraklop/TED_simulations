
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

lowerRate <- c(0.5, 0.25, 0.1) # number of observed states

tipLabelVec <- c(paste("t", 1:nbTaxa, sep = ""), paste("f", 1:nbFossils, sep = ""))
Colours <- brewer.pal(11,"Spectral")
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
# We assume that the underlying characters is an ordered, equal-frequency characters

for (b in 0:(length(morphoRate) -1) ){
    for (i in 1 : (length(lowerRate)) ){
        for (j in 1:nbReps){
            dataFile <- paste("Data_files/MorphoData_b", b, "_", nbTaxa, "Taxa_",i, "_", j, ".phy", sep = "")
            tree <- trees[[j]]
            nbVarSites <- 0
            morphoDat <- NULL
            while (nbVarSites < nbMorpho ){
                mD <- simSeq(tree, l = nbMorpho, type = "DNA", bf=c(lowerRate[i], 1-lowerRate[i], 0, 0), rate = morphoRate[b + 1])
                mD <- as.data.frame(mD)
                mD[mD == "c"] <- "0"
                mD[mD == "a"] <- "1"
                
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