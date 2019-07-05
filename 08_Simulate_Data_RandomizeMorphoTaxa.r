
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
morphoRate <- 0.25/treeAge
propMorpho <- 0.5
TreeFolder <- "../../1_Priors_TreeAge_CR/a_MakeTrees-D-simulateSeqs/"

propTaxRand <- c(0, 0.25, 0.5, 0.75, 1) # remove phylogenetic signal (in morph chars only)

tipLabelVec <- c(paste("t", 1:nbTaxa, sep = ""), paste("f", 1:nbFossils, sep = ""))
Colours <- brewer.pal(length(propTaxRand),"Spectral")
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

for (i in 1 : (length(propTaxRand)) ){
    for (j in 1:nbReps){
        tree <- trees[[j]]

        nbVarSites <- 0
        morphoDat <- NULL
        while (nbVarSites < nbMorpho ){
            mD <- simSeq(tree, l = nbMorpho, type = "USER", Q = rbind(c(0,1), c(1,0)), levels = c("0", "1"), rate = morphoRate)
            mD <- as.data.frame(mD)

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

        #b0: randomize a portion of the fossil taxon labels for morpho
        b <- 0
        dataFile <- paste("Data_files/MorphoData_b", b, "_", nbTaxa, "Taxa_",i, "_", j, ".phy", sep = "")
        newMD <- morphoDat
        if(propTaxRand[i] > 0) {
            fossil.indices <- grep("f", names(newMD))
            toRand <- sample(fossil.indices, size = propTaxRand[i] * nbFossils, replace = F)
            names(newMD)[toRand] <- sample(names(morphoDat)[toRand], length(toRand), replace = F)
        }
        #remove part of the morpho data for the fossils
        for ( fossil in 1: nbFossils ){
            fossilID <- which(names(newMD) == paste("f", fossil, sep = "") )  
            newMD[[fossilID]][sample(1:nbMorpho, size = nbMorpho*(1-propMorpho), replace = F )] <- "?"
        }
                
        newMD <- newMD[, tipLabelVec]
        newMD<-as.phyDat(newMD, levels=unique(unlist(newMD)), type="USER")
        write.phyDat(newMD, file=dataFile, nbcol=-1, append=T)
        write("", file=dataFile, sep="\n", append=T)
        
        #b1: randomize a portion of the extant taxon labels for morpho
        b <- 1
        dataFile <- paste("Data_files/MorphoData_b", b, "_", nbTaxa, "Taxa_",i, "_", j, ".phy", sep = "")
        newMD <- morphoDat
        if(propTaxRand[i] > 0) {
            extant.indices <- grep("t", names(newMD))
            toRand <- sample(extant.indices, size = propTaxRand[i] * nbTaxa, replace = F)
            names(newMD)[toRand] <- sample(names(morphoDat)[toRand], length(toRand), replace = F)
        }
        #remove part of the morpho data for the fossils
        for ( fossil in 1: nbFossils ){
            fossilID <- which(names(newMD) == paste("f", fossil, sep = "") )  
            newMD[[fossilID]][sample(1:nbMorpho, size = nbMorpho*(1-propMorpho), replace = F )] <- "?"
        }
        ext.rand.MD <- newMD #save here for later
        newMD <- newMD[, tipLabelVec]
        newMD<-as.phyDat(newMD, levels=unique(unlist(newMD)), type="USER")
        write.phyDat(newMD, file=dataFile, nbcol=-1, append=T)
        write("", file=dataFile, sep="\n", append=T)
        
        #b2: randomize both for morpho - use 'newMD' from before which already has extant taxa
        # randomized and proportion of fossil characters removed
        b <- 2
        dataFile <- paste("Data_files/MorphoData_b", b, "_", nbTaxa, "Taxa_",i, "_", j, ".phy", sep = "")
        newMD <- ext.rand.MD
        if(propTaxRand[i] > 0) {
            fossil.indices <- grep("f", names(newMD))
            toRand <- sample(fossil.indices, size = propTaxRand[i] * nbFossils, replace = F)
            names(newMD)[toRand] <- sample(names(morphoDat)[toRand], length(toRand), replace = F)
        }

        newMD <- newMD[, tipLabelVec]
        newMD<-as.phyDat(newMD, levels=unique(unlist(newMD)), type="USER")
        write.phyDat(newMD, file=dataFile, nbcol=-1, append=T)
        write("", file=dataFile, sep="\n", append=T)        
    }
}
#------------------------------------------------------------------------------- 