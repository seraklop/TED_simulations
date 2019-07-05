
library(phangorn)
#-----------------------------------
#settings
nbReps <- 10
treeAge <- 100
nbTaxa <- 48
nbFossils <- nbTaxa/3
nbSites <- 1000
nbMorpho <- 100
propMorpho <- 0.5
molecRate <- 0.25/treeAge
TreeFolder <- "../../1_Priors_TreeAge_CR/a_MakeTrees-D-simulateSeqs/"

propMorpho <- c( 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.99 )
Rate <- 0.25/treeAge

tipLabelVec <- c(paste("t", 1:nbTaxa, sep = ""), paste("f", 1:nbFossils, sep = ""))
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
#b0: simulate morpho datasets on each tree, under different sampling proportions for the fossils (morpho chars sampled randomly)
for (t in 1:length(trees)) {
    tree <- trees[[t]]
    dataFile <- paste("Data_files/b0_MorphoData_Random_", t, ".phy", sep = "")
    for (i in 1 : (length(propMorpho)) ){
            nbVarSites <- 0
            morphoDat <- NULL
            while (nbVarSites < nbMorpho ){
                mD <- simSeq(tree, l = nbMorpho, type = "USER", Q = cbind(c(0,1),c(1,0)), levels = c(0,1), rate = Rate)
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
            fossil.indices <- grep("f", names(morphoDat))
            for (f in 1:length(fossil.indices)){
                morphoDat[sample(1:nbMorpho, size = nbMorpho*(1-propMorpho[i]), replace = F ), fossil.indices[f]] <- "?"
            }

            morphoDat <- morphoDat[, tipLabelVec]
            morphoDat<-as.phyDat(morphoDat, levels=c("1","0"), type="USER")
            write.phyDat(morphoDat, file=dataFile, nbcol=-1, append=T)
            write("", file=dataFile, sep="\n", append=T)   
    }
}
#------------------------------------------------------------------------------- 
#b1  different sampling proportions for the fossils with half non-random, half random
for (t in 1:length(trees)) {
    tree <- trees[[t]]
    dataFile <- paste("Data_files/b1_MorphoData_HalfRandom_", t, ".phy", sep = "")
    for (i in 1 : (length(propMorpho)) ){
        nbVarSites <- 0
        morphoDat <- NULL
        while (nbVarSites < nbMorpho ){
            mD <- simSeq(tree, l = nbMorpho, type = "USER", Q = cbind(c(0,1),c(1,0)), levels = c(0,1), rate = Rate)
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
        
        #remove part of the morpho data for the fossils - half randomly, half always the same
        nbNonrandom <- ceiling( 0.5 *nbMorpho *(1-propMorpho[i]) )
        nbRandom <- floor( 0.5 *nbMorpho *(1-propMorpho[i]) )
        
        charsAllMissing <- sample(1:nbMorpho, size = nbNonrandom, replace = F )
        charsMaybeNotMissing <- setdiff(seq(from = 1, to = nbMorpho, by = 1), charsAllMissing)
        
        for ( fossil in 1: nbFossils ){
            fossilID <- which(names(morphoDat) == paste("f", fossil, sep = "") )
            morphoDat[[fossilID]][charsAllMissing] <- "?"
            charsRmRandom <- sample(charsMaybeNotMissing, size = nbRandom, replace = F )
            morphoDat[[fossilID]][charsRmRandom] <- "?"
        }
        
        morphoDat <- morphoDat[, tipLabelVec]
        morphoDat<-as.phyDat(morphoDat, levels=c("1","0"), type="USER")
        write.phyDat(morphoDat, file=dataFile, nbcol=-1, append=T)
        write("", file=dataFile, sep="\n", append=T)   
         
    }
}
#------------------------------------------------------------------------------- 
#b2: different sampling proportions for the fossils with all non-random
for (t in 1:length(trees)) {
    tree <- trees[[t]]
    dataFile <- paste("Data_files/b2_MorphoData_AllNonRandom_", t, ".phy", sep = "")
    for (i in 1 : (length(propMorpho)) ){
        nbVarSites <- 0
        morphoDat <- NULL
        while (nbVarSites < nbMorpho ){
            mD <- simSeq(tree, l = nbMorpho, type = "USER", Q = cbind(c(0,1),c(1,0)), levels = c(0,1), rate = Rate)
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
        
        #remove part of the morpho chars for all the fossils
        nbNonrandom <- nbMorpho *(1-propMorpho[i])
        charsAllMissing <- sample(1:nbMorpho, size = nbNonrandom, replace = F )
        
        for ( fossil in 1: nbFossils ){
            fossilID <- which(names(morphoDat) == paste("f", fossil, sep = "") )
            morphoDat[[fossilID]][charsAllMissing] <- "?"
        }
        
        morphoDat <- morphoDat[, tipLabelVec]
        morphoDat<-as.phyDat(morphoDat, levels=c("1","0"), type="USER")
        write.phyDat(morphoDat, file=dataFile, nbcol=-1, append=T)
        write("", file=dataFile, sep="\n", append=T)   
         
    }
}
#------------------------------------------------------------------------------- 
