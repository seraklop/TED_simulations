
library(phangorn)
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
TreeFolder <- "../../1_Priors_TreeAge_CR/a_MakeTrees-D-simulateSeqs/"
RateVector <- c(0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 2.0, 10.0)
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
#(2)simulate morpho datasets on each tree, under different rates

for (t in 1:length(trees)) {
    tree <- trees[[t]]
    dataFile <- paste("Data_files/MorphoData_", nbTaxa, "Taxa_", t, ".phy", sep = "")
    for (i in 1 : (length(RateVector)) ){
        nbParsInfSites <- 0
        morphoDat <- NULL
        while (nbParsInfSites < nbMorpho ){
            mD <- simSeq(tree, l = nbMorpho, type = "DNA", bf=c(0.5, 0.5, 0, 0), rate = RateVector[i]/treeAge)
            mD <- as.data.frame(mD)
            mD[mD == "c"] <- "0"
            mD[mD == "a"] <- "1"
            isParsInf <- rep(F, dim(mD)[1])
            for ( s in 1:dim(mD)[1] ){
                nb0s <- length(which(mD[s,] == "0"))
                nb1s <- length(which(mD[s,] == "1"))
                if ( nb0s >= 2 && nb1s >=2 ){
                    isParsInf[s] <- T
                }
            }
            if (length(morphoDat) == 0 ){
                morphoDat <- mD[isParsInf,]
            } else {
                morphoDat <- rbind(morphoDat, mD[isParsInf,])
            }
            nbParsInfSites <- nbParsInfSites + sum(isParsInf)
        }
        morphoDat <- morphoDat[1:nbMorpho,]
        
        #remove part of the morpho data for the fossils
        fossil.indices <- grep("f", names(morphoDat))
        for (f in 1:length(fossil.indices)){
            morphoDat[sample(1:nbMorpho, size = nbMorpho*(1-propMorpho), replace = F ), fossil.indices[f]] <- "?"
        }

        morphoDat <- morphoDat[, tipLabelVec]
        morphoDat <- as.phyDat(morphoDat, levels=c("1","0"), type="USER")
        write.phyDat(morphoDat, file=dataFile, nbcol=-1, append=T)
        write("", file=dataFile, sep="\n", append=T)   
          
    }
}
