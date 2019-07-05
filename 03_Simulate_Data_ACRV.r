
library(phangorn)
#-----------------------------------
#settings
nbReps <- 20
treeAge <- 100
nbTaxa <- 12
nbFossils <- nbTaxa/3
nbSites <- 1000
nbMorpho <- 100
propMorpho <- 0.5
molecRate <- 0.25/treeAge
TreeFolder <- "../../1_Priors_TreeAge_CR/a_MakeTrees-D-simulateSeqs/"
alphaVector <- c(1000000, 10, 1, 0.5, 0.1, 0.01)
tipLabelVec <- c(paste("t", 1:nbTaxa, sep = ""), paste("f", 1:nbFossils, sep = ""))
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
#(2) simulate morpho datasets on each tree, under different rate regimes

# A - discrete gamma rates, four categories
nbCats <- 4
for (t in 1:length(trees)) {
    tree <- trees[[t]]
    dataFile <- paste("Data_files/MorphoData_", nbTaxa, "Taxa_GammaDistr_", t, ".phy", sep = "")
    for (i in 1 : (length(alphaVector)) ){
        morphoDat <- NULL
        gammaRates <- discrete.gamma(alphaVector[i], nbCats)
        for (Cat in 1:nbCats){
            nbVarSites <- 0
            morphoDatTemp <- NULL  
            currRate <- molecRate * gammaRates[Cat]
            if (currRate < 0.00001){
                currRate <- 0.00001
            }
            while (nbVarSites < nbMorpho / nbCats){
                mD <- simSeq(tree, l = nbMorpho / nbCats, type = "DNA", bf=c(0.5, 0.5, 0, 0), rate = currRate)
                mD <- as.data.frame(mD)
                mD[mD == "c"] <- "0"
                mD[mD == "a"] <- "1"
                isVariable <- rep(F, dim(mD)[1])
                for ( s in 1:dim(mD)[1] ){
                    if (length(unique(as.numeric(mD[s,]))) == 2){
                        isVariable[s] <- T
                    }
                }
                if (length(morphoDatTemp) == 0 ){
                    morphoDatTemp <- mD[isVariable,]
                } else {
                    morphoDatTemp <- rbind(morphoDatTemp, mD[isVariable,])
                }
                nbVarSites <- nbVarSites + sum(isVariable)
            }
            morphoDatTemp <- morphoDatTemp[1:(nbMorpho / nbCats),]
            morphoDat <- rbind(morphoDat, morphoDatTemp)
        }
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

# B - two rate categories to mimic a bimodal distribution. Based as before on the discretized Gamma:
# rate 1 equals the lowest of the 4-category discrete gamma and contains half of the sites
# rate 2 equals the highest rate and contains the other half of the sites

nbCats <- 4
for (t in 1:length(trees)) {
  tree <- trees[[t]]
  dataFile <- paste("Data_files/MorphoData_", nbTaxa, "Taxa_BimodalDistr_", t, ".phy", sep = "")
  for (i in 1 : (length(alphaVector)) ){
    morphoDat <- NULL
    gammaRates <- discrete.gamma(alphaVector[i], nbCats)
    #lower rate, half the sites  
    nbVarSites <- 0
    morphoDatTemp <- NULL
    currRate <- molecRate * gammaRates[1]
    if (currRate < 0.00001){
        currRate <- 0.00001
    }
    while (nbVarSites < nbMorpho / 2){
        mD <- simSeq(tree, l = nbMorpho / 2, type = "DNA", bf=c(0.5, 0.5, 0, 0), rate = currRate)
        mD <- as.data.frame(mD)
        mD[mD == "c"] <- "0"
        mD[mD == "a"] <- "1"
        isVariable <- rep(F, dim(mD)[1])
        for ( s in 1:dim(mD)[1] ){
          if (length(unique(as.numeric(mD[s,]))) == 2){
            isVariable[s] <- T
          }
        }
        if (length(morphoDatTemp) == 0 ){
          morphoDatTemp <- mD[isVariable,]
        } else {
          morphoDatTemp <- rbind(morphoDatTemp, mD[isVariable,])
        }
        nbVarSites <- nbVarSites + sum(isVariable)
    }
    morphoDatTemp <- morphoDatTemp[1:(nbMorpho / 2),]
    morphoDat <- morphoDatTemp
    
    #higher rate, again half the sites  
    nbVarSites <- 0
    morphoDatTemp <- NULL
    higherRate <- molecRate * gammaRates[nbCats]
    if (higherRate < 0.00001){
      higherRate <- 0.00001
    }
    while (nbVarSites < nbMorpho / 2){
        mD <- simSeq(tree, l = nbMorpho / 2, type = "DNA", bf=c(0.5, 0.5, 0, 0), rate = higherRate)
        mD <- as.data.frame(mD)
        mD[mD == "c"] <- "0"
        mD[mD == "a"] <- "1"
        isVariable <- rep(F, dim(mD)[1])
        for ( s in 1:dim(mD)[1] ){
          if (length(unique(as.numeric(mD[s,]))) == 2){
            isVariable[s] <- T
          }
        }
        if (length(morphoDatTemp) == 0 ){
          morphoDatTemp <- mD[isVariable,]
        } else {
          morphoDatTemp <- rbind(morphoDatTemp, mD[isVariable,])
        }
        nbVarSites <- nbVarSites + sum(isVariable)
    }
    morphoDatTemp <- morphoDatTemp[1:(nbMorpho / 2),]
    morphoDat <- rbind(morphoDat, morphoDatTemp)
    
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