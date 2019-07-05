
library(phangorn)

#function calculating the proportion of internal nodes in the reftree which are recovered in a given tree
prop.recovered.reftree <- function(refTree, tree, rm.fossils = F) {
    
    if (rm.fossils) {
        fossils <- grep("f[0-9]*", refTree$tip.label)
        if (length(fossils) > 0 ) {
            refTree <- drop.tip(refTree, fossils)
        }
        fossils2 <- grep("f[0-9]*", tree$tip.label)
        if (length(fossils2) > 0 ) {
            tree <- drop.tip(tree, fossils2)
        }
    }
    
    if (length(setdiff(refTree$tip.label, tree$tip.label)) > 0) {
        warning("Tip labels do not match in tree comparison. Returning NA")
        return(NA)
    }
    
    nTips <- length(refTree$tip.label)
    parent<-refTree$edge[,1]
    child<-refTree$edge[,2]
    root<-as.integer(parent[!match(parent, child, 0)][1]) 
    tips<-child[is.na(!match(child,parent))]
    
    node.present <- rep(0, refTree$Nnode)
    for (n in 1:refTree$Nnode) {
        int.node <- n + nTips
        descs.ref <- refTree$tip.label[Descendants(refTree, int.node, type = "tips")[[1]]]
        MRCA <- getMRCA(tree, descs.ref)
        descs.tree <- tree$tip.label[Descendants(tree, MRCA, type = "tips")[[1]]]
        if ( length(setdiff(descs.tree, descs.ref)) == 0 ){
            node.present[n] <- 1 
        }
    }
    return(sum(node.present) / refTree$Nnode)
}

reftree <- read.tree("reftree.nwk")
contree <- read.nexus("con.tre")[[1]]
contree.woF <- read.nexus("con_woF.tre")[[1]]

prop.with.fossils <- prop.recovered.reftree(reftree, contree, rm.fossils = F)
prop.wo.fossils <- prop.recovered.reftree(reftree, contree.woF, rm.fossils = T)

write(paste("withFossils", "\t", prop.with.fossils, sep = ""), file = "nodes_correct.txt")
write(paste("woFossils", "\t", prop.wo.fossils, sep = ""), file = "nodes_correct.txt", append = T)

