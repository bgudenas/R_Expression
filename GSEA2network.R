## Written by Brian Gudenas
## 11/12/2019
##--Dependencies
# fgsea
# WGCNA
# dynamicTreeCut
# stringr
# dplyr
# reshape2


GSEA2network = function(fgseaRes, GO, outname, minMod = 3 ){
  ## fgseaRes is the output from the fgsea R package
  ## GO is a list of lists, containing the GO_terms and their respective gene lists
  ## minMod is the minimum cluster size
  ## outname is the prefix appended to the output files in ./GSEAnet/, the output is a node file and an edge file which can be imported into cytoscape
  
  library(dplyr)
  fgseaRes$Length =0
  for (i in 1:nrow(fgseaRes)){
    fgseaRes$Length[i] = length(unlist(GO[names(myGO) == fgseaRes$pathway[i] ]))
  }
  fgseaRes = fgseaRes[ ,colnames(fgseaRes) != "leadingEdge"] ## remove this column if present since fgsea stores as list
  
  adjMat = matrix(nrow=nrow(fgseaRes), ncol = nrow(fgseaRes), data = 0)
  rownames(adjMat) = fgseaRes$pathway
  colnames(adjMat) = rownames(adjMat)

  ## Create matrix of jaccard coefficients
  for ( i in 1:nrow(fgseaRes)) {
    genesI =  unlist(GO[names(GO) == fgseaRes$pathway[i] ])
    for (j in 1:nrow(fgseaRes)) {
      genesJ = unlist(GO[names(myGO) == fgseaRes$pathway[j] ])
      overlap = sum(!is.na(match(genesI, genesJ )))
      jacc_coeff = overlap / length(unique(c(genesI, genesJ) ))  ## num overlaps / unique genes in sets
      adjMat[i,j] = jacc_coeff
    }
  }
  
  rownames(adjMat) = fgseaRes$pathway
  colnames(adjMat) = rownames(adjMat)
  
  fgseaRes$pathway = stringr::str_replace_all(tolower(stringr::str_replace(fgseaRes$pathway, "GO_", "" )), "_", " ")
  colnames(adjMat) =  stringr::str_replace_all(tolower(stringr::str_replace(colnames(adjMat), "GO_", "" )), "_", " ")
  rownames(adjMat) =  stringr::str_replace_all(tolower(stringr::str_replace(rownames(adjMat), "GO_", "" )), "_", " ")
   ## melt to longform
  adjLong = reshape2::melt(adjMat)

  ## Filter self connections
  adjLong =adjLong[ adjLong$Var1 !=  adjLong$Var2, ] 
  
## Filter duplicate connections with different order
  vec = vector(length = nrow(adjLong), mode = "list")
  for ( i in 1:nrow(adjLong)) {
    neword = sort(adjLong[i,c(1,2) ] )
    vec[i] = paste(neword[1,1], neword[1,2] , sep="", collapse = "")
  }
  
  adjLong$Filt =vec
  adjLong = adjLong[!duplicated(adjLong$Filt), ]
  adjLong = adjLong[ ,-4]
  
  ### Keep Every nodes highest Edge
  keeps = adjLong %>% 
    group_by(Var1) %>% 
    top_n(1, value) %>% 
    as.data.frame()
  
  ## Select top 20% of edges for visualization
  adjLong = adjLong[adjLong$value > quantile(adjLong$value, 0.80), ]
  
  keepsUpper = adjLong %>% 
    group_by(Var1) %>% 
    top_n(3, value) %>% 
    as.data.frame()
  
  adjLong = rbind(keeps, keepsUpper) %>% 
    distinct()
  
  adjLong = adjLong[ ,1:3]
  colnames(adjLong) = c("SourceNode", "TargetNode", "weight")
  #length(unique(adjLong$SourceNode))
  dir.create("./GSEAnet/", showWarnings = FALSE)
  write.table(adjLong, paste0("./GSEAnet/", outname, "_network.txt") , sep = "\t", row.names = FALSE, quote = FALSE)
  
  fgseaRes$padj = -log10(fgseaRes$padj)
  
  ## Cluster nodes using dynamic tree cut, for colors
  distMat = as.dist(1 - adjMat)  ## 1 - jaccards coef = jaccards distance
  dendro = hclust(distMat, method = "average" )
  clust = dynamicTreeCut::cutreeDynamicTree( dendro, minModuleSize = minMod,  deepSplit = TRUE)
  print(table(clust) )
  fgseaRes$Colors = WGCNA::labels2colors(clust)
  
  ## Select node with highest NES for label
  keeps =   fgseaRes %>% 
    dplyr::group_by(Colors) %>%
    dplyr::top_n( 1, abs(NES) )
  
  fgseaRes$RName = fgseaRes$pathway
  fgseaRes$RName[ is.na(match(fgseaRes$pathway, keeps$pathway))] = ""
  fgseaRes$Colors[fgseaRes$Colors == "black"] = "gold" ## avoid black nodes for plotting
  fgseaRes$NES = abs(fgseaRes$NES)
  
write.table(fgseaRes,  paste0("./GSEAnet/", outname, "_Nodes.txt") , sep = "\t", row.names = FALSE, quote = FALSE)
return(fgseaRes)
}

## EXAMPLE
# fgseaRes = readRDS("./Results/fgsea_DGE_PDX_cohort.rds") ## output from fgsea function
# myGO = fgsea::gmtPathways("../Annots/GSEA/c5.bp.v6.2.symbols.gmt" )  ##from MsigDB should match set used in fgsea call

## Note GSEA2network was designed to work with  up- or down-regulated nodes one set at a time

##------------ Down-reg
# Gnet = GSEA2network(fgseaRes =  fgseaRes[ fgseaRes$NES < 0 & fgseaRes$padj <= 0.05, ],
# GO = myGO, 
# outname = "DOWN_reg")

##------------ UP-reg
# Gnet = GSEA2network( fgseaRes = fgseaRes[ fgseaRes$NES > 0 & fgseaRes$padj <= 0.05, ],
#        GO = myGO, 
#        outname = "UP_reg")
