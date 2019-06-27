
GSEA2network = function(fgseaRes, GO, outname, minMod = 2, grey = TRUE){
  library(dplyr)
  fgseaRes$Length =0
  for (i in 1:nrow(fgseaRes)){
    fgseaRes$Length[i] = length(unlist(GO[names(myGO) == fgseaRes$pathway[i] ]))
  }
  
  df = matrix(nrow=nrow(fgseaRes), ncol = nrow(fgseaRes), data = 0)
  rownames(df) = fgseaRes$pathway
  colnames(df) = rownames(df)

  for ( i in 1:nrow(fgseaRes)) {
    genesI =  unlist(GO[names(GO) == fgseaRes$pathway[i] ])
    for (j in 1:nrow(fgseaRes)) {
      genesJ = unlist(GO[names(myGO) == fgseaRes$pathway[j] ])
      overlap = sum(!is.na(match(genesI, genesJ )))
      shared = overlap / length(unique(c(genesI, genesJ) ))
      df[i,j] = shared
    }
  }
  # GOmat = matrix(nrow = length(genes), ncol = nrow(fgseaRes), data  = 0)
  # rownames(GOmat) = genes
  # colnames(GOmat) = fgseaRes$pathway
  # 
  # ## Create Gene x GOterm binary matrix
  # for (i in 1:nrow(fgseaRes)){
  #   GOmat[ ,i] = ifelse(is.na(match(genes, unlist(myGO[names(myGO) == fgseaRes$pathway[i] ]))), 0, 1)
  # }
  # 
  # GOsim = outer(1:ncol(GOmat), 1:ncol(GOmat), FUN = Vectorize( function(i,j) irr::kappa2( cbind(GOmat[ ,i], GOmat[ ,j])) [[5]]  ) )
  # df = GOsim
  rownames(df) = fgseaRes$pathway
  colnames(df) = rownames(df)
  
  fgseaRes$pathway = stringr::str_replace_all(tolower(stringr::str_replace(fgseaRes$pathway, "GO_", "" )), "_", " ")
  colnames(df) =  stringr::str_replace_all(tolower(stringr::str_replace(colnames(df), "GO_", "" )), "_", " ")
  rownames(df) =  stringr::str_replace_all(tolower(stringr::str_replace(rownames(df), "GO_", "" )), "_", " ")
  
  df2 = reshape2::melt(df)
 # df2 =df2[df2$value > 0, ]
  ## Filter self connections
  df2 =df2[ df2$Var1 !=  df2$Var2, ] 
  
## FIlter duplicate connections with different order
  vec = vector(length = nrow(df2), mode = "list")
  for ( i in 1:nrow(df2)) {
    neword = sort(df2[i,c(1,2) ] )
    vec[i] = paste(neword[1,1], neword[1,2] , sep="", collapse = "")
  }
  
  df2$Filt =vec
  df2 = df2[!duplicated(df2$Filt), ]
  df2 = df2[ ,-4]
  
  ### Keep Every nodes highest Edge
  keeps = df2 %>% 
    group_by(Var1) %>% 
    top_n(1, value) %>% 
    as.data.frame()
  
  ## Select top 10% of edges for visualization
  df2 = df2[df2$value > quantile(df2$value, 0.80), ]
  
  keepsUpper = df2 %>% 
    group_by(Var1) %>% 
    top_n(3, value) %>% 
    as.data.frame()
  
  df2 = rbind(keeps, keepsUpper) %>% 
    distinct()
  
  df2 = df2[ ,1:3]
  colnames(df2) = c("SourceNode", "TargetNode", "weight")
  length(unique(df2$SourceNode))
  write.table(df2, paste0("./Results/GSEAnet/", outname, "_network.txt") , sep = "\t", row.names = FALSE, quote = FALSE)
  
  #fgseaRes = fgseaRes[ ,c(1,3,5,9) ]
  fgseaRes$padj = -log10(fgseaRes$padj)
  
  ## Cluster nodes using dynamic tree cut, for colors
  distMat = as.dist(1-df)
  dendro = hclust(distMat, method = "average" )
  clust = dynamicTreeCut::cutreeDynamicTree( dendro, minModuleSize = minMod,  deepSplit = TRUE)
  print(table(clust) )
  fgseaRes$Colors = WGCNA::labels2colors(clust)
  
  
  keeps =   fgseaRes %>% 
    dplyr::group_by(Colors) %>%
    dplyr::top_n( 1, abs(NES) )
  
  fgseaRes$RName = fgseaRes$pathway
  fgseaRes$RName[ is.na(match(fgseaRes$pathway, keeps$pathway))] = ""
  fgseaRes$Colors[fgseaRes$Colors == "black"] = "gold"
  fgseaRes$NES = abs(fgseaRes$NES)
  
write.table(fgseaRes,  paste0("./Results/GSEAnet/", outname, "_Nodes.txt") , sep = "\t", row.names = FALSE, quote = FALSE)
return(fgseaRes)
}

##
# fgseaRes = readRDS("./Results/fgsea_DGE_ELP1cohort.rds")
# myGO = fgsea::gmtPathways("../Annots/GSEA/c5.bp.v6.2.symbols.gmt" )
# Gnet = GSEA2network(fgseaRes[ fgseaRes$NES < 0 & fgseaRes$padj <= 0.01, ], myGO, "DOWNreg")
# table(Gnet$Colors)
# Gnet[Gnet$Colors == "grey", ]
# Gnet[grepl("stress", Gnet$pathway), ]
# Gnet[grepl("protein", Gnet$pathway), ]
# Gnet[Gnet$Colors == "brown", ]

#Gnet = GSEA2network(fgseaRes, myGO, "UPreg" )
