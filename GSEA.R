GSEA = function(ranked, GO_file ) {
  set.seed(54321)
  library(gage)
  library(fgsea)
  
  rnks = ranked
  rnks =rnks[!is.na(names(rnks))]
  rnks = rnks[!duplicated(names(rnks))]
  rnks = sort(rnks, decreasing = TRUE)
  
  myGO = fgsea::gmtPathways(GO_file)
  
  fgseaRes <- fgsea::fgsea(pathways = myGO, 
                           stats = rnks,
                           minSize=15,
                           maxSize=400,
                           nperm=10000)
  
  fgseaRes = fgseaRes[fgseaRes$padj < 0.05, ]
  
  ## Filter FGSEA by using gage results. Must be significant and in same direction to keep 
  gobpres = gage(rnks, gsets=myGO, same.dir=TRUE, set.size =c(15,400))
  
  ups = as.data.frame(gobpres$greater) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < 0.05 ) %>%
    dplyr::select("Pathway")
  
  downs = as.data.frame(gobpres$less) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < 0.05 ) %>%
    dplyr::select("Pathway")
  
  tmpup = fgseaRes[fgseaRes$NES > 0 & !is.na(match(fgseaRes$pathway, ups$Pathway)), ]
  tmpdown = fgseaRes[fgseaRes$NES < 0 & !is.na(match(fgseaRes$pathway, downs$Pathway)), ]
  filtRes = rbind(tmpup, tmpdown)
  
  ### Collapse pathways
  keepsUp = fgsea::collapsePathways(filtRes[filtRes$NES > 1.5, ], myGO, rnks,  nperm = 500)
  keepsDown = fgsea::collapsePathways(filtRes[filtRes$NES < -1.5, ], myGO, rnks,  nperm = 500)
  keeps = c(keepsUp$mainPathways, keepsDown$mainPathways)
  filtRes = filtRes[ !is.na(match(filtRes$pathway, keeps)), ]
  
  filtResTidy = rbind(head(filtRes, n = 10), tail(filtRes, n = 10))
  filtResTidy$Enrichment = ifelse(filtResTidy$NES > 0, "Up-regulated", "Down-regulated")
  colvec =  c("blue3", "red3")
  names(colvec) = c("Down-regulated", "Up-regulated")
  
  g = ggplot(filtResTidy, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = Enrichment)) +
    scale_fill_manual(values= colvec) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="GSEA - Biological Processes (top 10)") + 
    theme_minimal()
  
  output = list("Results" = filtRes, "Plot" = g)
  return(output)
}
