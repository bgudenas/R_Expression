GSEA = function(ranked, GO_file, pval = 0.05, display_num = 10, UP_only = FALSE) {
  set.seed(54321)
  library(gage)
  library(fgsea)
  
  rnks = ranked
  rnks =rnks[!is.na(names(rnks))]
  rnks = rnks[!duplicated(names(rnks))]
  rnks = sort(rnks, decreasing = TRUE)
  
  if (is.character(GO_file)){
    myGO = fgsea::gmtPathways(GO_file)
  } else {
    myGO = GO_file 
  }
  
  fgseaRes <- fgsea::fgsea(pathways = myGO, 
                           stats = rnks,
                           minSize=15,
                           maxSize=400,
                           nperm=100000)
  
  fgseaRes = fgseaRes[fgseaRes$padj <= 0.05, ]
  
  ## Filter FGSEA by using gage results. Must be significant and in same direction to keep 
  gobpres = gage(rnks, gsets=myGO, same.dir=TRUE, set.size =c(15,400))
  
  ups = as.data.frame(gobpres$greater) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val <= 0.05 ) %>%
    dplyr::select("Pathway")
  
  downs = as.data.frame(gobpres$less) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val <= 0.05 ) %>%
    dplyr::select("Pathway")
  
  tmpup = fgseaRes[fgseaRes$NES > 0 & !is.na(match(fgseaRes$pathway, ups$Pathway)), ]
  tmpdown = fgseaRes[fgseaRes$NES < 0 & !is.na(match(fgseaRes$pathway, downs$Pathway)), ]
  filtRes = rbind(tmpup, tmpdown)
  
  ### Collapse pathways
  keepsUp = fgsea::collapsePathways(filtRes[filtRes$NES > 1.5, ], myGO, rnks,  nperm = 500)
  keepsDown = fgsea::collapsePathways(filtRes[filtRes$NES < -1.5, ], myGO, rnks,  nperm = 500)
  keeps = c(keepsUp$mainPathways, keepsDown$mainPathways)
  filtRes = filtRes[ !is.na(match(filtRes$pathway, keeps)), ] %>% 
    arrange(desc(NES))
  
  ## if all same sign select top 10 or as many as available if split sign select top 10 for each or top pos,neg available
  if (  all(filtRes$NES > 0 ) | all(filtRes$NES < 0 ) ){
    
    sel_sets = ifelse( nrow(filtRes) >= display_num, display_num, nrow(filtRes))
   filtResTidy  =  filtRes %>%
     arrange(desc(abs(NES))) %>% 
     top_n(n = sel_sets, wt = abs(NES))
  } else {
    NES_logvec = filtRes$NES > 0
    sel_pos = ifelse( nrow(filtRes[NES_logvec, ]) >= display_num, display_num, nrow(filtRes[NES_logvec, ]))
    sel_neg = ifelse( nrow(filtRes[!NES_logvec, ]) >= display_num, display_num, nrow(filtRes[!NES_logvec, ]))
    
    filtResTidy  = rbind(filtRes[1:sel_pos, ],
                         filtRes[(nrow(filtRes) - sel_neg + 1 ):nrow(filtRes),  ])
  }
  
  filtResTidy$Enrichment = ifelse(filtResTidy$NES > 0, "Up-regulated", "Down-regulated")
  if (UP_only == TRUE ){ filtResTidy = filtResTidy[ filtResTidy$NES > 0, ] }
  #colvec = setNames(object = c("blue3", "red3"), nm = c("Down-regulated", "Up-regulated"))
  upcols =  colorRampPalette(colors = c("red4", "red1", "lightpink"))( sum(filtResTidy$Enrichment == "Up-regulated"))
  downcols =  colorRampPalette(colors = c( "lightblue", "blue1", "blue4"))( sum(filtResTidy$Enrichment == "Down-regulated"))
  colos = c(upcols, downcols)
  names(colos) = 1:length(colos)
  filtResTidy$Index = as.factor(1:nrow(filtResTidy))
  filtResTidy$pathway = str_remove(filtResTidy$pathway, "GO_")
  filtResTidy$pathway = str_replace_all(filtResTidy$pathway, "_", " ")
  
  g = ggplot(filtResTidy,  aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = Index)) +
    scale_fill_manual(values= colos ) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="GSEA - Biological Processes") + 
    theme_minimal() +
    theme(legend.position = "none")
  
  output = list("Results" = filtRes, "Plot" = g)
  return(output)
}




### just fgsea

GSEA_single = function(ranked, GO_file, pval = 0.05, display_num = 10, UP_only = FALSE) {
  set.seed(54321)
  library(gage)
  library(fgsea)
  library(dplyr)
  library(stringr)
  
  rnks = ranked
  rnks =rnks[!is.na(names(rnks))]
  rnks = rnks[!duplicated(names(rnks))]
  rnks = sort(rnks, decreasing = TRUE)
  
  myGO = fgsea::gmtPathways(GO_file)
  
  fgseaRes <- fgsea::fgsea(pathways = myGO, 
                           stats = rnks,
                           minSize=15,
                           maxSize=400,
                           nperm=100000)
  
  filtRes = fgseaRes[fgseaRes$padj <= 0.05, ]

  
  ### Collapse pathways
  keepsUp = fgsea::collapsePathways(filtRes[filtRes$NES > 1.5, ], myGO, rnks,  nperm = 500)
  keepsDown = fgsea::collapsePathways(filtRes[filtRes$NES < -1.5, ], myGO, rnks,  nperm = 500)
  keeps = c(keepsUp$mainPathways, keepsDown$mainPathways)
  filtRes = filtRes[ !is.na(match(filtRes$pathway, keeps)), ] %>% 
    arrange(desc(NES))
  
  ## if all same sign select top 10 or as many as available if split sign select top 10 for each or top pos,neg available
  if (  all(filtRes$NES > 0 ) | all(filtRes$NES < 0 ) ){
    
    sel_sets = ifelse( nrow(filtRes) >= display_num, display_num, nrow(filtRes))
    filtResTidy  =  filtRes %>%
      arrange(desc(abs(NES))) %>% 
      top_n(n = sel_sets, wt = abs(NES))
  } else {
    NES_logvec = filtRes$NES > 0
    sel_pos = ifelse( nrow(filtRes[NES_logvec, ]) >= display_num, display_num, nrow(filtRes[NES_logvec, ]))
    sel_neg = ifelse( nrow(filtRes[!NES_logvec, ]) >= display_num, display_num, nrow(filtRes[!NES_logvec, ]))
    
    filtResTidy  = rbind(filtRes[1:sel_pos, ],
                         filtRes[(nrow(filtRes) - sel_neg + 1 ):nrow(filtRes),  ])
  }
  
  filtResTidy$Enrichment = ifelse(filtResTidy$NES > 0, "Up-regulated", "Down-regulated")
  if (UP_only == TRUE ){ filtResTidy = filtResTidy[ filtResTidy$NES > 0, ] }
  #colvec = setNames(object = c("blue3", "red3"), nm = c("Down-regulated", "Up-regulated"))
  upcols =  colorRampPalette(colors = c("red4", "red1", "lightpink"))( sum(filtResTidy$Enrichment == "Up-regulated"))
  downcols =  colorRampPalette(colors = c( "lightblue", "blue1", "blue4"))( sum(filtResTidy$Enrichment == "Down-regulated"))
  colos = c(upcols, downcols)
  names(colos) = 1:length(colos)
  filtResTidy$Index = as.factor(1:nrow(filtResTidy))
  filtResTidy$pathway = str_remove(filtResTidy$pathway, "GO_")
  filtResTidy$pathway = str_replace_all(filtResTidy$pathway, "_", " ")
  
  g = ggplot(filtResTidy,  aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = Index)) +
    scale_fill_manual(values= colos ) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="GSEA - Biological Processes") + 
    theme_minimal() +
    theme(legend.position = "none")
  
  output = list("Results" = filtRes, "Plot" = g)
  return(output)
}
