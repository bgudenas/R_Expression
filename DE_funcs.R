## Brian Gudenas
## 6/7/19
## Functions to help DESeq2 workflow

PlotVolcano = function(tab, fc, pval, Name, intGenes=NULL, topn=7){
  #tab = data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$padj), Gene = res$geneName)
  #g = plotVolcano(tab, 1.5, 0.1, "Het vs WT", intGenes="CRX")
  library(dplyr)
  library(ggplot2)
  library(ggrepel)

lfc = log2(fc)

colvec = rep("grey", nrow(tab))
colvec[tab$negLogPval > -log10(pval) &  tab$logFC > lfc ] = "firebrick"
colvec[tab$negLogPval > -log10(pval) &  tab$logFC < -lfc ] = "dodgerblue"

g = ggplot(tab, aes(y=negLogPval, logFC) ) +
    geom_point(col=colvec, size = 2) +
    theme_bw() +
    xlab(expression(Log[2]~fold~change)) +
    ylab(expression(-Log[10]~pvalue)) +
    ggtitle(Name) +
    theme(text = element_text(size=20)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_text_repel(data= tab %>% 
                      filter(abs(logFC) > lfc & negLogPval > -log10(pval)) %>% 
                      group_by(logFC > 0) %>% 
                      top_n( topn, (  abs(logFC) ) ),
                    aes(label=Gene, size = 2) ) +
                    theme(legend.position="none")

if ( is.null(intGenes) == FALSE ){
  
  g = g + geom_text_repel(data= tab %>% 
                        filter(Gene == intGenes), 
                     aes(label=Gene, size = 2), fontface = "bold")
}

return(g)
}


countsWrap = function(geneName, Group){
  ID = rownames(res)[res$geneName == geneName]
  d = data.frame(count = counts[ ,colnames(counts) == ID] )
  d$Group = metadata %>%  dplyr::select(Group)
  
  pval = res$padj[res$geneName == geneName]
  if (pval < 0.05 ){
    pval = "P-value < 0.05"
  } else { pval = "" }
  
  g = ggplot(d, aes(x=Group , y=count, fill = Group)) +
    geom_boxplot( outlier.shape=NA ) +
    geom_point(position=position_jitter(w=0.1,h=0), size = 3) +
    ggtitle(paste( geneName, ID) ) +
    #   scale_fill_manual(values=c( "#E69F00", "#56B4E9")) +
    theme_bw() +
    th +
    annotate("text", x = 1, y = max(d$count)*.9, label = pval, size = 6 )
  
  return(g)
}


DEannotate = function(res, map){
  # res should be result from DESeq2 
  # map should be annotable from biomaRt
  
  res = as.data.frame( res[!is.na(res$padj), ] )
  mapmatch = match(rownames(res), map$ensembl_gene_id)
  res$geneName =  map$external_gene_name[ mapmatch]
  res$description = map$description[ mapmatch]
  res$biotype =  map$gene_biotype[ mapmatch]
  res$chromosome_name =  map$chromosome_name[ mapmatch]
  res$start_position =  map$start_position[ mapmatch]
  res$band = map$band[mapmatch]
  
  if ( grepl("ENSMUSG", rownames(res)[1] )==1 ) {
    res$hsapiens_gene = map$hsapiens_homolog_ensembl_gene[mapmatch]
    res$hsapiens_ensembl = map$hsapiens_homolog_ensembl_gene[mapmatch]
  } else if ( grepl("ENSG", rownames(res)[1] )==1 ) {
    res$mouse_gene = map$mmusculus_homolog_associated_gene_name[mapmatch]
    res$mouse_ensembl = map$mmusculus_homolog_ensembl_gene[mapmatch]
  }
  
  res = res[!is.na(res$geneName), ]
  res= as.data.frame(res)
}

