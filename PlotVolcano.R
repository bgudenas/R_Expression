

PlotVolcano = function(tab, fc, pval, Name){
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  #tab =
 # tab = data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$padj), Gene = res$geneName)
  

lfc = log2(fc)

colvec = rep("grey", nrow(tab))
colvec[tab$negLogPval > -log10(pval) &  tab$logFC > lfc ] = "red3"
colvec[tab$negLogPval > -log10(pval) &  tab$logFC < -lfc ] = "blue3"

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
                      top_n( 10, ( negLogPval ) ),
                    aes(label=Gene, size = 1), fontface = "bold") +
                    theme(legend.position="none")

return(g)
}
# 
# g = plotVolcano(tab, 1.5, 0.1, "Het vs WT")
# g
