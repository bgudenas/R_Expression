# resOG = res
# res = as.data.frame(res)
# 
# chromosome = "X"
# rchrom = rchrom[rchrom$baseMean > 20, ]

chromPlot = function(res, chromosome ){
  library(dplyr)
  library(ggplot2)
  th <- theme(text = element_text(size=20, family = "Helvetica" ), panel.grid.major = element_blank(), panel.grid.minor = element_blank() )
  
  rchrom = res %>% 
    filter(chromosome_name == chromosome & baseMean > 20 ) %>% 
    arrange(start_position) %>% 
    mutate(Start = rank(start_position)) %>% 
    mutate(up = ifelse(log2FoldChange > 0.58 & padj < 0.05 , 1, 0)) %>% 
    mutate(down = ifelse(log2FoldChange < -0.58 & padj < 0.05 , -1, 0)) 
  
  rchrom$start_position = rchrom$start_position/1e+06
  rchrom$CNV = as.factor(rchrom$up + rchrom$down)
  
  df = windowSlide(rchrom$start_position, rchrom$log2FoldChange, rchrom$padj)
  ups = as.data.frame(df[1])
  downs = as.data.frame(df[2])
  
  g = ggplot(rchrom, aes(y = log2FoldChange, x = start_position )) +
    geom_point( aes(color = CNV), size = 1) +
    scale_color_manual(values = c('-1' = "dodgerblue",
                         '0'= "grey",
                         '1' = "red3")) +
    geom_hline( yintercept =  0 , col = "black", lwd =1 , lty =2, alpha = 0.5)  +
    theme_bw()+
    theme(legend.position = "none" ) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank()) +
    annotate(geom = "text", label =chromosome, x = 1, y = 0.8*(min(rchrom$log2FoldChange)), size =6)
  
  #   geom_rect(data=ups, inherit.aes=FALSE, aes(xmin=starts, xmax=ends, ymin=min(rchrom$log2FoldChange),
   #      ymax=max(rchrom$log2FoldChange), group=group), color="transparent", fill="red2", alpha=0.2) 
     # geom_rect(data=downs, inherit.aes=FALSE, aes(xmin=starts, xmax=ends, ymin=min(rchrom$log2FoldChange),
     #                                            ymax=max(rchrom$log2FoldChange), group=group), color="transparent", fill="blue2", alpha=0.2)
    if (nrow(ups) > 0 ){ 
      g = g +geom_rect(data=ups, inherit.aes=FALSE, aes(xmin=starts, xmax=ends, ymin=min(rchrom$log2FoldChange),
ymax=max(rchrom$log2FoldChange), group=group), color="transparent", fill="red2", alpha=0.3) 
    } 
  if (nrow(downs) > 0 ){ 
    g = g +geom_rect(data=downs, inherit.aes=FALSE, aes(xmin=starts, xmax=ends, ymin=min(rchrom$log2FoldChange),
                                                      ymax=max(rchrom$log2FoldChange), group=group), color="transparent", fill="blue2", alpha=0.3) 
  }      
        
    return(g)
}
   
  
  # 
  # starts = rchrom$start_position
  # L2FC = rchrom$log2FoldChange
slideMerge = function(df ){
  ## merges continuos ranges
  for ( j in 1:10) {
    i=1
    while(i < nrow(df) ) {
      end1 = df$ends[i]
      begin2 = df$starts[i + 1]
      
      if (end1 > begin2){ 
        df$ends[i] = df$ends[i + 1]
        df = df[-c(i+1), ]
      }
      i =i+1
    }
  }
  return(df)
}
  
  
  windowSlide = function(starts, L2FC, pvals, window = 1.5, slide = 0.5) {
    begin=0
    winscores = c()
    beginvec = c()
    endvec = c()
    for ( i in 1:round(max(starts)/slide) ){
      end = begin + window 
      beginvec = c(beginvec, begin )
      endvec = c(endvec, end)
      
      keeps = (starts >= begin ) & ( starts <= end )
      score = median( L2FC[keeps] )
      # if (!is.na(avepval)) {
      #   if (avepval > 0.1 ){ score = 0 }
      # }
      if ( sum( pvals[keeps] < 0.05) == 0 ) { score = 0}
      if (is.na(score) ){score = 0}
      if ( sum(keeps) < 3 ){score = 0}
      winscores = c(winscores, score )
      begin = begin + slide
     
    }
    
    df = data.frame(starts = beginvec, ends = endvec, scores = winscores)
    df$col = ""
    df$col[df$scores >= 1 ] = "Gain"
    df$col[df$scores <= -1 ] = "Loss"
    ups = df[df$col == "Gain", ]
    ups$group = seq_along(ups$starts)
    ups = slideMerge(ups)
    
    downs = df[df$col == "Loss", ]
    downs$group = seq_along(downs$starts)
    downs = slideMerge(downs)
    mega = list(ups, downs)
    return(mega)
  }
#df = windowSlide(starts, L2FC, pvals)      


plotAllChroms = function(res){
  chrs = as.character( c(1:22, "X"))
pvec = vector(length = 23, "list")
names(pvec) = chrs
for ( i in 1:length(chrs)){
  chrom = chrs[i]
  pvec[[i]] = chromPlot(as.data.frame(res), chrom )
}
mega = gridExtra::grid.arrange(nrow = 6, pvec[[1]], pvec[[2]],pvec[[3]], pvec[[4]], pvec[[5]],
 pvec[[6]], pvec[[7]], pvec[[8]], pvec[[9]],
pvec[[10]],pvec[[11]],pvec[[12]],pvec[[13]],pvec[[14]],
pvec[[15]],pvec[[16]],pvec[[17]],pvec[[18]],pvec[[19]],pvec[[20]],
pvec[[21]],pvec[[22]],pvec[[23]] )
return(mega)
}

#plotAllChroms(res)