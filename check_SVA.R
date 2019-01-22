## Func to do SVA on each cohort
#library(DESeq2)

    SVA_check = function(countMat, annots, Batch) {
      ##Assumes countmat is ordered to annots
      library(DESeq2)
      library(sva)
      
      countMat = countMat[  , annots$Batch == Batch]
      metadata = annots[annots$Batch == Batch, ]      

            DEmat =DESeqDataSetFromMatrix(countData = countMat, colData = metadata, design = ~ ELP1 ) 
            DE  = DESeq(DEmat)
            
      countNorm = counts(DE, normalized = TRUE)
      countFilt = countNorm[ rowSums(countNorm > 10) > 2 , ]
      
      mod = model.matrix(~ ELP1 , data =  metadata )
      mod0 = model.matrix(~ 1, data= metadata )
      
      n.sv = num.sv(countFilt, mod, method="leek")
      svseq = svaseq(countFilt , mod, mod0, n.sv = 1, B = 10)
      rownames(svseq$sv) = DE$IDs
      
      par(mfrow = c(1, 1), mar = c(3,5,3,1))
      for (i in 1:1) {
        stripchart(svseq$sv[, i] ~ DE$ELP1, vertical = TRUE, main = paste0("SV", i))
        abline(h = 0)
      }
      return(list(svseq$sv, DE$ELP1 ) )
    }
  
   # vals =  SVA_check(mega, annots, "SJYC07")
   # stripchart(vals[[1]] ~ vals[[2]], vertical = TRUE, main = paste0("SV", i))
   