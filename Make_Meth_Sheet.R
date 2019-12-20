
# 
 #dir = "./Methylation/IDATS/ETMRlikeDICER/"
 #outname = "./Methylation/IDATS/ETMRlikeDICER/newRefs_meta.csv"
 #group = "450k"

makeMethSheet = function(dir, group, outname){
  #dir = "./Methylation/IDATS/ETMRlikeDICER/"
  #outname = "./Methylation/IDATS/ETMRlikeDICER/newRefs_meta.csv"
  #group = "450k"
  idats = list.files(dir)
  idats = idats[grepl("idat", idats)]
  
  isplit = stringr::str_split(idats, "_")
  samps = paste0(unlist(lapply(isplit, "[[", 1)), "_",unlist(lapply(isplit, "[[", 2)))
  
  df = data.frame("Sample_Name" = samps, "Sample_Group" = group, status = "cancer", 
                  "Slide" = unlist(lapply(isplit, "[[", 1)),
                  "Array" = unlist(lapply(isplit, "[[", 2)))
  
  df =df[!duplicated(df$Sample_Name), ]
  
  write.csv(df, outname, quote = FALSE, row.names = FALSE)
}

#makeMethSheet(dir, group, outname)  
