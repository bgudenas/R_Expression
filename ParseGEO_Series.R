
# pathMBs = "../Datasets/Bulk_Expression/Microarray/Human/GSE37382_MB_1_1_ST/GSE37382_series_matrix.txt"
# MBs = ParseGEO_Microarray(series_matrix_path = pathMBs,
#                           platform = "affy_hugene_1_0_st_v1", 
#                           species_dataset = "hsapiens_gene_ensembl",
#                           outname = "GSE37382_MB_1_1_ST_Proccesed"
# )


ParseGEO_Microarray = function(series_matrix_path, platform = NULL, species_dataset = NULL, outname ) {
  options(stringsAsFactors = FALSE)
  all = read.table(series_matrix_path, fill = TRUE, sep = "\t",header = FALSE )
  summary = as.character(all[ which(all[ ,1] == "!Series_summary"), 2 ])
  design =  as.character(all[ which(all[ ,1] == "!Series_overall_design"), 2 ])
  markMeta = which(all[,1] == "!Sample_title") ## Line metadata tbale begins
  meta = read.table(series_matrix_path, fill = TRUE, sep = "\t", skip = markMeta, header = FALSE )
  markExpr = which(meta[ ,1] == "ID_REF") + markMeta - 1 ## Line expression matrix begins
  datExpr = read.table(series_matrix_path, fill = TRUE, sep = "\t", skip = markExpr, header = TRUE, row.names = 1) 
  
  which(meta[ ,1] == "!Sample_characteristics_ch1" )
  samp_chars = c()
  for ( i in  which(meta[ ,1] == "!Sample_characteristics_ch1" ) ) {
    
   samp_chars = cbind(samp_chars, as.character(unlist(meta[ i, ])) ) 
  }
  samp_chars = samp_chars[-1, ]
  if (!is.null(dim(samp_chars))){
  newnames = unique(unlist(lapply(lapply(samp_chars, FUN = function(x) lapply(unlist(stringr::str_split(x, ":")),"[[", 1)), "[[", 1 )))
  newnames = newnames[newnames != ""]
  colnames(samp_chars) =  newnames
  } else { names(samp_chars) = unique(unlist(lapply(stringr::str_split(samp_chars, ":"), "[[", 1)))
  
  }
  
  meta = data.frame("ID" = as.character(unlist(meta[ which(meta[ ,1] == "!Sample_geo_accession" ), ]))[-1],
                    "Group" =  as.character(unlist(meta[ which(meta[ ,1] == "!Sample_title" ), ]))[-1]
  )
  meta = cbind(meta, samp_chars)
  if ( is.null(dim(samp_chars)) ) {
   colnames(meta)[ncol(meta)] = unique(unlist(lapply(stringr::str_split(samp_chars, ":"), "[[", 1)))
  }
  for ( i in 1:ncol(meta)) {
   if ( any(grepl(": ", meta[ ,i] )) ) {
     
    drop_word = unique(unlist(lapply(stringr::str_split(meta[ ,i], ":") , "[[", 1 )))
    drop_word = drop_word[drop_word != ""]
    meta[ ,i] = stringr::str_remove_all(string = meta[ ,i], pattern = drop_word)
    meta[ ,i] = stringr::str_remove_all(string = meta[ ,i], pattern = ": ")
   }
  }
  
  if ( !is.null(platform)) {
    
    fil = paste0("../Annots/Microarray/", platform, ".rds")
    if (!file.exists(fil)){
      mart = biomaRt::useMart("ensembl", dataset = species_dataset )
      map = biomaRt::getBM(mart = mart, attributes=c(platform, "external_gene_name" ), filters = platform, values = rownames(datExpr) )
      saveRDS(map, file = fil)
    } else { map = readRDS( fil )
    print(paste("PROBE ID's =", fil ))
    }
  }
  
  ## map[ ,1] should be probe ID and map[ ,2] is gene Name
  genes = map[ ,2][ match(rownames(datExpr), map[ ,1])]
  datExpr = datExpr[!is.na(genes), ]
  genes =genes[!is.na(genes)]
  
  gMat = matrix(nrow = length(unique(genes)), ncol = ncol(datExpr), data = 0)
  rownames(gMat) = unique(genes)
  colnames(gMat) = colnames(datExpr)
  
  for ( i in rownames(gMat)){
    probeID = map[ ,1][ !is.na(match(map[ ,2], i )) ]
    
   gMat[rownames(gMat) == i, ] = colMeans(datExpr[rownames(datExpr) %in% probeID, ])
     
  }

  outfile = list("Expr" = gMat, "Metadata" = meta, "Summary" = summary, "Design" = design)
  print(outname)
  
  print(lapply(outfile, dim ))
  
  saveRDS(outfile, file.path( dirname(series_matrix_path), paste0(outname, ".rds" )))
  
  
}



# Merge multiple GEO series -----------------------------------------------

merge_GEO = function( Series ){ 
  ## slow but easy to do
  meta = Series[[1]]$Metadata[ ,1:3]
  colnames(meta) = c("GEO_ID", "Series_ID", "Group")
  meta$Batch = 1
  
  Expr = Series[[1]]$Expr %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column( "Gene")
  
  for ( i in 2:length(Series)){
    newMeta =  Series[[i]]$Metadata[ ,1:3]
    colnames(newMeta) = c("GEO_ID", "Series_ID", "Group")
    newMeta$Batch = i
    meta = rbind(meta, newMeta[ ,c(1:3, ncol(newMeta))] )
  }
  
  
  for ( i in 2:length(Series)){
    datExpr =  Series[[i]]$Expr %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column( "Gene")
    
    Expr = inner_join(Expr, datExpr, by = "Gene" )
  }
  
  
  Expr = Expr %>% 
    tibble::column_to_rownames("Gene")
  
  return(list("datExpr" = Expr, "Annots" = meta ))  
}

