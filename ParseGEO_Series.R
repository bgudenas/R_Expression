

#series_matrix_path = "/Users/bgudenas/Desktop/Projects/Datasets/Human/Expression/Microarray/MB/GSE37382/GSE37382_series_matrix.txt"

ParseGEO_Microarray = function(series_matrix_path ) {
  options(stringsAsFactors = FALSE)
  all = read.table(series_matrix_path, fill = TRUE, sep = "\t",header = FALSE )
  summary = as.character(all[ which(all[ ,1] == "!Series_summary"), 2 ])
  design =  as.character(all[ which(all[ ,1] == "!Series_overall_design"), 2 ])
  meta = read.table(series_matrix_path, fill = TRUE, sep = "\t", skip = 25, header = FALSE )
  markExpr = which(meta[ ,1] == "ID_REF") + 25
  datExpr = read.table(series_matrix_path, fill = TRUE, sep = "\t", skip = markExpr, header = TRUE, row.names = 1) 
  
  
  
  which(meta[ ,1] == "!Sample_characteristics_ch1" )
  samp_chars = c()
  for ( i in  which(meta[ ,1] == "!Sample_characteristics_ch1" ) ) {
    
   samp_chars = cbind(samp_chars, as.character(unlist(meta[ i, ])) ) 
  }
  samp_chars = samp_chars[-1, ]
  if (!is.null(dim(samp_chars))){
  colnames(samp_chars) =  unique(unlist(lapply(lapply(samp_chars, FUN = function(x) lapply(unlist(stringr::str_split(x, ":")),"[[", 1)), "[[", 1 )))
  } else { names(samp_chars) = unique(unlist(lapply(stringr::str_split(samp_chars, ":"), "[[", 1)))
  
  }
  
  meta = data.frame("ID" = as.character(unlist(meta[ which(meta[ ,1] == "!Sample_geo_accession" ), ]))[-1],
                    "Group" =  as.character(unlist(meta[ which(meta[ ,1] == "!Sample_title" ), ]))[-1]
  )
  meta = cbind(meta, samp_chars)
  if ( is.null(dim(samp_chars)) ) {
   colnames(meta)[ncol(meta)] = unique(unlist(lapply(stringr::str_split(samp_chars, ":"), "[[", 1)))
  }
  
  
  outfile = list("Expr" = datExpr, "Metadata" = meta, "Summary" = summary, "Design" = design)
  
  
}


mega = ParseGEO_Microarray(series_matrix_path = series_matrix_path)

mega = ParseGEO_Microarray("./Public_RBmodel/GSE124537-GPL6246_series_matrix.txt")
