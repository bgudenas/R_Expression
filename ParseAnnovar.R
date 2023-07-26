### Function to parse annovar output

ParseAnnovar = function( samp, exac_filt = 0.0001, af_filt = 0.1, ad_filt=3, dp_filt = 10, GERMLINE = FALSE, RNA = FALSE) {
  library(dplyr)
  library(stringr)
  options(stringsAsFactors = FALSE)
  if (!file.exists(samp)){
    stop(paste(samp, " does not exist"))
  }
  df =  read.table( samp, sep = "\t", fill = TRUE, check.names = FALSE )
  sampID = stringr::str_split( basename(samp), "\\." )[[1]][1]
  varcol = colnames(df) == sampID
  if (sum(varcol) == 0 ){
   # sampID = str_replace_all(sampID, "-", ".") 
    varcol = colnames(df) == str_replace_all(sampID, "-", ".")
  }
  if (sum(varcol) == 0 ){
    varcol = as.logical(adist(colnames(df),sampID) <= 3)
  }
  AF = unlist(lapply(stringr::str_split(df[ ,varcol], ":"), "[[", 3))
  
  AD = unlist(lapply(stringr::str_split(df[ ,varcol], ":"), "[[", 2))
  AD = as.numeric(unlist(lapply(str_split(AD, ","), "[[", 2)))
  ##TODO extract AF for indels/multiallelic
  AF[grepl(",", AF) ] = 99
  AF= as.numeric(AF)
  DP = as.numeric( unlist(lapply(stringr::str_split(df[ ,varcol], ":"), "[[", 4)))
  if ( RNA == TRUE ) {
    DP = as.numeric( unlist(lapply(stringr::str_split(df[ ,varcol], ":"), "[[", 3)))
    AF= rep(99, length(AF))
    }
  
  if ( GERMLINE ==  TRUE){
    AF = unlist(lapply(stringr::str_split(df$INFO, ";"), "[[", 2))
    AF = unlist(lapply(stringr::str_split(AF, "="), "[[", 2))
    AF[grepl(",", AF) ] = 1.1
    AF= as.numeric(AF)
    
    AD = unlist(lapply(stringr::str_split(df[, ncol(df) ], ":"), "[[", 2))
    AD = as.numeric(unlist(lapply(stringr::str_split(AD, ","), "[[", 2)))
    
    DP =as.numeric( unlist(lapply(stringr::str_split(df[, ncol(df) ], ":"), "[[", 3)) )
  }
  
  df$AF=AF
  df$Alt_AD = AD
  df$DP = DP

  df$ExAC_ALL[ df$ExAC_ALL == "."] = 0
  df$ExAC_ALL = as.numeric(as.character( df$ExAC_ALL))
  exdrop = sum(df$ExAC_ALL > exac_filt )
 # print(paste("ExAC Filters", exdrop, "Variants", "(", round(exdrop/nrow(df) *100, 2), "% of variants )"))
    afdrop = sum(df$AF < af_filt )
  #  print(paste("Allele Frequency Filters", afdrop, "Variants", "(", round(afdrop/nrow(df) *100, 2), "% of variants )"))
    addrop = sum(df$Alt_AD < ad_filt )
   # print(paste("AD Filters", dpdrop, "Variants", "(", round(dpdrop/nrow(df) *100, 2), "% of variants )"))
    
    df = df %>% 
    dplyr::filter(ExAC_ALL < exac_filt & AF >= af_filt & Alt_AD >= ad_filt & DP >= dp_filt) %>% 
      dplyr::filter(ExonicFunc.ensGene != "synonymous SNV" )
    
    df$ExonicFunc.ensGene[df$ExonicFunc.ensGene == "." & df$Func.ensGene == "splicing"] = "splicing"
    
    df = df[ ,c(1:56,  (ncol(df)-2):ncol(df) )]
  
  return(df)
}
