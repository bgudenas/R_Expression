# 6/7/19
#Brian Gudenas
## Make gene annotation tables for different species containing gene-level info

MakeHgAnno = function(path){
  library(biomaRt)
  dataset = "hsapiens_gene_ensembl"
  mart = biomaRt::useMart("ensembl", dataset = dataset )
  map = biomaRt::getBM(mart = mart, attributes=c("ensembl_gene_id","external_gene_name","hgnc_symbol", "gene_biotype", 
                                                 "chromosome_name", "start_position", "end_position","description",
                                                 "transcript_count", "strand", "band" ))
  
  map2 = biomaRt::getBM(mart = mart, attributes=c("ensembl_gene_id", 'mmusculus_homolog_ensembl_gene', "mmusculus_homolog_associated_gene_name" ) )
  map3 = biomaRt::getBM(mart = mart, attributes=c("ensembl_gene_id","entrezgene_id"))
  
  map$mmusculus_homolog_associated_gene_name = map2$mmusculus_homolog_associated_gene_name[ match(map$ensembl_gene_id, map2$ensembl_gene_id) ]
  map$mmusculus_homolog_ensembl_gene = map2$mmusculus_homolog_ensembl_gene[ match(map$ensembl_gene_id, map2$ensembl_gene_id) ]
  map$entrezgene_id = map3$entrezgene_id[ match(map$ensembl_gene_id, map3$ensembl_gene_id) ]

  saveRDS(map,  path )
}
#MakeHgAnno( "/Users/bgudenas/Desktop/Projects/Annots/Annotables/hg38.rds" )


MakeMmAnno = function(path){
  dataset = "mmusculus_gene_ensembl"
  mart = biomaRt::useMart("ensembl", dataset = dataset )
  map = biomaRt::getBM(mart = mart, attributes=c("ensembl_gene_id","entrezgene","external_gene_name","hgnc_symbol", "gene_biotype", 
                                                 "chromosome_name", "start_position", "end_position","description",
                                                 "transcript_count", "strand", "band" ))
  
  map2 = biomaRt::getBM(mart = mart, attributes=c("ensembl_gene_id",'hsapiens_homolog_ensembl_gene', "hsapiens_homolog_associated_gene_name" ) )
  
  map$hsapiens_homolog_associated_gene_name = map2$hsapiens_homolog_associated_gene_name[ match(map$ensembl_gene_id, map2$ensembl_gene_id) ]
  map$hsapiens_homolog_ensembl_gene = map2$hsapiens_homolog_ensembl_gene[ match(map$ensembl_gene_id, map2$ensembl_gene_id) ]
  
  saveRDS(map,  path )
}
#MakeMmAnno( "/Users/bgudenas/Desktop/Projects/Annots/Annotables/Mm10.rds" )


