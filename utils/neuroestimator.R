library(neuroestimator)


NeuroestimatorResults <- function(table, species){ 
  
  Sys.setenv(RETICULATE_MINICONDA_PATH="/opt/conda")
  reticulate::use_condaenv("neuroestimator")
  
  res <- read.table(table, header=T)
  
  rownames(res) <- res[,1]
  res <- res[3:ncol(res)]
  
  res <- round(res)
  
  res <- neuroestimator(res, species=species, id_type = "ensembl_gene_id") # rmbr to change species if needed
  write.table(res, "neuroestimator_results.txt", quote = FALSE, row.names = T, sep = "\t", col.names = T)
  
}
