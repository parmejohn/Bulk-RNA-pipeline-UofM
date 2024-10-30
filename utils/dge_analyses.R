set.seed(333)

##### MAIN #####
PerformDGETests <- function(kallisto.output, species){
  if (species == 'Homo sapiens') {
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    tx2gene <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_id"))
  } else if (species == 'Mus musculus'){
    mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    tx2gene <- transcripts(EnsDb.Mmusculus.v79, columns=c("tx_id", "gene_id"))
  }
  
  tx2gene <- as.data.frame(tx2gene[, c("tx_id", "gene_id")]) %>% dplyr::select(tx_id, gene_id)
  tx2gene <- tx2gene[!duplicated(tx2gene$tx_id), ]
  
	kallisto.output <- sort(kallisto.output)

  kallisto.abundance <- list.files(kallisto.output, pattern = "abundance.h5", full.names = T)

  sampleTable <- lapply(kallisto.output, extract_sample_info)
  sampleTable <- do.call(rbind, sampleTable)
  
  names(kallisto.abundance) <- sampleTable$sample
  
  txi <- tximport(kallisto.abundance, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)
  saveRDS(txi, "counts_object.rds")

  ## Perform DESeq2 analyses
  dds <- DESeqDataSetFromTximport(txi, colData = sampleTable, design = ~0 + condition)
  
  smallestGroupSize <- 3
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds <- dds[keep,]
  dds <- DESeq(dds)
  
  ### Perform every pairwise comparison possible
  m_df <- msigdbr(species = species, category = "C5", subcategory = "BP") # dont need to reload the dataset every time for GSEA
  fgsea.sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  group.pairs <- as.data.frame(combn(unique(sampleTable$condition), 2))
  group.pairs <- sapply(group.pairs, function(x) as.character(x), simplify = FALSE)
  res.df <- NA
  for (i in 1:length(group.pairs)){
    res <- results(dds, contrast = c("condition", group.pairs[[i]]))
    res$ensemblID <- rownames(res)
    if (species == 'Homo sapiens') {
      gene_mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                            filters = "ensembl_gene_id", 
                            values = res$ensemblID , 
                            mart = mart)
    } else if (species == 'Mus musculus'){
      gene_mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "mgi_id"),
                            filters = "ensembl_gene_id", 
                            values = res$ensemblID , 
                            mart = mart)
    }
    res.df <- merge(as.data.frame(res), gene_mapping, by.x = "ensemblID", by.y = "ensembl_gene_id", all.x = TRUE)
    if (length(unique(sampleTable$condition)) == 2){
      write.table(res.df, paste0("deseq2_", group.pairs[[i]][1], "_vs_", group.pairs[[i]][2], '_genewalk.txt'), row.names = F, quote = F)
    } else {
      write.table(res.df, paste0("deseq2_", group.pairs[[i]][1], "_vs_", group.pairs[[i]][2], '_res.txt'), row.names = F, quote = F)
    }
    p1 <- deg_volcano_plot(res.df, group.pairs[[i]][1], group.pairs[[i]][2])
    p2 <- GseaComparison(res.df, fgsea.sets, group.pairs[[i]][1], group.pairs[[i]][2])
    p3 <- deg_heatmap(res.df, as.data.frame(txi$abundance), group.pairs[[i]][1], group.pairs[[i]][2], sampleTable)
  }
  
  ### Perform 1 vs all comparison for gene walk results
  if (length(unique(sampleTable$condition)) > 2){
    for(i in 1:length(unique(sampleTable$condition))){
      res <- results(dds, 
                     contrast = list(paste0("condition", unique(sampleTable$condition)[i]),paste0("condition", unique(sampleTable$condition)[-i])),
                     listValues = c(1, -1/length(unique(sampleTable$condition)[-i]))
                     )
      res$ensemblID <- rownames(res)
      res.df <- merge(as.data.frame(res), gene_mapping, by.x = "ensemblID", by.y = "ensembl_gene_id", all.x = TRUE)
      write.table(res.df, paste0("deseq2_", unique(sampleTable$condition)[i], "_vs_all_genewalk.txt"), row.names = F, quote = F)
      
      p1 <- deg_volcano_plot(res.df, 
                             unique(sampleTable$condition)[i], 
                             "all"
                             )
      p2 <- GseaComparison(res.df, 
                           fgsea.sets, 
                           unique(sampleTable$condition)[i], 
                           "all"
                           )
      p3 <- deg_heatmap(res.df, 
                        as.data.frame(txi$abundance), 
                        unique(sampleTable$condition)[i], 
                        "all", 
                        sampleTable
                        )
      
    }
  }
}

##### FXNS #####
## helper fxns
extract_sample_info <- function(file_path) {
  # Extract the file name (without the path)
  file_name <- basename(file_path)
  
  # Use regex to capture the sample and condition from the file name
  sample <- sub(".*_[^_]+_([^_]+)_[^_]+_[^_]+$", "\\1", file_name)
  condition <- sub(".*_([^_]+)_[^_]+$", "\\1", file_name)
  
  # Create a data frame with the desired structure
  df <- data.frame(
    sample = sample,
    condition = condition,
    stringsAsFactors = FALSE
  )
  
  return(df)
}

deg_volcano_plot <- function(de_markers, ident.1, ident.2){
  p <- ggplot(de_markers, aes(log2FoldChange, -log10(pvalue))) + 
    geom_point(size = 0.5, alpha = 0.5) + 
    theme_bw() +
    ylab("-log10(unadjusted p-value)") + 
    geom_text_repel(aes(label = ifelse(((padj < 0.05 & log2FoldChange >= 2)|(padj < 0.05 & log2FoldChange <= -2)), hgnc_symbol,
                                       "")), colour = "red", size = 3) + 
    ggtitle(paste0("DESeq2 Volcano Plot: ", ident.1, " vs ", ident.2))
  ggsave(filename = paste0("deseq2_volcano_", ident.1, "_vs_", ident.2, '.pdf'), plot = p, width=8, height=8)
  
}

GseaComparison <- function(de.markers, fgsea.sets, ident.1, ident.2){
  cluster.genes <- de.markers %>%
    arrange(desc(log2FoldChange)) %>% 
    dplyr::select(hgnc_symbol, log2FoldChange) # use avg_log2FC as ranking for now; https://www.biostars.org/p/9526168/
  
  ranks <- deframe(cluster.genes)
  
  fgseaRes <- fgsea(fgsea.sets, stats = ranks)
  fgseaRes$leadingEdge <- as.character(fgseaRes$leadingEdge)
  write.table(fgseaRes, paste0("gsea_", ident.1, "_vs_", ident.2,".txt"), quote = FALSE,row.names = T, sep = "\t", col.names = T)
  fgseaRes <- dplyr::filter(fgseaRes, padj <= 0.05 & size >= 3) %>% arrange(desc(NES))
  fgseaRes$Enrichment = ifelse(fgseaRes$NES > 0, "red", "blue") 
  
  filtRes <-  rbind(head(fgseaRes, n = 10),
                    tail(fgseaRes, n = 10 ))
  
  p <- ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point(aes(fill = Enrichment, size = size), shape=21) + # size is equal to the amount of genes found in the given gene set
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    theme_bw() +
    labs(x="Pathway", y="Normalized Enrichment Score") + 
    scale_fill_identity() + 
    ggtitle(paste0("GSEA: ", ident.1, " vs ", ident.2))
  
  ggsave(filename = paste0("gsea_", ident.1, "_vs_", ident.2, '.pdf'), plot = p, width = 305, height = 152, units = "mm")
}

deg_heatmap <- function(de.markers, counts.df, ident.1, ident.2, sample.table){
  
  sample.table <- sample.table[order(sample.table$condition),]
  counts.df <- counts.df[, sample.table$sample]
  
  if (ident.2 == "all") {
    sub.samples <- sample.table
  } else {
    sub.samples <- subset(sample.table, condition == ident.1 | condition == ident.2)
  }
  
  ann <- data.frame(sub.samples$condition)
  colnames(ann) <- c('Condition')
  col = list(c=structure(brewer.pal(length(unique(sub.samples$condition)), "Set2"), 
                         names = unique(sub.samples$condition)))
  colAnn <- HeatmapAnnotation(df = ann,
                              which = 'column',
                              col = col,
                              annotation_width = unit(c(1, 4), 'cm'),
                              gap = unit(1, 'mm'))
  
  de.markers.filtered <- subset(de.markers, padj < 0.05 & abs(log2FoldChange) >= 2)
  counts.df.filtered <- counts.df %>%
    dplyr::filter(row.names(counts.df) %in% de.markers.filtered$ensemblID)
  row_labels = structure(de.markers.filtered$hgnc_symbol, names = de.markers.filtered$ensemblID)
  row_labels[is.na(row_labels) | row_labels == ""] <- names(row_labels)[is.na(row_labels) | row_labels == ""]
  if (ident.2 != "all") {
    counts.df.filtered <- counts.df.filtered %>% dplyr::select(sub.samples$sample)
  }
  heat <- t(scale(t(counts.df.filtered)))
  
  
  hmap <- Heatmap(
    heat,
    name = "expression",
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_column_dend = TRUE,
    show_row_dend = TRUE,
    row_dend_reorder = TRUE,
    column_dend_reorder = TRUE,
    width = unit(100, "mm"),
    top_annotation=colAnn,
    row_labels = row_labels[rownames(heat)],
    row_names_gp = gpar(fontsize = 4)
  )
  p <- draw(hmap, heatmap_legend_side="left", annotation_legend_side="right")
  
  pdf(paste0("deg_heatmap_", ident.1, "_vs_", ident.2, ".pdf"), width = 10, height = 10)
  print(p)
  graphics.off()
  
  return(p)
}
