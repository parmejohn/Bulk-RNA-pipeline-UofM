set.seed(333)

##### MAIN #####
PerformDGETests <- function(cnts, 
                            normalized.cnts, 
                            species, 
                            sample.table, 
                            filter.chr,
                            filter.gsea.genes){
  
  manual.gene.mapping <- cnts[,1:2]
  print(filter.chr)  
  print(filter.gsea.genes)

  if (length(filter.chr) == 1 & filter.chr[1] == "none"){
	print("no chr filter")	
  } else {
    if (species == 'Homo sapiens') {
      mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      gene_mapping <- getBM(attributes = c("ensembl_gene_id", "chromosome_name"),
                            filters = "ensembl_gene_id", 
                            values = cnts$gene_id, 
                            mart = mart)
    } else if (species == 'Mus musculus'){
      mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
      gene_mapping <- getBM(attributes = c("ensembl_gene_id", "chromosome_name"),
                            filters = "ensembl_gene_id", 
                            values = cnts$gene_id, 
                            mart = mart)
    }
    names(gene_mapping)[names(gene_mapping) == 'ensembl_gene_id'] <- 'gene_id'
    cnts <- merge(cnts, gene_mapping, by = "gene_id")
    cnts <- cnts[cnts$chromosome_name != filter.chr, ] %>% dplyr::select(-c("chromosome_name"))
    
    normalized.cnts <- merge(normalized.cnts, gene_mapping, by = "gene_id")
    normalized.cnts <- normalized.cnts[normalized.cnts$chromosome_name != filter.chr, ] %>% dplyr::select(-c("chromosome_name"))
  }
  
  rownames(cnts) <- NULL
  rownames(normalized.cnts) <- NULL
  cnts <- column_to_rownames(cnts, "gene_id") %>% dplyr::select(-1)
  normalized.cnts <- column_to_rownames(normalized.cnts, "gene_id") %>% dplyr::select(-1)
  
  ## Perform DESeq2 analyses
  dds <- DESeqDataSetFromMatrix(countData = round(cnts),
                                colData = sample.table,
                                design = ~0 + condition
                                )
  
  smallestGroupSize <- 3
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds <- dds[keep,]
  dds <- DESeq(dds)
  
  ### Perform every pairwise comparison possible
  m_df <- msigdbr(species = species, category = "C5", subcategory = "BP") # dont need to reload the dataset every time for GSEA
  fgsea.sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  group.pairs <- as.data.frame(combn(unique(sample.table$condition), 2))
  group.pairs <- sapply(group.pairs, function(x) as.character(x), simplify = FALSE)
  res.df <- NA
  res.df.list <- list()
  j = 0
  for (i in 1:length(group.pairs)){
    res <- results(dds, contrast = c("condition", group.pairs[[i]]))
    res$ensemblID <- rownames(res)
    if (species == 'Homo sapiens') {
      mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      gene_mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                            filters = "ensembl_gene_id", 
                            values = res$ensemblID, 
                            mart = mart)
    } else if (species == 'Mus musculus'){
      mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
      gene_mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "mgi_id"),
                            filters = "ensembl_gene_id", 
                            values = res$ensemblID, 
                            mart = mart)
    }
    names(gene_mapping)[names(gene_mapping) == 'ensembl_gene_id'] <- 'ensemblID'
    names(manual.gene.mapping)[names(manual.gene.mapping) == 'gene_id'] <- 'ensemblID'
    
    res.df <- merge(as.data.frame(res), manual.gene.mapping, by = "ensemblID")
    res.df <- merge(as.data.frame(res.df), gene_mapping, by = "ensemblID")
    
    res.df <- subset(res.df, !is.na(padj) |
                       !is.na(pvalue) |
                       !is.na(log2FoldChange))
    
    write.table(res.df, paste0("deseq2_", group.pairs[[i]][1], "_vs_", group.pairs[[i]][2], '_res.txt'), row.names = F, quote = F)

    p1 <- deg_volcano_plot(res.df, group.pairs[[i]][1], group.pairs[[i]][2])
    p2 <- GseaComparison(res.df, fgsea.sets, group.pairs[[i]][1], group.pairs[[i]][2], filter.gsea = filter.gsea.genes)
    p3 <- deg_heatmap(res.df, normalized.cnts, group.pairs[[i]][1], group.pairs[[i]][2], sample.table)
    
    rownames(res.df) <- NULL
    res.df.list <- append(res.df.list, list(res.df))
    j = j + 1
    names(res.df.list)[i] <- paste0(group.pairs[[i]][1], "_vs_", group.pairs[[i]][2])
  }
  
  ## creating upset plot from all comparisons and intersections
  res.df.list.up <- lapply(res.df.list, subset, padj < 0.05 & log2FoldChange >= 2)
  res.df.list.dn <- lapply(res.df.list, subset, padj < 0.05 & log2FoldChange <= -2)
  upset.plot.up <- deg_upset(res.df.list.up)
  upset.plot.dn <- deg_upset(res.df.list.dn)
  
  pdf("deg_upset_upreg.pdf", width = 10, height = 6, onefile=FALSE)
  print(upset.plot.up)
  grid.text("Upregulated genes\n padj < 0.05 & log2FC >= 2",x = 0.20, y=0.85, gp=gpar(fontsize=16))
  graphics.off()
  
  pdf("deg_upset_dnreg.pdf", width = 10, height = 6, onefile=FALSE)
  print(upset.plot.dn)
  grid.text("Downregulated genes\n padj < 0.05 & log2FC <= -2",x = 0.20, y=0.85, gp=gpar(fontsize=16))
  graphics.off()
  
  
  ## create PCA and heatmap to determine sample similarity
  # variance stabilizing transformation of the count matrix
  vsd <- vst(dds, blind = FALSE)

  sampleDists <- dist(t(assay(vsd)))
  
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- paste(sample.table$sample, sample.table$condition, sep ="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  # create hiearchical clustering heatmap from the vsd matrix
  pdf("sample_similarity_heatmap.pdf", width=8, height=8)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
  graphics.off()

  pcaPlot <- plotPCA(vsd)
  ggsave("sample_similarity_pca.pdf", pcaPlot, width = 8, height = 8)
}

##### FXNS #####
## helper fxns
deg_volcano_plot <- function(de_markers, ident.1, ident.2, p_val_adj_cutoff = 0.05, avg_log2FC_cutoff = 2){
  
  keyvals <- ifelse(
    de_markers$log2FoldChange <= -avg_log2FC_cutoff & de_markers$padj < p_val_adj_cutoff, 'blue',
    ifelse(de_markers$log2FoldChange >= avg_log2FC_cutoff & de_markers$padj < p_val_adj_cutoff, 'red',
           'grey'))
  keyvals[is.na(keyvals)] <- 'grey'
  names(keyvals)[keyvals == 'red'] <- paste0('Upregulated in ', ident.1)
  names(keyvals)[keyvals == 'grey'] <- 'Not significant'
  names(keyvals)[keyvals == 'blue'] <- paste0('Upregulated in ', ident.2)
  
  p <- EnhancedVolcano(de_markers,
                  lab = de_markers$gene_name,
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  FCcutoff=avg_log2FC_cutoff,
                  pCutoff=p_val_adj_cutoff,
                  pCutoffCol="padj",
                  pointSize = 1.0,
                  labSize = 3.0,
                  drawConnectors = TRUE,
                  arrowheads = F, 
                  colCustom = keyvals,
                  legendPosition = 'right',
                  title = "DESeq2: Volcano Plot",
                  subtitle = paste0(ident.1, " vs ", ident.2),
                  caption = paste0("total = ", nrow(de_markers), " genes"),
                  legendIconSize = 2.5,
                  legendLabSize = 11)
  
  ggsave(filename = paste0("deseq2_volcano_", ident.1, "_vs_", ident.2, '.pdf'), plot = p, width=8, height=8)
}

GseaComparison <- function(de.markers, fgsea.sets, ident.1, ident.2, filter.gsea){
  cluster.genes <- de.markers %>%
    arrange(desc(log2FoldChange)) %>% 
    dplyr::select(hgnc_symbol, log2FoldChange) # use avg_log2FC as ranking for now; https://www.biostars.org/p/9526168/

	print(filter.gsea)

  if (length(filter.gsea) == 1 & filter.gsea[1] == "none"){
	print("no filter applied")
	} else { 
	  cluster.genes <- cluster.genes[cluster.genes$hgnc_symbol != filter.gsea, ]
  }
  
  ranks <- deframe(cluster.genes)
  
  fgseaRes <- fgsea(fgsea.sets, stats = ranks)
  fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ", "))
  write.table(fgseaRes, paste0("gsea_", ident.1, "_vs_", ident.2,".txt"), quote = FALSE, row.names = F, sep = "\t", col.names = T)
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
  
  sub.samples <- subset(sample.table, condition == ident.1 | condition == ident.2)

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
  row_labels = structure(de.markers.filtered$gene_name, names = de.markers.filtered$ensemblID)
  row_labels[is.na(row_labels) | row_labels == ""] <- names(row_labels)[is.na(row_labels) | row_labels == ""]
  
  counts.df.filtered <- counts.df.filtered %>% dplyr::select(sub.samples$sample)

  
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
    row_names_gp = gpar(fontsize = 3),
    column_title = paste0(nrow(counts.df.filtered), " DEGs for ", ident.1, " vs ", ident.2), 
    column_title_gp = gpar(fontsize = 15, fontface = "bold")
  )
  p <- draw(hmap, heatmap_legend_side="left", annotation_legend_side="right")
  
  pdf(paste0("deg_heatmap_", ident.1, "_vs_", ident.2, ".pdf"), width = 10, height = 10)
  print(p)
  graphics.off()
  
  return(p)
}

deg_upset <- function(res.df.list) {
  for (i in 1:length(res.df.list)) {
    res.df.list[[i]] <- (res.df.list[[i]]$ensemblID)
  }
  
  upset.df <- fromList(res.df.list)
  
  p <- upset(upset.df,
        nsets = length(res.df.list),
        order.by = "freq",
        text.scale = 1)
}
