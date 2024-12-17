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
    
    write.table(res.df, paste0("deseq2_", group.pairs[[i]][1], "_vs_", group.pairs[[i]][2], '_res.txt'), row.names = F, quote = F)

    p1 <- deg_volcano_plot(res.df, group.pairs[[i]][1], group.pairs[[i]][2])
    p2 <- GseaComparison(res.df, fgsea.sets, group.pairs[[i]][1], group.pairs[[i]][2], filter.gsea = filter.gsea.genes)
    p3 <- deg_heatmap(res.df, normalized.cnts, group.pairs[[i]][1], group.pairs[[i]][2], sample.table)
    
    rownames(res.df) <- NULL
    #res.df <- column_to_rownames(res.df, "ensemblID")
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
}

##### FXNS #####
## helper fxns
deg_volcano_plot <- function(de_markers, ident.1, ident.2, p_val_adj_cutoff = 0.05, avg_log2FC_cutoff = 2){
  
  unadjusted.pval.cutoff <- -log10(max(de_markers$pvalue[!is.na(de_markers$padj) & 
                                                           de_markers$padj <= p_val_adj_cutoff]))
  
  de_markers$significance <- "Not Significant"
  de_markers$significance[de_markers[["log2FoldChange"]] >= avg_log2FC_cutoff & 
                            de_markers[["padj"]] < p_val_adj_cutoff] <- paste0("Upregulated in ", ident.1)
  de_markers$significance[de_markers[["log2FoldChange"]] < -(avg_log2FC_cutoff) & 
                            de_markers[["padj"]] <= p_val_adj_cutoff] <- paste0("Upregulated in ", ident.2)
  de_markers$significance <- factor(de_markers$sig, levels = c(paste0("Upregulated in ", ident.2), "Not Significant", paste0("Upregulated in ", ident.1)))
  
  colors <- setNames(
    c("red", "grey", "blue"),
    c(paste0("Upregulated in ", ident.1), "Not Significant", paste0("Upregulated in ", ident.2))
  )
  
  p <- ggplot(de_markers, aes(log2FoldChange, -log10(pvalue), color = significance, fill = significance)) + 
    geom_point(size = 0.5, alpha = 0.5) + 
    theme_bw() +
    scale_color_manual(values = colors) +
    ylab("-log10(unadjusted p-value)") + 
    geom_text_repel(aes(label = ifelse(((padj < p_val_adj_cutoff & log2FoldChange >= avg_log2FC_cutoff)|
                                          (padj < p_val_adj_cutoff & log2FoldChange <= -avg_log2FC_cutoff)), 
                                       gene_name, "")), colour = "black", size = 3) + 
    geom_hline(yintercept=unadjusted.pval.cutoff, linetype="dashed", color = "black") +
    geom_vline(xintercept=c(-avg_log2FC_cutoff, avg_log2FC_cutoff), linetype="dashed", color = "black") +
    ggtitle(paste0("DESeq2 Volcano Plot: ", ident.1, " vs ", ident.2))
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
  row_labels = structure(de.markers.filtered$gene_name, names = de.markers.filtered$ensemblID)
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
