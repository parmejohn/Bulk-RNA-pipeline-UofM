NeuroestimatorPlot <- function(countsTable, sampleTable){
  res <- read.table(countsTable)
  sampleTableData <- read.csv(sampleTable)
  res <- rownames_to_column(res, "sample")
  
  sampleTableData <- sampleTableData[order(sampleTableData$sample, decreasing = F),]
  res <- res[order(res$sample, decreasing = F),]
  
  res$group <- sampleTableData$condition
  
  # Perform pairwise Kolmogorov-Smirnov tests
  # The Kolmogorov Smirnov test (KS test or K-S test) is used to compare two distributions to determine if they are pulling from the same underlying distribution
  # Find if there is significant differense in distributions
  stat.test <- combn(unique(res$group), 2, simplify = FALSE) %>%
    map_df(~{
      group1 <- .x[1]
      group2 <- .x[2]
      ks_result <- ks.test(
        res$predicted_activity[res$group == group1],
        res$predicted_activity[res$group == group2]
      )
      tibble(
        group1 = group1,
        group2 = group2,
        p = ks_result$p.value
      )
    }) %>%
    mutate(
      p.adj = p.adjust(p, method = "BH"),
      p.adj.signif = case_when(
        p.adj < 0.001 ~ "***",
        p.adj < 0.01 ~ "**",
        p.adj < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      y.position = seq(max(res$predicted_activity) * 1.1, 
                       by = 0.1 * max(res$predicted_activity), 
                       length.out = n())
    )
  
  # Boxplot with significance annotations
   p <- ggplot(res, aes(x = group, y = predicted_activity)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.7) +
    stat_pvalue_manual(
      stat.test,
      label = "p.adj.signif",
      tip.length = 0.01,
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position"
    ) +
    labs(
      title = "Predicted Activity by Condition",
      x = "Condition",
      y = "Predicted Activity"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  ggsave("neuroestimator_results.pdf", p, width = 8, height = 8)
}