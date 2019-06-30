# Abstarcting PC1 for genes.
  g_pca_abs <- abs(g_pca$rotation[,1:5])
  g_PC1_ranked <- sort(g_pca_abs[,1], decreasing = TRUE)
  g_PC1_ranked_df <- data.frame(g_PC1_ranked)
  remove(g_PC1_ranked)
# Removing the lines with same values for each patient from the g_Mvalues, because they can't be anlysed with PCA.
  g_T_Mvalues[ , apply(g_T_Mvalues, 2, var) == 0]
 
    # Output: ENSG00000201801 ENSG00000269971 ENSG00000202354 ENSG00000238923 ENSG00000273667

  for_kmeans1 <- g_Mvalues[!rownames(g_Mvalues) %in% "ENSG00000201801" , ]
  for_kmeans2 <- for_kmeans1[!rownames(for_kmeans1) %in% "ENSG00000269971" , ]
  for_kmeans3 <- for_kmeans2[!rownames(for_kmeans2) %in% "ENSG00000202354" , ]
  for_kmeans4 <- for_kmeans3[!rownames(for_kmeans3) %in% "ENSG00000238923" , ]
  for_kmeans <- for_kmeans4[!rownames(for_kmeans4) %in% "ENSG00000273667" , ]
  remove(for_kmeans1, for_kmeans2, for_kmeans3, for_kmeans4)
# Binding PC1 ranked dataframe with the M values for clustering. 
  g_PC1_M <- merge(g_PC1_ranked_df, for_kmeans, by="row.names")
  g_PC1_M[order(g_PC1_M$g_PC1_ranked, decreasing = TRUE),]
  rownames(g_PC1_M)=g_PC1_M$Row.names
# Crealing dataframes from 5-100% in 5% steps of g_PC1_M. 