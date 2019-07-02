# Abstarcting PC1 for genes.
  g_pca_abs <- abs(g_pca$rotation[,1:5])
  g_PC1_ranked <- sort(g_pca_abs[,1], decreasing = TRUE)
  g_PC1_ranked_df <- data.frame(g_PC1_ranked)
  remove(g_PC1_ranked)
  
# Abstarcting PC1 for promoters.
  p_pca_abs <- abs(p_pca$rotation[,1:5])
  p_PC1_ranked <- sort(p_pca_abs[,1], decreasing = TRUE)
  p_PC1_ranked_df <- data.frame(p_PC1_ranked)
  remove(p_PC1_ranked)
# Removing the genes with 0 standard deviation from the g_Mvalues, because they can't be anlysed with PCA.
  # Code: g_T_Mvalues[ , apply(g_T_Mvalues, 2, var) == 0]
   # Output: ENSG00000201801 ENSG00000269971 ENSG00000202354 ENSG00000238923 ENSG00000273667

  for_kmeans1 <- g_Mvalues[!rownames(g_Mvalues) %in% "ENSG00000201801" , ]
  for_kmeans2 <- for_kmeans1[!rownames(for_kmeans1) %in% "ENSG00000269971" , ]
  for_kmeans3 <- for_kmeans2[!rownames(for_kmeans2) %in% "ENSG00000202354" , ]
  for_kmeans4 <- for_kmeans3[!rownames(for_kmeans3) %in% "ENSG00000238923" , ]
  for_kmeans <- for_kmeans4[!rownames(for_kmeans4) %in% "ENSG00000273667" , ]
  remove(for_kmeans1, for_kmeans2, for_kmeans3, for_kmeans4)
  
# Binding PC1 ranked dataframe with the M values for clustering for genes. 
  g_PC1_M <- merge(g_PC1_ranked_df, for_kmeans, by="row.names")
  remove(for_kmeans, g_PC1_ranked_df)
  g_PC1_M <- g_PC1_M[order(g_PC1_M$g_PC1_ranked, decreasing = TRUE),]
  rownames(g_PC1_M)=g_PC1_M$Row.names
  g_PC1_M <- g_PC1_M[,-1]
  g_PC1_M <- g_PC1_M[,-1]
  g_PC1_MTrans <- t(g_PC1_M)
  rownames(g_PC1_MTrans) <- c(rep("AML", times = 10), rep("Mono", times = 10))
  g_PC1_MTrans <- data.frame(g_PC1_MTrans)
  
# Binding PC1 ranked dataframe with the M values for clustering for promoters. 
  p_PC1_M <- merge(p_PC1_ranked_df, p_Mvalues, by="row.names")
  p_PC1_M <- p_PC1_M[order(p_PC1_M$p_PC1_ranked, decreasing = TRUE),]
  rownames(p_PC1_M)=p_PC1_M$Row.names
  p_PC1_M <- p_PC1_M[,-1]
  p_PC1_M <- p_PC1_M[,-1]
  p_PC1_MTrans <- t(p_PC1_M) 
  rownames(p_PC1_MTrans) <- c(rep("AML", times = 10), rep("Mono", times = 10))
  p_PC1_MTrans <- data.frame(p_PC1_MTrans)
  
# Crealing dataframes 5-100% in 5% steps of g_PC1_M for clustering.  
  g_PC1_5 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.05)])
  g_PC1_10 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.1)])
  g_PC1_15 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.15)])
  g_PC1_20 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.2)])
  g_PC1_25 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.25)])
  g_PC1_30 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.3)])
  g_PC1_35 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.35)])
  g_PC1_40 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.4)])
  g_PC1_45 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.45)])
  g_PC1_50 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.5)])
  g_PC1_55 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.55)])
  g_PC1_60 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.6)])
  g_PC1_65 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.65)])
  g_PC1_70 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.7)])
  g_PC1_75 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.75)])
  g_PC1_80 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.8)])
  g_PC1_85 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.85)])
  g_PC1_90 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.9)])
  g_PC1_95 <- as.data.frame(g_PC1_MTrans[,1:(48961*0.95)])
  g_PC1_100 <- as.data.frame(g_PC1_MTrans[,1:(48961*1)])
  
# Creating vectrors with only the clusters for each data frame. 
  km5 <- kmeans(g_PC1_5, centers = 2)
  km5 <- km5$cluster
  
  km10 <- kmeans(g_PC1_10, centers = 2)
  km10 <- km10$cluster
  
  km15 <- kmeans(g_PC1_15, centers = 2)
  km15 <- km15$cluster
  
  km20 <- kmeans(g_PC1_20, centers = 2)
  km20 <- km20$cluster
  
  km25 <- kmeans(g_PC1_25, centers = 2)
  km25 <- km25$cluster
  
  km30 <- kmeans(g_PC1_30, centers = 2)
  km30 <- km30$cluster
  
  km35 <- kmeans(g_PC1_35, centers = 2)
  km35 <- km35$cluster
  
  km40 <- kmeans(g_PC1_40, centers = 2)
  km40 <- km40$cluster
  
  km45 <- kmeans(g_PC1_45, centers = 2)
  km45 <- km45$cluster
  
  km50 <- kmeans(g_PC1_50, centers = 2)
  km50 <- km50$cluster
  
  km55 <- kmeans(g_PC1_55, centers = 2)
  km55 <- km55$cluster
  
  km60 <- kmeans(g_PC1_60, centers = 2)
  km60 <- km60$cluster
  
  km65 <- kmeans(g_PC1_65, centers = 2)
  km65 <- km65$cluster
  
  km70 <- kmeans(g_PC1_70, centers = 2)
  km70 <- km70$cluster
  
  km75 <- kmeans(g_PC1_75, centers = 2)
  km75 <- km75$cluster
  
  km80 <- kmeans(g_PC1_80, centers = 2)
  km80 <- km80$cluster
  
  km85 <- kmeans(g_PC1_85, centers = 2)
  km85 <- km85$cluster
  
  km90 <- kmeans(g_PC1_90, centers = 2)
  km90 <- km90$cluster
  
  km95 <- kmeans(g_PC1_95, centers = 2)
  km95 <- km95$cluster
  
  km100 <- kmeans(g_PC1_100, centers = 2)
  km100 <- km100$cluster
  
  g_kmeans_all <- data.frame(cbind(km5, km10, km15, km20, km25, km30, km35, km40, km45, km50, km55, km60, km65, km70, km75, km80, km85, km90, km95, km100))
  
  remove(g_PC1_5, g_PC1_10, g_PC1_15, g_PC1_20, g_PC1_25, g_PC1_30, g_PC1_35, g_PC1_40, g_PC1_45, g_PC1_50
         , g_PC1_55, g_PC1_60, g_PC1_65, g_PC1_70, g_PC1_75, g_PC1_80, g_PC1_85, g_PC1_90, g_PC1_95, g_PC1_100)
  remove(km5, km10, km15, km20, km25, km30, km35, km40, km45, km50, km55, km60, km65, km70, km75, km80, km85, km90, km95, km100)
  