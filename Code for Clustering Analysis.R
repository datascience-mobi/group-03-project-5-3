# setting up loadings data frames
g_PC1_df <- data.frame(abs(g_pca$rotation[,1]))
colnames(g_PC1_df) <- c("g_PC1_df")

# Binding PC1 ranked dataframe with the M values for clustering for genes. 
g_PC1_M <- merge(g_PC1_df, g_Mvalues, by="row.names")
remove(g_PC1_df)


g_PC1_M <- g_PC1_M[order(g_PC1_M$g_PC1_df, decreasing = TRUE),]
rownames(g_PC1_M)=g_PC1_M$Row.names
g_PC1_M <- g_PC1_M[, 3:22]

fg_cluster <- function(x) {
  
  if(x %% jcount == 0) {
    cat("|")
  }
  
  x <- x/100
  
  g_PC1_M_reduced <- g_PC1_M[1:floor(nrow(g_PC1_M)*x), ]
  g_cluster <- t(g_PC1_M_reduced)
  rownames(g_cluster) <- c(rep("AML", times = 10), rep("Mono", times = 10))
  
  set.seed(x)
  g_c <- kmeans(g_cluster, centers = 2, nstart = 100)
  g_c <- g_c$cluster
  
  return(g_c)
}

j <- data.frame(seq(1, 100, 1))
jcount <- floor(nrow(j)/100)
g_c <- apply(j, 1, fg_cluster)
g_c <- t(g_c)


fgc_accuracy <- function(x) {
  
  workingrow <- g_c[x, ]
  cluster1_assA <- length(which(workingrow[1:10] == 1))
  cluster1_assM <- length(which(workingrow[11:20] == 1))
  
  if(cluster1_assA >= cluster1_assM) {
    cluster1_ass <- c("A")
    workingrow[workingrow == 1] <- cluster1_ass
    cluster2_ass <- c("M")
    workingrow[workingrow == 2] <- cluster2_ass
  } else {
    cluster1_ass <- c("M")
    workingrow[workingrow == 1] <- cluster1_ass
    cluster2_ass <- c("A")
    workingrow[workingrow == 2] <- cluster2_ass
  }
  
  workingrow_AML <- workingrow[1:10]
  workingrow_AML[workingrow_AML == c("A")] <- TRUE 
  workingrow_AML[workingrow_AML == c("M")] <- FALSE
  
  workingrow_mon <- workingrow[11:20]
  workingrow_mon[workingrow_mon == c("M")] <- TRUE
  workingrow_mon[workingrow_mon == c("A")] <- FALSE 
  
  workingrow <- as.logical(c(workingrow_AML, workingrow_mon))
  
  return(mean(workingrow))
}

j <- data.frame(seq(1, nrow(g_c), 1))
gc_acc <- apply(j, 1, fgc_accuracy)
View(gc_acc)

remove(fg_cluster, fgc_accuracy, g_c, g_PC1_M, j, jcount)