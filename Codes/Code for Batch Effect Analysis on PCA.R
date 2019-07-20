##################################################################################################
# Batch effect analysis of genes and Promoters
###################################################################################################
# requires PC_batch data frames and surveillance package to function
library(surveillance)

###################################################################################################
# tests for genes
  # Kruskal-Wallis for Provider, Date and Age
  g_pvalues_prov <- c()
  g_pvalues_date <- c()
  g_pvalues_age <- c()
  
  for (i in 1:5) {
    p <- (kruskal.test(g_PC_batch[,i]~g_PC_batch[,6]))$p.value
    d <- (kruskal.test(g_PC_batch[,i]~g_PC_batch[,7]))$p.value
    a <- (kruskal.test(g_PC_batch[,i]~g_PC_batch[,9]))$p.value
    g_pvalues_prov <- append(g_pvalues_prov, p)
    g_pvalues_date <- append(g_pvalues_date, d)
    g_pvalues_age <- append(g_pvalues_age, a)
  }
  
  # Monte carlo for runs
  g_pvalues_runs <- c()
  
  for (i in 1:5) {
    r <- (permutationTest(g_PC_batch[,i], g_PC_batch[,8], nPermutation = 9999))$pVal.permut
    g_pvalues_runs <- append(g_pvalues_runs, r)
  }
  
  # Willcoxon for cell type and sex
  g_pvalues_CT <- c()
  g_pvalues_sex <- c()
  
  for (i in 1:5) {
    c <- (wilcox.test(g_PC_batch[,i] ~ g_PC_batch[,10]))$p.value
    s <- (wilcox.test(g_PC_batch[,i] ~ g_PC_batch[,11]))$p.value    
    g_pvalues_CT <- append(g_pvalues_CT, c)
    g_pvalues_sex <- append(g_pvalues_sex, s)
  }
  
  # removing uneccesary data sets
  remove(p, d, a, r, c, s, i)

# tests for promoters
  # Kruskal-Wallis for Provider (_provider), Date (_date) and Age (_age)
    p_pvalues_prov <- c()
    p_pvalues_date <- c()
    p_pvalues_age <- c()
    
    for (i in 1:5) {
      p <- (kruskal.test(p_PC_batch[,i]~p_PC_batch[,6]))$p.value
      d <- (kruskal.test(p_PC_batch[,i]~p_PC_batch[,7]))$p.value
      a <- (kruskal.test(p_PC_batch[,i]~p_PC_batch[,9]))$p.value
      p_pvalues_prov <- append(p_pvalues_prov, p)
      p_pvalues_date <- append(p_pvalues_date, d)
      p_pvalues_age <- append(p_pvalues_age, a)
    }
  
  # Monte carlo for runs
    p_pvalues_runs <- c()
  
    for (i in 1:5) {
      r <- (permutationTest(p_PC_batch[,i], p_PC_batch[,8], nPermutation = 9999))$pVal.permut
      p_pvalues_runs <- append(p_pvalues_runs, r)
    }
  
  # Willcoxon for cell type and sex
    p_pvalues_CT <- c()
    p_pvalues_sex <- c()
    
    for (i in 1:5) {
        c <- (wilcox.test(p_PC_batch[,i] ~ p_PC_batch[,10]))$p.value
        s <- (wilcox.test(p_PC_batch[,i] ~ p_PC_batch[,11]))$p.value    
        p_pvalues_CT <- append(p_pvalues_CT, c)
        p_pvalues_sex <- append(p_pvalues_sex, s)
    }

  # removing uneccesary data sets
  remove(p, d, a, r, c, s, i)
######################################################################################################
# Summarizing values in a data frame
  # combining and formatting data set
  g_PC_pvalues <- data.frame(g_pvalues_prov, g_pvalues_date, g_pvalues_age, g_pvalues_runs, g_pvalues_CT, g_pvalues_sex)
  rownames(g_PC_pvalues) = c("PC1", "PC2", "PC3", "PC4", "PC5")
  colnames(g_PC_pvalues) <- c("Provider", "Date", "Age", "Runs", "Cell type", "Sex")
  
  p_PC_pvalues <- data.frame(p_pvalues_prov, p_pvalues_date, p_pvalues_age, p_pvalues_runs, p_pvalues_CT, p_pvalues_sex)
  rownames(p_PC_pvalues) = c("PC1", "PC2", "PC3", "PC4", "PC5")
  colnames(p_PC_pvalues) <- c("Provider", "Date", "Age", "Runs", "Cell type", "Sex")
  
  g_PC_pvalues <- t(g_PC_pvalues)
  p_PC_pvalues <- t(p_PC_pvalues)

  # removing uneccesary data sets
  remove(g_pvalues_prov, g_pvalues_date, g_pvalues_age, g_pvalues_runs, g_pvalues_CT, g_pvalues_sex)
  remove(p_pvalues_prov, p_pvalues_date, p_pvalues_age, p_pvalues_runs, p_pvalues_CT, p_pvalues_sex)

#####################################################################################################
# Creating heatmaps
  # Binary Heatmap
    # creating data set
    g_PC_pvalues_bin <- ifelse(g_PC_pvalues <= 0.05, 0, 1)
    p_PC_pvalues_bin <- ifelse(p_PC_pvalues <= 0.05, 0, 1)
  
    # plotting heatmap
    batch_colors <- colorRampPalette(c("#FF9500", "#FFFF66"))
    
    heatmap(g_PC_pvalues_bin
            , scale = "none"
            , col = batch_colors(2)
            , main = "Heatmap of Genes Batch Effects (Significance <= 0.05)"
            , xlab = "Kruskal (P, D, A), Permutation (R), Wilcoxon (CT, S)"
            , Rowv = NA
            , Colv = NA
    )
    
    heatmap(p_PC_pvalues_bin
            , scale = "none"
            , col = batch_colors(2)
            , main = "Heatmap of Promoters Batch Effects (Significance <= 0.05)"
            , xlab = "Kruskal (P, D, A), Permutation (R), Wilcoxon (CT, S)"
            , Rowv = NA
            , Colv = NA
    )

    # Heatmap with continous colors
    heatmap(as.matrix(g_PC_pvalues)
            , scale = "none"
            , col = batch_colors(256)
            , main = "Heatmap of Batch Effect p-values (Genes)"
            , xlab = "Kruskal (P, D, A), Permutation (R), Wilcoxon (CT, S)"
            , Rowv = NA
            , Colv = NA
    )
    
    heatmap(as.matrix(p_PC_pvalues)
            , scale = "none"
            , col = batch_colors(256)
            , main = "Heatmap of Batch Effect p-values (Promoters)"
            , xlab = "Kruskal (P, D, A), Permutation (R), Wilcoxon (CT, S)"
            , Rowv = NA
            , Colv = NA
    )

    # removing uneccesary data sets
    remove(g_PC_pvalues_bin ,batch_colors)
    remove(p_PC_pvalues_bin)