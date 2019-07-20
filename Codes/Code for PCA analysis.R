###################################################################################################
# Code for PCA analysis
  # Transposing of g_Mvalues, renaming of rownames for grouping in ggbiplot, removing of columns 
  # containing constant values, running PCA and coercing g_pca to a data frame
  g_T_Mvalues <- t(g_Mvalues)
  rownames(g_T_Mvalues) <- c(rep("AML", times = 10), rep("Mono", times = 10))
  g_T_Mvalues_clean <- g_T_Mvalues[ , apply(g_T_Mvalues, 2, var) != 0]
  g_pca <- prcomp(g_T_Mvalues_clean, center = TRUE, scale. = TRUE)
  g_pca_df <- as.data.frame(g_pca$x)
  
  # Transposing of p_Mvalues, renaming of rownames for grouping in ggbiplot, running PCA for genes 
  # and coercing p_pca to a data frame.
  p_T_Mvalues <- t(p_Mvalues)
  rownames(p_T_Mvalues) <- c(rep("AML", times = 10), rep("Mono", times = 10))
  p_pca <- prcomp(p_T_Mvalues, center = TRUE, scale. = TRUE)
  p_pca_df <- as.data.frame(p_pca$x)
  
  # Creating data frames for batch effect detection in genes and promoters. They contain the first 
  # five principal components and three columns  with technical parameters (biomaterial provider, 
  # first submission date and sequence runs) and three biological parametrs(Age,Sex and Cell type),
  # which are checked for batch effects.
  g_PC_batch <- cbind(g_pca_df[,1:5],
                      c(rep("Groningen", times=10), rep("Cambridge", times=4), rep("Nijmegen", times=6)),
                      c(1032,1005,1032,rep(1005, times=7),1184,1184,1032,1032,788,0,788,80,0,0),
                      c(9,10,9,rep(11, times=5),15,11,38,38,38,20,14,4,14,15,17,18),
                      c(32.5,57.5,62.5,67.5,72.5,67.5,67.5,42.5,47.5,67.5,67.5,rep(62.5, times=3),rep(47.5, times=4),42.5,67.5),
                      c(rep("myeloid", times=10), rep("monocyte", times=10)), 
                      c(rep("Female", times=5),"Male",rep("Female", times=4),"Male","Female","Male","Female",rep("Male", times=5),"Female"))
  names(g_PC_batch)[6] <- "Provider"
  names(g_PC_batch)[7] <- "Date"
  names(g_PC_batch)[8] <- "Runs"
  names(g_PC_batch)[9] <- "Age"
  names(g_PC_batch)[10] <- "Cell_type"
  names(g_PC_batch)[11] <- "Sex"
  
  p_PC_batch <- cbind(p_pca_df[,1:5],
                      c(rep("Groningen", times=10), rep("Cambridge", times=4), rep("Nijmegen", times=6)),
                      c(1032,1005,1032,rep(1005, times=7),1184,1184,1032,1032,788,0,788,80,0,0),
                      c(9,10,9,rep(11, times=5),15,11,38,38,38,20,14,4,14,15,17,18),
                      c(32.5,57.5,62.5,67.5,72.5,67.5,67.5,42.5,47.5,67.5,67.5,rep(62.5, times=3),rep(47.5, times=4),42.5,67.5),
                      c(rep("myeloid", times=10), rep("monocyte", times=10)), 
                      c(rep("Female", times=5),"Male",rep("Female", times=4),"Male","Female","Male","Female",rep("Male", times=5),"Female"))
  names(p_PC_batch)[6] <- "Provider"
  names(p_PC_batch)[7] <- "Date"
  names(p_PC_batch)[8] <- "Runs"
  names(p_PC_batch)[9] <- "Age"
  names(p_PC_batch)[10] <- "Cell_type"
  names(p_PC_batch)[11] <- "Sex"