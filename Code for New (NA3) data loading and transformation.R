#initial formatting
input_data <- AML_Mono_list
genes_data_framexy <- input_data$genes
promoters_data_framexy <- input_data$promoters

# Removing chromosome X, because of hypermethylation. And chromosome Y because of lack of male patients.
genes_data_frame <- data.frame(genes_data_framexy[!(genes_data_framexy$Chromosome == "chrX" | genes_data_framexy$Chromosome == "chrY"),])
promoters_data_frame <- data.frame(promoters_data_framexy[!(promoters_data_framexy$Chromosome == "chrX" | promoters_data_framexy$Chromosome == "chrY"),])

# Creating datasets with only patients data.
g_patients <- data.frame(genes_data_frame[11:50])
p_patients <- data.frame(promoters_data_frame[11:50])

# Changing all beta values of genes with corresponding coverage value <25 into NA (due to paper) and > 98525 into NA (due to kneedle analysis)
for (j in 21:40) {
  for (i in 1:nrow(g_patients)) {
    
    if(g_patients[i,j] < 25 | g_patients[i,j] >= 98525){
      g_patients[i,j-20] <- NA
    }
  }
}

# Changing all beta values of genes with corresponding coverage value <25 into NA (due to paper) and > 14175 into NA (due to kneedle analysis)
for (j in 21:40) {
  for (i in 1:nrow(p_patients)) {
    
    if(p_patients[i,j] < 25 | p_patients[i,j] >= 14175){
      p_patients[i,j-20] <- NA
    }
  }
}

# Cleaning up newly gained data frame from rows with more then 3 NA´s
genes_clean <- g_patients[!(rowSums(is.na(g_patients[1:10])) > 3 | 
                              rowSums(is.na(g_patients[11:20])) > 3),  ]
promoters_clean <- p_patients[!(rowSums(is.na(p_patients[1:10])) > 3 | 
                              rowSums(is.na(p_patients[11:20])) > 3),  ]

# Converting beta values into M values for genes
genes_clean_reduced <- genes_clean[,1:20]
promoters_clean_reduced <- promoters_clean[,1:20]


f_BetaToM <- function(x) {
  
  if(!(is.na(x))){
    
    if(x == 1){
      x <- 0.9999999999
    }
    
    if(x == 0) {
      x <- 0.0000000001
    }
    
    x = log2(x / (1 - x))
  }
  
  return(x)
}

g_M_NA <- apply(genes_clean_reduced, c(1,2), f_BetaToM)
g_M_NA <- data.frame(g_M_NA)

p_M_NA <- apply(promoters_clean_reduced, c(1,2), f_BetaToM)
p_M_NA <- data.frame(p_M_NA)

#Replacing all Na's with rnorm10 imputation
fg_rnorm10_imputation <- function(x) {
  if(x %% 10000 == 0.0000000000000000) {
    print(x/10000)
  }
  workingrow <- g_M_NA[x, ]
  
  if(rowSums(is.na(workingrow) > 0)) {
    replacement_positions <- which(is.na(workingrow))
    rowvalues <- workingrow[!(is.na(workingrow))]
    
    winnercandidates <- c(rep(0, (1+length(replacement_positions))))
    
    for(i in 1:10) {
      set.seed(x + nrow(g_M_NA)*i)
      rowextend <- c(rnorm(mean = mean(rowvalues), sd = sd(rowvalues), length(replacement_positions)))
      rowvalues_ext <- c(rowvalues, rowextend)
      
      preextendstats <- c(mean(rowvalues), sd(rowvalues))
      postextendstats <- c(mean(rowvalues_ext), sd(rowvalues_ext))
      impdiff <- (postextendstats-preextendstats)
      impdiff <- (sum(abs(impdiff)))
      
      candidate <- c(impdiff, rowextend)
      winnercandidates <- cbind(winnercandidates, candidate)
      
    }
    
    winnercandidates <- data.frame(winnercandidates)
    winnercandidates <- winnercandidates[, 2:ncol(winnercandidates)]
    winnercol <- which.min(winnercandidates[1, ])
    winner <- winnercandidates[2:nrow(winnercandidates), winnercol]
    
    workingrow[replacement_positions] <- winner
  }
  return(workingrow)
}

fp_rnorm10_imputation <- function(x) {
  if(x %% 10000 == 0.0000000000000000) {
    print(x/10000)
  }
  workingrow <- p_M_NA[x, ]
  
  if(rowSums(is.na(workingrow) > 0)) {
    replacement_positions <- which(is.na(workingrow))
    rowvalues <- workingrow[!(is.na(workingrow))]
    
    winnercandidates <- c(rep(0, (1+length(replacement_positions))))
    
    for(i in 1:10) {
      set.seed(x + nrow(p_M_NA)*i)
      rowextend <- c(rnorm(mean = mean(rowvalues), sd = sd(rowvalues), length(replacement_positions)))
      rowvalues_ext <- c(rowvalues, rowextend)
      
      preextendstats <- c(mean(rowvalues), sd(rowvalues))
      postextendstats <- c(mean(rowvalues_ext), sd(rowvalues_ext))
      impdiff <- (postextendstats-preextendstats)
      impdiff <- (sum(abs(impdiff)))
      
      candidate <- c(impdiff, rowextend)
      winnercandidates <- cbind(winnercandidates, candidate)
      
    }
    
    winnercandidates <- data.frame(winnercandidates)
    winnercandidates <- winnercandidates[, 2:ncol(winnercandidates)]
    winnercol <- which.min(winnercandidates[1, ])
    winner <- winnercandidates[2:nrow(winnercandidates), winnercol]
    
    workingrow[replacement_positions] <- winner
  }
  return(workingrow)
}

j <- data.frame(seq(1, nrow(g_M_NA), 1))
g_Mvalues <- apply(j, c(1,2) , fg_rnorm10_imputation)
g_Mvalues <- matrix(unlist(g_Mvalues), ncol = nrow(g_M_NA), nrow = 20)
g_Mvalues <- data.frame(t(g_Mvalues))

j <- data.frame(seq(1, nrow(p_M_NA), 1))
p_Mvalues <- apply(j, c(1,2) , fp_rnorm10_imputation)
p_Mvalues <- matrix(unlist(p_Mvalues), ncol = nrow(p_M_NA), nrow = 20)
p_Mvalues <- data.frame(t(p_Mvalues))

# recreating colnames and rownames
colnames(g_Mvalues) <- colnames(g_M_NA)
rownames(g_Mvalues) <- rownames(g_M_NA)

colnames(p_Mvalues) <- colnames(p_M_NA)
rownames(p_Mvalues) <- rownames(p_M_NA)

remove(j, i, input_data, g_M_NA, g_patients, genes_clean, genes_clean_reduced, genes_data_frame, genes_data_framexy, p_M_NA, p_patients, promoters_clean, promoters_clean_reduced, promoters_data_frame, promoters_data_framexy)

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