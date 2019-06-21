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
