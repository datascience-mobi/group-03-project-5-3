###################################################################################################
# (1) initial formatting
input_data <- AML_Mono_list
genes_data_framexy <- input_data$genes
promoters_data_framexy <- input_data$promoters

# Removing chromosome X, because of hypermethylation. And chromosome Y because of lack of male patients.
genes_data_frame <- data.frame(genes_data_framexy[!(genes_data_framexy$Chromosome == "chrX" | genes_data_framexy$Chromosome == "chrY"),])
promoters_data_frame <- data.frame(promoters_data_framexy[!(promoters_data_framexy$Chromosome == "chrX" | promoters_data_framexy$Chromosome == "chrY"),])

# Creating datasets with only patients data.
g_patients <- data.frame(genes_data_frame[11:50])
p_patients <- data.frame(promoters_data_frame[11:50])

remove(input_data, genes_data_framexy, promoters_data_framexy)

###################################################################################################
# (2) Cleaning out unreliable rows with unreliable coverage and too many NAs
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

remove(g_patients, p_patients, i, j)

###################################################################################################
# (3) M value transformation and imputation  
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
  
  if(x %% jcount == 0) {
    cat("|")
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
  
  if(x %% jcount == 0) {
    cat("|")
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
jcount <- floor(nrow(j)/100)

g_Mvalues <- apply(j, c(1,2) , fg_rnorm10_imputation)
g_Mvalues <- matrix(unlist(g_Mvalues), ncol = nrow(g_M_NA), nrow = 20)
g_Mvalues <- data.frame(t(g_Mvalues))

j <- data.frame(seq(1, nrow(p_M_NA), 1))
jcount <- floor(nrow(j)/100)

p_Mvalues <- apply(j, c(1,2) , fp_rnorm10_imputation)
p_Mvalues <- matrix(unlist(p_Mvalues), ncol = nrow(p_M_NA), nrow = 20)
p_Mvalues <- data.frame(t(p_Mvalues))

# recreating colnames and rownames
colnames(g_Mvalues) <- colnames(g_M_NA)
rownames(g_Mvalues) <- rownames(g_M_NA)

colnames(p_Mvalues) <- colnames(p_M_NA)
rownames(p_Mvalues) <- rownames(p_M_NA)

remove(j, jcount, genes_clean, genes_clean_reduced, promoters_clean, promoters_clean_reduced)

###################################################################################################
# (4) g_T and p_T formatting
#splitting datasets to AML and Mon
g_M_NA_AML <- g_M_NA[, 1:10]
g_M_NA_mon <- g_M_NA[, 11:20]

p_M_NA_AML <- p_M_NA[, 1:10]
p_M_NA_mon <- p_M_NA[, 11:20]

#defining a function that wil determine the values required for t-tests later
f_MtoT <- function(x) {
  
  
  m <- mean(x, na.rm = TRUE)
  sd  <- sd(x, na.rm = TRUE)
  n <- sum(!(is.na(x)))
  
  return(c(m,sd,n))
}

# application of function to genes and promoters
g_T_AML <- t(apply(g_M_NA_AML, 1, f_MtoT))
g_T_mon <- t(apply(g_M_NA_mon, 1, f_MtoT))

p_T_AML <- t(apply(p_M_NA_AML, 1, f_MtoT))
p_T_mon <- t(apply(p_M_NA_mon, 1, f_MtoT))

#naming, formatting, rremoving unwanted data sets
g_T <- cbind(g_T_AML, g_T_mon)
colnames(g_T) <- c("Mean AML", "SD AML", "N AML", "Mean Mono", "SD Mono", "N Mono")
g_T <- data.frame(g_T)

p_T <- cbind(p_T_AML, p_T_mon)
colnames(p_T) <- c("Mean AML", "SD AML", "N AML", "Mean Mono", "SD Mono", "N Mono")
p_T <- data.frame(p_T)

remove(g_T_AML, g_T_mon)
remove(p_T_AML, p_T_mon)
remove(f_BetaToM, fg_rnorm10_imputation, fp_rnorm10_imputation, f_MtoT)

###################################################################################################
# (5) generating the resource data set

# collecting the mean of the m values, as per the m distribution, the mean represents the most common i.e. normal value
# (from each cohort in the promoters and genes data set)
g_RM_AML_m <- data.frame(apply(g_M_NA_AML, 1, mean, na.rm = TRUE))
g_RM_mon_m <- data.frame(apply(g_M_NA_mon, 1, mean, na.rm = TRUE))

p_RM_AML_m <- data.frame(apply(p_M_NA_AML, 1, mean, na.rm = TRUE))
p_RM_mon_m <- data.frame(apply(p_M_NA_mon, 1, mean, na.rm = TRUE))


# defining a function that converts the "normal" M value for each sequence's cohort to the corresponding "normal" beta value
f_MtoBeta <- function(x) {
  
  x = (2^x)/(1 + 2^x)
  return(x)
  
}

# oddities of the data set resulting from application of the above funtion: if the original gene was 1 for all samples in beta, the result here is 1, even though we changed those 1s to 0.999...9, so r has rounded here.
# however, genes that were 0 in bet are returned as 0.00...01 here, because we also changed them along the way. R has not rounded here. Very strange but coincidentally convenient for future formatting
# (otherwise we wouldve had to set 0.9...9 to 1, and leave 0.0...1 because we can't divide by zero )

g_RB_AML_m <- apply(g_RM_AML_m, c(1,2), f_MtoBeta)
g_RB_mon_m <- apply(g_RM_mon_m, c(1,2), f_MtoBeta)
p_RB_AML_m <- apply(p_RM_AML_m, c(1,2), f_MtoBeta)
p_RB_mon_m <- apply(p_RM_mon_m, c(1,2), f_MtoBeta)

# calculating log2 foldchange from Mon to AML (so that the "normal" variety / 0 fold change is the healthy one)
g_RB_foldchange <- log2(g_RB_AML_m/g_RB_mon_m)
p_RB_foldchange <- log2(p_RB_AML_m/p_RB_mon_m)

# extracting symbols for later from genes and promoters
g_R_Symbols <- data.frame(genes_data_frame$symbol)
rownames(g_R_Symbols) <- rownames(genes_data_frame)

p_R_Symbols <- data.frame(promoters_data_frame$symbol)
rownames(p_R_Symbols) <- rownames(promoters_data_frame)

# formatting the g_resource data set
g_Resource <- merge(g_RB_foldchange, g_R_Symbols , by = 0, all = FALSE)
g_Resource <- g_Resource[, c(1,3,2)]
rownames(g_Resource) <- g_Resource$Row.names
colnames(g_Resource) <- c("Ensign_ID", "Symbols", "Foldchange_Beta")

# formatting the p_resource data set
p_Resource <- merge(p_RB_foldchange, p_R_Symbols , by = 0, all = FALSE)
p_Resource <- p_Resource[, c(1,3,2)]
rownames(p_Resource) <- p_Resource$Row.names
colnames(p_Resource) <- c("Ensign_ID", "Symbols", "Foldchange_Beta")

# removal of uneccesary data sets
remove(g_M_NA_AML, g_M_NA_mon, g_R_Symbols, g_RB_AML_m, g_RB_mon_m, g_RB_foldchange, g_RM_AML_m, g_RM_mon_m, genes_data_frame, promoters_data_frame, f_MtoBeta)
remove(p_M_NA_AML, p_M_NA_mon, p_R_Symbols, p_RB_AML_m, p_RB_mon_m, p_RB_foldchange, p_RM_AML_m, p_RM_mon_m)
