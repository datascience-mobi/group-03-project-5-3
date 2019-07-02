###################################################################################################
# loading packages
library(qvalue)
###################################################################################################
# creating p values via manually calculated student t-tests
  # defining a general function that applies a student t-test based on given m1/2, s1/2 and n1/2 
  f_ttest_alt <- function(m1, m2, s1, s2, n1, n2) {
    
    if(m1 == m2) {
      p_value = 1
      
    } else {
      
      df = (n1 + n2 - 2)
      s12 = sqrt((((n1 - 1)*s1^2)+((n2 - 1)*s2^2))/df)
      
      t = ((m1 - m2)/(s12*sqrt((1/n1) + (1/n2))))
      
      p_value = 2*pt(-abs(t),df)
      
    }
  }

  #defining appplication functions that will comb through g_t and p_t respectively and calculate t-tests for the given inputs
  fg_ttest_apl <- function(x) {
    m1 <- g_T[x, 1]
    m2 <- g_T[x, 4]
    
    s1 <- g_T[x, 2]
    s2 <- g_T[x, 5]
    
    n1 <- g_T[x, 3]
    n2 <- g_T[x, 6]
    
    p_value <- f_ttest_alt(m1, m2, s1, s2, n1, n2)
    return(p_value)
  }
  
  fp_ttest_apl <- function(x) {
    m1 <- p_T[x, 1]
    m2 <- p_T[x, 4]
    
    s1 <- p_T[x, 2]
    s2 <- p_T[x, 5]
    
    n1 <- p_T[x, 3]
    n2 <- p_T[x, 6]
    
    p_value <- f_ttest_alt(m1, m2, s1, s2, n1, n2)
    return(p_value)
  }

  # applying our application functions to g_t and p_t
  j <- data.frame(seq(1, nrow(g_T), 1))
  g_p_values <- apply(j, 1, fg_ttest_apl)
  
  j <- data.frame(seq(1, nrow(p_T), 1))
  p_p_values <- apply(j, 1, fp_ttest_apl)
  
  # formatting
  g_T <- cbind(g_T, g_p_values)
  colnames(g_T) <- c(colnames(g_T)[1:(ncol(g_T)-1)], "p_values")
  
  p_T <- cbind(p_T, p_p_values)
  colnames(p_T) <- c(colnames(p_T)[1:(ncol(p_T)-1)], "p_values")
  
  # removal of uneccesary data
  remove(g_p_values, fg_ttest_apl, f_ttest_alt, j)
  remove(p_p_values, fp_ttest_apl)
  
  # checking to make sure the p value distribution is well behaved
  hist(g_T$p_values, 100)
  hist(p_T$p_values, 100)

###################################################################################################
# q values 
  # ordering data sets so that the q value output can be cbound directly to the ordered data set
  g_T <- g_T[order(g_T$p_values, decreasing = TRUE), ]
  p_T <- p_T[order(p_T$p_values, decreasing = TRUE), ]
   
  # creating the q-objects, that cointain all the q value statistics and are saved for later
  g_q <- qvalue(p = g_T$p_values)
  p_q <- qvalue(p = p_T$p_values)
  
  # formatting
  g_T_q <- cbind(g_T, g_q$qvalue)
  colnames(g_T_q)[ncol(g_T_q)] <- c("q_values")
  
  p_T_q <- cbind(p_T, p_q$qvalue)
  colnames(p_T_q)[ncol(p_T_q)] <- c("q_values")

###################################################################################################
# extracting genes deemed relevant after PCA analysis from the g_T data frame
  g_T_reduced <- g_T_q[which(rownames(g_T_q) %in% g_PC1_names), ]
  p_T_reduced <- p_T_q[which(rownames(p_T_q) %in% p_PC1_names), ]
###################################################################################################
# final fromatting and volcano plot elements
  # merging our reduced data sets with the resource data sets to have symbol and searchable ensign id for every sequence
  g_Finale <- merge(g_T_reduced, g_Resource, by = 0, all = FALSE)
  rownames(g_Finale) <- g_Finale$Row.names
  g_Finale <- g_Finale[, c(8, 9, 12, 11, 10)]
  
  p_Finale <- merge(p_T_reduced, p_Resource, by = 0, all = FALSE)
  rownames(p_Finale) <- p_Finale$Row.names
  p_Finale <- p_Finale[, c(8, 9, 12, 11, 10)]
  
  remove(g_T_reduced, p_T_reduced)
  
  # selection for high fold changes and significance
  g_Finale <- g_Finale[g_Finale$q_values <= 0.05, ]
  g_Finale <- g_Finale[g_Finale$Foldchange_Beta < log2(3/4) | g_Finale$Foldchange_Beta > log2(4/3), ]
  
  p_Finale <- p_Finale[p_Finale$q_values <= 0.05, ]
  p_Finale <- p_Finale[p_Finale$Foldchange_Beta < log2(3/4) | p_Finale$Foldchange_Beta > log2(4/3), ]
  
  #ordering Finale data sets for fold change and significance
  g_sig <- rank(g_Finale$p_values)
  g_fc_abs <- abs(g_Finale$Foldchange_Beta)
  g_fc <- rank(-g_fc_abs)
  g_rank <- (g_fc + g_sig)/2
  g_Finale <- g_Finale[order(-g_rank, decreasing = TRUE), ]
  
  p_sig <- rank(p_Finale$p_values)
  p_fc_abs <- abs(p_Finale$Foldchange_Beta)
  p_fc <- rank(-p_fc_abs)
  p_rank <- (p_fc + p_sig)/2
  p_Finale <- p_Finale[order(-p_rank, decreasing = TRUE), ]
  
  # removal of uneccesary data sets
  remove(g_T_q, g_PC1_names)
  remove(p_T_q, p_PC1_names)
  remove(g_sig, g_fc_abs, g_fc, g_rank)
  remove(p_sig, p_fc_abs, p_fc, p_rank)
  
  # final comparative results and plotting for genes
  g_FinaleComparison <- merge(g_T, g_Resource, by = 0, all = FALSE)
  rownames(g_FinaleComparison) <- g_FinaleComparison$Row.names
  g_FinaleComparison <- g_FinaleComparison[, c(8,11, 10, 9)]
  
  smoothScatter(g_FinaleComparison$Foldchange_Beta, -log10(g_FinaleComparison$p_values), nbin = 1000, colramp = colorRampPalette(c("white", "red")))
  abline(h = -log10(max(g_Finale$p_values)))
  abline(v = log2(3/4))
  abline(v = log2(4/3))
  
  plot(g_FinaleComparison$Foldchange_Beta, -log10(g_FinaleComparison$p_values), type = "n")
  points(g_Finale$Foldchange_Beta, -log10(g_Finale$p_values))
  abline(h = -log10(max(g_Finale$p_values)))
  
  # final comparative results and plotting for promoters
  p_FinaleComparison <- merge(p_T, p_Resource, by = 0, all = FALSE)
  rownames(p_FinaleComparison) <- p_FinaleComparison$Row.names
  p_FinaleComparison <- p_FinaleComparison[, c(8,11, 10, 9)]
  
  smoothScatter(p_FinaleComparison$Foldchange_Beta, -log10(p_FinaleComparison$p_values), nbin = 1000, colramp = colorRampPalette(c("white", "red")))
  abline(h = -log10(max(p_Finale$p_values)))
  abline(v = log2(3/4))
  abline(v = log2(4/3))
  
  plot(p_FinaleComparison$Foldchange_Beta, -log10(p_FinaleComparison$p_values), type = "n")
  points(p_Finale$Foldchange_Beta, -log10(p_Finale$p_values))
  abline(h = -log10(max(p_Finale$p_values)))
  
  remove(g_T, g_Resource, g_M_NA)
  remove(p_T, p_Resource, p_M_NA)

