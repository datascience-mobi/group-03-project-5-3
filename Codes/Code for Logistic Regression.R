# logistical regression
  # removing redundant items from genes data frame
  lg_Finale_reduced <- g_Finale[!(rownames(g_Finale) %in% rownames(p_Finale)), ]
  
  # extracting mvalues corresponding to the finale esen´mbl ids
  lg_Mvalues <- g_Mvalues[rownames(g_Mvalues) %in% rownames(lg_Finale_reduced), ]
  lp_Mvalues <- p_Mvalues[rownames(p_Mvalues) %in% rownames(p_Finale), ]
  
  # formatting and combining pMvalues and gMvalues
  lg_Mvalues <- cbind(lg_Mvalues, lg_Finale_reduced$p_values, lg_Finale_reduced$Foldchange_Beta)
  lp_Mvalues <- cbind(lp_Mvalues, p_Finale$p_values, p_Finale$Foldchange_Beta)
  colnames(lg_Mvalues)[(ncol(lg_Mvalues)-1):ncol(lg_Mvalues)] <- c("p_values", "Foldchange_Beta")
  colnames(lp_Mvalues)[(ncol(lp_Mvalues)-1):ncol(lp_Mvalues)] <- c("p_values", "Foldchange_Beta")
  
  l_Finale_Compound <- rbind(lg_Mvalues, lp_Mvalues)
  
  # ranking the data set
  l_sig <- rank(l_Finale_Compound$p_values)
  l_fc_abs <- abs(l_Finale_Compound$Foldchange_Beta)
  l_fc <- rank(-l_fc_abs)
  l_rank <- max(l_sig, l_fc)
  l_Finale_Compound <- l_Finale_Compound[order(l_rank), 1:20]
  
  # removing uneccesary datasets
  remove(l_sig, l_fc_abs, l_fc, l_rank, lg_Finale_reduced, lp_Mvalues, lg_Mvalues)
  
  # adding row for phentype
  l_Finale_Compound <- rbind(c(rep(0, 20)), l_Finale_Compound)
  rownames(l_Finale_Compound)[1] <- c("Phenotype")
  
  l_Finale_Compound <- data.frame(t(l_Finale_Compound))
  l_Finale_Compound[, 1] <- c(rep(1, 10), rep(0, 10))
  
  # function for creating a model on random permutations of test and training data sets
  fl_test <- function(x) {
    set.seed(x)
    l_trainseq <- c(sample(seq(1, 10, 1), 7), sample(seq(11, 20, 1), 7))
    l_testseq <- seq(1,20,1)
    l_testseq <- l_testseq[!(l_testseq %in% l_trainseq)]
    
    l_Finale_Compound_train <- l_Finale_Compound[l_trainseq, ]
    l_Finale_Compound_test <- l_Finale_Compound[l_testseq, ]
    
    l_Model <- glm(Phenotype ~., family = binomial(link = 'logit'), data = l_Finale_Compound_test)
    
    l_results <- predict(l_Model,newdata=l_Finale_Compound_test,type='response')
    l_results <- ifelse(l_results > 0.5, 1, 0)
    l_accuracy <- mean(l_results == c(1, 1, 1, 0, 0, 0))
    return(l_accuracy)
  }
  
  j <- data.frame(seq(1, 20, 1))
  l_results <- apply(j, 1, fl_test)
  print(paste(mean(l_results)))
  
  # generating one model for further analysis
  set.seed(123)
  l_trainseq <- c(sample(seq(1, 10, 1), 7), sample(seq(11, 20, 1), 7))
  l_testseq <- seq(1,20,1)
  l_testseq <- l_testseq[!(l_testseq %in% l_trainseq)]
  
  l_Finale_Compound_train <- l_Finale_Compound[l_trainseq, ]
  l_Finale_Compound_test <- l_Finale_Compound[l_testseq, ]
  
  l_Model <- glm(Phenotype ~., family = binomial(link = 'logit'), data = l_Finale_Compound_test)
