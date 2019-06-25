# HYPERLOOP ANALYSIS FOR NA <= 3

###################################################################################################  
# Creating the appropriate Data frames basing on AML_Mono_list

  #initial formatting
  input_data <- AML_Mono_list
  genes_data_framexy <- input_data$genes
  promoters_data_framexy <- input_data$promoters
  
  # Removing chromosome X, because of hypermethylation. And chromosome Y because of lack of male patients.
  genes_data_frame <- data.frame(genes_data_framexy[!(genes_data_framexy$Chromosome == "chrX" | genes_data_framexy$Chromosome == "chrY"),])
  promoters_data_frame <- data.frame(promoters_data_framexy[!(promoters_data_framexy$Chromosome == "chrX" | promoters_data_framexy$Chromosome == "chrY"),])
  
  # creating data sets of all coverage values in the data set, for later sequence generation
  # All coverage values in one dataframe, for better plotting.
  g_coverage_all <- data.frame(Coverage = c(t(genes_data_frame[,31:50])))
  p_coverage_all <- data.frame(Coverage = c(t(promoters_data_frame[,31:50])))
  
  # grabbing some numbers for later tracking of what percentage of genes were cut in qc
  g_qcprecount <- nrow(genes_data_frame)
  p_qcprecount <- nrow(promoters_data_frame)
  
  # Creating datasets with only patients data.
  p_patients <- data.frame(promoters_data_frame[11:50])
  g_patients <- data.frame(genes_data_frame[11:50])
  

  # Changing all beta values of genes with corresponding coverage value <25 into NA (due to paper)
  for (j in 21:40) {
    for (i in 1:nrow(g_patients)) {
      
      if(g_patients[i,j] < 25){
        g_patients[i,j-20] <- NA
        }
      }
    }

    # Changing all beta values of promoterss with corresponding coverage value <25 into NA
  for (j in 21:40) {
    for (i in 1:nrow(p_patients)) {
      
      if(p_patients[i,j] < 25){
        p_patients[i,j-20] <- NA
      }
    }
  }
  
  # Cleaning up newly gained data frame from rows with more then 3 NA´s
  promoters_clean <- p_patients[!(rowSums(is.na(p_patients[1:10])) > 3  |
                                   rowSums(is.na(p_patients[11:20])) > 3 ), ]
  genes_clean <- g_patients[!(rowSums(is.na(g_patients[1:10])) > 3 | 
                                rowSums(is.na(g_patients[11:20])) > 3),  ]
  
  # removing uneccesary data sets
  remove(input_data, genes_data_framexy, promoters_data_framexy, genes_data_frame, promoters_data_frame, g_patients, p_patients, i, j)


  
###################################################################################################  
# performing two detailed hyperloop analyses
  
  #Creating a sequence for g hyperloop
  g_coverage_all_seq <- sort(g_coverage_all$Coverage, decreasing = TRUE)
  g_coverage_all_seq <- g_coverage_all_seq[1:200]
  g_coverage_all_seq <- sort(g_coverage_all_seq)
  
  g <- c(seq(25, min(g_coverage_all_seq), 100), g_coverage_all_seq)
  remove(g_coverage_all_seq)
  
  plot(g, type = "l")
  
  #Creating a sequence for p hyperloop
  p_coverage_all_seq <- sort(p_coverage_all$Coverage, decreasing = TRUE)
  p_coverage_all_seq <- p_coverage_all_seq[1:200]
  p_coverage_all_seq <- sort(p_coverage_all_seq)
  
  p <- c(seq(25, min(p_coverage_all_seq), 50), p_coverage_all_seq)
  remove(p_coverage_all_seq)
  
  plot(p, type = "l")
  
  # The hyperloop for finding the upper threshold for genes.
  
  analyse <- genes_clean
  
  for (j in 1:20) {
    for (i in 1:nrow(analyse)) {
      if (is.na(analyse[i, j]))
      {
        analyse[i, j + 20] <- NA
      }
    }
  }
  
  for (j in 1:20) {
    for (i in 1:nrow(analyse)) {
      if (is.na(analyse[i, j]))
      {
        analyse[i, j] <- -1
      }
    }
  }
  
  
  g_hyperloop <- analyse
  
  vector.remaining.genes.per.coverage <- 0
  
  nrow_g <- nrow(g_hyperloop)
  
  for (y in g) {
    for (j in 21:40) {
      for (i in 1:nrow(g_hyperloop)) {
        if (!(is.na(g_hyperloop[i, j]))) {
          if (g_hyperloop[i, j] < y)        {
            g_hyperloop[i, j - 20] <- NA
          }
        }
      }
    }
    print(nrow(analyse))
    g_hyperloop <-
      g_hyperloop[!(rowSums(is.na(g_hyperloop[1:10])) > 6 &
                      rowSums(is.na(g_hyperloop[11:20])) > 6) , ]
    nrow.cut.genes <- (nrow_g - nrow(g_hyperloop))
    print(nrow.cut.genes)
    vector.remaining.genes.per.coverage <-
      c(vector.remaining.genes.per.coverage, nrow.cut.genes)
  }
  remove(i, j, g_hyperloop, nrow_g, nrow.cut.genes, y, analyse)
  
  ## The hyperloop for finding the upper threshold for promoters.
  
  analyse <- promoters_clean
  
  for (j in 1:20) {
    for (i in 1:nrow(analyse)) {
      if (is.na(analyse[i, j]))
      {
        analyse[i, j + 20] <- NA
      }
    }
  }
  
  for (j in 1:20) {
    for (i in 1:nrow(analyse)) {
      if (is.na(analyse[i, j]))
      {
        analyse[i, j] <- -1
      }
    }
  }
  
  
  p_hyperloop <- analyse
  
  vector.remaining.promoters.per.coverage <- 0
  
  nrow_p <- nrow(p_hyperloop)
  
  for (y in p) {
    for (j in 21:40) {
      for (i in 1:nrow(p_hyperloop)) {
        if (!(is.na(p_hyperloop[i, j]))) {
          if (p_hyperloop[i, j] < y)        {
            p_hyperloop[i, j - 20] <- NA
          }
        }
      }
    }
    print(nrow(analyse))
    p_hyperloop <-
      p_hyperloop[!(rowSums(is.na(p_hyperloop[1:10])) > 6 &
                      rowSums(is.na(p_hyperloop[11:20])) > 6) , ]
    nrow.cut.promoters <- (nrow_p - nrow(p_hyperloop))
    print(nrow.cut.promoters)
    vector.remaining.promoters.per.coverage <-
      c(vector.remaining.promoters.per.coverage, nrow.cut.promoters)
  }
  
  remove(i, j, p_hyperloop, nrow_p, nrow.cut.promoters, y, analyse)

###################################################################################################  
#Doing Kneedle evaluatins for both hyperloop outputs
  #Doing a Kneedle evaluation of genes Hyperloop Output
      #Formatting
      
      g_HyperResults_NA3 <- data.frame(vector.remaining.genes.per.coverage)
      g_HyperResults_NA3 <- g_HyperResults_NA3[2:(nrow(g_HyperResults_NA3)),]
      g <- g[1:length(g_HyperResults_NA3)]
      g_HyperResults_NA3 <- cbind(g_HyperResults_NA3, g)
      g_HyperResults_NA3 <- data.frame(g_HyperResults_NA3)
      
      colnames(g_HyperResults_NA3) <- c("Remaining_Genes", "Coverage_Value")
      
      #Determining Kneelde m
      X2 <- max(g_HyperResults_NA3$Coverage_Value)
      X1 <- min(g_HyperResults_NA3$Coverage_Value)
      Y2 <- max(g_HyperResults_NA3$Remaining_Genes)
      Y1 <- min(g_HyperResults_NA3$Remaining_Genes) 
      
      kneedle.m <- ((Y2-Y1)/(X2-X1))
      
      #Determining Kneedle Line, Kneedle Line to remaining gene count difference and cbinding
      kneedle.line <- ((g_HyperResults_NA3$Coverage_Value-X1)*kneedle.m)
      g_HyperResults_NA3 <- cbind(g_HyperResults_NA3, kneedle.line)
      
      kneedle.difference <- (g_HyperResults_NA3$Remaining_Genes - g_HyperResults_NA3$kneedle.line)
      g_HyperResults_NA3 <- cbind(g_HyperResults_NA3, kneedle.difference)
      
      #Plotting
      plot(g_HyperResults_NA3$Coverage_Value, g_HyperResults_NA3$Remaining_Genes, type = "n")
      lines(g_HyperResults_NA3$Coverage_Value, g_HyperResults_NA3$Remaining_Genes)
      lines(g_HyperResults_NA3$Coverage_Value, g_HyperResults_NA3$kneedle.line)
      lines(g_HyperResults_NA3$Coverage_Value, g_HyperResults_NA3$kneedle.difference)
      
      remove(X1, X2, Y1, Y2, kneedle.difference, kneedle.line, g)
      
      # getting coverage at the knee
      g_uppercov_row <- which.max(g_HyperResults_NA3[, 4])
      g_uppercov <- g_HyperResults_NA3[g_uppercov_row,2]
      
      remove(g_uppercov_row)
  
  
  #Doing a Kneedle evaluation of promoters Hyperloop Output
      #Formatting
      p_HyperResults_NA3 <- data.frame(vector.remaining.promoters.per.coverage)
      p_HyperResults_NA3 <- p_HyperResults_NA3[2:(nrow(p_HyperResults_NA3)),]
      p <- p[1:length(p_HyperResults_NA3)]
      p_HyperResults_NA3 <- cbind(p_HyperResults_NA3, p)
      p_HyperResults_NA3 <- data.frame(p_HyperResults_NA3)
      
      colnames(p_HyperResults_NA3) <- c("Remaining_Promoters", "Coverage_Value")
      
      #Determining Kneelde m
      X2 <- max(p_HyperResults_NA3$Coverage_Value)
      X1 <- min(p_HyperResults_NA3$Coverage_Value)
      Y2 <- max(p_HyperResults_NA3$Remaining_Promoters)
      Y1 <- min(p_HyperResults_NA3$Remaining_Promoters) 
      
      kneedle.m <- ((Y2-Y1)/(X2-X1))
      
      #Determining Kneedle Line, Kneedle Line to remaining promoters count difference and cbinding
      kneedle.line <- ((p_HyperResults_NA3$Coverage_Value-X1)*kneedle.m)
      p_HyperResults_NA3 <- cbind(p_HyperResults_NA3, kneedle.line)
      
      kneedle.difference <- (p_HyperResults_NA3$Remaining_Promoters - p_HyperResults_NA3$kneedle.line)
      p_HyperResults_NA3 <- cbind(p_HyperResults_NA3, kneedle.difference)
      
      #Plotting
      plot(p_HyperResults_NA3$Coverage_Value, p_HyperResults_NA3$Remaining_Promoters, type = "n")
      lines(p_HyperResults_NA3$Coverage_Value, p_HyperResults_NA3$Remaining_Promoters)
      lines(p_HyperResults_NA3$Coverage_Value, p_HyperResults_NA3$kneedle.line)
      lines(p_HyperResults_NA3$Coverage_Value, p_HyperResults_NA3$kneedle.difference)
      
      remove(X1, X2, Y1, Y2, kneedle.difference, kneedle.line, p)
      
      # getting coverage at the knee
      p_uppercov_row <- which.max(p_HyperResults_NA3[, 4])
      p_uppercov <- p_HyperResults_NA3[p_uppercov_row,2]
      
      remove(p_uppercov_row)