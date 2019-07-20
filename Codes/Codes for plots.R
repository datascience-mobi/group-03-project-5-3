# Distribution of all Coverage Values in log10: 
    plot(density(log10(coverage_all$Coverage)), main = "Distribution of coverage values")

# All coverage values in one dataframe, for better plotting.
    g_coverage_all <- data.frame(Coverage = c(t(genes_data_frame[,31:50])))
    p_coverage_all <- data.frame(Coverage = c(t(promoters_data_frame[,31:50])))
    
# Creating new matrixes with only one group of patients with beta values withoun NAs. 
    g_AML_bedall <- data.frame(Beta = c(t(na.omit(genes_data_frame[,11:20]))))
    g_Mono_bedall <- data.frame(Beta = c(t(na.omit(genes_data_frame[,21:30]))))
    p_AML_bedall <- data.frame(Beta = c(t(na.omit(promoters_data_frame[,11:20]))))
    p_Mono_bedall <- data.frame(Beta = c(t(na.omit(promoters_data_frame[,21:30]))))
    
# Distribution of all beta values by AML or Mono patiens in the promoters and genes data set.
    par(mfrow=c(1, 2))
    hist(g_Mono_bedall$Beta, main = "Distribution of beta values by Mono patients for genes", xlab = "beta value", col = "forestgreen", breaks = 50)
    hist(g_AML_bedall$Beta, main = "Distribution of beta values by AML patients for genes", xlab = "beta value", col = "red", breaks = 50)
    par(mfrow=c(1, 2))
    hist(p_Mono_bedall$Beta, main = "Distribution of beta values by Mono patients for promoters", xlab = "beta value", col = "forestgreen", breaks = 50)
    hist(p_AML_bedall$Beta, main = "Distribution of beta values by AML patients for genes", xlab = "beta value", col = "red", breaks = 50)
    
# Drawing an "if we make this the maximum coverage value we keep this % of objects" plot for promoters and genes: 
    # Genes
    Y = seq(0, 200000, 1)
    Q = ecdf(g_coverage_all$Coverage) (Y)
    plot(Y, Q, type = "n", main = "quantiles for coverage by genes")
    lines(Y, Q)
  # For promoters:   
    Y = seq(0, 200000, 1)
    Q = ecdf(p_coverage_all$Coverage) (Y)
    plot(Y, Q, type = "n", main = "quantiles for coverage by promoters")
    lines(Y, Q)
    
# Printing how many genes will remain in the data set if only all values above knee_cut_cov are cut according to similar criteria as were used in "Data loading and transformation" (that code must be executed first)
    for (j in 21:40) {
      for (i in 1:nrow(Z)) {
        
        if(Z[i,j] > knee_cut_cov){
          Z[i,j-20] <- NA
        }}}
    
    Z.nasum <- data.frame(cbind(Z, rowSums(is.na(Z[,1:10])), rowSums(is.na(Z[,11:20]))))
    Z_clean <- Z.nasum[!(Z.nasum$rowSums.is.na.Z...1.10... > 4 | Z.nasum$rowSums.is.na.Z...11.20... > 4),]
    Z_clean <- Z_clean[,1:40]
    remove(Z.nasum, Z, i, j)
    
    percent_genes_remain = 100*nrow(Z_clean)/nrow(genes_data_frame)
    print(percent_genes_remain)
    if(percent_genes_remain >= 0.85) {
      print(paste("yay!"))
    }
    
# The hyperloop for finding the upper threshold for genes.
    analyse <- Z_clean
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
    
    
    a <- seq(42000, 48000, 60)
    Z.for.counting.rows <- analyse
    vector.remaining.genes.per.coverage <- 0
    nrow_Z <- nrow(Z.for.counting.rows)
    for (y in a) {
      for (j in 21:40) {
        for (i in 1:nrow(Z.for.counting.rows)) {
          if (!(is.na(Z.for.counting.rows[i, j]))) {
            if (Z.for.counting.rows[i, j] < y)        {
              Z.for.counting.rows[i, j - 20] <- NA
            }
          }
        }
      }
      print(nrow(analyse))
      Z.for.counting.rows <-
        Z.for.counting.rows[!(rowSums(is.na(Z.for.counting.rows[1:10])) > 5 &
                                rowSums(is.na(Z.for.counting.rows[11:20])) > 5) , ]
      nrow.cut.genes <- (nrow_Z - nrow(Z.for.counting.rows))
      print(nrow.cut.genes)
      vector.remaining.genes.per.coverage <-
        c(vector.remaining.genes.per.coverage, nrow.cut.genes)
    }
    
#Doing a Kneedle evaluation of a Hyperloop Output
    #Determining Kneedle m
    HyperloopResults_Number <- data.frame(vector.remaining.genes.per.coverage)
    
    HyperloopResults_Number <- HyperloopResults_Number[2:(nrow(HyperloopResults_Number)),]
    HyperloopResults_Number <- cbind(HyperloopResults_Number, a)
    HyperloopResults_Number <- data.frame(HyperloopResults_Number)
    
    names(HyperloopResults_Number)[names(HyperloopResults_Number) == "HyperloopResults_Number"] <- "Remaining.count.of.genes"
    names(HyperloopResults_Number)[names(HyperloopResults_Number) == "a"] <- "Coverage.value"
    HyperloopResults_NumberBackup <- HyperloopResults_Number
    
    
    X2 <- max(HyperloopResults_Number$Coverage.value)
    X1 <- min(HyperloopResults_Number$Coverage.value)
    Y2 <- max(HyperloopResults_Number$Remaining.count.of.genes)
    Y1 <- min(HyperloopResults_Number$Remaining.count.of.genes) #CAREFUL, you may need to shorten HyperloopResults_Number to the last increment at which all genes are included
    
    kneedle.m <- ((Y2-Y1)/(X2-X1))
    
    #Determining Kneedle Line, Kneedle Line to remaining gene count difference and cbinding
    kneedle.line <- ((HyperloopResults_Number$Coverage.value-X1)*kneedle.m)
    HyperloopResults_Number <- cbind(HyperloopResults_Number, kneedle.line)
    
    kneedle.difference <- (HyperloopResults_Number$Remaining.count.of.genes - HyperloopResults_Number$kneedle.line)
    HyperloopResults_Number <- cbind(HyperloopResults_Number, kneedle.difference)
    
    #Plotting
    plot(HyperloopResults_Number$Coverage.value, HyperloopResults_Number$Remaining.count.of.genes, type = "n")
    lines(HyperloopResults_Number$Coverage.value, HyperloopResults_Number$Remaining.count.of.genes)
    lines(HyperloopResults_Number$Coverage.value, HyperloopResults_Number$kneedle.line)
    lines(HyperloopResults_Number$Coverage.value, HyperloopResults_Number$kneedle.difference)
    #Make sure to manually check the maximum of HyperloopResults_Number$kneedle.difference and check the corresponding coverage value!
    
#Doing a Kneedle Evaluation of a HyRes Hyperloop Output for already exisitng Hyperloop Results
    # To use this code, make sure to replace "_Number" with something to differentiate it from the original Hyperloop Results
    HyperloopResults_Number <- data.frame(vector.remaining.genes.per.coverage)
    
    HyperloopResults_Number <- HyperloopResults_Number[2:(nrow(HyperloopResults_Number)),]
    HyperloopResults_Number <- cbind(HyperloopResults_Number, a)
    HyperloopResults_Number <- data.frame(HyperloopResults_Number)
    
    names(HyperloopResults_Number)[names(HyperloopResults_Number) == "HyperloopResults_Number"] <- "Remaining.count.of.genes"
    names(HyperloopResults_Number)[names(HyperloopResults_Number) == "a"] <- "Coverage.value"
    HyperloopResults_NumberBackup <- HyperloopResults_Number
    
    #Determining Kneedle Line, Kneedle Line to remaining gene count difference and cbinding
    kneedle.line <- ((HyperloopResults_Number$Coverage.value-X1)*kneedle.m)
    HyperloopResults_Number <- cbind(HyperloopResults_Number, kneedle.line)
    
    kneedle.difference <- (HyperloopResults_Number$Remaining.count.of.genes - HyperloopResults_Number$kneedle.line)
    HyperloopResults_Number <- cbind(HyperloopResults_Number, kneedle.difference)
    
    #Plotting, make sure to plot in the range of a. xlim = c(min(a),max(a)), similar for ylim
    plot(HyperloopResults_Number$Coverage.value, HyperloopResults_Number$kneedle.difference, type = "n")
    lines(HyperloopResults_Number$Coverage.value, HyperloopResults_Number$kneedle.difference)
    #Make sure to manually check the maximum of HyperloopResults_Number$kneedle.difference and check the corresponding coverage value!
    #Feel free to plot the old kneedle difference too by copyng the line above this one and changing the referenced Data Frame
    
# Calculating the percentage of lost genes per coverage value - Gigaloop
    Z_clean.backup <- Z_clean
    vector.remaining.genes.per.coverage <- 0
    for (y in c(70000, 71000, 80000)) {
      for (j in 21:40) {
        for (i in 1:nrow(Z_clean)) {
          
          if(Z_clean[i,j] > y){
            Z_clean[i,j-20] <- NA
          }}}
      print(nrow(Z_clean))
      Z.nasum <- data.frame(cbind(Z_clean, rowSums(is.na(Z_clean[,1:10])), rowSums(is.na(Z_clean[,11:20]))))
      Z.for.counting.rows <- Z.nasum[!(Z.nasum$rowSums.is.na.Z_clean...1.10... > 4 | Z.nasum$rowSums.is.na.Z_clean...11.20... > 4),]
      print(nrow(Z.for.counting.rows))
      vector.remaining.genes.per.coverage <- c(vector.remaining.genes.per.coverage, nrow(Z.for.counting.rows))
      remove(Z.for.counting.rows, Z.nasum)
      Z_clean <- Z_clean.backup
    } 
    remove(y, i, j)
    vector.remaining.genes.per.coverage <- vector.remaining.genes.per.coverage[-(1:1)]
    vector.remaining.genes.percentage.coverage <- vector.remaining.genes.per.coverage/nrow(Z_clean)
    # Naming vector values
    remaining.genes.corresponding.to.coverage <- cbind(c(70000, 71000, 80000), vector.remaining.genes.percentage.coverage, vector.remaining.genes.per.coverage)
    colnames(remaining.genes.corresponding.to.coverage) <- c("Coverage value", "Remaining percentage genes", "Remaining count of genes")
    View(remaining.genes.corresponding.to.coverage)
    