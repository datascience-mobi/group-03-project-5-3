# Distribution of all Coverage Values in log10: 
    plot(density(log10(coverage_all$Coverage)), main = "Distribution of coverage values")
# Drawing a "if we make this the maximum coverage value we keep this % of objects" plot: 
    P = seq(0,1, 0.001)
    Y = seq(0, 200000, 1)
    Q = ecdf(g_AML.covmean$rowMeans.g_AMLpat.cov.) (Y)
    plot(Y, Q, type = "n", main = "quantiles for coverage values")
    lines(Y, Q)
# Drawing the Kneedle algorithm inspired transformation, D is the difference between the kneedle slope line and the second plot
    slope_kneedle_line = max(coverage_all$Coverage)
    D <- (Q - (Y/slope_kneedle_line))
    lines(Y, D)
    abline(h = max(D), col= "red")
    knee_cut_cov = max(D)
    knee_cut_cov
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