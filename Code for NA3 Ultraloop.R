#recquires coverage all and patients recquired
###################################################################################################
#Creating a sequence for g Ultraloop
g_coverage_all_seq <- sort(g_coverage_all$Coverage, decreasing = TRUE)
g_coverage_all_seq <- g_coverage_all_seq[1:200]
g_coverage_all_seq <- sort(g_coverage_all_seq)

g <- c(seq(25, min(g_coverage_all_seq), 100), g_coverage_all_seq)
g_count <- floor(length(g)/100)
remove(g_coverage_all_seq)

plot(g, type = "l")

#Creating a sequence for p Ultraloop
p_coverage_all_seq <- sort(p_coverage_all$Coverage, decreasing = TRUE)
p_coverage_all_seq <- p_coverage_all_seq[1:200]
p_coverage_all_seq <- sort(p_coverage_all_seq)

p <- c(seq(25, min(p_coverage_all_seq), 50), p_coverage_all_seq)
p_count <- floor(length(p)/100)
remove(p_coverage_all_seq)

plot(p, type = "l")

###################################################################################################
#g_ULTRALOOP
###################################################################################################
# creating the analyse data set and defining variables
analyse <- g_patients
g_lowercov <- 25

#cleaning out rows with low coverage
  # separating data set
  analyse_bet <- analyse[ , 1:20]
  analyse_cov <- analyse[, 21:40]
  
  #using the ifelse functions to record which positions in coverage should be NA and setting those positions in beta to NA
  cov_NAs <- ifelse(analyse_cov < g_lowercov, 1, 0)
  analyse_bet[cov_NAs == 1] <- NA
  
  # recombining and cleaning data set
  analyse <- cbind(analyse_bet, analyse_cov)
  analyse <- analyse[!(rowSums(is.na(analyse[1:10])) > 3 | 
                         rowSums(is.na(analyse[11:20])) > 3),  ]

#setting NAs to -1 and coverage NAs
  # separating data set again
  analyse_bet <- analyse[ , 1:20]
  analyse_cov <- analyse[, 21:40]
  
  #using ifelse to set NAs in bet to -1 and set NAs in bet to NA in coverage
  bet_NAs <- ifelse(is.na(analyse_bet), 1, 0 )
  analyse_bet[bet_NAs == 1] <- -1
  analyse_cov[bet_NAs == 1] <- NA

  # recombining and cleaning data set
  analyse <- cbind(analyse_bet, analyse_cov)

# removing uneccesary data sets
remove(analyse_cov, analyse_bet, cov_NAs, bet_NAs)

###################################################################################################
#Ultraloop Code for genes
g_Ultraloop <- analyse
vector_remaining_genes_per_coverage <- c()
nrow_g <- nrow(g_Ultraloop)

for (y in g) {
  
  if(which(g == y) %% g_count == 0) {
    cat("|")
  }
  
  #defining separate dat sets and uppercov
  g_Ultraloop_bet <- g_Ultraloop[ , 1:20]
  g_Ultraloop_cov <- g_Ultraloop[ , 21:40]
  g_uppercov <- y
  
  #using the ifelse functions to record which positions in coverage are trustworthy and turning them into NA
  cov_NAs <- ifelse(g_Ultraloop_cov < g_uppercov, 1, 0 )
  g_Ultraloop_bet[cov_NAs == 1] <- NA
  
  #recombining the data set and cleaning out of trustworthy rows (>6 NAs), because these will be trustworthy on the next run through too
  g_Ultraloop <- cbind(g_Ultraloop_bet, g_Ultraloop_cov)
  g_Ultraloop <-
    g_Ultraloop[!(rowSums(is.na(g_Ultraloop[1:10])) > 6 &
                    rowSums(is.na(g_Ultraloop[11:20])) > 6) , ]

  #calculating from the deficit of rows how many would've remained and storing them
  g_cut <- (nrow_g - nrow(g_Ultraloop))
  vector_remaining_genes_per_coverage <-
    c(vector_remaining_genes_per_coverage, g_cut)
  
  if(g_cut == nrow_g) {
    break()
  }
}

#removing uneccesary data sets
remove(y, g_Ultraloop, nrow_g, analyse, cov_NAs, g_cut, g_Ultraloop_bet, g_Ultraloop_cov, g_lowercov, g_uppercov)

###################################################################################################
#p_ULTRALOOP
###################################################################################################
# creating the analyse data set and defining variables
analyse <- p_patients
p_lowercov <- 25

# Cleaning out rows with low coverage
# separating data set
analyse_bet <- analyse[ , 1:20]
analyse_cov <- analyse[, 21:40]

# Using the ifelse functions to record which positions in coverage should be NA and setting those positions in beta to NA
cov_NAs <- ifelse(analyse_cov < p_lowercov, 1, 0)
analyse_bet[cov_NAs == 1] <- NA

# Recombining and cleaning data set
analyse <- cbind(analyse_bet, analyse_cov)
analyse <- analyse[!(rowSums(is.na(analyse[1:10])) > 3 | 
                       rowSums(is.na(analyse[11:20])) > 3),  ]

# Setting NAs to -1 and coverage NAs
# Separating data set again
analyse_bet <- analyse[ , 1:20]
analyse_cov <- analyse[, 21:40]

# Using ifelse to set NAs in bet to -1 and set NAs in bet to NA in coverage
bet_NAs <- ifelse(is.na(analyse_bet), 1, 0 )
analyse_bet[bet_NAs == 1] <- -1
analyse_cov[bet_NAs == 1] <- NA

# Recombining and cleaning data set
analyse <- cbind(analyse_bet, analyse_cov)

# Removing uneccesary data sets
remove(analyse_cov, analyse_bet, cov_NAs, bet_NAs)

###################################################################################################
#Ultraloop Code for promoters
p_Ultraloop <- analyse
vector_remaining_promoters_per_coverage <- c()
nrow_p <- nrow(p_Ultraloop)

for (y in p) {
  
  if(which(p == y) %% p_count == 0) {
    cat("|")
  }
  
  # Defining separate dat sets and uppercov
  p_Ultraloop_bet <- p_Ultraloop[ , 1:20]
  p_Ultraloop_cov <- p_Ultraloop[ , 21:40]
  p_uppercov <- y
  
  # Using the if-else functions to record which positions in coverage are trustworthy and turning them into NA
  cov_NAs <- ifelse(p_Ultraloop_cov < p_uppercov, 1, 0 )
  p_Ultraloop_bet[cov_NAs == 1] <- NA
  
  # Recombining the data set and cleaning out of trustworthy rows (>6 NAs), because these will be trustworthy on the next run through too
  p_Ultraloop <- cbind(p_Ultraloop_bet, p_Ultraloop_cov)
  p_Ultraloop <-
    p_Ultraloop[!(rowSums(is.na(p_Ultraloop[1:10])) > 6 &
                    rowSums(is.na(p_Ultraloop[11:20])) > 6) , ]
  
  #calculating from the deficit of rows how many would've remained and storing them
  p_cut <- (nrow_p - nrow(p_Ultraloop))
  vector_remaining_promoters_per_coverage <-
    c(vector_remaining_promoters_per_coverage, p_cut)
  
  if(p_cut == nrow_p) {
    break()
  }
}

#removing uneccesary data sets
remove(y, p_Ultraloop, nrow_p, analyse, cov_NAs, p_cut, p_Ultraloop_bet, p_Ultraloop_cov, p_lowercov, p_uppercov)

###################################################################################################
#Doing Kneedle evaluatins for both Ultraloop outputs
#Doing a Kneedle evaluation of genes Ultraloop Output
#Formatting

g_UltraResults_NA3 <- data.frame(vector_remaining_genes_per_coverage)
g <- g[1:nrow(g_UltraResults_NA3)]
g_UltraResults_NA3 <- cbind(g_UltraResults_NA3, g)
g_UltraResults_NA3 <- data.frame(g_UltraResults_NA3)

colnames(g_UltraResults_NA3) <- c("Remaining_Genes", "Coverage_Value")

#Determining Kneelde m
X2 <- max(g_UltraResults_NA3$Coverage_Value)
X1 <- min(g_UltraResults_NA3$Coverage_Value)
Y2 <- max(g_UltraResults_NA3$Remaining_Genes)
Y1 <- min(g_UltraResults_NA3$Remaining_Genes) 

kneedle.m <- ((Y2-Y1)/(X2-X1))

#Determining Kneedle Line, Kneedle Line to remaining gene count difference and cbinding
kneedle.line <- ((g_UltraResults_NA3$Coverage_Value-X1)*kneedle.m)
g_UltraResults_NA3 <- cbind(g_UltraResults_NA3, kneedle.line)

kneedle.difference <- (g_UltraResults_NA3$Remaining_Genes - g_UltraResults_NA3$kneedle.line)
g_UltraResults_NA3 <- cbind(g_UltraResults_NA3, kneedle.difference)

#Plotting
plot(g_UltraResults_NA3$Coverage_Value, g_UltraResults_NA3$Remaining_Genes, type = "n")
lines(g_UltraResults_NA3$Coverage_Value, g_UltraResults_NA3$Remaining_Genes)
lines(g_UltraResults_NA3$Coverage_Value, g_UltraResults_NA3$kneedle.line)
lines(g_UltraResults_NA3$Coverage_Value, g_UltraResults_NA3$kneedle.difference)
#Make sure to manually check the maximum of g_UltraResults_NA3$kneedle.difference and check the corresponding Coverage_Value!

remove(X1, X2, Y1, Y2, kneedle.difference, kneedle.line, kneedle.m, g)

# getting coverage at the knee
g_uppercov_row <- which.max(g_UltraResults_NA3[, 4])
g_uppercov <- g_UltraResults_NA3[g_uppercov_row,2]

remove(g_uppercov_row, g_count, kneedle.m)

###################################################################################################
#Doing a Kneedle evaluation of promoters Ultraloop Output
#Formatting
p_UltraResults_NA3 <- data.frame(vector_remaining_promoters_per_coverage)
p <- p[1:nrow(p_UltraResults_NA3)]
p_UltraResults_NA3 <- cbind(p_UltraResults_NA3, p)
p_UltraResults_NA3 <- data.frame(p_UltraResults_NA3)

colnames(p_UltraResults_NA3) <- c("Remaining_Promoters", "Coverage_Value")

#Determining Kneelde m
X2 <- max(p_UltraResults_NA3$Coverage_Value)
X1 <- min(p_UltraResults_NA3$Coverage_Value)
Y2 <- max(p_UltraResults_NA3$Remaining_Promoters)
Y1 <- min(p_UltraResults_NA3$Remaining_Promoters) 

kneedle.m <- ((Y2-Y1)/(X2-X1))

#Determining Kneedle Line, Kneedle Line to remaining promoters count difference and cbinding
kneedle.line <- ((p_UltraResults_NA3$Coverage_Value-X1)*kneedle.m)
p_UltraResults_NA3 <- cbind(p_UltraResults_NA3, kneedle.line)

kneedle.difference <- (p_UltraResults_NA3$Remaining_Promoters - p_UltraResults_NA3$kneedle.line)
p_UltraResults_NA3 <- cbind(p_UltraResults_NA3, kneedle.difference)

#Plotting
plot(p_UltraResults_NA3$Coverage_Value, p_UltraResults_NA3$Remaining_Promoters, type = "n")
lines(p_UltraResults_NA3$Coverage_Value, p_UltraResults_NA3$Remaining_Promoters)
lines(p_UltraResults_NA3$Coverage_Value, p_UltraResults_NA3$kneedle.line)
lines(p_UltraResults_NA3$Coverage_Value, p_UltraResults_NA3$kneedle.difference)
#Make sure to manually check the maximum of p_UltraResults_NA3$kneedle.difference and check the corresponding Coverage_Value!

remove(X1, X2, Y1, Y2, kneedle.difference, kneedle.line, p)

# getting coverage at the knee
p_uppercov_row <- which.max(p_UltraResults_NA3[, 4])
p_uppercov <- p_UltraResults_NA3[p_uppercov_row,2]

remove(p_uppercov_row, p_count)