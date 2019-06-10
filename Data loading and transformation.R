# Loading the data and naming the frames.
AML_Mono_list <- readRDS("~/Heidelberg/Uni/FS4 2019/Bioinfo/AMLvsMono/AML_Mono_list.RDS")
input_data <- AML_Mono_list
genes_data_framexy <- input_data$genes
promoters_data_framexy <- input_data$promoters
cpgislands_data_frame <- input_data$cpgislands
tiling_data_frame <- input_data$tiling
# Removing chromosome X, because of hypermethylation. And chromosome Y because of lack of male patients.
genes_data_frame <- data.frame(genes_data_framexy[!(genes_data_framexy$Chromosome == "chrX" | genes_data_framexy$Chromosome == "chrY"),])
promoters_data_frame <- data.frame(promoters_data_framexy[!(promoters_data_framexy$Chromosome == "chrX" | promoters_data_framexy$Chromosome == "chrY"),])
# Creating new matrixes with only one group of patients, either beta values or coverage.
g_AMLpat.bed <- genes_data_frame[,11:20]
g_Monopat.bed <- genes_data_frame[,21:30]
g_AMLpat.cov <- genes_data_frame[,31:40]
g_Monopat.cov <- genes_data_frame[,41:50]

# All coverage values in one dataframe, for better plotting.
g_covall <- genes_data_frame[,31:50]
p_covall <- promoters_data_frame[,31:50]
g_coverage_all <- data.frame(Coverage = c(t(g_covall)))
p_coverage_all <- data.frame(Coverage = c(t(p_covall)))
# Calculating the mean coverage and beta values for AML and Mono patients for each gene.
# Might be handy later, but the NAs still have to be handled.
g_Mono.bedmean <- data.frame(rowMeans(g_Monopat.bed))
g_AML.bedmean <- data.frame(rowMeans(g_Monopat.bed))
g_AML.covmean <- data.frame(rowMeans(g_AMLpat.cov))
g_Mono.covmean <- data.frame(rowMeans(g_Monopat.cov))
# Conjoining datasets with different patients, with average coverage value for each gene. 
g_AML <- cbind(g_AMLpat.bed, g_AML.covmean)
g_Mono <- cbind(g_Monopat.bed, g_Mono.covmean)
p_patients <- data.frame(promoters_data_frame[11:50])

#Removing lines with more than 4 NAs from the beta values in both patient groups.
g_Monona <- data.frame(rowSums(is.na(g_Monopat.bed)))
g_AMLna <- data.frame(rowSums(is.na(g_AMLpat.bed)))
g_bednasum <- data.frame(cbind(genes_data_frame[,11:50], g_AMLna, g_Monona))
g_cleanna <- g_bednasum[!(g_bednasum$rowSums.is.na.g_AMLpat.bed.. > 4 | g_bednasum$rowSums.is.na.g_Monopat.bed.. > 4),]
Z <- g_cleanna[, c(1:40)]
p_cleanna <-  p_patients[!(rowSums(is.na(p_patients[1:10]) > 4)  |
                             rowSums(is.na(p_patients[11:20])) > 4 ) , ]
# Changing all beta values with corresponding coverage value <30 into NA
for (j in 21:40) {
  for (i in 1:nrow(Z)) {
    
    if(Z[i,j] < 30){
      Z[i,j-20] <- NA
    }}}
# Removing coverage threshold of 30 by promoters.
for (j in 21:40) {
  for (i in 1:nrow(p_cleanna)) {
    
    if(p_cleanna[i,j] < 30){
      p_cleanna[i,j-20] <- NA
    }}}

Z_promoters.clean <- p_cleanna[!(rowSums(is.na(p_cleanna[1:10]) > 4)  |
                                   rowSums(is.na(p_cleanna[11:20])) > 4 ) , ]

remove(i, j)

# Cleaning up newly gained Z data frame from rows with more then 4 NA´s
Z.nasum <- data.frame(cbind(Z, rowSums(is.na(Z[,1:10])), rowSums(is.na(Z[,11:20]))))
Z_clean <- Z.nasum[!(Z.nasum$rowSums.is.na.Z...1.10... > 4 | Z.nasum$rowSums.is.na.Z...11.20... > 4),]
Z_clean <- Z_clean[,1:40]
remove(Z.nasum, Z, i, j)