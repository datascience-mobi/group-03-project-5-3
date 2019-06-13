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

# Conjoining datasets with different patients, with average coverage value for each gene. 
g_AML <- cbind(g_AMLpat.bed, g_AML.covmean)
g_Mono <- cbind(g_Monopat.bed, g_Mono.covmean)
p_patients <- data.frame(promoters_data_frame[11:50])

#Removing lines with more than 4 NAs from the beta values in both patient groups.
g_Monona <- data.frame(rowSums(is.na(g_Monopat.bed)))
g_AMLna <- data.frame(rowSums(is.na(g_AMLpat.bed)))
g_bednasum <- data.frame(cbind(genes_data_frame[,11:50], g_AMLna, g_Monona))
g_cleanna <- g_bednasum[!(g_bednasum$rowSums.is.na.g_AMLpat.bed.. > 4 | g_bednasum$rowSums.is.na.g_Monopat.bed.. > 4),]
g_cleanna <- g_cleanna[, c(1:40)]
p_cleanna <-  p_patients[!(rowSums(is.na(p_patients[1:10]) > 4)  |
                             rowSums(is.na(p_patients[11:20])) > 4 ) , ]
# Changing all beta values of genes with corresponding coverage value <30 and >95158 into NA
for (j in 21:40) {
  for (i in 1:nrow(g_cleanna)) {
    
    if(g_cleanna[i,j] < 30){
      g_cleanna[i,j-20] <- NA
    }
    if(g_cleanna[i,j] > 95158){
      g_cleanna[i,j-20] <- NA
    }}}
# Changing all beta values of promoterss with corresponding coverage value <30 and >14140 into NA
for (j in 21:40) {
  for (i in 1:nrow(p_cleanna)) {
    
    if(p_cleanna[i,j] < 30){
      p_cleanna[i,j-20] <- NA
    }
    if(p_cleanna[i,j] > 14140) {
      p_cleanna[i,j-20] <- NA
    }
  }}

# Cleaning up newly gained data frame from rows with more then 4 NA´s
promoters_clean <- p_cleanna[!(rowSums(is.na(p_cleanna[1:10])) > 4  |
                                 rowSums(is.na(p_cleanna[11:20])) > 4 ), ]
genes_clean <- g_patients[!(rowSums(is.na(g_patients[1:10])) > 4 | 
                              rowSums(is.na(g_patients[11:20])) > 4),  ]

# Converting beta values into M values for promoters
g_normalization <- genes_clean[,1:20]
g_normalization_backup <- g_normalization

for(j in 1:20){
  for (i in 1:nrow(g_normalization)){
    
    if(is.na(g_normalization[i, j])){
      
    }
    
    else{
      
      if(g_normalization[i, j] == 1){
        g_normalization[i, j] <- 0.9999999999
      }
      if(g_normalization[i, j] == 0) {
        g_normalization[i, j] <- 0.0000000001
      }
      
      g_normalization[i, j] = log2(g_normalization[i, j] / (1 - g_normalization[i, j]))
    }
  }}
g_Mvalues <- g_normalization
g_normalization <- g_normalization_backup
remove(g_normalization_backup, i, j)

#Converting beta values into M values for promoters

p_normalization <- promoters_clean[,1:20]
p_normalization_backup <- p_normalization

for(j in 1:20){
  for (i in 1:nrow(p_normalization)){
    
    if(is.na(p_normalization[i, j])){
      
    }
    
    else{
      
      if(p_normalization[i, j] == 1){
        p_normalization[i, j] <- 0.9999999999
      }
      if(p_normalization[i, j] == 0) {
        p_normalization[i, j] <- 0.0000000001
      }
      
      p_normalization[i, j] = log2(p_normalization[i, j] / (1 - p_normalization[i, j]))
    }
  }}
p_Mvalues <- p_normalization
p_normalization <- p_normalization_backup
remove(p_normalization_backup, i, j)

#Replacing all Na's with beta mean values of each row for g_Mvalues

for(j in 1:10){
  for(i in 1:nrow(g_Mvalues)){
    if(is.na(g_Mvalues[i, j])){
      g_Mvalues[i, j] <- rowMeans(g_Mvalues[i, 1:10], na.rm = TRUE)
    }}}

for(j in 11:20){
  for(i in 1:nrow(g_Mvalues)){
    if(is.na(g_Mvalues[i, j])){
      g_Mvalues[i, j] <- rowMeans(g_Mvalues[i, 11:20], na.rm = TRUE)
    }}}

#Replacing all Na's with beta mean values of each row for p_Mvalues

for(j in 1:10){
  for(i in 1:nrow(p_Mvalues)){
    if(is.na(p_Mvalues[i, j])){
      p_Mvalues[i, j] <- rowMeans(p_Mvalues[i, 1:10], na.rm = TRUE)
    }}}

for(j in 11:20){
  for(i in 1:nrow(p_Mvalues)){
    if(is.na(p_Mvalues[i, j])){
      p_Mvalues[i, j] <- rowMeans(p_Mvalues[i, 11:20], na.rm = TRUE)
    }}}
remove(i,j)
