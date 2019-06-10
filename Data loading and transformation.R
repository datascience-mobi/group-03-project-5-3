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
g_cleanna <- g_cleanna[, c(1:40)]
p_cleanna <-  p_patients[!(rowSums(is.na(p_patients[1:10]) > 4)  |
                             rowSums(is.na(p_patients[11:20])) > 4 ) , ]
# Changing all beta values with corresponding coverage value <30 and >95158 into NA
for (j in 21:40) {
  for (i in 1:nrow(g_cleanna)) {
    
    if(g_cleanna[i,j] < 30){
      g_cleanna[i,j-20] <- NA
    }
    if(g_cleanna[i,j] > 95158){
      g_cleanna[i,j-20] <- NA
    }}}
# Removing coverage threshold of 30 by promoters.
for (j in 21:40) {
  for (i in 1:nrow(p_cleanna)) {
    
    if(p_cleanna[i,j] < 30){
      p_cleanna[i,j-20] <- NA
    }
    if(p_cleanna[i,j] > 14140) {
      p_cleanna[i,j-20] <- NA
    }
  }}

promoters.clean <- p_cleanna[!(rowSums(is.na(p_cleanna[1:10])) > 4  |
                                 rowSums(is.na(p_cleanna[11:20])) > 4 ), ]


# Cleaning up newly gained Z data frame from rows with more then 4 NA´s
genes.clean <- g_cleanna[!(rowSums(is.na(g_cleanna[1:10])) > 4 | rowSums(is.na(g_cleanna[11:20])) > 4),  ]

#Converting beta values into M values for promoters

g.normalization <- genes.clean[,1:20]
g.normalization.backup <- g.normalization

for(j in 1:20){
  for (i in 1:nrow(g.normalization)){
    
    if(is.na(g.normalization[i, j])){
      
    }
    
    else{
      
      if(g.normalization[i, j] == 1){
        g.normalization[i, j] <- 0.9999999999
      }
      if(g.normalization[i, j] == 0) {
        g.normalization[i, j] <- 0.0000000001
      }
      
      g.normalization[i, j] = log2(g.normalization[i, j] / (1 - g.normalization[i, j]))
    }
  }}
g.Mvalues <- g.normalization
g.normalization <- g.normalization.backup
remove(g.normalization.backup, i, j)

#Converting beta values into M values for promoters

p.normalization <- promoters.clean[,1:20]
p.normalization.backup <- p.normalization

for(j in 1:20){
  for (i in 1:nrow(p.normalization)){
    
    if(is.na(p.normalization[i, j])){
      
    }
    
    else{
      
      if(p.normalization[i, j] == 1){
        p.normalization[i, j] <- 0.9999999999
      }
      if(p.normalization[i, j] == 0) {
        p.normalization[i, j] <- 0.0000000001
      }
      
      p.normalization[i, j] = log2(p.normalization[i, j] / (1 - p.normalization[i, j]))
    }
  }}
p.Mvalues <- p.normalization
p.normalization <- p.normalization.backup
remove(p.normalization.backup, i, j)

#Replacing all Na's with beta mean values of each row for g.Mvalues

for(j in 1:10){
  for(i in 1:nrow(g.Mvalues)){
    if(is.na(g.Mvalues[i, j])){
      g.Mvalues[i, j] <- rowMeans(g.Mvalues[i, 1:10], na.rm = TRUE)
    }}}

for(j in 11:20){
  for(i in 1:nrow(g.Mvalues)){
    if(is.na(g.Mvalues[i, j])){
      g.Mvalues[i, j] <- rowMeans(g.Mvalues[i, 11:20], na.rm = TRUE)
    }}}

#Replacing all Na's with beta mean values of each row for p.Mvalues

for(j in 1:10){
  for(i in 1:nrow(p.Mvalues)){
    if(is.na(p.Mvalues[i, j])){
      p.Mvalues[i, j] <- rowMeans(p.Mvalues[i, 1:10], na.rm = TRUE)
    }}}

for(j in 11:20){
  for(i in 1:nrow(p.Mvalues)){
    if(is.na(p.Mvalues[i, j])){
      p.Mvalues[i, j] <- rowMeans(p.Mvalues[i, 11:20], na.rm = TRUE)
    }}}
remove(i,j)