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

# Creating datasets with only patients data.
p_patients <- data.frame(promoters_data_frame[11:50])
g_patients <- data.frame(genes_data_frame[11:50])

#Removing lines with more than 4 NAs from the beta values in both patient groups.
g_cleanna <- g_patients[!(rowSums(is.na(g_patients[1:10]) > 4)  |
                            rowSums(is.na(g_patients[11:20])) > 4 ) , ]
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
g_Mvalues <- genes_clean[,1:20]

for(j in 1:20){
  for (i in 1:nrow(g_Mvalues)){
    
    if(is.na(g_Mvalues[i, j])){
      
    }
    
    else{
      
      if(g_Mvalues[i, j] == 1){
        g_Mvalues[i, j] <- 0.9999999999
      }
      if(g_Mvalues[i, j] == 0) {
        g_Mvalues[i, j] <- 0.0000000001
      }
      
      g_Mvalues[i, j] = log2(g_Mvalues[i, j] / (1 - g_Mvalues[i, j]))
    }
  }}

remove(i, j)

#Converting beta values into M values for promoters

p_Mvalues <- promoters_clean[,1:20]

for(j in 1:20){
  for (i in 1:nrow(p_Mvalues)){
    
    if(is.na(p_Mvalues[i, j])){
      
    }
    
    else{
      
      if(p_Mvalues[i, j] == 1){
        p_Mvalues[i, j] <- 0.9999999999
      }
      if(p_Mvalues[i, j] == 0) {
        p_Mvalues[i, j] <- 0.0000000001
      }
      
      p_Mvalues[i, j] = log2(p_Mvalues[i, j] / (1 - p_Mvalues[i, j]))
    }
  }}

remove(i, j)

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
remove(i,j, cpgislands_data_frame, g_cleanna, p_cleanna, g_patients, p_patients, genes_clean, promoters_clean, genes_data_frame, genes_data_framexy, promoters_data_frame, promoters_data_framexy, tiling_data_frame, input_data)
