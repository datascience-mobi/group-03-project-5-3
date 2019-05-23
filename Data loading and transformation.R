# Loading the data and naming the frames.
AML_Mono_list <- readRDS("~/Heidelberg/Uni/FS4 2019/Bioinfo/AMLvsMono/AML_Mono_list.RDS")
input_data <- AML_Mono_list
genes_data_frame <- input_data$genes
promoters_data_frame <- input_data$promoters
cpgislands_data_frame <- input_data$cpgislands
tiling_data_frame <- input_data$tiling
# Calculating the length of each gene.
genes_length <- data.frame(genes_data_frame$End - genes_data_frame$Start)
# Creating new matrixes with only one group of patients, either beta values of coverage.
g_AMLpat.bed <- genes_data_frame[,11:20]
g_Monopat.bed <- genes_data_frame[,21:30]
g_AMLpat.cov <- genes_data_frame[,31:40]
g_Monopat.cov <- genes_data_frame[,41:50]
# Calculating the mean coverage and beta values for AML and Mono patients for each gene.
# Might be handy later, but the NAs still have to be removed.
g_Mono.bedmean <- data.frame(rowMeans(g_Monopat.bed))
g_AML.bedmean <- data.frame(rowMeans(g_Monopat.bed))
g_AML.covmean <- data.frame(rowMeans(g_AMLpat.cov))
g_Mono.covmean <- data.frame(rowMeans(g_Monopat.cov))
# Conjoining datasets with different patients, with average coverage value for each gene. 
g_AML <- cbind(g_AMLpat.bed, g_AML.covmean)
g_Mono <- cbind(g_Monopat.bed, g_Mono.covmean)

