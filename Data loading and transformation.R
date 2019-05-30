# Loading the data and naming the frames.
AML_Mono_list <- readRDS("~/Heidelberg/Uni/FS4 2019/Bioinfo/AMLvsMono/AML_Mono_list.RDS")
input_data <- AML_Mono_list
genes_data_framexy <- input_data$genes
promoters_data_framexy <- input_data$promoters
cpgislands_data_frame <- input_data$cpgislands
tiling_data_frame <- input_data$tiling
# Removing chromosome X, because of hypermethylation. And chromosome Y because of lack of male patients.
genes_data_framey <- data.frame(genes_data_framexy[!(genes_data_framexy$Chromosome == "chrX"),])
promoters_data_framey <- data.frame(promoters_data_framexy[!(promoters_data_framexy$Chromosome == "chrX"),])
genes_data_frame <- data.frame(genes_data_framey[!(genes_data_framey$Chromosome == "chrY"),])
promoters_data_frame <- data.frame(promoters_data_framey[!(promoters_data_framey$Chromosome == "chrY"),])
# Creating new matrixes with only one group of patients, either beta values or coverage.
g_AMLpat.bed <- genes_data_frame[,11:20]
g_Monopat.bed <- genes_data_frame[,21:30]
g_AMLpat.cov <- genes_data_frame[,31:40]
g_Monopat.cov <- genes_data_frame[,41:50]
# All coverage values in one dataframe for better plotting.
g_covall <- genes_data_frame[,31:50]
coverage_all <- data.frame(Coverage = c(t(g_covall)))
# Calculating the mean coverage and beta values for AML and Mono patients for each gene.
# Might be handy later, but the NAs still have to be handled.
g_Mono.bedmean <- data.frame(rowMeans(g_Monopat.bed))
g_AML.bedmean <- data.frame(rowMeans(g_Monopat.bed))
g_AML.covmean <- data.frame(rowMeans(g_AMLpat.cov))
g_Mono.covmean <- data.frame(rowMeans(g_Monopat.cov))
# Conjoining datasets with different patients, with average coverage value for each gene. 
g_AML <- cbind(g_AMLpat.bed, g_AML.covmean)
g_Mono <- cbind(g_Monopat.bed, g_Mono.covmean)
#Removing lines with more than 4 NAs from the beta values in both patient groups.
g_Monona <- data.frame(rowSums(is.na(g_Monopat.bed)))
g_AMLna <- data.frame(rowSums(is.na(g_AMLpat.bed)))
g_bednasum <- data.frame(cbind(genes_data_frame[,11:50], g_AMLna, g_Monona))
g_cleanna <- g_bednasum[!(g_bednasum$rowSums.is.na.g_AMLpat.bed.. > 4 & g_bednasum$rowSums.is.na.g_Monopat.bed.. > 4),]
Z <- g_cleanna[, c(1:40)]
