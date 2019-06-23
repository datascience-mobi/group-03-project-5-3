p_PC_batch_promotors_leo <- read.csv("C:/Users/pierr/Documents/USB-Backup/Uni Heidelberg/MoBi_Bsc/Semester IV/Modul_Einführung in die Informatik/Modulelement2_Anwendung bioinformatischer Systeme/Projekt/Daten/bearbeites Set von Konsti/p_PC_batch_leo.csv")

#############################################################################################################################

# Kruskal-Wallis Test (for more than two categories; Provider (_prov))
g_pvalues_prov <- c()
for (i in 1:5) {
  g <- (kruskal.test(p_PC_batch_promotors_leo[,i] ~ p_PC_batch_promotors_leo[,6]))$p.value
  g_pvalues_prov <- c(g_pvalues_prov, g)
  remove(g)
}

#############################################################################################################################

# Wilcox Test (two logical categories; Sex (_sex) and Cell_type (_Ct)
for (j in 10:11) {
  g_pvalues_Ct_sex <- c()
  for (i in 1:5) {
    g <- (wilcox.test(p_PC_batch_promotors_leo[,i] ~ p_PC_batch_promotors_leo[,j]))$p.value
    g_pvalues_Ct_sex <- c(g_pvalues_Ct_sex, g)
    remove(g)
  }
  
  if(j ==10){
    g_pvalues_Ct <- g_pvalues_Ct_sex
  }
  
  else{
    g_pvalues_sex <- g_pvalues_Ct_sex
  }
  remove(g_pvalues_Ct_sex)
}

##########################################################################################################################

# Monte Carlo / Permutation test(for numerical categories; Date (_date), Runs (_runs), Age (_age))
# First installation of package surveillance is required

for (j in 7:9) {
  g_pvalues <- c()
  for (i in 1:5) {
    g <- (permutationTest(p_PC_batch_promotors_leo[,i], p_PC_batch_promotors_leo[,j], nPermutation = 9999))$pVal.permut
    g_pvalues <- c(g_pvalues, g)
    remove(g)
  }
  if(j ==7){
    g_pvalues_date <- g_pvalues
  }
  else{
    if(j == 8){
      g_pvalues_runs <- g_pvalues
    }
    else{
      g_pvalues_age <- g_pvalues
    }
  }
  remove(g_pvalues)
}

############################################################################################################################

# Reporting values in a data frame

pvalues_test1_pro_le <- data.frame(g_pvalues_prov, g_pvalues_date, g_pvalues_runs, g_pvalues_age, g_pvalues_Ct, g_pvalues_sex, row.names = c("PC1", "PC2", "PC3", "PC4", "PC5"))
colnames(pvalues_test1_pro_le) <- c("Provider", "Date", "Runs", "Age", "Cell type", "Sex")
pvalues_test_trans1_pro_le <- t(pvalues_test1_pro_le)

############################################################################################################################

# Creating Heatmap

norming_pvalues <- function(x){
  if(x<=0.01){
    x = 1
  }
  else{
    x=10
  }
}

pvalues_test_norm1_pro_le <- apply(pvalues_test_trans1_pro_le
                            , c(1,2)
                            , norming_pvalues)
remove(norming_pvalues)
####################################################################################################

# Heatmap with two colors
heatmap(as.matrix(pvalues_test_norm1_pro_le)
        , scale = "none"
        , col = heat.colors(2)
        , main = "K: P; Mon: D, R, A; Will: Ct, S; red < 0.01; [leo_promotors], 1"
        , Rowv = NA
        , Colv = NA
)

# Heatmap with continous colors

heatmap(as.matrix(pvalues_test_trans1_pro_le)
        , scale = "none"
        , col = heat.colors(256)
        , main = "K: P; Mon: D, R, A; Will: Ct, S; red < 0.01; [leo_promotors], 1"
        , Rowv = NA
        , Colv = NA
)
remove(i, j, g_pvalues_age, g_pvalues_date, g_pvalues_prov, g_pvalues_runs, g_pvalues_sex, g_pvalues_Ct)










#############################################################################################################################

#############################################################################################################################

# Kruskal-Wallis Test (for more than two logical categories; Provider (_prov), Date (_date))
for (j in 6:7) {
  g_pvalues_prov_date <- c()
  for (i in 1:5) {
    g <- (kruskal.test(p_PC_batch_promotors_leo[,i] ~ p_PC_batch_promotors_leo[,j]))$p.value
    g_pvalues_prov_date <- c(g_pvalues_prov_date, g)
    remove(g)
  }
  if(j == 6){
    g_pvalues_prov <- g_pvalues_prov_date
  }
  else{
    g_pvalues_date <- g_pvalues_prov_date
  }
  remove(g_pvalues_prov_date)
}

#############################################################################################################################

# Wilcox Test (two logical categories; Sex (_sex) and Cell_type (_Ct)
for (j in 10:11) {
  g_pvalues_Ct_sex <- c()
  for (i in 1:5) {
    g <- (wilcox.test(p_PC_batch_promotors_leo[,i] ~ p_PC_batch_promotors_leo[,j]))$p.value
    g_pvalues_Ct_sex <- c(g_pvalues_Ct_sex, g)
    remove(g)
  }
  
  if(j ==10){
    g_pvalues_Ct <- g_pvalues_Ct_sex
  }
  
  else{
    g_pvalues_sex <- g_pvalues_Ct_sex
  }
  remove(g_pvalues_Ct_sex)
}

##########################################################################################################################

# Monte Carlo / Permutation test(for numerical categories; Runs (_runs), Age (_age))
# First installation of package surveillance is required
library(surveillance)
for (j in 8:9) {
  g_pvalues <- c()
  for (i in 1:5) {
    g <- (permutationTest(p_PC_batch_promotors_leo[,i], p_PC_batch_promotors_leo[,j], nPermutation = 9999))$pVal.permut
    g_pvalues <- c(g_pvalues, g)
    remove(g)
  }
  if(j == 8){
    g_pvalues_runs <- g_pvalues
  }
  else{
    g_pvalues_age <- g_pvalues
  }
  remove(g_pvalues)
}

############################################################################################################################

# Reporting values in a data frame

pvalues_test2_pro_le <- data.frame(g_pvalues_prov, g_pvalues_date, g_pvalues_runs, g_pvalues_age, g_pvalues_Ct, g_pvalues_sex, row.names = c("PC1", "PC2", "PC3", "PC4", "PC5"))
colnames(pvalues_test2_pro_le) <- c("Provider", "Date", "Runs", "Age", "Cell type", "Sex")
pvalues_test_trans2_pro_le <- t(pvalues_test2_pro_le)

############################################################################################################################

# Creating Heatmap

norming_pvalues <- function(x){
  if(x<=0.01){
    x = 1
  }
  else{
    x=10
  }
}

pvalues_test_norm2_pro_le <- apply(pvalues_test_trans2_pro_le
                            , c(1,2)
                            , norming_pvalues)
remove(norming_pvalues)

##############################################################################################################

# Creating Heatmaps

# Heatmap with two colors
heatmap(as.matrix(pvalues_test_norm2_pro_le)
        , scale = "none"
        , col = heat.colors(2)
        , main = "K: P, D; Mon: R, A; Will: Ct, S; red < 0.01; [leo_promotors],2"
        , Rowv = NA
        , Colv = NA
)


# Heatmap with continous colors
heatmap(as.matrix(pvalues_test_trans2_pro_le)
        , scale = "none"
        , col = heat.colors(256)
        , main = "K: P, D; Mon: R, A; Will: Ct, S; red < 0.01; [leo_promotors],2"
        , Rowv = NA
        , Colv = NA
)

remove(i, j, g_pvalues_age, g_pvalues_date, g_pvalues_prov, g_pvalues_runs, g_pvalues_sex, g_pvalues_Ct)









####################################################################################################

####################################################################################################

# Kruskal-Wallis Test (for more than two logical categories; Provider (_prov), Date (_date), Age (_age))
for (j in 6:7) {
  g_pvalues_prov_date <- c()
  for (i in 1:5) {
    g <- (kruskal.test(p_PC_batch_promotors_leo[,i] ~ p_PC_batch_promotors_leo[,j]))$p.value
    g_pvalues_prov_date <- c(g_pvalues_prov_date, g)
    remove(g)
  }
  if(j == 6){
    g_pvalues_prov <- g_pvalues_prov_date
  }
  else{
    g_pvalues_date <- g_pvalues_prov_date
  }
  remove(g_pvalues_prov_date)
}

g_pvalues_age <- c()
for (i in 1:5) {
  g <- (kruskal.test(p_PC_batch_promotors_leo[,i] ~ p_PC_batch_promotors_leo[,9]))$p.value
  g_pvalues_age <- c(g_pvalues_age, g)
  remove(g)
}

#############################################################################################################################

# Wilcox Test (two logical categories; Sex (_sex) and Cell_type (_Ct)
for (j in 10:11) {
  g_pvalues_Ct_sex <- c()
  for (i in 1:5) {
    g <- (wilcox.test(p_PC_batch_promotors_leo[,i] ~ p_PC_batch_promotors_leo[,j]))$p.value
    g_pvalues_Ct_sex <- c(g_pvalues_Ct_sex, g)
    remove(g)
  }
  
  if(j ==10){
    g_pvalues_Ct <- g_pvalues_Ct_sex
  }
  
  else{
    g_pvalues_sex <- g_pvalues_Ct_sex
  }
  remove(g_pvalues_Ct_sex)
}

##########################################################################################################################

# Monte Carlo / Permutation test(for numerical categories; Runs (_runs))
# First installation of package surveillance is required
library(surveillance)

g_pvalues_runs <- c()
for (i in 1:5) {
  g <- (permutationTest(p_PC_batch_promotors_leo[,i], p_PC_batch_promotors_leo[,8], nPermutation = 9999))$pVal.permut
  g_pvalues_runs <- c(g_pvalues_runs, g)
  remove(g)
}

############################################################################################################################

# Reporting values in a data frame

pvalues_test3_pro_le <- data.frame(g_pvalues_prov, g_pvalues_date, g_pvalues_runs, g_pvalues_age, g_pvalues_Ct, g_pvalues_sex, row.names = c("PC1", "PC2", "PC3", "PC4", "PC5"))
colnames(pvalues_test3_pro_le) <- c("Provider", "Date", "Runs", "Age", "Cell type", "Sex")
pvalues_test_trans3_pro_le <- t(pvalues_test3_pro_le)

############################################################################################################################

# Creating Heatmap

norming_pvalues <- function(x){
  if(x<=0.01){
    x = 1
  }
  else{
    x=10
  }
}

pvalues_test_norm3_pro_le <- apply(pvalues_test_trans3_pro_le
                            , c(1,2)
                            , norming_pvalues)


##############################################################################################################

# Creating Heatmaps

# Heatmap with two colors
heatmap(as.matrix(pvalues_test_norm3_pro_le)
        , scale = "none"
        , col = heat.colors(2)
        , main = "K: P, D, A; Mon: R; Will: Ct, S; red < 0.01; [leo_promotors],3"
        , Rowv = NA
        , Colv = NA
)
remove(norming_pvalues)


# Heatmap with continous colors

heatmap(as.matrix(pvalues_test_trans3_pro_le)
        , scale = "none"
        , col = heat.colors(256)
        , main = "K: P, D, A; Mon: R; Will: Ct, S; red < 0.01; [leo_promotors],3"
        , Rowv = NA
        , Colv = NA
)

remove(i, j, g_pvalues_age, g_pvalues_date, g_pvalues_prov, g_pvalues_runs, g_pvalues_sex, g_pvalues_Ct)










########################################################################################################

########################################################################################################

# All kruskal

for (j in 6:11) {
  g_pvalues_prov_date_runs_age_ct_sex <- c()
  for (i in 1:5) {
    g <- (kruskal.test(p_PC_batch_promotors_leo[,i] ~ p_PC_batch_promotors_leo[,j]))$p.value
    g_pvalues_prov_date_runs_age_ct_sex <- c(g_pvalues_prov_date_runs_age_ct_sex, g)
    remove(g)
  }
  if(j == 6){
    g_pvalues_prov <- g_pvalues_prov_date_runs_age_ct_sex
  }
  else{
    if(j==7){
      g_pvalues_date <- g_pvalues_prov_date_runs_age_ct_sex
    }
    else{
      if(j == 8){
        g_pvalues_runs <- g_pvalues_prov_date_runs_age_ct_sex
      }
      else{
        if(j == 9){
          g_pvalues_age <- g_pvalues_prov_date_runs_age_ct_sex
        }
        else{
          if(j == 10){
            g_pvalues_Ct <- g_pvalues_prov_date_runs_age_ct_sex
          }
          else{
            g_pvalues_sex <- g_pvalues_prov_date_runs_age_ct_sex
          }
        }
      }
    }
    
  }
  remove(g_pvalues_prov_date_runs_age_ct_sex)
}

############################################################################################################################

# Reporting values in a data frame

pvalues_test4_pro_le <- data.frame(g_pvalues_prov, g_pvalues_date, g_pvalues_runs, g_pvalues_age, g_pvalues_Ct, g_pvalues_sex, row.names = c("PC1", "PC2", "PC3", "PC4", "PC5"))
colnames(pvalues_test4_pro_le) <- c("Provider", "Date", "Runs", "Age", "Cell type", "Sex")
pvalues_test_trans4_pro_le <- t(pvalues_test4_pro_le)

############################################################################################################################

# Creating Heatmap

norming_pvalues <- function(x){
  if(x<=0.01){
    x = 1
  }
  else{
    x=10
  }
}

pvalues_test_norm4_pro_le <- apply(pvalues_test_trans4_pro_le
                            , c(1,2)
                            , norming_pvalues)


##############################################################################################################

# Creating Heatmaps

# Heatmap with two colors
heatmap(as.matrix(pvalues_test_norm4_pro_le)
        , scale = "none"
        , col = heat.colors(2)
        , main = "K: P, D, R, A, Ct, S; red < 0.01; [leo_promotors],4"
        , Rowv = NA
        , Colv = NA
)
remove(norming_pvalues)


# Heatmap with continous colors

heatmap(as.matrix(pvalues_test_trans4_pro_le)
        , scale = "none"
        , col = heat.colors(256)
        , main = "K: P, D, R, A, Ct, S; red < 0.01; [leo_promotors],4"
        , Rowv = NA
        , Colv = NA
)
remove(i, j, g_pvalues_age, g_pvalues_date, g_pvalues_prov, g_pvalues_runs, g_pvalues_sex, g_pvalues_Ct)
