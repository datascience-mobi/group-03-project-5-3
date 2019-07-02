######################################################################################################

######################################################################################################

g_PC_batch_genes_leo <- read.csv("C:/Users/pierr/Documents/USB-Backup/Uni Heidelberg/MoBi_Bsc/Semester IV/Modul_Einführung in die Informatik/Modulelement2_Anwendung bioinformatischer Systeme/Projekt/Daten/bearbeites Set von Konsti/g_PC_batch_genes_leo.csv")

# New updated code with the different patients cohorts and both cohorts combined
# Kruskal-Wallis for Provider (_provider), Date (_date) and Age (_age)

g_pvalues_prov <- c()
g_pvalues_date <- c()
g_pvalues_age <- c()

for (i in 1:5) {
  p <- (kruskal.test(g_PC_batch_genes_leo[,i]~g_PC_batch_genes_leo[,6]))$p.value
  d <- (kruskal.test(g_PC_batch_genes_leo[,i]~g_PC_batch_genes_leo[,7]))$p.value
  a <- (kruskal.test(g_PC_batch_genes_leo[,i]~g_PC_batch_genes_leo[,9]))$p.value
  g_pvalues_prov <- append(g_pvalues_prov, p)
  g_pvalues_date <- append(g_pvalues_date, d)
  g_pvalues_age <- append(g_pvalues_age, a)
}

######################################################################################################

# Monte carlo for runs (_runs)
# First installation of package surveillance is required
library(surveillance)

g_pvalues_runs <- c()
for (i in 1:5) {
  g <- (permutationTest(g_PC_batch_genes_leo[,i], g_PC_batch_genes_leo[,8], nPermutation = 9999))$pVal.permut
  g_pvalues_runs <- append(g_pvalues_runs, g)
}

######################################################################################################

# Willcoxon for cell type (_CT) and sex (_sex)

for (j in 10:11) {
  g_pvalues_Ct_sex <- c()
  for (i in 1:5) {
    g <- (wilcox.test(g_PC_batch_genes_leo[,i] ~ g_PC_batch_genes_leo[,j]))$p.value
    g_pvalues_Ct_sex <- append(g_pvalues_Ct_sex, g)
  }
  
  if(j ==10){
    g_pvalues_Ct <- g_pvalues_Ct_sex
  }
  
  else{
    g_pvalues_sex <- g_pvalues_Ct_sex
  }
  remove(g_pvalues_Ct_sex)
}

######################################################################################################

# Reporting values in a data frame

pvalues_genes_both <- data.frame(g_pvalues_prov, g_pvalues_date, g_pvalues_runs, g_pvalues_age, g_pvalues_Ct, g_pvalues_sex, row.names = c("PC1", "PC2", "PC3", "PC4", "PC5"))
colnames(pvalues_genes_both) <- c("Provider", "Date", "Runs", "Age", "Cell type", "Sex")
pvalues_genes_both_trans <- t(pvalues_genes_both)

#######################################################################################################

# Norming data frame for two colored heatmap

norming_pvalues <- function(x){
  if(x<=0.01){
    x = 1
  }
  else{
    x=10
  }
}

pvalues_genes_both_norm<- apply(pvalues_genes_both_trans
                                   , c(1,2)
                                   , norming_pvalues)

#####################################################################################################

# Creating heatmaps

# Heatmap with two colors
heatmap(as.matrix(pvalues_genes_both_norm)
        , scale = "none"
        , col = heat.colors(2)
        , main = "K: P, D, A; Mon: R; Will: Ct, S; red < 0.01; [leo_genes_both]"
        , Rowv = NA
        , Colv = NA
)


# Heatmap with continous colors

heatmap(as.matrix(pvalues_genes_both_trans)
        , scale = "none"
        , col = heat.colors(256)
        , main = "K: P, D, A; Mon: R; Will: Ct, S; red < 0.01; [leo_genes_both]"
        , Rowv = NA
        , Colv = NA
)
remove(a,d,p,i,j,g,g_pvalues_age, g_pvalues_Ct, g_pvalues_date, g_pvalues_prov, g_pvalues_runs, g_pvalues_sex, norming_pvalues)



##################################################################################################################

#########################################################################################################

# Plots for the two cohorts alone

# AML
g_PC_batch_genes_leo_AML <- g_PC_batch_genes_leo[1:10,]

# New updated code with only the AML cohort
# Kruskal-Wallis for Date (_date) and Age (_age)

g_pvalues_date <- c()
g_pvalues_age <- c()

for (i in 1:5) {
  d <- (kruskal.test(g_PC_batch_genes_leo_AML[,i]~g_PC_batch_genes_leo_AML[,7]))$p.value
  a <- (kruskal.test(g_PC_batch_genes_leo_AML[,i]~g_PC_batch_genes_leo_AML[,9]))$p.value
  g_pvalues_date <- append(g_pvalues_date, d)
  g_pvalues_age <- append(g_pvalues_age, a)
}

######################################################################################################

# Monte carlo for runs (_runs)
# First installation of package surveillance is required
library(surveillance)

g_pvalues_runs <- c()
for (i in 1:5) {
  g <- (permutationTest(g_PC_batch_genes_leo_AML[,i], g_PC_batch_genes_leo_AML[,8], nPermutation = 9999))$pVal.permut
  g_pvalues_runs <- append(g_pvalues_runs, g)
}

######################################################################################################

# Willcoxon for sex (_sex)

g_pvalues_sex <- c()
  for (i in 1:5) {
    g <- (wilcox.test(g_PC_batch_genes_leo_AML[,i] ~ g_PC_batch_genes_leo_AML[,11]))$p.value
    g_pvalues_sex <- append(g_pvalues_sex, g)
}
  

######################################################################################################

# Reporting values in a data frame

pvalues_genes_AML <- data.frame(g_pvalues_date, g_pvalues_runs, g_pvalues_age, g_pvalues_sex, row.names = c("PC1", "PC2", "PC3", "PC4", "PC5"))
colnames(pvalues_genes_AML) <- c("Date", "Runs", "Age", "Sex")
pvalues_genes_AML_trans <- t(pvalues_genes_AML)

#######################################################################################################

# Norming data frame for two colored heatmap

norming_pvalues <- function(x){
  if(x<=0.01){
    x = 1
  }
  else{
    x=10
  }
}

pvalues_genes_AML_norm<- apply(pvalues_genes_AML_trans
                                , c(1,2)
                                , norming_pvalues)

#####################################################################################################

# Creating heatmaps

# Heatmap with two colors
heatmap(as.matrix(pvalues_genes_AML_norm)
        , scale = "none"
        , col = heat.colors(2)
        , main = "K: D, A; Mon: R; Will: S; red < 0.01; [leo_genes_AML]"
        , Rowv = NA
        , Colv = NA
)


# Heatmap with continous colors

heatmap(as.matrix(pvalues_genes_AML_trans)
        , scale = "none"
        , col = heat.colors(256)
        , main = "K: D, A; Mon: R; Will: S; red < 0.01; [leo_genes_AML]"
        , Rowv = NA
        , Colv = NA
)
remove(a,d,i,g,g_pvalues_age, g_pvalues_date, g_pvalues_runs, g_pvalues_sex, norming_pvalues)



# Mono
g_PC_batch_genes_leo <- read.csv("C:/Users/pierr/Documents/USB-Backup/Uni Heidelberg/MoBi_Bsc/Semester IV/Modul_Einführung in die Informatik/Modulelement2_Anwendung bioinformatischer Systeme/Projekt/Daten/bearbeites Set von Konsti/g_PC_batch_genes_leo.csv")

g_PC_batch_genes_leo_mono <- g_PC_batch_genes_leo[11:20,]

# New updated code with only Mono
# Kruskal-Wallis for Date (_date) and Age (_age)

g_pvalues_date <- c()
g_pvalues_age <- c()

for (i in 1:5) {
  d <- (kruskal.test(g_PC_batch_genes_leo_mono[,i]~g_PC_batch_genes_leo_mono[,7]))$p.value
  a <- (kruskal.test(g_PC_batch_genes_leo_mono[,i]~g_PC_batch_genes_leo_mono[,9]))$p.value
  g_pvalues_date <- append(g_pvalues_date, d)
  g_pvalues_age <- append(g_pvalues_age, a)
}

######################################################################################################

# Monte carlo for runs (_runs)
# First installation of package surveillance is required
library(surveillance)

g_pvalues_runs <- c()
for (i in 1:5) {
  g <- (permutationTest(g_PC_batch_genes_leo_mono[,i], g_PC_batch_genes_leo_mono[,8], nPermutation = 9999))$pVal.permut
  g_pvalues_runs <- append(g_pvalues_runs, g)
}

######################################################################################################

# Willcoxon for Provider (_prov) and sex (_sex)

g_pvalues_prov <- c()
g_pvalues_sex <- c()

for (i in 1:5) {
  p <- (wilcox.test(g_PC_batch_genes_leo_mono[,i]~g_PC_batch_genes_leo_mono[,6]))$p.value
  s <- (wilcox.test(g_PC_batch_genes_leo_mono[,i]~g_PC_batch_genes_leo_mono[,11]))$p.value
  g_pvalues_prov <- append(g_pvalues_prov, p)
  g_pvalues_sex <- append(g_pvalues_sex, a)
}

######################################################################################################

# Reporting values in a data frame

pvalues_genes_mono <- data.frame(g_pvalues_prov, g_pvalues_date, g_pvalues_runs, g_pvalues_age, g_pvalues_sex, row.names = c("PC1", "PC2", "PC3", "PC4", "PC5"))
colnames(pvalues_genes_mono) <- c("Provider", "Date", "Runs", "Age", "Sex")
pvalues_genes_mono_trans <- t(pvalues_genes_mono)

#######################################################################################################

# Norming data frame for two colored heatmap

norming_pvalues <- function(x){
  if(x<=0.01){
    x = 1
  }
  else{
    x=10
  }
}

pvalues_genes_mono_norm<- apply(pvalues_genes_mono_trans
                                , c(1,2)
                                , norming_pvalues)

#####################################################################################################

# Creating heatmaps

# Heatmap with two colors
heatmap(as.matrix(pvalues_genes_mono_norm)
        , scale = "none"
        , col = heat.colors(2)
        , main = "K: P, D, A; Mon: R; Will: S; red < 0.01; [leo_genes_mono]"
        , Rowv = NA
        , Colv = NA
)


# Heatmap with continous colors

heatmap(as.matrix(pvalues_genes_mono_trans)
        , scale = "none"
        , col = heat.colors(256)
        , main = "K: P, D, A; Mon: R; Will: S; red < 0.01; [leo_genes_mono]"
        , Rowv = NA
        , Colv = NA
)
remove(s,d,p,i,g,g_pvalues_age, g_pvalues_date, g_pvalues_prov, g_pvalues_runs, g_pvalues_sex, norming_pvalues)
