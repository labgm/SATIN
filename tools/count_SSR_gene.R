#!/home/cwillian/mambaforge/envs/SATIN/bin/Rscript

#####!/usr/bin/Rscript

#setwd("/Your/Working/Directory")
setwd("/home/cwillian/Documents/SATIN/SATIN-main/tools/")


library(dplyr)
library(tidyr)

SSR_contagem1 <- read.delim("SSR_counting.txt", header = TRUE, stringsAsFactors = FALSE)
SSR_Groups <- read.delim("ecoli_types.txt", header = FALSE, stringsAsFactors = FALSE)


rownames(SSR_Groups) <- SSR_Groups$V1
tmp <- select(SSR_Groups, V2)
SSR_Groups <- tmp
remove(tmp)

dir.create('STAT')
#####STAT/SSR_AOV_sig.csv#####
test <- c("", "","AOV", "AOV","AOV", "AOV", "AOV", "SHAPIRO-WILK", "SHAPIRO-WILK", "KRUSKAL-WALLIS", "KRUSKAL-WALLIS", "SSR-IN-GENOMES")
cat(test, "\n", file = "STAT/SSR_AOV_sig.csv", sep = "\t")
test <- c("SSR","Gene", "Pr(>F)", "F_value", "Mean_Sq", "Sum_Sq", "Df", "W", "P-value", "chi_squared", "p_value", "SSR_sum")
cat(test, "\n", file = "STAT/SSR_AOV_sig.csv", sep = "\t", append = TRUE)

######STAT/SSR_Tukey_sig.csv####
test <- c("SSR", "Gene", "diff", "lwr", "upr", "p.adj", "interactions")
cat(test, "\n", file = "STAT/SSR_Tukey_sig.csv", sep = "\t")
remove(test)

not_significant_genes <- list()

for (gene in unique(SSR_contagem1$Gene)) {
    
  SSR_contagem <- filter(SSR_contagem1, SSR_contagem1$Gene == gene)
  SSR_names <- SSR_contagem$SSR
  tmp = setNames(data.frame(t(SSR_contagem[,-1:-2])), SSR_contagem[,2])
  SSR_contagem <- tmp
  remove(tmp)

  
  SSR_contagem_Grouped <- transform(merge(SSR_contagem,SSR_Groups,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
  SSR_names <- append(SSR_names, "V2")
  names(SSR_contagem_Grouped) <- SSR_names
  SSR_contagem_Grouped[is.na(SSR_contagem_Grouped)] <- 0
  SSR_contagem_Grouped <- arrange(SSR_contagem_Grouped, desc(SSR_contagem_Grouped$V2))


  #########KRUSKAL WALLIS############################

  kruskal_result <- lapply(SSR_contagem_Grouped[-(ncol(SSR_contagem_Grouped))], FUN = function(x) kruskal.test(x~SSR_contagem_Grouped$V2))
  
  if (exists("kruskal_result") == FALSE ) {
    print("kruskal_result error in gene =")
    print(gene)
    next}

  #########AOV AND SHAPIRO WILK############################
  #
  ###############
  try(aov_result <- lapply(SSR_contagem_Grouped[-(ncol(SSR_contagem_Grouped))], FUN = function(x) aov(x~SSR_contagem_Grouped$V2)))

  try(shapiro_result <- lapply(SSR_contagem_Grouped[-(ncol(SSR_contagem_Grouped))], FUN = function(x) shapiro.test(x)))
  
  if (exists("shapiro_result") == FALSE ) {
    print("shapiro_result error in gene =")
    print(gene)
    shapiro_result <- list()
    }


  # significante_shapiro <- list()
  # significante_index <- list()
  # for (i in 1:length(aov_result)) {if (summary.aov(aov_result[[i]])[[1]][["Pr(>F)"]][1] < 0.05) {
  #   significante_index <- append(significante_index, aov_result[i])
  #   significante_shapiro <- append(significante_shapiro, shapiro_result[i])
  # }}

  significante_shapiro <- list()
  significante_index <- list()
  significante_kruskal <- list()
  
  
  for (i in 1:length(kruskal_result)) {if (is.na(kruskal_result[i][[1]][["p.value"]])) {
    print('Missing p-value from kruskal_result item')
    kruskal_result[i][[1]][["p.value"]] <- "NF"} #create a NF (Not Found) object
    else if (kruskal_result[i][[1]][["p.value"]] < 0.001) {
    significante_index <- append(significante_index, aov_result[i])
    significante_shapiro <- append(significante_shapiro, shapiro_result[i])
    significante_kruskal <- append(significante_kruskal, kruskal_result[i])
  }}

  
  if (length(significante_index) == 0) {
     print("Not significant SSR found! \n gene =")
     print(gene)
     not_significant_genes <- append(not_significant_genes, gene)
     next}


  #########POST HOC############################

  Tukey_result <- list()
  for (i in significante_index) {Tukey_result <- append(Tukey_result, TukeyHSD(i))}
  Tukey_names <- list()
  for (i in names(significante_index)) {Tukey_names <- append(Tukey_names, i)}

  Tukey_SSR <- NULL
  for (i in 1:length(Tukey_names)) {
    SSR = as.character(Tukey_names[i])
    a = as.data.frame(Tukey_result[i], row.names = NULL)
    d = gene
    b = rownames(Tukey_result$`SSR_contagem_Grouped$V2`)
    Tukey_SSR = rbind(Tukey_SSR, data.frame(SSR, d, a, b))}

  colnames(Tukey_SSR) <- c("SSR", "Gene", "diff", "lwr", "upr", "p.adj", "interactions")
  remove(a, d)

  #######################################################
  ###########ORGANIZAR OUTPUTS PARA SALVAR#################################
  SSR_Variacao_sig<- data.frame()

  significante_SSR_Grouped <- select(SSR_contagem_Grouped, names(significante_kruskal))

  tmp <- NULL
  for (i in 1:length(significante_index)) {
    GENE = gene
    SSR = names(significante_index[i])
    b = format(summary.aov(significante_index[[i]])[[1]][["Pr(>F)"]][1], scientific = TRUE, digits = 20) # Valor de P-value para aov() #"Pr(>F)"
    c = format(summary.aov(significante_index[[i]])[[1]][["F value"]][1], scientific = FALSE, digits = 4)  # Valor de F-value para aov()
    d = format(summary.aov(significante_index[[i]])[[1]][["Mean Sq"]][1], scientific = FALSE, digits = 4)  # Valor de "Mean Sq" para aov()
    e = format(summary.aov(significante_index[[i]])[[1]][["Sum Sq"]][1], scientific = FALSE, digits = 4)  # Valor de "Sum Sq" para aov()
    f = summary.aov(significante_index[[i]])[[1]][["Df"]][1]     # Valor de "Df" para aov()
    if (is.null(significante_shapiro[[i]])) {g = h = "NF"} else { #create a NF (Not Found) object
    g = format(significante_shapiro[[i]]$statistic[[1]], scientific = FALSE, digits = 4) # Valor de W para shapiro.test() #"W"
    h = format(significante_shapiro[[i]]$p.value, scientific = TRUE, digits = 20)  # Valor de P-value para shapiro.test()
    }
    j = summary.factor(significante_kruskal[i][[1]][[1]]) # Valor de Kruskal-Wallis chi-squared para kruskal.test()
    k = format(significante_kruskal[i][[1]][["p.value"]], scientific = TRUE, digits = 20) # Valor de P-value para kruskal.test()
    l = as.numeric(sum(significante_SSR_Grouped[,i]))
    tmp = rbind(tmp, data.frame(SSR, gene, b, c, d, e, f, g, h,  j, k, l))}
  
  #   tmp = rbind(tmp, data.frame(SSR, gene, b, c, d, e, f, g, h,  j, k))}
  # SSR_Variacao_sig <- merge(SSR_Variacao_sig, tmp, by = "SSR")

  SSR_Variacao_sig <- as.data.frame(tmp)
  names(SSR_Variacao_sig) <- c("SSR","Gene", "Pr(>F)", "F value", "Mean Sq", "Sum Sq", "Df", "W", "P-value", "chi_squared", "p_value")
  rownames(SSR_Variacao_sig) <- {}
  #SSR_Variacao_sig <- arrange(SSR_Variacao_sig, desc(SSR))
  remove(tmp)
  remove(significante_index)
  remove(significante_shapiro)
  remove(significante_kruskal)
  remove(Tukey_names)
  remove(Tukey_result)
  remove(SSR_contagem_Grouped)
  remove(SSR_contagem)
  remove(SSR_names)
  remove(aov_result)
  remove(shapiro_result)
  remove(b, c, d, e, f, g, gene, GENE, h, i, j, k, SSR)
  gc()

  
  write.table(SSR_Variacao_sig, "STAT/SSR_AOV_sig.csv", row.names=FALSE, col.names = FALSE, append = TRUE, sep = "\t")
  write.table(Tukey_SSR, "STAT/SSR_Tukey_sig.csv", row.names=FALSE, col.names = FALSE, append = TRUE, sep = "\t")
  #remove(Tukey_SSR, SSR_Variacao_sig)
}

###SALVAR RESULTADOS

