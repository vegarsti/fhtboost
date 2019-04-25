library(devtools)
load_all()
seeds <- 1:99
seed_strings <- c(formatC(seeds, width=2, flag="0"), '100')

directory <- "../dataset/oberthuer/oberthur_both/"
full_filenames <- paste0(directory, seed_strings, "_", "test_result.csv")
full_df <- c()
for (filename in full_filenames) {
  test_result <- read.csv(filename)
  full_df <- rbind(test_result, full_df)
}
full_deviance <- full_df$deviance
# full_filename <- '../dataset/oberthuer/deviance_boxplot.pdf'
# pdf(full_filename, width=6, height=6)
# boxplot(full_deviance, ylab='Difference in deviance')
# abline(h=0)
# dev.off()


full_filenames <- paste0(directory, seed_strings, "_", "beta.csv")
full_df <- c()
genes <- rep(0, 9987)
for (filename in full_filenames) {
  test_result <- read.csv(filename)
  chosen <- test_result$j[-1]
  genes[chosen] <- genes[chosen] + 1
}

chosen_genes_df <- data.frame(cbind(non_null_parameters(genes), genes[non_null_parameters(genes)]))
names(chosen_genes_df) <- c("j", "times")
chosen_genes_df <- chosen_genes_df[order(-chosen_genes_df$times), ]
top_20_genes <- head(chosen_genes_df, n=20)
full_filename <- '../dataset/oberthuer/top20_genes.csv'
write.csv(top_20_genes, file=full_filename, row.names=FALSE)






# CLINICAL
directory <- "../dataset/oberthuer/oberthur_clinical/"
full_filenames <- paste0(directory, seed_strings, "_", "clinical_test_result.csv")
full_df <- c()
for (filename in full_filenames) {
  test_result <- read.csv(filename)
  full_df <- rbind(test_result, full_df)
}
clinical_deviance <- full_df$deviance
full_filename <- '../dataset/oberthuer/deviance_clinical_boxplot.pdf'
# pdf(full_filename, width=6, height=6)
# boxplot(full_df$deviance, ylab='Difference in deviance')
# abline(h=0)
# dev.off()



directory <- "../dataset/oberthuer/oberthur_genetic/"
full_filenames <- paste0(directory, seed_strings, "_", "genetic", "_", "test_result.csv")
full_df <- c()
for (filename in full_filenames) {
  test_result <- read.csv(filename)
  full_df <- rbind(test_result, full_df)
}
genetic_deviance <- full_df$deviance
# full_filename <- '../dataset/oberthuer/deviance_boxplot.pdf'
# pdf(full_filename, width=6, height=6)
# boxplot(full_deviance, ylab='Difference in deviance')
# abline(h=0)
# dev.off()


full_filenames <- paste0(directory, seed_strings, "_", "genetic", "_", "beta.csv")
full_df <- c()
genes <- rep(0, 9987)
for (filename in full_filenames) {
  test_result <- read.csv(filename)
  chosen <- test_result$j[-1]
  genes[chosen] <- genes[chosen] + 1
}

chosen_genes_df <- data.frame(cbind(non_null_parameters(genes), genes[non_null_parameters(genes)]))
names(chosen_genes_df) <- c("j", "times")
chosen_genes_df <- chosen_genes_df[order(-chosen_genes_df$times), ]
top_20_genes <- head(chosen_genes_df, n=20)
full_filename <- '../dataset/oberthuer/top20_genes_genetic_only.csv'
write.csv(top_20_genes, file=full_filename, row.names=FALSE)






full_filename <- '../dataset/oberthuer/deviance_all_boxplot_w_title.pdf'
pdf(full_filename, width=12, height=6)
boxplot(full_deviance, genetic_deviance, clinical_deviance, xlab='Difference in deviance', horizontal=TRUE,
        main="Difference in deviance for all three explanatory models")
axis(2, labels=c("Full", "Genetic", "Clinical"), at=1:3, las=2)
abline(v=0, lty=3)
dev.off()

full_filename <- '../dataset/oberthuer/deviance_all_boxplot_no_title.pdf'
pdf(full_filename, width=12, height=6)
boxplot(full_deviance, genetic_deviance, clinical_deviance, xlab='Difference in deviance', horizontal=TRUE)
axis(2, labels=c("Full", "Genetic", "Clinical"), at=1:3, las=2)
abline(v=0, lty=3)
dev.off()

