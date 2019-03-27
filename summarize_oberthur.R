seeds <- 1:99
directory <- "../dataset/oberthuer/oberthur/"
seed_strings <- c(formatC(seeds, width=2, flag="0"), '100')
full_filenames <- paste0(directory, seed_strings, "_", "test_result.csv")

full_df <- c()
for (filename in full_filenames) {
  test_result <- read.csv(filename)
  full_df <- rbind(test_result, full_df)
}
full_filename <- '../dataset/oberthuer/deviance_boxplot2.pdf'
pdf(full_filename, width=6, height=6)
boxplot(full_df$deviance, ylab='Difference in deviance')
abline(h=0)
dev.off()
