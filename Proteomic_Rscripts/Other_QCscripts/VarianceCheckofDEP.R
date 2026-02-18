mat <- assay(data_filt)

protein_var <- apply(mat, 1, var, na.rm = TRUE)

summary(protein_var)
hist(protein_var, breaks = 50)


protein_mean <- rowMeans(mat, na.rm = TRUE)

plot(
  protein_mean,
  protein_var,
  pch = 16, cex = 0.4,
  xlab = "Mean intensity",
  ylab = "Variance",
  main = "Mean–variance relationship"
)

# Fit linear model
mv_lm <- lm(protein_var ~ protein_mean)

# Extract slope
slope <- coef(mv_lm)[2]
slope
