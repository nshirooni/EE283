library(dplyr)
library(qqman)

layout(matrix(c(1,3,2,3), nrow = 2, byrow = TRUE), widths = c(3, 2))

manhattan(merged_results_qqman_clean, 
          main = "Model 1: Treat + Founder + Interaction")

manhattan(merged_results_qqman_clean_nest, 
          main = "Model 2: Founder + Treat (Nested)")


scatter_data <- data.frame(
  model1 = merged_results_qqman_clean$P,
  model2 = merged_results_qqman_clean_nest$P
)

plot(scatter_data$model1, scatter_data$model2,
     xlab = "-log10(p) Model 1",
     ylab = "-log10(p) Model 2",
     main = "Scatter plot of -log10(p)'s",
     pch = 20, col = "darkgreen")
abline(0, 1, col = "red", lty = 2)
