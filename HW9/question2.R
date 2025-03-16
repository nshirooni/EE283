library(tiff)
library(ggplot2)
library(gridExtra) 
library(grid)        

file_paths <- list.files(path = "/dfs6/pub/nshiroon/EE283/DATA/ATACseq/hw6", pattern = "\\.tiff$", full.names = TRUE)
img_list <- lapply(file_paths, readTIFF)

ggplots <- lapply(img_list, function(img) {
  ggplot() +
    annotation_raster(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    theme_void() + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
})

combined_plot <- grid.arrange(grobs = ggplots, ncol = 3, nrow = 2)
ggsave("question2.tiff", plot = combined_plot, device = "tiff", width = 12, height = 8, dpi = 300)