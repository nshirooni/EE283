library(ggplot2)
library(gridExtra)
library(grid)
library(nycflights13)
library(dplyr)

P1 <- ggplot(flights, aes(x = distance, y = arr_delay)) +
  geom_point(size = 0.3, alpha = 0.3, color = "blue") +
  geom_smooth(color = "red", se = FALSE) +
  labs(x = "distance", y = "arr_delay") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6),
        plot.margin = margin(2, 2, 2, 2))

temp_flights <- flights %>%
  group_by(carrier) %>%
  summarize(m_arr_delay = mean(arr_delay, na.rm = TRUE))

P2 <- ggplot(temp_flights, aes(x = carrier, y = m_arr_delay, fill = carrier)) +
  geom_bar(stat = "identity") +
  labs(x = "carrier", y = "m_arr_delay") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 6),
        legend.position = "none",
        plot.margin = margin(2, 2, 2, 2))

P3 <- ggplot(flights, aes(x = carrier, y = arr_delay, fill = carrier)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3) +
  labs(x = "carrier", y = "arr_delay") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 6),
        legend.position = "none",
        plot.margin = margin(2, 2, 2, 2))

P4 <- ggplot(flights, aes(x = arr_delay)) +
  geom_histogram(binwidth = 10, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(x = "arr_delay", y = "count") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        plot.margin = margin(2, 2, 2, 2))

lay <- rbind(c(1, 1, 2),
             c(1, 1, 3),
             c(1, 1, 4))

tiff("figure1.tiff", width = 7, height = 6, units = "in", res = 600)
grid.arrange(P1, P2, P3, P4, layout_matrix = lay)
dev.off()

