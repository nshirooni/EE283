library(tidyverse)

mal = read_tsv("~/allhaps.malathion.200kb.txt.gz")

selected_positions <- mal %>%
  distinct(chr, pos) %>%
  slice_sample(n = 10)  

mal_subset <- mal %>%
  semi_join(selected_positions, by = c("chr", "pos"))

run_anova <- function(df) {
  df <- df %>%
    mutate(
      treat = as.factor(str_sub(pool, 2, 2)),  
      founder = as.factor(founder)  
    )
  
  if (nlevels(df$treat) < 2 | nlevels(df$founder) < 2) {
    return(tibble(founder_p = NA, nested_treat_p = NA))
  }
  
  model <- lm(asin(sqrt(freq)) ~ founder + treat %in% founder, data = df)
  anova_table <- anova(model)
  
  tibble(
    founder_p = anova_table$`Pr(>F)`[1],  
    nested_treat_p = anova_table$`Pr(>F)`[2]  
  )
}

mal_results <- mal_subset %>%
  group_by(chr, pos) %>%
  nest() %>%
  mutate(anova_results = map(data, safely(run_anova))) %>%
  mutate(anova_results = map(anova_results, ~ .x$result %||% tibble(founder_p = NA, nested_treat_p = NA))) %>%
  unnest(anova_results) %>%
  mutate(
    founder_logp = -log10(founder_p),
    nested_treat_logp = -log10(nested_treat_p)
  ) %>%
  select(chr, pos, founder_logp, nested_treat_logp)

print(mal_results)
write.csv(mal_results,file = 'problem2_results.txt')
