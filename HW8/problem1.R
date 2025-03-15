library(tidyverse)

# Read the dataset
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
    return(tibble(treat_p = NA, founder_p = NA, interaction_p = NA))
  }
  
  model <- lm(asin(sqrt(freq)) ~ treat + founder + treat:founder, data = df)
  anova_table <- anova(model)
  
  tibble(
    treat_p = anova_table$`Pr(>F)`[1],
    founder_p = anova_table$`Pr(>F)`[2],
    interaction_p = anova_table$`Pr(>F)`[3]
  )
}

mal_results <- mal_subset %>%
  group_by(chr, pos) %>%
  nest() %>%
  mutate(anova_results = map(data, safely(run_anova))) %>%
  mutate(anova_results = map(anova_results, ~ .x$result %||% tibble(treat_p = NA, founder_p = NA, interaction_p = NA))) %>%
  unnest(anova_results) %>%
  mutate(
    treat_logp = -log10(treat_p),
    founder_logp = -log10(founder_p),
    interaction_logp = -log10(interaction_p)
  ) %>%
  select(chr, pos, treat_logp, founder_logp, interaction_logp)

print(mal_results)
write.csv(mal_results,file = 'problem1_results.txt')

