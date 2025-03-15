library(tidyverse)

mal = read_tsv("~/allhaps.malathion.200kb.txt.gz")

selected_positions <- mal %>%
  distinct(chr, pos) %>%
  slice_sample(n = 10) 

mal_subset <- mal %>%
  semi_join(selected_positions, by = c("chr", "pos"))


run_anova_model1 <- function(df) {
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

run_anova_model2 <- function(df) {
  df <- df %>%
    mutate(
      treat = as.factor(str_sub(pool, 2, 2)),  
      founder = as.factor(founder)
    )
  
  if (nlevels(df$treat) < 2 | nlevels(df$founder) < 2) {
    return(tibble(founder_p2 = NA, nested_treat_p = NA))
  }
  
  model <- lm(asin(sqrt(freq)) ~ founder + treat %in% founder, data = df)
  anova_table <- anova(model)
  
  tibble(
    founder_p2 = anova_table$`Pr(>F)`[1], 
    nested_treat_p = anova_table$`Pr(>F)`[2] 
  )
}


mal_results1 <- mal_subset %>%
  group_by(chr, pos) %>%
  nest() %>%
  mutate(anova_results = map(data, safely(run_anova_model1))) %>%
  mutate(anova_results = map(anova_results, ~ .x$result %||% tibble(treat_p = NA, founder_p = NA, interaction_p = NA))) %>%
  unnest(anova_results) %>%
  mutate(
    treat_logp = -log10(treat_p),
    founder_logp = -log10(founder_p),
    interaction_logp = -log10(interaction_p)
  ) %>%
  select(chr, pos, treat_logp, founder_logp, interaction_logp)

mal_results2 <- mal_subset %>%
  group_by(chr, pos) %>%
  nest() %>%
  mutate(anova_results = map(data, safely(run_anova_model2))) %>%
  mutate(anova_results = map(anova_results, ~ .x$result %||% tibble(founder_p2 = NA, nested_treat_p = NA))) %>%
  unnest(anova_results) %>%
  mutate(
    founder_logp2 = -log10(founder_p2),
    nested_treat_logp = -log10(nested_treat_p)
  ) %>%
  select(chr, pos, founder_logp2, nested_treat_logp)

merged_results <- full_join(mal_results1, mal_results2, by = c("chr", "pos"))

print(merged_results)
write.csv(merged_results,file = 'problem3_results.txt')
