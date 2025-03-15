library(dplyr)
library(qqman)

merged_results_qqman_clean <- merged_results %>%
  ungroup() %>%                                   
  mutate(CHR = factor(chr, levels = unique(chr))) %>%  
  mutate(CHR = as.numeric(CHR)) %>%              
  select(CHR, pos, interaction_logp) %>%        
  rename(BP = pos, P = interaction_logp) %>%      
  mutate(P = P, SNP = merged_results$pos)            

merged_results_qqman_clean_nest <- merged_results %>%
  ungroup() %>%                                   
  mutate(CHR = factor(chr, levels = unique(chr))) %>%  
  mutate(CHR = as.numeric(CHR)) %>%             
  select(CHR, pos, nested_treat_logp) %>%        
  rename(BP = pos, P = nested_treat_logp) %>%       
  mutate(P = P, SNP = merged_results$pos)            

par(mfrow = c(2, 1))  

manhattan(merged_results_qqman_clean, 
          main = "Model 1: Treat + Founder + Interaction")

manhattan(merged_results_qqman_clean_nest, 
          main = "Model 2: Founder + Treat (Nested)")



