
# Sequencing depth is reflected in TCR repertoire diversity i.e. to the number of individual TCR sequences in the sample.
# Thus we selected a cut-off of 10 000 reads,
# removed samples below this (n=6) and resampled the of our samples to 10 000 reads, and we ended
# in (CND-RA n=6, SN-RA n=7, SP-RA n=46 and HC n=28).

age <- fread("data/clinical.txt")
over10k_resampled <- fread("results/diversity/unsampled_over10k.diversity.aa.resampled.txt") %>%
  mutate(clonality = 1 - normalizedShannonWienerIndex_mean) %>%
  mutate(cohort = substr(sample_id, 1, 2)) %>%
  mutate(cohort = gsub("HR", "SP", cohort)) %>%
  mutate(cohort = plyr::revalue(as.factor(cohort), c("DR" = "CND-RA", "SN" = "SN-RA", "SP" = "SP-RA")))  %>%
  left_join(select(age, sample_id, age))


## Select only cols we need
div_columns <- c( "cohort", "sample_id", "age","reads", "diversity", "clonality", grep("mean", colnames(over10k_resampled), value = T))
over10k_resampled <- over10k_resampled %>% select(div_columns) %>% select(-chaoE_mean)


## Visualize the diversity indices as well as their correlation with age
over10k_resampled %>%
  melt(id = c("sample_id", "cohort", "age")) %>%
  ggplot(aes(cohort,value,fill=cohort)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) +
    facet_wrap(~variable, scales = "free", ncol = 5) +
    theme_bw() +
    scale_fill_manual(values = getPalette(4)) + theme(legend.position = "none") + ggsignif::geom_signif(comparisons = list(c("CND-RA", "HC"), c("CND-RA", "SN-RA"), c("CND-RA", "SP-RA")), step_increase = 0.05) +
  facets_nice
ggsave("results/diversity/plots/box_resampled_10k_diversity_full.pdf", width = 12, height = 8)

over10k_resampled %>%
  select("sample_id", "cohort", "diversity", "clonality", "inverseSimpsonIndex_mean") %>%
  dplyr::rename("inverse Simpson index" = inverseSimpsonIndex_mean) %>%
  melt(id = c("sample_id", "cohort")) %>%
  ggplot(aes(cohort,value,fill=cohort)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) +
  facet_wrap(~variable, scales = "free", ncol = 5) +
  theme_bw() +
  scale_fill_manual(values = getPalette3(4)) + theme(legend.position = "none") + ggsignif::geom_signif(comparisons = list(c("CND-RA", "HC"), c("CND-RA", "SN-RA"), c("CND-RA", "SP-RA")), step_increase = 0.05, map_signif_level = T) +
  facets_nice + labs(x = "")
ggsave("results/diversity/plots/box_resampled_10k_diversity_manu.pdf", width = 10, height = 4)

over10k_resampled %>%
  select("sample_id", "cohort", "diversity", "clonality", "inverseSimpsonIndex_mean", "age") %>%
  dplyr::rename("inverse Simpson index" = inverseSimpsonIndex_mean) %>%
  melt(id = c("sample_id", "cohort")) %>%
  ggplot(aes(cohort,value,fill=cohort)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) +
  facet_wrap(~variable, scales = "free", ncol = 5) +
  theme_bw() +
  scale_fill_manual(values = getPalette3(4)) + theme(legend.position = "none") + ggsignif::geom_signif(comparisons = list(c("CND-RA", "HC"), c("CND-RA", "SN-RA"), c("CND-RA", "SP-RA")), step_increase = 0.07, map_signif_level = T) +
  facets_nice + labs(x = "")
ggsave("results/diversity/plots/box_resampled_10k_diversity_manu.pdf", width = 10, height = 4)

over10k_resampled %>%
  select("sample_id", "cohort", "diversity", "clonality", "inverseSimpsonIndex_mean", "age") %>%
  dplyr::rename("inverse Simpson index" = inverseSimpsonIndex_mean) %>%
  melt(id = c("sample_id", "cohort", "age")) %>%
  ggplot(aes(age,value,fill=cohort)) + geom_point(shape = 21, size = 3) +
  facet_wrap(~variable, scales = "free", ncol = 5) +
  theme_bw() + geom_smooth(method = "lm", aes(color = cohort), fill = NA) +
  scale_color_manual(values = getPalette3(4)) +
  scale_fill_manual(values = getPalette3(4)) +
  #theme(legend.position = "none") +
  #ggsignif::geom_signif(comparisons = list(c("CND-RA", "HC"), c("CND-RA", "SN-RA"), c("CND-RA", "SP-RA")), step_increase = 0.05, map_signif_level = T) +
  facets_nice
ggsave("results/diversity/plots/box_resampled_10k_diversity_manu_age.pdf", width = 15, height = 5)
