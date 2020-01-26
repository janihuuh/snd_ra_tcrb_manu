
## As many autoimmune disorders have been associated with common viral infections, we hard-matched the amino acid sequences within these clusters against VDJdb,
##  currently the biggest repository of TCRs with known epitope-specificities, mostly against viral epitopes [ref Bagaev NAR 2019

## Read in the files produced with VDJtools
dr_vdjdb <- fread("results/vdjdb/dr_.dr_tcrgp_tot.annot.txt")  %>% filterVDJdb %>% preprocessMultiVDJdb
sn_vdjdb <- fread("results/vdjdb/sn_.sn_tcrgp_tot.annot.txt")  %>% filterVDJdb %>% preprocessMultiVDJdb
sp_vdjdb <- fread("results/vdjdb/sp_.sp_tcrgp_tot.annot.txt")  %>% filterVDJdb %>% preprocessMultiVDJdb
hc_vdjdb <- fread("results/vdjdb/hc_.hc_tcrgp_tot.annot.txt")  %>% filterVDJdb %>% preprocessMultiVDJdb

tot_vdjdb <- rbind(dr_vdjdb, sn_vdjdb, sp_vdjdb, hc_vdjdb) %>%
  mutate(cohort = substr(name,1,2)) %>%
  mutate(cohort = plyr::revalue(as.factor(cohort), c("DR" = "CND-RA", "SN" = "SN-RA", "SP" = "SP-RA")))



## Visualize results, by species
df <- tot_vdjdb %>% filter(target == "Single") %>% filter(count > 1) %>%
  group_by(name, antigen.species, cohort, .drop = T) %>% summarise(tot.freq = sum(freq))

df <- df %>% ungroup() %>%
  mutate(name = as.factor(name), antigen.species = as.factor(antigen.species)) %>%
  tidyr::complete(name, antigen.species, fill = list(tot.freq = 0)) %>%
  mutate(cohort =  plyr::revalue(as.factor(substr(name, 1, 2)), c("DR" = "CND-RA", "SN" = "SN-RA", "SP" = "SP-RA")) ) %>%
  filter(antigen.species %in% c("CMV", "EBV", "InfluenzaA"))

df %>% filter(antigen.species %in% c("CMV", "EBV", "InfluenzaA")) %>%
  ggplot(aes(cohort,tot.freq,fill=cohort)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) +
    facet_wrap(~antigen.species, scales = "free") +
    scale_fill_manual(values = getPalette3(4)) +
    ggsignif::geom_signif(comparisons = list(c("CND-RA", "HC"), c("CND-RA", "SN-RA"), c("CND-RA", "SP-RA")),  map_signif_level = T, y_position = c(0.055, 0.05, 0.045)) + theme_bw() + facets_nice + theme(legend.position = "none") + labs(x = "", y = "freq of CD8+ repertoire") +
  ylim(values = c(0, 0.058))
ggsave("results/vdjdb/plots/box_antigen_species_min1.pdf", width = 8, height = 3.5)



## Visualize results, by epitope
df <- tot_vdjdb %>% filter(target == "Single") %>% filter(count > 1) %>%
  group_by(name, antigen.epitope, antigen.species, cohort, .drop = T) %>% summarise(tot.freq = sum(freq))

df <- df %>% ungroup() %>% filter(antigen.species %in% c("CMV", "EBV", "InfluenzaA")) %>%
  mutate(name = as.factor(name), antigen.epitope = as.factor(antigen.epitope)) %>%
  tidyr::complete(name, antigen.epitope, fill = list(tot.freq = 0)) %>%
  mutate(cohort =  plyr::revalue(as.factor(substr(name, 1, 2)), c("DR" = "CND-RA", "SN" = "SN-RA", "SP" = "SP-RA")) )

df %>%
  ggplot(aes(cohort,tot.freq,fill=cohort)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) +
    facet_wrap(~antigen.epitope, scales = "free") +
    scale_fill_manual(values = getPalette3(4)) + theme(legend.position = "none") +
    ggsignif::geom_signif(comparisons = list(c("CND-RA", "HC"), c("CND-RA", "SN-RA"), c("CND-RA", "SP-RA")), step_increase = 0.02, map_signif_level = T) + theme_bw() + facets_nice + theme(legend.position = "none") + labs(x = "", y = "freq of CD8+ repertoire")
ggsave("results/vdjdb/plots/box_antigen_epitope_min1.pdf", width = 12, height = 12)
