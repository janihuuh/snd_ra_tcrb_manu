
## Recent work by us [35] and others [36,37] have highlighted that structurally similar TCRs can bind similar antigens.
## Thus, to find T-cell specificities against previously unknown antigens, we grouped TCRs with one amino acid mismatch into potentially epitope-specific clusters.



##############################

## Prepare for gliph

## Including naive
min_count <- 0

dr_gliph <- list.files("data/tcr_filtered_10k_downsampled/", full.names = T) %>% grep(pattern = "DR", value = T) %>%
 lapply(FUN = function(x){fread(x) %>% filter(count > min_count) %>% mutate(name = extractFileName(x)) %>% vdjToGliph()}) %>% rbindlist()

sn_gliph <- list.files("data/tcr_filtered_10k_downsampled/", full.names = T) %>% grep(pattern = "SN", value = T) %>%
  lapply(FUN = function(x){fread(x) %>% filter(count > min_count) %>% mutate(name = extractFileName(x)) %>% vdjToGliph()}) %>% rbindlist()

sp_gliph <- list.files("data/tcr_filtered_10k_downsampled/", full.names = T) %>% grep(pattern = "SP", value = T) %>%
  lapply(FUN = function(x){fread(x) %>% filter(count > min_count) %>% mutate(name = extractFileName(x)) %>% vdjToGliph()}) %>% rbindlist()

hc_gliph <- list.files("data/tcr_filtered_10k_downsampled/", full.names = T) %>% grep(pattern = "HC", value = T) %>%
  lapply(FUN = function(x){fread(x) %>% filter(count > min_count) %>% mutate(name = extractFileName(x)) %>% vdjToGliph()}) %>% rbindlist()

dir.create("results/gliph/input_files/", showWarnings = F)
fwrite(dr_gliph, "results/gliph/input_files/dr_gliph_tot.txt", sep = "\t", quote = F, row.names = F)
fwrite(sn_gliph, "results/gliph/input_files/sn_gliph_tot.txt", sep = "\t", quote = F, row.names = F)
fwrite(sp_gliph, "results/gliph/input_files/sp_gliph_tot.txt", sep = "\t", quote = F, row.names = F)
fwrite(hc_gliph, "results/gliph/input_files/hc_gliph_tot.txt", sep = "\t", quote = F, row.names = F)


## Non-naive
min_count <- 1

dr_gliph <- list.files("data/tcr_filtered_10k_downsampled/", full.names = T) %>% grep(pattern = "DR", value = T) %>%
  lapply(FUN = function(x){fread(x) %>% filter(count > min_count) %>% mutate(name = extractFileName(x)) %>% vdjToGliph()}) %>% rbindlist()

sn_gliph <- list.files("data/tcr_filtered_10k_downsampled/", full.names = T) %>% grep(pattern = "SN", value = T) %>%
  lapply(FUN = function(x){fread(x) %>% filter(count > min_count) %>% mutate(name = extractFileName(x)) %>% vdjToGliph()}) %>% rbindlist()

sp_gliph <- list.files("data/tcr_filtered_10k_downsampled/", full.names = T) %>% grep(pattern = "SP", value = T) %>%
  lapply(FUN = function(x){fread(x) %>% filter(count > min_count) %>% mutate(name = extractFileName(x)) %>% vdjToGliph()}) %>% rbindlist()

hc_gliph <- list.files("data/tcr_filtered_10k_downsampled/", full.names = T) %>% grep(pattern = "HC", value = T) %>%
  lapply(FUN = function(x){fread(x) %>% filter(count > min_count) %>% mutate(name = extractFileName(x)) %>% vdjToGliph()}) %>% rbindlist()

fwrite(dr_gliph, "results/gliph/input_files/dr_gliph_min1.txt", sep = "\t", quote = F, row.names = F)
fwrite(sn_gliph, "results/gliph/input_files/sn_gliph_min1.txt", sep = "\t", quote = F, row.names = F)
fwrite(sp_gliph, "results/gliph/input_files/sp_gliph_min1.txt", sep = "\t", quote = F, row.names = F)
fwrite(hc_gliph, "results/gliph/input_files/hc_gliph_min1.txt", sep = "\t", quote = F, row.names = F)

tot_gliph <- rbind(dr_gliph, sn_gliph, sp_gliph, hc_gliph)
fwrite(tot_gliph, "results/gliph/input_files/tot_gliph_min1.txt", sep = "\t", quote = F, row.names = F)



###########################################################################



## Read in the gliph-produced clustering (with Hamming distance=1 on non-naive clones)
tot_el    <- read.delim("results/gliph/input_files/tot_gliph_min1-clone-network.txt", header = F) %>% filter(V3 != "singleton")
tot_am    <- read.delim("results/gliph/input_files/tot_gliph_min1-convergence-groups.txt", header = F) %>% filter(V1 > 2) %>% arrange(desc(V1))
tot_gliph <- read.delim("results/gliph/input_files/tot_gliph_min1.txt")

tot_interaction_map <- data.frame(CDR3b       = do.call("c", strsplit(as.character(tot_am$V3), split = " ")),
                                  cluster     = rep(paste0("cluster", 1:nrow(tot_am)), tot_am$V1),
                                  clustersize = rep(tot_am$V1, tot_am$V1))
tot_gliph <- merge(tot_interaction_map, tot_gliph, by = "CDR3b")

tot_gliph %>% group_by(Patient) %>% summarise(tot.freq = sum(Counts)) %>% mutate(cohort = substr(Patient, 1, 2)) %>%
  mutate(cohort = plyr::revalue(as.factor(cohort), c("DR" = "CND-RA","SN" = "SN-RA", "SP" = "SP-RA"))) %>%

  ggplot(aes(cohort,tot.freq,fill=cohort)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) +
    theme_bw() +
    scale_fill_manual(values = getPalette3(4)) + theme(legend.position = "none") + ggsignif::geom_signif(comparisons = list(c("CND-RA", "HC"), c("CND-RA", "SN-RA"), c("CND-RA", "SP-RA")), step_increase = 0.1, map_signif_level = T) +
    facets_nice + labs(y = "freq of clustered \nCD8+ cells")
ggsave("results/gliph/input_files/box_clustered_min1.pdf", width = 4, height = 4)

## Show one of these clusters as a logo plot
tot_interaction_map %>% filter(cluster %in% cdn_ra_clusters$cluster) %>%
  # mutate(CDR3b = substr(CDR3b, 4, nchar(as.character(CDR3b)) - 3)) %>%
  filter(cluster == "cluster13") %>% select(CDR3b) %>%
  write.table("results/cluster13.txt", quote = F, row.names = F)
