

## Public clonotypes are TCR-sequences that are shared between different individuals, while private sequences are unique for each individual.
## Such clonotypes could be public either because these clonotypes are more readily generated during thymic maturation or because these clonotypes are important for the disease pathogenesis.
## To find the clonotypes that are more important to disease pathogenesis,
## we subgroups by calculating the amount of public clonotypes that are enriched to a given subgroup in comparison to other subgroup.
## Low amount of such enriched public clonotypes would mean high similarity between the tested groups, and high amount of such clonotypes would mean dissimilarity between the tested groups


## Read in the VDJtools produced files
readVDJpublic <- function(df){

  df_info <- df[,c(1:14)]
  df      <- df[,-c(1:14)]
  df[df > 0] <- 1
  rownames(df) <- df_info$cdr3aa
  return(df)

}

dr_hc_public <- fread("results/public/dr_hc.join.aa.table.txt") %>% readVDJpublic
sn_hc_public <- fread("results/public/sn_hc.join.aa.table.txt") %>% readVDJpublic
sn_hc_public <- fread("results/public/sn_hc.join.aa.table.txt") %>% readVDJpublic
sp_hc_public <- fread("results/public/sp_hc.join.aa.table.txt") %>% readVDJpublic

dr_sn_public <- fread("results/public/dr_sn.join.aa.table.txt") %>% readVDJpublic
dr_sp_public <- fread("results/public/dr_sp.join.aa.table.txt") %>% readVDJpublic
sp_sn_public <- fread("results/public/sp_sn.join.aa.summary.txt") %>% readVDJpublic

sn_tot_sp_hc_public <- fread("results/public/sn_dr_sp.join.aa.table.txt") %>% readVDJpublic
sn_tot_hc_public    <- fread("results/public/sn_dr_hc.join.aa.table.txt") %>% readVDJpublic



## Calculate enrichment
# i) dr vs hc
dr_hc_public_fisher <- fisher_meta(dr_hc_public, "DR", "HC",  substr(colnames(dr_hc_public), 1, 2))
write.table(dr_hc_public_fisher, "results/public/dr_hc_public_fisher.txt", sep = "\t", quote = F)

# ii) sn vs hc
sn_hc_public_fisher <- fisher_meta(sn_hc_public, "SN", "HC",  substr(colnames(sn_hc_public), 1, 2))
write.table(sn_hc_public_fisher, "results/public/sn_hc_public_fisher.txt", sep = "\t", quote = F)

# iii) sp vs hc
sp_hc_public_fisher <- fisher_meta(sp_hc_public, "SP", "HC",  substr(colnames(sp_hc_public), 1, 2))
write.table(sp_hc_public_fisher, "results/public/sp_hc_public_fisher.txt", sep = "\t", quote = F)

# iv) dr vs hc
hc_dr_public_fisher <- fisher_meta(dr_hc_public, "HC", "DR",  substr(colnames(dr_hc_public), 1, 2))
write.table(hc_dr_public_fisher, "results/public/hc_dr_public_fisher.txt", sep = "\t", quote = F)

# v) sn vs hc
hc_sn_public_fisher <- fisher_meta(sn_hc_public, "HC", "SN",  substr(colnames(sn_hc_public), 1, 2))
write.table(hc_sn_public_fisher, "results/public/hc_sn_public_fisher.txt", sep = "\t", quote = F)

# vi) sp vs hc
hc_sp_public_fisher <- fisher_meta(sp_hc_public, "HC", "SP",  substr(colnames(sp_hc_public), 1, 2))
write.table(hc_sp_public_fisher, "results/public/hc_sp_public_fisher.txt", sep = "\t", quote = F)

# vii) dr vs sn
dr_sn_public_fisher <- fisher_meta(dr_sn_public, "DR", "SN",  substr(colnames(dr_sn_public), 1, 2))
write.table(dr_sn_public_fisher, "results/public/dr_sn_public_fisher.txt", sep = "\t", quote = F)

# viii) dr vs sp
dr_sp_public_fisher <- fisher_meta(dr_sp_public, "DR", "SP",  substr(colnames(dr_sp_public), 1, 2))
write.table(dr_sp_public_fisher, "results/public/dr_sp_public_fisher.txt", sep = "\t", quote = F)

# ix) sp vs sn
sp_sn_public_fisher <- fisher_meta(sp_sn_public, "SP", "SN",  substr(colnames(sp_sn_public), 1, 2))
write.table(sp_sn_public_fisher, "results/public/sp_sn_public_fisher.txt", sep = "\t", quote = F)

# x) sn vs dr
sn_dr_public_fisher <- fisher_meta(dr_sn_public, "SN", "DR",  substr(colnames(dr_sn_public), 1, 2))
write.table(sn_dr_public_fisher, "results/public/sn_dr_public_fisher.txt", sep = "\t", quote = F)

# xi) sp vs dr
sp_dr_public_fisher <- fisher_meta(dr_sp_public, "SP", "DR",  substr(colnames(dr_sp_public), 1, 2))
write.table(sp_dr_public_fisher, "results/public/sp_dr_public_fisher.txt", sep = "\t", quote = F)

# xii) sn vs sp
sn_sp_public_fisher <- fisher_meta(sp_sn_public, "SN", "SP",  substr(colnames(sp_sn_public), 1, 2))
write.table(sn_sp_public_fisher, "results/public/sn_sp_public_fisher.txt", sep = "\t", quote = F)

# xiii) sn + dr vs sp
sn_tot_sp_hc_public_fisher <- fisher_meta(sn_tot_sp_hc_public, "SN", "SP", class_vector = ifelse(substr(colnames(sn_tot_sp_hc_public), 1, 2) != "SP", "SN", "SP"))
write.table(sn_tot_sp_hc_public_fisher, "results/public/sn_tot_sp_hc_public_fisher.txt", sep = "\t", quote = F)

# xiv) sn + dr vs hc
sn_tot_hc_public_fisher <- fisher_meta(sn_tot_hc_public, "SN", "HC", class_vector = ifelse(substr(colnames(sn_tot_hc_public), 1, 2) != "HC", "SN", "HC"))
write.table(sn_tot_hc_public_fisher, "results/public/sn_tot_hc_public_fisher.txt", sep = "\t", quote = F)

# xv) sn + dr vs sp
sp_sn_tot_public_fisher <- fisher_meta(sn_tot_sp_hc_public, "SP", "SN", class_vector = ifelse(substr(colnames(sn_tot_sp_hc_public), 1, 2) != "SP", "SN", "SP"))
write.table(sp_sn_tot_public_fisher, "results/public/sp_sn_tot_public_fisher.txt", sep = "\t", quote = F)

# xvi) dr + sp vs sn
dr_sp_tot_public_fisher <- fisher_meta(sn_tot_sp_hc_public, "DRSP", "SN", class_vector = ifelse(substr(colnames(sn_tot_sp_hc_public), 1, 2) != "SN", "DRSP", "SN"))
write.table(dr_sp_tot_public_fisher, "results/public/dr_sp_tot_public_fisher.txt", sep = "\t", quote = F)



## Visualise the amounts of enriched clonotypes
total_public <- rbind(dr_hc_public_fisher        %>% filter(p_val < 0.05) %>% mutate(type = "CND-RA \nvs \nHC", class = "single"),
                      sn_hc_public_fisher        %>% filter(p_val < 0.05) %>% mutate(type = "SN-RA \nvs \nHC", class = "single"),
                      sp_hc_public_fisher        %>% filter(p_val < 0.05) %>% mutate(type = "SP-RA \nvs \nHC", class = "single"),

                      hc_dr_public_fisher        %>% filter(p_val < 0.05) %>% mutate(type = "HC \nvs \nCND-RA", class = "single"),
                      hc_sn_public_fisher        %>% filter(p_val < 0.05) %>% mutate(type = "HC \nvs \nSN-RA", class = "single"),
                      hc_sp_public_fisher        %>% filter(p_val < 0.05) %>% mutate(type = "HC \nvs \nSP-RA", class = "single"),


                      dr_sn_public_fisher        %>% filter(p_val < 0.05) %>% mutate(type = "CND-RA \nvs \nSN-RA", class = "single"),
                      dr_sp_public_fisher        %>% filter(p_val < 0.05) %>% mutate(type = "CND-RA \nvs \nSP-RA", class = "single"),
                      sp_sn_public_fisher        %>% filter(p_val < 0.05) %>% mutate(type = "SP-RA \nvs \nSN-RA", class = "single"),

                      sn_dr_public_fisher        %>% filter(p_val < 0.05) %>% mutate(type = "SN-RA \nvs \nCND-RA", class = "single"),
                      sp_dr_public_fisher        %>% filter(p_val < 0.05) %>% mutate(type = "SP-RA \nvs \nCND-RA", class = "single"),
                      sn_sp_public_fisher        %>% filter(p_val < 0.05) %>% mutate(type = "SN-RA \nvs \nSP-RA", class = "single"),



                      sn_tot_sp_hc_public_fisher %>% filter(p_val < 0.05) %>% mutate(type = "CND-RA + SN-RA \nvs \nSP-RA", class = "multi"),
                      dr_sp_tot_public_fisher    %>% filter(p_val < 0.05) %>% mutate(type = "CND-RA + SP-RA \nvs \nSN-RA", class = "multi"),
                      # sp_sn_tot_public_fisher %>% filter(p_val < 0.05) %>% mutate(type = "SP vs SN tot"),
                      sn_tot_hc_public_fisher    %>% filter(p_val < 0.05) %>% mutate(type = "CND-RA + SN-RA \nvs \nHC", class = "multi"),
                      dr_sp_tot_public_fisher    %>% filter(p_val < 0.05)  %>% mutate(type = "SP-RA + SN-RA \nvs \nCND+RA", class = "multi"))

write.table(total_public, "results/public/total_public.txt", sep = "\t", quote = F)



public_df <- total_public %>% mutate(class2 = factor(as.character(class), levels = c("single", "multi"))) %>%
  mutate(type2 = as.factor(as.character(type))) %>%
  group_by(type2, class2) %>% summarise(n = n())

rbind(as.data.frame(public_df), data.frame(type2 = "CND-RA + SP-RA \nvs \nSN-RA", class2 = "multi", n = 0)) %>%
  ggplot(aes(reorder(type2, n), n, fill = type2, label = n)) + geom_bar(stat = "identity") + scale_fill_brewer(palette = "Pastel1") + theme(legend.position = "none") +
    labs(x = "", y = "# public clonotypes") + facet_wrap(~class2, nrow = 2, scales = "free_x") + geom_text() + theme_bw() + facets_nice
ggsave("results/public/amounts.pdf", width = 4, height = 5)


rbind(as.data.frame(public_df), data.frame(type2 = "CND-RA + SP-RA \nvs \nSN-RA", class2 = "multi", n = 0)) %>%
  filter(class2 == "single") %>%
  ggplot(aes(reorder(type2, n), n, fill = type2, label = n)) + geom_bar(stat = "identity") + scale_fill_brewer(palette = "Pastel1") + theme(legend.position = "none") +
  labs(x = "", y = "# statistically significant \npublic clonotypes") +
  # facet_wrap(~class2, nrow = 2, scales = "free_x") +
  geom_text() + theme_bw() + facets_nice + theme(legend.position = "none")
ggsave("results/public/bar_single.png", width = 4, height = 3)

rbind(as.data.frame(public_df), data.frame(type2 = "CND-RA + SP-RA \nvs \nSN-RA", class2 = "multi", n = 0)) %>%
  filter(class2 == "multi") %>%
  ggplot(aes(reorder(type2, n), n, fill = type2, label = n)) + geom_bar(stat = "identity") + scale_fill_brewer(palette = "Set3") + theme(legend.position = "none") +
  labs(x = "", y = "# statistically significant \npublic clonotypes") +
  # facet_wrap(~class2, nrow = 2, scales = "free_x") +
  geom_text() + theme_bw() + facets_nice + theme(legend.position = "none")
ggsave("results/public/bar_multi.png", width = 4, height = 3)



## Hand coded matrix as heatmap
mtx <-  c(NA, 0, 59, 49,
          3, NA, 56, 248,
          1, 3, NA, 15,
          1, 0, 25, NA) %>% matrix(ncol = 4, byrow = T)

colnames(mtx) <- c("CND-RA", "SN-RA", "SP-RA", "HC")
rownames(mtx) <- c("CND-RA", "SN-RA", "SP-RA", "HC")

pdf("results/public/heatmap.pdf", width = 3.5, height = 3)
pheatmap::pheatmap(mtx, cluster_cols = F, cluster_rows = F, display_numbers = T, annotation_legend = F) %>% print
dev.off()

mtx <- c(NA, 0, 59,
         0, NA, 56,
         59, 56, NA) %>% matrix(ncol = 3, byrow = T)

colnames(mtx) <- c("CND-RA", "SN-RA", "SP-RA")
rownames(mtx) <- c("CND-RA", "SN-RA", "SP-RA")

pdf("results/public/heatmap2.pdf", width = 3.5, height = 3)
pheatmap::pheatmap(mtx[1:3,1:3], cluster_cols = F, cluster_rows = F, display_numbers = T, annotation_legend = F) %>% print
dev.off()
