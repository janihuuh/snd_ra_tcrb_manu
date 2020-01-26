
extractFileName <- function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub(".*\\/", "", str1)
}

vdjToGliph <- function(vdj_df){

  ## Write gliph files to the vdj files
  # @ param
  # input: df from vdj

  df <- vdj_df %>% select(cdr3aa, v, j, name, freq) %>%
    dplyr::rename("CDR3b"   = "cdr3aa",
                  "TRBV"    = "v",
                  "TRBJ"    = "j",
                  "Patient" =  name,
                  "Counts"  = freq)

  return(df)

}








fisher_meta <- function(bin_df, class_a, class_b, class_vector){

  require(pbapply)

  row_fisher <- function(row, alt = 'greater', cnf = 0.95) {
    # Set the function to count Fisher's exact for each row

    f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)

    return(c(row,
             p_val = f$p.value))
  }


  # Count total number of samples in each phenotype
  class_a = which(class_vector == class_a)
  class_b = which(class_vector == class_b)

  class_a_df      <- rowSums(bin_df[,..class_a])
  class_b_df      <- rowSums(bin_df[,..class_b])

  no_class_a   <- length(class_a)
  no_class_b   <- length(class_b)

  no_samples   <- ncol(bin_df)

  tot_df        <- cbind("Phenotype_positives" = class_a_df, "Ctrl_positives" = class_b_df,
                         "Phenotype_negatives" = no_class_a - class_a_df, "Ctrl_negatives" = no_class_b - class_b_df)
  tot_df <- data.frame(tot_df)

  # Run the code to count Fisher's exact test
  p  <- t(pbapply(tot_df, 1, row_fisher))
  p  <- as.data.frame(p)
  p$cdr3aa <- rownames(bin_df)

  p <- p[order(p$p_val), ]


  ## Calculate p-adj
  p <- merge(p, data.frame(p_val = unique(p$p_val), p_val_adj = p.adjust(unique(p$p_val), method = "BH")), by = "p_val")

  return(p)

}
