library(data.table)

# --------------- Core helper functions ---------------

all_repertoires <- function(func,
                            dirs = c("A","B","C","D","E","G"),
                            files = c("H_F","H_T","K_F","K_T","L_F","L_T")) {
  sample_data_tables <- list()
  for (d in dirs) for (f in files) {
    path <- file.path(d, paste0(f, ".tsv")) 
    processed <- func(fread(path))
    sample_name <- paste(d, f, sep = "_") #i.e. A_H_T
    sample_data_table <- as.data.table(as.list(processed))
    sample_data_table[, sample := sample_name] #set 'sample' for all values in the table
    sample_data_tables[[sample_name]] <- sample_data
  }
  combined_data <- rbindlist(sample_data_tables, fill = TRUE)
  setcolorder(combined_data, c("sample", setdiff(names(combined_data), "sample")))
  combined_data
}

# Pull out only the data for a given chain and stage
get_group <- function(df, chain, stage) {
  chain_and_stage_matching_pattern <- paste0("_", chain, "_", stage, "$")# i.e. H_F
  hits <- grepl(chain_and_stage_matching_pattern, df$sample)
  df[hits, !("sample"), with = FALSE]# those rows' data, but without "sample" col
}

# Calculate paired Wilcoxon p-values between two matching tables
calc_pvals <- function(df1, df2) {
  common_cols <- intersect(names(df1), names(df2))
  p_vec <- sapply(common_cols, function(col) {
    wilcox.test(df1[[col]], df2[[col]], paired = TRUE)$p.value
  })
  data.table(feature = common_cols, p_value = unname(p_vec))[order(p_value)]
}

# --------------- Task definitions ---------------

# 1. k-mer counts (lengths 3–9)

count_kmers <- function(data) table(factor(nchar(data$cdr3_aa), levels = 3:9))

# 2. average CDR3 length
mean_cdr3    <- function(data) mean(nchar(data$cdr3_aa), na.rm = TRUE)

# 3. amino-acid frequencies
aa_freq      <- function(data) table(unlist(strsplit(paste(data$cdr3_aa, collapse = ""), "")))

# Put them in a list so you can add more without touching the loop
tasks <- list(
  kmer_counts = count_kmers,
  cdr3_means  = mean_cdr3,
  aa_freqs    = aa_freq
)

# --------------- Run tasks automatically ---------------
# for (task in names(tasks)) {
#   result <- do_for_all(tasks[[task]])
#   fwrite(result, paste0(task, ".csv"))
# }

# --------------- Automated comparisons ---------------
# Load the amino-acid frequencies (or whichever table you want)
aa_table <- fread("aa_freqs.csv")

# Detect unique chains and stages from the sample names
parts <- tstrsplit(aa_table$sample, "_")
chains <- unique(parts[[2]])
stages <- unique(parts[[3]])

# Compare every chain between every pair of stages (e.g. F vs T)
for (ch in chains) for (s1 in stages) for (s2 in stages) {
  if (s1 >= s2) next
  df1 <- get_group(aa_table, ch, s1)
  df2 <- get_group(aa_table, ch, s2)
  ptab <- calc_pvals(df1, df2)
  out_file <- sprintf("pvals_%s_%s_vs_%s.csv", ch, s1, s2)
  fwrite(ptab, out_file)
}
