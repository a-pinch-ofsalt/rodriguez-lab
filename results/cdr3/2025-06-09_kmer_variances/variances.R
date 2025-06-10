library(data.table)
library(stringr)

input_dir <- "stage_kmer_counts"
all_files <- list.files(input_dir, full.names = TRUE)
output_dir <- "kmer_variance_tables"
dir.create(output_dir, showWarnings = FALSE)

# Extract metadata: individual, chain, stage
file_info <- data.table(
  path = all_files,
  filename = basename(all_files)
)
file_info[, c("individual", "chain", "stage") := tstrsplit(filename, "_", fixed = TRUE)[1:3]]
file_info[, stage := as.integer(str_match(filename, "stage_(\\d+)_")[,2])]

# Group by chain and stage
grouped <- split(file_info, by = c("chain", "stage"), drop = TRUE)

for (group_key in names(grouped)) {
  group <- grouped[[group_key]]
  chain <- unique(group$chain)
  stage <- unique(group$stage)
  
  
  message("Processing chain: ", chain, " stage: ", stage)
  
  # Load top 10000 kmers from each file
  top_kmers_list <- lapply(group$path, function(f) {
    dt <- fread(f)
    dt <- dt[order(-count)][1:min(10000, .N)]
    individual <- str_match(basename(f), "^([A-Z])_")[,2]
    dt[, .(kmer, count, individual)]
  })
  
  combined_dt <- rbindlist(top_kmers_list)
  if (nrow(combined_dt) == 0) next
  
  # Get global top 10000 kmers by total count
  summed <- combined_dt[, .(total_count = sum(count)), by = kmer]
  top10000_kmers <- summed[order(-total_count)][1:min(10000, .N), kmer]
  
  # Get individual counts for top kmers
  counts_per_indiv <- combined_dt[kmer %in% top10000_kmers, .(count = sum(count)), by = .(kmer, individual)]
  counts_wide <- dcast(counts_per_indiv, kmer ~ individual, value.var = "count", fill = 0)
  
  # Total kmer count for scaling
  total_kmer_count <- sum(sapply(group$path, function(f) sum(fread(f)$count)))
  
  # Scaled IQR calculation
  counts_matrix <- as.matrix(counts_wide[, -1, with = FALSE])
  scaled_iqr <- apply(counts_matrix, 1, function(x) {
    (quantile(x, 0.75) - quantile(x, 0.25)) / total_kmer_count
  })
  
  # Output
  output_dt <- data.table(kmer = counts_wide$kmer, scaled_iqr = scaled_iqr)
  fwrite(output_dt, file = file.path(output_dir, paste0("chain_", chain, "_stage_", stage, "_kmer_variance.csv")))
}

v <- fread(paste0("kmer_variance_tables/chain_H_stage_1_kmer_variance.csv"))
hist(x = v$iqr, xlim = c(0, 0.00005), breaks = 120, main = "H1")
?hist
library(data.table)

# Load all variance files
files <- list.files("kmer_variance_tables", full.names = TRUE, pattern = "\\.csv$")

# Extract metadata
meta <- data.table(
  path = files,
  fname = basename(files),
  chain = sub("chain_(.)_stage_.*", "\\1", basename(files)),
  stage = as.integer(sub(".*_stage_(\\d+)_.*", "\\1", basename(files)))
)

# For each chain, plot all stages in one image
unique_chains <- unique(meta$chain)
unique_chains
for (chain in unique_chains) {
  subset_meta <- meta[meta$chain == chain, ][order(stage)]
  n_stages <- nrow(subset_meta)
  
  # Save to PNG
  png(filename = paste0("kmer_variance_tables/hist_chain_", chain, "_all_stages.png"),
      width = 1600, height = 300 * ceiling(n_stages / 3))
  
  # Setup grid layout: 3 per row
  par(mfrow = c(ceiling(n_stages / 3), 3), mar = c(4, 4, 2, 1))
  
  for (i in 1:n_stages) {
    v <- fread(subset_meta$path[i])
    stage_label <- subset_meta$stage[i]
    
    fixed_breaks <- seq(0, 0.000051, length.out = 41)  # pad a bit past 0.00005
    # 120 bins between 0 and 0.00005
    clipped_vals <- v$scaled_iqr[v$scaled_iqr >= 0 & v$scaled_iqr <= 0.00005]
    hist(
      clipped_vals,
      breaks = fixed_breaks,
      xlim = c(0, 0.00005),
      main = paste0("Stage ", stage_label),
      xlab = "Scaled IQR"
    )
    
    
  }
  dev.off()
  
}
