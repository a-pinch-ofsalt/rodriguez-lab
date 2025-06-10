m <- fread("CDR3_length_means.csv")
m[, c("individual", "chain", "survived", "stage") := tstrsplit(sample,split='_')]
means <- m[survived=="T",]
means
means[, sample := NULL]
means[, survived := NULL]

setcolorder(means, c("individual", "chain", "stage", "V1"))
setorder(means, stage)
setorder(means, chain)
means <- dcast(means, chain + stage ~ individual, value.var = "V1")
means[, row_mean := rowMeans(.SD), .SDcols = c("A", "B", "C", "D", "E", "G")]
means[, c("A", "B", "C", "D", "E", "G") := NULL]
means <- dcast(means, chain ~ stage, value.var = "row_mean")
fwrite(means, 'cdr3_means.csv')
means
means
png("myplot_H.png", width=800, height=800)
hist(means[chain=="H"][['V1']])
dev.off()
png("myplot_K.png", width=800, height=800)
hist(means[chain=="K"][['V1']])
dev.off()
png("myplot_L.png", width=800, height=800)
hist(means[chain=="L"][['V1']])
dev.off()

H <- 17
K <- 8
L <- 8
library(data.table)

dt <- fread("A/H_F.tsv")

counts <- list()
counts_normalized <- list()

for (i in 1:17) {
  # Get the i-th character of each CDR3 string
  pos <- substr(dt$cdr3_aa, i, i)
  
  # Make a frequency table
  pos_counts <- table(pos)
  
  # Convert to data.table and label with position
  pos_dt <- as.data.table(pos_counts)
  pos_dt[, position := i]
  pos_dt <- pos_dt[pos != "*" & pos != ""]
  print(pos_dt)
  
  
  # Normalize
  
  pos_dt_norm <- pos_dt[, .(pos, freq = N / sum(N), position)]
  
  # Store
  counts[[i]] <- pos_dt
  counts_normalized[[i]] <- pos_dt_norm
}

# Combine all
all_counts <- rbindlist(counts)
all_counts_normalized <- rbindlist(counts_normalized)

library(data.table)

aa_matrix <- dcast(all_counts_normalized, 
                   position ~ pos, 
                   value.var = "freq", 
                   fill = 0)
aa_matrix
fwrite(aa_matrix, 'amino_acid_frequencies_at_cdr3_positions_1-17.csv')


aa_matrix
all_counts
all_counts_normalized
aa_soup <- paste0(dt$cdr3_aa, collapse="")
counts <- table(strsplit(aa_soup, ""))

counts_df <- as.data.frame(as.list(counts))  # one row, named columns

rowname <- paste(dir, fname, sep = "_")

counts_df$File <- rowname

all_counts_list[[rowname]] <- counts_df



library(data.table)

# Define max CDR3 lengths per chain
cdr3_lengths <- list(H = 17, K = 8, L = 8)

# Directories and chain files to process
dirs <- c("A", "B", "C", "D", "E", "G")
chains <- c("H_T", "K_T", "L_T")

# For each dataset, process and save its normalized amino acid frequency matrix
for (dir in dirs) {
  for (chain_file in chains) {
    
    filepath <- file.path(dir, paste0(chain_file, ".tsv"))
    if (!file.exists(filepath)) next  # Skip if file is missing
    
    dt <- fread(filepath)
    
    # Get chain type ("H", "K", or "L") and its CDR3 length
    chain_type <- substr(chain_file, 1, 1)
    max_pos <- cdr3_lengths[[chain_type]]
    
    counts <- list()
    counts_normalized <- list()
    
    for (i in 1:max_pos) {
      pos <- substr(dt$cdr3_aa, i, i)
      pos_counts <- table(pos)
      pos_dt <- as.data.table(pos_counts)
      pos_dt[, position := i]
      pos_dt <- pos_dt[pos != "*" & pos != ""]
      pos_dt_norm <- pos_dt[, .(pos, freq = N / sum(N), position)]
      counts[[i]] <- pos_dt
      counts_normalized[[i]] <- pos_dt_norm
    }
    
    all_counts_normalized <- rbindlist(counts_normalized)
    
    # Pivot to wide format
    aa_matrix <- dcast(all_counts_normalized,
                       position ~ pos,
                       value.var = "freq",
                       fill = 0)
    
    # Save
    outfile <- paste0("AA_freq_", dir, "_", chain_file, ".csv")
    fwrite(aa_matrix, outfile)
  }
}















library(data.table)

cdr3_lengths <- list(H = 17, K = 8, L = 8)
dirs <- c("A", "B", "C", "D", "E", "G")
chains <- c("H_T", "K_T", "L_T")

master_dt <- list()

for (dir in dirs) {
  for (chain_file in chains) {
    
    filepath <- file.path(dir, paste0(chain_file, ".tsv"))
    if (!file.exists(filepath)) next
    
    dt <- fread(filepath)
    chain_type <- substr(chain_file, 1, 1)
    max_pos <- cdr3_lengths[[chain_type]]
    
    for (i in 1:max_pos) {
      pos <- substr(dt$cdr3_aa, i, i)
      pos_counts <- table(pos)
      pos_dt <- as.data.table(pos_counts)
      pos_dt <- pos_dt[pos != "*" & pos != ""]
      pos_dt[, `:=`(
        position = i,
        chain = chain_type,
        person = dir,
        freq = N / sum(N)
      )]
      setnames(pos_dt, "pos", "amino_acid")
      master_dt[[length(master_dt) + 1]] <- pos_dt[, .(chain, position, person, amino_acid, freq)]
    }
  }
}

# Combine and write the full table
final_dt <- rbindlist(master_dt)
yes <- dcast(final_dt, chain + position + person ~ amino_acid, value.var="freq")
fwrite(yes, "AA_freq_full_long.csv")






library(data.table)

cdr3_lengths <- list(H = 17, K = 8, L = 8)
dirs <- c("A", "B", "C", "D", "E", "G")
chains <- c("H_T", "K_T", "L_T")

master_dt <- list()

for (dir in dirs) {
  for (chain_file in chains) {
    
    filepath <- file.path(dir, paste0(chain_file, ".tsv"))
    if (!file.exists(filepath)) next
    
    dt <- fread(filepath)
    chain_type <- substr(chain_file, 1, 1)
    max_pos <- cdr3_lengths[[chain_type]]
    
    for (stage_val in unique(dt$stage)) {
      dt_stage <- dt[stage == stage_val]
      
      for (i in 1:max_pos) {
        pos <- substr(dt_stage$cdr3_aa, i, i)
        pos_counts <- table(pos)
        pos_dt <- as.data.table(pos_counts)
        pos_dt <- pos_dt[pos != "*" & pos != ""]
        pos_dt[, `:=`(
          position = i,
          chain = chain_type,
          person = dir,
          stage = stage_val,
          freq = N / sum(N)
        )]
        setnames(pos_dt, "pos", "amino_acid")
        master_dt[[length(master_dt) + 1]] <- pos_dt[, .(chain, position, person, stage, amino_acid, freq)]
      }
    }
  }
}

final_dt <- rbindlist(master_dt)
yes <- dcast(final_dt, chain + position + person + stage ~ amino_acid, value.var="freq")
head(final_dt)
fwrite(final_dt, "AA_freq_full_long_by_stage.csv")






library(data.table)
library(ggplot2)
library(gridExtra)

# Read input CSV
df <- fread("AA_freq_full_long.csv")

# Reshape to long format
long <- melt(df, id.vars = c("chain", "position", "person"),
             variable.name = "AA", value.name = "freq")
colnames(df)
# Add averaged row (across people)
avg <- long[, .(freq = mean(freq, na.rm = TRUE)), by = .(chain, position, AA)]
avg[, person := "AVG"]
long_all <- rbind(long, avg)

# Define number of stages per chain
n_stages <- list(H = 4, K = 3, L = 3)

# Main plotting function
plot_chain_grid <- function(chain_type) {
  subset <- long_all[chain == chain_type]
  
  # Split by stage: we'll assume positions map to stages like 1-4 = stage 1, 5-8 = stage 2, etc.
  # Modify this logic if stage info is separate
  pos_per_stage <- max(subset$position) / n_stages[[chain_type]]
  subset[, stage := ceiling(position / pos_per_stage)]
  
  # Generate one heatmap per (person, stage)
  top_aa_table <- subset[, .SD[which.max(freq)], by = .(person, stage, position)]

  people <- unique(subset$person)
  plot_list <- list()
  for (p in people) {
    for (s in 1:n_stages[[chain_type]]) {
      dt <- subset[person == p & stage == s]
      g <- ggplot(dt, aes(x = factor(position), y = AA, fill = freq)) +
        geom_tile(color = "white") +
        scale_fill_gradient(low = "white", high = "red") +
        labs(title = paste(p, "- Stage", s), x = "Position", y = "AA") +
        theme_minimal(base_size = 9) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 5),
              plot.title = element_text(size = 8),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
      plot_list[[paste(p, s)]] <- g
    }
  }
  
  grid <- marrangeGrob(plot_list, nrow = n_stages[[chain_type]], ncol = 7,
                       top = paste("AA Frequencies – Chain", chain_type))
  ggsave(paste0("AA_freq_grid_chain_", chain_type, ".png"), grid,
         width = 20, height = 10)
}

# Generate and save all
for (ch in c("H", "K", "L")) {
  plot_chain_grid(ch)
}
# Make sure freq is numeric
long_all[, freq := as.numeric(freq)]

# Get top AA per group
top_aa_table <- long_all[
  , .SD[which.max(freq)],
  by = .(chain, person, stage, position)
][
  , .SD[1], by = .(chain, person, stage, position)
]

# Build two tables: one for AA, one for freq
aa_table <- top_aa_table[, .(chain, position, stage, person, value = amino_acid)]
freq_table <- top_aa_table[, .(chain, position, stage, person, value = freq)]

# Cast wide by person
aa_wide <- dcast(aa_table, chain + position + stage ~ person, value.var = "value")
freq_wide <- dcast(freq_table, chain + position + stage ~ person, value.var = "value")

# Rename columns to show AA/freq
aa_cols <- setdiff(names(aa_wide), c("chain", "position", "stage"))
freq_cols <- setdiff(names(freq_wide), c("chain", "position", "stage"))

setnames(aa_wide, aa_cols, paste0(aa_cols, "_AA"))
setnames(freq_wide, freq_cols, paste0(freq_cols, "_freq"))

# Merge them on chain + position + stage
final <- merge(aa_wide, freq_wide, by = c("chain", "position", "stage"))

# Save it
fwrite(final, "top_AA_per_position_chain_stage_wide.csv")



library(data.table)
library(ggplot2)
library(gridExtra)

# Read the file
df <- fread("AA_freq_full_long_by_stage.csv")
colnames(df)
df
# Make sure these columns exist: chain, position, person, stage, A-Y...
if (!"stage" %in% colnames(df)) stop("Missing 'stage' column!")

long <- copy(df)  # or just use df directly
avg <- long[, .(freq = mean(freq, na.rm = TRUE)), by = .(chain, position, stage, amino_acid)]
avg[, person := "AVG"]
setcolorder(avg, c("chain", "position", "person", "stage", "amino_acid", "freq"))

long_all <- rbind(long, avg)


# Number of stages per chain (just for safety)
n_stages <- list(H = 4, K = 3, L = 3)

# Function to plot each chain
plot_chain_grid <- function(chain_type) {
  subset <- long_all[chain == chain_type]
  people <- unique(subset$person)
  stages <- sort(unique(subset$stage))
  subset[, freq := as.numeric(freq)]
  
  top_aa_table <- subset[, .SD[which.max(freq)], by = .(person, stage, position)]
  
  
  plots <- list()
  for (p in people) {
    for (s in stages) {
      dt <- subset[person == p & stage == s]
      dt[, rank := frank(-freq, ties.method = "min"), by = position]
      
      top_labels <- top_aa_table[person == p & stage == s]
      
      g <- ggplot(dt, aes(x = factor(position), y = amino_acid, fill = freq)) +
        geom_tile(color = "white") +
        geom_text(aes(label = rank), size = 2, color = "black") +
        geom_text(data = top_labels, aes(x = factor(position), y = max(dt$AA), label = AA),
                  inherit.aes = FALSE, vjust = -1.5, size = 3.5, fontface = "bold") +
        scale_fill_gradient(low = "white", high = "red") +
        labs(title = paste(p, "- Stage", s), x = "Position", y = "AA") +
        theme_minimal(base_size = 9) +
        theme(axis.text.x = element_text(angle = 90, size = 6),
              axis.text.y = element_text(size = 6),
              plot.title = element_text(size = 8),
              axis.title = element_blank())
      
      plots[[paste(p, s)]] <- g
    }
  }
  
  
  grid <- marrangeGrob(plots, nrow = length(stages), ncol = 7,
                       top = paste("AA Frequencies – Chain", chain_type))
  ggsave(paste0("AA_freq_grid_chain_", chain_type, ".png"), grid, width = 20, height = 10)
}

# Run for each chain
for (ch in c("H", "K", "L")) {
  plot_chain_grid(ch)
}






library(data.table)
library(ggplot2)
library(gridExtra)

# Load long-format data
df <- fread("AA_freq_full_long_by_stage.csv")

# Sanity check
stopifnot(all(c("chain", "position", "person", "stage") %in% colnames(df)))

# Average across people
avg <- df[, .(freq = mean(freq, na.rm = TRUE)), by = .(chain, position, stage, amino_acid)]
avg[, person := "AVG"]
df_all <- rbind(df, avg)

# Define number of stages per chain
n_stages <- list(H = 4, K = 3, L = 3)

# Plot function using dcast + geom_tile
plot_chain_grid <- function(chain_type) {
  subset <- df_all[chain == chain_type]
  people <- unique(subset$person)
  stages <- sort(unique(subset$stage))
  subset[, freq := as.numeric(freq)]
  
  top_aa_table <- subset[, .SD[which.max(freq)], by = .(person, stage, position)]
  
  
  plots <- list()
  for (p in people) {
    for (s in stages) {
      dt <- subset[person == p & stage == s]
      dt[, rank := frank(-freq, ties.method = "min"), by = position]
      
      top_labels <- top_aa_table[person == p & stage == s]
      dt_with_top <- rbind(dt, data.table(
        position = unique(dt$position),
        amino_acid = "TOP",  # dummy row
        freq = NA_real_,
        rank = NA_integer_,
        person = p,
        stage = s,
        chain = chain_type
      ))
      
      aa_levels <- c(sort(unique(dt$amino_acid)), "TOP")
      dt_with_top[, amino_acid := factor(amino_acid, levels = aa_levels)]
      
      
      g <- ggplot(dt_with_top, aes(x = factor(position), y = amino_acid, fill = freq)) +
        geom_tile(color = "white", na.rm = TRUE) +
        geom_text(data = dt, aes(label = rank), size = 2, color = "black", family="mono", alpha=0.25) +
        scale_fill_gradient(low = "white", high = "red", na.value = "white") +
        labs(title = paste(p, "- Stage", s), x = "Position", y = "amino_acid") +
        theme_minimal(base_size = 9) +
        theme(axis.text.x = element_text(angle = 90, size = 6),
              axis.text.y = element_text(size = 6),
              plot.title = element_text(size = 8),
              axis.title = element_blank())
      
      
      
      
      plots[[paste(p, s)]] <- g
    }
  }
  
  
  grid <- marrangeGrob(plots, nrow = n_stages[[chain_type]], ncol = 7,
                       top = paste("AA Frequencies – Chain", chain_type))
  ggsave(paste0("AA_freq_grid_chain_", chain_type, ".png"), grid, width = 20, height = 10)
}

# Run it
for (ch in c("H", "K", "L")) {
  plot_chain_grid(ch)
}






#DONT USE CHAIN_DATA
freqs <- c(fread("combined_AA_freq_H_T.csv"), fread("combined_AA_freq_K_T.csv"), fread("combined_AA_freq_L_T.csv"))
h_data <- fread("combined_AA_freq_H_T.csv")
k_data <- fread("combined_AA_freq_K_T.csv")
l_data <- fread("combined_AA_freq_L_T.csv")
head(h_data)

h_data
person_freqs <- split(h_data, h_data$person)[[1]]

for freq_set in person_freqs) {
  freqs <- 
}




ggplot(h_data, aes(x = pos, y = factor(position), fill = freq)) +
  geom_tile() +
  facet_wrap(~ person, ncol = 3) +
  scale_fill_viridis_c(option = "C") +
  labs(title = "CDR3 AA Frequencies by Position per Person (H chain)",
       x = "Amino Acid", y = "Position") +
  theme_minimal()
class(chain_data)




for (x in chain_data) {
  y <- rbindlist(x)
  
  print(y[pos=='A',])
  
  
}
library(ggplot2)
h_data <- rbindlist(chain_data[["H"]])
k_data <- rbindlist(chain_data[["K"]])
l_data <- rbindlist(chain_data[["L"]])

split(h_data, h_data$stage)

s1 <- h_data[stage=='4' & position=='17']
s1

stage_4 <- h_data[stage=='4']
stage_3 <- h_data[stage=='3']
aas <- unique(h_data$pos)



lapply(h_data[stage=='4'], h_data)
rbindlist(lapply(1:17, function(cdr3_position){
  s1 <- stage_4[position==cdr3_position]
  s2 <- stage_3[position==cdr3_position]
  aa_pvals_this_position <- lapply(aas, function(aa){
    data.frame(
      pos = cdr3_position,
      aa = aa,
      pval = wilcox.test(s1[pos==aa]$freq, s2[pos==aa]$freq, paired=TRUE)$p.value
    )
    
  })
  
  dcast(rbindlist(aa_pvals_this_position), pos ~ aa, value.var = "pval")
}))


unique(h_data$stage)

h_data <- rbindlist(chain_data[["H"]])
h_data

h_data <- fread("AA_freq_A_H_T.csv")
for (stg in sort(unique(h_data$stage))) {
  subset_stage <- h_data[stage == stg]
  
  p <- ggplot(subset_stage, aes(x = position, y = freq, color = pos, group = interaction(person, pos))) +
    geom_line(alpha = 0.3) +
    stat_summary(aes(group = pos), fun = mean, geom = "line", size = 1.2) +
    labs(title = paste("H Chain - Stage", stg),
         x = "CDR3 Position", y = "Normalized Frequency", color = "AA") +
    theme_minimal() +
    xlim(4, 17) +
    ylim(0, .3)
  
  print(p)
  ggsave(paste0("H_chain_stage_", stg, "_lineplot.png"), p, width = 10, height = 6)
}

ggplot(h_data, aes(x = pos, y = factor(position), fill = freq)) +
  geom_tile() +
  facet_wrap(~ person, ncol = 3) +
  scale_fill_viridis_c(option = "C") +
  labs(title = "CDR3 AA Frequencies by Position per Person (H chain)",
       x = "Amino Acid", y = "Position") +
  theme_minimal()
class(chain_data)




#STATISTICALLY SIGNIFICANT------------------------------------------------------------------


# First, make sure your input data is all combined
h_data <- rbindlist(chain_data[["H"]])  # has columns like: position, pos (amino acid), stage, freq

# Make sure 'stage' is a character or factor first so comparisons work
h_data[, stage := as.character(stage)]

# Get the list of unique CDR3 positions (e.g., 1 to 17)
cdr3_positions <- sort(unique(h_data$position))

# Get the list of all amino acids
aa_list <- unique(h_data$pos)

# Now run the comparison
aa_position_pvals <- rbindlist(lapply(cdr3_positions, function(cdr3_pos) {
  s1 <- h_data[stage == "3" & position == cdr3_pos]
  s2 <- h_data[stage == "4" & position == cdr3_pos]
  
  lapply(aa_list, function(aa) {
    v1 <- s1[pos == aa]$freq
    v2 <- s2[pos == aa]$freq
    
    # Only run test if there's data
    if (length(v1) < 2 || length(v2) < 2 || all(is.na(v1)) || all(is.na(v2))) return(NULL)
    
    test <- wilcox.test(v1, v2, paired = TRUE, exact = FALSE)
    data.table(
      position = cdr3_pos,
      aa = aa,
      pval = signif(test$p.value, 3),
      significant = test$p.value < 0.05
    )
  }) |> rbindlist()
}), fill = TRUE)
ss
# Optional: turn it into a wide table with one row per position, one column per AA
aa_wide <- dcast(aa_position_pvals, position ~ aa, value.var = "pval")



library(data.table)

# Assuming chain_data contains a list with H, K, L each being a list of data.tables per file
# and each one has columns: position, pos (AA), stage, freq
fread()
# Combine each chain's data into one data.table
h_data <- rbindlist(chain_data[["H"]])
k_data <- rbindlist(chain_data[["K"]])
l_data <- rbindlist(chain_data[["L"]])


k_data
all_chains <- list(H = h_data, K = k_data, L = l_data)
names(all_chains)
# Make sure 'stage' is numeric so you can do arithmetic comparisons
for (chain in names(all_chains)) {
  all_chains[[chain]][, stage := as.numeric(as.character(stage))]
}
all_chains

# Now do stage comparisons for each chain
all_pvals <- rbindlist(lapply(names(all_chains), function(chain) {
  dt <- all_chains[[chain]]
  aa_list <- unique(dt$pos)
  pos_list <- unique(dt$position)
  stage_list <- sort(unique(dt$stage))
  print(stage_list)
  
  # Compare each adjacent pair of stages (1-2, 2-3, 3-4, etc.)
  rbindlist(lapply(seq_len(length(stage_list) - 1), function(i) {
    s1 <- stage_list[i]
    s2 <- stage_list[i + 1]
    
    # For each CDR3 position
    rbindlist(lapply(pos_list, function(pos) {
      dat1 <- dt[stage == s1 & position == pos]
      dat2 <- dt[stage == s2 & position == pos]
      
      rbindlist(lapply(aa_list, function(aa) {
        v1 <- dat1[pos == aa]$freq
        v2 <- dat2[pos == aa]$freq
        
        if (length(v1) < 2 || length(v2) < 2 || all(is.na(v1)) || all(is.na(v2))) return(NULL)
        
        test <- wilcox.test(v1, v2, paired = TRUE, exact = FALSE)
        
        data.table(
          chain = chain,
          position = pos,
          aa = aa,
          stage_from = s1,
          stage_to = s2,
          pval = signif(test$p.value, 3),
          significant = test$p.value < 0.05
        )
      }))
    }))
  }))
}), fill = TRUE)

