# run_analysis.R
library(data.table)
library(ggplot2)
library(tidytext)
library(patchwork)
library(viridisLite) # For colors

# --- Settings ---
cdr3_lengths <- list(H = 17, K = 8, L = 8)
dirs <- c("A", "B", "C", "D", "E", "G")
chains <- c("H_T", "K_T", "L_T")

# --- 1. Load Data & Calculate Frequencies ---
chain_data <- list(H = list(), K = list(), L = list())

for (chain_file in chains) {
  chain_type <- substr(chain_file, 1, 1)
  max_pos <- cdr3_lengths[[chain_type]]
  
  for (dir in dirs) {
    filepath <- file.path(dir, paste0(chain_file, ".tsv"))
    if (!file.exists(filepath)) next
    dt <- fread(filepath)
    
    for (stg in unique(dt$stage)) {
      dt_stage <- dt[stage == stg]
      for (i in 1:max_pos) {
        positions <- substr(dt_stage$cdr3_aa, i, i)
        aa_dt <- as.data.table(table(positions))
        setnames(aa_dt, c("aa", "N"))
        aa_dt <- aa_dt[aa != "*" & aa != ""]
        if (sum(aa_dt$N) > 0) {
            aa_dt[, freq := N / sum(N)]
        } else {
            aa_dt[, freq := 0]
        }
        aa_dt[, `:=`(position = i, person = dir, stage = stg, chain = chain_type)]
        chain_data[[chain_type]][[length(chain_data[[chain_type]]) + 1]] <- aa_dt[, .(aa, freq, position, person, stage, chain)]
      }
    }
  }
}

# Save combined frequency files
combined_aa_freq_files <- list()
for (chain_type_key in names(chain_data)) {
  if (length(chain_data[[chain_type_key]]) == 0) next
  all_counts <- rbindlist(chain_data[[chain_type_key]], fill = TRUE)
  aa_matrix <- dcast(all_counts, stage + position + person ~ aa, value.var = "freq", fill = 0)
  output_filename <- paste0("combined_AA_freq_", chain_type_key, "_T.csv")
  fwrite(aa_matrix, output_filename)
  combined_aa_freq_files[[chain_type_key]] <- output_filename
}

# --- 2. Run Wilcoxon Tests ---

# Function to compare stages across all people
compare_stages_all_people <- function(dt) {
  results <- list()
  stages <- sort(unique(dt$stage))
  positions <- sort(unique(dt$position))
  aa_cols <- setdiff(colnames(dt), c("stage", "position", "person"))

  for (pos in positions) {
    for (aa in aa_cols) {
      for (i in 1:(length(stages) - 1)) {
        stg1 <- stages[i]
        stg2 <- stages[i + 1]
        group1 <- dt[stage == stg1 & position == pos, .(person, value = get(aa))]
        group2 <- dt[stage == stg2 & position == pos, .(person, value = get(aa))]
        merged <- merge(group1, group2, by = "person", suffixes = c("_s1", "_s2"))
        
        if (nrow(merged) >= 3) {
          test <- tryCatch({ wilcox.test(merged$value_s1, merged$value_s2, paired = TRUE) },
                           error = function(e) { list(p.value = NA) })
          results[[length(results) + 1]] <- data.table(
            position = pos, aa = aa, stage1 = stg1, stage2 = stg2, pval = signif(test$p.value, 3), n = nrow(merged)
          )
        }
      }
    }
  }
  if (length(results) > 0) return(rbindlist(results)) else return(data.table())
}

# Load freq data and run tests for each chain
chains_freq_data <- list()
if ("H" %in% names(combined_aa_freq_files)) chains_freq_data[["H"]] <- fread(combined_aa_freq_files[["H"]])
if ("K" %in% names(combined_aa_freq_files)) chains_freq_data[["K"]] <- fread(combined_aa_freq_files[["K"]])
if ("L" %in% names(combined_aa_freq_files)) chains_freq_data[["L"]] <- fread(combined_aa_freq_files[["L"]])

wilcox_results_list <- list()
for (chain_name in names(chains_freq_data)) {
  wilcox_dt <- compare_stages_all_people(chains_freq_data[[chain_name]])
  if (nrow(wilcox_dt) > 0) {
    wilcox_dt_wide <- dcast(wilcox_dt, position + stage1 + stage2 ~ aa, value.var = "pval")
  } else {
    wilcox_dt_wide <- data.table(position = numeric(), stage1 = numeric(), stage2 = numeric())
  }
  wilcox_results_list[[chain_name]] <- wilcox_dt_wide
  fwrite(wilcox_dt_wide, paste0('wilcox_pos_', chain_name, '.csv'))
}

wilcox_H <- wilcox_results_list[["H"]]
wilcox_K <- wilcox_results_list[["K"]]
wilcox_L <- wilcox_results_list[["L"]]

# --- 3. Plotting ---

# Common themes and color scales
light_theme <- theme_minimal(base_size = 11) + theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"), axis.text = element_text(size = 8), legend.position = "bottom")
blue_fill <- scale_fill_gradient(low = "#e0f2ff", high = "#045a8d", name = "1 / p-value")

# Plotting function
generate_heatmap <- function(dt, chain_name, title_prefix, fill_scale) {
  if (is.null(dt) || nrow(dt) == 0) { message(paste0("No data for ", chain_name, " plot.")); return(NULL) }
  
  melted <- melt(dt, id.vars = c("position", "stage1", "stage2"), measure.vars = patterns("^[A-Z]$"), variable.name = "aa", value.name = "pval")
  melted[, aa := as.character(aa)]
  melted[, inv_pval := pmin(1 / pval, 100)]
  melted[, transition_label := paste0("Stage ", stage1, " → ", stage2)]
  melted[, transition_label := factor(transition_label, levels = unique(melted[order(stage2)]$transition_label))]
  
  intensity <- melted[, .(mean_inv_pval = mean(inv_pval, na.rm = TRUE)), by = transition_label]
  ordered_transitions <- intensity[order(-mean_inv_pval), transition_label]
  melted[, transition_label := factor(transition_label, levels = ordered_transitions)]
  melted[, aa_ordered := reorder_within(aa, inv_pval, transition_label)]

  p1 <- ggplot(melted, aes(x = position, y = aa, fill = inv_pval)) + geom_tile() + fill_scale + facet_wrap(~ transition_label, ncol = 3) + 
        labs(title = "Alphabetical Order", x = "CDR3 Position", y = "Amino Acid") + light_theme
  
  p2 <- ggplot(melted, aes(x = position, y = aa_ordered, fill = inv_pval)) + geom_tile() + fill_scale + facet_wrap(~ transition_label, ncol = 3, scales = "free_y") + 
        scale_y_reordered() + labs(title = "By Significance", x = "CDR3 Position", y = "Amino Acid (per transition)") + light_theme

  final_plot <- (p1 / p2) + plot_layout(heights = c(1, 1.1)) + plot_annotation(title = paste0(title_prefix, " (", chain_name, " Chain)"), theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
  return(final_plot)
}

# Make and save plots for H, K, L chains
plot_H <- generate_heatmap(wilcox_H, "H", "Wilcoxon p-value Heatmaps Across Stage Transitions", blue_fill)
if (!is.null(plot_H)) ggsave("wilcoxon_H_heatmap.png", plot_H, width = 11, height = 8, dpi = 300)

plot_K <- generate_heatmap(wilcox_K, "K", "Wilcoxon p-value Heatmaps Across Stage Transitions", blue_fill)
if (!is.null(plot_K)) ggsave("wilcoxon_K_heatmap.png", plot_K, width = 11, height = 8, dpi = 300)

plot_L <- generate_heatmap(wilcox_L, "L", "Wilcoxon p-value Heatmaps Across Stage Transitions", blue_fill)
if (!is.null(plot_L)) ggsave("wilcoxon_L_heatmap.png", plot_L, width = 11, height = 8, dpi = 300)

cat("Analysis complete. Check your directory for output files and plots.\n")