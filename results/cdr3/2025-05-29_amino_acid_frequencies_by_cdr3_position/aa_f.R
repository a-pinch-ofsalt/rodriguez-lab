library(data.table)
library(ggplot2)
library(gridExtra)

# Load long-format amino_acid frequency data
df <- fread("AA_freq_full_long_by_stage.csv")

# Add average across people
avg <- df[, .(freq = mean(freq, na.rm = TRUE)), by = .(chain, position, stage, amino_acid)]
avg[, person := "AVG"]
df_all <- rbind(df, avg)

# Compute frequency diffs between consecutive stages
get_diffs <- function(dt) {
  setorder(dt, stage)
  dt[, freq_diff := freq - shift(freq), by = .(chain, position, amino_acid, person)]
  dt[!is.na(freq_diff)]
}

# Get global scale range
delta_all <- get_diffs(df_all)
# Compute max for just this chain
# Compute max for just this chain

# Plot function per chain
plot_deltas <- function(chain_type) {
  subset <- delta_all[chain == chain_type]
  people <- unique(subset$person)
  
  # These are stages that contain diffs (i.e., stage 2 holds diff 1→2)
  stage_transitions <- sort(unique(subset[!is.na(freq_diff)]$stage))
  
  scale_limit <- max(abs(subset$freq_diff), na.rm = TRUE)
  
  plots <- list()
  for (p in people) {
    for (s in stage_transitions) {
      dt <- subset[person == p & stage == s]
      mat <- dcast(dt, position ~ amino_acid, value.var = "freq_diff", fill = 0)
      molten <- melt(mat, id.vars = "position", variable.name = "amino_acid", value.name = "freq_diff")
      
      g <- ggplot(molten, aes(x = amino_acid, y = factor(position), fill = freq_diff)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(
          low = "blue", mid = "white", high = "red", midpoint = 0,
          limits = c(-scale_limit, scale_limit)
        ) +
        labs(title = paste(p, "- ΔStage", s - 1, "→", s), x = "amino_acid", y = "Position") +
        theme_minimal(base_size = 9) +
        theme(axis.text.x = element_text(angle = 90, size = 6),
              axis.text.y = element_text(size = 6),
              plot.title = element_text(size = 8),
              axis.title = element_blank())
      
      plots[[paste(p, s)]] <- g
    }
  }
  
  n_rows <- length(stage_transitions)
  grid <- marrangeGrob(plots, nrow = n_rows, ncol = 7,
                       top = paste("Δ amino_acid Frequencies – Chain", chain_type))
  ggsave(paste0("delta_amino_acid_freq_grid_chain_", chain_type, ".png"), grid,
         width = 20, height = 7)
}


# Generate plots per chain
for (ch in c("H", "K", "L")) {
  plot_deltas(ch)
}





plot_deltas <- function(chain_type) {
  subset <- delta_all[chain == chain_type & person == "AVG"]
  
  # Just include stage transitions that actually have diffs
  stage_transitions <- sort(unique(subset[!is.na(freq_diff)]$stage))
  
  scale_limit <- max(abs(subset$freq_diff), na.rm = TRUE)
  
  plots <- list()
  for (s in stage_transitions) {
    dt <- subset[stage == s]
    mat <- dcast(dt, position ~ amino_acid, value.var = "freq_diff", fill = 0)
    molten <- melt(mat, id.vars = "position", variable.name = "amino_acid", value.name = "freq_diff")
    
    g <- ggplot(molten, aes(x = amino_acid, y = factor(position), fill = freq_diff)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(
        low = "blue", mid = "white", high = "red", midpoint = 0,
        limits = c(-scale_limit, scale_limit)
      ) +
      labs(title = paste("AVG - ΔStage", s - 1, "→", s), x = "amino_acid", y = "Position") +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 90, size = 6),
            axis.text.y = element_text(size = 6),
            plot.title = element_text(size = 8),
            axis.title = element_blank())
    
    plots[[paste("AVG", s)]] <- g
  }
  
  n_rows <- length(stage_transitions)
  grid <- marrangeGrob(plots, nrow = n_rows, ncol = 1,
                       top = paste("Δ amino_acid Frequencies (AVG) – Chain", chain_type))
  ggsave(paste0("delta_amino_acid_freq_grid_chain_", chain_type, "_AVG_only.png"),
         grid, width = 4, height = 10)
}
for (ch in c("H", "K", "L")) {
  plot_deltas(ch)
}

