library(data.table)
library(ggplot2)
library(ggrepel)

# ── Settings ─────────────────────────────────────────────
input_dir <- "stage_kmer_counts"
output_dir <- "kmer_stage_wilcoxon"
dir.create(output_dir, showWarnings = FALSE)

chains <- c("L", "K", "H")
top_k <- 1000

# ── Main loop per chain ─────────────────────────────────
for (chain in chains) {
  cat("Processing chain:", chain, "\n")
  
  files <- list.files(input_dir, pattern = paste0("_", chain, "_T.tsv_cdr3_only_stage_\\d+_kmer_counts\\.csv$"), full.names = TRUE)
  meta <- data.table(file = files)
  meta[, fname := basename(file)]
  meta[, individual := sub("_.*", "", fname)]
  meta[, stage := as.integer(sub(".*stage_(\\d+)_kmer_counts\\.csv$", "\\1", fname))]
  
  # Get top kmers from each file
  top_list <- rbindlist(lapply(1:nrow(meta), function(i) {
    dt <- fread(meta[i]$file, select = c("kmer", "count"), nrows = top_k)
    dt[, `:=`(individual = meta[i]$individual, stage = meta[i]$stage)]
    return(dt)
  }))
  
  # Sum across to get overall top 1000 kmers
  top_summary <- top_list[, .(total_count = sum(count)), by = kmer][order(-total_count)][1:top_k]
  top_kmers <- top_summary$kmer
  top_list <- top_list[kmer %in% top_kmers]
  
  # Pivot to wide format: kmer × (individual_stage)
  top_wide <- dcast(top_list, kmer ~ paste0(individual, "_stage", stage), value.var = "count", fill = 0)
  
  # Get list of stages present
  stages <- sort(unique(meta$stage))
  
  # Prepare results table
  result_table <- copy(top_wide[, .(kmer)])
  
  # Add log2FC and p-values for every stage pair
  for (i in 1:(length(stages) - 1)) {
    for (j in (i + 1):length(stages)) {
      s1 <- stages[i]
      s2 <- stages[j]
      cols1 <- grep(paste0("_stage", s1, "$"), names(top_wide), value = TRUE)
      cols2 <- grep(paste0("_stage", s2, "$"), names(top_wide), value = TRUE)
      
      # Wilcoxon per row
      test_res <- apply(top_wide[, -1, with = FALSE], 1, function(row) {
        v1 <- as.numeric(row[cols1])
        v2 <- as.numeric(row[cols2])
        if (length(v1) < 2 || length(v2) < 2) return(c(NA, NA))
        pval <- tryCatch(wilcox.test(v1, v2)$p.value, error = function(e) NA)
        log2fc <- log2(mean(v2) + 1e-9) - log2(mean(v1) + 1e-9)
        c(log2fc, pval)
      })
      
      test_res <- t(test_res)
      result_table[[paste0("log2fc_", s1, "_vs_", s2)]] <- test_res[, 1]
      result_table[[paste0("pval_", s1, "_vs_", s2)]] <- test_res[, 2]
    }
  }
  
  # Save main result
  fwrite(result_table, file.path(output_dir, paste0("kmer_wilcoxon_table_", chain, ".csv")))
  
  # ── Volcano plots for each stage pair ────────────────
  stage_pairs <- grep("^log2fc_", names(result_table), value = TRUE)
  
  for (fc_col in stage_pairs) {
    pval_col <- gsub("log2fc", "pval", fc_col)
    tag <- gsub("log2fc_", "", fc_col)
    
    volcano <- result_table[, .(kmer, log2fc = get(fc_col), pval = get(pval_col))]
    volcano <- volcano[!is.na(pval) & !is.na(log2fc)]
    volcano[, `:=`(
      neglog10p = -log10(pval + 1e-20),
      sqrt_fc = sign(log2fc) * sqrt(abs(log2fc)),
      dot_id = .I
    )]
    
    p <- ggplot(volcano, aes(x = sqrt_fc, y = neglog10p, color = dot_id)) +
      geom_point(size = 2, alpha = 0.85) +
      scale_color_gradient(low = "black", high = "red", name = "Row Index") +
      geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "darkred") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      theme_minimal() +
      labs(
        title = paste("Volcano:", chain, tag),
        x = "Sqrt(Log2 Fold Change)", y = "-log10 P-Value"
      )
    
    ggsave(file.path(output_dir, paste0("volcano_", chain, "_", tag, ".png")),
           p, width = 8, height = 6, dpi = 300)
  }
}

