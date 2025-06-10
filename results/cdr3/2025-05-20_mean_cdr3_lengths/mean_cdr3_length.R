#gpt did the write-csv part

library(data.table)

dirs <- c("A", "B", "C", "D", "E", "G")
files <- c("H_T", "K_T", "L_T")

# Prepare an empty list to collect results
results <- list()

for (dir in dirs) {
  for (file in files) {
    path <- file.path(dir, paste0(file, ".tsv"))
    if (!file.exists(path)) {
      warning(paste("File not found:", path))
      next
    }
    dat <- fread(path)
    mean_length <- mean(nchar(dat$cdr3_aa))
    print()
    # print(file)
    # 
    # mean_length <- mean(dat$cdr3_end[valid_rows] - dat$cdr3_start[valid_rows])
    # print(mean_length)
    
    print()
    results[[length(results) + 1]] <- data.table(
      individual = dir,
      c("chain", stage)
      Directory = dir,
      
      File = file,
      MeanCDR3Length = mean_length
    )
  }
}

# Combine all results into one data table
all_results <- rbindlist(results)

# Write to CSV
fwrite(all_results, "cdr3_means.csv")

print("Saved means to cdr3_means.csv")
